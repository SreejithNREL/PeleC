#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_Print.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_ParmParse.H>

#include "mechanism.H"

#include "PeleC.H"
#include "IndexDefines.H"
#include "PelePhysics.H"
#include "Tagging.H"
#include "ProblemDerive.H"
#include "prob_parm.H"
#include "Forcing.H"
#include "Utilities.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
  amrex::Real u[3] = {0.0};
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  amrex::Real rho, energy, T, pres, e;
  amrex::Real pert = 0.0;

  amrex::Real kernel_location[3];
  const amrex::Real kernel_diameter = 100.0e-04;
  const amrex::Real kernel_y_loc    = 1.92;
  const amrex::Real kernel_height   = 1.01e-01;
  double r;

#ifndef AMREX_USE_GPU
  kernel_location[0] = 0.0;
  kernel_location[1] = kernel_y_loc-kernel_height/2.0;
  kernel_location[2] = 0.0;

  if(prob_parm.init_kernel){

    if (prob_parm.pertmag > 0.0) {
      pert =
        prob_parm.pertmag *
        (1.0 * std::sin(2 * constants::PI() * 4 * x / 0.1) *
           std::sin(2 * constants::PI() * 5 * z / 0.1) +
         1.023 *
           std::sin(2 * constants::PI() * 2 * (x - 0.4598) / 0.1) *
           std::sin(2 * constants::PI() * 4 * (z - 0.53765) / 0.1) +
         0.945 *
           std::sin(2 * constants::PI() * 3 * (x - 0.712435) / 0.1) *
           std::sin(2 * constants::PI() * 3 * (z - 2.137) / 0.1) +
         1.017 *
           std::sin(2 * constants::PI() * 5 * (x - 0.33) / 0.1) *
           std::sin(2 * constants::PI() * 6 * (z - 1.8) / 0.1) +
         0.982 *
           std::sin(2 * constants::PI() * 5 * (x - 1.4234) / 0.1));
    }
    
    r = sqrt(pow((kernel_location[0]-x),2)+pow((kernel_location[2]-z),2));
    if(y > prob_parm.d_yarray[0] and y < prob_parm.d_yarray[prob_parm.ny-1] 
      and r < prob_parm.d_xarray[prob_parm.nx-1]){

      // Fill in the velocities and interpolated quantities.
      amrex::Real u[3] = {0.0};
      amrex::Real uinterp[prob_parm.nscal];

      // Interpolation factors
      amrex::Real mod[3] = {0.0};
      int idx   = 0;
      int idxp1 = 0;
      int idy   = 0;
      int idyp1 = 0;
      amrex::Real slp[3] = {0.0};
      
      amrex::Real radius_pert = r+pert;
      if(radius_pert >= prob_parm.d_xarray[prob_parm.nx-1]) //this avoids extrapolation
        radius_pert = prob_parm.d_xarray[prob_parm.nx-1];

      mod[0] = std::fmod(radius_pert, prob_parm.Lxinput);
      locate(prob_parm.d_xarray, prob_parm.nx, mod[0], idx);
      idxp1 = (idx + 1) % prob_parm.nx;
      //idxp1 = (idx + 1);
      slp[0] = (mod[0] - prob_parm.d_xarray[idx]) / prob_parm.d_xdiff[idx];

      mod[1] = std::fmod(y, prob_parm.Lyinput);
      locate(prob_parm.d_yarray, prob_parm.ny, mod[1], idy);
      idyp1 = (idy + 1) % prob_parm.ny;
      //idyp1 = (idy + 1);
      slp[1] = (mod[1] - prob_parm.d_yarray[idy]) / prob_parm.d_ydiff[idy];

      const amrex::Real f0 = (1 - slp[0]) * (1 - slp[1]);
      const amrex::Real f1 = slp[0] * (1 - slp[1]);
      const amrex::Real f2 = (1 - slp[0]) * slp[1];
      const amrex::Real f3 = slp[0] * slp[1];

      // Interpolate data
      for (int iscal=0;iscal<prob_parm.nscal;iscal++){
        uinterp[iscal] = prob_parm.d_data_ic[idx  +idy  *prob_parm.nx+iscal*prob_parm.nx*prob_parm.ny]*f0 +
                         prob_parm.d_data_ic[idxp1+idy  *prob_parm.nx+iscal*prob_parm.nx*prob_parm.ny]*f1 +
                         prob_parm.d_data_ic[idx  +idyp1*prob_parm.nx+iscal*prob_parm.nx*prob_parm.ny]*f2 +
                         prob_parm.d_data_ic[idxp1+idyp1*prob_parm.nx+iscal*prob_parm.nx*prob_parm.ny]*f3;
      
      }

      // // temperature, pressure, internal energy and composition
      double p ;//= uinterp[1]*10.0; //Pa to dyn/cm2
      rho = uinterp[0]/1.0e+3;
      T = uinterp[2];
      double massfrac[NUM_SPECIES]={0.0};

     for (int n = 0; n < NUM_SPECIES; n++)
       massfrac[n] = uinterp[n+6];

     double sum = 0.0;
     for (int n = 0; n < NUM_SPECIES; n++)
       sum = sum + massfrac[n];
      // printf("Sum of all species = % e \n", sum);
      massfrac[N2_ID] = massfrac[N2_ID] + 1.0 - sum; //making sure that sum of ys is 1
      // massfrac[N2_ID] = 0.767;
      // massfrac[O2_ID] = 0.233;

      auto eos = pele::physics::PhysicsType::eos();
      eos.RTY2P(rho, T, massfrac, p);
      eos.PYT2RE(p, massfrac, T, rho, e);

      double beta  = atan(z/x);
      double alpha = atan(x/z);

      //Velocity in each direction
      if(x >= 0.0){
        u[0] = uinterp[5]*cos(beta)*100.0 + prob_parm.vn_in;
      }
      else{
        u[0] = -uinterp[5]*cos(beta)*100.0 + prob_parm.vn_in;
      }

      u[1] = uinterp[4]*100.0;

      if(z >= 0){
        u[2] = uinterp[5]*cos(alpha)*100.0;
      }
      else{
        u[2] = -uinterp[5]*cos(alpha)*100.0;
      }
      for (int n = 0; n < NUM_SPECIES; n++)
        state(i, j, k, UFS + n) = rho * massfrac[n];

      state(i, j, k, URHO)  = rho;
      state(i, j, k, UEINT) = rho * e;
      state(i, j, k, UTEMP) = T;
      state(i, j, k, UMX)   = rho * u[0];
      state(i, j, k, UMY)   = rho * u[1];
      state(i, j, k, UMZ)   = rho * u[2];
      state(i, j, k, UEDEN)  = rho * (e + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2])); 
    }
    else{ 
      for (int n = 0; n < NVAR; n++)
        state(i, j, k, n) = prob_parm.fuel_state[n];
    }
  }
  //if init_kernel is false: initialize premixed composition
  else{ 
    for (int n = 0; n < NVAR; n++)
      state(i, j, k, n) = prob_parm.fuel_state[n];
  }

  if(prob_parm.turbulence){

    amrex::Real uinterp[3] = {0.0};

    // Interpolation factors
    amrex::Real mod[3] = {0.0};
    int idx[3] = {0};
    int idxp1[3] = {0};
    amrex::Real slp[3] = {0.0};


    // mod[0] = std::fmod(x+1.0, prob_parm.Lxturb); //this line is to start from S3D
    mod[0] = std::fmod(x+0.5, prob_parm.Lxturb);
    locate(prob_parm.d_xarray_turb, prob_parm.ires_x, mod[0], idx[0]);
    idxp1[0] = (idx[0] + 1) % prob_parm.ires_x;
    slp[0] = (mod[0] - prob_parm.d_xarray_turb[idx[0]]) / prob_parm.d_xdiff_turb[idx[0]];

    // mod[1] = std::fmod(y-0.92, prob_parm.Lxturb); //this line is to start from S3D
    mod[1] = std::fmod(y, prob_parm.Lyturb);
    locate(prob_parm.d_yarray_turb, prob_parm.ires_y, mod[1], idx[1]);
    idxp1[1] = (idx[1] + 1) % prob_parm.ires_y;
    slp[1] = (mod[1] - prob_parm.d_yarray_turb[idx[1]]) / prob_parm.d_ydiff_turb[idx[1]];

    // mod[2] = std::fmod(z+1.0, prob_parm.Lxturb); //this line is to start from S3D
    mod[2] = std::fmod(z+0.4, prob_parm.Lzturb);
    locate(prob_parm.d_zarray_turb, prob_parm.ires_z, mod[2], idx[2]);
    idxp1[2] = (idx[2] + 1) % prob_parm.ires_z;
    slp[2] = (mod[2] - prob_parm.d_zarray_turb[idx[2]]) / prob_parm.d_zdiff_turb[idx[2]];

    const amrex::Real f0 = (1 - slp[0]) * (1 - slp[1]) * (1 - slp[2]);
    const amrex::Real f1 = slp[0] * (1 - slp[1]) * (1 - slp[2]);
    const amrex::Real f2 = (1 - slp[0]) * slp[1] * (1 - slp[2]);
    const amrex::Real f3 = (1 - slp[0]) * (1 - slp[1]) * slp[2];
    const amrex::Real f4 = slp[0] * (1 - slp[1]) * slp[2];
    const amrex::Real f5 = (1 - slp[0]) * slp[1] * slp[2];
    const amrex::Real f6 = slp[0] * slp[1] * (1 - slp[2]);
    const amrex::Real f7 = slp[0] * slp[1] * slp[2];
    int nx = prob_parm.ires_x;
    int ny = prob_parm.ires_y;
    int nz = prob_parm.ires_z;


    // Interpolate data
    for (int iscal=0;iscal<3;iscal++){

      uinterp[iscal] = prob_parm.d_data_ic_turb[idx[0]   + nx*idx[1]   + nx*ny*idx[2]   + iscal*nx*ny*nz]*f0 +
                       prob_parm.d_data_ic_turb[idxp1[0] + nx*idx[1]   + nx*ny*idx[2]   + iscal*nx*ny*nz]*f1 +
                       prob_parm.d_data_ic_turb[idx[0]   + nx*idxp1[1] + nx*ny*idx[2]   + iscal*nx*ny*nz]*f2 +
                       prob_parm.d_data_ic_turb[idx[0]   + nx*idx[1]   + nx*ny*idxp1[2] + iscal*nx*ny*nz]*f3 +
                       prob_parm.d_data_ic_turb[idxp1[0] + nx*idx[1]   + nx*ny*idxp1[2] + iscal*nx*ny*nz]*f4 +
                       prob_parm.d_data_ic_turb[idx[0]   + nx*idxp1[1] + nx*ny*idxp1[2] + iscal*nx*ny*nz]*f5 +
                       prob_parm.d_data_ic_turb[idxp1[0] + nx*idxp1[1] + nx*ny*idx[2]   + iscal*nx*ny*nz]*f6 +
                       prob_parm.d_data_ic_turb[idxp1[0] + nx*idxp1[1] + nx*ny*idxp1[2] + iscal*nx*ny*nz]*f7 ;
    }

    rho = state(i, j, k, URHO);
    e   = state(i, j, k, UEINT) / rho;

    amrex::Real velx = state(i, j, k, UMX)/rho + uinterp[0];
    amrex::Real vely = state(i, j, k, UMY)/rho + uinterp[1];
    amrex::Real velz = state(i, j, k, UMZ)/rho + uinterp[2];

    state(i, j, k, UMX) = rho * velx;
    state(i, j, k, UMY) = rho * vely;
    state(i, j, k, UMZ) = rho * velz;
    state(i, j, k, UEDEN)  = rho * (e + 0.5 * (velx * velx + vely * vely + velz * velz)); 
  }
#endif
}


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  amrex::Real u[3] = {0.0};
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  amrex::Real rho, energy, T, pres;
  amrex::Real l_fluct = 1.24E-02*2.0; //corresponds to lt for the turbulence level we are aiming at 
  amrex::Real pert = prob_parm.vn_in * 0.1; //velocity perturbation at the inlet. For now it is 10 % the mean value

  s_ext[URHO]  = prob_parm.fuel_state[URHO];
 // s_ext[UMX]   = prob_parm.fuel_state[URHO]*(sin(x[1]/l_fluct*3.14156)*pert + prob_parm.vn_in);
  s_ext[UMX]   = prob_parm.fuel_state[UMX];
  s_ext[UMY]   = prob_parm.fuel_state[UMY];
  s_ext[UMZ]   = prob_parm.fuel_state[UMZ];
  s_ext[UEINT] = prob_parm.fuel_state[UEINT];
  s_ext[UEDEN] = prob_parm.fuel_state[UEDEN];
  s_ext[UTEMP] = prob_parm.fuel_state[UTEMP];
  for (int n = 0; n < NUM_SPECIES; n++)
    s_ext[UFS + n] = prob_parm.fuel_state[UFS + n];
}

void read_input_file(
  const std::string iname,
  int z_coord,
  int nx,
  int ny,
  int nz,
  int nscal,
  amrex::Vector<amrex::Real>& gridx_input,
  amrex::Vector<amrex::Real>& gridy_input,
  amrex::Vector<amrex::Real>& gridz_input,
  amrex::Vector<amrex::Real>& data);

void get_inputs(
  const std::string& iname,
   int z_coord);

void pc_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
