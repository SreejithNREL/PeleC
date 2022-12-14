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

  const amrex::Real kernel_y_loc    = 1.92;
  const amrex::Real kernel_height   = 1.01e-01;
  amrex::Real kernel_location[3];
  amrex::Real r;

  amrex::Real y_centered = y - prob_parm.y_center;
  amrex::Real u[3] = {0.0};
  amrex::Real rho;

  u[0] = state(i, j, k, UMX);
  u[1] = state(i, j, k, UMY);
  u[2] = state(i, j, k, UMZ);

  kernel_location[0] = 0.0;
  kernel_location[1] = kernel_y_loc-kernel_height/2.0;
  kernel_location[2] = 0.0;

  r = sqrt(pow((kernel_location[0]-x),2)+pow((kernel_location[1]-y),2)+pow((kernel_location[2]-z),2));

  // if(y_centered > -prob_parm.pure_h2_width and y_centered < prob_parm.pure_h2_width and
  //             z > -prob_parm.pure_h2_width and          z < prob_parm.pure_h2_width){

  // if(r < prob_parm.pure_h2_width){
  //   // printf("Initializing pure H2 composition\n");

  //   for (int n = 0; n < NUM_SPECIES; n++)
  //    state(i, j, k, UFS + n) = prob_parm.H2_state[UFS + n];

  //   rho = prob_parm.H2_state[URHO];

  //   //initial velocity is coming from plt file

  //   state(i, j, k, URHO)  = prob_parm.H2_state[URHO];
  //   state(i, j, k, UEINT) = prob_parm.H2_state[UEINT];
  //   state(i, j, k, UTEMP) = prob_parm.H2_state[UTEMP];
  //   state(i, j, k, UEDEN) = rho * (prob_parm.H2_state[UEINT]/rho + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2])); 
  // }
  // else
  // {
    // printf("Initializing H2/ammonia composition\n");

    for (int n = 0; n < NUM_SPECIES; n++)
     state(i, j, k, UFS + n) = prob_parm.fuel_state[UFS + n];

    rho = prob_parm.fuel_state[URHO];

    //initial velocity is coming from plt file

    state(i, j, k, UMX) = u[0]*rho;
    state(i, j, k, UMY) = u[1]*rho;
    state(i, j, k, UMZ) = u[2]*rho;

    state(i, j, k, URHO)  = prob_parm.fuel_state[URHO];
    state(i, j, k, UEINT) = prob_parm.fuel_state[UEINT];
    state(i, j, k, UTEMP) = prob_parm.fuel_state[UTEMP];
    state(i, j, k, UEDEN) = rho * (prob_parm.fuel_state[UEINT]/rho + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));     
    state(i, j, k, UFA)   = prob_parm.fuel_state[UFA];
  // }
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
  const amrex::Real* dx = geomdata.CellSize();
  amrex::Real velx, vely, velz;
  amrex::Real rho, T, pres;

  amrex::Real y_centered = x[1] - prob_parm.y_center;

  // if(sgn == 1){
    if(prob_parm.turbulence){

      velx = s_ext[UMX] + prob_parm.vn_in;
      vely = s_ext[UMY];
      velz = s_ext[UMZ];      
    }
    else{
      velx = 0.0 + prob_parm.vn_in;
      vely = 0.0;
      velz = 0.0;
    }

    // if(y_centered > -prob_parm.pure_h2_width and y_centered < prob_parm.pure_h2_width and
    //          x[2] > -prob_parm.pure_h2_width and       x[2] < prob_parm.pure_h2_width and
    //           time < prob_parm.time_h2_injection){

      // for (int n = 0; n < NUM_SPECIES; n++)
      //  s_ext[UFS+ n] = prob_parm.H2_state[UFS + n];

      // rho = prob_parm.H2_state[URHO];

      // //initial velocity is coming from plt file
      // s_ext[UMX]   = velx*rho;
      // s_ext[UMY]   = vely*rho;
      // s_ext[UMZ]   = velz*rho;
      // s_ext[URHO]  = prob_parm.H2_state[URHO];
      // s_ext[UEINT] = prob_parm.H2_state[UEINT];
      // s_ext[UTEMP] = prob_parm.H2_state[UTEMP];
      // s_ext[UEDEN]  = rho * (prob_parm.H2_state[UEINT]/rho + 0.5 * (velx * velx + vely * vely + velz * velz)); 

    // }
    // else
    // {
    
      for (int n = 0; n < NUM_SPECIES; n++)
       s_ext[UFS+ n] = prob_parm.fuel_state[UFS + n];

      rho = prob_parm.fuel_state[URHO];

      //initial velocity is coming from plt file
      s_ext[UMX]   = velx*rho;
      s_ext[UMY]   = vely*rho;
      s_ext[UMZ]   = velz*rho;
      s_ext[URHO]  = prob_parm.fuel_state[URHO];
      s_ext[UEINT] = prob_parm.fuel_state[UEINT];
      s_ext[UTEMP] = prob_parm.fuel_state[UTEMP];
      s_ext[UEDEN]  = rho * (prob_parm.fuel_state[UEINT]/rho + 0.5 * (velx * velx + vely * vely + velz * velz));     

    // }
  // }
  // else{
  //   amrex::Real velx = 0.0;
  //   amrex::Real vely = 0.0;
  //   amrex::Real velz = 0.0;
  //   s_ext[URHO]  = prob_parm.fuel_state[URHO];
  //   s_ext[UMX]   = prob_parm.fuel_state[URHO]*velx;
  //   s_ext[UMY]   = prob_parm.fuel_state[URHO]*vely;
  //   s_ext[UMZ]   = prob_parm.fuel_state[URHO]*velz;
  //   s_ext[UEINT] = prob_parm.fuel_state[UEINT];
  //   s_ext[UEDEN] = prob_parm.fuel_state[UEINT] + prob_parm.fuel_state[URHO] * (0.5 * (velx * velx + vely * vely + velz * velz));
  //   s_ext[UTEMP] = prob_parm.fuel_state[UTEMP];
  //   for (int n = 0; n < NUM_SPECIES; n++)
  //     s_ext[UFS + n] = prob_parm.fuel_state[UFS + n];
  //     // s_ext[UFS + n] = prob_parm.H2_state[UFS + n];
  // }
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

void
init_composition_ammonia(amrex::Real molefrac[NUM_SPECIES]);

void
init_composition_air(amrex::Real molefrac[NUM_SPECIES]);

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
