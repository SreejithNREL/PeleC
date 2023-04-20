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

  for (int n = 0; n < NVAR; n++)
   state(i, j, k, n) = prob_parm.fuel_state[n];

  state(i, j, k, UFA) = 0.0;

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

  if(sgn == 0){
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

    rho = prob_parm.fuel_state[URHO];

    s_ext[URHO]  = rho;
    s_ext[UMX]   = rho*velx;
    s_ext[UMY]   = rho*vely;
    s_ext[UMZ]   = rho*velz;
    s_ext[UEINT] = prob_parm.fuel_state[UEINT];
    s_ext[UEDEN] = prob_parm.fuel_state[UEINT] + rho * (0.5 * (velx * velx + vely * vely + velz * velz));
    s_ext[UTEMP] = prob_parm.fuel_state[UTEMP];
    for (int n = 0; n < NUM_SPECIES; n++)
      s_ext[UFS + n] = prob_parm.fuel_state[UFS + n];
    
    // printf("End of turbulent BC \n");
  }
  else{
    amrex::Real velx = 0.0;
    amrex::Real vely = 0.0;
    amrex::Real velz = 0.0;
    s_ext[URHO]  = prob_parm.fuel_state[URHO];
    s_ext[UMX]   = prob_parm.fuel_state[URHO]*velx;
    s_ext[UMY]   = prob_parm.fuel_state[URHO]*vely;
    s_ext[UMZ]   = prob_parm.fuel_state[URHO]*velz;
    s_ext[UEINT] = prob_parm.fuel_state[UEINT];
    s_ext[UEDEN] = prob_parm.fuel_state[UEINT] + prob_parm.fuel_state[URHO] * (0.5 * (velx * velx + vely * vely + velz * velz));
    s_ext[UTEMP] = prob_parm.fuel_state[UTEMP];
    for (int n = 0; n < NUM_SPECIES; n++)
      s_ext[UFS + n] = prob_parm.fuel_state[UFS + n];
      // s_ext[UFS + n] = prob_parm.fuel_state[UFS + n];
  }
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
init_composition(amrex::Real molefrac[NUM_SPECIES]);

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
