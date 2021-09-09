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

  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  // const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  amrex::Real x[3] = {
    prob_lo[0] + static_cast<amrex::Real>(i + 0.5) * dx[0],
    prob_lo[1] + static_cast<amrex::Real>(j + 0.5) * dx[1],
    prob_lo[2] + static_cast<amrex::Real>(k + 0.5) * dx[2]};

  amrex::Real p = PeleC::h_prob_parm_device->pamb;
  amrex::Real T = PeleC::h_prob_parm_device->T_in;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  amrex::Real rho, e, vt;
  vt = 0.0;
  amrex::Real ek = 0.5 * (vt * vt);

  for (int n = 0; n < NUM_SPECIES; n++)
    massfrac[n] = 1.0;


  amrex::Real mean          = 0.0;
  amrex::Real sigma         = 0.01;
  amrex::Real peak_pressure = p*0.2;
  amrex::Real pi            = 3.14156;
  amrex::Real gaussian;
  
  // gaussian_function(mean,mean,sigma,gaussian);
  amrex::Real peak_gaussian = 1/sigma/(2.*pi)*exp(-0.5*pow((mean - mean)/sigma, 2));

  // gaussian_function(x,mean,sigma,gaussian);
  gaussian =  1/sigma/(2.*pi)*exp(-0.5*pow((x[0] - mean)/sigma, 2));
  amrex::Real P_local = p + gaussian/peak_gaussian * peak_pressure;  
  // printf("Pressure %f \n", P_local);


  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(P_local, massfrac, T, rho, e);


  state(i, j, k, URHO)  = rho;
  state(i, j, k, UMX)   = 0.0;
  state(i, j, k, UMY)   = 0.0;
  state(i, j, k, UMZ)   = 0.0;
  state(i, j, k, UEINT) = rho * e;
  state(i, j, k, UEDEN) = rho * (e + ek);
  state(i, j, k, UTEMP) = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k,UFS + n) = rho * massfrac[n];
  }
}


void
gaussian_function(amrex::Real x,amrex::Real mean,amrex::Real sigma, amrex::Real gaussian);


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const amrex::Real* /*s_int[NVAR]*/,
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real /*time*/,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm,
  int bc_type = 0,
  amrex::Real* bc_params = nullptr,
  amrex::Real* bc_target = nullptr)
{
  if(idir == 0){
    for (int n = 0; n < NVAR; n++){
      s_ext[n] = prob_parm.air_state[n];
      // printf("s_ext[%i] = %f\n", n, s_ext[n]);
    }

    if(bc_type = -1){
      bc_type = 8; //outflow

      amrex::Real massfrac[NUM_SPECIES] = {0.0};
      amrex::Real P, rho;
      for (int n = 0; n < NUM_SPECIES; n++)
        massfrac[n] = prob_parm.air_state[UFS + n];

      auto eos = pele::physics::PhysicsType::eos();
      rho = prob_parm.air_state[URHO];
      eos.RTY2P(rho, prob_parm.air_state[UTEMP], massfrac, P);
      // printf("Target pressure: %f [bar]\n",P/1.e+6);

      amrex::Real TARGET_VX = prob_parm.air_state[UMX]/rho;
      amrex::Real TARGET_VY = prob_parm.air_state[UMY]/rho;
      amrex::Real TARGET_VZ = prob_parm.air_state[UMZ]/rho;
      amrex::Real TARGET_TEMPERATURE = prob_parm.air_state[UTEMP];
      amrex::Real TARGET_PRESSURE = P;

      bc_target[0] = TARGET_VX;
      bc_target[1] = TARGET_VY;
      bc_target[2] = TARGET_VZ;
      bc_target[3] = TARGET_TEMPERATURE;
      bc_target[4] = TARGET_PRESSURE;


      amrex::Real relax_U = 0.5;
      amrex::Real relax_V = 0.5;
      amrex::Real relax_W = 0.5;
      amrex::Real relax_T = -0.2;
      amrex::Real beta = -0.6;
      amrex::Real sigma = 0.3;

      // bc_params[0] = relax_U;
      // bc_params[1] = relax_V;
      // bc_params[2] = relax_W;
      // bc_params[3] = relax_T;
      // bc_params[4] = beta;
      // bc_params[5] = sigma;
    }

  }
}

void pc_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
