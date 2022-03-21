#ifndef _FORCING_H_
#define _FORCING_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>
#include "IndexDefines.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
LTP_ignition_source(
      amrex::Real time,
      amrex::Real dt,
      int i,
      int j,
      int k,
      amrex::GeometryData const& geomdata,
      ProbParmDevice const& prob_parm,
      amrex::Array4<const amrex::Real> const& state,
      amrex::Real& S_ltp);

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real heavyside(amrex::Real t);

#endif
