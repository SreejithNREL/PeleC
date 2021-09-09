#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  bool binfmt = false;
  bool restart = false;
  
  amrex::Real T_in  = 298.0;
  amrex::Real vn_in = 0.2;
  amrex::Real pamb  = 1.0e+6; //1 bar

  amrex::GpuArray<amrex::Real, NVAR> air_state = {{0.0}};

  amrex::Real bc_params_init[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  amrex::Real bc_target_init[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

};

struct ProbParmHost
{
  std::string iname;
  amrex::Vector<amrex::Real> h_xinput;
  amrex::Vector<amrex::Real> h_uinput;
  amrex::Vector<amrex::Real> h_vinput;
  amrex::Vector<amrex::Real> h_winput;
  amrex::Vector<amrex::Real> h_xarray;
  amrex::Vector<amrex::Real> h_xdiff;
  amrex::Gpu::DeviceVector<amrex::Real> xinput;
  amrex::Gpu::DeviceVector<amrex::Real> uinput;
  amrex::Gpu::DeviceVector<amrex::Real> vinput;
  amrex::Gpu::DeviceVector<amrex::Real> winput;
  amrex::Gpu::DeviceVector<amrex::Real> xarray;
  amrex::Gpu::DeviceVector<amrex::Real> xdiff;
  ProbParmHost()
    : xinput(0), uinput(0), vinput(0), winput(0), xarray(0), xdiff()
  {
  }
};

#endif
