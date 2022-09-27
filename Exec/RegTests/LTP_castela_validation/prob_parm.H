#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>
#include <turbinflow.H>

struct ProbParmDevice
{
  amrex::GpuArray<amrex::Real, NVAR> H2_state  = {{0.0}};
  amrex::GpuArray<amrex::Real, NVAR> fuel_state  = {{0.0}};
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> L = {{1.0}};

  amrex::Real forcing_u0 = 0.0;
  amrex::Real forcing_v0 = 0.0;
  amrex::Real forcing_w0 = 0.0;
  amrex::Real forcing_force = 0.0;

  //LTP model parameters
  amrex::Real base_power     = 0.0;
  amrex::Real dwell          = 0.0;
  amrex::Real npulses        = 0.0;
  amrex::Real pulse_duration = 0.0;
  amrex::Real ltp_start_time = 0.0;
  amrex::Real time_h2_injection = 0.0;

  amrex::Real pamb    = 1013250.0 * 1.0;
  amrex::Real phi_in  = 1.0;
  amrex::Real moles_NH3 = 0.0;
  amrex::Real moles_H2 = 0.0;
  amrex::Real moles_N2 = 0.0;
  amrex::Real pure_h2_width = 10.0;
  amrex::Real y_center = 1.92;
  amrex::Real T_in    = 298.0;
  amrex::Real vn_in   = 0.0;
  int nt = 0; //number of points to be read at a time in the turbulent file 
  bool restart = false;
  bool init_kernel = false;
  bool first_time = true;
  int x_index[2] = {0};


  amrex::Real* d_xarray       = nullptr;
  amrex::Real* d_xdiff        = nullptr;
  amrex::Real* d_yarray       = nullptr;
  amrex::Real* d_ydiff        = nullptr;
  amrex::Real* d_data_ic      = nullptr;
  amrex::Real* d_fuel_state   = nullptr;
  amrex::Real* d_kernel_state = nullptr;

 // TurbParm tp;
  int  meanFlowDir = 1;

  int nx    = 0;
  int ny    = 0;
  int nz    = 0;
  int nscal = 0;

  amrex::GpuArray<bool,2*AMREX_SPACEDIM> turb_ok = {{false}};

  // == variables to initialize turbulence ==
  amrex::Real* d_xarray_turb  = nullptr;
  amrex::Real* d_xdiff_turb   = nullptr;
  amrex::Real* d_yarray_turb  = nullptr;
  amrex::Real* d_ydiff_turb   = nullptr;
  amrex::Real* d_zarray_turb  = nullptr;
  amrex::Real* d_zdiff_turb   = nullptr;
  amrex::Real* d_data_ic_turb = nullptr;

  int ires_x = 0;
  int ires_y = 0;
  int ires_z = 0;

  bool turbulence = false;
  
  amrex::Real Lxturb = 0.0;
  amrex::Real Lyturb = 0.0;
  amrex::Real Lzturb = 0.0;

  //==================================
};

struct ProbParmHost
{
   bool do_turb = false;
// == variables to initialize turbulence ==
  amrex::Vector<amrex::Real> T_xinput;
  amrex::Vector<amrex::Real> T_xarray;
  amrex::Vector<amrex::Real> T_xdiff;

  amrex::Vector<amrex::Real> T_yinput;
  amrex::Vector<amrex::Real> T_yarray;
  amrex::Vector<amrex::Real> T_ydiff;

  amrex::Vector<amrex::Real> T_zinput;
  amrex::Vector<amrex::Real> T_zarray;
  amrex::Vector<amrex::Real> T_zdiff;

  amrex::Vector<amrex::Real> T_uinput;
  amrex::Vector<amrex::Real> T_vinput;
  amrex::Vector<amrex::Real> T_winput;

  amrex::Vector<amrex::Real> data_turb;


  amrex::Gpu::DeviceVector<amrex::Real> xarray_turb;
  amrex::Gpu::DeviceVector<amrex::Real> xdiff_turb;
  amrex::Gpu::DeviceVector<amrex::Real> yarray_turb;
  amrex::Gpu::DeviceVector<amrex::Real> ydiff_turb;
  amrex::Gpu::DeviceVector<amrex::Real> zarray_turb;
  amrex::Gpu::DeviceVector<amrex::Real> zdiff_turb;
  amrex::Gpu::DeviceVector<amrex::Real> data_ic_turb;

  int nx_tmp    = 0;
  int ny_tmp    = 0;
  int nz_tmp    = 0;
  int nscal_tmp = 0;
  std::string turb_ic;

  // ======= kernel related stuff ======
  amrex::Vector<amrex::Real> v_ydiff;
  amrex::Vector<amrex::Real> v_xdiff;
  amrex::Vector<amrex::Real> v_xinput;
  amrex::Vector<amrex::Real> v_yinput;
  amrex::Vector<amrex::Real> data;
  amrex::Vector<amrex::Real> data_ic;

  std::string iname = "";
  
  amrex::Gpu::DeviceVector<amrex::Real> xarray;
  amrex::Gpu::DeviceVector<amrex::Real> xdiff;
  amrex::Gpu::DeviceVector<amrex::Real> yarray;
  amrex::Gpu::DeviceVector<amrex::Real> ydiff;

  ProbParmHost()
    : xarray(0), 
      yarray(0), 
      xarray_turb(0),
      yarray_turb(0),
      zarray_turb(0),
      data_ic_turb(0),
      xdiff(), 
      ydiff(), 
      xdiff_turb(),
      ydiff_turb(),
      zdiff_turb()
  {
  }

};

#endif
