#ifndef _FORCING_H_
#define _FORCING_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>
#include "IndexDefines.H"

#ifdef PELEC_USE_FORCING
void change_cfl_during_ignition(
  amrex::Real time,
  amrex::Real dt,
  ProbParmDevice const& prob_parm);


amrex::Real heavyside(amrex::Real t);

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real heavyside_device(amrex::Real t){
  if(t < 0.0){
    amrex::Real x = 0.0;
    return x;
  }
  else if(t == 0.0){
    amrex::Real x = 0.5;
    return x;    
  }
  else{
    amrex::Real x = 1.0;
    return x;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void LTP_ignition_source(
  amrex::Real time,
  amrex::Real dt,
  int i,
  int j,
  int k,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm,
  amrex::Array4<const amrex::Real> const& state,
  amrex::Real& S_ltp){ 

  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real kernel_offset   = 1.92-0.15; //[cm] 0.15 comes from electrode gap divided by 2
  amrex::Real kernel_location[3];
  
  amrex::Real za      = 3.0e-3; //anode tip [m]
  amrex::Real zc      = 0.0;    //cathode tip
  amrex::Real zi      = zc + (za-zc)/4.0; // streamer collision location

  amrex::Real x = 1.e-2*(prob_lo[0] + (i + 0.5) * dx[0]); // [m]
  amrex::Real y = 1.e-2*(prob_lo[1] + (j + 0.5) * dx[1] - kernel_offset); //setting y to be zero at the cathode tip [m]
  amrex::Real z = 1.e-2*(prob_lo[2] + (k + 0.5) * dx[2]); // [m]


  amrex::Real alpha   = 1.05;
  amrex::Real beta    = 0.98;
  amrex::Real delta_g = 1.e-4; //streamer length scale [mm]
  amrex::Real delta_f = 3.*delta_g; //width of the streamer wake
  amrex::Real Ha      = 0.95;
  amrex::Real Hc      = 0.0 ;

  amrex::Real ign_duration = prob_parm.pulse_duration;         //duration of the ignition pulse [s]
  amrex::Real tau_p   = prob_parm.dwell;             //inter-pulse width [s]
  amrex::Real tau_hat = ign_duration/tau_p;   //pulse width
  amrex::Real N       = prob_parm.npulses;    //Number of pulses
  // amrex::Real A       = 1.5e+14*1.e-06*1.e+7;  //base power [J/m3/s] -> [erg/cm3/s]
  amrex::Real A       = prob_parm.base_power*1.e-06*1.e+7;  //base power [J/m3/s] -> [erg/cm3/s]
  // A = A * dx[0]*dx[1]*dx[2]; // [erg/s]

  amrex::Real t = time - 1.0e-3; //time used to account for local convection of the ignition source
  
  //parameter B controls ignition over time
  amrex::Real B = (1.+ floor(t/tau_p - 1) - floor(t/tau_p - tau_hat))*heavyside_device(t-floor(t/N/tau_p));
  

  // if(time > ign_duration and B > 0.0) t = t - tau_p; //set t=0 at the start of the second pulse to compute convection distance
  S_ltp = 0.0;
  if(B > 0.0 and x > -5e-02 and x < 5e-02){

    amrex::Real rho = state(i, j, k, URHO);

    //accounting for local cross-flow velocity
    amrex::Real Vel[3] = {0.0};
    Vel[0] = state(i, j, k, UMX)/rho;
    Vel[1] = state(i, j, k, UMY)/rho;
    Vel[2] = state(i, j, k, UMZ)/rho;
    // x = x - t*Vel[0];
    // y = y - t*Vel[1];
    // z = z - t*Vel[2];

    if(y > zc and y < za){
      amrex::Real ha;
      if(y >= zi and y <= za){
        amrex::Real da = zi*zi-2.0*zi*za+za*za;
        amrex::Real aa = 4.0*(1.0-Ha)/da;
        amrex::Real ba = 4.0*(Ha*(zi+za)-zi-za)/da;
        amrex::Real ca = (2.0*zi*za*(1.0-2.0*Ha)+zi*zi+za*za)/da;

        ha = aa*y*y + ba*y + ca;
      }
      else{
        ha = 0.0;
      }

      amrex::Real hc;
      if(y >= zc and y <= zi){
        // Ha = maximum value of ha between za and zi
        amrex::Real dc = zi*zi-2.0*zi*zc+zc*zc;
        amrex::Real ac = 4.0*(1.0-Hc)/dc;
        amrex::Real bc = 4.0*(Hc*(zi+zc)-zi-zc)/dc;
        amrex::Real cc = (2.0*zi*zc*(1.0-2.0*Hc)+zi*zi+zc*zc)/dc;

        hc = ac*y*y + bc*y + cc;
      }
      else{
        hc = 0.0;
      }

      amrex::Real g = exp(-(x*x + z*z)/delta_g/delta_g/2.);

      amrex::Real mz;
      if(y >= alpha*zi){
        mz = za;
      }
      else if(y <= beta*zi){
        mz = zc;
      }
      else{
        mz = zi;
      }

      amrex::Real f = exp(-((x*x + z*z)/delta_g/delta_g/2. + pow((y-mz),2.)/delta_f/delta_f/2.));

      S_ltp = A*B*(f+g*(ha+hc));
    }
  }
}

#endif
#endif
