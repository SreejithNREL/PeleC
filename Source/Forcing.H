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
    amrex::Real x = 0.0;
    return x;    
  }
  else{
    amrex::Real x = 1.0;
    return x;
  }
}


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real exp_correlation(
  amrex::Real c,
  amrex::Real pk,
  amrex::Real ak,
  amrex::Real T,
  amrex::Real bk){
  //experimental correlation for the relaxation time
  return c/pk*exp(ak*(pow(T,-1./3.)-bk)-18.42);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real evib_equil(amrex::Real Temp){
    amrex::Real theta_1 = 3396; //[K]  first quantum vibrational state of N2
    
    amrex::Real mw = 28.0; // N2 molecular weight [kg/kmol]
    amrex::Real r = pele::physics::Constants::RU/mw; //[J/kg/K]

    return (r*theta_1)/(exp(theta_1/Temp)-1);

}


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real tau_vt(amrex::Real massfrac[NUM_SPECIES],
                   amrex::Real Temp,
                   amrex::Real rho){

    auto eos = pele::physics::PhysicsType::eos();

    amrex::Real c = 1013250.0; // [dyne/cm2] as pressure in PeleC.   Castela et al. 2016 use c=1 [atm]
    amrex::Real molefrac[NUM_SPECIES] = {0.0};
    amrex::Real P;

    eos.RTY2P(rho, Temp, massfrac, P);
    // P = P * 0.1; //Convert Pele pressure to Pa

    eos.Y2X(massfrac,molefrac);
    amrex::Real mw[NUM_SPECIES];
    eos.molecular_weight(mw);
    
    amrex::Real pk,ak,bk;

    //relaxation time of atomic oxygen 
    amrex::Real tau_O = 1.0e+10; // This will set the term '1/tau_O' close to zero
    //Massfraction of atomic oxygen is zero at the beginning. Avoid division by zero
    if(molefrac[O_ID] > 0.0){
      pk = P*molefrac[O_ID]; //partial pressure obtained assuming ideal gas law
      ak = 72.4;
      bk = 0.015;
      tau_O = exp_correlation(c,pk,ak,Temp,bk);
    }

    //relaxation time of O2
    pk = P*molefrac[O2_ID]; //partial pressure obtained assuming ideal gas law
    ak = 229.0;
    bk = 0.0295;
    amrex::Real tau_O2 = exp_correlation(c,pk,ak,Temp,bk);

    //relaxation time of N2
    pk = P*molefrac[N2_ID]; //partial pressure obtained assuming ideal gas law
    ak = 221.0;
    bk = 0.029;
    amrex::Real tau_N2 = exp_correlation(c,pk,ak,Temp,bk);

    return pow((1/tau_O + 1/tau_O2 + 1/tau_N2),-1);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real R_vt(amrex::Real rho, 
                 amrex::Real e_vib,
                 amrex::Real massfrac[NUM_SPECIES],
                 amrex::Real Temp){

 return rho*(e_vib - evib_equil(Temp))/tau_vt(massfrac,Temp,rho);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real vyaas_energy_dist(
  amrex::Real B,
  amrex::Real x,
  amrex::Real y,
  amrex::Real z){

  const amrex::Real kernel_offset   = 1.92e-2-0.15e-2; //[cm] 0.15 comes from electrode gap divided by 2
  y = y - kernel_offset; //setting y to be zero at the cathode tip [m]

  amrex::Real za      = 3.0e-3; //anode tip [m]
  amrex::Real zc      = 0.0;    //cathode tip
  amrex::Real zi      = zc + (za-zc)/4.0; // streamer collision location

  amrex::Real alpha   = 1.05;
  amrex::Real beta    = 0.98;
  amrex::Real delta_g = 1.e-4; //streamer length scale [mm]
  amrex::Real delta_f = 3.*delta_g; //width of the streamer wake
  amrex::Real Ha      = 0.95;
  amrex::Real Hc      = 0.0 ;

  if(y > zc and y < za){
    // amrex::Print() << "In LTP_ignition_source: B = " << B << "  t = " << t << std::endl;
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

    return B*(f+g*(ha+hc));
  }
  else{
    return 0.0;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real castela_energy_dist(amrex::Real x,
                                amrex::Real y,
                                amrex::Real z){
  // amrex::Real r = 225.e-6; //ignition kernel radius in meters
  // amrex::Real a = 1.63242383e-04;
  // amrex::Real b = 10.;
  amrex::Real a = 690e-10;
  amrex::Real b = 2.5;

  amrex::Real r = sqrt(pow((0.0-x),2)+pow((0.0-y),2)+pow((0.0-z),2)); //radial distance from the ignition center

  return 1.0-pow(erf(pow(r,2)/a),b);
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
  amrex::Real& E_tot_source,
  amrex::Real& E_vib_source,
  amrex::GpuArray<amrex::Real, NUM_SPECIES>& omega){ 

  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  
  amrex::Real x = 1.e-2*(prob_lo[0] + (i + 0.5) * dx[0]); // [m]
  amrex::Real y = 1.e-2*(prob_lo[1] + (j + 0.5) * dx[1] ); // [m]
#if AMREX_SPACEDIM == 3
  amrex::Real z = 1.e-2*(prob_lo[2] + (k + 0.5) * dx[2]); // [m]
#else
  amrex::Real z = 0.0;
#endif

  // ------ Local variables -----
  amrex::Real ign_duration = prob_parm.pulse_duration;         //duration of the ignition pulse [s]
  amrex::Real tau_p   = prob_parm.dwell;             //inter-pulse width [s]
  amrex::Real tau_hat = ign_duration/tau_p;   //pulse width
  amrex::Real N       = prob_parm.npulses;    //Number of pulses

  // amrex::Real Ep      = prob_parm.base_power*1.e+1;  //base power [J/m3] -> [erg/cm3]
  amrex::Real Ep      = prob_parm.base_power*1.e+1/ign_duration;  //base power [J/m3] -> [erg/cm3/s]

  //variables related to Castela's model
  amrex::Real g_chem = 0.35; //fraction of energy going to ultrafast dissociation
  amrex::Real g_heat = 0.2; //fraction of energy going to ultrafast heating
  amrex::Real alpha_castela = g_chem+g_heat; // fraction of energy going to ultrafast dissociation and heating

  amrex::Real massfrac[NUM_SPECIES]; 
  amrex::Real rho = state(i, j, k, URHO);
  amrex::Real e_vib_local = state(i, j, k, UFA);
  amrex::Real Temp        = state(i, j, k, UTEMP);

  //local massfraction
  for (int n = 0; n < NUM_SPECIES; n++)
    massfrac[n] = state(i, j, k, UFS+n)/rho;

  //relaxation rate of the vibrational energy
  amrex::Real Rvt;
  Rvt = R_vt(rho,e_vib_local,massfrac,Temp);
  
  //reaction rates vector is initialized with zero. It is changed to a computed value during ignition (B > 0)
  for (int n = 0; n < NUM_SPECIES; n++)
    omega[n] = 0.0;
  
  // adding source term for the relaxation of the vibrational state
  E_tot_source =  Rvt;
  E_vib_source = -Rvt;


  // ---- Control of energy deposition in time -----
  amrex::Real t = 0.0;
  if(time >= prob_parm.ltp_start_time){
    t = time - prob_parm.ltp_start_time; //time used to compute parameter B
  }
  
  //parameter B controls ignition over time
  amrex::Real B = (1.+ floor(t/tau_p - 1) - floor(t/tau_p - tau_hat))*heavyside_device(t-floor(t/N/tau_p));
  // -----------------------------------------------
  
  if(B > 0.0){

    //accounting for local cross-flow velocity
    amrex::Real Vel[3] = {0.0};
    Vel[0] = state(i, j, k, UMX)/rho;
    Vel[1] = state(i, j, k, UMY)/rho;
    Vel[2] = state(i, j, k, UMZ)/rho;
    // x = x - t*Vel[0];
    // y = y - t*Vel[1];
    // z = z - t*Vel[2];

    // amrex::Real energy_dist = vyaas_energy_dist(B,x,y,z);
    amrex::Real energy_dist = castela_energy_dist(x,y,z);
    
    // Finished computing spatial distribution. 
    // Lets now move to the division of deposited energy into 
    // ultrafast heating and increase in vibrational energy and species dissociation
      
    amrex::Real Y_O2_U; //Oxygen massfraction in the oxidizer
    amrex::Real e[NUM_SPECIES]; //internal energy
    amrex::Real mw[NUM_SPECIES]; //molecular weights
    
    auto eos = pele::physics::PhysicsType::eos();
    eos.T2Ei(Temp, e);
    eos.molecular_weight(mw);

    //Oxygen massfraction in the oxidizer
    Y_O2_U = prob_parm.fuel_state[UFS + O2_ID]/prob_parm.fuel_state[URHO];
    
    //Reaction source terms for ultrafast dissociation. Units [g/cm3/s]
    omega[O_ID]  = g_chem * massfrac[O2_ID]/Y_O2_U * Ep/e[O_ID];
    omega[O2_ID] = - mw[O_ID]/mw[O2_ID] * omega[O_ID];

    //discharge energy consumed in the dissociation process
    amrex::Real E_chem = Ep * g_chem * massfrac[O2_ID]/Y_O2_U*(1 - e[O2_ID]/e[O_ID]);

    //discharge energy used to ultrafast gas heating
    amrex::Real E_heat = alpha_castela*Ep - E_chem;

    //discharge energy used to ultrafast increase in vibrational energy
    amrex::Real E_vib =  Ep*(1-g_heat-g_chem);

    //Taking into account energy distribution from Vyaas' model
    for (int n = 0; n < NUM_SPECIES; n++)
      omega[n] *= energy_dist;

    E_tot_source += (E_heat+E_chem)*energy_dist;
    E_vib_source += (E_vib)        *energy_dist;

  }
}


#endif
#endif
