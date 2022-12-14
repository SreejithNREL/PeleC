#include "prob.H"

void
init_bc()
{
  amrex::Real vt, ek, a, yl, yr, sumY, T, rho, e;
  amrex::Real molefrac[NUM_SPECIES], massfrac[NUM_SPECIES];

  init_composition_air(molefrac);

  double p = PeleC::h_prob_parm_device->pamb;
  auto eos = pele::physics::PhysicsType::eos();
  
  amrex::Real mw[NUM_SPECIES];
  eos.molecular_weight(mw);
  amrex::Real mw_N2 = mw[N2_ID];

  amrex::Real Ru = pele::physics::Constants::RU/mw_N2; // [erg/g/K]
  amrex::Real theta_v_s = 3371.; // [K]
  amrex::Real theta_dof = 5.;
  amrex::Real Evib = Ru*theta_v_s / (std::exp(theta_v_s/PeleC::h_prob_parm_device->T_vib)-1);
  amrex::Real Etr  = theta_dof*Ru*PeleC::h_prob_parm_device->T_tr/2.0; // [erg/g]

  eos.X2Y(molefrac, massfrac);
  eos.EY2T(Etr, massfrac, T); //get temperature
  eos.PYT2RE(p, massfrac, T, rho, Etr); //get rho

  vt = 0.0;
  ek = 0.5 * (vt * vt);
  //Evib = 0.0; 

  // (PeleC::h_prob_parm_device->H2_state)[URHO]  = rho;
  // (PeleC::h_prob_parm_device->H2_state)[UMX]   = rho * vt;
  // (PeleC::h_prob_parm_device->H2_state)[UMY]   = 0.0;
  // (PeleC::h_prob_parm_device->H2_state)[UMZ]   = 0.0;
  // (PeleC::h_prob_parm_device->H2_state)[UEINT] = rho * e;
  // (PeleC::h_prob_parm_device->H2_state)[UEDEN] = rho * (e + ek);
  // (PeleC::h_prob_parm_device->H2_state)[UTEMP] = T;
  // for (int n = 0; n < NUM_SPECIES; n++){
  //   (PeleC::h_prob_parm_device->H2_state)[UFS + n] = rho * massfrac[n];
  // }
  

  (PeleC::h_prob_parm_device->fuel_state)[URHO]  = rho;
  (PeleC::h_prob_parm_device->fuel_state)[UMX]   = 0.0;
  (PeleC::h_prob_parm_device->fuel_state)[UMY]   = 0.0;
  (PeleC::h_prob_parm_device->fuel_state)[UMZ]   = 0.0;
  (PeleC::h_prob_parm_device->fuel_state)[UEINT] = rho * Etr;
  (PeleC::h_prob_parm_device->fuel_state)[UEDEN] = rho * (Etr + ek);
  (PeleC::h_prob_parm_device->fuel_state)[UTEMP] = T;
  (PeleC::h_prob_parm_device->fuel_state)[UFA]   = rho*Evib; // 1000 [K] vibrational temperature
  for (int n = 0; n < NUM_SPECIES; n++){
    (PeleC::h_prob_parm_device->fuel_state)[UFS + n] = rho * massfrac[n];
  }
}

void
init_composition_air(amrex::Real molefrac[NUM_SPECIES])
{
  for (int n = 0; n < NUM_SPECIES; n++)
    molefrac[n] = 0.0;

  amrex::Real a;
  a = 0.5;

  molefrac[O2_ID] = 0.0;
  molefrac[N2_ID] = 1.0 -  molefrac[O2_ID];

  //molefrac[O_ID] = 0.061543*2.0;
  //molefrac[O2_ID] = 0.21 - molefrac[O_ID];
  //molefrac[N2_ID] = 1.0 -  molefrac[O2_ID] - molefrac[O_ID];

  amrex::Print() << "== pure air composition ==" << std::endl;
  amrex::Print() << "X_O2 = " << molefrac[O2_ID] << std::endl;
  amrex::Print() << "X_N2 = " << molefrac[N2_ID] << std::endl;

}


void
pc_prob_close()
{

}

extern "C" 
{
  void amrex_probinit(
              const int* init,
              const int* name,
              const int* namelen,
              const amrex_real* problo,
              const amrex_real* probhi)
    {

      amrex::ParmParse pp("prob");
      pp.query("pamb"           , PeleC::h_prob_parm_device->pamb);
      pp.query("phi_in"         , PeleC::h_prob_parm_device->phi_in);
      pp.query("base_power"     , PeleC::h_prob_parm_device->base_power);
      pp.query("npulses"        , PeleC::h_prob_parm_device->npulses);
      pp.query("pulse_duration" , PeleC::h_prob_parm_device->pulse_duration);
      pp.query("ltp_start_time" , PeleC::h_prob_parm_device->ltp_start_time);
      pp.query("time_h2_injection", PeleC::h_prob_parm_device->time_h2_injection);
      pp.query("dwell"          , PeleC::h_prob_parm_device->dwell);
      pp.query("T_tr"           , PeleC::h_prob_parm_device->T_tr);
      pp.query("T_vib"          , PeleC::h_prob_parm_device->T_vib);
      pp.query("vn_in"          , PeleC::h_prob_parm_device->vn_in);
      pp.query("turbulence"     , PeleC::h_prob_parm_device->turbulence);
      pp.query("restart"        , PeleC::h_prob_parm_device->restart);

      pp.query("moles_NH3"      , PeleC::h_prob_parm_device->moles_NH3);
      pp.query("moles_H2"       , PeleC::h_prob_parm_device->moles_H2);
      pp.query("moles_N2"       , PeleC::h_prob_parm_device->moles_N2);
      pp.query("pure_h2_width"  , PeleC::h_prob_parm_device->pure_h2_width);


      amrex::ParmParse pp2("geometry");

      amrex::Vector<amrex::Real> prob_lo(AMREX_SPACEDIM);
      amrex::Vector<amrex::Real> prob_hi(AMREX_SPACEDIM);

      pp2.queryarr("prob_lo",prob_lo,0,AMREX_SPACEDIM);
      pp2.queryarr("prob_hi",prob_hi,0,AMREX_SPACEDIM);

      PeleC::h_prob_parm_device->y_center = 1.92;

      init_bc();
    }
}




void PeleC::problem_post_timestep()
{
}

// void PeleC::problem_post_init(){}
void PeleC::problem_post_init(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata)
{
//
//  amrex::Real rho = h_prob_parm_device->fuel_state[URHO];
//  amrex::Real velx = state(i, j, k, UMX); //state() at this point is storing velocity, not momentum
//  amrex::Real vely = state(i, j, k, UMY); //state() at this point is storing velocity, not momentum
//  amrex::Real velz = state(i, j, k, UMZ); //state() at this point is storing velocity, not momentum
//
//  amrex::Real sum = 0.0;
//  for (int n = 0; n < NUM_SPECIES; n++){
//    state(i, j, k, UFS + n) = h_prob_parm_device->fuel_state[UFS + n];
//  }
//
//  amrex::Real e = h_prob_parm_device->fuel_state[UEINT];
//
//  state(i, j, k, URHO)  = rho;
//  state(i, j, k, UMX)   = rho * velx;
//  state(i, j, k, UMY)   = rho * vely;
//  state(i, j, k, UMZ)   = rho * velz;
//  state(i, j, k, UTEMP) = h_prob_parm_device->fuel_state[UTEMP];
//  state(i, j, k, UEINT) = e;
//  state(i, j, k, UEDEN) = e + rho * (0.5 * (velx * velx + vely * vely + velz * velz)); 
//
//
}

void PeleC::problem_post_restart(){}

