#include "PeleC.H"
#include "IndexDefines.H"
#include "Forcing.H"
#include "pelec_params.H"

void
PeleC::construct_old_ext_source(amrex::Real time, amrex::Real dt)
{
  const amrex::MultiFab& S_old = get_old_data(State_Type);

  int ng = 0; // None filled

  old_sources[ext_src]->setVal(0.0);

  if (!add_ext_src) {
    return;
  }

  fill_ext_source(time, dt, S_old, S_old, *old_sources[ext_src], ng);

  old_sources[ext_src]->FillBoundary(geom.periodicity());
}

void
PeleC::construct_new_ext_source(amrex::Real time, amrex::Real dt)
{
  const amrex::MultiFab& S_old = get_old_data(State_Type);
  const amrex::MultiFab& S_new = get_new_data(State_Type);

  int ng = 0;

  new_sources[ext_src]->setVal(0.0);

  if (!add_ext_src) {
    return;
  }

  fill_ext_source(time, dt, S_old, S_new, *new_sources[ext_src], ng);
}

void
PeleC::fill_ext_source(
  amrex::Real time,
  amrex::Real dt,
  const amrex::MultiFab& state_old
  /*unused*/,
  const amrex::MultiFab& state_new,
  amrex::MultiFab& ext_src,
  int ng)
{
  auto const& fact =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(state_old.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();

//  auto const& Sns = state_new.const_arrays();
//  auto const& Farrs = ext_src.arrays();
//  auto const& flagarrs = flags.const_arrays();
//  const amrex::IntVect ngs(ng);
//  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> ext_force = {0.0};
//  for (int i = 0; i < static_cast<int>(external_forcing.size()); i++) {
//    ext_force[i] = external_forcing[i];
//  }
//  amrex::ParallelFor(
//    ext_src, ngs, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
//      if (!flagarrs[nbx](i, j, k).isCovered()) {
//        amrex::Real e_force = 0.0;
//        for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
//          Farrs[nbx](i, j, k, UMX + idir) = ext_force[idir];
//          e_force += Sns[nbx](i, j, k, UMX + idir) / Sns[nbx](i, j, k, URHO) *
//                     ext_force[idir];
//        }
//        Farrs[nbx](i, j, k, UEDEN) = e_force;
//      }
//    });
//  amrex::Gpu::synchronize();

#ifdef PELEC_USE_FORCING
    const ProbParmDevice* lprobparm = d_prob_parm_device;
    change_cfl_during_ignition(time,dt,*lprobparm);
#endif

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(ext_src, amrex::TilingIfNotGPU()); mfi.isValid();
       ++mfi) {
    const amrex::Box& bx = mfi.growntilebox(ng);

    const auto& flag_fab = flags[mfi];
    amrex::FabType typ = flag_fab.getType(bx);
    if (typ == amrex::FabType::covered) {
      continue;
    }

    // auto const& So = state_old.array(mfi);
    auto const& Sn = state_new.array(mfi);
    auto const& Farr = ext_src.array(mfi);

#ifdef PELEC_USE_FORCING
    //compute term B in here and only call ignition_source() when B > 0.0
    const auto geomdata = geom.data();

#endif
    
    // Evaluate the external source
    amrex::ParallelFor(
      bx, NVAR, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
        Farr(i, j, k, n) = 0.0;

#ifdef PELEC_USE_FORCING
      amrex::Real E_tot_source;
      amrex::GpuArray<amrex::Real, NUM_SPECIES> omega = {{0.0}};
      amrex::Real E_vib_source;
      LTP_ignition_source(time,
                      dt,
                      i,
                      j,
                      k,
                      geomdata,
            		      *lprobparm,
                      Sn,
                      E_tot_source,
                      E_vib_source,
                      omega);

      // Farr(i, j, k, UEINT) = E_tot_source;
      Farr(i, j, k, UEDEN) = E_tot_source;
      Farr(i, j, k, UFA)   = E_vib_source;
  
      for (int n = 0; n < NUM_SPECIES; n++)
        Farr(i, j, k, UFS+n) = omega[n];

      // Farr(i, j, k, UEDEN) = Sn(i, j, k, URHO)*(E_tot_source + 0.5 * (vx*vx + vy*vy + vz*vz));
#endif
      });
  }
}


#ifdef PELEC_USE_FORCING

void change_cfl_during_ignition(
  amrex::Real time,
  amrex::Real dt,
  ProbParmDevice const& prob_parm){ 

  amrex::Real ign_duration = prob_parm.pulse_duration;         //duration of the ignition pulse [s]
  amrex::Real tau_p   = prob_parm.dwell;             //inter-pulse width [s]
  amrex::Real tau_hat = ign_duration/tau_p;   //pulse width
  amrex::Real N       = prob_parm.npulses;    //Number of pulses
  // amrex::Real A       = 1.5e+14*1.e-06*1.e+7;  //base power [J/m3/s] -> [erg/cm3/s]
  amrex::Real A       = prob_parm.base_power*1.e-06*1.e+7;  //base power [J/m3/s] -> [erg/cm3/s]
  // A = A * dx[0]*dx[1]*dx[2]; // [erg/s]

  amrex::Real t = 0.0;
  if(time >= prob_parm.ltp_start_time){
    t = time - prob_parm.ltp_start_time; //time used to account for local convection of the ignition source
  }

  //parameter B controls ignition over time
  amrex::Real B = (1.+ floor(t/tau_p - 1) - floor(t/tau_p - tau_hat))*heavyside(t-floor(t/N/tau_p));
  
  cfl = 0.3;
  if(B > 0.0){
    //Hack to reduce timestep automatically during ignition
    cfl = 0.01;  
  }
}


amrex::Real heavyside(amrex::Real t){
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
#endif
