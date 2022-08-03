#include "prob.H"

void
init_bc()
{
  amrex::Real vt;
  amrex::Real ek;
  amrex::Real T;
  amrex::Real rho;
  amrex::Real e;
  amrex::Real massfrac[NUM_SPECIES];

  for (int n = 0; n < NUM_SPECIES; n++)
    massfrac[n] = 1.0;

  T = PeleC::h_prob_parm_device->T_in;

  const amrex::Real p = PeleC::h_prob_parm_device->pamb;

  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(p, massfrac, T, rho, e);

  vt = 0.0; // PeleC::h_prob_parm_device->vn_in;
  ek = 0.5 * (vt * vt);

  PeleC::h_prob_parm_device->air_state[URHO] = rho;
  PeleC::h_prob_parm_device->air_state[UMX] = 0.0;
  PeleC::h_prob_parm_device->air_state[UMY] = 0.0;
  PeleC::h_prob_parm_device->air_state[UMZ] = 0.0;
  PeleC::h_prob_parm_device->air_state[UEINT] = rho * e;
  PeleC::h_prob_parm_device->air_state[UEDEN] = rho * (e + ek);
  PeleC::h_prob_parm_device->air_state[UTEMP] = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    PeleC::h_prob_parm_device->air_state[UFS + n] = rho * massfrac[n];
  }
}

void
pc_prob_close()
{
}

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex_real* problo,
  const amrex_real* probhi)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("pamb", PeleC::h_prob_parm_device->pamb);
  pp.query("T_in", PeleC::h_prob_parm_device->T_in);
  pp.query("vn_in", PeleC::h_prob_parm_device->vn_in);
  pp.query("direction", PeleC::h_prob_parm_device->direction);

  init_bc();
}
}

void
gaussian_function(
  amrex::Real x, amrex::Real mean, amrex::Real sigma, amrex::Real gaussian)
{
  amrex::Real pi = 3.1415;
  gaussian = 1 / sigma / (2. * pi) * exp(-0.5 * pow((x - mean) / sigma, 2));
}

void
PeleC::problem_post_timestep()
{
}

void
PeleC::problem_post_init()
{
}

void
PeleC::problem_post_restart()
{
}
