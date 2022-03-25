#include "prob.H"

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
  {
    amrex::ParmParse pp("prob");
    pp.query("reynolds", PeleC::h_prob_parm_device->reynolds);
    pp.query("mach", PeleC::h_prob_parm_device->mach);
    pp.query("prandtl", PeleC::h_prob_parm_device->prandtl);
    pp.query("rho_x_fact", PeleC::h_prob_parm_device->rho_x_fact);
    pp.query("rho_y_fact", PeleC::h_prob_parm_device->rho_y_fact);
    pp.query("rho_z_fact", PeleC::h_prob_parm_device->rho_z_fact);
    pp.query("u_0_fact", PeleC::h_prob_parm_device->u_0_fact);
    pp.query("v_0_fact", PeleC::h_prob_parm_device->v_0_fact);
    pp.query("w_0_fact", PeleC::h_prob_parm_device->w_0_fact);
    pp.query("u_r_fact", PeleC::h_prob_parm_device->u_r_fact);
    pp.query("v_r_fact", PeleC::h_prob_parm_device->v_r_fact);
    pp.query("w_r_fact", PeleC::h_prob_parm_device->w_r_fact);
    pp.query("p_r_fact", PeleC::h_prob_parm_device->p_r_fact);
    pp.query("a_rhox", PeleC::h_prob_parm_device->a_rhox);
    pp.query("a_rhoy", PeleC::h_prob_parm_device->a_rhoy);
    pp.query("a_rhoz", PeleC::h_prob_parm_device->a_rhoz);
    pp.query("a_ux", PeleC::h_prob_parm_device->a_ur);
    pp.query("a_ux", PeleC::h_prob_parm_device->a_vr);
    pp.query("a_ux", PeleC::h_prob_parm_device->a_wr);
    pp.query("a_ux", PeleC::h_prob_parm_device->a_pr);
  }

  // Define the length scale
  PeleC::h_prob_parm_device->L_x = probhi[0] - problo[0];
  PeleC::h_prob_parm_device->L_y = probhi[1] - problo[1];
  PeleC::h_prob_parm_device->L_z = probhi[2] - problo[2];

  // Initial density, velocity, and material properties
  amrex::Real eint;
  amrex::Real cs;
  amrex::Real cp;
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(
    PeleC::h_prob_parm_device->p0, massfrac, PeleC::h_prob_parm_device->T0,
    PeleC::h_prob_parm_device->rho0, eint);
  eos.RTY2Cs(
    PeleC::h_prob_parm_device->rho0, PeleC::h_prob_parm_device->T0, massfrac,
    cs);
  eos.TY2Cp(PeleC::h_prob_parm_device->T0, massfrac, cp);

  PeleC::h_prob_parm_device->u0 = PeleC::h_prob_parm_device->mach * cs;

  auto& trans_parm = PeleC::trans_parms.host_trans_parm();
  trans_parm.const_bulk_viscosity =
    PeleC::h_prob_parm_device->rho0 * PeleC::h_prob_parm_device->u0 *
    PeleC::h_prob_parm_device->L_x / PeleC::h_prob_parm_device->reynolds;
  trans_parm.const_diffusivity = 0.0;
  trans_parm.const_viscosity =
    PeleC::h_prob_parm_device->rho0 * PeleC::h_prob_parm_device->u0 *
    PeleC::h_prob_parm_device->L_x / PeleC::h_prob_parm_device->reynolds;
  trans_parm.const_conductivity =
    trans_parm.const_viscosity * cp / PeleC::h_prob_parm_device->prandtl;
  PeleC::trans_parms.sync_to_device();

  // clang-format off
  // MASA parameters for the following functions
  // rho = rho_0 + rho_x * cos(a_rhox * PI * x / L) * rho_y * cos(a_rhoy * PI * y / L) * rho_z * cos(a_rhoz * PI * z / L);
  // u = u_0 + u_r * cos(a_ux * PI * r / L)
  // v = v_0 + v_r * cos(a_vx * PI * r / L)
  // w = w_0 + w_r * cos(a_wx * PI * r / L)
  // p = p_0 + p_r * cos(a_px * PI * r / L)
  // clang-format on
  masa_set_param("L", PeleC::h_prob_parm_device->L_x);
  masa_set_param(
    "R", pele::physics::Constants::RU / pele::physics::Constants::AIRMW);
  masa_set_param("k", trans_parm.const_conductivity);
  masa_set_param("Gamma", eos.gamma);
  masa_set_param("mu", trans_parm.const_viscosity);
  masa_set_param("mu_bulk", trans_parm.const_bulk_viscosity);
  masa_set_param("rho_0", PeleC::h_prob_parm_device->rho0);
  masa_set_param(
    "rho_x",
    PeleC::h_prob_parm_device->rho_x_fact * PeleC::h_prob_parm_device->rho0);
  masa_set_param("rho_y", PeleC::h_prob_parm_device->rho_y_fact);
  masa_set_param("rho_z", PeleC::h_prob_parm_device->rho_z_fact);
  masa_set_param(
    "u_0", PeleC::h_prob_parm_device->u_0_fact * PeleC::h_prob_parm_device->u0);
  masa_set_param(
    "v_0", PeleC::h_prob_parm_device->v_0_fact * PeleC::h_prob_parm_device->u0);
  masa_set_param(
    "w_0", PeleC::h_prob_parm_device->w_0_fact * PeleC::h_prob_parm_device->u0);
  masa_set_param("p_0", PeleC::h_prob_parm_device->p0);
  masa_set_param(
    "u_r", PeleC::h_prob_parm_device->u_r_fact * PeleC::h_prob_parm_device->u0);
  masa_set_param(
    "v_r", PeleC::h_prob_parm_device->v_r_fact * PeleC::h_prob_parm_device->u0);
  masa_set_param(
    "w_r", PeleC::h_prob_parm_device->w_r_fact * PeleC::h_prob_parm_device->u0);
  masa_set_param(
    "p_r", PeleC::h_prob_parm_device->p_r_fact * PeleC::h_prob_parm_device->p0);
  masa_set_param("a_rhox", PeleC::h_prob_parm_device->a_rhox);
  masa_set_param("a_rhoy", PeleC::h_prob_parm_device->a_rhoy);
  masa_set_param("a_rhoz", PeleC::h_prob_parm_device->a_rhoz);
  masa_set_param("a_ur", PeleC::h_prob_parm_device->a_ur);
  masa_set_param("a_vr", PeleC::h_prob_parm_device->a_vr);
  masa_set_param("a_wr", PeleC::h_prob_parm_device->a_wr);
  masa_set_param("a_pr", PeleC::h_prob_parm_device->a_pr);

  // Display and check
  if (amrex::ParallelDescriptor::IOProcessor()) {
    masa_display_param();
  }
  masa_sanity_check();
}
}

void
PeleC::problem_post_timestep()
{
  if ((verbose <= 0) || (!do_mms)) {
    return;
  }

  int finest_level = parent->finestLevel();
  amrex::Real time = state[State_Type].curTime();
  amrex::Real rho_mms_err = 0.0;
  amrex::Real u_mms_err = 0.0;
  amrex::Real v_mms_err = 0.0;
  amrex::Real w_mms_err = 0.0;
  amrex::Real p_mms_err = 0.0;
  amrex::Real rho_residual = 0.0;
  amrex::Real rhou_residual = 0.0;
  amrex::Real rhov_residual = 0.0;
  amrex::Real rhow_residual = 0.0;
  amrex::Real rhoE_residual = 0.0;

#ifdef PELEC_USE_MASA
  if (level == 0) {
    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "... MMS problem post timestep" << std::endl;
    }

    // Calculate the errors and residuals
    for (int lev = 0; lev <= finest_level; lev++) {
      PeleC& pc_lev = getLevel(lev);

      const bool local_flag = true;
      rho_mms_err += pc_lev.volWgtSquaredSum("rhommserror", time, local_flag);
      u_mms_err += pc_lev.volWgtSquaredSum("ummserror", time, local_flag);
      v_mms_err += pc_lev.volWgtSquaredSum("vmmserror", time, local_flag);
      w_mms_err += pc_lev.volWgtSquaredSum("wmmserror", time, local_flag);
      p_mms_err += pc_lev.volWgtSquaredSum("pmmserror", time, local_flag);

      rho_residual += pc_lev.volWgtSquaredSumDiff(Density, time, local_flag);
      rhou_residual += pc_lev.volWgtSquaredSumDiff(Xmom, time, local_flag);
      rhov_residual += pc_lev.volWgtSquaredSumDiff(Ymom, time, local_flag);
      rhow_residual += pc_lev.volWgtSquaredSumDiff(Zmom, time, local_flag);
      rhoE_residual += pc_lev.volWgtSquaredSumDiff(Eden, time, local_flag);
    }

    // Reductions
    amrex::ParallelDescriptor::ReduceRealSum(
      &rho_mms_err, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
      &u_mms_err, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
      &v_mms_err, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
      &w_mms_err, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
      &p_mms_err, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
      &rho_residual, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
      &rhou_residual, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
      &rhov_residual, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
      &rhow_residual, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
      &rhoE_residual, 1, amrex::ParallelDescriptor::IOProcessorNumber());

    // Get the norm and normalize it
    amrex::MultiFab vol(grids, dmap, 1, 0);
    amrex::MultiFab::Copy(vol, volume, 0, 0, 1, 0);
    if (eb_in_domain) {
      amrex::MultiFab::Multiply(vol, vfrac, 0, 0, 1, 0);
    }
    amrex::Real V = vol.sum(0, false);
    rho_mms_err = std::sqrt(rho_mms_err / V);
    u_mms_err = std::sqrt(u_mms_err / V);
    v_mms_err = std::sqrt(v_mms_err / V);
    w_mms_err = std::sqrt(w_mms_err / V);
    p_mms_err = std::sqrt(p_mms_err / V);
    rho_residual = std::sqrt(rho_residual / V);
    rhou_residual = std::sqrt(rhou_residual / V);
    rhov_residual = std::sqrt(rhov_residual / V);
    rhow_residual = std::sqrt(rhow_residual / V);
    rhoE_residual = std::sqrt(rhoE_residual / V);

    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "TIME= " << time << " RHO MMS ERROR  = " << rho_mms_err
                     << '\n';
      amrex::Print() << "TIME= " << time << " U MMS ERROR    = " << u_mms_err
                     << '\n';
      amrex::Print() << "TIME= " << time << " V MMS ERROR    = " << v_mms_err
                     << '\n';
      amrex::Print() << "TIME= " << time << " W MMS ERROR    = " << w_mms_err
                     << '\n';
      amrex::Print() << "TIME= " << time << " P MMS ERROR    = " << p_mms_err
                     << '\n';
      amrex::Print() << "TIME= " << time << " RHO RESIDUAL   = " << rho_residual
                     << '\n';
      amrex::Print() << "TIME= " << time
                     << " RHO*U RESIDUAL = " << rhou_residual << '\n';
      amrex::Print() << "TIME= " << time
                     << " RHO*V RESIDUAL = " << rhov_residual << '\n';
      amrex::Print() << "TIME= " << time
                     << " RHO*W RESIDUAL = " << rhow_residual << '\n';
      amrex::Print() << "TIME= " << time
                     << " RHO*E RESIDUAL = " << rhoE_residual << '\n';

      const int log_index = find_datalog_index("mmslog");
      if (log_index >= 0) {

        std::ostream& data_log2 = parent->DataLog(log_index);

        // Write the quantities at this time
        const int datwidth = 14;
        const int datprecision = 6;
        data_log2 << std::setw(datwidth) << time;
        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << rho_mms_err;
        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << u_mms_err;
        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << v_mms_err;
        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << w_mms_err;
        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << p_mms_err;
        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << rho_residual;
        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << rhou_residual;
        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << rhov_residual;
        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << rhow_residual;
        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << rhoE_residual;
        data_log2 << std::endl;
      }
    }
  }

#else
  Error("MASA is not turned on. Turn on with USE_MASA=TRUE.");
#endif
}

void
PeleC::problem_post_init()
{
  if ((verbose <= 0) || (!do_mms)) {
    return;
  }

  amrex::Real time = state[State_Type].curTime();

  if (level == 0) {
    if (amrex::ParallelDescriptor::IOProcessor()) {
      const int log_index = find_datalog_index("mmslog");
      if (log_index >= 0) {
        std::ostream& data_log2 = parent->DataLog(log_index);
        if (time == 0.0) {
          const int datwidth = 14;
          data_log2 << std::setw(datwidth) << "          time";
          data_log2 << std::setw(datwidth) << "   rho_mms_err";
          data_log2 << std::setw(datwidth) << "     u_mms_err";
          data_log2 << std::setw(datwidth) << "     v_mms_err";
          data_log2 << std::setw(datwidth) << "     w_mms_err";
          data_log2 << std::setw(datwidth) << "     p_mms_err";
          data_log2 << std::setw(datwidth) << "  rho_residual";
          data_log2 << std::setw(datwidth) << " rhou_residual";
          data_log2 << std::setw(datwidth) << " rhov_residual";
          data_log2 << std::setw(datwidth) << " rhow_residual";
          data_log2 << std::setw(datwidth) << " rhoE_residual";
          data_log2 << std::endl;
        }
      }
    }
  }
}

void
PeleC::problem_post_restart()
{
}
