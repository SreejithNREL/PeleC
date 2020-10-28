#include <AMReX_DistributionMapping.H>

#include "PeleC.H"
#include "React.H"
#include "EOS.H"
#include "base_getrates.h"
#ifdef USE_SUNDIALS_PP
#include <reactor.h>
#endif

void
PeleC::react_state(
  amrex::Real time, amrex::Real dt, bool react_init, amrex::MultiFab* A_aux)
{
  /*
    Update I_R, and recompute S_new
   */
  BL_PROFILE("PeleC::react_state()");

  const amrex::Real strt_time = amrex::ParallelDescriptor::second();

  AMREX_ASSERT(do_react == 1);

  if (verbose && amrex::ParallelDescriptor::IOProcessor()) {
    if (react_init) {
      amrex::Print() << "... Initializing reactions, using interval dt = " << dt
                     << std::endl;
    } else {
      amrex::Print() << "... Computing reactions for dt = " << dt << std::endl;
    }
  }

  amrex::MultiFab& S_new = get_new_data(State_Type);
  const int ng = S_new.nGrow();
  prefetchToDevice(S_new);

  // Create a MultiFab with all of the non-reacting source terms.

  amrex::MultiFab Atmp, *Ap;

  if (A_aux == nullptr || react_init) {
    Atmp.define(grids, dmap, NVAR, ng, amrex::MFInfo(), Factory());
    Atmp.setVal(0);
    Ap = &Atmp;
  }

  if (!react_init) {
    // Build non-reacting source term, and an S_new that does not include
    // reactions

    if (A_aux == nullptr) {
      for (int n = 0; n < src_list.size(); ++n) {
        amrex::MultiFab::Saxpy(
          Atmp, 0.5, *new_sources[src_list[n]], 0, 0, NVAR, ng);
        amrex::MultiFab::Saxpy(
          Atmp, 0.5, *old_sources[src_list[n]], 0, 0, NVAR, ng);
      }
      if (do_hydro && !do_mol) {
        amrex::MultiFab::Add(Atmp, hydro_source, 0, 0, NVAR, ng);
      }
    } else {
      Ap = A_aux;
    }

    amrex::MultiFab& S_old = get_old_data(State_Type);
    amrex::MultiFab::Copy(S_new, S_old, 0, 0, NVAR, ng);
    amrex::MultiFab::Saxpy(S_new, dt, *Ap, 0, 0, NVAR, ng);
  }

  amrex::MultiFab& reactions = get_new_data(Reactions_Type);
  reactions.setVal(0.0);
  prefetchToDevice(reactions);
  if (use_reactions_work_estimate) {
    amrex::Abort("Need to implement redistribution of chemistry work");
  }

#ifdef PELEC_USE_EB
  auto const& fact =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(S_new.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();
#endif

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  {
    for (amrex::MFIter mfi(S_new, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {

      const amrex::Box& bx = mfi.growntilebox(ng);

      auto const& uold =
        react_init ? S_new.array(mfi) : get_old_data(State_Type).array(mfi);
      auto const& unew = S_new.array(mfi);
      auto const& a = Ap->array(mfi);
      amrex::FArrayBox w(bx, 1);
      amrex::Elixir w_eli = w.elixir();
      auto const& w_arr = w.array();
      auto const& I_R = reactions.array(mfi);
      const int do_update = react_init ? 0 : 1; // TODO: Update here? Or just get reaction source?

#ifdef PELEC_USE_EB
      const auto& flag_fab = flags[mfi];
      amrex::FabType typ = flag_fab.getType(bx);
      if (typ == amrex::FabType::covered) {
        continue;
      } else if (
        typ == amrex::FabType::singlevalued || typ == amrex::FabType::regular)
#endif
      {
        if (chem_integrator == 1) {
          const int nsubsteps_min = adaptrk_nsubsteps_min;
          const int nsubsteps_max = adaptrk_nsubsteps_max;
          const int nsubsteps_guess = adaptrk_nsubsteps_guess;
          const amrex::Real errtol = adaptrk_errtol;
          const auto len = amrex::length(bx); 
          const auto lo = amrex::lbound(bx);
          const auto uo = amrex::ubound(bx);
          const int ncells = len.x * len.y * len.z; //bx.numPts()
          
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
          
          amrex::Real dt_rk = dt / nsubsteps_guess;
          const amrex::Real dt_min = dt / nsubsteps_min;
          const amrex::Real dt_max = dt / nsubsteps_max;
          amrex::Real updt_time = 0.0;
          
          amrex::Real* rhoedot_ext;
          amrex::Real* rhoe_carryover;
          amrex::Real* rhoe_rk;
          
          cudaError_t cuda_status = cudaSuccess;
          cudaMallocManaged(&rhoedot_ext, ncells * sizeof(amrex::Real));
          cudaMallocManaged(&rhoe_carryover, ncells * sizeof(amrex::Real));
          cudaMallocManaged(&rhoe_rk, ncells * sizeof(amrex::Real));
          
          
          amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            amrex::Real rhou = uold(i, j, k, UMX);
            amrex::Real rhov = uold(i, j, k, UMY);
            amrex::Real rhow = uold(i, j, k, UMZ);
            amrex::Real rho_old = uold(i, j, k, URHO);
            amrex::Real rhoInv = 1.0 / rho_old;
            amrex::Real rho = 0.;

            for (int nsp = UFS; nsp < (UFS + NUM_SPECIES); nsp++) {
              rho += uold(i, j, k, nsp);
            }

            amrex::Real nrg = (uold(i, j, k, UEDEN) - (0.5 * (rhou * rhou + rhov * rhov + rhow * rhow) * rhoInv)) * rhoInv;

            rhou = unew(i, j, k, UMX);
            rhov = unew(i, j, k, UMY);
            rhow = unew(i, j, k, UMZ);
            rhoInv = 1.0 / unew(i, j, k, URHO);
            
            int offset = (k - lo.z) * len.x * len.y + (j - lo.y) * len.x + (i - lo.x);
            
            rhoedot_ext[offset] = ((unew(i, j, k, UEDEN) -
                                  (0.5 * (rhou * rhou + rhov * rhov + rhow * rhow) * rhoInv)) - rho * nrg) /  dt;
            
            rhoe_carryover[offset] = rho * nrg;
            rhoe_rk[offset] = rho * nrg;
            
          });
          
          cuda_status = cudaStreamSynchronize(amrex::Gpu::gpuStream());
          
          amrex::Array4<amrex::Real> urk = unew;          
          amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            for (int n = 0; n < NVAR; n++) {
              urk(i,j,k,n) = uold(i,j,k,n);
            }
          });
          
          cuda_status = cudaStreamSynchronize(amrex::Gpu::gpuStream());
          
	  amrex::Array4<amrex::Real> urk_carryover = unew;
          amrex::Array4<amrex::Real> urk_err = unew;
          
          amrex::Real* massfrac;
          amrex::Real* pressure;
          amrex::Real* temperature;
          amrex::Real* mixMW;
          amrex::Real* diffusion;
          amrex::Real* wdot;
          
          cudaMallocManaged(&massfrac, NUM_SPECIES  * ncells * sizeof(amrex::Real));
          cudaMallocManaged(&wdot, NUM_SPECIES  * ncells * sizeof(amrex::Real));
          cudaMallocManaged(&diffusion, NUM_SPECIES  * ncells * sizeof(amrex::Real));
          cudaMallocManaged(&pressure, ncells * sizeof(amrex::Real));
          cudaMallocManaged(&temperature, ncells * sizeof(amrex::Real));
          cudaMallocManaged(&mixMW, ncells * sizeof(amrex::Real));
            
          // Do the RK!
          int steps = 0;
          while (updt_time < dt) {
            
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              for (int n = 0; n < NVAR; n++) {
                urk_carryover(i,j,k,n) = urk(i,j,k,n);
                urk_err(i,j,k,n) = 0.0;
              }
              
              int offset = (k - lo.z) * len.x * len.y + (j - lo.y) * len.x + (i - lo.x);
              rhoe_carryover[offset] = rhoe_rk[offset];
            });
                      
            for (int stage = 0; stage < 6; stage++) {
              
              amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real Y[NUM_SPECIES] = {};
                int offset = (k - lo.z) * len.x * len.y + (j - lo.y) * len.x + (i - lo.x);
                for (int n = 0; n < NUM_SPECIES; ++n) {
                  massfrac[offset * NUM_SPECIES + n] = urk(i,j,k,UFS + n) / urk(i,j,k,URHO);
                  Y[n] = urk(i,j,k,UFS + n) / urk(i,j,k,URHO);
                } 
                amrex::Real P;
                amrex::Real R = urk(i,j,k,URHO);
                amrex::Real T = urk(i,j,k,UTEMP);
                EOS::RTY2P(R, T, Y, P);
                pressure[offset] = P;
                temperature[offset] = urk(i,j,k,UTEMP);
                
                amrex::Real wbar;
                CKMMWY(&Y[0], &wbar); 
                mixMW[offset] = 1.0/wbar;
                for (int n = 0; n < NUM_SPECIES; n++)  diffusion[offset * NUM_SPECIES + n] = 0.0;
                
              });
              
              cuda_status = cudaStreamSynchronize(amrex::Gpu::gpuStream());
              //base_getrates (pressure, temperature, mixMW, massfrac, diffusion, dt_rk, wdot);
              
              
            /*================ Adapt Time step! ======================== */
            } // end rk stages
            updt_time += dt_rk;
            steps += 1;
            //adapt_timestep(urk_err, dt_max, dt_rk, dt_min, errtol);
          } // end timestep loop

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!// 

          /*amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
              pc_expl_reactions(
                i, j, k, uold, unew, a, w_arr, I_R, dt, nsubsteps_min,
                nsubsteps_max, nsubsteps_guess, errtol, do_update);
            });
          //printf("\nncells = %d (%d, %d, %d)\n", ncells, uo.x, uo.y, uo.z);*/
 
        } else if (chem_integrator == 2) {
#ifdef USE_SUNDIALS_PP
          const auto len = amrex::length(bx);
          const auto lo = amrex::lbound(bx);
          const int ncells = len.x * len.y * len.z;
          int reactor_type = 1;
          amrex::Real fabcost;
          amrex::Real current_time = 0.0;

          amrex::Real* rY_in;
          amrex::Real* rY_src_in;
          amrex::Real* re_in;
          amrex::Real* re_src_in;

#ifdef USE_CUDA_SUNDIALS_PP
          cudaError_t cuda_status = cudaSuccess;
          cudaMallocManaged(
            &rY_in, (NUM_SPECIES + 1) * ncells * sizeof(amrex::Real));
          cudaMallocManaged(
            &rY_src_in, NUM_SPECIES * ncells * sizeof(amrex::Real));
          cudaMallocManaged(&re_in, ncells * sizeof(amrex::Real));
          cudaMallocManaged(&re_src_in, ncells * sizeof(amrex::Real));

          int ode_ncells = ncells;
#else
          rY_in = new amrex::Real[ncells * (NUM_SPECIES + 1)];
          rY_src_in = new amrex::Real[ncells * (NUM_SPECIES)];
          re_in = new amrex::Real[ncells];
          re_src_in = new amrex::Real[ncells];

          int ode_ncells = 1;
#endif
          amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              amrex::Real rhou = uold(i, j, k, UMX);
              amrex::Real rhov = uold(i, j, k, UMY);
              amrex::Real rhow = uold(i, j, k, UMZ);
              amrex::Real rho_old = uold(i, j, k, URHO);
              amrex::Real rhoInv = 1.0 / rho_old;
              amrex::Real rho = 0.;

              for (int nsp = UFS; nsp < (UFS + NUM_SPECIES); nsp++) {
                rho += uold(i, j, k, nsp);
              }

              amrex::Real nrg =
                (uold(i, j, k, UEDEN) -
                 (0.5 * (rhou * rhou + rhov * rhov + rhow * rhow) * rhoInv)) *
                rhoInv;

              rhou = unew(i, j, k, UMX);
              rhov = unew(i, j, k, UMY);
              rhow = unew(i, j, k, UMZ);
              rhoInv = 1.0 / unew(i, j, k, URHO);

              amrex::Real rhoedot_ext =
                ((unew(i, j, k, UEDEN) -
                  (0.5 * (rhou * rhou + rhov * rhov + rhow * rhow) * rhoInv)) -
                 rho * nrg) /
                dt;

              int offset =
                (k - lo.z) * len.x * len.y + (j - lo.y) * len.x + (i - lo.x);
              for (int nsp = 0; nsp < NUM_SPECIES; nsp++) {
                rY_in[offset * (NUM_SPECIES + 1) + nsp] =
                  uold(i, j, k, UFS + nsp);
                rY_src_in[offset * NUM_SPECIES + nsp] = a(i, j, k, UFS + nsp);
              }
              rY_in[offset * (NUM_SPECIES + 1) + NUM_SPECIES] =
                uold(i, j, k, UTEMP);
              re_in[offset] = uold(i, j, k, UEINT);
              re_src_in[offset] = rhoedot_ext;
            });

#ifdef USE_CUDA_SUNDIALS_PP
          cuda_status = cudaStreamSynchronize(amrex::Gpu::gpuStream());
#endif
          fabcost = 0.0;
          for (int i = 0; i < ncells; i += ode_ncells) {

#ifdef USE_CUDA_SUNDIALS_PP
            fabcost += react(
              rY_in + i * (NUM_SPECIES + 1), rY_src_in + i * NUM_SPECIES,
              re_in + i, re_src_in + i, &dt, &current_time, reactor_type,
              ode_ncells, amrex::Gpu::gpuStream());
#else
            fabcost += react(
              rY_in + i * (NUM_SPECIES + 1), rY_src_in + i * NUM_SPECIES,
              re_in + i, re_src_in + i, &dt, &current_time);
#endif
          }
          fabcost = fabcost / ncells;

          // unpack data
          amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              w_arr(i, j, k) = fabcost;
              amrex::Real rhou = uold(i, j, k, UMX);
              amrex::Real rhov = uold(i, j, k, UMY);
              amrex::Real rhow = uold(i, j, k, UMZ);
              amrex::Real rho_old = uold(i, j, k, URHO);
              amrex::Real rhoInv = 1.0 / rho_old;

              amrex::Real rho = 0.;
              for (int nsp = UFS; nsp < (UFS + NUM_SPECIES); nsp++) {
                rho += uold(i, j, k, nsp);
              }
              amrex::Real nrg =
                (uold(i, j, k, UEDEN) -
                 (0.5 * (rhou * rhou + rhov * rhov + rhow * rhow) * rhoInv)) *
                rhoInv;

              rhou = unew(i, j, k, UMX);
              rhov = unew(i, j, k, UMY);
              rhow = unew(i, j, k, UMZ);
              rhoInv = 1.0 / unew(i, j, k, URHO);

              amrex::Real rhoedot_ext =
                ((unew(i, j, k, UEDEN) -
                  (0.5 * (rhou * rhou + rhov * rhov + rhow * rhow) * rhoInv)) -
                 rho * nrg) /
                dt;

              amrex::Real umnew = uold(i, j, k, UMX) + dt * a(i, j, k, UMX);
              amrex::Real vmnew = uold(i, j, k, UMY) + dt * a(i, j, k, UMY);
              amrex::Real wmnew = uold(i, j, k, UMZ) + dt * a(i, j, k, UMZ);
              amrex::Real rhonew;

              int offset =
                (k - lo.z) * len.x * len.y + (j - lo.y) * len.x + (i - lo.x);
              for (int nsp = 0; nsp < NUM_SPECIES; nsp++) {
                rhonew += rY_in[offset * (NUM_SPECIES + 1) + nsp];
              }

              if (do_update) {
                unew(i, j, k, URHO) = rhonew;
                unew(i, j, k, UMX) = umnew;
                unew(i, j, k, UMY) = vmnew;
                unew(i, j, k, UMZ) = wmnew;
                for (int nsp = 0; nsp < NUM_SPECIES; nsp++) {
                  unew(i, j, k, UFS + nsp) =
                    rY_in[offset * (NUM_SPECIES + 1) + nsp];
                }
                unew(i, j, k, UTEMP) =
                  rY_in[offset * (NUM_SPECIES + 1) + NUM_SPECIES];
              }

              for (int nsp = 0; nsp < NUM_SPECIES; nsp++) {
                I_R(i, j, k, nsp) = (rY_in[offset * (NUM_SPECIES + 1) + nsp] -
                                     uold(i, j, k, UFS + nsp)) /
                                      dt -
                                    rY_src_in[offset * (NUM_SPECIES) + nsp];
              }
              I_R(i, j, k, NUM_SPECIES) =
                ((nrg * rho_old) + dt * rhoedot_ext +
                 0.5 * (umnew * umnew + vmnew * vmnew + wmnew * wmnew) /
                   rhonew -
                 uold(i, j, k, UEDEN)) /
                  dt -
                a(i, j, k, UEDEN);
            });

          if (do_react_load_balance || do_mol_load_balance) {
            get_new_data(Work_Estimate_Type)[mfi].plus<amrex::RunOn::Device>(w);
          }
#else
          amrex::Abort(
            "chem_integrator=2 which requires Sundials to be enabled");
#endif
        } else {
          amrex::Abort("chem_integrator must be equal to 1 or 2");
        }
      }
    }
  }

  if (ng > 0)
    S_new.FillBoundary(geom.periodicity());

  if (verbose > 1) {

    const int IOProc = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::Real run_time = amrex::ParallelDescriptor::second() - strt_time;

#ifdef AMREX_LAZY
    Lazy::QueueReduction([=]() mutable {
#endif
      amrex::ParallelDescriptor::ReduceRealMax(run_time, IOProc);

      if (amrex::ParallelDescriptor::IOProcessor())
        amrex::Print() << "PeleC::react_state() time = " << run_time << "\n";
#ifdef AMREX_LAZY
    });
#endif
  }
}
