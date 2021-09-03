//***** THIS IS NOT FINISHED FROM BEING PORTED FROM FORTRAN YET *******

#include "PeleC.H"
#include "NSCBC.H"
#include "prob.H"

//------------------------------
// Imposing Ghost-Cells Navier-Stokes Characteristic BCs.
// For the theory, see Motheau et al. AIAA J. Vol. 55, No. 10 : pp. 3399-3408,
// 2017.
//
// Note that for the corner treatment, we depart from the AIAA paper, because
// we found out that the corner coupling method was superfluous and that
// providing transverse terms computed from one-sided derivative do the job.
//
//------------------------------

void
PeleC::impose_NSCBC(
  const amrex::Box& bx,
  const amrex::Array4<amrex::Real>& uin,
  const amrex::Array4<amrex::Real>& q,
  const amrex::Array4<amrex::Real>& qaux,
  const amrex::Box& qbox,
  AMREX_D_DECL(
    const amrex::Array4<int>& x_bcMask,
    const amrex::Array4<int>& y_bcMask,
    const amrex::Array4<int>& z_bcMask),
  const int nscbc_isAnyPerio,
  const amrex::Real time,
  const amrex::Real dt)
{
  const amrex::Box dom = geom.Domain();
  const int* domlo = dom.loVect();
  const int* domhi = dom.hiVect();
  const int domlox = domlo[0];
  const int domhix = domhi[0];
  const int* q_lo = qbox.loVect();
  const int* q_hi = qbox.hiVect();
  const int domloy = domlo[1];
  const int domhiy = domhi[1];
  const int domloz = domlo[2];
  const int domhiz = domhi[2];
  const amrex::Real* prob_lo = geom.ProbLo();
  const amrex::Real* prob_hi = geom.ProbHi();
  const amrex::Real* dx = geom.CellSize();
  // BC params are relaxation factors for the NSCBC
  // TODO: Hard-coded for now, should be variable
  const amrex::Real relax_U = 0.5;
  const amrex::Real relax_V = 0.5;
  const amrex::Real relax_W = 0.5;
  const amrex::Real relax_T = -0.2;
  const amrex::Real beta = 1.;
  const amrex::Real sigma = -0.6;
  amrex::GpuArray<amrex::Real, 6> bc_params = {relax_U, relax_V, relax_W, relax_T, beta, sigma};
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> problen =
    {prob_hi[0] - prob_lo[0], prob_hi[1] - prob_lo[1], prob_hi[2] - prob_lo[2]};
  const auto& bcs = PeleC::phys_bc;
  const ProbParmDevice* lprobparm = d_prob_parm_device;

  // print*,"===== Inputs to impose_NSCBC() ===="
  // print*,"x_bcMask",x_bcMask(x_bcMask_l1:x_bcMask_h1,x_bcMask_l2:x_bcMask_h2,x_bcMask_l3:x_bcMask_h3)
  // print*,"y_bcMask",y_bcMask(y_bcMask_l1:y_bcMask_h1,y_bcMask_l2:y_bcMask_h2,y_bcMask_l3:y_bcMask_h3)
  // print*,"z_bcMask",z_bcMask(z_bcMask_l1:z_bcMask_h1,z_bcMask_l2:z_bcMask_h2,z_bcMask_l3:z_bcMask_h3)

  // print*,"q",q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
  // print*,"qaux",qaux(qa_l1:qa_h1,qa_l2:qa_h2,qa_l3:qa_h3,NQAUX)
  // print*,"uin",uin(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,NVAR)
  // print*,""

  // Note the BC indices are
  // 0: Interior (periodic)
  // 1: Hard
  // 2: FOExtrap
  // 3: Symmetry
  // 4: SlipWall
  // 5: NoSlipWal
  // 6: UserBC
  // TODO: bcnormal should allow for returning an Inflow and Outflow types for NSCBCs
  // For now we assume 7: Inflow and 8: Outflow
  if (nscbc_isAnyPerio == 0) {
    //--------------------------------------------------------------------------
    // corners
    //--------------------------------------------------------------------------
    if (((q_hi[0] > domhi[0]) || (q_lo[0] < domlo[0])) &&
        ((q_hi[1] > domhi[1]) || (q_lo[1] < domlo[1])) &&
        ((q_hi[2] > domhi[2]) || (q_lo[2] < domlo[2]))) {
      amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if (((i == domlo[0]) || (i == domhi[0])) &&
            ((j == domlo[1]) || (j == domhi[1])) &&
            ((k == domlo[2]) || (k == domhi[2]))) {
          int test_keyword_x, test_keyword_y, test_keyword_z;
          int x_isgn, y_isgn, z_isgn;
          int x_idx_Mask, y_idx_Mask, z_idx_Mask;
          if (i == domhi[0]) {
            x_isgn = -1;
            test_keyword_x = bcs.hi(0);
            x_idx_Mask = i + 1;
          } else {
            x_isgn = 1;
            test_keyword_x = bcs.lo(0);
            x_idx_Mask = i;
          }
          if (j == domhi[1]) {
            y_isgn = -1;
            test_keyword_y = bcs.hi(1);
            y_idx_Mask = j + 1;
          } else {
            y_isgn = 1;
            test_keyword_y = bcs.lo(1);
            y_idx_Mask = j;
          }
          // Normal derivative along x
          amrex::Real dpdx, dudx, dvdx, dwdx, drhodx;
          normal_derivative(i, j, k, 0, x_isgn, dx[0], dpdx, dudx, dvdx, dwdx, drhodx, q);
          // Normal derivative along y
          amrex::Real dpdy, dudy, dvdy, dwdy, drhody;
          normal_derivative(i, j, k, 1, y_isgn, dx[1], dpdy, dudy, dvdy, dwdy, drhody, q);
          // Normal derivative along x
          amrex::Real dpdz, dudz, dvdz, dwdz, drhodz;
          normal_derivative(i, j, k, 2, z_isgn, dx[2], dpdz, dudz, dvdz, dwdz, drhodz, q);

          amrex::GpuArray<amrex::Real, 5> Tx = {{0.0}};
          amrex::GpuArray<amrex::Real, 5> Ty = {{0.0}};
          amrex::GpuArray<amrex::Real, 5> Tz = {{0.0}};
          // Compute transverse terms for X
          compute_transverse_terms(
            i, j, k, 0, Tx.data(), dpdx, dudx, dvdx, dwdx, drhodx, dpdy, dudy, dvdy, dwdy, drhody,
            dpdz, dudz, dvdz, dwdz, drhodz, q, qaux);
          // Compute transverse terms for Y
          compute_transverse_terms(
            i, j, k, 1, Ty.data(), dpdx, dudx, dvdx, dwdx, drhodx, dpdy, dudy, dvdy, dwdy, drhody,
            dpdz, dudz, dvdz, dwdz, drhodz, q, qaux);
          // Compute transverse terms for Z
          compute_transverse_terms(
            i, j, k, 2, Tz.data(), dpdx, dudx, dvdx, dwdx, drhodx, dpdy, dudy, dvdy, dwdy, drhody,
            dpdz, dudz, dvdz, dwdz, drhodz, q, qaux);
          amrex::GpuArray<amrex::Real, NVAR> x_bc_target;
          amrex::GpuArray<amrex::Real, NVAR> y_bc_target;
          amrex::GpuArray<amrex::Real, NVAR> z_bc_target;
          amrex::GpuArray<amrex::Real, NVAR> s_int;

          // LODI system waves for X
          amrex::GpuArray<amrex::Real, 5> Lx = {{0.0}};
          // LODI system waves for Y
          amrex::GpuArray<amrex::Real, 5> Ly = {{0.0}};
          // LODI system waves for Z
          amrex::GpuArray<amrex::Real, 5> Lz = {{0.0}};
          const amrex::Real x[AMREX_SPACEDIM] = {AMREX_D_DECL(
            prob_lo[0] + static_cast<amrex::Real>(i + 0.5) * dx[0],
            prob_lo[1] + static_cast<amrex::Real>(j + 0.5) * dx[1],
            prob_lo[2] + static_cast<amrex::Real>(k + 0.5) * dx[2])};
          int x_bc_type = test_keyword_x;
          int y_bc_type = test_keyword_y;
          int z_bc_type = test_keyword_z;
          // UserBC type is 6
          // TODO: I believe this won't work until bcnormal is designed to allow
          //       BC types to be returned based on current index, commenting out for now
          // Now just assumes low is inflow and high is outflow
          if (test_keyword_x == 6) {
            bcnormal(x, s_int.data(), x_bc_target.data(), 0, x_isgn, time, geom.data(), *lprobparm);
            // TODO: Hard-coded Inflow and Outflow, bcnormal should provide these
            if (x_isgn == 1)
              x_bc_type = 7; // Inflow
            else
              x_bc_type = 8; //Outflow
          }
          x_bcMask(i, j, k) = x_bc_type;
          if (test_keyword_y == 6) {
            bcnormal(x, s_int.data(), y_bc_target.data(), 1, y_isgn, time, geom.data(), *lprobparm);
            // TODO: Hard-coded Inflow and Outflow, bcnormal should provide these
            if (y_isgn == 1)
              y_bc_type = 7; // Inflow
            else
              y_bc_type = 8; //Outflow
          }
          y_bcMask(i, j, k) = y_bc_type;
          if (test_keyword_x == 6) {
            bcnormal(x, s_int.data(), z_bc_target.data(), 2, y_isgn, time, geom.data(), *lprobparm);
            // TODO: Hard-coded Inflow and Outflow, bcnormal should provide these
            if (z_isgn == 1)
              z_bc_type = 7; // Inflow
            else
              z_bc_type = 8; //Outflow
          }
          z_bcMask(i, j, k) = z_bc_type;
          compute_waves(i, j, k, 0, x_isgn, x_bc_type, problen.data(), bc_params.data(),
            x_bc_target.data(), Tx.data(), Lx.data(), dpdx, dudx, dvdx, dwdx, drhodx, q, qaux);
          compute_waves(i, j, k, 1, y_isgn, y_bc_type, problen.data(), bc_params.data(),
            y_bc_target.data(), Ty.data(), Ly.data(), dpdy, dudy, dvdy, dwdy, drhody, q, qaux);
          compute_waves(i, j, k, 2, z_isgn, z_bc_type, problen.data(), bc_params.data(),
            z_bc_target.data(), Tz.data(), Lz.data(), dpdz, dudz, dvdz, dwdz, drhodz, q, qaux);
          update_ghost_cells(i, j, k, x_bc_type, 0, x_isgn, dx[0], domlo, domhi, Lx.data(), uin, q, qaux);
          update_ghost_cells(i, j, k, y_bc_type, 1, y_isgn, dx[1], domlo, domhi, Ly.data(), uin, q, qaux);
          update_ghost_cells(i, j, k, z_bc_type, 2, z_isgn, dx[2], domlo, domhi, Lz.data(), uin, q, qaux);
        }
       });
    }
  }
// }
//     if (
//       ((q_hi[0] > domhi[0]) || (q_lo[0] < domlo[0])) &&
//       ((q_hi[1] > domhi[1]) || (q_lo[1] < domlo[1])) &&
//       ((q_hi[2] > domhi[2]) || (q_lo[2] < domlo[2]))) {

//       if (q_hi[0] > domhi[1]) {
//         test_keyword_x = PeleC::phys_bc.hi()[0];
//         i = domhi[0];
//         x_isign = -1;
//         x_idx_Mask = i + 1;
//       } else if (q_lo[0] < domlo[0]) {
//         test_keyword_x = PeleC::phys_bc.lo()[0];
//         i = domlo[0];
//         x_isign = 1;
//         x_idx_Mask = i;
//       }

//       if (q_hi[1] > domhi[1]) {
//         test_keyword_y = PeleC::phys_bc.hi()[1];
//         j = domhi[1];
//         y_isign = -1;
//         y_idx_Mask = j + 1;
//       } else if (q_lo[1] < domlo[1]) {
//         test_keyword_y = PeleC::phys_bc.lo()[1];
//         j = domlo[1];
//         y_isign = 1;
//         y_idx_Mask = j;
//       }

//       if (q_hi[2] > domhi[2]) {
//         test_keyword_z = PeleC::phys_bc.hi()[2];
//         k = domhi[2];
//         z_isign = -1;
//         z_idx_Mask = k + 1;
//       } else if (q_lo[2] < domlo[2]) {
//         test_keyword_z = PeleC::phys_bc.lo()[2];
//         k = domlo[2];
//         z_isign = 1;
//         z_idx_Mask = k;
//       }

//       x = (static_cast<amrex::Real>(i) + 0.5) * dx;
//       y = (static_cast<amrex::Real>(j) + 0.5) * dx[1];
//       z = (static_cast<amrex::Real>(k) + 0.5) * dx[2];

//       // Normal derivative along x
//       normal_derivative(
//         i, j, k, 1, x_isign, dx, dpdx, dudx, dvdx, dwdx, drhodx, q, q_l1, q_l2,
//         q_l3, q_h1, q_h2, q_h3);

//       // Normal derivative along y
//       normal_derivative(
//         i, j, k, 2, y_isign, dx[1], dpdy, dudy, dvdy, dwdy, drhody, q, q_l1, q_l2,
//         q_l3, q_h1, q_h2, q_h3);

//       // Normal derivative along z
//       normal_derivative(
//         i, j, k, 3, z_isign, dx[2], dpdz, dudz, dvdz, dwdz, drhodz, q, q_l1, q_l2,
//         q_l3, q_h1, q_h2, q_h3);

//       // Compute transverse terms for X
//       compute_transverse_terms(
//         i, j, k, 1, T1_X, T2_X, T3_X, T4_X, T5_X, dpdx, dudx, dvdx, dwdx,
//         drhodx, dpdy, dudy, dvdy, dwdy, drhody, dpdz, dudz, dvdz, dwdz, drhodz,
//         q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1,
//         qa_h2, qa_h3);

//       // Compute transverse terms for X
//       compute_transverse_terms(
//         i, j, k, 2, T1_Y, T2_Y, T3_Y, T4_Y, T5_Y, dpdx, dudx, dvdx, dwdx,
//         drhodx, dpdy, dudy, dvdy, dwdy, drhody, dpdz, dudz, dvdz, dwdz, drhodz,
//         q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1,
//         qa_h2, qa_h3);

//       // Compute transverse terms for X
//       compute_transverse_terms(
//         i, j, k, 3, T1_Z, T2_Z, T3_Z, T4_Z, T5_Z, dpdx, dudx, dvdx, dwdx,
//         drhodx, dpdy, dudy, dvdy, dwdy, drhody, dpdz, dudz, dvdz, dwdz, drhodz,
//         q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1,
//         qa_h2, qa_h3);

//       // Calling user target BC values
//       // x face
//       if (test_keyword_x == UserBC) {
//         // bcnormal([x,y,z],U_dummy,U_ext,1,x_isign,time,x_bc_type,x_bc_params,x_bc_target);
//       } else {
//         x_bc_type = test_keyword_x;
//       }
//       x_bcMask(x_idx_Mask, j, k) = x_bc_type;

//       // y face
//       if (test_keyword_y == UserBC) {
//         // bcnormal([x,y,z],U_dummy,U_ext,2,y_isign,time,y_bc_type,y_bc_params,y_bc_target);
//       } else {
//         y_bc_type = test_keyword_y;
//       }
//       y_bcMask(i, y_idx_Mask, k) = y_bc_type;

//       // z face
//       if (test_keyword_z == UserBC) {
//         // bcnormal([x,y,z],U_dummy,U_ext,3,z_isign,time,z_bc_type,z_bc_params,z_bc_target);
//       } else {
//         z_bc_type = test_keyword_z;
//       }
//       z_bcMask(i, j, z_idx_Mask) = z_bc_type;

//       // Computing the LODI system waves along X
//       compute_waves(
//         i, j, k, 1, x_isign, x_bc_type, x_bc_params, x_bc_target, T1_X, T2_X,
//         T3_X, T4_X, T5_X, L1, L2, L3, L4, L5, dpdx, dudx, dvdx, dwdx, drhodx, q,
//         q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1,
//         qa_h2, qa_h3);

//       // Computing the LODI system waves along Y
//       compute_waves(
//         i, j, k, 2, y_isign, y_bc_type, y_bc_params, y_bc_target, T1_Y, T2_Y,
//         T3_Y, T4_Y, T5_Y, M1, M2, M3, M4, M5, dpdy, dudy, dvdy, dwdy, drhody, q,
//         q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1,
//         qa_h2, qa_h3);

//       // Computing the LODI system waves along Z
//       compute_waves(
//         i, j, k, 3, z_isign, z_bc_type, z_bc_params, z_bc_target, T1_Z, T2_Z,
//         T3_Z, T4_Z, T5_Z, N1, N2, N3, N4, N5, dpdz, dudz, dvdz, dwdz, drhodz, q,
//         q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1,
//         qa_h2, qa_h3);

//       // Recomputing ghost-cells values with the LODI waves along X
//       update_ghost_cells(
//         i, j, k, x_bc_type, 1, x_isign, dx, domlo, domhi, L1, L2, L3, L4, L5,
//         uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, q, q_l1, q_l2,
//         q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3);

//       // Recomputing ghost-cells values with the LODI waves along Y
//       update_ghost_cells(
//         i, j, k, y_bc_type, 2, y_isign, dx[1], domlo, domhi, M1, M2, M3, M4, M5,
//         uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, q, q_l1, q_l2,
//         q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3);

//       // Recomputing ghost-cells values with the LODI waves along Y
//       update_ghost_cells(
//         i, j, k, z_bc_type, 3, z_isign, dx[2], domlo, domhi, N1, N2, N3, N4, N5,
//         uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, q, q_l1, q_l2,
//         q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3);
//     }
//   }
  //--------------------------------------------------------------------------
  // lower X
  //--------------------------------------------------------------------------

  if ((q_lo[0] < domlo[0]) && (PeleC::phys_bc.lo()[0] == UserBC)) {
    int i = domlo[0];
    amrex::Real  x = (static_cast<amrex::Real>(i) + 0.5) * dx[0];
    for (int j = q_lo[1] + 1; j < q_hi[1] - 1; j++) {
      amrex::Real y = (static_cast<amrex::Real>(j) + 0.5) * dx[1];
      if (nscbc_isAnyPerio == 0) {
        if ((j == domlo[1]) || (j == domhi[1])) {
          continue; // Doing that to avoid ghost cells already filled by corners
        }
      }
      for (int k = q_lo[2] + 1; k < q_hi[2] - 1; k++) {
        amrex::Real z = (static_cast<amrex::Real>(k) + 0.5) * dx[2];
        if (nscbc_isAnyPerio == 0) {
          if ((k == domlo[2]) || (k == domhi[2])) {
            continue; // Doing that to avoid ghost cells already filled by
                      // corners
          }
        }
        const amrex::Real x_array[AMREX_SPACEDIM] = {AMREX_D_DECL(x,y,z)};

        amrex::Real dpdx, dudx, dvdx, dwdx, drhodx;
        // Normal derivative along x
        normal_derivative(i, j, k, 0, 1, dx[0], dpdx, dudx, dvdx, dwdx, drhodx, q);

        amrex::Real dpdy, dudy, dvdy, dwdy, drhody;
        // Tangential derivative along y
        tangential_derivative(i, j, k, 1, dx[1], dpdy, dudy, dvdy, dwdy, drhody, q);

        amrex::Real dpdz, dudz, dvdz, dwdz, drhodz;
        // Tangential derivative along z
        tangential_derivative(i, j, k, 2, dx[2], dpdz, dudz, dvdz, dwdz, drhodz, q);

        amrex::GpuArray<amrex::Real, 5> Tx = {{0.0}};
        // Compute transverse terms
        compute_transverse_terms(
          i, j, k, 0, Tx.data(), dpdx, dudx, dvdx, dwdx,
          drhodx, dpdy, dudy, dvdy, dwdy, drhody, dpdz, dudz, dvdz, dwdz,
          drhodz, q, qaux);

        amrex::GpuArray<amrex::Real, NVAR> x_bc_target;
        amrex::GpuArray<amrex::Real, NVAR> s_int;

        // Calling user target BC values
        bcnormal(x_array, s_int.data(), x_bc_target.data(), 0, 1, time, geom.data(), *lprobparm);
        int x_bc_type = 7; // Hard-coded inflow. This variable should be updated in bcnormal()

        // Filling bcMask with specific user defined BC type
        if (
          (j < q_lo[1] + 3) || (j > q_hi[1] - 3) || (k < q_lo[2] + 3) ||
          (k > q_hi[2] - 3)) {
          // do nothing.
          // There is just 1 ghost-cell with bcMask because of the Riemann
          // solver
        } else {
          x_bcMask(i, j, k) = x_bc_type;
        }

        // LODI system waves for X
        amrex::GpuArray<amrex::Real, 5> Lx = {{0.0}};

        // Computing the LODI system waves
        compute_waves(i, j, k, 0, 1, x_bc_type, problen.data(), bc_params.data(),
            x_bc_target.data(), Tx.data(), Lx.data(), dpdx, dudx, dvdx, dwdx, drhodx, q, qaux);

        // Recomputing ghost-cells values with the LODI waves
        update_ghost_cells(i, j, k, x_bc_type, 0, 1, dx[0], domlo, domhi, Lx.data(), uin, q, qaux);
      }
    }
  }

  //--------------------------------------------------------------------------
  // upper X
  //--------------------------------------------------------------------------

  if ((q_hi[0] > domhi[0]) && (PeleC::phys_bc.hi()[0] == UserBC)) {
    int i = domhi[0];
    amrex::Real  x = (static_cast<amrex::Real>(i) + 0.5) * dx[0];
    for (int j = q_lo[1] + 1; j < q_hi[1] - 1; j++) {
      amrex::Real  y = (static_cast<amrex::Real>(j) + 0.5) * dx[1];
      if (nscbc_isAnyPerio == 0) {
        if ((j == domlo[1]) || (j == domhi[1])) {
          continue; // Doing that to avoid ghost cells already filled by corners
        }
      }
      for (int k = q_lo[2] + 1; k < q_hi[2] - 1; k++) {
        amrex::Real  z = (static_cast<amrex::Real>(k) + 0.5) * dx[2];
        if (nscbc_isAnyPerio == 0) {
          if ((k == domlo[2]) || (k == domhi[2])) {
            continue; // Doing that to avoid ghost cells already filled by
                      // corners
          }
        }
        const amrex::Real x_array[AMREX_SPACEDIM] = {AMREX_D_DECL(x,y,z)};

        amrex::Real dpdx, dudx, dvdx, dwdx, drhodx;
        // Normal derivative along x
        normal_derivative(i, j, k, 0, 1, dx[0], dpdx, dudx, dvdx, dwdx, drhodx, q);
        
        amrex::Real dpdy, dudy, dvdy, dwdy, drhody;
        // Tangential derivative along y
        tangential_derivative(i, j, k, 1, dx[1], dpdy, dudy, dvdy, dwdy, drhody, q);

        amrex::Real dpdz, dudz, dvdz, dwdz, drhodz;
        // Tangential derivative along z
        tangential_derivative(i, j, k, 2, dx[2], dpdz, dudz, dvdz, dwdz, drhodz, q);

        amrex::GpuArray<amrex::Real, 5> Tx = {{0.0}};
        // Compute transverse terms
        compute_transverse_terms(
          i, j, k, 0, Tx.data(), dpdx, dudx, dvdx, dwdx,
          drhodx, dpdy, dudy, dvdy, dwdy, drhody, dpdz, dudz, dvdz, dwdz,
          drhodz, q, qaux);

        amrex::GpuArray<amrex::Real, NVAR> x_bc_target;
        amrex::GpuArray<amrex::Real, NVAR> s_int;

        // Calling user target BC values
        bcnormal(x_array, s_int.data(), x_bc_target.data(), 0, -1, time, geom.data(), *lprobparm);
        int x_bc_type = 8; // Hard-coded outflow. This variable should be updated in bcnormal()

        // Filling bcMask with specific user defined BC type
        // bcnormal([x,y,z],U_dummy,U_ext,1,-1,time,bc_type,bc_params,bc_target);
        if (
          (j < q_lo[1] + 3) || (j > q_hi[1] - 3) || (k < q_lo[2] + 3) ||
          (k > q_hi[2] - 3)) {
          // do nothing
          // There is just 1 ghost-cell with bcMask because of the Riemann
          // solver
        } else {
          x_bcMask(i + 1, j, k) = x_bc_type;
        }

        // LODI system waves for X
        amrex::GpuArray<amrex::Real, 5> Lx = {{0.0}};

        // Computing the LODI system waves
        compute_waves(i, j, k, 0, -1, x_bc_type, problen.data(), bc_params.data(),
            x_bc_target.data(), Tx.data(), Lx.data(), dpdx, dudx, dvdx, dwdx, drhodx, q, qaux);

        // Recomputing ghost-cells values with the LODI waves
        update_ghost_cells(i, j, k, x_bc_type, 0, -1, dx[0], domlo, domhi, Lx.data(), uin, q, qaux);
      }
    }
  }

  //--------------------------------------------------------------------------
  // lower Y
  //--------------------------------------------------------------------------

  if ((q_lo[1] < domlo[1]) && (PeleC::phys_bc.lo()[1] == UserBC)) {
    int j = domlo[1];
    amrex::Real  y = (static_cast<amrex::Real>(j) + 0.5) * dx[1];
    for (int i = q_lo[0] + 1; i < q_hi[0] - 1; i++) {
      amrex::Real  x = (static_cast<amrex::Real>(i) + 0.5) * dx[0];
      if (nscbc_isAnyPerio == 0) {
        if ((i == domlo[0]) || (i == domhi[0])) {
          continue; // Doing that to avoid ghost cells already filled by corners
        }
      }
      for (int k = q_lo[2] + 1; k < q_hi[2] - 1; k++) {
        amrex::Real  z = (static_cast<amrex::Real>(k) + 0.5) * dx[2];
        if (nscbc_isAnyPerio == 0) {
          if ((k == domlo[2]) || (k == domhi[2])) {
            continue; // Doing that to avoid ghost cells already filled by
                      // corners
          }
        }
        const amrex::Real x_array[AMREX_SPACEDIM] = {AMREX_D_DECL(x,y,z)};
        
        // Normal derivative along y
        amrex::Real dpdy, dudy, dvdy, dwdy, drhody;
        normal_derivative(i, j, k, 1, 1, dx[1], dpdy, dudy, dvdy, dwdy, drhody, q);

        amrex::Real dpdx, dudx, dvdx, dwdx, drhodx;
        // Tangential derivative along x
        tangential_derivative(i, j, k, 0, dx[0], dpdx, dudx, dvdx, dwdx, drhodx, q);

        amrex::Real dpdz, dudz, dvdz, dwdz, drhodz;
        // Tangential derivative along z
        tangential_derivative(i, j, k, 2, dx[2], dpdz, dudz, dvdz, dwdz, drhodz, q);


        amrex::GpuArray<amrex::Real, 5> Ty = {{0.0}};
        // Compute transverse terms
        compute_transverse_terms(
          i, j, k, 1, Ty.data(), dpdx, dudx, dvdx, dwdx,
          drhodx, dpdy, dudy, dvdy, dwdy, drhody, dpdz, dudz, dvdz, dwdz,
          drhodz, q, qaux);

        amrex::GpuArray<amrex::Real, NVAR> y_bc_target;
        amrex::GpuArray<amrex::Real, NVAR> s_int;

        // Filling bcMask with specific user defined BC type
        // bcnormal([x,y,z],U_dummy,U_ext,1,1,time,bc_type,bc_params,bc_target);
        bcnormal(x_array, s_int.data(), y_bc_target.data(), 1, 1, time, geom.data(), *lprobparm);
        int y_bc_type = 7; // Hard-coded inflow. This variable should be updated in bcnormal()

        if (
          (i < q_lo[0] + 3) || (i > q_hi[0] - 3) || (k < q_lo[2] + 3) ||
          (k > q_hi[2] - 3)) {
          // do nothing
          // There is just 1 ghost-cell with bcMask because of the Riemann
          // solver
        } else {
          y_bcMask(i, j, k) = y_bc_type;
        }

        // LODI system waves for y
        amrex::GpuArray<amrex::Real, 5> Ly = {{0.0}};

        // Computing the LODI system waves
        compute_waves(i, j, k, 1, 1, y_bc_type, problen.data(), bc_params.data(),
                    y_bc_target.data(), Ty.data(), Ly.data(), dpdy, dudy, dvdy, dwdy, drhody, q, qaux);

        // Recomputing ghost-cells values with the LODI waves
        update_ghost_cells(i, j, k, y_bc_type, 1, 1, dx[1], domlo, domhi, Ly.data(), uin, q, qaux);
      }
    }
  }

  //--------------------------------------------------------------------------
  // upper Y
  //--------------------------------------------------------------------------

  if ((q_hi[1] > domhi[1]) && (PeleC::phys_bc.hi()[1] == UserBC)) {
    int j = domhi[1];
    amrex::Real  y = (static_cast<amrex::Real>(j) + 0.5) * dx[1];
    for (int i = q_lo[0] + 1; i < q_hi[0] - 1; i++) {
      amrex::Real  x = (static_cast<amrex::Real>(i) + 0.5) * dx[0];
      if (nscbc_isAnyPerio == 0) {
        if ((i == domlo[0]) || (i == domhi[0])) {
          continue; // Doing that to avoid ghost cells already filled by corners
        }
      }
      for (int k = q_lo[2] + 1; k < q_hi[2] - 1; k++) {
        amrex::Real  z = (static_cast<amrex::Real>(k) + 0.5) * dx[2];
        if (nscbc_isAnyPerio == 0) {
          if ((k == domlo[2]) || (k == domhi[2])) {
            continue; // Doing that to avoid ghost cells already filled by
                      // corners
          }
        }
        const amrex::Real x_array[AMREX_SPACEDIM] = {AMREX_D_DECL(x,y,z)};
        
        // Normal derivative along y
        amrex::Real dpdy, dudy, dvdy, dwdy, drhody;
        normal_derivative(i, j, k, 1, -1, dx[1], dpdy, dudy, dvdy, dwdy, drhody, q);

        amrex::Real dpdx, dudx, dvdx, dwdx, drhodx;
        // Tangential derivative along x
        tangential_derivative(i, j, k, 0, dx[0], dpdx, dudx, dvdx, dwdx, drhodx, q);

        amrex::Real dpdz, dudz, dvdz, dwdz, drhodz;
        // Tangential derivative along z
        tangential_derivative(i, j, k, 2, dx[2], dpdz, dudz, dvdz, dwdz, drhodz, q);

        amrex::GpuArray<amrex::Real, 5> Ty = {{0.0}};
        // Compute transverse terms
        compute_transverse_terms(
          i, j, k, 1, Ty.data(), dpdx, dudx, dvdx, dwdx,
          drhodx, dpdy, dudy, dvdy, dwdy, drhody, dpdz, dudz, dvdz, dwdz,
          drhodz, q, qaux);

        amrex::GpuArray<amrex::Real, NVAR> y_bc_target;
        amrex::GpuArray<amrex::Real, NVAR> s_int;

        // Filling bcMask with specific user defined BC type
        // bcnormal([x,y,z],U_dummy,U_ext,2,-1,time,bc_type,bc_params,bc_target);
        bcnormal(x_array, s_int.data(), y_bc_target.data(), 1, -1, time, geom.data(), *lprobparm);
        int y_bc_type = 8; // Hard-coded outlow. This variable should be updated in bcnormal()

        if (
          (i < q_lo[0] + 3) || (i > q_hi[0] - 3) || (k < q_lo[2] + 3) ||
          (k > q_hi[2] - 3)) {
          // do nothing
          // There is just 1 ghost-cell with bcMask because of the Riemann
          // solver
        } else {
          y_bcMask(i, j + 1, k) = y_bc_type;
        }

        // LODI system waves for y
        amrex::GpuArray<amrex::Real, 5> Ly = {{0.0}};

        // Computing the LODI system waves
        compute_waves(i, j, k, 1, -1, y_bc_type, problen.data(), bc_params.data(),
                    y_bc_target.data(), Ty.data(), Ly.data(), dpdy, dudy, dvdy, dwdy, drhody, q, qaux);

        // Recomputing ghost-cells values with the LODI waves
        update_ghost_cells(i, j, k, y_bc_type, 1, -1, dx[1], domlo, domhi, Ly.data(), uin, q, qaux);

      }
    }
  }

  //--------------------------------------------------------------------------
  // lower Z
  //--------------------------------------------------------------------------

  if ((q_lo[2] < domlo[2]) && (PeleC::phys_bc.lo()[2] == UserBC)) {
    int k = domlo[2];
    amrex::Real z = (static_cast<amrex::Real>(k) + 0.5) * dx[2];
    for (int i = q_lo[0] + 1; i < q_hi[0] - 1; i++) {
      amrex::Real x = (static_cast<amrex::Real>(i) + 0.5) * dx[0];
      if (nscbc_isAnyPerio == 0) {
        if ((i == domlo[0]) || (i == domhi[0])) {
          continue; // Doing that to avoid ghost cells already filled by corners
        }
      }
      for (int j = q_lo[1] + 1; j < q_hi[1] - 1; j++) {
        amrex::Real  y = (static_cast<amrex::Real>(j) + 0.5) * dx[1];
        if (nscbc_isAnyPerio == 0) {
          if ((j == domlo[1]) || (j == domhi[1])) {
            continue; // Doing that to avoid ghost cells already filled by
                      // corners
          }
        }
        const amrex::Real x_array[AMREX_SPACEDIM] = {AMREX_D_DECL(x,y,z)};

        // Normal derivative along z
        amrex::Real dpdz, dudz, dvdz, dwdz, drhodz;
        normal_derivative(i, j, k, 2, 1, dx[1], dpdz, dudz, dvdz, dwdz, drhodz, q);

        amrex::Real dpdx, dudx, dvdx, dwdx, drhodx;
        // Tangential derivative along x
        tangential_derivative(i, j, k, 0, dx[0], dpdx, dudx, dvdx, dwdx, drhodx, q);

        amrex::Real dpdy, dudy, dvdy, dwdy, drhody;
        // Tangential derivative along y
        tangential_derivative(i, j, k, 1, dx[1], dpdy, dudy, dvdy, dwdy, drhody, q);

        amrex::GpuArray<amrex::Real, 5> Tz = {{0.0}};
        // Compute transverse terms
        compute_transverse_terms(
          i, j, k, 2, Tz.data(), dpdx, dudx, dvdx, dwdx, drhodx, dpdy, dudy, dvdy, dwdy, drhody,
          dpdz, dudz, dvdz, dwdz, drhodz, q, qaux); 

        amrex::GpuArray<amrex::Real, NVAR> z_bc_target;
        amrex::GpuArray<amrex::Real, NVAR> s_int;

        // Filling bcMask with specific user defined BC type
        // bcnormal([x,y,z],U_dummy,U_ext,3,1,time,bc_type,bc_params,bc_target);
        bcnormal(x_array, s_int.data(), z_bc_target.data(), 2, 1, time, geom.data(), *lprobparm);
        int z_bc_type = 7; // Hard-coded inflow. This variable should be updated in bcnormal()

        if (
          (i < q_lo[0] + 3) || (i > q_hi[0] - 3) || (j < q_lo[1] + 3) ||
          (j > q_hi[1] - 3)) {
          // do nothing
          // There is just 1 ghost-cell with bcMask because of the Riemann
          // solver
        } else {
          z_bcMask(i, j, k) = z_bc_type;
        }

        // LODI system waves for z
        amrex::GpuArray<amrex::Real, 5> Lz = {{0.0}};

        // Computing the LODI system waves
        compute_waves(i, j, k, 2, 1, z_bc_type, problen.data(), bc_params.data(),
          z_bc_target.data(), Tz.data(), Lz.data(), dpdz, dudz, dvdz, dwdz, drhodz, q, qaux);

        // Recomputing ghost-cells values with the LODI waves
        update_ghost_cells(i, j, k, z_bc_type, 2, 1, dx[2], domlo, domhi, Lz.data(), uin, q, qaux);

      }
    }
  }

  //--------------------------------------------------------------------------
  // upper Z
  //--------------------------------------------------------------------------

  if ((q_hi[2] > domhi[2]) && (PeleC::phys_bc.hi()[2] == UserBC)) {
    int k = domhi[2];
    amrex::Real  z = (static_cast<amrex::Real>(k) + 0.5) * dx[2];
    for (int i = q_lo[0] + 1; i < q_hi[0] - 1; i++) {
      amrex::Real  x = (static_cast<amrex::Real>(i) + 0.5) * dx[0];
      if (nscbc_isAnyPerio == 0) {
        if ((i == domlo[0]) || (i == domhi[0])) {
          continue; // Doing that to avoid ghost cells already filled by corners
        }
      }
      for (int j = q_lo[1] + 1; j < q_hi[1] - 1; j++) {
        amrex::Real  y = (static_cast<amrex::Real>(j) + 0.5) * dx[1];
        if (nscbc_isAnyPerio == 0) {
          if ((j == domlo[1]) || (j == domhi[1])) {
            continue; // Doing that to avoid ghost cells already filled by
                      // corners
          }
        }
        const amrex::Real x_array[AMREX_SPACEDIM] = {AMREX_D_DECL(x,y,z)};

        // Normal derivative along z
        amrex::Real dpdz, dudz, dvdz, dwdz, drhodz;
        normal_derivative(i, j, k, 2, -1, dx[1], dpdz, dudz, dvdz, dwdz, drhodz, q);

        amrex::Real dpdx, dudx, dvdx, dwdx, drhodx;
        // Tangential derivative along x
        tangential_derivative(i, j, k, 0, dx[0], dpdx, dudx, dvdx, dwdx, drhodx, q);

        amrex::Real dpdy, dudy, dvdy, dwdy, drhody;
        // Tangential derivative along y
        tangential_derivative(i, j, k, 1, dx[1], dpdy, dudy, dvdy, dwdy, drhody, q);

        amrex::GpuArray<amrex::Real, 5> Tz = {{0.0}};
        // Compute transverse terms
        compute_transverse_terms(
          i, j, k, 2, Tz.data(), dpdx, dudx, dvdx, dwdx, drhodx, dpdy, dudy, dvdy, dwdy, drhody,
          dpdz, dudz, dvdz, dwdz, drhodz, q, qaux); 

        amrex::GpuArray<amrex::Real, NVAR> z_bc_target;
        amrex::GpuArray<amrex::Real, NVAR> s_int;

        // Filling bcMask with specific user defined BC type
        // bcnormal([x,y,z],U_dummy,U_ext,3,-1,time,bc_type,bc_params,bc_target)
        bcnormal(x_array, s_int.data(), z_bc_target.data(), 2, -1, time, geom.data(), *lprobparm);
        int z_bc_type = 8; // Hard-coded outflow. This variable should be updated in bcnormal()

        if (
          (i < q_lo[0] + 3) || (i > q_hi[0] - 3) || (j < q_lo[1] + 3) ||
          (j > q_hi[1] - 3)) {
          // do nothing
          // There is just 1 ghost-cell with bcMask because of the Riemann
          // solver
        } else {
          z_bcMask(i, j, k) = z_bc_type;
        }

        // LODI system waves for z
        amrex::GpuArray<amrex::Real, 5> Lz = {{0.0}};

        // Computing the LODI system waves
        compute_waves(i, j, k, 2, -1, z_bc_type, problen.data(), bc_params.data(),
          z_bc_target.data(), Tz.data(), Lz.data(), dpdz, dudz, dvdz, dwdz, drhodz, q, qaux);

        // Recomputing ghost-cells values with the LODI waves
        update_ghost_cells(i, j, k, z_bc_type, 2, -1, dx[2], domlo, domhi, Lz.data(), uin, q, qaux);
      }
    }
  }
}

void
PeleC::set_bc_mask(
  const amrex::Box& bx,
  const int nscbc_isAnyPerio,
  AMREX_D_DECL(
    const amrex::Array4<int>& x_bcMask,
    const amrex::Array4<int>& y_bcMask,
    const amrex::Array4<int>& z_bcMask),
  AMREX_D_DECL(
    const amrex::Box& x_bcMask_bx,
    const amrex::Box& y_bcMask_bx,
    const amrex::Box& z_bcMask_bx))
{
  const int* qlo = bx.loVect();
  const int* qhi = bx.hiVect();

  const amrex::Box dom = geom.Domain();
  const int* domlo = dom.loVect();
  const int* domhi = dom.hiVect();
  const int domlox = domlo[0];
  const int domhix = domhi[0];
  if(domlo[0] < qlo[0] && domlo[1] < qlo[1] && domhi[0] > qhi[0] &&
     domhi[1] > qhi[1] && domlo[2] < qlo[2] && domhi[2] > qhi[2])
    return;

  // Grab the BCs
  const auto& bcs = PeleC::phys_bc;

  amrex::Box dbox(x_bcMask_bx);
  if (dbox.loVect()[0] == domlo[0]) {
    amrex::Box cbox(amrex::IntVect(AMREX_D_DECL(domlo[0], dbox.loVect()[1], dbox.loVect()[2])),
                    amrex::IntVect(AMREX_D_DECL(domlo[0], dbox.hiVect()[1], dbox.hiVect()[2])));
    amrex::ParallelFor(cbox, [=] AMREX_GPU_DEVICE(int i, int j, int k)
      {
        x_bcMask(i, j, k) = bcs.lo(0);
      });
  }

  if (dbox.hiVect()[0] == domhi[0] + 1) {
    amrex::Box cbox(amrex::IntVect(AMREX_D_DECL(domhi[0] + 1, dbox.loVect()[1], dbox.loVect()[2])),
                    amrex::IntVect(AMREX_D_DECL(domhi[0] + 1, dbox.hiVect()[1], dbox.hiVect()[2])));
    amrex::ParallelFor(cbox, [=] AMREX_GPU_DEVICE(int i, int j, int k)
      {
        x_bcMask(i, j, k) = bcs.hi(0);
      });
  }

  dbox = y_bcMask_bx;
  if (dbox.loVect()[1] == domlo[1]) {
    amrex::Box cbox(amrex::IntVect(AMREX_D_DECL(dbox.loVect()[0], domlo[1], dbox.loVect()[2])),
                    amrex::IntVect(AMREX_D_DECL(dbox.hiVect()[0], domlo[1], dbox.hiVect()[2])));
    amrex::ParallelFor(cbox, [=] AMREX_GPU_DEVICE(int i, int j, int k)
      {
        y_bcMask(i, j, k) = bcs.lo(1);
      });
  }

  if (dbox.hiVect()[1] == domhi[1] + 1) {
    amrex::Box cbox(amrex::IntVect(AMREX_D_DECL(dbox.loVect()[0], domhi[1] + 1, dbox.loVect()[2])),
                    amrex::IntVect(AMREX_D_DECL(dbox.hiVect()[0], domhi[1] + 1, dbox.hiVect()[2])));
    amrex::ParallelFor(cbox, [=] AMREX_GPU_DEVICE(int i, int j, int k)
      {
        y_bcMask(i, j, k) = bcs.hi(1);
      });
  }

  dbox = z_bcMask_bx;
  if (dbox.loVect()[2] == domlo[2]) {
    amrex::Box cbox(amrex::IntVect(AMREX_D_DECL(dbox.loVect()[0], dbox.loVect()[1], domlo[2])),
                    amrex::IntVect(AMREX_D_DECL(dbox.hiVect()[0], dbox.hiVect()[1], domlo[2])));
    amrex::ParallelFor(cbox, [=] AMREX_GPU_DEVICE(int i, int j, int k)
      {
        z_bcMask(i, j, k) = bcs.lo(2);
      });
  }

  if (dbox.hiVect()[2] == domhi[2] + 1) {
    amrex::Box cbox(amrex::IntVect(AMREX_D_DECL(dbox.loVect()[0], dbox.loVect()[1], domhi[2] + 1)),
                    amrex::IntVect(AMREX_D_DECL(dbox.hiVect()[0], dbox.hiVect()[1], domhi[2] + 1)));
    amrex::ParallelFor(cbox, [=] AMREX_GPU_DEVICE(int i, int j, int k)
      {
        z_bcMask(i, j, k) = bcs.hi(2);
      });
  }
}
