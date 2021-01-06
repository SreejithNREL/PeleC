//***** THIS IS NOT FINISHED FROM BEING PORTED FROM FORTRAN YET *******

#include "PeleC.H"

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
  int lo[3],
  int hi[3],
  int domlo[3],
  int domhi[3],
  const amrex::Array4<amrex::Real>& uin,
  int uin_l1,
  int uin_l2,
  int uin_l3,
  int uin_h1,
  int uin_h2,
  int uin_h3,
  const amrex::Array4<amrex::Real>& q,
  int q_l1,
  int q_l2,
  int q_l3,
  int q_h1,
  int q_h2,
  int q_h3,
  const amrex::Array4<amrex::Real>& qaux,
  int qa_l1,
  int qa_l2,
  int qa_l3,
  int qa_h1,
  int qa_h2,
  int qa_h3,
  const amrex::Array4<int>& x_bcMask,
  int x_bcMask_l1,
  int x_bcMask_l2,
  int x_bcMask_l3,
  int x_bcMask_h1,
  int x_bcMask_h2,
  int x_bcMask_h3,
  const amrex::Array4<int>& y_bcMask,
  int y_bcMask_l1,
  int y_bcMask_l2,
  int y_bcMask_l3,
  int y_bcMask_h1,
  int y_bcMask_h2,
  int y_bcMask_h3,
  const amrex::Array4<int>& z_bcMask,
  int z_bcMask_l1,
  int z_bcMask_l2,
  int z_bcMask_l3,
  int z_bcMask_h1,
  int z_bcMask_h2,
  int z_bcMask_h3,
  int flag_nscbc_isAnyPerio,
  int flag_nscbc_perio,
  amrex::Real time,
  amrex::Real delta[3],
  amrex::Real dt,
  int verbose)
{
  // use amrex_constants_module
  // use prob_params_module, only : physbc_lo, physbc_hi, problo, probhi, UserBC
  // use meth_params_module, only : NVAR, NQAUX,QVAR
  // use bc_fill_module, only: bcnormal
  // implicit none
  // integer, intent(in) :: lo(3), hi(3), verbose
  // integer, intent(in) :: domlo(3), domhi(3)
  // integer, intent(in) :: q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
  // integer, intent(in) :: qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3
  // integer, intent(in) :: uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3
  // integer, intent(in) :: x_bcMask_l1, x_bcMask_l2, x_bcMask_l3, x_bcMask_h1,
  // x_bcMask_h2, x_bcMask_h3 integer, intent(in) :: y_bcMask_l1, y_bcMask_l2,
  // y_bcMask_l3, y_bcMask_h1, y_bcMask_h2, y_bcMask_h3 integer, intent(in) ::
  // z_bcMask_l1, z_bcMask_l2, z_bcMask_l3, z_bcMask_h1, z_bcMask_h2,
  // z_bcMask_h3 integer, intent(in) :: flag_nscbc_isAnyPerio integer,
  // intent(in) :: flag_nscbc_perio(3) amrex::Real, intent(inout) ::
  // q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR) amrex::Real, intent(inout) ::
  // qaux(qa_l1:qa_h1,qa_l2:qa_h2,qa_l3:qa_h3,NQAUX) amrex::Real, intent(inout)
  // :: uin(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,NVAR) integer,
  // intent(inout) ::
  // x_bcMask(x_bcMask_l1:x_bcMask_h1,x_bcMask_l2:x_bcMask_h2,x_bcMask_l3:x_bcMask_h3)
  // integer, intent(inout) ::
  // y_bcMask(y_bcMask_l1:y_bcMask_h1,y_bcMask_l2:y_bcMask_h2,y_bcMask_l3:y_bcMask_h3)
  // integer, intent(inout) ::
  // z_bcMask(z_bcMask_l1:z_bcMask_h1,z_bcMask_l2:z_bcMask_h2,z_bcMask_l3:z_bcMask_h3)
  // amrex::Real, intent(in) :: delta(3), dt, time

  // Local
  amrex::Real dx, dy, dz;
  amrex::Real x, y, z;
  amrex::Real drhodx, dudx, dvdx, dwdx, dpdx;
  amrex::Real dpdy, dudy, dvdy, dwdy, drhody;
  amrex::Real dpdz, dudz, dvdz, dwdz, drhodz;
  amrex::Real L1, L2, L3, L4, L5;
  amrex::Real M1, M2, M3, M4, M5;
  amrex::Real N1, N2, N3, N4, N5;
  amrex::Real T1_X, T2_X, T3_X, T4_X, T5_X;
  amrex::Real T1_Y, T2_Y, T3_Y, T4_Y, T5_Y;
  amrex::Real T1_Z, T2_Z, T3_Z, T4_Z, T5_Z;
  amrex::Real U_dummy[NVAR];
  amrex::Real U_ext[NVAR];
  const amrex::Real small = 1.0e-8;
  int q_lo[3] = {q_l1, q_l2, q_l3};
  int q_hi[3] = {q_h1, q_h2, q_h3};
  int uin_lo[3] = {uin_l1, uin_l2, uin_l3};
  int uin_hi[3] = {uin_h1, uin_h2, uin_h3};
  int i, j, k;
  int bc_type, x_bc_type, y_bc_type, z_bc_type;
  int x_isign, y_isign, z_isign, x_idx_Mask, y_idx_Mask, z_idx_Mask;
  int test_keyword_x, test_keyword_y, test_keyword_z;
  amrex::Real bc_params[6] = {0.0};
  amrex::Real x_bc_params[6] = {0.0};
  amrex::Real y_bc_params[6] = {0.0};
  amrex::Real z_bc_params[6] = {0.0};
  amrex::Real bc_target[5] = {0.0};
  amrex::Real x_bc_target[5] = {0.0};
  amrex::Real y_bc_target[5] = {0.0};
  amrex::Real z_bc_target[5] = {0.0};

  dx = delta[0];
  dy = delta[1];
  dz = delta[2];

  /*
      x_bcMask.setVal<amrex::RunOn::Device>(0);
      y_bcMask.setVal<amrex::RunOn::Device>(0);
      z_bcMask.setVal<amrex::RunOn::Device>(0);
  */
  if (flag_nscbc_isAnyPerio == 0) {

    //--------------------------------------------------------------------------
    // corners
    //--------------------------------------------------------------------------

    if (
      ((q_hi[0] > domhi[0]) || (q_lo[0] < domlo[0])) &&
      ((q_hi[1] > domhi[1]) || (q_lo[1] < domlo[1])) &&
      ((q_hi[2] > domhi[2]) || (q_lo[2] < domlo[2]))) {

      if (q_hi[0] > domhi[1]) {
        test_keyword_x = PeleC::phys_bc.hi()[0];
        i = domhi[0];
        x_isign = -1;
        x_idx_Mask = i + 1;
      } else if (q_lo[0] < domlo[0]) {
        test_keyword_x = PeleC::phys_bc.lo()[0];
        i = domlo[0];
        x_isign = 1;
        x_idx_Mask = i;
      }

      if (q_hi[1] > domhi[1]) {
        test_keyword_y = PeleC::phys_bc.hi()[1];
        j = domhi[1];
        y_isign = -1;
        y_idx_Mask = j + 1;
      } else if (q_lo[1] < domlo[1]) {
        test_keyword_y = PeleC::phys_bc.lo()[1];
        j = domlo[1];
        y_isign = 1;
        y_idx_Mask = j;
      }

      if (q_hi[2] > domhi[2]) {
        test_keyword_z = PeleC::phys_bc.hi()[2];
        k = domhi[2];
        z_isign = -1;
        z_idx_Mask = k + 1;
      } else if (q_lo[2] < domlo[2]) {
        test_keyword_z = PeleC::phys_bc.lo()[2];
        k = domlo[2];
        z_isign = 1;
        z_idx_Mask = k;
      }

      x = (static_cast<amrex::Real>(i) + 0.5) * dx;
      y = (static_cast<amrex::Real>(j) + 0.5) * dy;
      z = (static_cast<amrex::Real>(k) + 0.5) * dz;

      // Normal derivative along x
      normal_derivative(
        i, j, k, 1, x_isign, dx, dpdx, dudx, dvdx, dwdx, drhodx, q, q_l1, q_l2,
        q_l3, q_h1, q_h2, q_h3);

      // Normal derivative along y
      normal_derivative(
        i, j, k, 2, y_isign, dy, dpdy, dudy, dvdy, dwdy, drhody, q, q_l1, q_l2,
        q_l3, q_h1, q_h2, q_h3);

      // Normal derivative along z
      normal_derivative(
        i, j, k, 3, z_isign, dz, dpdz, dudz, dvdz, dwdz, drhodz, q, q_l1, q_l2,
        q_l3, q_h1, q_h2, q_h3);

      // Compute transverse terms for X
      compute_transverse_terms(
        i, j, k, 1, T1_X, T2_X, T3_X, T4_X, T5_X, dpdx, dudx, dvdx, dwdx,
        drhodx, dpdy, dudy, dvdy, dwdy, drhody, dpdz, dudz, dvdz, dwdz, drhodz,
        q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1,
        qa_h2, qa_h3);

      // Compute transverse terms for X
      compute_transverse_terms(
        i, j, k, 2, T1_Y, T2_Y, T3_Y, T4_Y, T5_Y, dpdx, dudx, dvdx, dwdx,
        drhodx, dpdy, dudy, dvdy, dwdy, drhody, dpdz, dudz, dvdz, dwdz, drhodz,
        q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1,
        qa_h2, qa_h3);

      // Compute transverse terms for X
      compute_transverse_terms(
        i, j, k, 3, T1_Z, T2_Z, T3_Z, T4_Z, T5_Z, dpdx, dudx, dvdx, dwdx,
        drhodx, dpdy, dudy, dvdy, dwdy, drhody, dpdz, dudz, dvdz, dwdz, drhodz,
        q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1,
        qa_h2, qa_h3);

      // Calling user target BC values
      // x face
      if (test_keyword_x == UserBC) {
        // bcnormal([x,y,z],U_dummy,U_ext,1,x_isign,time,x_bc_type,x_bc_params,x_bc_target);
      } else {
        x_bc_type = test_keyword_x;
      }
      x_bcMask(x_idx_Mask, j, k) = x_bc_type;

      // y face
      if (test_keyword_y == UserBC) {
        // bcnormal([x,y,z],U_dummy,U_ext,2,y_isign,time,y_bc_type,y_bc_params,y_bc_target);
      } else {
        y_bc_type = test_keyword_y;
      }
      y_bcMask(i, y_idx_Mask, k) = y_bc_type;

      // z face
      if (test_keyword_z == UserBC) {
        // bcnormal([x,y,z],U_dummy,U_ext,3,z_isign,time,z_bc_type,z_bc_params,z_bc_target);
      } else {
        z_bc_type = test_keyword_z;
      }
      z_bcMask(i, j, z_idx_Mask) = z_bc_type;

      // Computing the LODI system waves along X
      compute_waves(
        i, j, k, 1, x_isign, x_bc_type, x_bc_params, x_bc_target, T1_X, T2_X,
        T3_X, T4_X, T5_X, L1, L2, L3, L4, L5, dpdx, dudx, dvdx, dwdx, drhodx, q,
        q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1,
        qa_h2, qa_h3);

      // Computing the LODI system waves along Y
      compute_waves(
        i, j, k, 2, y_isign, y_bc_type, y_bc_params, y_bc_target, T1_Y, T2_Y,
        T3_Y, T4_Y, T5_Y, M1, M2, M3, M4, M5, dpdy, dudy, dvdy, dwdy, drhody, q,
        q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1,
        qa_h2, qa_h3);

      // Computing the LODI system waves along Z
      compute_waves(
        i, j, k, 3, z_isign, z_bc_type, z_bc_params, z_bc_target, T1_Z, T2_Z,
        T3_Z, T4_Z, T5_Z, N1, N2, N3, N4, N5, dpdz, dudz, dvdz, dwdz, drhodz, q,
        q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1,
        qa_h2, qa_h3);

      // Recomputing ghost-cells values with the LODI waves along X
      update_ghost_cells(
        i, j, k, x_bc_type, 1, x_isign, dx, domlo, domhi, L1, L2, L3, L4, L5,
        uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, q, q_l1, q_l2,
        q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3);

      // Recomputing ghost-cells values with the LODI waves along Y
      update_ghost_cells(
        i, j, k, y_bc_type, 2, y_isign, dy, domlo, domhi, M1, M2, M3, M4, M5,
        uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, q, q_l1, q_l2,
        q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3);

      // Recomputing ghost-cells values with the LODI waves along Y
      update_ghost_cells(
        i, j, k, z_bc_type, 3, z_isign, dz, domlo, domhi, N1, N2, N3, N4, N5,
        uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, q, q_l1, q_l2,
        q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3);
    }
  }
  //--------------------------------------------------------------------------
  // lower X
  //--------------------------------------------------------------------------

  if ((q_lo[0] < domlo[0]) && (PeleC::phys_bc.lo()[0] == UserBC)) {
    i = domlo[0];
    x = (static_cast<amrex::Real>(i) + 0.5) * dx;
    for (j = q_lo[1] + 1; j < q_hi[1] - 1; j++) {
      y = (static_cast<amrex::Real>(j) + 0.5) * dy;
      if (flag_nscbc_isAnyPerio == 0) {
        if ((j == domlo[1]) || (j == domhi[1])) {
          continue; // Doing that to avoid ghost cells already filled by corners
        }
      }
      for (k = q_lo[2] + 1; k < q_hi[2] - 1; k++) {
        z = (static_cast<amrex::Real>(k) + 0.5) * dz;
        if (flag_nscbc_isAnyPerio == 0) {
          if ((k == domlo[2]) || (k == domhi[2])) {
            continue; // Doing that to avoid ghost cells already filled by
                      // corners
          }
        }
        // Normal derivative along x
        normal_derivative(
          i, j, k, 1, 1, dx, dpdx, dudx, dvdx, dwdx, drhodx, q, q_l1, q_l2,
          q_l3, q_h1, q_h2, q_h3);

        // Tangential derivative along y
        tangential_derivative(
          i, j, k, 2, dy, dpdy, dudy, dvdy, dwdy, drhody, q, q_l1, q_l2, q_l3,
          q_h1, q_h2, q_h3);

        // Tangential derivative along z
        tangential_derivative(
          i, j, k, 3, dz, dpdz, dudz, dvdz, dwdz, drhodz, q, q_l1, q_l2, q_l3,
          q_h1, q_h2, q_h3);

        // Compute transverse terms
        compute_transverse_terms(
          i, j, k, 1, T1_X, T2_X, T3_X, T4_X, T5_X, dpdx, dudx, dvdx, dwdx,
          drhodx, dpdy, dudy, dvdy, dwdy, drhody, dpdz, dudz, dvdz, dwdz,
          drhodz, q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2,
          qa_l3, qa_h1, qa_h2, qa_h3);

        // Calling user target BC values
        // bcnormal([x,y,z],U_dummy,U_ext,1,1,time,bc_type,bc_params,bc_target);

        // Filling bcMask with specific user defined BC type
        if (
          (j < q_lo[1] + 3) || (j > q_hi[1] - 3) || (k < q_lo[2] + 3) ||
          (k > q_hi[2] - 3)) {
          // do nothing.
          // There is just 1 ghost-cell with bcMask because of the Riemann
          // solver
        } else {
          x_bcMask(i, j, k) = bc_type;
        }

        // Computing the LODI system waves
        compute_waves(
          i, j, k, 1, 1, bc_type, bc_params, bc_target, T1_X, T2_X, T3_X, T4_X,
          T5_X, L1, L2, L3, L4, L5, dpdx, dudx, dvdx, dwdx, drhodx, q, q_l1,
          q_l2, q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2,
          qa_h3);

        // Recomputing ghost-cells values with the LODI waves
        update_ghost_cells(
          i, j, k, bc_type, 1, 1, dx, domlo, domhi, L1, L2, L3, L4, L5, uin,
          uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, q, q_l1, q_l2, q_l3,
          q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3);
      }
    }
  }

  //--------------------------------------------------------------------------
  // upper X
  //--------------------------------------------------------------------------

  if ((q_hi[0] > domhi[0]) && (PeleC::phys_bc.hi()[0] == UserBC)) {
    i = domhi[0];
    x = (static_cast<amrex::Real>(i) + 0.5) * dx;
    for (j = q_lo[1] + 1; j < q_hi[1] - 1; j++) {
      y = (static_cast<amrex::Real>(j) + 0.5) * dy;
      if (flag_nscbc_isAnyPerio == 0) {
        if ((j == domlo[1]) || (j == domhi[1])) {
          continue; // Doing that to avoid ghost cells already filled by corners
        }
      }
      for (k = q_lo[2] + 1; k < q_hi[2] - 1; k++) {
        z = (static_cast<amrex::Real>(k) + 0.5) * dz;
        if (flag_nscbc_isAnyPerio == 0) {
          if ((k == domlo[2]) || (k == domhi[2])) {
            continue; // Doing that to avoid ghost cells already filled by
                      // corners
          }
        }

        // Normal derivative along x
        normal_derivative(
          i, j, k, 1, -1, dx, dpdx, dudx, dvdx, dwdx, drhodx, q, q_l1, q_l2,
          q_l3, q_h1, q_h2, q_h3);

        // Tangential derivative along y
        tangential_derivative(
          i, j, k, 2, dy, dpdy, dudy, dvdy, dwdy, drhody, q, q_l1, q_l2, q_l3,
          q_h1, q_h2, q_h3);

        // Tangential derivative along z
        tangential_derivative(
          i, j, k, 3, dz, dpdz, dudz, dvdz, dwdz, drhodz, q, q_l1, q_l2, q_l3,
          q_h1, q_h2, q_h3);

        // Compute transverse terms
        compute_transverse_terms(
          i, j, k, 1, T1_X, T2_X, T3_X, T4_X, T5_X, dpdx, dudx, dvdx, dwdx,
          drhodx, dpdy, dudy, dvdy, dwdy, drhody, dpdz, dudz, dvdz, dwdz,
          drhodz, q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2,
          qa_l3, qa_h1, qa_h2, qa_h3);

        // Filling bcMask with specific user defined BC type
        // bcnormal([x,y,z],U_dummy,U_ext,1,-1,time,bc_type,bc_params,bc_target);
        if (
          (j < q_lo[1] + 3) || (j > q_hi[1] - 3) || (k < q_lo[2] + 3) ||
          (k > q_hi[2] - 3)) {
          // do nothing
          // There is just 1 ghost-cell with bcMask because of the Riemann
          // solver
        } else {
          x_bcMask(i + 1, j, k) = bc_type;
        }

        // Computing the LODI system waves
        compute_waves(
          i, j, k, 1, -1, bc_type, bc_params, bc_target, T1_X, T2_X, T3_X, T4_X,
          T5_X, L1, L2, L3, L4, L5, dpdx, dudx, dvdx, dwdx, drhodx, q, q_l1,
          q_l2, q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2,
          qa_h3);

        // Recomputing ghost-cells values with the LODI waves
        update_ghost_cells(
          i, j, k, bc_type, 1, -1, dx, domlo, domhi, L1, L2, L3, L4, L5, uin,
          uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, q, q_l1, q_l2, q_l3,
          q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3);
      }
    }
  }

  //--------------------------------------------------------------------------
  // lower Y
  //--------------------------------------------------------------------------

  if ((q_lo[1] < domlo[1]) && (PeleC::phys_bc.lo()[1] == UserBC)) {
    j = domlo[1];
    y = (static_cast<amrex::Real>(j) + 0.5) * dy;
    for (i = q_lo[0] + 1; i < q_hi[0] - 1; i++) {
      x = (static_cast<amrex::Real>(i) + 0.5) * dx;
      if (flag_nscbc_isAnyPerio == 0) {
        if ((i == domlo[0]) || (i == domhi[0])) {
          continue; // Doing that to avoid ghost cells already filled by corners
        }
      }
      for (k = q_lo[2] + 1; k < q_hi[2] - 1; k++) {
        z = (static_cast<amrex::Real>(k) + 0.5) * dz;
        if (flag_nscbc_isAnyPerio == 0) {
          if ((k == domlo[2]) || (k == domhi[2])) {
            continue; // Doing that to avoid ghost cells already filled by
                      // corners
          }
        }

        // Normal derivative along y
        normal_derivative(
          i, j, k, 2, 1, dy, dpdy, dudy, dvdy, dwdy, drhody, q, q_l1, q_l2,
          q_l3, q_h1, q_h2, q_h3);

        // Tangential derivative along x
        tangential_derivative(
          i, j, k, 1, dx, dpdx, dudx, dvdx, dwdx, drhodx, q, q_l1, q_l2, q_l3,
          q_h1, q_h2, q_h3);

        // Tangential derivative along z
        tangential_derivative(
          i, j, k, 3, dz, dpdz, dudz, dvdz, dwdz, drhodz, q, q_l1, q_l2, q_l3,
          q_h1, q_h2, q_h3);

        // Compute transverse terms
        compute_transverse_terms(
          i, j, k, 2, T1_Y, T2_Y, T3_Y, T4_Y, T5_Y, dpdx, dudx, dvdx, dwdx,
          drhodx, dpdy, dudy, dvdy, dwdy, drhody, dpdz, dudz, dvdz, dwdz,
          drhodz, q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2,
          qa_l3, qa_h1, qa_h2, qa_h3);

        // Filling bcMask with specific user defined BC type
        // bcnormal([x,y,z],U_dummy,U_ext,2,1,time,bc_type,bc_params,bc_target);
        if (
          (i < q_lo[0] + 3) || (i > q_hi[0] - 3) || (k < q_lo[2] + 3) ||
          (k > q_hi[2] - 3)) {
          // do nothing
          // There is just 1 ghost-cell with bcMask because of the Riemann
          // solver
        } else {
          y_bcMask(i, j, k) = bc_type;
        }

        // Computing the LODI system waves
        compute_waves(
          i, j, k, 2, 1, bc_type, bc_params, bc_target, T1_Y, T2_Y, T3_Y, T4_Y,
          T5_Y, L1, L2, L3, L4, L5, dpdy, dudy, dvdy, dwdy, drhody, q, q_l1,
          q_l2, q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2,
          qa_h3);

        // Recomputing ghost-cells values with the LODI waves
        update_ghost_cells(
          i, j, k, bc_type, 2, 1, dy, domlo, domhi, L1, L2, L3, L4, L5, uin,
          uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, q, q_l1, q_l2, q_l3,
          q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3);
      }
    }
  }

  //--------------------------------------------------------------------------
  // upper Y
  //--------------------------------------------------------------------------

  if ((q_hi[1] > domhi[1]) && (PeleC::phys_bc.hi()[1] == UserBC)) {
    j = domhi[1];
    y = (static_cast<amrex::Real>(j) + 0.5) * dy;
    for (i = q_lo[0] + 1; i < q_hi[0] - 1; i++) {
      x = (static_cast<amrex::Real>(i) + 0.5) * dx;
      if (flag_nscbc_isAnyPerio == 0) {
        if ((i == domlo[0]) || (i == domhi[0])) {
          continue; // Doing that to avoid ghost cells already filled by corners
        }
      }
      for (k = q_lo[2] + 1; k < q_hi[2] - 1; k++) {
        z = (static_cast<amrex::Real>(k) + 0.5) * dz;
        if (flag_nscbc_isAnyPerio == 0) {
          if ((k == domlo[2]) || (k == domhi[2])) {
            continue; // Doing that to avoid ghost cells already filled by
                      // corners
          }
        }

        // Normal derivative along y
        normal_derivative(
          i, j, k, 2, -1, dy, dpdy, dudy, dvdy, dwdy, drhody, q, q_l1, q_l2,
          q_l3, q_h1, q_h2, q_h3);

        // Tangential derivative along x
        tangential_derivative(
          i, j, k, 1, dx, dpdx, dudx, dvdx, dwdx, drhodx, q, q_l1, q_l2, q_l3,
          q_h1, q_h2, q_h3);

        // Tangential derivative along z
        tangential_derivative(
          i, j, k, 3, dz, dpdz, dudz, dvdz, dwdz, drhodz, q, q_l1, q_l2, q_l3,
          q_h1, q_h2, q_h3);

        // Compute transverse terms
        compute_transverse_terms(
          i, j, k, 2, T1_Y, T2_Y, T3_Y, T4_Y, T5_Y, dpdx, dudx, dvdx, dwdx,
          drhodx, dpdy, dudy, dvdy, dwdy, drhody, dpdz, dudz, dvdz, dwdz,
          drhodz, q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2,
          qa_l3, qa_h1, qa_h2, qa_h3);

        // Filling bcMask with specific user defined BC type
        // bcnormal([x,y,z],U_dummy,U_ext,2,-1,time,bc_type,bc_params,bc_target);

        if (
          (i < q_lo[0] + 3) || (i > q_hi[0] - 3) || (k < q_lo[2] + 3) ||
          (k > q_hi[2] - 3)) {
          // do nothing
          // There is just 1 ghost-cell with bcMask because of the Riemann
          // solver
        } else {
          y_bcMask(i, j + 1, k) = bc_type;
        }

        // Computing the LODI system waves
        compute_waves(
          i, j, k, 2, -1, bc_type, bc_params, bc_target, T1_Y, T2_Y, T3_Y, T4_Y,
          T5_Y, L1, L2, L3, L4, L5, dpdy, dudy, dvdy, dwdy, drhody, q, q_l1,
          q_l2, q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2,
          qa_h3);

        // Recomputing ghost-cells values with the LODI waves
        update_ghost_cells(
          i, j, k, bc_type, 2, -1, dy, domlo, domhi, L1, L2, L3, L4, L5, uin,
          uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, q, q_l1, q_l2, q_l3,
          q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3);
      }
    }
  }

  //--------------------------------------------------------------------------
  // lower Z
  //--------------------------------------------------------------------------

  if ((q_lo[2] < domlo[2]) && (PeleC::phys_bc.lo()[2] == UserBC)) {
    k = domlo[2];
    z = (static_cast<amrex::Real>(k) + 0.5) * dz;
    for (i = q_lo[0] + 1; i < q_hi[0] - 1; i++) {
      x = (static_cast<amrex::Real>(i) + 0.5) * dx;
      if (flag_nscbc_isAnyPerio == 0) {
        if ((i == domlo[0]) || (i == domhi[0])) {
          continue; // Doing that to avoid ghost cells already filled by corners
        }
      }
      for (j = q_lo[1] + 1; j < q_hi[1] - 1; j++) {
        y = (static_cast<amrex::Real>(j) + 0.5) * dy;
        if (flag_nscbc_isAnyPerio == 0) {
          if ((j == domlo[1]) || (j == domhi[1])) {
            continue; // Doing that to avoid ghost cells already filled by
                      // corners
          }
        }

        // Normal derivative along y
        normal_derivative(
          i, j, k, 3, 1, dz, dpdz, dudz, dvdz, dwdz, drhodz, q, q_l1, q_l2,
          q_l3, q_h1, q_h2, q_h3);

        // Tangential derivative along x
        tangential_derivative(
          i, j, k, 1, dx, dpdx, dudx, dvdx, dwdx, drhodx, q, q_l1, q_l2, q_l3,
          q_h1, q_h2, q_h3);

        // Tangential derivative along z
        tangential_derivative(
          i, j, k, 2, dy, dpdy, dudy, dvdy, dwdy, drhody, q, q_l1, q_l2, q_l3,
          q_h1, q_h2, q_h3);

        // Compute transverse terms
        compute_transverse_terms(
          i, j, k, 3, T1_Z, T2_Z, T3_Z, T4_Z, T5_Z, dpdx, dudx, dvdx, dwdx,
          drhodx, dpdy, dudy, dvdy, dwdy, drhody, dpdz, dudz, dvdz, dwdz,
          drhodz, q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2,
          qa_l3, qa_h1, qa_h2, qa_h3);

        // Filling bcMask with specific user defined BC type
        // bcnormal([x,y,z],U_dummy,U_ext,3,1,time,bc_type,bc_params,bc_target);

        if (
          (i < q_lo[0] + 3) || (i > q_hi[0] - 3) || (j < q_lo[1] + 3) ||
          (j > q_hi[1] - 3)) {
          // do nothing
          // There is just 1 ghost-cell with bcMask because of the Riemann
          // solver
        } else {
          z_bcMask(i, j, k) = bc_type;
        }

        // Computing the LODI system waves
        compute_waves(
          i, j, k, 3, 1, bc_type, bc_params, bc_target, T1_Z, T2_Z, T3_Z, T4_Z,
          T5_Z, L1, L2, L3, L4, L5, dpdz, dudz, dvdz, dwdz, drhodz, q, q_l1,
          q_l2, q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2,
          qa_h3);

        // Recomputing ghost-cells values with the LODI waves
        update_ghost_cells(
          i, j, k, bc_type, 3, 1, dz, domlo, domhi, L1, L2, L3, L4, L5, uin,
          uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, q, q_l1, q_l2, q_l3,
          q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3);
      }
    }
  }

  //--------------------------------------------------------------------------
  // upper Z
  //--------------------------------------------------------------------------

  if ((q_hi[2] > domhi[2]) && (PeleC::phys_bc.hi()[2] == UserBC)) {
    k = domhi[2];
    z = (static_cast<amrex::Real>(k) + 0.5) * dz;
    for (i = q_lo[0] + 1; i < q_hi[0] - 1; i++) {
      x = (static_cast<amrex::Real>(i) + 0.5) * dx;
      if (flag_nscbc_isAnyPerio == 0) {
        if ((i == domlo[0]) || (i == domhi[0])) {
          continue; // Doing that to avoid ghost cells already filled by corners
        }
      }
      for (j = q_lo[1] + 1; j < q_hi[1] - 1; j++) {
        y = (static_cast<amrex::Real>(j) + 0.5) * dy;
        if (flag_nscbc_isAnyPerio == 0) {
          if ((j == domlo[1]) || (j == domhi[1])) {
            continue; // Doing that to avoid ghost cells already filled by
                      // corners
          }
        }

        // Normal derivative along y
        normal_derivative(
          i, j, k, 3, -1, dz, dpdz, dudz, dvdz, dwdz, drhodz, q, q_l1, q_l2,
          q_l3, q_h1, q_h2, q_h3);

        // Tangential derivative along x
        tangential_derivative(
          i, j, k, 1, dx, dpdx, dudx, dvdx, dwdx, drhodx, q, q_l1, q_l2, q_l3,
          q_h1, q_h2, q_h3);

        // Tangential derivative along z
        tangential_derivative(
          i, j, k, 2, dy, dpdy, dudy, dvdy, dwdy, drhody, q, q_l1, q_l2, q_l3,
          q_h1, q_h2, q_h3);

        // Compute transverse terms
        compute_transverse_terms(
          i, j, k, 3, T1_Z, T2_Z, T3_Z, T4_Z, T5_Z, dpdx, dudx, dvdx, dwdx,
          drhodx, dpdy, dudy, dvdy, dwdy, drhody, dpdz, dudz, dvdz, dwdz,
          drhodz, q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2,
          qa_l3, qa_h1, qa_h2, qa_h3);

        // Filling bcMask with specific user defined BC type
        // bcnormal([x,y,z],U_dummy,U_ext,3,-1,time,bc_type,bc_params,bc_target)
        if (
          (i < q_lo[0] + 3) || (i > q_hi[0] - 3) || (j < q_lo[1] + 3) ||
          (j > q_hi[1] - 3)) {
          // do nothing
          // There is just 1 ghost-cell with bcMask because of the Riemann
          // solver
        } else {
          z_bcMask(i, j, k) = bc_type;
        }

        // Computing the LODI system waves
        compute_waves(
          i, j, k, 3, -1, bc_type, bc_params, bc_target, T1_Z, T2_Z, T3_Z, T4_Z,
          T5_Z, L1, L2, L3, L4, L5, dpdz, dudz, dvdz, dwdz, drhodz, q, q_l1,
          q_l2, q_l3, q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2,
          qa_h3);

        // Recomputing ghost-cells values with the LODI waves
        update_ghost_cells(
          i, j, k, bc_type, 3, -1, dz, domlo, domhi, L1, L2, L3, L4, L5, uin,
          uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, q, q_l1, q_l2, q_l3,
          q_h1, q_h2, q_h3, qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3);
      }
    }
  }
}

//-------------------------------------------------
// Generic routines below
//-------------------------------------------------

void
PeleC::normal_derivative(
  int i,
  int j,
  int k,
  int idir,
  int isign,
  amrex::Real delta,
  amrex::Real dp,
  amrex::Real du,
  amrex::Real dv,
  amrex::Real dw,
  amrex::Real drho,
  const amrex::Array4<amrex::Real>& q,
  int q_l1,
  int q_l2,
  int q_l3,
  int q_h1,
  int q_h2,
  int q_h3)
{
  // use meth_params_module, only : QVAR, QPRES, QU, QV, QW, QRHO
  // implicit none
  // integer, intent(in) :: i,j,k,idir,isign
  // integer, intent(in) :: q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
  // amrex::Real, intent(in) :: delta
  // amrex::Real, intent(in) :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
  // amrex::Real, intent(out) :: dp, du, dv, dw, drho

  if (idir == 1) {
    if (isign == 1) {
      // 2nd order
      dp = ((-3.0 / 2.0) * q(i, j, k, QPRES) + 2.0 * q(i + 1, j, k, QPRES) -
            0.5 * q(i + 2, j, k, QPRES)) /
           delta;
      du = ((-3.0 / 2.0) * q(i, j, k, QU) + 2.0 * q(i + 1, j, k, QU) -
            0.5 * q(i + 2, j, k, QU)) /
           delta;
      dv = ((-3.0 / 2.0) * q(i, j, k, QV) + 2.0 * q(i + 1, j, k, QV) -
            0.5 * q(i + 2, j, k, QV)) /
           delta;
      dw = ((-3.0 / 2.0) * q(i, j, k, QW) + 2.0 * q(i + 1, j, k, QW) -
            0.5 * q(i + 2, j, k, QW)) /
           delta;
      drho = ((-3.0 / 2.0) * q(i, j, k, QRHO) + 2.0 * q(i + 1, j, k, QRHO) -
              0.5 * q(i + 2, j, k, QRHO)) /
             delta;
    } else if (isign == -1) {
      // 2nd order
      dp = ((3.0 / 2.0) * q(i, j, k, QPRES) - 2.0 * q(i - 1, j, k, QPRES) +
            0.5 * q(i - 2, j, k, QPRES)) /
           delta;
      du = ((3.0 / 2.0) * q(i, j, k, QU) - 2.0 * q(i - 1, j, k, QU) +
            0.5 * q(i - 2, j, k, QU)) /
           delta;
      dv = ((3.0 / 2.0) * q(i, j, k, QV) - 2.0 * q(i - 1, j, k, QV) +
            0.5 * q(i - 2, j, k, QV)) /
           delta;
      dw = ((3.0 / 2.0) * q(i, j, k, QW) - 2.0 * q(i - 1, j, k, QW) +
            0.5 * q(i - 2, j, k, QW)) /
           delta;
      drho = ((3.0 / 2.0) * q(i, j, k, QRHO) - 2.0 * q(i - 1, j, k, QRHO) +
              0.5 * q(i - 2, j, k, QRHO)) /
             delta;
    } else {
      amrex::Abort("Problem of isign in impose_NSCBC_3d:normal_derivative");
    }
  } else if (idir == 2) {
    if (isign == 1) {
      // 2nd order
      dp = ((-3.0 / 2.0) * q(i, j, k, QPRES) + 2.0 * q(i, j + 1, k, QPRES) -
            0.5 * q(i, j + 2, k, QPRES)) /
           delta;
      du = ((-3.0 / 2.0) * q(i, j, k, QU) + 2.0 * q(i, j + 1, k, QU) -
            0.5 * q(i, j + 2, k, QU)) /
           delta;
      dv = ((-3.0 / 2.0) * q(i, j, k, QV) + 2.0 * q(i, j + 1, k, QV) -
            0.5 * q(i, j + 2, k, QV)) /
           delta;
      dw = ((-3.0 / 2.0) * q(i, j, k, QW) + 2.0 * q(i, j + 1, k, QW) -
            0.5 * q(i, j + 2, k, QW)) /
           delta;
      drho = ((-3.0 / 2.0) * q(i, j, k, QRHO) + 2.0 * q(i, j + 1, k, QRHO) -
              0.5 * q(i, j + 2, k, QRHO)) /
             delta;
    } else if (isign == -1) {
      // 2nd order
      dp = ((3.0 / 2.0) * q(i, j, k, QPRES) - 2.0 * q(i, j - 1, k, QPRES) +
            0.5 * q(i, j - 2, k, QPRES)) /
           delta;
      du = ((3.0 / 2.0) * q(i, j, k, QU) - 2.0 * q(i, j - 1, k, QU) +
            0.5 * q(i, j - 2, k, QU)) /
           delta;
      dv = ((3.0 / 2.0) * q(i, j, k, QV) - 2.0 * q(i, j - 1, k, QV) +
            0.5 * q(i, j - 2, k, QV)) /
           delta;
      dw = ((3.0 / 2.0) * q(i, j, k, QW) - 2.0 * q(i, j - 1, k, QW) +
            0.5 * q(i, j - 2, k, QW)) /
           delta;
      drho = ((3.0 / 2.0) * q(i, j, k, QRHO) - 2.0 * q(i, j - 1, k, QRHO) +
              0.5 * q(i, j - 2, k, QRHO)) /
             delta;
    } else {
      amrex::Abort("Problem of isign in impose_NSCBC_3d:normal_derivative");
    }
  } else if (idir == 3) {
    if (isign == 1) {
      // 2nd order
      dp = ((-3.0 / 2.0) * q(i, j, k, QPRES) + 2.0 * q(i, j, k + 1, QPRES) -
            0.5 * q(i, j, k + 2, QPRES)) /
           delta;
      du = ((-3.0 / 2.0) * q(i, j, k, QU) + 2.0 * q(i, j, k + 1, QU) -
            0.5 * q(i, j, k + 2, QU)) /
           delta;
      dv = ((-3.0 / 2.0) * q(i, j, k, QV) + 2.0 * q(i, j, k + 1, QV) -
            0.5 * q(i, j, k + 2, QV)) /
           delta;
      dw = ((-3.0 / 2.0) * q(i, j, k, QW) + 2.0 * q(i, j, k + 1, QW) -
            0.5 * q(i, j, k + 2, QW)) /
           delta;
      drho = ((-3.0 / 2.0) * q(i, j, k, QRHO) + 2.0 * q(i, j, k + 1, QRHO) -
              0.5 * q(i, j, k + 2, QRHO)) /
             delta;
    } else if (isign == -1) {
      // 2nd order
      dp = ((3.0 / 2.0) * q(i, j, k, QPRES) - 2.0 * q(i, j, k - 1, QPRES) +
            0.5 * q(i, j, k - 2, QPRES)) /
           delta;
      du = ((3.0 / 2.0) * q(i, j, k, QU) - 2.0 * q(i, j, k - 1, QU) +
            0.5 * q(i, j, k - 2, QU)) /
           delta;
      dv = ((3.0 / 2.0) * q(i, j, k, QV) - 2.0 * q(i, j, k - 1, QV) +
            0.5 * q(i, j, k - 2, QV)) /
           delta;
      dw = ((3.0 / 2.0) * q(i, j, k, QW) - 2.0 * q(i, j, k - 1, QW) +
            0.5 * q(i, j, k - 2, QW)) /
           delta;
      drho = ((3.0 / 2.0) * q(i, j, k, QRHO) - 2.0 * q(i, j, k - 1, QRHO) +
              0.5 * q(i, j, k - 2, QRHO)) /
             delta;
    } else {
      amrex::Abort("Problem of isign in impose_NSCBC_3d:normal_derivative");
    }
  } else {
    amrex::Abort("Problem of idir in impose_NSCBC_3d:normal_derivative");
  }
}

void
PeleC::tangential_derivative(
  int i,
  int j,
  int k,
  int idir,
  amrex::Real delta,
  amrex::Real dp,
  amrex::Real du,
  amrex::Real dv,
  amrex::Real dw,
  amrex::Real drho,
  const amrex::Array4<amrex::Real>& q,
  int q_l1,
  int q_l2,
  int q_l3,
  int q_h1,
  int q_h2,
  int q_h3)
{
  // use meth_params_module, only : QVAR, QPRES, QU, QV, QW, QRHO
  // implicit none
  // integer, intent(in) :: i,j,k,idir
  // integer, intent(in) :: q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
  // amrex::Real, intent(in) :: delta
  // amrex::Real, intent(in) :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
  // amrex::Real, intent(out) :: dp, du, dv, dw, drho

  // Warning, idir means the tangential direction, this is different from 2D
  // (sorry)
  if (idir == 1) {
    // 2nd order Central
    dp = (q(i + 1, j, k, QPRES) - q(i - 1, j, k, QPRES)) / (2.0 * delta);
    du = (q(i + 1, j, k, QU) - q(i - 1, j, k, QU)) / (2.0 * delta);
    dv = (q(i + 1, j, k, QV) - q(i - 1, j, k, QV)) / (2.0 * delta);
    dw = (q(i + 1, j, k, QW) - q(i - 1, j, k, QW)) / (2.0 * delta);
    drho = (q(i + 1, j, k, QRHO) - q(i - 1, j, k, QRHO)) / (2.0 * delta);
  } else if (idir == 2) {
    // 2nd order Central
    dp = (q(i, j + 1, k, QPRES) - q(i, j - 1, k, QPRES)) / (2.0 * delta);
    du = (q(i, j + 1, k, QU) - q(i, j - 1, k, QU)) / (2.0 * delta);
    dv = (q(i, j + 1, k, QV) - q(i, j - 1, k, QV)) / (2.0 * delta);
    dw = (q(i, j + 1, k, QW) - q(i, j - 1, k, QW)) / (2.0 * delta);
    drho = (q(i, j + 1, k, QRHO) - q(i, j - 1, k, QRHO)) / (2.0 * delta);
  } else if (idir == 3) {
    // 2nd order Central
    dp = (q(i, j, k + 1, QPRES) - q(i, j, k - 1, QPRES)) / (2.0 * delta);
    du = (q(i, j, k + 1, QU) - q(i, j, k - 1, QU)) / (2.0 * delta);
    dv = (q(i, j, k + 1, QV) - q(i, j, k - 1, QV)) / (2.0 * delta);
    dw = (q(i, j, k + 1, QW) - q(i, j, k - 1, QW)) / (2.0 * delta);
    drho = (q(i, j, k + 1, QRHO) - q(i, j, k - 1, QRHO)) / (2.0 * delta);
  } else {
    amrex::Abort("Problem of idir in impose_NSCBC_3d:tangential_derivative");
  }
}

void
PeleC::compute_transverse_terms(
  int i,
  int j,
  int k,
  int idir,
  amrex::Real T1,
  amrex::Real T2,
  amrex::Real T3,
  amrex::Real T4,
  amrex::Real T5,
  amrex::Real dpdx,
  amrex::Real dudx,
  amrex::Real dvdx,
  amrex::Real dwdx,
  amrex::Real drhodx,
  amrex::Real dpdy,
  amrex::Real dudy,
  amrex::Real dvdy,
  amrex::Real dwdy,
  amrex::Real drhody,
  amrex::Real dpdz,
  amrex::Real dudz,
  amrex::Real dvdz,
  amrex::Real dwdz,
  amrex::Real drhodz,
  const amrex::Array4<amrex::Real>& q,
  int q_l1,
  int q_l2,
  int q_l3,
  int q_h1,
  int q_h2,
  int q_h3,
  const amrex::Array4<amrex::Real>& qaux,
  int qa_l1,
  int qa_l2,
  int qa_l3,
  int qa_h1,
  int qa_h2,
  int qa_h3)
{
  // use meth_params_module, only : QVAR, QPRES, QU, QV, QW, QRHO, NQAUX, QC,
  // QGAMC implicit none

  // integer, intent(in) :: i,j,k,idir
  // integer, intent(in) :: q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
  // integer, intent(in) :: qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3
  // amrex::Real, intent(in) :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
  // amrex::Real, intent(in) :: qaux(qa_l1:qa_h1,qa_l2:qa_h2,qa_l3:qa_h3,NQAUX)
  // amrex::Real, intent(in) :: dpdx, dudx, dvdx, dwdx, drhodx
  // amrex::Real, intent(in) :: dpdy, dudy, dvdy, dwdy, drhody
  // amrex::Real, intent(in) :: dpdz, dudz, dvdz, dwdz, drhodz
  // amrex::Real, intent(out) :: T1, T2, T3, T4, T5
  // amrex::Real :: inv_rho

  const amrex::Real inv_rho = 1.0 / q(i, j, k, QRHO);

  if (idir == 1) {
    T1 =
      (q(i, j, k, QV) * (dpdy - q(i, j, k, QRHO) * qaux(i, j, k, QC) * dudy)) +
      (q(i, j, k, QW) * (dpdz - q(i, j, k, QRHO) * qaux(i, j, k, QC) * dudz)) +
      (qaux(i, j, k, QGAMC) * q(i, j, k, QPRES) * (dvdy + dwdz));
    T2 = (q(i, j, k, QV) *
          ((qaux(i, j, k, QC) * qaux(i, j, k, QC) * drhody) - dpdy)) +
         (q(i, j, k, QW) *
          ((qaux(i, j, k, QC) * qaux(i, j, k, QC) * drhodz) - dpdz));
    T3 = q(i, j, k, QV) * dvdy + q(i, j, k, QW) * dvdz + dpdy * inv_rho;
    T4 = q(i, j, k, QV) * dwdy + q(i, j, k, QW) * dwdz + dpdz * inv_rho;
    T5 =
      (q(i, j, k, QV) * (dpdy + q(i, j, k, QRHO) * qaux(i, j, k, QC) * dudy)) +
      (q(i, j, k, QW) * (dpdz + q(i, j, k, QRHO) * qaux(i, j, k, QC) * dudz)) +
      (qaux(i, j, k, QGAMC) * q(i, j, k, QPRES) * (dvdy + dwdz));
  } else if (idir == 2) {
    T1 =
      (q(i, j, k, QU) * (dpdx - q(i, j, k, QRHO) * qaux(i, j, k, QC) * dvdx)) +
      (q(i, j, k, QW) * (dpdz - q(i, j, k, QRHO) * qaux(i, j, k, QC) * dvdz)) +
      (qaux(i, j, k, QGAMC) * q(i, j, k, QPRES) * (dudx + dwdz));
    T2 = q(i, j, k, QU) * dudx + q(i, j, k, QW) * dudz + dpdx * inv_rho;
    T3 = (q(i, j, k, QU) *
          ((qaux(i, j, k, QC) * qaux(i, j, k, QC) * drhodx) - dpdx)) +
         (q(i, j, k, QW) *
          ((qaux(i, j, k, QC) * qaux(i, j, k, QC) * drhodz) - dpdz));

    T4 = q(i, j, k, QU) * dwdx + q(i, j, k, QW) * dwdz + dpdz * inv_rho;
    T5 =
      (q(i, j, k, QU) * (dpdx + q(i, j, k, QRHO) * qaux(i, j, k, QC) * dvdx)) +
      (q(i, j, k, QW) * (dpdz + q(i, j, k, QRHO) * qaux(i, j, k, QC) * dvdz)) +
      (qaux(i, j, k, QGAMC) * q(i, j, k, QPRES) * (dudx + dwdz));
  } else if (idir == 3) {
    T1 =
      (q(i, j, k, QU) * (dpdx - q(i, j, k, QRHO) * qaux(i, j, k, QC) * dwdx)) +
      (q(i, j, k, QV) * (dpdy - q(i, j, k, QRHO) * qaux(i, j, k, QC) * dwdy)) +
      (qaux(i, j, k, QGAMC) * q(i, j, k, QPRES) * (dudx + dvdy));
    T2 = q(i, j, k, QU) * dudx + q(i, j, k, QV) * dudy + dpdx * inv_rho;
    T3 = q(i, j, k, QU) * dvdx + q(i, j, k, QV) * dvdy + dpdy * inv_rho;
    T4 = (q(i, j, k, QU) *
          ((qaux(i, j, k, QC) * qaux(i, j, k, QC) * drhodx) - dpdx)) +
         (q(i, j, k, QV) *
          ((qaux(i, j, k, QC) * qaux(i, j, k, QC) * drhody) - dpdy));
    T5 =
      (q(i, j, k, QU) * (dpdx + q(i, j, k, QRHO) * qaux(i, j, k, QC) * dwdx)) +
      (q(i, j, k, QV) * (dpdy + q(i, j, k, QRHO) * qaux(i, j, k, QC) * dwdy)) +
      (qaux(i, j, k, QGAMC) * q(i, j, k, QPRES) * (dudx + dvdy));
  } else {
    amrex::Abort("Problem of idir in impose_NSCBC_2d:compute_transverse_terms");
  }
}

void
PeleC::update_ghost_cells(
  int i,
  int j,
  int k,
  int bc_type,
  int idir,
  int isign,
  amrex::Real delta,
  int domlo[3],
  int domhi[3],
  amrex::Real L1,
  amrex::Real L2,
  amrex::Real L3,
  amrex::Real L4,
  amrex::Real L5,
  const amrex::Array4<amrex::Real>& uin,
  int uin_l1,
  int uin_l2,
  int uin_l3,
  int uin_h1,
  int uin_h2,
  int uin_h3,
  const amrex::Array4<amrex::Real>& q,
  int q_l1,
  int q_l2,
  int q_l3,
  int q_h1,
  int q_h2,
  int q_h3,
  const amrex::Array4<amrex::Real>& qaux,
  int qa_l1,
  int qa_l2,
  int qa_l3,
  int qa_h1,
  int qa_h2,
  int qa_h3)
{
  // use eos_module
  // use amrex_constants_module, only : ONE
  // use fuego_chemistry, only : NUM_SPECIES
  // use prob_params_module, only : SlipWall, NoSlipWall

  // integer, intent(in) :: i,j,k,idir,isign,bc_type
  // integer, intent(in) :: domlo(3), domhi(3)
  // integer, intent(in) :: q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
  // integer, intent(in) :: qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3
  // integer, intent(in) :: uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3
  // amrex::Real, intent(in) :: L1, L2, L3, L4, L5
  // amrex::Real, intent(in) :: delta
  // amrex::Real, intent(inout) :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
  // amrex::Real, intent(inout) ::
  // qaux(qa_l1:qa_h1,qa_l2:qa_h2,qa_l3:qa_h3,NQAUX) amrex::Real,
  // intent(inout) :: uin(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,NVAR)

  int idx_gc1, idx_gc2, idx_gc3, idx_gc4, idx_int1, idx_int2, idx_int3;
  int idx_start, idx_end, hop, n, local_index;
  amrex::Real drho, du, dv, dw, dp, wall_sign;
  const amrex::Real small = 1.e-8;

  // if ((idir == 1) || (idir == 2) || (idir == 3)) then
  //  do nothing
  // else
  //  bl_abort("Problem of idir in impose_NSCBC_3d:update_ghost_cells")
  // end if

  // if ((isign == 1) || (isign == -1)) then
  //  do nothing
  // else
  //  bl_abort("Problem of isign in impose_NSCBC_3d:update_ghost_cells")
  // end if

  // Compute new spatial derivative
  if (idir == 1) {
    local_index = i;
    drho = (L2 + 0.5 * (L1 + L5)) / (qaux(i, j, k, QC) * qaux(i, j, k, QC));
    du = (L5 - L1) / (2.0 * qaux(i, j, k, QC) * q(i, j, k, QRHO));
    dv = L3;
    dw = L4;
    dp = 0.5 * (L1 + L5);
  } else if (idir == 2) {
    local_index = j;
    drho = (L3 + 0.5 * (L1 + L5)) / (qaux(i, j, k, QC) * qaux(i, j, k, QC));
    du = L2;
    dv = (L5 - L1) / (2.0 * qaux(i, j, k, QC) * q(i, j, k, QRHO));
    dw = L4;
    dp = 0.5 * (L1 + L5);
  } else if (idir == 3) {
    local_index = k;
    drho = (L4 + 0.5 * (L1 + L5)) / (qaux(i, j, k, QC) * qaux(i, j, k, QC));
    du = L2;
    dv = L3;
    dw = (L5 - L1) / (2.0 * qaux(i, j, k, QC) * q(i, j, k, QRHO));
    dp = 0.5 * (L1 + L5);
  }

  if (isign == 1) {
    idx_gc1 = local_index - 1;
    idx_gc2 = local_index - 2;
    idx_gc3 = local_index - 3;
    idx_gc4 = local_index - 4;
    idx_int1 = local_index + 1;
    idx_int2 = local_index + 2;
    idx_int3 = local_index + 3;
    idx_start = domlo[idir] - 1;
    idx_end = domlo[idir] - 4;
  } else if (isign == -1) {
    idx_gc1 = local_index + 1;
    idx_gc2 = local_index + 2;
    idx_gc3 = local_index + 3;
    idx_gc4 = local_index + 4;
    idx_int1 = local_index - 1;
    idx_int2 = local_index - 2;
    idx_int3 = local_index - 3;
    idx_start = domhi[idir] + 1;
    idx_end = domhi[idir] + 4;
  }

  if (idir == 1) {
    // Update ghost cells
    // 2nd order
    q(idx_gc1, j, k, QU) = q(idx_int1, j, k, QU) - 2.0 * delta * du * isign;
    q(idx_gc1, j, k, QV) = q(idx_int1, j, k, QV) - 2.0 * delta * dv * isign;
    q(idx_gc1, j, k, QW) = q(idx_int1, j, k, QW) - 2.0 * delta * dw * isign;
    q(idx_gc1, j, k, QRHO) =
      q(idx_int1, j, k, QRHO) - 2.0 * delta * drho * isign;
    q(idx_gc1, j, k, QPRES) =
      q(idx_int1, j, k, QPRES) - 2.0 * delta * dp * isign;

    q(idx_gc2, j, k, QU) = -2.0 * q(idx_int1, j, k, QU) - 3.0 * q(i, j, k, QU) +
                           6.0 * q(idx_gc1, j, k, QU) +
                           6.0 * delta * du * isign;
    q(idx_gc2, j, k, QV) = -2.0 * q(idx_int1, j, k, QV) - 3.0 * q(i, j, k, QV) +
                           6.0 * q(idx_gc1, j, k, QV) +
                           6.0 * delta * dv * isign;
    q(idx_gc2, j, k, QW) = -2.0 * q(idx_int1, j, k, QW) - 3.0 * q(i, j, k, QW) +
                           6.0 * q(idx_gc1, j, k, QW) +
                           6.0 * delta * dw * isign;
    q(idx_gc2, j, k, QRHO) =
      -2.0 * q(idx_int1, j, k, QRHO) - 3.0 * q(i, j, k, QRHO) +
      6.0 * q(idx_gc1, j, k, QRHO) + 6.0 * delta * drho * isign;
    q(idx_gc2, j, k, QPRES) =
      -2.0 * q(idx_int1, j, k, QPRES) - 3.0 * q(i, j, k, QPRES) +
      6.0 * q(idx_gc1, j, k, QPRES) + 6.0 * delta * dp * isign;

    q(idx_gc3, j, k, QU) = 3.0 * q(idx_int1, j, k, QU) + 10.0 * q(i, j, k, QU) -
                           18.0 * q(idx_gc1, j, k, QU) +
                           6.0 * q(idx_gc2, j, k, QU) -
                           12.0 * delta * du * isign;
    q(idx_gc3, j, k, QV) = 3.0 * q(idx_int1, j, k, QV) + 10.0 * q(i, j, k, QV) -
                           18.0 * q(idx_gc1, j, k, QV) +
                           6.0 * q(idx_gc2, j, k, QV) -
                           12.0 * delta * dv * isign;
    q(idx_gc3, j, k, QW) = 3.0 * q(idx_int1, j, k, QW) + 10.0 * q(i, j, k, QW) -
                           18.0 * q(idx_gc1, j, k, QW) +
                           6.0 * q(idx_gc2, j, k, QW) -
                           12.0 * delta * dw * isign;
    q(idx_gc3, j, k, QRHO) =
      3.0 * q(idx_int1, j, k, QRHO) + 10.0 * q(i, j, k, QRHO) -
      18.0 * q(idx_gc1, j, k, QRHO) + 6.0 * q(idx_gc2, j, k, QRHO) -
      12.0 * delta * drho * isign;
    q(idx_gc3, j, k, QPRES) =
      3.0 * q(idx_int1, j, k, QPRES) + 10.0 * q(i, j, k, QPRES) -
      18.0 * q(idx_gc1, j, k, QPRES) + 6.0 * q(idx_gc2, j, k, QPRES) -
      12.0 * delta * dp * isign;

    q(idx_gc4, j, k, QU) =
      -2.0 * q(idx_int1, j, k, QU) - 13.0 * q(i, j, k, QU) +
      24.0 * q(idx_gc1, j, k, QU) - 12.0 * q(idx_gc2, j, k, QU) +
      4.0 * q(idx_gc3, j, k, QU) + 12.0 * delta * du * isign;
    q(idx_gc4, j, k, QV) =
      -2.0 * q(idx_int1, j, k, QV) - 13.0 * q(i, j, k, QV) +
      24.0 * q(idx_gc1, j, k, QV) - 12.0 * q(idx_gc2, j, k, QV) +
      4.0 * q(idx_gc3, j, k, QV) + 12.0 * delta * dv * isign;
    q(idx_gc4, j, k, QW) =
      -2.0 * q(idx_int1, j, k, QW) - 13.0 * q(i, j, k, QW) +
      24.0 * q(idx_gc1, j, k, QW) - 12.0 * q(idx_gc2, j, k, QW) +
      4.0 * q(idx_gc3, j, k, QW) + 12.0 * delta * dw * isign;
    q(idx_gc4, j, k, QRHO) =
      -2.0 * q(idx_int1, j, k, QRHO) - 13.0 * q(i, j, k, QRHO) +
      24.0 * q(idx_gc1, j, k, QRHO) - 12.0 * q(idx_gc2, j, k, QRHO) +
      4.0 * q(idx_gc3, j, k, QRHO) + 12.0 * delta * drho * isign;
    q(idx_gc4, j, k, QPRES) =
      -2.0 * q(idx_int1, j, k, QPRES) - 13.0 * q(i, j, k, QPRES) +
      24.0 * q(idx_gc1, j, k, QPRES) - 12.0 * q(idx_gc2, j, k, QPRES) +
      4.0 * q(idx_gc3, j, k, QPRES) + 12.0 * delta * dp * isign;

    if ((bc_type == NoSlipWall) || (bc_type == SlipWall)) {
      if (bc_type == NoSlipWall) {
        wall_sign = -1.0;
      } else if (bc_type == SlipWall) {
        wall_sign = 1.0;
      }

      q(idx_gc1, j, k, QU) = -q(i, j, k, QU);
      q(idx_gc2, j, k, QU) = -q(idx_int1, j, k, QU);
      q(idx_gc3, j, k, QU) = -q(idx_int2, j, k, QU);
      q(idx_gc4, j, k, QU) = -q(idx_int3, j, k, QU);

      q(idx_gc1, j, k, QV) = wall_sign * q(i, j, k, QV);
      q(idx_gc2, j, k, QV) = wall_sign * q(idx_int1, j, k, QV);
      q(idx_gc3, j, k, QV) = wall_sign * q(idx_int2, j, k, QV);
      q(idx_gc4, j, k, QV) = wall_sign * q(idx_int3, j, k, QV);

      q(idx_gc1, j, k, QW) = wall_sign * q(i, j, k, QW);
      q(idx_gc2, j, k, QW) = wall_sign * q(idx_int1, j, k, QW);
      q(idx_gc3, j, k, QW) = wall_sign * q(idx_int2, j, k, QW);
      q(idx_gc4, j, k, QW) = wall_sign * q(idx_int3, j, k, QW);

      q(idx_gc1, j, k, QRHO) = q(i, j, k, QRHO);
      q(idx_gc2, j, k, QRHO) = q(idx_int1, j, k, QRHO);
      q(idx_gc3, j, k, QRHO) = q(idx_int2, j, k, QRHO);
      q(idx_gc4, j, k, QRHO) = q(idx_int3, j, k, QRHO);

      q(idx_gc1, j, k, QPRES) = q(i, j, k, QPRES);
      q(idx_gc2, j, k, QPRES) = q(idx_int1, j, k, QPRES);
      q(idx_gc3, j, k, QPRES) = q(idx_int2, j, k, QPRES);
      q(idx_gc4, j, k, QPRES) = q(idx_int3, j, k, QPRES);
    }

    // Recompute missing values thanks to EOS
    for (int hop = idx_start; hop < idx_end; hop = hop - isign) {
      amrex::Real eos_state_p = q(hop, j, k, QPRES);
      amrex::Real eos_state_rho = q(hop, j, k, QRHO);
      amrex::Real eos_state_massfrac[NUM_SPECIES];
      for (int n = 0; n < NUM_SPECIES; n++) {
        eos_state_massfrac[n] = q(hop, j, k, QFS + n - 1);
      }
      amrex::Real eos_state_aux[NUM_AUX];
      for (int n = 0; n < NUM_AUX; n++) {
        eos_state_aux[n] = q(hop, j, k, QFX + n - 1);
      }

      amrex::Real eos_state_T;
      amrex::Real eos_state_e;
      // EOS::eos_rp(eos_state_p, eos_state_rho, eos_state_massfrac,
      // eos_state_T, eos_state_e);
      q(hop, j, k, QTEMP) = eos_state_T;
      q(hop, j, k, QREINT) = eos_state_e * q(hop, j, k, QRHO);
      q(hop, j, k, QGAME) = q(hop, j, k, QPRES) / q(hop, j, k, QREINT) + 1.0;

      // qaux(hop, j, k, QDPDR) = eos_state_dpdr_e;
      // qaux(hop, j, k, QDPDE) = eos_state_dpde;
      // qaux(hop, j, k, QGAMC) = eos_state_gam1;
      // qaux(hop, j, k, QC) = eos_state_cs;
      qaux(hop, j, k, QCSML) =
        amrex::max<amrex::Real>(small, small * qaux(hop, j, k, QC));

      // Here the update of the conservative variables uin seems to have only an
      // impact on the application of artificial viscosity difmag.
      uin(hop, j, k, URHO) = eos_state_rho;
      uin(hop, j, k, UMX) = q(hop, j, k, QU) * eos_state_rho;
      uin(hop, j, k, UMY) = q(hop, j, k, QV) * eos_state_rho;
      uin(hop, j, k, UMZ) = q(hop, j, k, QW) * eos_state_rho;
      uin(hop, j, k, UEINT) = eos_state_rho * eos_state_e;
      uin(hop, j, k, UEDEN) =
        eos_state_rho *
        (eos_state_e + 0.5 * (q(hop, j, k, QU) * q(hop, j, k, QU) +
                              q(hop, j, k, QV) * q(hop, j, k, QV) +
                              q(hop, j, k, QW) * q(hop, j, k, QW)));
      uin(hop, j, k, UTEMP) = eos_state_T;
      for (int n = 1; n < NUM_SPECIES; n++) {
        uin(hop, j, k, UFS + n - 1) = eos_state_rho * eos_state_massfrac[n];
      }
    }
  } else if (idir == 2) {
    // Update ghost cells
    // 2nd order
    q(i, idx_gc1, k, QU) = q(i, idx_int1, k, QU) - 2.0 * delta * du * isign;
    q(i, idx_gc1, k, QV) = q(i, idx_int1, k, QV) - 2.0 * delta * dv * isign;
    q(i, idx_gc1, k, QW) = q(i, idx_int1, k, QW) - 2.0 * delta * dw * isign;
    q(i, idx_gc1, k, QRHO) =
      q(i, idx_int1, k, QRHO) - 2.0 * delta * drho * isign;
    q(i, idx_gc1, k, QPRES) =
      q(i, idx_int1, k, QPRES) - 2.0 * delta * dp * isign;

    q(i, idx_gc2, k, QU) = -2.0 * q(i, idx_int1, k, QU) - 3.0 * q(i, j, k, QU) +
                           6.0 * q(i, idx_gc1, k, QU) +
                           6.0 * delta * du * isign;
    q(i, idx_gc2, k, QV) = -2.0 * q(i, idx_int1, k, QV) - 3.0 * q(i, j, k, QV) +
                           6.0 * q(i, idx_gc1, k, QV) +
                           6.0 * delta * dv * isign;
    q(i, idx_gc2, k, QW) = -2.0 * q(i, idx_int1, k, QW) - 3.0 * q(i, j, k, QW) +
                           6.0 * q(i, idx_gc1, k, QW) +
                           6.0 * delta * dw * isign;
    q(i, idx_gc2, k, QRHO) =
      -2.0 * q(i, idx_int1, k, QRHO) - 3.0 * q(i, j, k, QRHO) +
      6.0 * q(i, idx_gc1, k, QRHO) + 6.0 * delta * drho * isign;
    q(i, idx_gc2, k, QPRES) =
      -2.0 * q(i, idx_int1, k, QPRES) - 3.0 * q(i, j, k, QPRES) +
      6.0 * q(i, idx_gc1, k, QPRES) + 6.0 * delta * dp * isign;

    q(i, idx_gc3, k, QU) = 3.0 * q(i, idx_int1, k, QU) + 10.0 * q(i, j, k, QU) -
                           18.0 * q(i, idx_gc1, k, QU) +
                           6.0 * q(i, idx_gc2, k, QU) -
                           12.0 * delta * du * isign;
    q(i, idx_gc3, k, QV) = 3.0 * q(i, idx_int1, k, QV) + 10.0 * q(i, j, k, QV) -
                           18.0 * q(i, idx_gc1, k, QV) +
                           6.0 * q(i, idx_gc2, k, QV) -
                           12.0 * delta * dv * isign;
    q(i, idx_gc3, k, QW) = 3.0 * q(i, idx_int1, k, QW) + 10.0 * q(i, j, k, QW) -
                           18.0 * q(i, idx_gc1, k, QW) +
                           6.0 * q(i, idx_gc2, k, QW) -
                           12.0 * delta * dw * isign;
    q(i, idx_gc3, k, QRHO) =
      3.0 * q(i, idx_int1, k, QRHO) + 10.0 * q(i, j, k, QRHO) -
      18.0 * q(i, idx_gc1, k, QRHO) + 6.0 * q(i, idx_gc2, k, QRHO) -
      12.0 * delta * drho * isign;
    q(i, idx_gc3, k, QPRES) =
      3.0 * q(i, idx_int1, k, QPRES) + 10.0 * q(i, j, k, QPRES) -
      18.0 * q(i, idx_gc1, k, QPRES) + 6.0 * q(i, idx_gc2, k, QPRES) -
      12.0 * delta * dp * isign;

    q(i, idx_gc4, k, QU) =
      -2.0 * q(i, idx_int1, k, QU) - 13.0 * q(i, j, k, QU) +
      24.0 * q(i, idx_gc1, k, QU) - 12.0 * q(i, idx_gc2, k, QU) +
      4.0 * q(i, idx_gc3, k, QU) + 12.0 * delta * du * isign;
    q(i, idx_gc4, k, QV) =
      -2.0 * q(i, idx_int1, k, QV) - 13.0 * q(i, j, k, QV) +
      24.0 * q(i, idx_gc1, k, QV) - 12.0 * q(i, idx_gc2, k, QV) +
      4.0 * q(i, idx_gc3, k, QV) + 12.0 * delta * dv * isign;
    q(i, idx_gc4, k, QW) =
      -2.0 * q(i, idx_int1, k, QW) - 13.0 * q(i, j, k, QW) +
      24.0 * q(i, idx_gc1, k, QW) - 12.0 * q(i, idx_gc2, k, QW) +
      4.0 * q(i, idx_gc3, k, QW) + 12.0 * delta * dw * isign;
    q(i, idx_gc4, k, QRHO) =
      -2.0 * q(i, idx_int1, k, QRHO) - 13.0 * q(i, j, k, QRHO) +
      24.0 * q(i, idx_gc1, k, QRHO) - 12.0 * q(i, idx_gc2, k, QRHO) +
      4.0 * q(i, idx_gc3, k, QRHO) + 12.0 * delta * drho * isign;
    q(i, idx_gc4, k, QPRES) =
      -2.0 * q(i, idx_int1, k, QPRES) - 13.0 * q(i, j, k, QPRES) +
      24.0 * q(i, idx_gc1, k, QPRES) - 12.0 * q(i, idx_gc2, k, QPRES) +
      4.0 * q(i, idx_gc3, k, QPRES) + 12.0 * delta * dp * isign;

    if ((bc_type == NoSlipWall) || (bc_type == SlipWall)) {
      if (bc_type == NoSlipWall) {
        wall_sign = -1.0;
      } else if (bc_type == SlipWall) {
        wall_sign = 1.0;
      }

      q(i, idx_gc1, k, QU) = wall_sign * q(i, j, k, QU);
      q(i, idx_gc2, k, QU) = wall_sign * q(i, idx_int1, k, QU);
      q(i, idx_gc3, k, QU) = wall_sign * q(i, idx_int2, k, QU);
      q(i, idx_gc4, k, QU) = wall_sign * q(i, idx_int3, k, QU);

      q(i, idx_gc1, k, QV) = -q(i, j, k, QV);
      q(i, idx_gc2, k, QV) = -q(i, idx_int1, k, QV);
      q(i, idx_gc3, k, QV) = -q(i, idx_int2, k, QV);
      q(i, idx_gc4, k, QV) = -q(i, idx_int3, k, QV);

      q(i, idx_gc1, k, QW) = wall_sign * q(i, j, k, QW);
      q(i, idx_gc2, k, QW) = wall_sign * q(i, idx_int1, k, QW);
      q(i, idx_gc3, k, QW) = wall_sign * q(i, idx_int2, k, QW);
      q(i, idx_gc4, k, QW) = wall_sign * q(i, idx_int3, k, QW);

      q(i, idx_gc1, k, QRHO) = q(i, j, k, QRHO);
      q(i, idx_gc2, k, QRHO) = q(i, idx_int1, k, QRHO);
      q(i, idx_gc3, k, QRHO) = q(i, idx_int2, k, QRHO);
      q(i, idx_gc4, k, QRHO) = q(i, idx_int3, k, QRHO);

      q(i, idx_gc1, k, QPRES) = q(i, j, k, QPRES);
      q(i, idx_gc2, k, QPRES) = q(i, idx_int1, k, QPRES);
      q(i, idx_gc3, k, QPRES) = q(i, idx_int2, k, QPRES);
      q(i, idx_gc4, k, QPRES) = q(i, idx_int3, k, QPRES);
    }

    // Recompute missing values thanks to EOS
    for (int hop = idx_start; hop < idx_end; hop = hop - isign) {
      amrex::Real eos_state_p = q(i, hop, k, QPRES);
      amrex::Real eos_state_rho = q(i, hop, k, QRHO);
      amrex::Real eos_state_massfrac[NUM_SPECIES];
      amrex::Real eos_state_aux[NUM_AUX];
      for (int n = 0; n < NUM_SPECIES; n++) {
        eos_state_massfrac[n] = q(i, hop, k, QFS + NUM_SPECIES - 1);
      }
      for (int n = 0; n < NUM_AUX; n++) {
        eos_state_aux[n] = q(i, hop, k, QFX + NUM_AUX - 1);
      }

      amrex::Real eos_state_T;
      amrex::Real eos_state_e;
      // EOS::eos_rp(eos_state_p, eos_state_rho, eos_state_massfrac,
      // eos_state_T, eos_state_e);
      q(i, hop, k, QTEMP) = eos_state_T;
      q(i, hop, k, QREINT) = eos_state_e * q(i, hop, k, QRHO);
      q(i, hop, k, QGAME) = q(i, hop, k, QPRES) / q(i, hop, k, QREINT) + 1.0;

      // qaux(i, hop, k, QDPDR) = eos_state_dpdr_e;
      // qaux(i, hop, k, QDPDE) = eos_state_dpde;
      // qaux(i, hop, k, QGAMC) = eos_state_gam1;
      // qaux(i, hop, k, QC) = eos_state_cs;
      qaux(i, hop, k, QCSML) =
        amrex::max<amrex::Real>(small, small * qaux(i, hop, k, QC));

      // Here the update of the conservative variables uin seems to have only an
      // impact on the application of artificial viscosity difmag.
      uin(i, hop, k, URHO) = eos_state_rho;
      uin(i, hop, k, UMX) = q(i, hop, k, QU) * eos_state_rho;
      uin(i, hop, k, UMY) = q(i, hop, k, QV) * eos_state_rho;
      uin(i, hop, k, UMZ) = q(i, hop, k, QW) * eos_state_rho;
      uin(i, hop, k, UEINT) = eos_state_rho * eos_state_e;
      uin(i, hop, k, UEDEN) =
        eos_state_rho *
        (eos_state_e + 0.5 * (q(i, hop, k, QU) * q(i, hop, k, QU) +
                              q(i, hop, k, QV) * q(i, hop, k, QV) +
                              q(i, hop, k, QW) * q(i, hop, k, QW)));
      uin(i, hop, k, UTEMP) = eos_state_T;
      for (int n = 1; n < NUM_SPECIES; n++) {
        uin(i, hop, k, UFS + n - 1) = eos_state_rho * eos_state_massfrac[n];
      }
    }
  } else if (idir == 3) {
    // Update ghost cells
    // 2nd order
    q(i, j, idx_gc1, QU) = q(i, j, idx_int1, QU) - 2.0 * delta * du * isign;
    q(i, j, idx_gc1, QV) = q(i, j, idx_int1, QV) - 2.0 * delta * dv * isign;
    q(i, j, idx_gc1, QW) = q(i, j, idx_int1, QW) - 2.0 * delta * dw * isign;
    q(i, j, idx_gc1, QRHO) =
      q(i, j, idx_int1, QRHO) - 2.0 * delta * drho * isign;
    q(i, j, idx_gc1, QPRES) =
      q(i, j, idx_int1, QPRES) - 2.0 * delta * dp * isign;

    q(i, j, idx_gc2, QU) = -2.0 * q(i, j, idx_int1, QU) - 3.0 * q(i, j, k, QU) +
                           6.0 * q(i, j, idx_gc1, QU) +
                           6.0 * delta * du * isign;
    q(i, j, idx_gc2, QV) = -2.0 * q(i, j, idx_int1, QV) - 3.0 * q(i, j, k, QV) +
                           6.0 * q(i, j, idx_gc1, QV) +
                           6.0 * delta * dv * isign;
    q(i, j, idx_gc2, QW) = -2.0 * q(i, j, idx_int1, QW) - 3.0 * q(i, j, k, QW) +
                           6.0 * q(i, j, idx_gc1, QW) +
                           6.0 * delta * dw * isign;
    q(i, j, idx_gc2, QRHO) =
      -2.0 * q(i, j, idx_int1, QRHO) - 3.0 * q(i, j, k, QRHO) +
      6.0 * q(i, j, idx_gc1, QRHO) + 6.0 * delta * drho * isign;
    q(i, j, idx_gc2, QPRES) =
      -2.0 * q(i, j, idx_int1, QPRES) - 3.0 * q(i, j, k, QPRES) +
      6.0 * q(i, j, idx_gc1, QPRES) + 6.0 * delta * dp * isign;

    q(i, j, idx_gc3, QU) = 3.0 * q(i, j, idx_int1, QU) + 10.0 * q(i, j, k, QU) -
                           18.0 * q(i, j, idx_gc1, QU) +
                           6.0 * q(i, j, idx_gc2, QU) -
                           12.0 * delta * du * isign;
    q(i, j, idx_gc3, QV) = 3.0 * q(i, j, idx_int1, QV) + 10.0 * q(i, j, k, QV) -
                           18.0 * q(i, j, idx_gc1, QV) +
                           6.0 * q(i, j, idx_gc2, QV) -
                           12.0 * delta * dv * isign;
    q(i, j, idx_gc3, QW) = 3.0 * q(i, j, idx_int1, QW) + 10.0 * q(i, j, k, QW) -
                           18.0 * q(i, j, idx_gc1, QW) +
                           6.0 * q(i, j, idx_gc2, QW) -
                           12.0 * delta * dw * isign;
    q(i, j, idx_gc3, QRHO) =
      3.0 * q(i, j, idx_int1, QRHO) + 10.0 * q(i, j, k, QRHO) -
      18.0 * q(i, j, idx_gc1, QRHO) + 6.0 * q(i, j, idx_gc2, QRHO) -
      12.0 * delta * drho * isign;
    q(i, j, idx_gc3, QPRES) =
      3.0 * q(i, j, idx_int1, QPRES) + 10.0 * q(i, j, k, QPRES) -
      18.0 * q(i, j, idx_gc1, QPRES) + 6.0 * q(i, j, idx_gc2, QPRES) -
      12.0 * delta * dp * isign;

    q(i, j, idx_gc4, QU) =
      -2.0 * q(i, j, idx_int1, QU) - 13.0 * q(i, j, k, QU) +
      24.0 * q(i, j, idx_gc1, QU) - 12.0 * q(i, j, idx_gc2, QU) +
      4.0 * q(i, j, idx_gc3, QU) + 12.0 * delta * du * isign;
    q(i, j, idx_gc4, QV) =
      -2.0 * q(i, j, idx_int1, QV) - 13.0 * q(i, j, k, QV) +
      24.0 * q(i, j, idx_gc1, QV) - 12.0 * q(i, j, idx_gc2, QV) +
      4.0 * q(i, j, idx_gc3, QV) + 12.0 * delta * dv * isign;
    q(i, j, idx_gc4, QW) =
      -2.0 * q(i, j, idx_int1, QW) - 13.0 * q(i, j, k, QW) +
      24.0 * q(i, j, idx_gc1, QW) - 12.0 * q(i, j, idx_gc2, QW) +
      4.0 * q(i, j, idx_gc3, QW) + 12.0 * delta * dw * isign;
    q(i, j, idx_gc4, QRHO) =
      -2.0 * q(i, j, idx_int1, QRHO) - 13.0 * q(i, j, k, QRHO) +
      24.0 * q(i, j, idx_gc1, QRHO) - 12.0 * q(i, j, idx_gc2, QRHO) +
      4.0 * q(i, j, idx_gc3, QRHO) + 12.0 * delta * drho * isign;
    q(i, j, idx_gc4, QPRES) =
      -2.0 * q(i, j, idx_int1, QPRES) - 13.0 * q(i, j, k, QPRES) +
      24.0 * q(i, j, idx_gc1, QPRES) - 12.0 * q(i, j, idx_gc2, QPRES) +
      4.0 * q(i, j, idx_gc3, QPRES) + 12.0 * delta * dp * isign;

    if ((bc_type == NoSlipWall) || (bc_type == SlipWall)) {
      if (bc_type == NoSlipWall) {
        wall_sign = -1.0;
      } else if (bc_type == SlipWall) {
        wall_sign = 1.0;
      }

      q(i, j, idx_gc1, QU) = wall_sign * q(i, j, k, QU);
      q(i, j, idx_gc2, QU) = wall_sign * q(i, j, idx_int1, QU);
      q(i, j, idx_gc3, QU) = wall_sign * q(i, j, idx_int2, QU);
      q(i, j, idx_gc4, QU) = wall_sign * q(i, j, idx_int3, QU);

      q(i, j, idx_gc1, QV) = wall_sign * q(i, j, k, QV);
      q(i, j, idx_gc2, QV) = wall_sign * q(i, j, idx_int1, QV);
      q(i, j, idx_gc3, QV) = wall_sign * q(i, j, idx_int2, QV);
      q(i, j, idx_gc4, QV) = wall_sign * q(i, j, idx_int3, QV);

      q(i, j, idx_gc1, QW) = -q(i, j, k, QW);
      q(i, j, idx_gc2, QW) = -q(i, j, idx_int1, QW);
      q(i, j, idx_gc3, QW) = -q(i, j, idx_int2, QW);
      q(i, j, idx_gc4, QW) = -q(i, j, idx_int3, QW);

      q(i, j, idx_gc1, QRHO) = q(i, j, k, QRHO);
      q(i, j, idx_gc2, QRHO) = q(i, j, idx_int1, QRHO);
      q(i, j, idx_gc3, QRHO) = q(i, j, idx_int2, QRHO);
      q(i, j, idx_gc4, QRHO) = q(i, j, idx_int3, QRHO);

      q(i, j, idx_gc1, QPRES) = q(i, j, k, QPRES);
      q(i, j, idx_gc2, QPRES) = q(i, j, idx_int1, QPRES);
      q(i, j, idx_gc3, QPRES) = q(i, j, idx_int2, QPRES);
      q(i, j, idx_gc4, QPRES) = q(i, j, idx_int3, QPRES);
    }

    // Recompute missing values thanks to EOS
    for (int hop = idx_start; hop < idx_end; hop = hop - isign) {
      amrex::Real eos_state_p = q(i, j, hop, QPRES);
      amrex::Real eos_state_rho = q(i, j, hop, QRHO);
      amrex::Real eos_state_massfrac[NUM_SPECIES];
      for (int n = 0; n < NUM_SPECIES; n++) {
        eos_state_massfrac[n] = q(i, j, hop, QFS + NUM_SPECIES - 1);
      }
      amrex::Real eos_state_aux[NUM_AUX];
      for (int n = 0; n < NUM_SPECIES; n++) {
        eos_state_aux[n] = q(i, j, hop, QFX + NUM_AUX - 1);
      }

      amrex::Real eos_state_T;
      amrex::Real eos_state_e;
      // eos_rp(eos_state);
      q(i, j, hop, QTEMP) = eos_state_T;
      q(i, j, hop, QREINT) = eos_state_e * q(i, j, hop, QRHO);
      q(i, j, hop, QGAME) = q(i, j, hop, QPRES) / q(i, j, hop, QREINT) + 1.0;

      // qaux(i, j, hop, QDPDR) = eos_state_dpdr_e;
      // qaux(i, j, hop, QDPDE) = eos_state_dpde;
      // qaux(i, j, hop, QGAMC) = eos_state_gam1;
      // qaux(i, j, hop, QC) = eos_state_cs;
      qaux(i, j, hop, QCSML) =
        amrex::max<amrex::Real>(small, small * qaux(i, j, hop, QC));

      // Here the update of the conservative variables uin seems to have only an
      // impact on the application of artificial viscosity difmag.
      uin(i, j, hop, URHO) = eos_state_rho;
      uin(i, j, hop, UMX) = q(i, j, hop, QU) * eos_state_rho;
      uin(i, j, hop, UMY) = q(i, j, hop, QV) * eos_state_rho;
      uin(i, j, hop, UMZ) = q(i, j, hop, QW) * eos_state_rho;
      uin(i, j, hop, UEINT) = eos_state_rho * eos_state_e;
      uin(i, j, hop, UEDEN) =
        eos_state_rho *
        (eos_state_e + 0.5 * (q(i, j, hop, QU) * q(i, j, hop, QU) +
                              q(i, j, hop, QV) * q(i, j, hop, QV) +
                              q(i, j, hop, QW) * q(i, j, hop, QW)));
      uin(i, j, hop, UTEMP) = eos_state_T;
      for (int n = 1; n < NUM_SPECIES; n++) {
        uin(i, j, hop, UFS + n - 1) = eos_state_rho * eos_state_massfrac[n];
      }
    }
  }
  // destroy(eos_state);
}

void
PeleC::compute_waves(
  int i,
  int j,
  int k,
  int idir,
  int isign,
  int bc_type,
  amrex::Real bc_params[6],
  amrex::Real bc_target[5],
  amrex::Real T1,
  amrex::Real T2,
  amrex::Real T3,
  amrex::Real T4,
  amrex::Real T5,
  amrex::Real L1,
  amrex::Real L2,
  amrex::Real L3,
  amrex::Real L4,
  amrex::Real L5,
  amrex::Real dp,
  amrex::Real du,
  amrex::Real dv,
  amrex::Real dw,
  amrex::Real drho,
  const amrex::Array4<amrex::Real>& q,
  int q_l1,
  int q_l2,
  int q_l3,
  int q_h1,
  int q_h2,
  int q_h3,
  const amrex::Array4<amrex::Real>& qaux,
  int qa_l1,
  int qa_l2,
  int qa_l3,
  int qa_h1,
  int qa_h2,
  int qa_h3)
{
  // use prob_params_module, only : probhi, Inflow, Outflow, SlipWall,
  // NoSlipWall

  // Note: for convenience, all waves are called L because we have just 1
  // direction here Waves M and N will be employed for corners

  // Local
  amrex::Real mach_local, TARGET_VX, TARGET_VY, TARGET_VZ, TARGET_TEMPERATURE,
    TARGET_PRESSURE;
  amrex::Real beta, relax_T, relax_U, relax_V, relax_W, sigma_out, Kout;

  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> probhi =
    geom.ProbHiArray();

  // if ((idir == 1) || (idir == 2) || (idir == 3)) {
  //  do nothing
  //} else {
  //  amrex::Abort("Problem of idir in impose_NSCBC_3d:compute_waves");
  //}

  // if ((isign == 1) || (isign == -1)) {
  //  do nothing
  //} else {
  //  amrex::Abort("Problem of isign in impose_NSCBC_3d:compute_waves");
  //}

  mach_local =
    std::sqrt(
      q(i, j, k, QU) * q(i, j, k, QU) + q(i, j, k, QV) * q(i, j, k, QV) +
      q(i, j, k, QW) * q(i, j, k, QW)) /
    qaux(i, j, k, QC);

  // Recasting targets values and numerical parameters
  TARGET_VX = bc_target[0];
  TARGET_VY = bc_target[1];
  TARGET_VZ = bc_target[2];
  TARGET_TEMPERATURE = bc_target[3];
  TARGET_PRESSURE = bc_target[4];

  relax_T = bc_params[0];
  relax_U = bc_params[1];
  relax_V = bc_params[2];
  relax_W = bc_params[3];
  // Here we have the abilities to set beta=local Mach number
  // it may works better for outflow BCs
  if (bc_params[4] < 0.0) {
    beta = mach_local;
  } else {
    beta = bc_params[4];
  }
  sigma_out = bc_params[5];

  // Computing known numerical LODI waves
  if (idir == 1) {
    // Numerical LODI waves along X
    L1 = (q(i, j, k, QU) - qaux(i, j, k, QC)) *
         (dp - (q(i, j, k, QRHO) * qaux(i, j, k, QC)) * du);
    L2 =
      q(i, j, k, QU) * (((qaux(i, j, k, QC) * qaux(i, j, k, QC)) * drho) - dp);
    L3 = q(i, j, k, QU) * dv;
    L4 = q(i, j, k, QU) * dw;
    L5 = (q(i, j, k, QU) + qaux(i, j, k, QC)) *
         (dp + (q(i, j, k, QRHO) * qaux(i, j, k, QC)) * du);
  } else if (idir == 2) {
    // Numerical LODI waves along Y
    L1 = (q(i, j, k, QV) - qaux(i, j, k, QC)) *
         (dp - (q(i, j, k, QRHO) * qaux(i, j, k, QC)) * dv);
    L2 = q(i, j, k, QV) * du;
    L3 =
      q(i, j, k, QV) * (((qaux(i, j, k, QC) * qaux(i, j, k, QC)) * drho) - dp);
    L4 = q(i, j, k, QV) * dw;
    L5 = (q(i, j, k, QV) + qaux(i, j, k, QC)) *
         (dp + (q(i, j, k, QRHO) * qaux(i, j, k, QC)) * dv);
  } else if (idir == 3) {
    // Numerical LODI waves along Z
    L1 = (q(i, j, k, QW) - qaux(i, j, k, QC)) *
         (dp - (q(i, j, k, QRHO) * qaux(i, j, k, QC)) * dw);
    L2 = q(i, j, k, QW) * du;
    L3 = q(i, j, k, QW) * dv;
    L4 =
      q(i, j, k, QW) * (((qaux(i, j, k, QC) * qaux(i, j, k, QC)) * drho) - dp);
    L5 = (q(i, j, k, QW) + qaux(i, j, k, QC)) *
         (dp + (q(i, j, k, QRHO) * qaux(i, j, k, QC)) * dw);
  }

  // Computing missing LODI waves from BC model
  if (bc_type == Inflow) {
    if (idir == 1) {
      if (isign == 1) {
        L5 = relax_U *
               ((q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QC)) *
                (1.0 - mach_local * mach_local) / probhi[idir]) *
               (q(i, j, k, QU) - TARGET_VX) -
             ((1.0 - beta) * T5);
      } else if (isign == -1) {
        L1 = relax_U *
               ((q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QC)) *
                (1.0 - mach_local * mach_local) / probhi[idir]) *
               (q(i, j, k, QU) - TARGET_VX) -
             ((1.0 - beta) * T1);
      }
      L2 = relax_T *
             (q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QRSPEC) /
              probhi[idir]) *
             (q(i, j, k, QTEMP) - TARGET_TEMPERATURE) -
           ((1.0 - beta) * T2);
      L3 = relax_V * (qaux(i, j, k, QC) / probhi[idir]) *
             (q(i, j, k, QV) - TARGET_VY) -
           ((1.0 - beta) * T3);
      L4 = relax_W * (qaux(i, j, k, QC) / probhi[idir]) *
             (q(i, j, k, QW) - TARGET_VZ) -
           ((1.0 - beta) * T4);
    } else if (idir == 2) {
      if (isign == 1) {
        L5 = relax_V *
               ((q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QC)) *
                (1.0 - mach_local * mach_local) / probhi[idir]) *
               (q(i, j, k, QV) - TARGET_VY) -
             ((1.0 - beta) * T5);
      } else if (isign == -1) {
        L1 = relax_V *
               ((q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QC)) *
                (1.0 - mach_local * mach_local) / probhi[idir]) *
               (q(i, j, k, QV) - TARGET_VY) -
             ((1.0 - beta) * T1);
      }
      L2 = relax_U * (qaux(i, j, k, QC) / probhi[idir]) *
             (q(i, j, k, QU) - TARGET_VX) -
           ((1.0 - beta) * T2);
      L3 = relax_T *
             (q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QRSPEC) /
              probhi[idir]) *
             (q(i, j, k, QTEMP) - TARGET_TEMPERATURE) -
           ((1.0 - beta) * T3);
      L4 = relax_W * (qaux(i, j, k, QC) / probhi[idir]) *
             (q(i, j, k, QW) - TARGET_VZ) -
           ((1.0 - beta) * T4);
    } else if (idir == 3) {
      if (isign == 1) {
        L5 = relax_W *
               ((q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QC)) *
                (1.0 - mach_local * mach_local) / probhi[idir]) *
               (q(i, j, k, QW) - TARGET_VZ) -
             ((1.0 - beta) * T5);
      } else if (isign == -1) {
        L1 = relax_W *
               ((q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QC)) *
                (1.0 - mach_local * mach_local) / probhi[idir]) *
               (q(i, j, k, QW) - TARGET_VZ) -
             ((1.0 - beta) * T1);
      }
      L2 = relax_W * (qaux(i, j, k, QC) / probhi[idir]) *
             (q(i, j, k, QU) - TARGET_VX) -
           ((1.0 - beta) * T2);
      L3 = relax_W * (qaux(i, j, k, QC) / probhi[idir]) *
             (q(i, j, k, QV) - TARGET_VY) -
           ((1.0 - beta) * T3);
      L4 = relax_T *
             (q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QRSPEC) /
              probhi[idir]) *
             (q(i, j, k, QTEMP) - TARGET_TEMPERATURE) -
           ((1.0 - beta) * T4);
    } else {
      amrex::Abort("Error:: Wait, is this the fourth dimension?");
    }
  } else if ((bc_type == SlipWall) || (bc_type == NoSlipWall)) {
    // Values long Y will be computed by mirror functions below
    // but we set waves values to 0 to avoid undefined variables
    L1 = 0.0;
    L2 = 0.0;
    L3 = 0.0;
    L4 = 0.0;
    L5 = 0.0;
  } else if (bc_type == Outflow) {
    Kout = sigma_out * (1.0 - (mach_local * mach_local)) *
           (qaux(i, j, k, QC) / probhi[idir]);
    if (isign == 1) {
      L5 = (Kout * (q(i, j, k, QPRES) - TARGET_PRESSURE)) - ((1.0 - beta) * T5);
    } else if (isign == -1) {
      L1 = (Kout * (q(i, j, k, QPRES) - TARGET_PRESSURE)) - ((1.0 - beta) * T1);
    }
  } else {
    amrex::Abort("Error:: This BC is not yet implemented for x dir in "
                 "characteristic form");
  }

  // Shaping the waves to be at the good dimension
  if (idir == 1) {
    L1 = L1 / (q(i, j, k, QU) - qaux(i, j, k, QC));
    L5 = L5 / (q(i, j, k, QU) + qaux(i, j, k, QC));
    if (q(i, j, k, QU) == 0.0) {
      L2 = 0.0;
      L3 = 0.0;
      L4 = 0.0;
    } else {
      L2 = L2 / q(i, j, k, QU);
      L3 = L3 / q(i, j, k, QU);
      L4 = L4 / q(i, j, k, QU);
    }
  } else if (idir == 2) {
    L1 = L1 / (q(i, j, k, QV) - qaux(i, j, k, QC));
    L5 = L5 / (q(i, j, k, QV) + qaux(i, j, k, QC));
    if (q(i, j, k, QV) == 0.0) {
      L2 = 0.0;
      L3 = 0.0;
      L4 = 0.0;
    } else {
      L2 = L2 / q(i, j, k, QV);
      L3 = L3 / q(i, j, k, QV);
      L4 = L4 / q(i, j, k, QV);
    }
  } else if (idir == 3) {
    L1 = L1 / (q(i, j, k, QW) - qaux(i, j, k, QC));
    L5 = L5 / (q(i, j, k, QW) + qaux(i, j, k, QC));
    if (q(i, j, k, QW) == 0.0) {
      L2 = 0.0;
      L3 = 0.0;
      L4 = 0.0;
    } else {
      L2 = L2 / q(i, j, k, QW);
      L3 = L3 / q(i, j, k, QW);
      L4 = L4 / q(i, j, k, QW);
    }
  }
}

void
PeleC::set_bc_mask(
  const int lo[3],
  const int hi[3],
  const int domlo[3],
  const int domhi[3],
  AMREX_D_DECL(
    const amrex::Array4<int>& x_bcMask,
    const amrex::Array4<int>& y_bcMask,
    const amrex::Array4<int>& z_bcMask))
{
  const int x_bcMask_l1 = lo[0];
  const int x_bcMask_l2 = lo[1];
  const int x_bcMask_l3 = lo[2];
  const int x_bcMask_h1 = hi[0];
  const int x_bcMask_h2 = hi[1];
  const int x_bcMask_h3 = hi[2];
  const int y_bcMask_l1 = lo[0];
  const int y_bcMask_l2 = lo[1];
  const int y_bcMask_l3 = lo[2];
  const int y_bcMask_h1 = hi[0];
  const int y_bcMask_h2 = hi[1];
  const int y_bcMask_h3 = hi[2];
  const int z_bcMask_l1 = lo[0];
  const int z_bcMask_l2 = lo[1];
  const int z_bcMask_l3 = lo[2];
  const int z_bcMask_h1 = hi[0];
  const int z_bcMask_h2 = hi[1];
  const int z_bcMask_h3 = hi[2];

  if (x_bcMask_l1 == domlo[0]) {
    for (int i = x_bcMask_l2; i < x_bcMask_h2; i++) {
      for (int j = x_bcMask_l3; j < x_bcMask_h3; j++) {
        x_bcMask(domlo[0], i, j) = PeleC::phys_bc.lo()[0];
      }
    }
  }

  if (x_bcMask_h1 == domhi[0] + 1) {
    for (int i = x_bcMask_l2; i < x_bcMask_h2; i++) {
      for (int j = x_bcMask_l3; j < x_bcMask_h3; j++) {
        x_bcMask(domhi[0] + 1, i, j) = PeleC::phys_bc.hi()[0];
      }
    }
  }

  if (y_bcMask_l2 == domlo[1]) {
    for (int i = y_bcMask_l1; i < y_bcMask_h1; i++) {
      for (int j = y_bcMask_l3; j < y_bcMask_h3; j++) {
        y_bcMask(i, domlo[1], j) = PeleC::phys_bc.lo()[1];
      }
    }
  }

  if (y_bcMask_h2 == domhi[1] + 1) {
    for (int i = y_bcMask_l1; i < y_bcMask_h1; i++) {
      for (int j = y_bcMask_l3; j < y_bcMask_h3; j++) {
        y_bcMask(i, domhi[1] + 1, j) = PeleC::phys_bc.hi()[1];
      }
    }
  }

  if (z_bcMask_l3 == domlo[2]) {
    for (int i = z_bcMask_l1; i < z_bcMask_h1; i++) {
      for (int j = z_bcMask_l2; j < z_bcMask_h2; j++) {
        z_bcMask(i, j, domlo[2]) = PeleC::phys_bc.lo()[2];
      }
    }
  }

  if (z_bcMask_h3 == domhi[2] + 1) {
    for (int i = z_bcMask_l1; i < z_bcMask_h1; i++) {
      for (int j = z_bcMask_l2; j < z_bcMask_h2; j++) {
        z_bcMask(i, j, domhi[2] + 1) = PeleC::phys_bc.hi()[2];
      }
    }
  }
}
