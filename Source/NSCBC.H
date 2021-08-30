#ifndef _NSCBC_H_
#define _NSCBC_H_

#include "PeleC.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void normal_derivative(
  int i,
  int j,
  int k,
  int idir,
  int isign,
  amrex::Real delta,
  amrex::Real& dp,
  amrex::Real& du,
  amrex::Real& dv,
  amrex::Real& dw,
  amrex::Real& drho,
  const amrex::Array4<amrex::Real>& q)
{
  if (idir == 0) {
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
  } else if (idir == 1) {
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
  } else if (idir == 2) {
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


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void compute_transverse_terms(
  int i,
  int j,
  int k,
  int idir,
  amrex::Real* T,
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
  const amrex::Array4<amrex::Real>& qaux)
{
  const amrex::Real inv_rho = 1.0 / q(i, j, k, QRHO);
  amrex::Real T1, T2, T3, T4, T5;
  if (idir == 0) {
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
  } else if (idir == 1) {
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
  } else if (idir == 2) {
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
  T[0] = T1;
  T[1] = T2;
  T[2] = T3;
  T[3] = T4;
  T[5] = T5;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void tangential_derivative(
  int i,
  int j,
  int k,
  int idir,
  amrex::Real delta,
  amrex::Real& dp,
  amrex::Real& du,
  amrex::Real& dv,
  amrex::Real& dw,
  amrex::Real& drho,
  const amrex::Array4<amrex::Real>& q)
{
  // Warning, idir means the tangential direction, this is different from 2D
  // (sorry)
  if (idir == 0) {
    // 2nd order Central
    dp = (q(i + 1, j, k, QPRES) - q(i - 1, j, k, QPRES)) / (2.0 * delta);
    du = (q(i + 1, j, k, QU) - q(i - 1, j, k, QU)) / (2.0 * delta);
    dv = (q(i + 1, j, k, QV) - q(i - 1, j, k, QV)) / (2.0 * delta);
    dw = (q(i + 1, j, k, QW) - q(i - 1, j, k, QW)) / (2.0 * delta);
    drho = (q(i + 1, j, k, QRHO) - q(i - 1, j, k, QRHO)) / (2.0 * delta);
  } else if (idir == 1) {
    // 2nd order Central
    dp = (q(i, j + 1, k, QPRES) - q(i, j - 1, k, QPRES)) / (2.0 * delta);
    du = (q(i, j + 1, k, QU) - q(i, j - 1, k, QU)) / (2.0 * delta);
    dv = (q(i, j + 1, k, QV) - q(i, j - 1, k, QV)) / (2.0 * delta);
    dw = (q(i, j + 1, k, QW) - q(i, j - 1, k, QW)) / (2.0 * delta);
    drho = (q(i, j + 1, k, QRHO) - q(i, j - 1, k, QRHO)) / (2.0 * delta);
  } else if (idir == 2) {
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

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void update_ghost_cells(
  int i,
  int j,
  int k,
  int bc_type,
  int idir,
  int isign,
  amrex::Real delta,
  const int* domlo,
  const int* domhi,
  const amrex::Real* L,
  const amrex::Array4<amrex::Real>& uin,
  const amrex::Array4<amrex::Real>& q,
  const amrex::Array4<amrex::Real>& qaux)
{
  int idx_gc1, idx_gc2, idx_gc3, idx_gc4, idx_int1, idx_int2, idx_int3;
  int idx_start, idx_end, hop, n, local_index;
  amrex::Real drho, du, dv, dw, dp, wall_sign;
  const amrex::Real small = 1.e-8;
  AMREX_ASSERT((idir == 0) || (idir == 1) || (idir == 2));
  AMREX_ASSERT((isign == 1) || (isign == -1));
  amrex::Real L1 = L[0];
  amrex::Real L2 = L[1];
  amrex::Real L3 = L[2];
  amrex::Real L4 = L[3];
  amrex::Real L5 = L[4];

  // Compute new spatial derivative
  if (idir == 0) {
    local_index = i;
    drho = (L2 + 0.5 * (L1 + L5)) / (qaux(i, j, k, QC) * qaux(i, j, k, QC));
    du = (L5 - L1) / (2.0 * qaux(i, j, k, QC) * q(i, j, k, QRHO));
    dv = L3;
    dw = L4;
    dp = 0.5 * (L1 + L5);
  } else if (idir == 1) {
    local_index = j;
    drho = (L3 + 0.5 * (L1 + L5)) / (qaux(i, j, k, QC) * qaux(i, j, k, QC));
    du = L2;
    dv = (L5 - L1) / (2.0 * qaux(i, j, k, QC) * q(i, j, k, QRHO));
    dw = L4;
    dp = 0.5 * (L1 + L5);
  } else if (idir == 2) {
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

  if (idir == 0) {
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

    // If BC is SlipWall or NoSlipWall
    if ((bc_type == 4) || (bc_type == 5)) {
      if (bc_type == 5) {
        wall_sign = -1.0;
      } else if (bc_type == 4) {
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


  } else if (idir == 1) {
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

    if ((bc_type == 4) || (bc_type == 5)) {
      if (bc_type == 5) {
        wall_sign = -1.0;
      } else if (bc_type == 4) {
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
  } else if (idir == 2) {
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

    if ((bc_type == 4) || (bc_type == 5)) {
      if (bc_type == 5) {
        wall_sign = -1.0;
      } else if (bc_type == 4) {
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
  }
  auto eos = pele::physics::PhysicsType::eos();
  // Recompute missing values thanks to EOS
  int ii = i;
  int jj = j;
  int kk = k;
  for (int hop = idx_start; hop < idx_end; hop = hop - isign) {
    if (idir == 0) {
      ii = hop;
    } else if (idir == 1) {
      jj = hop;
    } else if (idir == 2) {
      kk = hop;
    }
    amrex::Real eos_state_p = q(ii, jj, kk, QPRES);
    amrex::Real eos_state_rho = q(ii, jj, kk, QRHO);
    amrex::Real eos_state_massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; n++) {
      eos_state_massfrac[n] = q(ii, jj, kk, QFS + n - 1);
    }
    amrex::Real eos_state_aux[NUM_AUX];
    for (int n = 0; n < NUM_AUX; n++) {
      eos_state_aux[n] = q(ii, jj, kk, QFX + n - 1);
    }

    amrex::Real eos_state_T, eos_state_e, eos_state_gamma, eos_state_cs;
    amrex::Real eos_state_dpdr_e, eos_state_dpde, eos_state_gam1, wbar;
    eos.Y2WBAR(eos_state_massfrac, wbar);
    eos.RYP2E(eos_state_rho, eos_state_massfrac, eos_state_p, eos_state_e);
    eos.EY2T(eos_state_e, eos_state_massfrac, eos_state_T);
    eos.TY2G(eos_state_T, eos_state_massfrac, eos_state_gamma);
    eos.RTY2Cs(eos_state_rho, eos_state_T, eos_state_massfrac, eos_state_cs);
    eos.RTY2dpde_dpdre(eos_state_rho, eos_state_T, eos_state_massfrac,
      eos_state_dpde, eos_state_dpdr_e);
    q(ii, jj, kk, QTEMP) = eos_state_T;
    q(ii, jj, kk, QREINT) = eos_state_e * q(ii, jj, kk, QRHO);
    q(ii, jj, kk, QGAME) = q(ii, jj, kk, QPRES) / q(ii, jj, kk, QREINT) + 1.0;

    qaux(ii, jj, kk, QDPDR) = eos_state_dpdr_e;
    qaux(ii, jj, kk, QDPDE) = eos_state_dpde;
    qaux(ii, jj, kk, QGAMC) = eos_state_gam1;
    qaux(ii, jj, kk, QC) = eos_state_cs;
    qaux(ii, jj, kk, QCSML) = amrex::max<amrex::Real>(
      constants::small_num(), constants::small_num() * eos_state_cs);
    qaux(ii, jj, kk, QRSPEC) = pele::physics::Constants::RU / wbar;

    // Here the update of the conservative variables uin seems to have only an
    // impact on the application of artificial viscosity difmag.
    uin(ii, jj, kk, URHO) = eos_state_rho;
    uin(ii, jj, kk, UMX) = q(ii, jj, kk, QU) * eos_state_rho;
    uin(ii, jj, kk, UMY) = q(ii, jj, kk, QV) * eos_state_rho;
    uin(ii, jj, kk, UMZ) = q(ii, jj, kk, QW) * eos_state_rho;
    uin(ii, jj, kk, UEINT) = eos_state_rho * eos_state_e;
    uin(ii, jj, kk, UEDEN) =
      eos_state_rho *
      (eos_state_e + 0.5 * (q(ii, jj, kk, QU) * q(ii, jj, kk, QU) +
                            q(ii, jj, kk, QV) * q(ii, jj, kk, QV) +
                            q(ii, jj, kk, QW) * q(ii, jj, kk, QW)));
    uin(ii, jj, kk, UTEMP) = eos_state_T;
    for (int n = 1; n < NUM_SPECIES; n++) {
      uin(ii, jj, kk, UFS + n - 1) = eos_state_rho * eos_state_massfrac[n];
    }
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
compute_waves(
  int i,
  int j,
  int k,
  int idir,
  int isign,
  int bc_type,
  const amrex::Real* problen,
  const amrex::Real bc_params[6],
  const amrex::Real bc_target[5],
  const amrex::Real T[5],
  amrex::Real* L,
  amrex::Real dp,
  amrex::Real du,
  amrex::Real dv,
  amrex::Real dw,
  amrex::Real drho,
  const amrex::Array4<amrex::Real>& q,
  const amrex::Array4<amrex::Real>& qaux)
{
  // Note: for convenience, all waves are called L because we have just 1
  // direction here Waves M and N will be employed for corners

  AMREX_ASSERT((idir == 0) || (idir == 1) || (idir == 2));

  amrex::Real mach_local =
    std::sqrt(
      q(i, j, k, QU) * q(i, j, k, QU) + q(i, j, k, QV) * q(i, j, k, QV) +
      q(i, j, k, QW) * q(i, j, k, QW)) /
    qaux(i, j, k, QC);

  // Recasting targets values and numerical parameters
  amrex::Real TARGET_VX = bc_target[0];
  amrex::Real TARGET_VY = bc_target[1];
  amrex::Real TARGET_VZ = bc_target[2];
  amrex::Real TARGET_TEMPERATURE = bc_target[3];
  amrex::Real TARGET_PRESSURE = bc_target[4];

  amrex::Real relax_T = bc_params[0];
  amrex::Real relax_U = bc_params[1];
  amrex::Real relax_V = bc_params[2];
  amrex::Real relax_W = bc_params[3];
  // Here we have the abilities to set beta=local Mach number
  // it may works better for outflow BCs
  amrex::Real beta = bc_params[4];
  if (bc_params[4] < 0.0) {
    beta = mach_local;
  }
  amrex::Real sigma_out = bc_params[5];

  amrex::Real L1, L2, L3, L4, L5, Kout;
  const amrex::Real T1 = T[0];
  const amrex::Real T2 = T[1];
  const amrex::Real T3 = T[2];
  const amrex::Real T4 = T[3];
  const amrex::Real T5 = T[4];
  // Computing known numerical LODI waves
  if (idir == 0) {
    // Numerical LODI waves along X
    L1 = (q(i, j, k, QU) - qaux(i, j, k, QC)) *
         (dp - (q(i, j, k, QRHO) * qaux(i, j, k, QC)) * du);
    L2 =
      q(i, j, k, QU) * (((qaux(i, j, k, QC) * qaux(i, j, k, QC)) * drho) - dp);
    L3 = q(i, j, k, QU) * dv;
    L4 = q(i, j, k, QU) * dw;
    L5 = (q(i, j, k, QU) + qaux(i, j, k, QC)) *
         (dp + (q(i, j, k, QRHO) * qaux(i, j, k, QC)) * du);
  } else if (idir == 1) {
    // Numerical LODI waves along Y
    L1 = (q(i, j, k, QV) - qaux(i, j, k, QC)) *
         (dp - (q(i, j, k, QRHO) * qaux(i, j, k, QC)) * dv);
    L2 = q(i, j, k, QV) * du;
    L3 =
      q(i, j, k, QV) * (((qaux(i, j, k, QC) * qaux(i, j, k, QC)) * drho) - dp);
    L4 = q(i, j, k, QV) * dw;
    L5 = (q(i, j, k, QV) + qaux(i, j, k, QC)) *
         (dp + (q(i, j, k, QRHO) * qaux(i, j, k, QC)) * dv);
  } else if (idir == 2) {
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
  if (bc_type == 7) {
    // Inflow
    if (idir == 0) {
      if (isign == 1) {
        L5 = relax_U *
               ((q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QC)) *
                (1.0 - mach_local * mach_local) / problen[idir]) *
               (q(i, j, k, QU) - TARGET_VX) -
             ((1.0 - beta) * T5);
      } else if (isign == -1) {
        L1 = relax_U *
               ((q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QC)) *
                (1.0 - mach_local * mach_local) / problen[idir]) *
               (q(i, j, k, QU) - TARGET_VX) -
             ((1.0 - beta) * T1);
      }
      L2 = relax_T *
             (q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QRSPEC) /
              problen[idir]) *
             (q(i, j, k, QTEMP) - TARGET_TEMPERATURE) -
           ((1.0 - beta) * T2);
      L3 = relax_V * (qaux(i, j, k, QC) / problen[idir]) *
             (q(i, j, k, QV) - TARGET_VY) -
           ((1.0 - beta) * T3);
      L4 = relax_W * (qaux(i, j, k, QC) / problen[idir]) *
             (q(i, j, k, QW) - TARGET_VZ) -
           ((1.0 - beta) * T4);
    } else if (idir == 1) {
      if (isign == 1) {
        L5 = relax_V *
               ((q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QC)) *
                (1.0 - mach_local * mach_local) / problen[idir]) *
               (q(i, j, k, QV) - TARGET_VY) -
             ((1.0 - beta) * T5);
      } else if (isign == -1) {
        L1 = relax_V *
               ((q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QC)) *
                (1.0 - mach_local * mach_local) / problen[idir]) *
               (q(i, j, k, QV) - TARGET_VY) -
             ((1.0 - beta) * T1);
      }
      L2 = relax_U * (qaux(i, j, k, QC) / problen[idir]) *
             (q(i, j, k, QU) - TARGET_VX) -
           ((1.0 - beta) * T2);
      L3 = relax_T *
             (q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QRSPEC) /
              problen[idir]) *
             (q(i, j, k, QTEMP) - TARGET_TEMPERATURE) -
           ((1.0 - beta) * T3);
      L4 = relax_W * (qaux(i, j, k, QC) / problen[idir]) *
             (q(i, j, k, QW) - TARGET_VZ) -
           ((1.0 - beta) * T4);
    } else if (idir == 2) {
      if (isign == 1) {
        L5 = relax_W *
               ((q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QC)) *
                (1.0 - mach_local * mach_local) / problen[idir]) *
               (q(i, j, k, QW) - TARGET_VZ) -
             ((1.0 - beta) * T5);
      } else if (isign == -1) {
        L1 = relax_W *
               ((q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QC)) *
                (1.0 - mach_local * mach_local) / problen[idir]) *
               (q(i, j, k, QW) - TARGET_VZ) -
             ((1.0 - beta) * T1);
      }
      L2 = relax_W * (qaux(i, j, k, QC) / problen[idir]) *
             (q(i, j, k, QU) - TARGET_VX) -
           ((1.0 - beta) * T2);
      L3 = relax_W * (qaux(i, j, k, QC) / problen[idir]) *
             (q(i, j, k, QV) - TARGET_VY) -
           ((1.0 - beta) * T3);
      L4 = relax_T *
             (q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QRSPEC) /
              problen[idir]) *
             (q(i, j, k, QTEMP) - TARGET_TEMPERATURE) -
           ((1.0 - beta) * T4);
    } else {
      amrex::Abort("Error:: Wait, is this the fourth dimension?");
    }
  } else if ((bc_type == 4) || (bc_type == 5)) {
    // Slipwall and NoSlipWall
    // Values long Y will be computed by mirror functions below
    // but we set waves values to 0 to avoid undefined variables
    L1 = 0.0;
    L2 = 0.0;
    L3 = 0.0;
    L4 = 0.0;
    L5 = 0.0;
  } else if (bc_type == 8) {
    // Outflow
    Kout = sigma_out * (1.0 - (mach_local * mach_local)) *
           (qaux(i, j, k, QC) / problen[idir]);
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
  if (idir == 0) {
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
  } else if (idir == 1) {
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
  } else if (idir == 2) {
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
  L[0] = L1;
  L[1] = L2;
  L[2] = L3;
  L[3] = L4;
  L[4] = L5;
}

#endif
