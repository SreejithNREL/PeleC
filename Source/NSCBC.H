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
  T[0] = T1;
  T[1] = T2;
  T[2] = T3;
  T[3] = T4;
  T[5] = T5;
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
  if (bc_type == 7) {
    // Inflow
    if (idir == 1) {
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
    } else if (idir == 2) {
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
    } else if (idir == 3) {
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
  L[0] = L1;
  L[1] = L2;
  L[2] = L3;
  L[3] = L4;
  L[4] = L5;
}

#endif
