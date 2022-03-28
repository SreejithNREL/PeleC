#ifndef _TAGGING_H_
#define _TAGGING_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_TagBox.H>
#include "prob_parm.H"

struct TaggingParm
{
  amrex::Real denerr = 1.0e10;
  amrex::Real dengrad = 1.0e10;
  amrex::Real presserr = 1.0e10;
  amrex::Real pressgrad = 1.0e10;
  amrex::Real velerr = 1.0e10;
  amrex::Real velgrad = 1.0e10;
  amrex::Real vorterr = 1.0e10;
  amrex::Real temperr = 1.0e10;
  amrex::Real tempgrad = 1.0e10;
  amrex::Real ftracerr = 1.0e10;
  amrex::Real ftracgrad = 1.0e10;
  amrex::Real vfracerr = 1.0e10;

  std::string eb_refine_type = "static";

  int max_denerr_lev = 10;
  int max_dengrad_lev = 10;
  int max_presserr_lev = 10;
  int max_pressgrad_lev = 10;
  int max_velerr_lev = 10;
  int max_velgrad_lev = 10;
  int max_vorterr_lev = 10;
  int max_temperr_lev = 10;
  int max_tempgrad_lev = 10;
  int max_ftracerr_lev = 10;
  int max_ftracgrad_lev = 10;
  int max_vfracerr_lev = 10;


  amrex::Real dom_lbnd_x = 1.0e10;
  amrex::Real dom_hbnd_x = 1.0e10;
  amrex::Real dom_lbnd_y = 1.0e10;
  amrex::Real dom_hbnd_y = 1.0e10;
  amrex::Real dom_lbnd_z = 1.0e10;
  amrex::Real dom_hbnd_z = 1.0e10;
  amrex::Vector<amrex::Real> dom_lbnd[AMREX_SPACEDIM] = {{1.0e10}};
 // amrex::Vector<amrex::Real> dom_hbnd[AMREX_SPACEDIM] = {{1.0e10}};
  int max_dom_bnd_lev = 10;

  int max_eb_refine_lev = 10;
  int min_eb_refine_lev = -1;
  int adapt_eb_refined_lev = -1;

};

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
tag_error_region(
  const int i,
  const int j,
  const int k,
  amrex::GeometryData const& geomdata,
  amrex::Real x_lo,
  amrex::Real x_hi,
  amrex::Real y_lo,
  amrex::Real y_hi,
  amrex::Real z_lo,
  amrex::Real z_hi,
  amrex::Array4<char> const& tag,
  char tagval) noexcept
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  // const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];


  // Tag on regions of high field values
  if ((x_lo < x) && (x < x_hi) && (y_lo < y) && (y < y_hi)
#if AMREX_SPACEDIM == 3
      && (z_lo < z) && (z < z_hi)
  //    (dom_lbnd[2] < z) && (z < dom_hbnd[2])) {
#endif
     ){
        tag(i, j, k) = tagval;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
tag_error(
  const int i,
  const int j,
  const int k,
  amrex::Array4<char> const& tag,
  amrex::Array4<amrex::Real const> const& field,
  const amrex::Real fielderr,
  char tagval) noexcept
{
  // Tag on regions of high field values
  if (field(i, j, k) >= fielderr) {
    tag(i, j, k) = tagval;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
tag_graderror(
  const int i,
  const int j,
  const int k,
  amrex::Array4<char> const& tag,
  amrex::Array4<amrex::Real const> const& field,
  const amrex::Real fieldgrad,
  char tagval) noexcept
{
  // Tag on regions of high field gradient
  AMREX_D_TERM(
    amrex::Real ax = amrex::Math::abs(field(i + 1, j, k) - field(i, j, k));
    ax = amrex::max<amrex::Real>(
      ax, amrex::Math::abs(field(i, j, k) - field(i - 1, j, k)));
    , amrex::Real ay = amrex::Math::abs(field(i, j + 1, k) - field(i, j, k));
    ay = amrex::max<amrex::Real>(
      ay, amrex::Math::abs(field(i, j, k) - field(i, j - 1, k)));
    , amrex::Real az = amrex::Math::abs(field(i, j, k + 1) - field(i, j, k));
    az = amrex::max<amrex::Real>(
      az, amrex::Math::abs(field(i, j, k) - field(i, j, k - 1)));)
#if AMREX_SPACEDIM > 1
  if (amrex::max<amrex::Real>(AMREX_D_DECL(ax, ay, az)) >= fieldgrad)
#elif AMREX_SPACEDIM == 1
  if (ax >= fieldgrad)
#endif
  {
    tag(i, j, k) = tagval;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
tag_abserror(
  const int i,
  const int j,
  const int k,
  amrex::Array4<char> const& tag,
  amrex::Array4<amrex::Real const> const& field,
  const amrex::Real fielderr,
  char tagval) noexcept
{
  // Tag on regions of high field values
  if (amrex::Math::abs(field(i, j, k)) >= fielderr) {
    tag(i, j, k) = tagval;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
tag_error_bounds(
  const int i,
  const int j,
  const int k,
  amrex::Array4<char> const& tag,
  amrex::Array4<amrex::Real const> const& field,
  const amrex::Real lbnd,
  const amrex::Real ubnd,
  char tagval) noexcept
{
  // Tag on regions inside bounds
  if ((lbnd < field(i, j, k)) && (field(i, j, k) < ubnd)) {
    tag(i, j, k) = tagval;
  }
}

struct EmptyProbTagStruct
{
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void set_problem_tags(
    const int /*i*/,
    const int /*j*/,
    const int /*k*/,
    amrex::Array4<char> const& /*tag*/,
    amrex::Array4<amrex::Real const> const& /*field*/,
    char /*tagval*/,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> /*dx*/,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> /*prob_lo*/,
    const amrex::Real /*time*/,
    const int /*level*/,
    ProbParmDevice const& /*d_prob_parm_device*/) noexcept
  {
  }
};

template <typename ProbTagStruct, typename... Args>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
set_problem_tags(Args&&... args)
{
  ProbTagStruct::set_problem_tags(std::forward<Args>(args)...);
}

#endif
