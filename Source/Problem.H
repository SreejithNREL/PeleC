//#ifndef _PROBLEM_H_
//#define _PROBLEM_H_

// problem-specific PeleC:: declarations go here

void problem_post_timestep();

void problem_post_restart();

void problem_post_init();
void problem_post_init( int i,
 int j,
 int k,
 amrex::Array4<amrex::Real> const& state,
 amrex::GeometryData const& geomdata);

//#endif
