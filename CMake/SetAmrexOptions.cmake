set(USE_XSDK_DEFAULTS OFF)
set(AMReX_SPACEDIM "${PELEC_DIM}" CACHE STRING "Number of physical dimensions" FORCE)
set(AMReX_MPI ${PELEC_ENABLE_MPI})
set(AMReX_OMP ${PELEC_ENABLE_OPENMP})
set(AMReX_EB ${PELEC_ENABLE_AMREX_EB})
set(AMReX_PARTICLES ${PELEC_ENABLE_PARTICLES})
set(AMReX_CUDA ${PELEC_ENABLE_CUDA})
set(AMReX_DPCPP ${PELEC_ENABLE_DPCPP})
set(AMReX_HIP ${PELEC_ENABLE_HIP})
set(AMReX_PLOTFILE_TOOLS ${PELEC_ENABLE_FCOMPARE})
set(AMReX_SUNDIALS OFF)
set(AMReX_FORTRAN OFF)
set(AMReX_FORTRAN_INTERFACES OFF)
set(AMReX_PIC OFF)
set(AMReX_PRECISION "${PELEC_PRECISION}" CACHE STRING "Floating point precision" FORCE)
set(AMReX_LINEAR_SOLVERS OFF)
set(AMReX_AMRDATA OFF)
set(AMReX_ASCENT OFF)
set(AMReX_SENSEI OFF)
set(AMReX_CONDUIT OFF)
set(AMReX_HYPRE OFF)
set(AMReX_FPE OFF)
set(AMReX_ASSERTIONS OFF)
set(AMReX_BASE_PROFILE OFF)
set(AMReX_TINY_PROFILE ${PELEC_ENABLE_TINY_PROFILE})
set(AMReX_TRACE_PROFILE OFF)
set(AMReX_MEM_PROFILE OFF)
set(AMReX_COMM_PROFILE OFF)
set(AMReX_BACKTRACE OFF)
set(AMReX_PROFPARSER OFF)
set(AMReX_ACC OFF)
set(AMReX_INSTALL OFF)

if(PELEC_ENABLE_CUDA)
  set(AMReX_GPU_BACKEND CUDA CACHE STRING "AMReX GPU type" FORCE)
elseif(PELEC_ENABLE_HIP)
  set(AMReX_GPU_BACKEND HIP CACHE STRING "AMReX GPU type" FORCE)
elseif(PELEC_ENABLE_DPCPP)
  set(AMReX_GPU_BACKEND SYCL CACHE STRING "AMReX GPU type" FORCE)
else()
  set(AMReX_GPU_BACKEND NONE CACHE STRING "AMReX GPU type" FORCE)
endif()
