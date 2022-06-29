set(SUNDIALS_PRECISION "${PELEC_PRECISION}" CACHE STRING "Floating point precision" FORCE)
set(BUILD_ARKODE ON)
set(BUILD_CVODE ON)
set(BUILD_SUNMATRIX_CUSPARSE ${PELEC_ENABLE_CUDA})
set(BUILD_SUNLINSOL_CUSOLVERSP ${PELEC_ENABLE_CUDA})
set(SUNDIALS_INDEX_SIZE 32)
set(BUILD_CVODES OFF)
set(BUILD_IDA OFF)
set(BUILD_IDAS OFF)
set(BUILD_KINSOL OFF)
set(BUILD_CPODES OFF)
set(BUILD_EXAMPLES OFF)
set(_BUILD_EXAMPLES OFF)
set(BUILD_TESTING OFF)
set(ENABLE_CUDA ${PELEC_ENABLE_CUDA})
set(ENABLE_HIP ${PELEC_ENABLE_HIP})
set(ENABLE_SYCL ${PELEC_ENABLE_DPCPP})
set(EXAMPLES_ENABLE_C OFF)
set(EXAMPLES_ENABLE_CXX OFF)
set(EXAMPLES_ENABLE_CUDA OFF)
