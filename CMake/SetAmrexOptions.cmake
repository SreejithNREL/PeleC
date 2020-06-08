#Set amrex options
set(USE_XSDK_DEFAULTS OFF)
set(DIM ${PELEC_DIM})
set(ENABLE_PIC OFF)
set(ENABLE_MPI ${PELEC_ENABLE_MPI})
set(ENABLE_OMP ${PELEC_ENABLE_OPENMP})
set(ENABLE_DP ON)
set(ENABLE_EB ${PELEC_ENABLE_EB})
set(ENABLE_FORTRAN_INTERFACES OFF)
set(ENABLE_LINEAR_SOLVERS OFF)
set(ENABLE_AMRDATA OFF)
set(ENABLE_PARTICLES ${PELEC_ENABLE_PARTICLES})
set(ENABLE_SENSEI_INSITU OFF)
set(ENABLE_CONDUIT OFF)
set(ENABLE_SUNDIALS OFF)
set(ENABLE_FPE OFF)
set(ENABLE_ASSERTIONS OFF)
set(ENABLE_BASE_PROFILE OFF)
set(ENABLE_TINY_PROFILE OFF)
set(ENABLE_TRACE_PROFILE OFF)
set(ENABLE_MEM_PROFILE OFF)
set(ENABLE_COMM_PROFILE OFF)
set(ENABLE_BACKTRACE OFF)
set(ENABLE_PROFPARSER OFF)
set(ENABLE_CUDA ${PELEC_ENABLE_CUDA})
set(ENABLE_ACC OFF)
set(ENABLE_PLOTFILE_TOOLS ${PELEC_ENABLE_FCOMPARE})
#set(ENABLE_FORTRAN OFF)