add_library( getm_parallel
    parallel.F90
)
if(GETM_USE_PARALLEL)
  target_link_libraries(getm_parallel PUBLIC MPI::MPI_Fortran)
endif()
set_property( TARGET getm_parallel PROPERTY FOLDER parallel )
