#KBadd_subdirectory(parallel)

add_subdirectory(domain)
add_subdirectory(operators)
add_subdirectory(airsea)
add_subdirectory(dynamics)
add_subdirectory(physics)
add_subdirectory(input)
add_subdirectory(output)

add_executable( getm_exe
    main.F90
    getm.F90
)
set_property( TARGET getm_exe PROPERTY FOLDER getm )
target_include_directories( getm_exe PRIVATE ${CMAKE_BINARY_DIR}/modules )
target_link_libraries( getm_exe PRIVATE getm_domain getm_operators getm_physics getm_dynamics getm_input getm_output flom )
