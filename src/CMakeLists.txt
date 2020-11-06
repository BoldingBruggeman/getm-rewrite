add_subdirectory(parallel)

add_subdirectory(domain)
add_subdirectory(operators)
add_subdirectory(airsea)
add_subdirectory(dynamics)

set_property( TARGET getm_physics PROPERTY FOLDER physics )
target_sources( getm_physics PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/physics/physics.F90
    ${CMAKE_CURRENT_LIST_DIR}/physics/radiation.F90
    ${CMAKE_CURRENT_LIST_DIR}/physics/salinity.F90
    ${CMAKE_CURRENT_LIST_DIR}/physics/temperature.F90
    ${CMAKE_CURRENT_LIST_DIR}/physics/density.F90
    ${CMAKE_CURRENT_LIST_DIR}/physics/mixing.F90
#    ${CMAKE_CURRENT_LIST_DIR}/physics/diffusion.F90
)

set_property( TARGET getm_input PROPERTY FOLDER input )
target_sources( getm_input PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/input/input.F90
    ${CMAKE_CURRENT_LIST_DIR}/input/bathymetry.F90
)

set_property( TARGET getm_output PROPERTY FOLDER output )
target_sources( getm_output PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/output/output.F90
)

set_property( TARGET getm_exe PROPERTY FOLDER getm )
target_sources( getm_exe PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/main.F90
    ${CMAKE_CURRENT_LIST_DIR}/getm.F90
)
