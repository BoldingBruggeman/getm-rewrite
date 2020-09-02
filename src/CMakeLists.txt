
set_property( TARGET getm_domain PROPERTY FOLDER domain )
target_sources( getm_domain PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/domain/domain.F90
    ${CMAKE_CURRENT_LIST_DIR}/domain/grids.F90
    ${CMAKE_CURRENT_LIST_DIR}/domain/metrics.F90
    ${CMAKE_CURRENT_LIST_DIR}/domain/register.F90
    ${CMAKE_CURRENT_LIST_DIR}/domain/uv_depths.F90
    ${CMAKE_CURRENT_LIST_DIR}/domain/cfl_check.F90
    ${CMAKE_CURRENT_LIST_DIR}/domain/depth_update.F90
    ${CMAKE_CURRENT_LIST_DIR}/domain/vertical_coordinates.F90
    ${CMAKE_CURRENT_LIST_DIR}/domain/vertical_gvc.F90
    ${CMAKE_CURRENT_LIST_DIR}/domain/vertical_sigma.F90
)

set_property( TARGET getm_airsea PROPERTY FOLDER airsea )
target_sources( getm_airsea PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/airsea/airsea.F90
)

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

set_property( TARGET getm_dynamics PROPERTY FOLDER dynamics )
target_sources( getm_dynamics PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/dynamics/dynamics.F90
    ${CMAKE_CURRENT_LIST_DIR}/dynamics/sealevel.F90
    ${CMAKE_CURRENT_LIST_DIR}/dynamics/pressure.F90
    ${CMAKE_CURRENT_LIST_DIR}/dynamics/pressure_surface.F90
    ${CMAKE_CURRENT_LIST_DIR}/dynamics/pressure_internal.F90
    ${CMAKE_CURRENT_LIST_DIR}/dynamics/momentum.F90
    ${CMAKE_CURRENT_LIST_DIR}/dynamics/momentum_2d.F90
    ${CMAKE_CURRENT_LIST_DIR}/dynamics/momentum_3d.F90
    ${CMAKE_CURRENT_LIST_DIR}/dynamics/velocities.F90
    #${CMAKE_CURRENT_LIST_DIR}/dynamics/velocity_diffusion.F90
    #${CMAKE_CURRENT_LIST_DIR}/dynamics/shear.F90
    ${CMAKE_CURRENT_LIST_DIR}/dynamics/slow_terms.F90
    ${CMAKE_CURRENT_LIST_DIR}/dynamics/momentum_register.F90
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

set_property( TARGET getm_operators PROPERTY FOLDER operators )
target_sources( getm_operators PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/operators/operators.F90
    ${CMAKE_CURRENT_LIST_DIR}/operators/advection.F90
    ${CMAKE_CURRENT_LIST_DIR}/operators/advection_upstream.F90
    ${CMAKE_CURRENT_LIST_DIR}/operators/diffusion.F90
)

#add_subdirectory(operators)

set_property( TARGET getm_exe PROPERTY FOLDER getm )
target_sources( getm_exe PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/main.F90
    ${CMAKE_CURRENT_LIST_DIR}/getm.F90
)
