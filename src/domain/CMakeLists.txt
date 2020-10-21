add_library( getm_domain
    domain.F90
    grids.F90
    metrics.F90
    register.F90
    uv_depths.F90
    cfl_check.F90
    depth_update.F90
    vertical_coordinates.F90
    vertical_gvc.F90
    vertical_sigma.F90
)
set_property( TARGET getm_domain PROPERTY FOLDER domain )
target_link_libraries(getm_domain flom field_manager)