add_library( getm_domain
    domain.F90
    grids.F90
    metrics.F90
    register.F90
    uvx_depths.F90
    mirror_bdys.F90
    cfl_check.F90
    depth_update.F90
    vertical_coordinates.F90
    vertical_gvc.F90
    vertical_sigma.F90
)
set_property( TARGET getm_domain PROPERTY FOLDER getm )
target_link_libraries(getm_domain flom field_manager)
