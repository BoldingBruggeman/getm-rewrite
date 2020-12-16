add_library( getm_input
    input.F90
    bathymetry.F90
)
set_property( TARGET getm_input PROPERTY FOLDER input )
target_link_libraries( getm_input PRIVATE getm_domain getm_airsea getm_dynamics flom )

