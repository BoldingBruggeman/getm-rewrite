add_library( getm_airsea
    airsea.F90
)
set_property( TARGET getm_airsea PROPERTY FOLDER getm )
target_link_libraries( getm_airsea PRIVATE getm_domain )
