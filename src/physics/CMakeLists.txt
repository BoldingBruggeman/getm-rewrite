add_library( getm_physics
    physics.F90
    radiation.F90
    salinity.F90
    temperature.F90
    density.F90
    mixing.F90
)
set_property( TARGET getm_physics PROPERTY FOLDER physics )
target_link_libraries( getm_physics PRIVATE getm_airsea getm_operators getm_domain flom turbulence )
