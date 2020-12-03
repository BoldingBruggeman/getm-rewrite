add_library(getm_dynamics
    dynamics.F90
    sealevel.F90
    pressure.F90
    pressure_surface.F90
    pressure_internal.F90
    ip_blumberg_mellor.F90
    ip_shchepetkin_mcwilliams.F90
    momentum.F90
    momentum_2d.F90
    momentum_3d.F90
    uv_advection.F90
    uv_diffusion.F90
    velocities.F90
    #shear.F90
    slow_terms.F90
    momentum_register.F90
)
set_property( TARGET getm_dynamics PROPERTY FOLDER getm )
target_link_libraries( getm_dynamics PRIVATE getm_operators getm_domain flom )
