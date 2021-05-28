add_library( getm_operators
    operators.F90
    advection.F90
    advection_schemes/advection_superbee.F90
    advection_schemes/advection_p2_pdm.F90
    advection_schemes/advection_splmax13.F90
    advection_schemes/advection_hsimt.F90
    advection_schemes/advection_muscl.F90
    advection_schemes/advection_upstream.F90
    diffusion.F90
)
target_link_libraries(getm_operators getm_domain)
set_property( TARGET getm_operators PROPERTY FOLDER getm )
