set_property( TARGET getm_operators PROPERTY FOLDER operators )
target_sources( getm_operators PRIVATE
    getm_operators.F90
    advection.F90
#    advection_upstream.F90
    diffusion.F90
)

