set_property( TARGET getm_operators PROPERTY FOLDER operators )
set(USE_KB on)
if (USE_KB)
target_sources( getm_operators PRIVATE
    getm_operators_kb.F90
    advection_kb.F90
#    advection_upstream.F90
    diffusion.F90
)
else()
target_sources( getm_operators PRIVATE
    getm_operators.F90
    advection.F90
#    advection_upstream.F90
    diffusion.F90
)
endif()

