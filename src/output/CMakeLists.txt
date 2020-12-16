add_library( getm_output
    output.F90
)
set_property( TARGET getm_output PROPERTY FOLDER output )
target_link_libraries( getm_output PRIVATE flom )
