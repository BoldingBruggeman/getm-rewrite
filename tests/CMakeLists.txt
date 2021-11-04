add_compile_options(-pg -g)

#add_executable(test_uv_depths EXCLUDE_FROM_ALL
#               ${CMAKE_CURRENT_LIST_DIR}/test_uv_depths.F90
#              )
#target_link_libraries(test_uv_depths PRIVATE getm_grid)
#set_property(TARGET test_uv_depths PROPERTY FOLDER grid_tests)


add_custom_target( test_all_getm )

add_executable(test_vertical_coordinates EXCLUDE_FROM_ALL
               ${CMAKE_CURRENT_LIST_DIR}/test_vertical_coordinates.F90
              )
target_link_libraries(test_vertical_coordinates PRIVATE getm_domain getm_dynamics flom )

add_executable(test_getm_domain EXCLUDE_FROM_ALL
               ${CMAKE_CURRENT_LIST_DIR}/test_getm_domain.F90
              )
target_link_libraries(test_getm_domain PRIVATE getm_input getm_domain flom )
#KB add_dependencies( test_all_getm ${test_all_getm}  test_getm_domain )
#set_property(TARGET test_uv_depths PROPERTY FOLDER test_grids)

add_executable(test_bathymetry EXCLUDE_FROM_ALL
               ${CMAKE_CURRENT_LIST_DIR}/test_bathymetry.F90
              )
target_link_libraries(test_bathymetry PRIVATE getm_input getm_domain flom )
add_dependencies( test_all_getm ${test_all_getm} test_bathymetry  )

add_executable(test_density EXCLUDE_FROM_ALL
               ${CMAKE_CURRENT_LIST_DIR}/test_density.F90
              )
target_link_libraries(test_density PRIVATE getm_physics getm_domain flom)
add_dependencies( test_all_getm ${test_all_getm}  test_density )

add_executable(test_uvx_depths EXCLUDE_FROM_ALL
               ${CMAKE_CURRENT_LIST_DIR}/test_uvx_depths.F90
              )
target_link_libraries(test_uvx_depths PRIVATE getm_domain flom)
add_dependencies( test_all_getm ${test_all_getm}  test_uvx_depths )

add_executable(test_advection EXCLUDE_FROM_ALL
               ${CMAKE_CURRENT_LIST_DIR}/test_advection.F90
              )
target_link_libraries(test_advection PRIVATE getm_output getm_operators getm_domain flom)
add_dependencies( test_all_getm ${test_all_getm}  test_advection )

add_executable(test_vertical_diffusion EXCLUDE_FROM_ALL
               ${CMAKE_CURRENT_LIST_DIR}/test_vertical_diffusion.F90
              )
target_link_libraries(test_vertical_diffusion PRIVATE getm_output getm_operators  getm_physics flom)
add_dependencies( test_all_getm ${test_all_getm}  test_vertical_diffusion )

add_executable(test_parallel EXCLUDE_FROM_ALL
               ${CMAKE_CURRENT_LIST_DIR}/test_parallel.F90
              )
target_link_libraries(test_parallel PRIVATE getm_domain getm_parallel flom)
#KBadd_dependencies( test_all_getm ${test_all_getm}  test_parallel )

