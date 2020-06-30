add_compile_options(-pg -g)

#add_executable(test_uv_depths EXCLUDE_FROM_ALL
#               ${CMAKE_CURRENT_LIST_DIR}/test_uv_depths.F90
#              )
#target_link_libraries(test_uv_depths PRIVATE getm_grid)
#set_property(TARGET test_uv_depths PROPERTY FOLDER grid_tests)


add_custom_target( test_all_loms )

add_executable(test_getm_domain EXCLUDE_FROM_ALL
               ${CMAKE_CURRENT_LIST_DIR}/test_getm_domain.F90
              )
target_link_libraries(test_getm_domain PRIVATE getm_input getm_domain flom )
add_dependencies( test_all_loms ${test_all_loms}  test_getm_domain )
#set_property(TARGET test_uv_depths PROPERTY FOLDER test_grids)

add_executable(test_bathymetry EXCLUDE_FROM_ALL
               ${CMAKE_CURRENT_LIST_DIR}/test_bathymetry.F90
              )
target_link_libraries(test_bathymetry PRIVATE getm_input getm_domain flom )
add_dependencies( test_all_loms ${test_all_loms} test_bathymetry  )

add_executable(test_density EXCLUDE_FROM_ALL
               ${CMAKE_CURRENT_LIST_DIR}/test_density.F90
              )
target_link_libraries(test_density PRIVATE getm_physics getm_domain flom)
add_dependencies( test_all_loms ${test_all_loms}  test_density )

add_executable(test_uv_depths EXCLUDE_FROM_ALL
               ${CMAKE_CURRENT_LIST_DIR}/test_uv_depths.F90
              )
target_link_libraries(test_uv_depths PRIVATE getm_domain flom)
add_dependencies( test_all_loms ${test_all_loms}  test_uv_depths )
