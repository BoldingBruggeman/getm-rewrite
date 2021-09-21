cmake_minimum_required(VERSION 3.1.3)

# Create main project.
project(pygetm VERSION 0.1.0 LANGUAGES Fortran C)

# Specify position-independent code since we will build a shared library.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

if (MSVC)
  add_compile_options(/libs:static)
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/modules")

add_subdirectory(../extern/python-otps2 python-otps2)

add_subdirectory(../extern/pygsw pygsw)

set(GOTM_BUILD_LIBRARIES_ONLY ON)
set(GOTM_USE_FABM OFF CACHE BOOL "Use FABM" FORCE)
set(GOTM_USE_NetCDF OFF CACHE BOOL "Use NetCDF" FORCE)
add_subdirectory(../extern/gotm gotm EXCLUDE_FROM_ALL)

add_library( flom 
	../extern/flom/src/grid.F90
	../extern/flom/src/logging.F90
	../extern/flom/src/memory_manager.F90
	)
set_property( TARGET flom PROPERTY FOLDER flom )

add_subdirectory(../src/domain domain)
add_subdirectory(../src/airsea airsea)
add_subdirectory(../src/operators operators)
add_subdirectory(../src/dynamics dynamics)

add_library(cgetm
  wrapper.F90)

target_link_libraries(cgetm getm_domain getm_operators getm_airsea getm_dynamics)

add_library(cgotm
  pygotm_wrapper.F90)
target_link_libraries(cgotm turbulence util yaml)

# Find target Python interpreter.
find_package(Python3 COMPONENTS Interpreter Development NumPy)

add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/_pygetm.c
  COMMAND ${Python3_EXECUTABLE} -m cython _pygetm.pyx
  MAIN_DEPENDENCY _pygetm.pyx
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/_pygotm.c
  COMMAND ${Python3_EXECUTABLE} -m cython _pygotm.pyx
  MAIN_DEPENDENCY _pygotm.pyx
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

Python3_add_library(_pygetm MODULE
  _pygetm.c)
target_link_libraries(_pygetm PRIVATE cgetm Python3::NumPy)


Python3_add_library(_pygotm MODULE
  _pygotm.c)
target_link_libraries(_pygotm PRIVATE cgotm Python3::NumPy)

# Build wheel
add_custom_command(OUTPUT wheel.stamp
                   COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/pygetm" pygetm
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:_pygetm>" pygetm/
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:_pygotm>" pygetm/
                   COMMAND ${Python3_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/setup.py" bdist_wheel
                   VERBATIM)
add_custom_target(wheel ALL DEPENDS wheel.stamp)

# Let pip install the wheel
install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install --no-index --find-links=${CMAKE_CURRENT_BINARY_DIR}/dist --user --force-reinstall --upgrade pygetm WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dist)")
