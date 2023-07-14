# Ensure FindPython3 supports NumPy
cmake_minimum_required(VERSION 3.14)

# Create main project.
project(pygetm VERSION 0.1.0 LANGUAGES Fortran C)

# Specify position-independent code since we will build a shared library.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

if (MSVC)
  add_compile_options(/libs:static)
endif()

# Find target Python interpreter.
set(Python_FIND_STRATEGY LOCATION)
set(Python3_FIND_STRATEGY LOCATION)
find_package(Python3 COMPONENTS Interpreter Development NumPy)

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/modules")

add_subdirectory(../extern/python-otps2 python-otps2 EXCLUDE_FROM_ALL)

add_subdirectory(../extern/pygsw pygsw EXCLUDE_FROM_ALL)
add_subdirectory(airsea pyairsea EXCLUDE_FROM_ALL)

set(FLEXOUT_USE_NetCDF OFF CACHE BOOL "Enable output in NetCDF format" FORCE)
mark_as_advanced(FLEXOUT_USE_NetCDF)
add_subdirectory(../extern/gotm/extern/flexout flexout EXCLUDE_FROM_ALL)

add_subdirectory(../extern/gotm/src/util gotm/util EXCLUDE_FROM_ALL)
add_subdirectory(../extern/gotm/src/turbulence gotm/turbulence EXCLUDE_FROM_ALL)
target_include_directories(turbulence PRIVATE ../extern/gotm/include)
target_include_directories(util PRIVATE ../extern/gotm/include)
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  target_compile_options(turbulence PRIVATE -ffree-line-length-none)
  target_compile_options(util PRIVATE -ffree-line-length-none)
endif()

set(PYFABM_DIR ${CMAKE_CURRENT_BINARY_DIR}/pyfabm CACHE STRING "Directory to write pyfabm library to" FORCE)
set(PYFABM_NAME "fabm" CACHE STRING "Base name of library" FORCE)
set(PYFABM_DIM_COUNT 3 CACHE STRING "Number of FABM dimensions" FORCE)
set(FABM_FORCED_REAL_KIND c_double)
add_compile_definitions(_FABM_DEPTH_DIMENSION_INDEX_=3 _FABM_MASK_TYPE_=integer _FABM_MASKED_VALUE_=0 _FABM_HORIZONTAL_MASK_ _FABM_VERTICAL_BOTTOM_TO_SURFACE_)
add_subdirectory(../extern/fabm/src/drivers/python fabm EXCLUDE_FROM_ALL)

set(GET_GIT_INFO_SCRIPT "../extern/fabm/src/cmake/Modules/GetGitInfo.cmake")
add_custom_target(pygetm_version
   ${CMAKE_COMMAND} -DINFILE=version.pxi.in -DOUTFILE=${CMAKE_CURRENT_BINARY_DIR}/version.pxi -P ${GET_GIT_INFO_SCRIPT}
   WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
   COMMENT "Retrieving description of last pygetm commit..."
   SOURCES version.pxi.in ${GET_GIT_INFO_SCRIPT}
   VERBATIM
)

add_library( flom 
  ../extern/flom/src/grid.F90
  ../extern/flom/src/logging.F90
  ../extern/flom/src/memory_manager.F90
)
set_property( TARGET flom PROPERTY FOLDER flom )

add_subdirectory(../src/domain domain)
add_subdirectory(../src/operators operators)
add_subdirectory(../src/dynamics dynamics)

add_library(cgetm
  src/wrapper.F90
  src/momentum.F90
  src/radiation.F90
  src/interpolation.F90
  src/internal_pressure.F90
)

target_link_libraries(cgetm getm_domain getm_operators getm_dynamics)

add_library(cgotm
  pygotm_wrapper.F90)
target_link_libraries(cgotm turbulence util yaml)

add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/_pygetm.c
  COMMAND cython -I ${CMAKE_CURRENT_BINARY_DIR} _pygetm.pyx
  MAIN_DEPENDENCY _pygetm.pyx
  DEPENDS pygetm_version ${CMAKE_CURRENT_BINARY_DIR}/version.pxi
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/_pygotm.c
  COMMAND cython _pygotm.pyx
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
                   COMMAND ${CMAKE_COMMAND} -E make_directory pygetm/pyfabm
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:fabm>" pygetm/pyfabm/
                   COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/../extern/fabm/src/pyfabm/__init__.py" pygetm/pyfabm/
                   COMMAND ${Python3_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/setup.py" bdist_wheel
                   VERBATIM)
add_custom_target(wheel ALL DEPENDS wheel.stamp)

add_dependencies(wheel otps2_wheel pygsw_wheel airsea_wheel)

# Let pip install the wheel
option(PYGETM_USER_INSTALL "Install with pip --user" OFF)
if(PYGETM_USER_INSTALL)
  SET(EXTRA_PIP_ARGS "--user")
endif()
install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install --no-index --find-links=${CMAKE_CURRENT_BINARY_DIR}/dist --find-links=${CMAKE_CURRENT_BINARY_DIR}/pyairsea/dist --find-links=${CMAKE_CURRENT_BINARY_DIR}/pygsw/dist --find-links=${CMAKE_CURRENT_BINARY_DIR}/python-otps2/dist ${EXTRA_PIP_ARGS} --force-reinstall --upgrade pygetm pyairsea pygsw otps2 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dist)")
