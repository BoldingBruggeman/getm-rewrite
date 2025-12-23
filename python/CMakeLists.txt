# Ensure FindPython3 supports NumPy and CMP0094 (3.15) and WITH_SOABI (3.17)
cmake_minimum_required(VERSION 3.17)

# Create main project.
project(pygetm VERSION 0.9.18 LANGUAGES Fortran C)

# Specify position-independent code since we will build a shared library.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

if (WIN32 AND "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
  add_compile_options(/libs:static)
endif()

# Find target Python interpreter.
find_package(Python3 COMPONENTS Interpreter Development NumPy)

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/modules")

set(FLEXOUT_USE_NetCDF OFF CACHE BOOL "Enable output in NetCDF format" FORCE)
mark_as_advanced(FLEXOUT_USE_NetCDF)
add_subdirectory(../extern/flexout flexout EXCLUDE_FROM_ALL)

set(GET_GIT_INFO_SCRIPT "../extern/fabm/src/cmake/Modules/GetGitInfo.cmake")
add_custom_target(pygetm_version
   ${CMAKE_COMMAND} -DINFILE=version.pxi.in -DOUTFILE=${CMAKE_CURRENT_BINARY_DIR}/version.pxi -P ${GET_GIT_INFO_SCRIPT}
   WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
   COMMENT "Retrieving description of last pygetm commit..."
   SOURCES version.pxi.in ${GET_GIT_INFO_SCRIPT}
   VERBATIM
)

add_subdirectory(../extern/python-otps2 python-otps2 EXCLUDE_FROM_ALL)
add_subdirectory(../extern/pygsw pygsw EXCLUDE_FROM_ALL)

add_custom_target(pip_install pip install -v "${CMAKE_SOURCE_DIR}" VERBATIM)
add_dependencies(pip_install otps2_wheel pygsw_wheel)

add_subdirectory(../src/memory memory)
add_subdirectory(../src/domain domain)
add_subdirectory(../src/operators operators)

add_library(cgetm
  src/wrapper.F90
  src/momentum.F90
  src/radiation.F90
  src/interpolation.F90
  src/gradient.F90
  src/internal_pressure.F90
)

target_link_libraries(cgetm getm_domain getm_operators)

add_custom_command(OUTPUT _pygetm.c
  COMMENT "Creating ${CMAKE_CURRENT_BINARY_DIR}/_pygetm.c from ${CMAKE_CURRENT_SOURCE_DIR}/_pygetm.pyx"
  COMMAND Python3::Interpreter -m cython -I "${CMAKE_CURRENT_BINARY_DIR}" _pygetm.pyx -o "${CMAKE_CURRENT_BINARY_DIR}/_pygetm.c"
  MAIN_DEPENDENCY _pygetm.pyx
  DEPENDS pygetm_version
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  VERBATIM
)

Python3_add_library(_pygetm MODULE WITH_SOABI
  "${CMAKE_CURRENT_BINARY_DIR}/_pygetm.c")
target_link_libraries(_pygetm PRIVATE cgetm Python3::NumPy)

install(TARGETS _pygetm DESTINATION .)
