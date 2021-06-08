cmake_minimum_required(VERSION 3.1.3)

# Create main project.
project(pygetm VERSION 0.1.0 LANGUAGES Fortran)

# Specify position-independent code since we will build a shared library.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/modules")

add_subdirectory(../extern/flom/extern/flexout/src flexout EXCLUDE_FROM_ALL)

add_subdirectory(../extern/python-otps2 python-otps2)

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

add_library(cgetm SHARED
  wrapper.F90
  pygetm/__init__.py
)
target_link_libraries(cgetm getm_domain getm_operators getm_airsea getm_dynamics)

# Find target Python interpreter.
find_package(PythonInterp)

# Build wheel
add_custom_command(OUTPUT wheel.stamp
                   COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/pygetm" pygetm
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_LINKER_FILE:cgetm>" pygetm/
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:cgetm>" pygetm/
                   COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/setup.py" build_ext --library-dirs=pygetm --libraries=$<TARGET_FILE_BASE_NAME:cgetm> --rpath=$ORIGIN bdist_wheel
                   VERBATIM)
add_custom_target(wheel ALL DEPENDS wheel.stamp)

# Let pip install the wheel
install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install --no-index --find-links=${CMAKE_CURRENT_BINARY_DIR}/dist --user --force-reinstall --upgrade pygetm WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dist)")
