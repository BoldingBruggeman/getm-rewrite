cmake_minimum_required(VERSION 3.1.3)

# Create main project.
project(airsea VERSION 0.1.0 LANGUAGES Fortran C)


# Specify position-independent code since we will build a shared library.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

if (MSVC)
  add_compile_options(/libs:static)
endif()

add_library(cairsea
   src/airsea_fluxes.F90
   src/airsea_variables.F90
   src/humidity.F90
   src/longwave_radiation.F90
   src/fairall.F90
   src/t10m.F90
   src/kara.F90
   src/kondo.F90
   src/solar_zenith_angle.F90
   src/albedo_water.F90
   src/shortwave_radiation.F90
   src/airsea_2d.F90
)
target_compile_definitions(cairsea PRIVATE GUSTINESS)

# Find target Python interpreter.
set(Python_FIND_STRATEGY LOCATION)
find_package(Python3 COMPONENTS Interpreter Development NumPy)

add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/_pyairsea.c
  COMMAND ${Python3_EXECUTABLE} -m cython _pyairsea.pyx
  MAIN_DEPENDENCY _pyairsea.pyx
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

Python3_add_library(_pyairsea MODULE
  _pyairsea.c)
target_link_libraries(_pyairsea PRIVATE cairsea Python3::NumPy)

# Build wheel
add_custom_command(OUTPUT airsea_wheel.stamp
                   COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/pyairsea" pyairsea
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:_pyairsea>" pyairsea/
                   COMMAND ${Python3_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/setup.py" bdist_wheel
                   VERBATIM)
add_custom_target(airsea_wheel ALL DEPENDS airsea_wheel.stamp)

# Let pip install the wheel
install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install --no-index --find-links=${CMAKE_CURRENT_BINARY_DIR}/dist --user --force-reinstall --upgrade pyairsea WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dist)")
