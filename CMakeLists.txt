cmake_minimum_required(VERSION 3.1.3)

# Create main project.
project(getm VERSION 0.1.0 LANGUAGES Fortran)

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Create placeholder file with version information
# used when STIM_EMBED_VERSION is off.
#set(GIT_COMMIT_ID ${VERSION})
#set(GIT_BRANCH_NAME unknown)
#configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/stim_version.F90.in" ${CMAKE_CURRENT_BINARY_DIR}/src/stim_version.F90)

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")

# Use use position-independent code (-fPIC) everywhere if building shared libraries
if(BUILD_SHARED_LIBS)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

option(GETM_USE_PARALLEL "Enable MPI parallelization" OFF)
if(GETM_USE_PARALLEL)
   find_package(MPI REQUIRED)
   if(MPI_Fortran_INCLUDE_PATH)
     include_directories("${MPI_Fortran_INCLUDE_PATH}")
   endif()
   add_definitions(-D_PARALLEL_)
   set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${MPI_Fortran_LIBRARIES})
endif(GETM_USE_PARALLEL)


add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/extern/flom" extern/flom)

# Various compilation options
if( NOT GETM_USE_STATIC )
option(GETM_USE_STATIC "Static compilation (domain size is set at compile time)" OFF)
endif()

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  add_compile_options(-ffree-line-length-none)
endif()
add_compile_options(-pg -g)
#target_link_options(add PRIVATE -pg -g)

#add_link_options(-pg -g)

include_directories("${CMAKE_CURRENT_LIST_DIR}/include")

# Libraries for the GETM 3D model
#add_library(getm_variables "")
#target_link_libraries( getm_variables PUBLIC getm_domain )
add_library(getm_airsea "")
target_link_libraries( getm_airsea PRIVATE getm_domain )
add_library(getm_physics "")
target_link_libraries( getm_physics PRIVATE getm_airsea getm_operators getm_domain flom )
add_library(getm_dynamics "")
target_link_libraries( getm_dynamics PRIVATE getm_operators getm_domain flom )
add_library(getm_input "")
target_link_libraries( getm_input PRIVATE getm_domain getm_airsea getm_dynamics flom )
add_library(getm_output "")
target_link_libraries( getm_output PRIVATE flom )
# .. and now the executable
add_executable(getm_exe)

#include(src/CMakeLists.txt)
add_subdirectory(src)
include(tests/CMakeLists.txt)

target_include_directories( getm_exe PRIVATE ${CMAKE_BINARY_DIR}/modules )
target_link_libraries( getm_exe PRIVATE getm_domain getm_parallel getm_operators getm_physics getm_dynamics getm_input getm_output flom )


#install(EXPORT stimConfig DESTINATION cmake)
