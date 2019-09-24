#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(TGF_propa C CXX Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/src/cmake/modules/) # path to custom .cmake files to find packages (e.g. matio)

execute_process(COMMAND "source ~/.bashrc")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if (WITH_GEANT4_UIVIS)
    find_package(Geant4 REQUIRED ui_all vis_all)
else ()
    find_package(Geant4 REQUIRED)
endif ()

#-----------------MATIO-----------------------------------------------------
#find_package(MATIO REQUIRED) # call to the .cmake script located at ${CMAKE_SOURCE_DIR}/src/cmake/modules/
#file(GLOB matio_headers ${MATIO_INCLUDE_DIR}/*.h)

#set(CMAKE_BUILD_TYPE Release)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#-----------------MATLAB-----------------------------------------------------
#set(Matlab_ROOT_DIR "/scratch/Matlab/R2017b/")
#find_package(Matlab REQUIRED MAT_LIBRARY)
#file(GLOB matlab_headers ${Matlab_INCLUDE_DIRS}/*.h)
#${Matlab_LIBRARIES}

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/src/include ${Geant4_INCLUDE_DIR}/src ${Geant4_INCLUDE_DIR}/src/src  ${PROJECT_SOURCE_DIR}/src/src/C/WMM/src/ ${PROJECT_SOURCE_DIR}/src/src/C/MSIS/)
#include_directories(${Matlab_INCLUDE_DIRS})

file(GLOB main ${PROJECT_SOURCE_DIR}/src/tgf_propa.cc)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/src/*.cc)
file(GLOB sourceshh ${PROJECT_SOURCE_DIR}/src/include/*.hh)

file(GLOB sourcesc ${PROJECT_SOURCE_DIR}/src/src/C/*.c)

if (UNIX)
    file(GLOB sourcesf ${PROJECT_SOURCE_DIR}/src/src/fortran/*.f)
endif (UNIX)

file(GLOB sources_mag_c ${PROJECT_SOURCE_DIR}/src/src/C/WMM/src/*.c)
file(GLOB sources_mag_h ${PROJECT_SOURCE_DIR}/src/src/C/WMM/src/*.h)

file(GLOB sources_msis_c ${PROJECT_SOURCE_DIR}/src/src/C/MSIS/*.c)
file(GLOB sources_msis_h ${PROJECT_SOURCE_DIR}/src/src/C/MSIS/*.h)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries

if (UNIX)
    add_executable(TGF_Propa ${main} ${sources_mag_h} ${sourcesf} ${sources_msis_c} ${sources_msis_h} ${sources_mag_c} ${sourcesc} ${sourceshh} ${sources})
endif (UNIX)

if (WIN32)
    add_executable(TGF_Propa ${main} ${sources_msis_c} ${sources_msis_h} ${sources_mag_c} ${sources_mag_h} ${sourcesc} ${sourceshh} ${sources})
endif (WIN32)

target_link_libraries(TGF_Propa ${Geant4_LIBRARIES})

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # using GCC
    target_compile_options(TGF_Propa PRIVATE -std=c++11 -mfpmath=sse -msse2 -ffp-contract=off)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    # using Intel C++
    target_compile_options(TGF_Propa PRIVATE -fp-model=extended)
endif ()

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(TGF_propa DEPENDS TGF_Propa)
