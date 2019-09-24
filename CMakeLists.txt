#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.6 FATAL_ERROR)

if (UNIX AND NOT APPLE)
    project(TGF_propa C CXX Fortran)
endif (UNIX AND NOT APPLE)

if (APPLE)
    project(TGF_propa C CXX)
endif (APPLE)

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
find_package(PkgConfig REQUIRED)
find_package(Threads REQUIRED)
pkg_search_module(UUID REQUIRED uuid)

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
include_directories(${PROJECT_SOURCE_DIR}/src/libsgp4 ${PROJECT_SOURCE_DIR}/src/src/C/UUID/ ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/src/include ${Geant4_INCLUDE_DIR}/src ${Geant4_INCLUDE_DIR}/src/src ${PROJECT_SOURCE_DIR}/src/src/C/WMM/src/ ${PROJECT_SOURCE_DIR}/src/src/C/MSIS/)
#include_directories(${Matlab_INCLUDE_DIRS})

file(GLOB main ${PROJECT_SOURCE_DIR}/src/tgf_propa.cc)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/src/*.cc)
file(GLOB sourceshh ${PROJECT_SOURCE_DIR}/src/include/*.hh)

file(GLOB sourcesc ${PROJECT_SOURCE_DIR}/src/src/C/*.c)

if (UNIX AND NOT APPLE)
    file(GLOB sourcesf ${PROJECT_SOURCE_DIR}/src/src/fortran/*.f)
endif (UNIX AND NOT APPLE)

file(GLOB sources_mag_c ${PROJECT_SOURCE_DIR}/src/src/C/WMM/src/*.c)
file(GLOB sources_mag_h ${PROJECT_SOURCE_DIR}/src/src/C/WMM/src/*.h)

file(GLOB sources_msis_c ${PROJECT_SOURCE_DIR}/src/src/C/MSIS/*.c)
file(GLOB sources_msis_h ${PROJECT_SOURCE_DIR}/src/src/C/MSIS/*.h)

file(GLOB SGP4_cpp ${PROJECT_SOURCE_DIR}/src/libsgp4/*.cc)
file(GLOB SGP4_h ${PROJECT_SOURCE_DIR}/src/libsgp4/*.h)

file(GLOB geolib_cpp ${PROJECT_SOURCE_DIR}/src/GeographicLib/*.cpp)
file(GLOB geolib_hpp ${PROJECT_SOURCE_DIR}/src/GeographicLib/*.hpp)
file(GLOB geolib_h ${PROJECT_SOURCE_DIR}/src/GeographicLib/*.h)

file(GLOB sources_UUID ${PROJECT_SOURCE_DIR}/src/src/C/UUID/*.cc ${PROJECT_SOURCE_DIR}/src/src/C/UUID/*.hpp ${PROJECT_SOURCE_DIR}/src/src/C/UUID/*.c ${PROJECT_SOURCE_DIR}/src/src/C/UUID/*.h)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries

if (UNIX AND NOT APPLE)
    add_executable(TGF_Propa ${main} ${sources_UUID} ${geolib_cpp} ${geolib_h} ${geolib_hpp} ${SGP4_cpp} ${SGP4_h} ${sources_mag_h} ${sourcesf} ${sources_msis_c} ${sources_msis_h} ${sources_mag_c} ${sourcesc} ${sourceshh} ${sources})
endif (UNIX AND NOT APPLE)

if (APPLE)
    add_executable(TGF_Propa ${main} ${sources_UUID} ${geolib_cpp} ${geolib_h} ${geolib_hpp} ${SGP4_cpp} ${SGP4_h} ${sources_msis_c} ${sources_msis_h} ${sources_mag_c} ${sources_mag_h} ${sourcesc} ${sourceshh} ${sources})
endif (APPLE)

target_link_libraries(TGF_Propa ${Geant4_LIBRARIES})

target_include_directories(TGF_Propa PUBLIC ${UUID_INCLUDE_DIRS})
target_link_libraries(TGF_Propa ${UUID_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

#if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#    # using GCC
#    target_compile_options(TGF_Propa PRIVATE -std=c++11 -mfpmath=sse -msse2 -ffp-contract=off)
#elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
#    # using Intel C++
#    target_compile_options(TGF_Propa PRIVATE -fp-model=extended)
#elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
#    # using Intel C++
#    target_compile_options(TGF_Propa PRIVATE -std=c++1y)
#endif ()

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(TGF_p DEPENDS TGF_Propa)

