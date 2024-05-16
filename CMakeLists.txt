#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.6 FATAL_ERROR)
project(TGF_propa LANGUAGES C CXX Fortran)

# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

# Define project-wide C++ standard and apply compiler flags
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#----------------------------------------------------------------------------
# Setup module paths and include directories
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/src/cmake/modules/") # Path to custom .cmake files to find packages

# Configure packages based on UI and Visual drivers option
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
find_package(Geant4 REQUIRED COMPONENTS ui_all vis_all)
find_package(PkgConfig REQUIRED)
find_package(Threads REQUIRED)
pkg_search_module(UUID REQUIRED uuid)

#----------------------------------------------------------------------------
# Include directories
include_directories(
    ${PROJECT_SOURCE_DIR}/src/libsgp4 
    ${PROJECT_SOURCE_DIR}/src/src/C/UUID/
    ${PROJECT_SOURCE_DIR}/src 
    ${PROJECT_SOURCE_DIR}/src/include 
    ${Geant4_INCLUDE_DIR}/src
    ${Geant4_INCLUDE_DIR}/src/src
    ${PROJECT_SOURCE_DIR}/src/src/C/WMM/src/
    ${PROJECT_SOURCE_DIR}/src/src/C/MSIS/
)

#----------------------------------------------------------------------------
# Locate sources
file(GLOB_RECURSE SOURCES
     "${PROJECT_SOURCE_DIR}/src/*.cc"
     "${PROJECT_SOURCE_DIR}/src/src/C/*.c"
     "${PROJECT_SOURCE_DIR}/src/src/C/WMM/src/*.c"
     "${PROJECT_SOURCE_DIR}/src/src/C/MSIS/*.c"
     "${PROJECT_SOURCE_DIR}/src/GeographicLib/*.cpp"
)

if (UNIX AND NOT APPLE)
    file(GLOB_RECURSE FORTRAN_SOURCES "${PROJECT_SOURCE_DIR}/src/src/fortran/*.f")
    list(APPEND SOURCES ${FORTRAN_SOURCES})
endif()

#----------------------------------------------------------------------------
# Add the executable
add_executable(TGF_Propa ${SOURCES})

# Link libraries
target_link_libraries(TGF_Propa PUBLIC ${Geant4_LIBRARIES} ${UUID_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

# Include and link directories
target_include_directories(TGF_Propa PUBLIC ${UUID_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Add custom target
add_custom_target(TGF_p DEPENDS TGF_Propa)
