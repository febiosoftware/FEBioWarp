cmake_minimum_required(VERSION 3.15)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
mark_as_advanced(FORCE CMAKE_INSTALL_PREFIX)

project(FEBioWarp)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Set a default build type if none was specified
set(default_build_type "Release")
 
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

find_package(FEBio REQUIRED)

##### Set appropriate defines and includes #####

if(WIN32)
    add_definitions(-DWIN32 /MP /openmp)
elseif(APPLE)
    add_definitions(-D__APPLE__)
    set(CMAKE_OSX_DEPLOYMENT_TARGET 10.15)
else()
	add_definitions(-DLINUX)
    add_compile_options(-static-libstdc++ -static-libgcc -w -Wall -fopenmp)
    
    set(CMAKE_BUILD_RPATH_USE_LINK_PATH FALSE)
    set(CMAKE_BUILD_RPATH $ORIGIN/../lib/)
endif()

##### Add library #####

file(GLOB HDR_FEBioWarp "FEWarp/*.h")
file(GLOB SRC_FEBioWarp "FEWarp/*.cpp")

add_library(FEBioWarp SHARED ${HDR_FEBioWarp} ${SRC_FEBioWarp})
set_property(TARGET FEBioWarp PROPERTY AUTOGEN_BUILD_DIR ${CMAKE_BINARY_DIR}/CMakeFiles/AutoGen/FEBioWarp_autogen)

##### Link Libraries #####

target_link_libraries(FEBioWarp FEBio::FECore FEBio::FEBioMech FEBio::FEImgLib FEBio::FEBioXML)
