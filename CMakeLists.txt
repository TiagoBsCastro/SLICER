cmake_minimum_required(VERSION 3.10)
project(SLICER VERSION 1.0 DESCRIPTION "SLICER - Simulation LIght conE buildeR")

# Set the CMake module path to include your "cmake" directory
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler and linker options
set(CMAKE_CXX_COMPILER mpic++)
add_compile_options(-O3 -g -fPIC -Wall -Wpedantic)

# Option for enabling custom features
#option(ENABLE_CUSTOM_FEATURES "Enable custom features in SLICER" ON)

# Check if the custom features are enabled and add definitions accordingly
if(USE_REPLICATION)
	add_definitions(-DReplicationOnPerpendicularPlane)
endif()

if(USE_FIXED_PLC_VERTEX)
	add_definitions(-DFixedPLCVertex)
endif()


# Define executable output path
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

# Add the subdirectory containing the source files
add_subdirectory(SLICER)
