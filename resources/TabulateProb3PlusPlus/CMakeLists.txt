#
# GUNDAM Example: Tabulate Oscillation Weights
#

# cmake_minimum_required() should be called prior to this top-level project()
# 3.5 minimum a priori. Taking a lower version as min will make recent CMake
# version complain about the deprecation version older than 3.5.
# Might require higher version for specific features.
cmake_minimum_required( VERSION 3.5 FATAL_ERROR )

# Define the project
project( GUNDAM-Example )

message("System name: ${CMAKE_SYSTEM_NAME}")
message("System processor: ${CMAKE_SYSTEM_PROCESSOR}")

# Get the Prob3plusplus package and build it.
include(ExternalProject)
ExternalProject_add(Prob3plusplus
  GIT_REPOSITORY https://github.com/rogerwendell/Prob3plusplus.git
  GIT_TAG main)

ExternalProject_Get_Property(Prob3plusplus BINARY_DIR)
set(Prob3plusplus_DIR "${BINARY_DIR}/${CMAKE_SYSTEM_NAME}")

link_directories(${Prob3plusplus_DIR}/lib)
include_directories(${Prob3plusplus_DIR}/include)

# Build a shared library that can be included into GUNDAM
add_library(TabulatedProb3PlusPlus SHARED TabulatedProb3PlusPlus.cc)
target_link_libraries(TabulatedProb3PlusPlus PUBLIC libProb3plusplus.so)

# A test executable
add_executable(testTabulatedProb3PlusPlus testTabulatedProb3PlusPlus.cc)
add_dependencies(TabulatedProb3PlusPlus Prob3plusplus)

# A command that works to compile gundam from any sub-directory.  You
# can get this aliased to "gundam-build" by sourcing
#
# . "$(git rev-parse --show-toplevel)/cmake/scripts/gundam-setup.sh"
#
# That runs:
# Local Variables:
# compile-command:"$(git rev-parse --show-toplevel)/cmake/gundam-build.sh"
# End:
