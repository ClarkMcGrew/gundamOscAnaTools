#
# GUNDAM Example: Tabulate Oscillation Weights
#

# cmake_minimum_required() should be called prior to this top-level project()
# 3.5 minimum a priori. Taking a lower version as min will make recent CMake
# version complain about the deprecation version older than 3.5.
# Might require higher version for specific features.
cmake_minimum_required( VERSION 3.5 FATAL_ERROR )

# Define the project
project( GUNDAM-NuOscillator )

message("System name: ${CMAKE_SYSTEM_NAME}")
message("System processor: ${CMAKE_SYSTEM_PROCESSOR}")

set(UseGPU 1 CACHE STRING "NuOscillator Option -- Use the GPU")
set(UseDoubles 1 CACHE STRING "NuOscillator Option -- Use doubles")
set(UseCUDAProb3 1 CACHE STRING  "NuOscillator Option -- Enable CUDAProb3")
set(UseCUDAProb3Linear 0 CACHE STRING  "NuOscillator Option -- Enable CUDAProb3Linear")
set(UseOscProb 1 CACHE STRING  "NuOscillator Option -- Enable OscProb")
set(UseProb3ppLinear 1 CACHE STRING  "NuOscillator Option -- Enable Prob3ppLinear")
set(UseNuFASTLinear 1 CACHE STRING  "NuOscillator Option -- Enable NuFASTLinear")
set(UseProbGPULinear 1 CACHE STRING  "NuOscillator Option -- Enable ProbGPULinear")
set(NuOscillatorTag "main" CACHE STRING "NuOscillator tag to build")

# Get the NuOscillator package and build it.
include(FetchContent)
FetchContent_Declare(nuoscillator
  GIT_REPOSITORY https://github.com/dbarrow257/NuOscillator.git
  GIT_TAG ${NuOscillatorTag}
)
FetchContent_MakeAvailable(nuoscillator)

link_directories(${nuoscillator_BINARY_DIR}/lib)
include_directories("${nuoscillator_SOURCE_DIR}")

cmessage(STATUS "nuOscillator is populated:     ${nuoscillator_POPULATED}")
cmessage(STATUS "nuOscillator source directory: ${nuoscillator_SOURCE_DIR}")
cmessage(STATUS "nuOscillator binary directory: ${nuoscillator_BINARY_DIR}")

# Build a shared library that can be included into GUNDAM
add_library(TabulatedNuOscillator SHARED TabulatedNuOscillator.cc)
target_link_libraries(TabulatedNuOscillator PUBLIC NuOscillator::All)
install(TARGETS TabulatedNuOscillator DESTINATION lib)

# Add the local configuration files.
install(FILES Configs/GUNDAM_NuFASTLinear DESTINATION Configs)
install(FILES Configs/GUNDAM_OscProb DESTINATION Configs)

# A test executable
add_executable(testTabulatedNuOscillator testTabulatedNuOscillator.cc)
install(TARGETS testTabulatedNuOscillator DESTINATION bin)
if (UseNuFASTLinear)
  target_compile_definitions(testTabulatedNuOscillator PUBLIC TestNUFAST)
endif(UseNuFASTLinear)
if (UseOscProb)
  target_compile_definitions(testTabulatedNuOscillator PUBLIC TestOscProb)
endif(UseOscProb)
if (UseCUDAProb3)
  target_compile_definitions(testTabulatedNuOscillator PUBLIC TestCUDAProb3)
endif(UseCUDAProb3)
if (UseOscProb)
  target_compile_definitions(testTabulatedNuOscillator PUBLIC TestOscProb)
endif(UseOscProb)
if (UseCUDAProb3)
  target_compile_definitions(testTabulatedNuOscillator PUBLIC TestCUDAProb3)
endif(UseCUDAProb3)


# An executable to draw the atmospheric oscillation probabilites
add_executable(dumpAtmOscProb dumpAtmOscProb.cc)
install(TARGETS dumpAtmOscProb DESTINATION bin)
target_link_libraries(dumpAtmOscProb PUBLIC NuOscillator::All)
if (UseOscProb)
  target_compile_definitions(dumpAtmOscProb PUBLIC TestOscProb)
endif(UseOscProb)
if (UseCUDAProb3)
  target_compile_definitions(dumpAtmOscProb PUBLIC TestCUDAProb3)
endif(UseCUDAProb3)
