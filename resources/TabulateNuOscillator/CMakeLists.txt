#
# GUNDAM Example: Tabulate Oscillation Weights
#

# cmake_minimum_required() should be called prior to this top-level project()
# 3.5 minimum a priori. Taking a lower version as min will make recent CMake
# version complain about the deprecation version older than 3.5.
# Might require higher version for specific features.
cmake_minimum_required( VERSION 3.12 FATAL_ERROR )

# Define the project
project( GUNDAM-NuOscillator )

message("System name: ${CMAKE_SYSTEM_NAME}")
message("System processor: ${CMAKE_SYSTEM_PROCESSOR}")

include( CheckLanguage )
check_language( CUDA )
if( CMAKE_CUDA_COMPILER )
  message( STATUS "CUDA support enabled." )
  enable_language(CUDA)
  if( NOT DEFINED CMAKE_CUDA_ARCHITECTURES )
    # The default is taken from the CUDAARCHS environment
    # variable.  If it isn't set, then set it to the earliest
    # non-deprecated architecture.
    #   2022: architectures before 52 are deprecated.
    if( ${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.23 )
      # After cmake 3.23, this can be set to all or all-major
      set( CMAKE_CUDA_ARCHITECTURES all )
    else()
      set( CMAKE_CUDA_ARCHITECTURES 52 )
    endif()
  endif()
  message( STATUS "CUDA compilation architectures: \"${CMAKE_CUDA_ARCHITECTURES}\"")
  set(UseGPU 1 CACHE STRING "NuOscillator Option -- Use the GPU")
else( CMAKE_CUDA_COMPILER )
  message( WARNING "CUDA support is not present." )
  set(UseGPU 0 CACHE STRING "NuOscillator Option -- Use the GPU")
endif( CMAKE_CUDA_COMPILER )

set(UseDoubles 1 CACHE STRING "NuOscillator Option -- Use doubles")
set(UseCUDAProb3 1 CACHE STRING  "NuOscillator Option -- Enable CUDAProb3")
set(UseCUDAProb3Linear 0 CACHE STRING  "NuOscillator Option -- Enable CUDAProb3Linear")
set(UseProb3ppLinear 1 CACHE STRING "NuOscillator Option -- Enable Prob3ppLinear")
set(UseProbGPULinear 1 CACHE STRING  "NuOscillator Option -- Enable ProbGPULinear")
set(UseNuFASTLinear 1 CACHE STRING  "NuOscillator Option -- Enable NuFASTLinear")
set(UseOscProb 1 CACHE STRING  "NuOscillator Option -- Enable OscProb")

set(NuOscillatorTag "v1.3.1" CACHE STRING "NuOscillator tag to build")

# Get the NuOscillator package and build it.
include(FetchContent)
FetchContent_Declare(nuoscillator
  GIT_REPOSITORY https://github.com/dbarrow257/NuOscillator.git
  GIT_TAG ${NuOscillatorTag}
)
FetchContent_MakeAvailable(nuoscillator)

link_directories(${nuoscillator_BINARY_DIR}/lib)
include_directories("${nuoscillator_SOURCE_DIR}")

message(STATUS "nuOscillator is populated:     ${nuoscillator_POPULATED}")
message(STATUS "nuOscillator source directory: ${nuoscillator_SOURCE_DIR}")
message(STATUS "nuOscillator binary directory: ${nuoscillator_BINARY_DIR}")

# Build a shared library that can be included into GUNDAM
add_library(TabulatedNuOscillator SHARED TabulatedNuOscillator.cc)
target_link_libraries(TabulatedNuOscillator PUBLIC NuOscillator::All)
target_link_libraries(TabulatedNuOscillator PRIVATE ${CMAKE_DL_LIBS})
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

# An executable to draw the long baseline probabilities
add_executable(dumpLblOscProb dumpLblOscProb.cc)
install(TARGETS dumpLblOscProb DESTINATION bin)
target_link_libraries(dumpLblOscProb PUBLIC NuOscillator::All)
