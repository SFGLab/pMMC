cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
project(3dgnome LANGUAGES C CXX CUDA)

set(SOURCES
  src/main.cpp
  src/Anchor.cpp
  src/chromosomesset.cpp
  src/HierarchicalChromosome.cu
  src/LooperSolver.cpp
  src/BedRegion.cpp
  src/cluster.cpp
  src/InteractionArc.cpp
  src/ParallelMarkArcs.cu
  src/BedRegions.cpp
  src/density.cpp
  src/InteractionArcs.cpp
  src/ParallelMonteCarloHeatmap.cu
  src/ParallelMonteCarloArcs.cu
  src/ParallelMonteCarloSmooth.cu
  src/chromosome.cpp
  src/heatmap.cpp
  src/Settings.cpp
  src/SyntheticGenerator.cpp
  src/BenchmarkRunner.cpp
  src/MultiscaleEnergy.cpp
  src/SimulationLogger.cpp
  src/DistanceMapGenerator.cpp
  src/MetricsFramework.cpp
  src/CifWriter.cpp
  src/PdbWriter.cpp
  src/MdsInitializer.cpp
  thirdparty/common.cpp
  thirdparty/ini.c
  thirdparty/INIReader.cpp
  thirdparty/mtxlib.cpp
  thirdparty/rmsd.cpp
)

# Add portable getopt on Windows
if(WIN32)
  list(APPEND SOURCES thirdparty/getopt.c)
endif()

set(HEADERS
  include/Anchor.h
  include/Chromosome.h
  include/Density.h
  include/InteractionArc.h
  include/Settings.h
  include/BedRegion.h
  include/ChromosomesSet.h
  include/Heatmap.h
  include/InteractionArcs.h
  include/BedRegions.h
  include/Cluster.h
  include/HierarchicalChromosome.h
  include/LooperSolver.h
  include/SyntheticGenerator.h
  include/BenchmarkRunner.h
  include/MultiscaleEnergy.h
  include/SimulationLogger.h
  include/DistanceMapGenerator.h
  include/MetricsFramework.h
  include/platform.h
  include/CifWriter.h
  include/PdbWriter.h
  include/MdsInitializer.h
  thirdparty/common.h
  thirdparty/ini.h
  thirdparty/INIReader.h
  thirdparty/mtxlib.h
  thirdparty/rmsd.h
)

add_executable(pMMC)

# Platform-specific compiler options
if(MSVC)
  target_compile_options(pMMC PUBLIC /W3)
  target_compile_definitions(pMMC PUBLIC _CRT_SECURE_NO_WARNINGS)
else()
  target_compile_options(pMMC PUBLIC -Wall)
endif()

target_compile_features(pMMC PUBLIC cxx_std_17)
set_target_properties(pMMC PROPERTIES CXX_STANDARD 17)
set_target_properties(pMMC PROPERTIES CMAKE_CUDA_STANDARD 11)
set_target_properties(pMMC PROPERTIES CMAKE_CUDA_STANDARD_REQUIRED ON)

if(NOT CUDA_ARCH)
  message(SEND_ERROR "YOU NEED TO SPECIFY CUDA ARCHITECTURE E.G. -DCUDA_ARCH=\"80\"")
endif()

set_target_properties(pMMC PROPERTIES CUDA_ARCHITECTURES "${CUDA_ARCH}")

set_target_properties(pMMC PROPERTIES CUDA_RUNTIME_LIBRARY static)

target_include_directories(pMMC PUBLIC
  ${CMAKE_SOURCE_DIR}/thirdparty/
  ${CMAKE_SOURCE_DIR}/include/
)

target_sources(pMMC PUBLIC ${SOURCES} ${HEADERS})

set_target_properties(pMMC PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

set_property(TARGET pMMC PROPERTY CUDA_SEPARABLE_COMPILATION ON)

# Link psapi on Windows for GetProcessMemoryInfo (memory instrumentation)
if(WIN32)
  target_link_libraries(pMMC PUBLIC psapi)
endif()
