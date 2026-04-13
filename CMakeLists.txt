cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
project(pMMC LANGUAGES C CXX CUDA)

# ── CUDA architecture ───────────────────────────────────────────────
if(NOT CUDA_ARCH)
  message(SEND_ERROR "YOU NEED TO SPECIFY CUDA ARCHITECTURE E.G. -DCUDA_ARCH=\"89\"")
endif()

# ── pMMC-lib  (static library: CUDA kernels + thirdparty) ──────────
set(LIB_CUDA_SOURCES
  src/HierarchicalChromosome.cu
  src/ParallelMarkArcs.cu
  src/ParallelMonteCarloHeatmap.cu
  src/ParallelMonteCarloArcs.cu
  src/ParallelMonteCarloSmooth.cu
)

set(LIB_CPP_SOURCES
  src/thirdparty/common.cpp
  src/thirdparty/ini.c
  src/thirdparty/INIReader.cpp
  src/thirdparty/mtxlib.cpp
  src/thirdparty/rmsd.cpp
)

if(WIN32)
  list(APPEND LIB_CPP_SOURCES src/thirdparty/getopt.c)
endif()

add_library(pMMC-lib STATIC ${LIB_CUDA_SOURCES} ${LIB_CPP_SOURCES})

target_include_directories(pMMC-lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/include
  ${CMAKE_SOURCE_DIR}/src/thirdparty
)

target_compile_features(pMMC-lib PUBLIC cxx_std_17)

if(MSVC)
  target_compile_options(pMMC-lib PUBLIC $<$<COMPILE_LANGUAGE:CXX>:/W3> $<$<COMPILE_LANGUAGE:C>:/W3>)
  target_compile_definitions(pMMC-lib PUBLIC _CRT_SECURE_NO_WARNINGS WIN32 _WIN64)
else()
  target_compile_options(pMMC-lib PUBLIC $<$<COMPILE_LANGUAGE:CXX>:-Wall> $<$<COMPILE_LANGUAGE:C>:-Wall>)
endif()

set_target_properties(pMMC-lib PROPERTIES
  CUDA_ARCHITECTURES "${CUDA_ARCH}"
  CUDA_SEPARABLE_COMPILATION ON
  CUDA_RESOLVE_DEVICE_SYMBOLS ON
  CUDA_RUNTIME_LIBRARY Static
)

# OpenMP for CPU-side parallel execution (CXX only, not CUDA)
include(CheckCXXCompilerFlag)
if(MSVC)
  check_cxx_compiler_flag(/openmp HAS_OPENMP)
  if(HAS_OPENMP)
    target_compile_options(pMMC-lib PUBLIC $<$<COMPILE_LANGUAGE:CXX>:/openmp>)
    message(STATUS "OpenMP found — CPU multi-threading enabled")
  else()
    message(STATUS "OpenMP not found — CPU multi-threading disabled")
  endif()
else()
  find_package(OpenMP COMPONENTS CXX)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(pMMC-lib PUBLIC OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP found — CPU multi-threading enabled")
  else()
    message(STATUS "OpenMP not found — CPU multi-threading disabled")
  endif()
endif()

# ── Helper: add an executable that links pMMC-lib ──────────────────
function(add_pmmc_exe TARGET_NAME SOURCE_FILE)
  add_executable(${TARGET_NAME} ${SOURCE_FILE})
  target_link_libraries(${TARGET_NAME} PRIVATE pMMC-lib)
  if(WIN32)
    target_link_libraries(${TARGET_NAME} PRIVATE psapi)
  endif()
  if(MSVC)
    set_target_properties(${TARGET_NAME} PROPERTIES LINK_FLAGS "/NODEFAULTLIB:LIBCMT")
  endif()
  set_target_properties(${TARGET_NAME} PROPERTIES
    CUDA_ARCHITECTURES "${CUDA_ARCH}"
    CUDA_RUNTIME_LIBRARY Static
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
  )
endfunction()

# ── Executables ────────────────────────────────────────────────────
add_pmmc_exe(pMMC-reconstruct  src/apps/pMMC-reconstruct/main_reconstruct.cpp)
add_pmmc_exe(pMMC-analyze      src/apps/pMMC-analyze/main_analyze.cpp)
add_pmmc_exe(pMMC-convert      src/apps/pMMC-convert/main_convert.cpp)
add_pmmc_exe(pMMC-generate     src/apps/pMMC-generate/main_generate.cpp)
add_pmmc_exe(pMMC-benchmark    src/apps/pMMC-benchmark/main_benchmark.cpp)
add_pmmc_exe(pMMC-flowchart    src/apps/pMMC-flowchart/main_flowchart.cpp)
