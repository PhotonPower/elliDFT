cmake_minimum_required(VERSION 3.20)
project(elliDFT VERSION 0.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(ELLIDFT_BUILD_TESTS "Build unit and regression tests" ON)
option(ELLIDFT_WERROR "Treat compiler warnings as errors" OFF)

# JSON input/output only (approved dependency, header-only). No numerics.
find_package(nlohmann_json 3.10 REQUIRED)

# Linear algebra (default of the specification). Used from the basis module on.
find_package(Eigen3 REQUIRED NO_MODULE)

# Density functionals (default of the specification). Only src/xc/libxc_evaluator.cpp
# includes its header; everything else goes through XcEvaluator.
find_package(Libxc REQUIRED CONFIG)

# Threads (specification: OpenMP from M3 on). Apple clang needs the runtime
# from Homebrew (libomp); without OpenMP the code runs in one thread.
option(ELLIDFT_OPENMP "Use OpenMP threads" ON)
if(ELLIDFT_OPENMP)
  if(APPLE AND NOT DEFINED OpenMP_ROOT)
    execute_process(COMMAND brew --prefix libomp OUTPUT_VARIABLE _ellidft_libomp
                    OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(_ellidft_libomp)
      set(OpenMP_ROOT "${_ellidft_libomp}")
    endif()
  endif()
  find_package(OpenMP COMPONENTS CXX)
  if(NOT OpenMP_CXX_FOUND)
    message(WARNING "OpenMP not found: elliDFT will run in one thread")
  endif()
endif()

# External BLAS for Eigen's matrix products (OpenBLAS, MKL, Apple Accelerate).
# Off by default: it is to be switched on only after measurements on larger
# molecules show a gain (decision of the client, 2026-11-14).
option(ELLIDFT_USE_BLAS "Let Eigen call an external BLAS" OFF)
if(ELLIDFT_USE_BLAS)
  find_package(BLAS REQUIRED)
endif()

add_library(ellidft_warnings INTERFACE)
target_compile_options(ellidft_warnings INTERFACE -Wall -Wextra -Wpedantic -Wshadow)
if(ELLIDFT_WERROR)
  target_compile_options(ellidft_warnings INTERFACE -Werror)
endif()

# Module order follows the architecture: each module depends only on the
# modules listed before it (core, basis, integrals, grid, xc, df, scf,
# response, realtime, optics, gpu). The module ellipsometry (evaluation of measurements)
# uses only optics.
add_library(ellidft STATIC
  src/core/elements.cpp
  src/core/molecule.cpp
  src/core/parallel.cpp
  src/core/xyz.cpp
  src/basis/basis_library.cpp
  src/basis/basis_set.cpp
  src/basis/solid_harmonics.cpp
  src/grid/basis_on_grid.cpp
  src/grid/lebedev.cpp
  src/grid/molecular_grid.cpp
  src/integrals/boys.cpp
  src/integrals/hermite.cpp
  src/integrals/one_electron.cpp
  src/integrals/shell_pair.cpp
  src/integrals/two_electron.cpp
  src/xc/functional.cpp
  src/xc/libxc_evaluator.cpp
  src/xc/own_lda.cpp
  src/df/df_integrals.cpp
  src/df/df_jk.cpp
  src/scf/diis.cpp
  src/scf/jk.cpp
  src/scf/rhf.cpp
  src/scf/xc_matrix.cpp
  src/response/excitations.cpp
  src/response/polarizability.cpp
  src/response/response_operator.cpp
  src/realtime/nonlinear.cpp
  src/realtime/propagator.cpp
  src/realtime/spectrum.cpp
  src/optics/optics.cpp
  src/ellipsometry/analysis.cpp
  src/ellipsometry/dispersion_fit.cpp
  src/ellipsometry/least_squares.cpp
  src/ellipsometry/materials.cpp
  src/ellipsometry/stack.cpp
)
target_include_directories(ellidft PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(ellidft PUBLIC Eigen3::Eigen PRIVATE Libxc::xc ellidft_warnings)
if(ELLIDFT_USE_BLAS)
  target_compile_definitions(ellidft PUBLIC EIGEN_USE_BLAS)
  target_link_libraries(ellidft PUBLIC BLAS::BLAS)
endif()
if(ELLIDFT_OPENMP AND OpenMP_CXX_FOUND)
  target_link_libraries(ellidft PRIVATE OpenMP::OpenMP_CXX)
else()
  target_compile_options(ellidft PRIVATE -Wno-unknown-pragmas)
endif()
# Default location of the basis set files; ELLIDFT_BASIS_DIR overrides at run time.
target_compile_definitions(ellidft PRIVATE
  ELLIDFT_DEFAULT_BASIS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data/basis")

# Main program: control file in, results and protocol out.
execute_process(COMMAND git rev-parse --short HEAD WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                OUTPUT_VARIABLE ELLIDFT_GIT_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(NOT ELLIDFT_GIT_COMMIT)
  set(ELLIDFT_GIT_COMMIT "unknown")
endif()
add_executable(ellidft-run apps/ellidft.cpp)
target_link_libraries(ellidft-run PRIVATE ellidft nlohmann_json::nlohmann_json ellidft_warnings)
set_target_properties(ellidft-run PROPERTIES OUTPUT_NAME ellidft)
target_compile_definitions(ellidft-run PRIVATE
  ELLIDFT_VERSION="${PROJECT_VERSION}" ELLIDFT_GIT_COMMIT="${ELLIDFT_GIT_COMMIT}"
  ELLIDFT_BUILD_TYPE="${CMAKE_BUILD_TYPE}" ELLIDFT_COMPILER="${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

add_executable(ellidft-basisinfo apps/basisinfo.cpp)
target_link_libraries(ellidft-basisinfo PRIVATE ellidft nlohmann_json::nlohmann_json ellidft_warnings)

add_executable(ellidft-scf apps/scf.cpp)
target_link_libraries(ellidft-scf PRIVATE ellidft nlohmann_json::nlohmann_json ellidft_warnings)

add_executable(ellidft-response apps/response.cpp)
target_link_libraries(ellidft-response PRIVATE ellidft nlohmann_json::nlohmann_json ellidft_warnings)

add_executable(ellidft-realtime apps/realtime.cpp)
target_link_libraries(ellidft-realtime PRIVATE ellidft nlohmann_json::nlohmann_json ellidft_warnings)

add_executable(ellidft-nonlinear apps/nonlinear.cpp)
target_link_libraries(ellidft-nonlinear PRIVATE ellidft nlohmann_json::nlohmann_json ellidft_warnings)

add_executable(ellidft-kerr-film apps/kerr_film.cpp)
target_link_libraries(ellidft-kerr-film PRIVATE ellidft nlohmann_json::nlohmann_json ellidft_warnings)

add_executable(ellidft-optics apps/optics.cpp)
target_link_libraries(ellidft-optics PRIVATE ellidft nlohmann_json::nlohmann_json ellidft_warnings)

add_executable(ellidft-ellipsometry apps/ellipsometry.cpp)
target_link_libraries(ellidft-ellipsometry PRIVATE ellidft nlohmann_json::nlohmann_json ellidft_warnings)

if(ELLIDFT_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
