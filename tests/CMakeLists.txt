cmake_minimum_required(VERSION 3.24)
project(AS3Tests LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build configuration" FORCE)
endif()

include(CTest)
if(NOT BUILD_TESTING)
  return()
endif()

find_package(GTest CONFIG QUIET)
if(NOT GTest_FOUND)
  include(FetchContent)
  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(googletest
    URL https://codeload.github.com/google/googletest/tar.gz/refs/tags/v1.17.0
    URL_HASH SHA256=65fab701d9829d38cb77c14acdc431d2108bfdbf8979e40eb8ae567edf10b27c
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  )
  FetchContent_MakeAvailable(googletest)
endif()

# These numerical routines use the LP64 BLAS/LAPACK interface.
set(BLA_SIZEOF_INTEGER 4)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

set(AS3_SOLVER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../solver")
add_library(as3_test_numerics STATIC
  "${AS3_SOLVER_DIR}/src/error_structure.cpp"
  "${AS3_SOLVER_DIR}/src/integration_rule_structure.cpp"
  "${AS3_SOLVER_DIR}/src/linear_algebra.cpp"
  "${AS3_SOLVER_DIR}/src/polynomial_structure.cpp"
  "${AS3_SOLVER_DIR}/src/quadrature_structure.cpp"
  "${AS3_SOLVER_DIR}/src/standard_element_structure.cpp"
  "${AS3_SOLVER_DIR}/src/tensor_structure.cpp"
)
target_include_directories(as3_test_numerics PUBLIC
  "${AS3_SOLVER_DIR}/include"
  "${AS3_SOLVER_DIR}/include/quadrature"
)
target_compile_definitions(as3_test_numerics PRIVATE DEBUG=1)
target_link_libraries(as3_test_numerics PUBLIC LAPACK::LAPACK BLAS::BLAS)

add_executable(as3_unit_tests
  unit/quadrature_test.cpp
  unit/linear_algebra_test.cpp
  unit/polynomial_test.cpp
  unit/tensor_product_test.cpp
)
target_link_libraries(as3_unit_tests PRIVATE as3_test_numerics GTest::gtest_main)
target_compile_definitions(as3_unit_tests PRIVATE DEBUG=1)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU|IntelLLVM")
  target_compile_options(as3_unit_tests PRIVATE -Wall -Wextra -Wconversion)
endif()
include(GoogleTest)
gtest_discover_tests(as3_unit_tests DISCOVERY_TIMEOUT 30)

# Preserve the existing standalone regressions while new tests use GoogleTest.
add_executable(integration_rule_tables_regression regression/integration_rule_tables_test.cpp)
target_link_libraries(integration_rule_tables_regression PRIVATE as3_test_numerics)
add_test(NAME Regression.IntegrationRuleTables COMMAND integration_rule_tables_regression)

add_executable(element_dimensions_regression regression/element_dimensions_test.cpp)
target_link_libraries(element_dimensions_regression PRIVATE as3_test_numerics)
add_test(NAME Regression.ElementDimensions COMMAND element_dimensions_regression)
