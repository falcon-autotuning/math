cmake_minimum_required(VERSION 3.20)
project(
  falcon-math
  VERSION 1.0.0
  DESCRIPTION "Mathematical models and nonlinear curve fitting routines for Falcon"
  HOMEPAGE_URL "https://github.com/falcon-autotuning/math"
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(BUILD_TESTS "Build tests" ON)

# Ensure RPATH includes library directories of linked targets
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_BUILD_RPATH_USE_ORIGIN TRUE)
if(UNIX AND NOT APPLE)
  add_link_options("-Wl,--disable-new-dtags")
endif()

# Find dependencies via vcpkg
find_package(falcon-typing CONFIG REQUIRED)
find_package(falcon-core   CONFIG REQUIRED)
find_package(ensmallen     CONFIG REQUIRED)
find_package(Armadillo     CONFIG REQUIRED)

# Core C++ shared library: falcon-math
add_library(falcon-math SHARED
  src/math.cpp
)

target_include_directories(falcon-math
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(falcon-math
  PUBLIC
    falcon-core::falcon-core
    ensmallen::ensmallen
    armadillo
)

target_compile_definitions(falcon-math PRIVATE falcon_math_EXPORTS)
target_compile_options(falcon-math PRIVATE -Wall -Wextra -O3 -fPIC)

# FFI wrapper shared library: math-wrapper
add_library(math-wrapper SHARED
  math-wrapper.cpp
)

set_target_properties(math-wrapper PROPERTIES
  PREFIX ""
  OUTPUT_NAME "math-wrapper"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build"
)

target_include_directories(math-wrapper
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(math-wrapper
  PRIVATE
    falcon-math
    falcon-core::falcon-core
    falcon::falcon-typing
)

target_compile_options(math-wrapper PRIVATE -Wall -Wextra -O3 -fPIC)

# Tests
if(BUILD_TESTS)
  enable_testing()
  find_package(GTest CONFIG REQUIRED)

  add_executable(falcon_math_unit_tests
    tests/unit/test_math.cpp
  )

  target_include_directories(falcon_math_unit_tests
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/include
  )

  target_link_libraries(falcon_math_unit_tests
    PRIVATE
      falcon-math
      GTest::gtest
      GTest::gtest_main
  )

  add_test(NAME falcon_math_unit_tests COMMAND falcon_math_unit_tests)
endif()

# Install rules
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS falcon-math
  EXPORT falcon-math-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY include/falcon-math
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING
    PATTERN "*.hpp"
    PATTERN "*.h"
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/falcon-math-config-version.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)

install(EXPORT falcon-math-targets
  FILE falcon-math-targets.cmake
  NAMESPACE falcon::
  DESTINATION ${CMAKE_INSTALL_DATADIR}/falcon-math
)
