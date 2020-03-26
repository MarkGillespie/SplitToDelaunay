cmake_minimum_required(VERSION 2.8.9)

# Maybe stop from CMAKEing in the wrong place
if (CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(FATAL_ERROR "Source and build directories cannot be the same. Go use the /build directory.")
endif()

add_library(
  my_assert
  my_assert.cpp
  my_assert.h
  )

add_library(
  splitter
  delaunay.h
  delaunay.cpp
  )
target_link_libraries(splitter geometry-central polyscope my_assert)

add_executable(
  run
  main.cpp
)
target_include_directories(run PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../include")
target_link_libraries(run polyscope geometry-central splitter)

option(PACKAGE_TESTS "Build the tests" ON)

if(PACKAGE_TESTS)
  enable_testing()
  add_subdirectory(
    tests
  )
endif()
