cmake_minimum_required(VERSION 3.10)

project(cerebra VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB_RECURSE PROJECT_SOURCES
    "${PROJECT_SOURCE_DIR}/src/*.c"
    "${PROJECT_SOURCE_DIR}/src/linear/cblas_linear.c"
)

if(PROJECT_SOURCES)
  add_executable(${PROJECT_NAME} ${PROJECT_SOURCES})
  target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/include)
endif()

find_package(BLAS REQUIRED)
if(TARGET ${PROJECT_NAME})
  target_link_libraries(${PROJECT_NAME} PUBLIC ${BLAS_LIBRARIES})
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/bench/CMakeLists.txt")
  add_subdirectory(bench)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/tests/CMakeLists.txt")
  enable_testing()
  add_subdirectory(tests)
endif()
