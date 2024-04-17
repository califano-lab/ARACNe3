# ---- Project Definitions ----
cmake_minimum_required(VERSION 3.1...3.15)

project(ARACNe3	VERSION 0.2.0
	DESCRIPTION "Algorithm for the Reconstruction of Accurate Cellular Networks"
		LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)  # Must be copied later

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

message(STATUS "CMake version ${CMAKE_VERSION}")

# Set the default build type if not specified with -DCMAKE_BUILD_TYPE=<type>
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Set compiler flags for Debug and Release builds
if(MSVC)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Zi /W4")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Ox")
else()
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -Wall -Wextra -Wpedantic")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
endif()

# ifdef _DEBUG for conditional compilation
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
	add_definitions(-D_DEBUG)
endif()

if(WIN32)
    set(CACHE_PATH "$ENV{LOCALAPPDATA}/aracne3")
elseif(APPLE)
    set(CACHE_PATH "$ENV{HOME}/Library/Caches/aracne3")
else() # Linux/Unix
    set(CACHE_PATH "$ENV{HOME}/.cache/aracne3")
endif()

# ---- Generate Config.h from template header ----
# Config.h will pass the caching path and project version
configure_file(
  "${PROJECT_SOURCE_DIR}/include/ARACNe3/generated/config.h.in"
  "${PROJECT_BINARY_DIR}/include/ARACNe3/config.h"
  )

include_directories("${PROJECT_BINARY_DIR}/include/ARACNe3")

# ---- Required Libraries ----
add_subdirectory(deps/boost EXCLUDE_FROM_ALL)
add_subdirectory(deps/eigen EXCLUDE_FROM_ALL)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# ---- Options ----
option(USE_OPENMP "Compile with multithreading support" ON)

# ---- Add Subdirectories ----
add_subdirectory(src)

# ---- Add testing suite ----
# 'gtest' exclude from all build configurations to prevent double linking
add_subdirectory(deps/gtest EXCLUDE_FROM_ALL)

enable_testing()
include(GoogleTest)  # enable ctest to discover tests from test targets

add_subdirectory(tests)

# ---- Post-configuration actions ----

# copy compile commands
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                        "${CMAKE_BINARY_DIR}/compile_commands.json"
                        "${CMAKE_SOURCE_DIR}/compile_commands.json")
