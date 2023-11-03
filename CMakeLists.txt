cmake_minimum_required(VERSION 3.1...3.15)

project(ARACNe3	VERSION 1.0.0
	DESCRIPTION "Algorithm for the Reconstruction of Accurate Cellular Networks 3"
		LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

message(STATUS "CMake version ${CMAKE_VERSION}")

# Set the default build type to Release if not specified.  RelWithDebInfo
# minimizes performance costs while including debugging information.
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

# #ifdef _DEBUG used for conditional compilation
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
	add_definitions(-D_DEBUG)
endif()

# Set compiler flags for Debug and Release builds
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")


##### REQUIRED LIBRARIES #####
add_subdirectory(deps/boost EXCLUDE_FROM_ALL)


##### OPTIONS #####

## Dev
option(ENABLE_ALL_WARNINGS "Enable all compiler warnings. Only effective for
GCC/Clang" OFF)

## Sanitizers
option(USE_SANITIZER "Use santizer flags" OFF)
option(SANITIZER_PATH "Path to sanitizers.")
set(ENABLED_SANITIZERS "address" "leak" CACHE STRING
	  "Semicolon separated list of sanitizer names. E.g 'address;leak'.
	  Supported sanitizers are
	  address, leak, undefined and thread.")

if(USE_SANITIZER)
	  include(cmake/Sanitizer.cmake)
	  enable_sanitizers("${ENABLED_SANITIZERS}")
endif(USE_SANITIZER)

## OpenMP
option(USE_OPENMP "Compile with multithreading support" ON)

##### ADD SUBDIRECTORIES #####

add_subdirectory(src)


##### TESTING #####
option(BUILD_TESTS "Build the tests" OFF)

if(BUILD_TESTS)
    add_subdirectory(${PROJECT_SOURCE_DIR}/deps/gtest)
    enable_testing()
    add_subdirectory(tests)
endif()
