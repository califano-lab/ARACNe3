# Include header files in root/include/ARACNe3 directory
include_directories(${PROJECT_SOURCE_DIR}/include/ARACNe3)
include_directories(${OpenMP_CXX_INCLUDE_DIRS})

# Set the source files for ARACNe3
set(SOURCES
	ARACNe3.cpp 
	cmdline_parser.cpp
	stopwatch.cpp
	io.cpp
	algorithms.cpp
	NullModel.cpp
	subnet_operations.cpp
)

# Name of executable & target, source file list (h files ignored)
add_executable(ARACNe3_app ${SOURCES})

set_target_properties(ARACNe3_app PROPERTIES
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED ON
  DEBUG_POSTFIX "_debug"
  RELEASE_POSTFIX "_release"
)

if (WIN32 AND MINGW)
	target_compile_options(ARACNe3_app PUBLIC -static-libstdc++)
endif (WIN32 AND MINGW)

if(USE_OPENMP)
	target_link_libraries(ARACNe3_app PUBLIC OpenMP::OpenMP_CXX)
endif()

target_include_directories(ARACNe3_app PUBLIC "${CMAKE_SOURCE_DIR}/deps/boost/libs/math/include/boost/math")
target_link_libraries(ARACNe3_app PUBLIC Boost::math)

# The -lstdc++fs flag is needed for the std::filesystem library on some platforms, like CentOS 6.5 with older GCC versions. However, on most systems, the std::filesystem library is already included in the standard library.
if(UNIX AND NOT APPLE)
  target_compile_options(ARACNe3_app PRIVATE -lstdc++fs)
  target_link_libraries(ARACNe3_app PRIVATE stdc++fs)
endif()

# For testing suite
add_library(ARACNe3_lib ${SOURCES})

target_include_directories(ARACNe3_lib
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include/ARACNe3
)
