# Include header files in root/include/MINDy3 directory
include_directories(${PROJECT_SOURCE_DIR}/include/MINDy3)
include_directories(${OpenMP_CXX_INCLUDE_DIRS})

# Set the source files for MINDy3
set(SOURCES
	MINDy3.cpp 
	algorithms.cpp 
	io.cpp
	apcmi_nullmodel.cpp
)

# Name of executable & target, source file list (h files ignored)
add_executable(MINDy3_app ${SOURCES})

set_target_properties(MINDy3_app PROPERTIES
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED ON
  DEBUG_POSTFIX "_debug"
  RELEASE_POSTFIX "_release"
)

if (WIN32 AND MINGW)
	target_compile_options(MINDy3_app PUBLIC -static-libstdc++)
endif (WIN32 AND MINGW)

if(USE_OPENMP)
	target_link_libraries(MINDy3_app PUBLIC OpenMP::OpenMP_CXX)
endif()

target_include_directories(MINDy3_app PUBLIC "${CMAKE_SOURCE_DIR}/deps/boost/libs/math/include/boost/math")
target_link_libraries(MINDy3_app PUBLIC Boost::math)

# The -lstdc++fs flag is needed for the std::filesystem library on some platforms, like CentOS 6.5 with older GCC versions. However, on most systems, the std::filesystem library is already included in the standard library.
if(UNIX AND NOT APPLE)
  target_compile_options(MINDy3_app PRIVATE -lstdc++fs)
  target_link_libraries(MINDy3_app PRIVATE stdc++fs)
endif()

add_library(MINDy3_lib ${SOURCES})

target_include_directories(MINDy3_lib
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include/MINDy3
)
