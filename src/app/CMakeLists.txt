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
	apmi_nullmodel.cpp
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
	include(FindOpenMP)
	if(OpenMP_CXX_FOUND)
		add_definitions(-DOPENMP_ENABLED=1)
		message(STATUS "OpenMP package found by include(FindOpenMP)")

		target_link_libraries(ARACNe3_app PUBLIC OpenMP::OpenMP_CXX)
	else()
		message(STATUS "OpenMP package not found by include(FindOpenMP). Attempting fallback methods.")
		# Check for the macOS platform
		if(APPLE)
			# Check for x86 architecture
			execute_process(
				COMMAND uname -m
				OUTPUT_VARIABLE ARCHITECTURE
				OUTPUT_STRIP_TRAILING_WHITESPACE
				)

			if("${ARCHITECTURE}" STREQUAL "x86_64")
				# Set the flags if on x86 macOS
				set(LIBOMP_PATH "/usr/local/opt/libomp")

				# Create libomp as an imported target
				message(STATUS "Attempting libomp_imported target creation")
				add_library(libomp_imported STATIC IMPORTED)
				set_target_properties(libomp_imported PROPERTIES
					IMPORTED_LOCATION "${LIBOMP_PATH}/lib/libomp.a"
					INTERFACE_INCLUDE_DIRECTORIES "${LIBOMP_PATH}/include"
					)
				
				message(STATUS "Attempting to link libomp_imported")
				target_link_libraries(ARACNe3_app PRIVATE libomp_imported)
			endif()
		endif()
	endif()
endif(USE_OPENMP)

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
