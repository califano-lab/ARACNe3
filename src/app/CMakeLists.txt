# Set the source files for ARACNe3
set(SOURCES
	cmdline_parser.cpp
	stopwatch.cpp
	io.cpp
	algorithms.cpp
	apmi_nullmodel.cpp
	subnet_operations.cpp
)

# Mainly for testing suite, but also so ARACNe3_app can easily add includes
add_library(ARACNe3_lib ${SOURCES})

# Public includes so ARACNe3_app can see them, yet sources are private
target_include_directories(ARACNe3_lib
	PUBLIC
        ${PROJECT_SOURCE_DIR}/include/ARACNe3
	${OpenMP_CSS_INCLUDE_DIRS}
	"${CMAKE_SOURCE_DIR}/deps/boost/libs/math/include/boost/math"
)

if(USE_OPENMP)
  include(FindOpenMP)
	if(OpenMP_CXX_FOUND)
		add_definitions(-DOPENMP_ENABLED=1)
		message(STATUS "OpenMP package found by include(FindOpenMP)")

		target_link_libraries(ARACNe3_lib PUBLIC OpenMP::OpenMP_CXX)
	else()
		message(STATUS "OpenMP package not found by include(FindOpenMP).")
		# Check for the macOS platform
		if(APPLE)
      message(STATUS "Attempting APPLE fallback method.")
			# Check for x86 architecture
			execute_process(
				COMMAND uname -m
				OUTPUT_VARIABLE ARCHITECTURE
				OUTPUT_STRIP_TRAILING_WHITESPACE
				)

			if("${ARCHITECTURE}" STREQUAL "x86_64")
				# Set the flags if on x86 macOS
				set(LIBOMP_PATH "/usr/local/opt/libomp")
      else()
				# Set the flags if on arm64 macOS
				set(LIBOMP_PATH "/opt/homebrew/opt/libomp")
      endif()

      # Create libomp as an imported target
      message(STATUS "Attempting libomp_imported target creation")

      add_library(libomp_imported STATIC IMPORTED)
      set_target_properties(libomp_imported PROPERTIES
        IMPORTED_LOCATION "${LIBOMP_PATH}/lib/libomp.a"
        INTERFACE_INCLUDE_DIRECTORIES "${LIBOMP_PATH}/include"
        )
      
      message(STATUS "Attempting to link libomp_imported")

      target_link_libraries(ARACNe3_lib PUBLIC libomp_imported)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xclang -fopenmp")
    endif(APPLE)
  endif()
endif(USE_OPENMP)

target_link_libraries(ARACNe3_lib PUBLIC Boost::math)

if (WIN32 AND MINGW)
	target_compile_options(ARACNe3_lib PUBLIC -static-libstdc++)
endif (WIN32 AND MINGW)

# The -lstdc++fs flag is needed for the std::filesystem library on some platforms, like CentOS 6.5 with older GCC versions. However, on most systems, the std::filesystem library is already included in the standard library.
if(UNIX AND NOT APPLE)
	target_compile_options(ARACNe3_lib PUBLIC -lstdc++fs)
	target_link_libraries(ARACNe3_lib PUBLIC stdc++fs)
endif()

add_executable(ARACNe3_app ARACNe3.cpp)

set_target_properties(ARACNe3_app PROPERTIES
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED ON
  DEBUG_POSTFIX "_debug"
  RELEASE_POSTFIX "_release"
)

# Links the public libraries and includes, including the headers for A3
target_link_libraries(ARACNe3_app PRIVATE ARACNe3_lib)
