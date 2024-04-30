# ---- Define aracne3tools target ----
add_library(aracne3tools
	cmdline_parser.cpp
	stopwatch.cpp
	algorithms.cpp
	apmi_nullmodel.cpp
  logger.cpp
  subnet_logger.cpp
  analysis.cpp
  )

# TODO: Make serialization link privately so not all libs need to link serialization
target_link_libraries(aracne3tools PUBLIC
  Boost::serialization
  )

target_include_directories(aracne3tools PRIVATE
  ${PROJECT_SOURCE_DIR}/include/ARACNe3
)

# TODO: Remove apmi_nullmodel dependency on linking fs (only fsio should)
if(UNIX AND NOT APPLE)
  target_link_libraries(aracne3tools PRIVATE stdc++fs)
endif()

# ---- Define the fsio target ----
add_library(fsio STATIC
	fsio.cpp
  )

target_include_directories(fsio PRIVATE
	${PROJECT_SOURCE_DIR}/include/ARACNe3
)

target_link_libraries(fsio PRIVATE
  aracne3tools
  )

if(UNIX AND NOT APPLE)
  target_link_libraries(fsio PRIVATE stdc++fs)
endif()

# ---- Build ARACNe3 ----
add_executable(ARACNe3_app main.cpp subnet_operations.cpp)

set_target_properties(ARACNe3_app PROPERTIES
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED ON
  DEBUG_POSTFIX "_debug"
  RELEASE_POSTFIX "_release"
)

target_link_libraries(ARACNe3_app PRIVATE 
  fsio
  aracne3tools
  )

target_include_directories(ARACNe3_app PRIVATE 
  ${PROJECT_SOURCE_DIR}/include/ARACNe3
  )

if(USE_OPENMP)
  include(FindOpenMP)
	if(OpenMP_CXX_FOUND)
		add_definitions(-DOPENMP_ENABLED=1)
		message(STATUS "OpenMP package found by include(FindOpenMP)")

		target_link_libraries(ARACNe3_app PUBLIC OpenMP::OpenMP_CXX)
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

      target_link_libraries(ARACNe3_app PUBLIC libomp_imported)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xclang -fopenmp")
    endif(APPLE)
  endif()
endif(USE_OPENMP)
