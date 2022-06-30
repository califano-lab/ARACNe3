cmake_minimum_required(VERSION 3.5)

if (MSVC)
	add_compile_options(/W4)
else()
	add_compile_options(-Wall -Wextra -Wpedantic)
endif()

add_executable(ARACNe3 ARACNe3.cpp NullModel.cpp MatrixReglistIO.cpp APMI.cpp AlphaPruning.cpp MaxEntPruning.cpp Consolidator.cpp RegWebFns.cpp)
target_compile_features(ARACNe3 gnu++20)
