include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(algorithms_test test_algorithms.cpp)
target_link_libraries(algorithms_test gtest gtest_main ARACNe3_lib)

# Create a test called "AlgorithmTests" based on the executable "algorithm_tests"
add_test(NAME AlgorithmsTest COMMAND algorithms_test)
