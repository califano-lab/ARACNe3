# Make sure to include the relevant header files in test files. 

# In the top-level CMakeLists.txt file, the add_subdirectory(src/app) command
# should come before add_subdirectory(tests) to ensure the MINDy3_lib library
# is built before the tests are.

add_executable(tests
    test_main.cpp
    test_algorithms.cpp
    test_apcmi_nullmodel.cpp
    test_io.cpp
    test_MINDy3.cpp
)

target_link_libraries(tests
    PRIVATE
        gtest
        gtest_main
	MINDy3_lib
	OpenMP::OpenMP_CXX
	Boost::math
)

add_test(
   NAME google_test
   COMMAND tests
)
