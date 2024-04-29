# ---- Define the test for core app functionality ----
add_executable(aracne3_test
  test_algorithms.cpp
  test_apmi_nullmodel.cpp
)

target_link_libraries(aracne3_test
  PRIVATE
    gtest_main
    aracne3tools
    Eigen3::Eigen
    Boost::random
    Boost::math
)

target_include_directories(aracne3_test
  PRIVATE
  ${PROJECT_SOURCE_DIR}/include/ARACNe3
  )

gtest_discover_tests(aracne3_test)

# ---- Define the test for filesystem IO ----
add_executable(fsio_test
  test_fsio.cpp
)

target_link_libraries(fsio_test
  PRIVATE
  gtest_main
  fsio
)

target_include_directories(fsio_test
  PRIVATE
  ${PROJECT_SOURCE_DIR}/include/ARACNe3
  )

gtest_discover_tests(fsio_test)
