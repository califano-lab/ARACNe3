# ---- Define the test for core app functionality ----
add_executable(aracne3tools_test
  algorithm_test.cpp
  apmi_nullmodel_test.cpp
)

target_link_libraries(aracne3tools_test
  PRIVATE
    gtest_main
    aracne3tools
    Eigen3::Eigen
    Boost::random
    Boost::math
)

target_include_directories(aracne3tools_test
  PRIVATE
  ${PROJECT_SOURCE_DIR}/include/ARACNe3
  )

gtest_discover_tests(aracne3tools_test)

# ---- Define the test for filesystem IO ----
add_executable(fsio_test
  fsio_test.cpp
)

target_link_libraries(fsio_test
  PRIVATE
  gtest_main
  aracne3tools  # for algorithm.cpp
  aracne3analysis  # for subnet logger tools
  fsio
)

target_include_directories(fsio_test
  PRIVATE
  ${PROJECT_SOURCE_DIR}/include/ARACNe3
  )

gtest_discover_tests(fsio_test)

# ---- Define the test for APMI null model ----
add_executable(apmi_nullmodel_test
  apmi_nullmodel_test.cpp
)

target_link_libraries(apmi_nullmodel_test
  PRIVATE
  gtest_main
  aracne3tools
)

target_include_directories(apmi_nullmodel_test
  PRIVATE
  ${PROJECT_SOURCE_DIR}/include/ARACNe3
  )

gtest_discover_tests(apmi_nullmodel_test)
