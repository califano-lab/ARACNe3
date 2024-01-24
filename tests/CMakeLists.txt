add_executable(aracne3_test
  test_algorithms.cpp
)

target_link_libraries(aracne3_test
  PRIVATE
    gtest_main
    aracne3tools
)

target_include_directories(aracne3_test
  PRIVATE
  ${PROJECT_SOURCE_DIR}/include/ARACNe3
  )

gtest_discover_tests(aracne3_test)
