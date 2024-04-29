#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <random>
#include <string>

#include "fsio.hpp"

class FilesystemIOTest : public ::testing::Test {
protected:
  std::filesystem::path temp_dir = std::filesystem::temp_directory_path() /
                                   ("test_fsio_" + std::to_string(std::rand()));
  std::filesystem::path test_output_dir = temp_dir / "aracne3";
  std::filesystem::path test_final_output_file =
      test_output_dir / "network_testid.tsv";
  std::filesystem::path test_subnets_dir = test_output_dir / "subnetworks";

  std::filesystem::path test_exp_mat_file = temp_dir / "test_exp_mat.txt";
  std::filesystem::path test_reg_list_file = temp_dir / "test_reg_list.txt";

  FilesystemIOHandler io = FilesystemIOHandler(
      test_exp_mat_file.string(), test_reg_list_file.string(),
      test_final_output_file.string(), test_subnets_dir.string(), "testid",
      '\t');

  void SetUp() override {

    // ---- Make the test input files ----

    // make the test dir root
    std::filesystem::create_directories(temp_dir);

    // make the simulated output dir
    std::filesystem::create_directories(test_output_dir);
    std::filesystem::create_directories(test_subnets_dir);

    // Make input files
    std::ofstream ofs(test_exp_mat_file);
    ofs << std::string("gene\tsamp1\t2\t3.3\tsamp4\n") +
          "Gene1\t0.1\t0.2\t0.3\t0.4\n" +
          "5\t0.2\t0.3\t0.1\t0.15\n" +
          "Three\t0.\t0.1\t0.\t0.1";
    ofs.close();

    ofs = std::ofstream(test_reg_list_file);
    ofs << std::string("5\n") +
          "gene\n" +
          "Gene1\n" +
          "Gene2";  // "gene" and "Gene2" should not be added in memory"
    ofs.close();
  }

  void TearDown() override {
    // Remove the temporary directory and all its contents
    std::filesystem::remove_all(temp_dir);
  }
};

TEST_F(FilesystemIOTest, ReadExpMatrixAndCopulaTransform) {
  std::mt19937 rnd(123); // Use a fixed seed for reproducibility

  auto [gexp_matrix, genes, compressor, decompressor] =
      io.readExpMatrixAndCopulaTransform(rnd);

  // Check the size of the matrix and maps
  ASSERT_EQ(gexp_matrix.size(), 3);
  ASSERT_EQ(gexp_matrix[0].size(), 4);
  ASSERT_EQ(genes.size(), 3);
  ASSERT_EQ(compressor.size(), 3);
  ASSERT_EQ(decompressor.size(), 3);

  // Check whether the compressor and decompressor fully work!
  for (auto [gene, val] : compressor)
    ASSERT_EQ(decompressor[val], gene);

  for (size_t val = 0u; val < decompressor.size(); ++val)
    ASSERT_EQ(compressor[decompressor[val]], val);

  /*
   * Matrix is the following:
   * ---------------------------
   * Gene1    1/4, 2/4, 3/4, 4/4
   * 5        3/4, 4/4, 1/4, 2/4
   * Three    1/4, 3/4, 1/4, 3/4
   * ---------------------------
   *
   * Where 1/4 and 1/4 should tie each time. Ideally we would check that
   * repeats do this, but that was tested separately. We just check that
   * either is 1/4 and 2/4.
   *
   * Also leaving out not testing all values.
   */
  // row 1 test
  EXPECT_NEAR(gexp_matrix[0][0], 1. / 4, .01);
  EXPECT_NEAR(gexp_matrix[0][2], 3. / 4, .01);

  // row 2 test
  EXPECT_NEAR(gexp_matrix[1][1], 4. / 4, .01);
  EXPECT_NEAR(gexp_matrix[1][3], 2. / 4, .01);

  // row 3 test for the equal values. They should be shuffled so we test both
  auto near = [](float matrix_val, float expected_1, float expected_2,
                 float tolerance) {
    bool is_near_1 = std::abs(matrix_val - expected_1) <= tolerance;
    bool is_near_2 = std::abs(matrix_val - expected_2) <= tolerance;
    return is_near_1 || is_near_2;
  };

  float tol = 0.01f;
  EXPECT_TRUE(near(gexp_matrix[2][0], 1. / 4, 2. / 4, tol));
  EXPECT_TRUE(near(gexp_matrix[2][1], 3. / 4, 4. / 4, tol));
  EXPECT_TRUE(near(gexp_matrix[2][2], 1. / 4, 2. / 4, tol));
  EXPECT_TRUE(near(gexp_matrix[2][3], 3. / 4, 4. / 4, tol));
}

TEST_F(FilesystemIOTest, ReadRegList) {
  std::mt19937 rnd(123); // Use a fixed seed for reproducibility

  compression_map defined_genes = {{"Gene1", 1}, {"5", 42}, {"Three", 2}}; 

  auto [regs, w_list] = io.readRegList(defined_genes);

  // There were 2 defined in the compressor, 2 not defined
  ASSERT_EQ(regs.size(), 2);
  ASSERT_EQ(w_list.size(), 2);
  ASSERT_TRUE(regs.find(42) != regs.end());  // "5" was in the regulators
  ASSERT_TRUE(regs.find(0) == regs.end());  // 0 should have no corresponding
}
