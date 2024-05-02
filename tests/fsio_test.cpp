#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <random>
#include <string>

#include "algorithm.hpp"
#include "fsio.hpp"
#include "subnet_logger.hpp" // much easier for generating pseudo logs

class FilesystemIOTest : public ::testing::Test {
protected:
  // ---- Define all expected resources for this runtime ----

  std::filesystem::path temp_dir = std::filesystem::temp_directory_path() /
                                   ("fsio_test_" + std::to_string(std::rand()));
  std::filesystem::path test_output_dir = temp_dir / "aracne3/";
  std::filesystem::path test_final_output_file =
      test_output_dir / "network_testid.tsv";
  std::filesystem::path test_exp_mat_file = temp_dir / "test_exp_mat.txt";
  std::filesystem::path test_reg_list_file = temp_dir / "test_reg_list.txt";

  std::filesystem::path test_subnets_dir = test_output_dir / "subnetworks/";
  std::filesystem::path test_subnets_log_dir =
      test_output_dir / "log-subnetworks/";

  std::filesystem::path test_subnet_file_1 = "subnetwork-1_testid.tsv";
  std::filesystem::path test_subnet_file_2 = "subnetwork-5_testid.tsv";

  std::filesystem::path test_subnet_log_file_1 = "log-subnetwork-1_testid.txt";
  std::filesystem::path test_subnet_log_file_2 = "log-subnetwork-5_testid.txt";

  compression_map defined_genes = {{"Gene1", 1}, {"5", 42}, {"Three", 2}};
  decompression_map decompressor;
  geneset regulators = {defined_genes["5"], defined_genes["Gene1"]};

  uint32_t tot_possible_edges = 4u, num_edges_after_threshold_pruning_1 = 4u,
           num_edges_after_threshold_pruning_2 = 2u,
           num_edges_after_MaxEnt_pruning_2 = 1u;

  FilesystemIOHandler io = FilesystemIOHandler(
      test_exp_mat_file.string(), test_reg_list_file.string(),
      test_final_output_file.string(), test_subnets_dir.string(), "testid",
      '\t');

  void SetUp() override {
    std::filesystem::remove_all(temp_dir); // sometimes artifacts from testing

    // --- Set up decompressor ----
    // Your decompressor is a niche case

    decompressor = std::vector<std::string>(43, "");
    decompressor[1] = std::string("Gene1");
    decompressor[42] = std::string("5");
    decompressor[2] = std::string("Three");

    // ---- Make the test directories ----

    // make the test dir root
    std::filesystem::create_directories(temp_dir);

    // make the simulated output dir
    std::filesystem::create_directories(test_output_dir);
    std::filesystem::create_directories(test_subnets_dir);
    std::filesystem::create_directories(test_subnets_log_dir);

    // ---- Make the test ARACNe inputs ----

    std::ofstream ofs(test_exp_mat_file);
    ofs << std::string("gene\tsamp1\t2\t3.3\tsamp4\n") +
               "Gene1\t0.1\t0.2\t0.3\t0.4\n" + "5\t0.2\t0.3\t0.1\t0.15\n" +
               "Three\t0.\t0.1\t0.\t0.1";
    ofs.close();

    ofs = std::ofstream(test_reg_list_file);
    ofs << std::string("5\n") + "gene\n" + "Gene1\n" +
               "Gene2"; // "gene" and "Gene2" should not be added in memory"
    ofs.close();

    // ---- Make the test consolidate inputs ----

    // Create log files

    uint32_t tot_possible_edges = 4u;

    SubnetLogger sl1(test_subnets_log_dir / test_subnet_log_file_1);
    sl1.initSubnetLog("testid", 1, 2u, 3u, 500u, 300u, "FDR", 0.05f,
                      false); // did not MaxEnt prune
    sl1.write("\nRaw subnetwork computation time: 1s\n");
    sl1.write("Size of subnetwork: " + std::to_string(tot_possible_edges) +
              " edges.\n");
    sl1.write("\nThreshold pruning time (FDR): 1s\n");
    sl1.write("Edges removed: " +
              std::to_string(tot_possible_edges -
                             num_edges_after_threshold_pruning_1) +
              " edges.\n");
    sl1.write("Size of subnetwork: " +
              std::to_string(num_edges_after_threshold_pruning_1) +
              " edges.\n");
    sl1.write("\nPrinting subnetwork...");
    sl1.write("1s \n");

    SubnetLogger sl2(test_subnets_log_dir / test_subnet_log_file_2);
    sl2.initSubnetLog("testid", 1, 2u, 3u, 500u, 300u, "FDR", 0.05f,
                      true); // did MaxEnt prune
    sl2.write("\nRaw subnetwork computation time: 1s\n");
    sl2.write("Size of subnetwork: " + std::to_string(tot_possible_edges) +
              " edges.\n");
    sl2.write("\nThreshold pruning time (FDR): 1s\n");
    sl2.write("Edges removed: " +
              std::to_string(tot_possible_edges -
                             num_edges_after_threshold_pruning_2) +
              " edges.\n");
    sl2.write("Size of subnetwork: " +
              std::to_string(num_edges_after_threshold_pruning_2) +
              " edges.\n");
    sl2.write("\nMaxEnt pruning time: 1s\n");
    sl2.write("Edges removed: " +
              std::to_string(num_edges_after_threshold_pruning_2 -
                             num_edges_after_MaxEnt_pruning_2) +
              " edges.\n");
    sl2.write("Size of subnetwork: " +
              std::to_string(num_edges_after_MaxEnt_pruning_2) + " edges.\n");
    sl2.write("\nPrinting subnetwork...");
    sl2.write("1s \n");

    // Create network files

    ofs = std::ofstream(test_subnets_dir / test_subnet_file_1);
    ofs << std::string("regulator.values\ttarget.values\tmi.values\n") +
               "Gene1\t5\t0.0\n" + "Gene1\tThree\t0.5\n" + "5\tGene1\t0.0\n" +
               "5\tThree\t10\n";
    ofs.close();

    ofs = std::ofstream(test_subnets_dir / test_subnet_file_2);
    ofs << std::string("regulator.values\ttarget.values\tmi.values\n") +
               "5\tThree\t1\n";
    ofs.close();
  }

  void TearDown() override { std::filesystem::remove_all(temp_dir); }
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

  auto [regs, w_list] = io.readRegList(defined_genes);

  // There were 2 defined in the compressor, 2 not defined
  ASSERT_EQ(regs.size(), 2);
  ASSERT_EQ(w_list.size(), 2);
  ASSERT_TRUE(regs.find(42) != regs.end()); // "5" was in the regulators
  ASSERT_TRUE(regs.find(0) == regs.end());  // 0 should have no corresponding
}

TEST_F(FilesystemIOTest, WriteNetworkRegTarMI) {
  gene_to_gene_to_float subnet = {{0, {{1, 0.5f}, {2, 0.3f}}}, {3, {{0, 0.f}}}};
  decompression_map decompressor = {"gene1", "Gene 2", "3", "gene"};

  io.writeNetworkRegTarMI(1, subnet, decompressor);

  std::filesystem::path expected_file =
      test_subnets_dir / "subnetwork-1_testid.tsv";
  ASSERT_TRUE(std::filesystem::exists(expected_file));

  std::ifstream result_file(expected_file);
  ASSERT_TRUE(result_file.is_open());

  std::string line;
  std::getline(result_file, line);
  ASSERT_EQ(line, "regulator.values\ttarget.values\tmi.values");

  // Prepare expected results set. Because of range-based iteration,
  // the output line order may not match the data order.
  std::set<std::string> expected_lines = {"gene1\tGene 2\t0.5", "gene1\t3\t0.3",
                                          "gene\tgene1\t0"};

  while (std::getline(result_file, line)) {
    auto it = expected_lines.find(line);
    ASSERT_TRUE(it != expected_lines.end());
    expected_lines.erase(it);
  }

  // Verify that all expected lines were found
  ASSERT_TRUE(expected_lines.empty());
}

TEST_F(FilesystemIOTest, WriteARACNe3DF) {
  std::vector<ARACNe3_df> aracne3_final_df = {{0, 1, 0.5f, 0.8f, 3, -2.0f},
                                              {1, 2, 0.3f, 0.6f, 5, -3.5f},
                                              {2, 0, 0.2f, 0.9f, 2, -1.2f}};
  decompression_map decompressor = {"Gene 1", "5", "gene"};

  io.writeARACNe3DF(aracne3_final_df, decompressor);
  std::filesystem::path output_file = test_final_output_file;

  ASSERT_TRUE(std::filesystem::exists(output_file));

  std::ifstream result_file(output_file);
  ASSERT_TRUE(result_file.is_open());

  std::string line;
  std::getline(result_file, line);
  ASSERT_EQ(line, "regulator.values\ttarget.values\tmi.values\tscc."
                  "values\tcount.values\tlog.p.values");

  // Prepare expected results set. Because of range-based iteration,
  // the output line order may not match the data order.
  std::set<std::string> expected_lines = {"Gene 1\t5\t0.5\t0.8\t3\t-2",
                                          "5\tgene\t0.3\t0.6\t5\t-3.5",
                                          "gene\tGene 1\t0.2\t0.9\t2\t-1.2"};

  while (std::getline(result_file, line)) {
    auto it = expected_lines.find(line);
    ASSERT_TRUE(it != expected_lines.end());
    expected_lines.erase(it);
  }

  // Verify that all expected lines were found
  ASSERT_TRUE(expected_lines.empty());
}

TEST_F(FilesystemIOTest, FindSubnetFilesAndSubnetLogFiles) {
  pair_string_vecs results = io.findSubnetFilesAndSubnetLogFiles(
      test_subnets_dir, test_subnets_log_dir);

  using strvec = std::vector<std::string>;
  using strset = std::unordered_set<std::string>;

  strvec sn_paths = results.first;
  strvec sn_log_paths = results.second;

  ASSERT_EQ(sn_paths.size(), 2);
  ASSERT_EQ(sn_log_paths.size(), 2);

  strset sn_paths_set(sn_paths.begin(), sn_paths.end());
  strset sn_log_paths_set(sn_log_paths.begin(), sn_log_paths.end());

  ASSERT_TRUE(sn_paths_set.find(test_subnet_file_1) != sn_paths_set.end());
  ASSERT_TRUE(sn_paths_set.find(test_subnet_file_2) != sn_paths_set.end());
  ASSERT_TRUE(sn_log_paths_set.find(test_subnet_log_file_1) !=
              sn_log_paths_set.end());
  ASSERT_TRUE(sn_log_paths_set.find(test_subnet_log_file_2) !=
              sn_log_paths_set.end());
}

TEST_F(FilesystemIOTest, LoadARACNe3SubnetAndUpdateFPRFromLog) {
  std::filesystem::path snfp1 = test_subnets_dir / test_subnet_file_1;
  std::filesystem::path snlfp1 = test_subnets_log_dir / test_subnet_log_file_1;
  std::filesystem::path snfp2 = test_subnets_dir / test_subnet_file_2;
  std::filesystem::path snlfp2 = test_subnets_log_dir / test_subnet_log_file_2;

  auto [subnet_1, fpr_1] = io.loadARACNe3SubnetAndUpdateFPRFromLog(
      snfp1, snlfp1, defined_genes, decompressor, regulators);

  auto [subnet_2, fpr_2] = io.loadARACNe3SubnetAndUpdateFPRFromLog(
      snfp2, snlfp2, defined_genes, decompressor, regulators);

  gene_to_gene_to_float exp_subnet_1 = {
      {defined_genes["Gene1"],
       {{defined_genes["5"], 0.}, {defined_genes["Three"], 0.5}}},
      {defined_genes["5"],
       {{defined_genes["Gene1"], 0.}, {defined_genes["Three"], 10.}}}};

  gene_to_gene_to_float exp_subnet_2 = {
      {defined_genes["5"], {{defined_genes["Three"], 1.}}}};

  ASSERT_EQ(subnet_1, exp_subnet_1);
  ASSERT_EQ(fpr_1,
            estimateFPRNoMaxEnt(.05, "FDR", num_edges_after_threshold_pruning_1,
                                tot_possible_edges));

  ASSERT_EQ(subnet_2, exp_subnet_2);
  ASSERT_EQ(fpr_2, estimateFPRWithMaxEnt(
                       .05, "FDR", num_edges_after_threshold_pruning_2,
                       num_edges_after_MaxEnt_pruning_2, tot_possible_edges));
}

// TODO: Add tests for exception conditions among files.
