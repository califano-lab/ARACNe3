#include "subnet_operations.hpp"
#include "algorithms.hpp"
#include "apmi_nullmodel.hpp"
#include "aracne3io.hpp"
#include "subnet_logger.hpp"
#include "stopwatch.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>

const vv_float sampleExpMatAndReCopulaTransform(const vv_float &exp_mat,
                                                const uint16_t n_subsamp,
                                                std::mt19937 &rnd) {
  std::vector<uint16_t> idxs(exp_mat[0].size());
  std::iota(idxs.begin(), idxs.end(), 0U);

  std::vector<uint16_t> fold(n_subsamp);
  std::sample(idxs.begin(), idxs.end(), fold.begin(), n_subsamp, rnd);

  vv_float subsample_exp_mat(exp_mat.size(),
                             std::vector<float>(n_subsamp, 0.f));
  for (gene_id gene = 0u; gene < exp_mat.size(); ++gene) {
    for (uint16_t i = 0u; i < n_subsamp; ++i)
      subsample_exp_mat[gene][i] = exp_mat[gene][fold[i]];

    subsample_exp_mat[gene] = copulaTransform(subsample_exp_mat[gene], rnd);
  }

  return subsample_exp_mat;
}

/*
 Prunes a network by control of alpha using the Benjamini-Hochberg Procedure if
 method = FDR, or FWER if method = FWER.
 */
std::tuple<gene_to_gene_to_float, uint32_t, gene_to_gene_to_float>
pruneAlpha(const vv_float &network, const std::vector<gene_id> &regs_c,
           const std::vector<gene_id> &genes_c, uint32_t size_of_network,
           const std::string &method, const float alpha,
           const APMINullModel &nullmodel, const geneset &regulators) {

  // Flatten network
  std::vector<std::tuple<gene_id, gene_id, float>> reg_tar_mi;
  reg_tar_mi.reserve(size_of_network);

  for (uint32_t reg_idx = 0u; reg_idx < regs_c.size(); ++reg_idx) {
    const gene_id reg = regs_c.at(reg_idx);
    for (uint32_t tar_idx = 0u; tar_idx < genes_c.size(); ++tar_idx) {
      const gene_id tar = genes_c.at(tar_idx);
      if (reg != tar)
        reg_tar_mi.emplace_back(reg, tar, network.at(reg_idx).at(tar_idx));
    }
  }

  // sort descending
  std::sort(reg_tar_mi.begin(), reg_tar_mi.end(),
            [](const auto &rtm1, const auto &rtm2) {
              return std::get<2>(rtm1) > std::get<2>(rtm2);
            });

  uint32_t lowest_idx_that_doesnt_pass_thresh = 0u;
  if (method == "FDR") {
    // Benjamini-Hochberg
    for (auto it = reg_tar_mi.cbegin(); it != reg_tar_mi.cend(); ++it) {
      const std::size_t k = it - reg_tar_mi.cbegin();
      const float p_k = nullmodel.getMIPVal(std::get<2>(*it));
      if (p_k <= (k + 1.f) / size_of_network * alpha)
        lowest_idx_that_doesnt_pass_thresh = k + 1u;
    }
  } else if (method == "FWER") {
    for (auto it = reg_tar_mi.begin(); it != reg_tar_mi.end(); ++it) {
      const std::size_t k = it - reg_tar_mi.cbegin();
      const float p_k = nullmodel.getMIPVal(std::get<2>(*it));
      if (p_k <= alpha / size_of_network)
        lowest_idx_that_doesnt_pass_thresh = k + 1u;
    }
  } else if (method == "FPR") {
    for (auto it = reg_tar_mi.begin(); it != reg_tar_mi.end(); ++it) {
      const std::size_t k = it - reg_tar_mi.cbegin();
      const float p_k = nullmodel.getMIPVal(std::get<2>(*it));
      if (p_k <= alpha)
        lowest_idx_that_doesnt_pass_thresh = k + 1u;
    }
  }

  reg_tar_mi.erase(reg_tar_mi.cbegin() + lowest_idx_that_doesnt_pass_thresh,
                   reg_tar_mi.cend());

  // TODO: Remove all .at()
  // rebuild network
  gene_to_gene_to_float pruned_net;
  gene_to_gene_to_float pruned_net_reg_reg_only;

  pruned_net.reserve(regs_c.size());
  pruned_net_reg_reg_only.reserve(regs_c.size());

  for (const auto &[reg, tar, mi] : reg_tar_mi) {
    pruned_net[reg].insert({tar, mi});
    if (regulators.find(tar) != regulators.end())
      pruned_net_reg_reg_only[reg].insert({tar, mi});
  }

  return {pruned_net, reg_tar_mi.size(), pruned_net_reg_reg_only};
}

/*
 Prune the network according to the MaxEnt weakest-edge reduction.
 */
std::pair<gene_to_gene_to_float, uint32_t>
pruneMaxEnt(gene_to_gene_to_float network, uint32_t size_of_network,
            const geneset &regulators,
            gene_to_gene_to_float network_reg_reg_only,
            const uint16_t nthreads) {

  gene_to_geneset edges_to_remove;
  edges_to_remove.reserve(regulators.size());

  // make unique data structure that stores redundant edges
  std::set<std::pair<gene_id, gene_id>>
      to_remove; // can't use unordered because hash fn not defined for pair;
                 // set is binary tree
  for (const auto &[reg1, reg2_mi] : network_reg_reg_only)
    for (const auto [reg2, mi_regs] : reg2_mi)
      if (to_remove.find({reg1, reg2}) == to_remove.end() &&
          to_remove.find({reg2, reg1}) == to_remove.end())
        to_remove.insert({reg2, reg1});

  // make into triangular matrix (remove all reg from all reg2 targets regulon)
  for (const auto &[reg1, reg2] : to_remove)
    network_reg_reg_only.at(reg1).erase(reg2);

#pragma omp parallel num_threads(nthreads)
  {
    // Local version for each thread
    gene_to_geneset local_edges_to_remove;

    // schedule in skips as opposed to chunks, faster now
#pragma omp for schedule(static, 1)
    for (int i = 0u; i < static_cast<int>(network_reg_reg_only.size()); ++i) {
      auto it = network_reg_reg_only.cbegin();
      std::advance(it, i);
      const uint16_t reg1 = it->first;
      const gene_to_float &reg2_mi = it->second;

      const gene_to_float &reg1_regulon = network.at(reg1);
      geneset &remove_from_reg1 = local_edges_to_remove[reg1];

      for (const auto [reg2, mi_regs] : reg2_mi) {
        // check if reg2 has regulon
        if (network.find(reg2) != network.end()) {
          const gene_to_float &reg2_regulon = network.at(reg2);
          geneset &remove_from_reg2 = local_edges_to_remove[reg2];

          for (const auto [tar, mi_reg1_tar] : reg1_regulon) {

            if (reg2_regulon.find(tar) != reg2_regulon.end()) {
              const float mi_reg2_tar = reg2_regulon.at(tar);
              if (mi_reg1_tar < mi_regs && mi_reg1_tar < mi_reg2_tar)
                remove_from_reg1.insert(tar);
              else if (mi_reg2_tar < mi_regs && mi_reg2_tar < mi_reg1_tar)
                remove_from_reg2.insert(tar);
              else {
                remove_from_reg1.insert(reg2);
                remove_from_reg2.insert(reg1);
              }
            }
          }
        }
      }
    }

// Merging local results into global result
#pragma omp critical
    {
      for (const auto &[reg, remove_set] : local_edges_to_remove) {
        edges_to_remove[reg].insert(remove_set.begin(), remove_set.end());
      }
    }
  }

  for (const auto &[reg, remove_from_reg] : edges_to_remove) {
    for (const gene_id tar : remove_from_reg)
      network[reg].erase(tar);
    size_of_network -= remove_from_reg.size();
  }

  return std::make_pair(network, size_of_network);
}

/*
 Generates an ARACNe3 subnet (called from main).
*/
std::tuple<gene_to_gene_to_float, float, uint32_t> createARACNe3Subnet(
    const vv_float &subsample_exp_mat, const geneset &regulators,
    const geneset &genes, const uint16_t tot_num_samps,
    const uint16_t tot_num_subsample, const uint16_t subnet_number,
    const bool prune_alpha, const APMINullModel &nullmodel,
    const std::string &method, const float alpha, const bool prune_MaxEnt,
    const std::string &subnets_log_dir, const uint16_t nthreads,
    const std::string &runid, const decompression_map &decompressor,
    const bool save_subnet, const ARACNe3IOHandler &io) {

  std::unique_ptr<SubnetLogger> subnet_logger;

  // ---- Quick macros ----

  auto qlog = [&](const std::string &cur_msg) {
    if (subnet_logger)
      subnet_logger->write(cur_msg);
    return;
  };

  if (save_subnet) {
    const std::string log_file_name = subnets_log_dir + "log-subnetwork-" +
                                      std::to_string(subnet_number) + "_" +
                                      runid + ".txt";
    subnet_logger = std::make_unique<SubnetLogger>(log_file_name);
    subnet_logger->initSubnetLog(runid, subnet_number, regulators, genes,
                                 tot_num_samps, tot_num_subsample, method,
                                 alpha, prune_MaxEnt);
  }

  // ---- Generate raw subnetwork ----

  qlog("\nRaw subnetwork computation time: ");
  Stopwatch watch1{};

  const uint32_t tot_possible_edges = regulators.size() * (genes.size() - 1);

  // vectorize sets and network for parallelism
  const std::vector<gene_id> regs_c(regulators.begin(), regulators.end()),
      genes_c(genes.begin(), genes.end());

  vv_float subnetwork_vec(regulators.size(),
                          std::vector<float>(genes.size(), 0.f));

#pragma omp parallel for num_threads(nthreads)
  for (int reg_idx = 0; reg_idx < static_cast<int>(regs_c.size()); ++reg_idx) {
    const gene_id reg = regs_c.at(reg_idx);
    for (uint32_t tar_idx = 0u; tar_idx < genes_c.size(); ++tar_idx) {
      const gene_id tar = genes_c.at(tar_idx);
      if (reg != tar)
        subnetwork_vec.at(reg_idx).at(tar_idx) =
            calcAPMI(subsample_exp_mat.at(reg), subsample_exp_mat.at(tar));
    }
  }

  qlog(watch1.getSeconds() + "\n");
  qlog("Size of subnetwork: " + std::to_string(tot_possible_edges) +
       " edges.\n");

  // ---- Prune by threshold ----

  qlog("\nThreshold pruning time (" + method + "): ");
  watch1.reset();

  // unpack tuple into objects
  uint32_t num_edges_after_threshold_pruning;
  gene_to_gene_to_float subnetwork, subnetwork_reg_reg_only;

  if (prune_alpha)
    std::tie(subnetwork, num_edges_after_threshold_pruning,
             subnetwork_reg_reg_only) =
        pruneAlpha(subnetwork_vec, regs_c, genes_c, tot_possible_edges, method,
                   alpha, nullmodel, regulators);
  else
    num_edges_after_threshold_pruning = tot_possible_edges;

  subnetwork_vec.clear();
  vv_float(subnetwork_vec).swap(subnetwork_vec);

  qlog(watch1.getSeconds() + "\n");
  qlog("Edges removed: " +
       std::to_string(tot_possible_edges - num_edges_after_threshold_pruning) +
       " edges.\n");
  qlog("Size of subnetwork: " +
       std::to_string(num_edges_after_threshold_pruning) + " edges.\n");

  // ---- Prune by maximizing entropy ----

  uint32_t size_of_subnet;
  float FPR_estimate_subnet;

  if (prune_MaxEnt) {
    qlog("\nMaxEnt pruning time: ");
    watch1.reset();

    uint32_t num_edges_after_MaxEnt_pruning;

    std::tie(subnetwork, num_edges_after_MaxEnt_pruning) =
        pruneMaxEnt(subnetwork, num_edges_after_threshold_pruning, regulators,
                    subnetwork_reg_reg_only, nthreads);

    qlog(watch1.getSeconds() + "\n");
    qlog("Edges removed: " +
         std::to_string(num_edges_after_threshold_pruning -
                        num_edges_after_MaxEnt_pruning) +
         " edges.\n");
    qlog("Size of subnetwork: " +
         std::to_string(num_edges_after_MaxEnt_pruning) + " edges.\n");

    FPR_estimate_subnet = estimateFPRWithMaxEnt(
        alpha, method, num_edges_after_threshold_pruning,
        num_edges_after_MaxEnt_pruning, tot_possible_edges);

    size_of_subnet = num_edges_after_MaxEnt_pruning;
  } else {
    estimateFPRNoMaxEnt(alpha, method, num_edges_after_threshold_pruning,
                        tot_possible_edges);

    size_of_subnet = num_edges_after_threshold_pruning;
  }

  // ---- Return and/or print subnetwork  ----

  if (save_subnet) {
    qlog("\nPrinting subnetwork...");
    watch1.reset();

    io.writeNetworkRegTarMI(subnet_number, subnetwork, decompressor);

    qlog(watch1.getSeconds() + "\n");
  }

  return {subnetwork, FPR_estimate_subnet, size_of_subnet};
}

const std::vector<ARACNe3_df>
consolidateSubnetsVec(const std::vector<gene_to_gene_to_float> &subnets,
                      const float FPR_estimate, const vv_float &exp_mat,
                      const geneset &regulators, const geneset &genes) {
  std::vector<ARACNe3_df> final_df;

  for (const gene_id reg : regulators) {
    for (const gene_id tar : genes) {
      uint16_t num_occurrences = 0u;
      for (uint16_t sn = 0u; sn < subnets.size(); ++sn) {
        if (subnets.at(sn).find(reg) != subnets[sn].end())
          if (subnets.at(sn).at(reg).find(tar) != subnets.at(sn).at(reg).end())
            ++num_occurrences;
      }
      if (num_occurrences > 0) {
        const float final_mi = calcAPMI(exp_mat.at(reg), exp_mat.at(tar));
        const float final_scc = pearsonsR(exp_mat.at(reg), exp_mat.at(tar));
        const double final_log_p =
            lRightTailBinomialP(subnets.size(), num_occurrences, FPR_estimate);
        final_df.emplace_back(reg, tar, final_mi, final_scc, num_occurrences,
                              final_log_p);
      }
    }
  }

  return final_df;
}
