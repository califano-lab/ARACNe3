#include "subnet_operations.hpp"
#include "ARACNe3.hpp"
#include "algorithms.hpp"
#include "apmi_nullmodel.hpp"
#include "io.hpp"
#include "stopwatch.hpp"
#include <boost/math/distributions/beta.hpp>
#include <fstream>
#include <iostream>
#include <omp.h>

/*
 Prunes a network by control of alpha using the Benjamini-Hochberg Procedure if
 method = FDR, or FWER if method = FWER.
 */
std::tuple<gene_to_gene_to_float, uint32_t, gene_to_gene_to_float>
pruneAlpha(const gene_to_gene_to_float &network, const geneset &regulators,
           uint32_t size_of_network, const std::string &method,
           const float alpha, const APMINullModel &nullmodel) {

  /* A vector that describes each regulator-mi-target interaction must be
   * initialized for sorting-based pruning */
  std::vector<std::tuple<gene_id, gene_id, float>> reg_tar_mi;
  reg_tar_mi.reserve(size_of_network);

  for (gene_id reg : regulators)
    for (const auto [tar, mi] : network.at(reg))
      reg_tar_mi.emplace_back(reg, tar, mi);

  // sort descending
  std::sort(reg_tar_mi.begin(), reg_tar_mi.end(),
            [](const std::tuple<gene_id, gene_id, float> &rtm1,
               const std::tuple<gene_id, gene_id, float> &rtm2) -> bool {
              return std::get<2>(rtm1) > std::get<2>(rtm2);
            });

  uint32_t argmax_k = 0;
  uint32_t m = size_of_network;
  if (method == "FDR") {
    // Benjamini-Hochberg
    for (auto it = reg_tar_mi.begin(); it != reg_tar_mi.end(); ++it) {
      const auto k = it - reg_tar_mi.begin();
      const float p_k = nullmodel.getMIPVal(std::get<2>(*it));
      if (p_k < k * alpha / m)
        argmax_k = static_cast<uint32_t>(k) + 1;
    }
  } else if (method == "FWER") {
    for (auto it = reg_tar_mi.begin(); it != reg_tar_mi.end(); ++it) {
      const auto k = it - reg_tar_mi.begin();
      const float p_k = nullmodel.getMIPVal(std::get<2>(*it));
      if (p_k < alpha / m)
        argmax_k = static_cast<uint32_t>(k) + 1;
    }
  } else if (method == "FPR") {
    // Only for benchmarking
    for (auto it = reg_tar_mi.begin(); it != reg_tar_mi.end(); ++it) {
      const auto k = it - reg_tar_mi.begin();
      const float p_k = nullmodel.getMIPVal(std::get<2>(*it));
      if (p_k < alpha)
        argmax_k = static_cast<uint32_t>(k) + 1;
    }
  }

  // create the new vector that is a pruned version of original
  std::vector<std::tuple<gene_id, gene_id, float>> pruned_vec(
      &reg_tar_mi[0], &reg_tar_mi[argmax_k]);

  // rebuild network
  size_of_network = pruned_vec.size();
  gene_to_gene_to_float pruned_net;
  gene_to_gene_to_float pruned_net_reg_reg_only;

  pruned_net.reserve(regulators.size());
  pruned_net_reg_reg_only.reserve(regulators.size());

  for (const auto [reg, tar, mi] : pruned_vec) {
    pruned_net[reg][tar] = mi;
    if (regulators.find(tar) != regulators.end())
      pruned_net_reg_reg_only[reg][tar] = mi;
  }

  return std::make_tuple(pruned_net, size_of_network, pruned_net_reg_reg_only);
}

/*
 Prune the network according to the MaxEnt weakest-edge reduction.
 */
std::pair<gene_to_gene_to_float, uint32_t>
pruneMaxEnt(gene_to_gene_to_float network, uint32_t size_of_network,
            const geneset &regulators,
            gene_to_gene_to_float network_reg_reg_only) {

  gene_to_geneset edges_to_remove;
  edges_to_remove.reserve(regulators.size());

  // triangular matrix (reg and reg+1)
  for (const auto [reg1, reg2_mi] : network_reg_reg_only) {
    const std::unordered_map<gene_id, float> &reg1_regulon = network.at(reg1);
    std::unordered_set<gene_id> &remove_from_reg1 = edges_to_remove[reg1];
    for (const auto [reg2, mi_regs] : reg2_mi) {
      network_reg_reg_only[reg2].erase(reg1);
      const std::unordered_map<gene_id, float> &reg2_regulon = network.at(reg2);
      std::unordered_set<gene_id> &remove_from_reg2 = edges_to_remove[reg2];
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
std::pair<gene_to_gene_to_float, float> createARACNe3Subnet(
    const gene_to_floats &subsample_exp_mat, const geneset &regulators,
    const geneset &genes, const uint16_t tot_num_samps,
    const uint16_t tot_num_subsample, const uint16_t cur_subnet_ct,
    const bool prune_alpha, const APMINullModel &nullmodel,
    const std::string &method, const float alpha, const bool prune_MaxEnt,
    const std::string &output_dir, const std::string &subnets_dir,
    const std::string &subnets_log_dir, const uint16_t nthreads,
    const std::string &runid) {

  float FPR_estimate_subnet;
  std::ofstream log_output(subnets_log_dir + "log_subnet" +
                           std::to_string(cur_subnet_ct + 1) + "_" + runid +
                           ".txt");
  std::time_t t = std::time(nullptr);

  log_output << "---------" << std::put_time(std::localtime(&t), "%c %Z")
             << "---------" << std::endl
             << std::endl;
  log_output << "Subnetwork #: " + std::to_string(cur_subnet_ct + 1) << std::endl;
  log_output << "Total # regulators (defined in gexp mat): " +
                    std::to_string(regulators.size())
             << std::endl;
  log_output << "Total # targets: " + std::to_string(genes.size()) << std::endl;
  log_output << "Total # samples: " + std::to_string(tot_num_samps)
             << std::endl;
  log_output << "Subsampled quantity: " + std::to_string(tot_num_subsample)
             << std::endl;
  log_output << "Total possible edges: " +
                    std::to_string(regulators.size() * (genes.size() - 1))
             << std::endl;
  log_output << "Method of first pruning step: " + method << std::endl;
  log_output << "Alpha: " + std::to_string(alpha) << std::endl;
  log_output << "MaxEnt Pruning: " +
                    std::string(prune_MaxEnt ? "true" : "false")
             << std::endl;
  log_output << "\n-----------Begin Network Generation-----------" << std::endl;

  // begin subnet computation

  //-------time module-------
  Watch watch1;
  log_output << "\nRaw subnetwork computation time: ";
  log_output.flush();
  watch1.reset();
  //-------------------------

  uint32_t size_of_subnetwork = 0U;

  // vectorize sets and network for parallelism
  const std::vector<gene_id> regs_vec(regulators.begin(), regulators.end()),
      genes_vec(genes.begin(), genes.end());

  std::vector<std::vector<float>> subnetwork_vec(
      regulators.size(), std::vector<float>(genes.size(), 0.f));

#pragma omp parallel for num_threads(nthreads)
  for (int reg_idx = 0; reg_idx < regulators.size(); ++reg_idx) {
    const gene_id reg = regs_vec[reg_idx];
    for (int tar_idx = 0; tar_idx < genes.size(); ++tar_idx) {
      const gene_id tar = genes_vec[tar_idx];
      if (reg != tar)
        subnetwork_vec[reg_idx][tar_idx] =
            calcAPMI(subsample_exp_mat.at(reg), subsample_exp_mat.at(tar));
    }
  }

  // transfer back to hash map structure
  gene_to_gene_to_float subnetwork;
  subnetwork.reserve(regulators.size());
  for (uint16_t reg_idx = 0U; reg_idx < regulators.size(); ++reg_idx) {
    const gene_id reg = regs_vec[reg_idx];
    for (uint16_t tar_idx = 0U; tar_idx < genes.size(); ++tar_idx) {
      const gene_id tar = genes_vec[tar_idx];
      if (reg != tar) {
        subnetwork[reg][tar] = subnetwork_vec[reg_idx][tar_idx];
        ++size_of_subnetwork;
      }
    }
  }

  //-------time module-------
  log_output << watch1.getSeconds() << std::endl;
  log_output << "Size of subnetwork: " << size_of_subnetwork << " edges."
             << std::endl;
  //-------------------------

  //-------time module-------
  log_output << "\nThreshold pruning time (" + method + "): ";
  log_output.flush();
  watch1.reset();
  //-------------------------

  uint32_t size_prev = size_of_subnetwork;

  // unpack tuple into objects
  gene_to_gene_to_float subnetwork_reg_reg_only;

  std::tie(subnetwork, size_of_subnetwork, subnetwork_reg_reg_only) =
      pruneAlpha(subnetwork, regulators, size_of_subnetwork, method, alpha,
                 nullmodel);

  //-------time module-------
  log_output << watch1.getSeconds() << std::endl;
  log_output << "Edges removed: " << size_prev - size_of_subnetwork << " edges."
             << std::endl;
  log_output << "Size of subnetwork: " << size_of_subnetwork << " edges."
             << std::endl;
  //-------------------------

  // save for binomial distribution parameter (theta)
  uint32_t num_edges_after_threshold_pruning = size_of_subnetwork;

  if (prune_MaxEnt) {
    //-------time module-------
    log_output << "\nMaxEnt pruning time: ";
    log_output.flush();
    watch1.reset();
    //-------------------------

    size_prev = size_of_subnetwork;
    std::tie(subnetwork, size_of_subnetwork) = pruneMaxEnt(
        subnetwork, size_of_subnetwork, regulators, subnetwork_reg_reg_only);

    //-------time module-------
    log_output << watch1.getSeconds() << std::endl;
    log_output << "Edges removed: " << size_prev - size_of_subnetwork
               << " edges." << std::endl;
    log_output << "Size of subnetwork: " << size_of_subnetwork << " edges."
               << std::endl;
    //-------------------------

    uint32_t num_edges_after_MaxEnt_pruning = size_of_subnetwork;
    if (method == "FDR")
      FPR_estimate_subnet = (alpha * num_edges_after_MaxEnt_pruning) /
                            (regulators.size() * (genes.size() - 1) -
                             (1 - alpha) * num_edges_after_threshold_pruning);
    else if (method == "FWER")
      FPR_estimate_subnet = (alpha / (regulators.size() * (genes.size() - 1))) *
                            (num_edges_after_MaxEnt_pruning) /
                            (num_edges_after_threshold_pruning);
    else if (method == "FPR")
      FPR_estimate_subnet = alpha * num_edges_after_MaxEnt_pruning /
                            num_edges_after_threshold_pruning;
  } else {
    if (method == "FDR")
      FPR_estimate_subnet = (alpha * num_edges_after_threshold_pruning) /
                            (regulators.size() * (genes.size() - 1) -
                             (1 - alpha) * num_edges_after_threshold_pruning);
    else if (method == "FWER")
      FPR_estimate_subnet = alpha / (regulators.size() * (genes.size() - 1));
    else if (method == "FPR")
      FPR_estimate_subnet = alpha;
  }

  //-------time module-------
  log_output << "\nPrinting subnetwork in directory \"" + subnets_dir + "\"...";
  log_output.flush();
  watch1.reset();
  //-------------------------

  // writes the individual subnet output
  writeNetworkRegTarMI(subnetwork, subnets_dir + "subnet" +
                                       std::to_string(cur_subnet_ct + 1) + "_" +
                                       runid + ".tsv");

  //-------time module-------
  log_output << watch1.getSeconds() << std::endl;
  //-------------------------

  std::cout << "...subnetwork " + std::to_string(cur_subnet_ct + 1) +
                   " completed = " + std::to_string(size_of_subnetwork) +
                   " edges returned..."
            << std::endl;

  return std::make_pair(subnetwork, FPR_estimate_subnet);
}

const std::vector<consolidated_df_row>
consolidateSubnetsVec(const std::vector<gene_to_gene_to_float> &subnets,
                      const float FPR_estimate, const gene_to_floats &exp_mat,
                      const geneset &regulators, const geneset &genes,
                      const gene_to_shorts &ranks_mat) {
  std::vector<consolidated_df_row> final_df;
  const uint32_t tot_poss_edgs = regulators.size() * (genes.size() - 1);

  for (const gene_id &reg : regulators) {
    for (const gene_id &tar : genes) {
      uint16_t num_occurrences = 0;
      for (uint16_t sn = 0U; sn < subnets.size(); ++sn) {
        if (subnets[sn].find(reg) != subnets[sn].end())
          if (subnets[sn].at(reg).find(tar) != subnets[sn].at(reg).end())
            ++num_occurrences;
      }
      if (num_occurrences > 0) {
        const float final_mi = calcAPMI(exp_mat.at(reg), exp_mat.at(tar));
        const float final_scc = calcSCC(ranks_mat.at(reg), ranks_mat.at(tar));
        const double final_log_p = std::log(
            rightTailBinomialP(subnets.size(), num_occurrences, FPR_estimate));
        final_df.emplace_back(reg, tar, final_mi, final_scc, num_occurrences,
                              final_log_p);
      }
    }
  }

  return final_df;
}
