#pragma once

#include <string>
#include "ARACNe3.hpp"


typedef struct consolidated_df_row {
	const gene_id regulator;
	const gene_id target;
	const float final_mi;
	const float final_scc;
	const uint16_t num_subnets_incident;
	const double final_p;
	consolidated_df_row(const gene_id& r, const gene_id& t, const float& mi, const float& scc, const uint16_t& n, const double& p) : regulator(r), target(t), num_subnets_incident(n), final_mi(mi), final_scc(scc), final_p(p) {};
} consolidated_df_row;

std::string makeUnixDirectoryNameUniversal(std::string &dir_name);
std::string makeUnixDirectoryNameUniversal(std::string &&dir_name);
void makeDir(std::string &dir_name);
std::vector<uint16_t> rank_indexes(const std::vector<float>& vec, std::mt19937 &rand);

std::unordered_set<gene_id> readRegList(const std::string &filename);
std::pair<gene_to_floats, gene_to_shorts> readExpMatrixAndCopulaTransform(const std::string &filename, std::mt19937 &rand);

gene_to_floats sampleExpMatAndReCopulaTransform(gene_to_floats &exp_mat, std::mt19937 &rand);

void writeNetworkRegTarMI(gene_to_gene_to_float &network, const std::string &output_dir, const std::string &output_suffix);

void writeConsolidatedNetwork(const std::vector<consolidated_df_row>& final_df, std::string filename);

class TooManySubnetsRequested : public std::exception {	
public:
	const char * what() {
		return "TRIED TO CONSOLIDATE TOO MANY NETWORKS.  LOG/SUBNET FILE DOES NOT EXIST UP TO SPECIFIED NUMBER OF REQUESTED SUBNETS TO CONSOLIDATE.";
	}
};
gene_to_edge_tars readSubNetAndUpdateFPRFromLog(const std::string &output_dir, const uint16_t subnet_num);

