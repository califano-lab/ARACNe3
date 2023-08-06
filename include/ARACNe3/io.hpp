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
std::vector<uint16_t> rank_indexes(const std::vector<float>& vec);

/*
 Does not return a list of regulators as a string vector, as we are using compression
 */
void readRegList(std::string &filename);

/*
 Returns a map of gene identifier -> gene expression.
 */
void readExpMatrix(std::string &filename);

genemap sampleFromGlobalGenemap();

void writeNetworkRegTarMI(const reg_web &network, const std::string &output_dir, const std::string &output_suffix);

void writeConsolidatedNetwork(const std::vector<consolidated_df_row>& final_df, std::string filename);

class TooManySubnetsRequested : public std::exception {	
public:
	const char * what() {
		return "TRIED TO CONSOLIDATE TOO MANY NETWORKS.  LOG/SUBNET FILE DOES NOT EXIST UP TO SPECIFIED NUMBER OF REQUESTED SUBNETS TO CONSOLIDATE.";
	}
};
reg_web readSubNetAndUpdateFPRFromLog(const std::string &output_dir, const uint16_t subnet_num);

