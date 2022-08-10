#include <vector>
#include <string>
#include <tuple>
#include <unordered_set>
#include <set>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <random>
#include <thread>
#include <numeric>
#include <math.h>
#include <filesystem>
#include <chrono>
#include <regex>

#ifndef _ARACNE3_HPP_
#define _ARACNE3_HPP_

/*
 I assume we will need conditional compilation for multithreading and stuff later on.
 */
#if defined __linux__ || defined __APPLE__
const std::string hiddenfpre = ".";
const char directory_slash = '/';
#elif defined _WIN32 // in windows you must set a hidden file via properties
const char directory_slash = '\\';
const std::string hiddenfpre = "";
#else
const std::string hiddenfpre = "";
#endif /* __linux__ || __APPLE__ */


//--------------------- Overhead Data Structures	 -----------------------
/*
 Switching to this type alias to reduce confusion in other parts of the program
 */
typedef uint16_t gene_id_t;

/*
 * The "ball on stick" data structure for edges in network.  Will be used with
 * dictionary so that a regulator can have a group of these associated.
 */
typedef struct edge_tar{
	gene_id_t target;
	float mi;
	edge_tar(const gene_id_t &t, const float &mi) : target(t), mi(mi) {};
	
	/*
	 Note, the edge_tar structs are compared based on 'target identifier' and NOT MI
	 */
	bool operator<(const edge_tar& et) const
	{
	    return target < et.target;
	}
	
	bool operator>(const edge_tar& et) const
	{
	    return target > et.target;
	}
	
	bool operator==(const edge_tar& et) const
	{
	    return target == et.target;
	}
} edge_tar;

typedef struct consolidated_df {
	const gene_id_t regulator;
	const gene_id_t target;
	const float final_mi;
	const float final_scc;
	const uint16_t num_subnets_incident;
	const double final_p;
	consolidated_df(const gene_id_t& r, const gene_id_t& t, const float& mi, const float& scc, const uint16_t& n, const double& p) : regulator(r), target(t), num_subnets_incident(n), final_mi(mi), final_scc(scc), final_p(p) {};
} consolidated_df;

/*
 * The regulator and target list is represented by this data type, an unordered
 * hash map between regulator names (string or number in uncompressed/compressed
 * implementation) as well as a vector of "ball-on-stick" structs, or all the
 * targets and their MIs.
 */
typedef std::unordered_map<gene_id_t, std::vector<edge_tar>> reg_web;
typedef std::unordered_map<gene_id_t, std::unordered_map<gene_id_t, float>> map_map;


//--------------------- IO.cpp 		-----------------------
std::string makeUnixDirectoryNameUniversal(std::string &dir_name);
std::string makeUnixDirectoryNameUniversal(std::string &&dir_name);
void makeDir(std::string &dir_name);
std::vector<uint16_t> rank_indexes(const std::vector<float>& vec);

/*
 * Maps gene identifiers to gene expression matrices
 */
typedef std::unordered_map<gene_id_t, std::vector<float>> genemap;

/*
 Same as genemap, except mapped to ranked vectors
 */
typedef std::unordered_map<gene_id_t, std::vector<uint16_t>> genemap_r;

/*
 Does not return a list of regulators as a string vector, as we are using compression
 */
void readRegList(std::string &filename);

/*
 Returns a map of gene identifier -> gene expression.
 */
void readExpMatrix(std::string &filename);

genemap sampleFromGlobalGenemap();

void writeNetworkRegTarMI(const reg_web &network, std::string &output_dir, const std::string &output_suffix);

void writeConsolidatedNetwork(const std::vector<consolidated_df>& final_df, std::string filename);


//--------------------- APMI.cpp	 		-----------------------

float APMI(const std::vector<float>& vec_x, const std::vector<float>& vec_y, const float& q_thresh = 7.815, 
		const uint16_t& size_thresh = 4);

std::vector<edge_tar> genemapAPMI(genemap &matrix, const gene_id_t& identifier, const float& q_thresh = 7.815, const uint16_t& size_thresh = 4);

const std::vector<float> permuteAPMI(const std::vector<float> &ref_vec,
		const std::vector<std::vector<float>> &target_vec, const float &q_thresh = 7.815,
		const uint16_t &size_thresh = 4);


//--------------------- NullModel.cpp	 		-----------------------
// sets a file static variable of an ordered null distribution
const std::vector<float> initNullMIs(const uint16_t& tot_num_subsample);

// the function below requires that initNullMIs has been called
const float getMIPVal(const float& mi, const float& p_precise = 0.001f);

// the function below requires that initNullMIs has been called
const std::vector<float> getMIPVals(const std::vector<float>& mis, const float& p_precise = 0.001f);


//--------------------- AlphaPruning.cpp	 	-----------------------
std::pair<reg_web, map_map> pruneAlpha(reg_web &network, uint32_t& size_of_network);

//--------------------- MaxEntPruning.cpp	 	-----------------------
reg_web pruneMaxEnt(reg_web &network, map_map &tftfNetwork, uint32_t &size_of_network);

//--------------------- RegWebFns.cpp	 		-----------------------
reg_web sort_edge_tars(reg_web &regweb);	

/*
 Used for MaxEnt pruning to copy ARACNe-AP; we should probably either overhaul our own data structure or find a way to use our own, as any conversion is essentially a memory and runtime cost.
 */
std::unordered_map<gene_id_t, std::unordered_map<gene_id_t, float>> regweb_to_mapmap(reg_web &network);
//--------------------- Consolidator.cpp	 	-----------------------
std::vector<consolidated_df> consolidate_subnets_vec(std::vector<reg_web> &subnets);


//--------------------- ARACNe3.cpp	 		-----------------------


#endif /* #ifndef _ARACNE3_HPP_ */
