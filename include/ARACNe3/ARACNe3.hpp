#pragma once

#include <vector>
#include <string>
#include <sstream>
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
#include <iomanip>
#include <omp.h>

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
typedef uint16_t gene_id;

/*
 * The "ball on stick" data structure for edges in network.  Will be used with
 * dictionary so that a regulator can have a group of these associated.
 */
typedef struct edge_tar {
	gene_id target;
	float mi;
	edge_tar(const gene_id &t, const float &mi) : target(t), mi(mi) {};
	
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

// Maps gene to regulon
typedef std::unordered_map<gene_id, std::unordered_set<gene_id>>
    gene_to_geneset;
typedef std::unordered_map<gene_id, std::vector<edge_tar>> gene_to_edge_tars;

// Used for gexp and ranked gexp matrix storage
typedef std::unordered_map<gene_id, std::vector<float>> gene_to_floats;
typedef std::unordered_map<gene_id, std::vector<uint16_t>> gene_to_shorts;

// used for edge strength
typedef std::unordered_map<gene_id, std::unordered_map<gene_id, float>>
    gene_to_gene_to_float;
typedef std::unordered_map<gene_id, std::unordered_map<gene_id, float>>
    gene_to_gene_to_float;

//--------------------- APMI.cpp	 		-----------------------

float APMI(const std::vector<float>& vec_x, const std::vector<float>& vec_y, const float& q_thresh = 7.815, 
		const uint16_t& size_thresh = 4);

std::vector<edge_tar> gene_to_floatsAPMI(gene_to_floats &matrix, const gene_id& identifier, const float& q_thresh = 7.815, const uint16_t& size_thresh = 4);

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
std::pair<gene_to_edge_tars, gene_to_gene_to_float> pruneAlpha(gene_to_edge_tars &network, uint32_t& size_of_network);

//--------------------- MaxEntPruning.cpp	 	-----------------------
gene_to_edge_tars pruneMaxEnt(gene_to_edge_tars &network, gene_to_gene_to_float tftfNetwork, uint32_t &size_of_network);

//--------------------- RegWebFns.cpp	 		-----------------------
gene_to_edge_tars sort_edge_tars(gene_to_edge_tars &regweb);	

/*
 Used for MaxEnt pruning to copy ARACNe-AP; we should probably either overhaul our own data structure or find a way to use our own, as any conversion is essentially a memory and runtime cost.
 */
std::unordered_map<gene_id, std::unordered_map<gene_id, float>> regweb_to_mapmap(gene_to_edge_tars &network);

//--------------------- ARACNe3.cpp	 		-----------------------
