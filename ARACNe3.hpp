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

#ifndef _ARACNE3_HPP_
#define _ARACNE3_HPP_

/*
 I assume we will need conditional compilation for multithreading and stuff later on.
 */
#if defined __linux__ || defined __APPLE__
const std::string hiddenfpre = ".";
#elif defined _WIN32 // in windows you must set a hidden file via properties
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
	gene_id_t regulator;
	gene_id_t target;
	uint16_t num_subnets_incident;
	float final_mi;
	float final_scc;
	float final_p;
	consolidated_df(const gene_id_t& r, const gene_id_t& t, const uint16_t& n, const float& mi, const float& scc, const float& p) : regulator(r), target(t), num_subnets_incident(n), final_mi(mi), final_scc(scc), final_p(p) {};
} consolidated_df;

/*
 * The regulator and target list is represented by this data type, an unordered
 * hash map between regulator names (string or number in uncompressed/compressed
 * implementation) as well as a vector of "ball-on-stick" structs, or all the
 * targets and their MIs.
 */
typedef std::unordered_map<gene_id_t, std::vector<edge_tar>> reg_web;
typedef std::unordered_map<gene_id_t, std::unordered_map<gene_id_t, float>> map_map;


//--------------------- MatrixReglistIO.cpp 		-----------------------
void makeDir(const std::string &dir_name);

/*
 * Maps gene identifiers to gene expression matrices
 */
typedef std::unordered_map<gene_id_t, std::vector<float>> genemap;

/*
 Does not return a list of regulators as a string vector, as we are using compression
 */
void readRegList(std::string filename = "regulators.txt");

/*
 Returns a map of gene identifier -> gene expression.
 */
std::vector<genemap> readExpMatrix(std::string filename = "exp_mat.txt");

void writeNetworkRegTarMI(const reg_web &network, const std::string &output_dir = "output", const std::string &output_suffix = "0");


//--------------------- APMI.cpp	 		-----------------------
/*
 * Square struct for APMI estimator.
 * 'x_bound1' is the x-coordinate of the bottom left of square
 * 'y_bound1' is the y-coordinate of the bottom left of square
 * 'width' is the width of the square
 * 'pts' is an array of indices
 * 'num_pts' is the size of 'pts'
 */
// Can we somehow make uint16_t *pts into const uint16_t *&pts?
typedef struct {const float &x_bound1, &y_bound1, &width; 
	uint16_t *pts, &num_pts;} square;

float APMI(const std::vector<float>& vec_x, const std::vector<float>& vec_y, const float q_thresh = 7.815, 
		const uint16_t size_thresh = 4);

std::vector<edge_tar> genemapAPMI(genemap &matrix, const gene_id_t& identifier, const float& q_thresh = 7.815, const uint16_t& size_thresh = 4);

const std::vector<float> permuteAPMI(std::vector<float> &ref_vec,
		std::vector<std::vector<float>> &target_vec, const float q_thresh = 7.815,
		const uint16_t size_thresh = 4);


//--------------------- NullModel.cpp	 		-----------------------
// sets a file static variable of an ordered null distribution
const std::vector<float> initNullMIs(const uint16_t& tot_num_samps);

// the function below requires that initNullMIs has been called
const float getMIPVal(const float& mi, const float& p_precise = 0.001f);

// the function below requires that initNullMIs has been called
const std::vector<float> getMIPVals(const std::vector<float>& mis, const float& p_precise = 0.001f);


//--------------------- AlphaPruning.cpp	 	-----------------------
std::pair<reg_web, map_map> pruneAlpha(reg_web &network, uint32_t& size_of_network);

//--------------------- MaxEntPruning.cpp	 	-----------------------
reg_web pruneMaxEnt(reg_web &network, map_map& tftfNetwork, uint32_t &size_of_network);

//--------------------- RegWebFns.cpp	 		-----------------------
reg_web sort_edge_tars(reg_web &regweb);	

/*
 Used for MaxEnt pruning to copy ARACNe-AP; we should probably either overhaul our own data structure or find a way to use our own, as any conversion is essentially a memory and runtime cost.
 */
std::unordered_map<gene_id_t, std::unordered_map<gene_id_t, float>> regweb_to_mapmap(reg_web &network);
//--------------------- Consolidator.cpp	 	-----------------------
std::vector<consolidated_df> consolidate(std::vector<reg_web> &subnets);


//--------------------- ARACNe3.cpp	 		-----------------------


#endif /* #ifndef _ARACNE3_HPP_ */
