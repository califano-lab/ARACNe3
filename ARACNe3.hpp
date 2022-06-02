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
#include <numeric>
#include <math.h>
#include <filesystem>
#include <chrono>

#ifndef _ARACNE3_HPP_
#define _ARACNE3_HPP_


//--------------------- Overhead Data Structures	 -----------------------
/*
 Switching to this type alias to reduce confusion in other parts of the program
 */
typedef uint16_t reg_id_t;

/*
 * The "ball on stick" data structure for edges in network.  Will be used with
 * dictionary so that a regulator can have a group of these associated.
 */
typedef struct edge_tar{
	reg_id_t target;
	float mi;
	edge_tar(const reg_id_t &t, const float &mi) : target(t), mi(mi) {};
	
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

/*
 * The regulator and target list is represented by this data type, an unordered
 * hash map between regulator names (string or number in uncompressed/compressed
 * implementation) as well as a vector of "ball-on-stick" structs, or all the
 * targets and their MIs.
 */
typedef std::unordered_map<reg_id_t, std::vector<edge_tar>> reg_web;
typedef std::unordered_map<reg_id_t, std::unordered_map<reg_id_t, float>> map_map;


//--------------------- MatrixReglistIO.cpp 		-----------------------
/*
 * Maps gene identifiers to gene expression matrices
 */
typedef std::unordered_map<reg_id_t, std::vector<float>> genemap;

/*
 Does not return a list of regulators as a string vector, as we are using compression
 */
void readRegList(std::string filename = "regulators.txt");

/*
 Returns a map of gene identifier -> gene expression.
 */
genemap readTransformedGexpMatrix(std::string filename = "exp_mat.txt");

void printNetworkRegTarMI(const reg_web &network, const std::string &filename = "output.txt");


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

float APMI(std::vector<float>, std::vector<float>, const float q_thresh = 7.815, 
		const uint16_t size_thresh = 4);

std::vector<edge_tar> genemapAPMI(genemap &matrix, reg_id_t identifier, const float q_thresh = 7.815, const uint16_t size_thresh = 4);

const std::vector<float> permuteAPMI(std::vector<float> &ref_vec,
		std::vector<std::vector<float>> &target_vec, const float q_thresh = 7.815,
		const uint16_t size_thresh = 4);


//--------------------- NullModel.cpp	 		-----------------------
void checkInitNullMIs();

// sets a file static variable of an ordered null distribution
const std::vector<float> initNullMIs(uint16_t);

// the function below requires that initNullMIs has been called
const float getMIPVal(const float &);

// the function below requires that initNullMIs has been called
const std::vector<float> getMIPVals(const std::vector<float> &);


//--------------------- FDRPruning.cpp	 		-----------------------
reg_web pruneFDR(reg_web &network, uint32_t network_size, float FDR = 0.05);

//--------------------- MaxEntPruning.cpp	 		-----------------------
reg_web pruneMaxEnt(reg_web &network);

//--------------------- RegWebFns.cpp	 		-----------------------
reg_web sort_edge_tars(reg_web &regweb);

/*
 Used for MaxEnt pruning to copy ARACNe-AP; we should probably either overhaul our own data structure or find a way to use our own, as any conversion is essentially a memory and runtime cost.
 */
std::unordered_map<reg_id_t, std::unordered_map<reg_id_t, float>> regweb_to_mapmap(reg_web &network);

//--------------------- ARACNe3.cpp	 		-----------------------

#endif /* #ifndef _ARACNE3_HPP_ */
