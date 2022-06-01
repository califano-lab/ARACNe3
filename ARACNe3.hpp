#include <vector>
#include <string>
#include <tuple>
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
 * The "ball on stick" data structure for edges in network.  Will be used with
 * dictionary so that a regulator can have a group of these associated.
 */
typedef struct edge_tar{
	uint16_t target;
	float mi;
	edge_tar(const uint16_t &t, const float &mi) : target(t), mi(mi) {};
} edge_tar;

/*
 * The regulator and target list is represented by this data type, an unordered
 * hash map between regulator names (string or number in uncompressed/compressed
 * implementation) as well as a vector of "ball-on-stick" structs, or all the
 * targets and their MIs.
 */
typedef std::unordered_map<uint16_t, std::vector<edge_tar>> reg_web;

//--------------------- MatrixReglistIO.cpp 		-----------------------
/*
 * Maps gene identifiers to gene expression matrices
 */
typedef std::unordered_map<uint16_t, std::vector<float>> genemap;

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

std::vector<edge_tar> genemapAPMI(genemap &matrix, uint16_t identifier, const float q_thresh = 7.815, const uint16_t size_thresh = 4);

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

//--------------------- DPIPruning.cpp	 		-----------------------
reg_web pruneDPI(reg_web &network);

//--------------------- ARACNe3.cpp	 		-----------------------


#endif /* #ifndef _ARACNE3_HPP_ */
