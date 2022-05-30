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

#ifndef _ARACNE3_HPP_
#define _ARACNE3_HPP_


//--------------------- Overhead Data Structures	 -----------------------
/*
 * The "ball on stick" data structure for edges in network.  Will be used with
 * dictionary so that a regulator can have a group of these associated.
 */
typedef struct edge_tar{
	const std::string &target;
	const float mi;
	edge_tar(const std::string &t, const float mi) : target(t), mi(mi) {};
} edge_tar;

typedef struct edge_tar_p : edge_tar {
	const float p_value;
	edge_tar_p(const std::string &t, const float mi, const float p) :
	edge_tar(t, mi), p_value(p) {};
} edge_tar_p;

/*
 * The regulator and target list is represented by this data type, an unordered
 * hash map between regulator names (string or number in uncompressed/compressed
 * implementation) as well as a vector of "ball-on-stick" structs, or all the
 * targets and their MIs.
 */
typedef std::unordered_map<std::string, std::vector<edge_tar>> reg_web;
typedef std::unordered_map<std::string, std::vector<edge_tar_p>> reg_web_p;

//--------------------- MatrixReglistIO.cpp 		-----------------------
/*
 * Maps strings to gene expression matrices
 */
typedef std::unordered_map<std::string, std::vector<float>> genemap;

std::vector<std::string> readRegList(std::string);

genemap readTransformedGexpMatrix(std::string);


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

float APMI(std::vector<float>, std::vector<float>, const float q_thresh, 
		const uint16_t);

std::vector<edge_tar> genemapAPMI(genemap &, const std::string &, const float q_thresh = 7.815, const uint16_t size_thresh = 4);

const std::vector<float> permuteAPMI(std::vector<float> &ref_vec, 
		std::vector<std::vector<float>> &target_vec, const float,
		const uint16_t);


//--------------------- NullModel.cpp	 		-----------------------
// sets a file static variable of an ordered null distribution
std::vector<float> initNullMIs(const unsigned int);

// the function below requires that initNullMIs has been called
float getMIPVal(const float);

// the function below requires that initNullMIs has been called
std::vector<float> getMIPVals(const std::vector<float>);



#endif /* #ifndef _ARACNE3_HPP_ */
