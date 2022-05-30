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
 * Compressed and non-compressed versions.
 */
typedef struct {const std::string &target; const float mi;} edge_tar;
typedef struct {edge_tar &etr; const float p_value;} edge_tar_p;
typedef struct {const uint16_t &target; const float mi;} edge_tar_compressed;
typedef struct {edge_tar_compressed &etr; const float p_value;}
edge_tar_p_compressed;

/*
 * The regulator and target list is represented by this data type, an unordered
 * hash map between regulator names (string or number in uncompressed/compressed
 * implementation) as well as a vector of "ball-on-stick" structs, or all the
 * targets and their MIs.
 */
typedef std::unordered_map<std::string, std::vector<edge_tar>> reg_web;
typedef std::unordered_map<std::string, std::vector<edge_tar_p>> reg_web_p;
typedef std::unordered_map<uint16_t, std::vector<edge_tar_compressed>>
reg_web_compressed;
typedef std::unordered_map<uint16_t, std::vector<edge_tar_p_compressed>>
reg_web_p_compressed;

//--------------------- MatrixReglistIO.cpp 		-----------------------
/*
 * Maps strings to gene expression matrices
 */
typedef std::unordered_map<std::string, std::vector<float>> genemap;

/*
 * Compressed version of the original genemap, limited to 65536 possible feature
 * names.  Same as genemap, but swapping strings with unsigned short.
 */
typedef std::unordered_map<uint16_t, std::vector<float>>
genemap_compressed;

/*
 * Maps strings to unsigned short, mainly used for speeding up compression
 * feature
 */
typedef std::unordered_map<std::string, uint16_t> string_map;

std::vector<std::string> readRegList(std::string);

genemap readTransformedGexpMatrix(std::string);

// NOTE TO SELF: const on a std::object will make the entire object const, so do
// not define const within the object.  Returning const stuff is okay and
// automatically cast.  const tuple can have const contents.
const std::tuple<const std::vector<std::string>, const string_map>
readRegList_compressed(std::string);

const std::tuple<const genemap_compressed, const std::vector<std::string>>
readTransformedGexpMatrix_compressed(const std::tuple<const
		std::vector<std::string>, const string_map>, std::string);



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

void genemapAPMI(genemap &, const std::string &, const float, const unsigned
		short);

const std::vector<const float> permuteAPMI(std::vector<float> &ref_vec, 
		std::vector<std::vector<float>> &target_vec, const float,
		const uint16_t);


//--------------------- NullModel.cpp	 		-----------------------
std::vector<float> initNullMis(const unsigned int);



#endif /* #ifndef _ARACNE3_HPP_ */
