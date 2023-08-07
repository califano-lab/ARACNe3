#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <tuple>
#include <set>
#include <unordered_set>
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

// Maps gene to regulon
typedef std::unordered_map<gene_id, std::unordered_set<gene_id>>
    gene_to_geneset;

// Used for gexp and ranked gexp matrix storage
typedef std::unordered_map<gene_id, std::vector<float>> gene_to_floats;
typedef std::unordered_map<gene_id, std::vector<uint16_t>> gene_to_shorts;

// used for edge strength
typedef std::unordered_map<gene_id, std::unordered_map<gene_id, float>>
    gene_to_gene_to_float;
typedef std::unordered_map<gene_id, std::unordered_map<gene_id, float>>
    gene_to_gene_to_float;

//--------------------- NullModel.cpp	 		-----------------------
// sets a file static variable of an ordered null distribution
const std::vector<float> initNullMIs(const uint16_t& tot_num_subsample);

// the function below requires that initNullMIs has been called
const float getMIPVal(const float& mi, const float& p_precise = 0.001f);

// the function below requires that initNullMIs has been called
const std::vector<float> getMIPVals(const std::vector<float>& mis, const float& p_precise = 0.001f);

