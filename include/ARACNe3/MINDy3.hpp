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
#include <numeric>
#include <math.h>
#include <filesystem>
#include <chrono>
#include <regex>
#include <iomanip> // necessary for CentOS
#include <omp.h>

#pragma once

#if defined __linux__ || defined __APPLE__
const std::string hiddenfpre = ".";
const char directory_slash = '/';
#elif defined _WIN32 // in windows you must set a hidden file via properties
const char directory_slash = '\\';
const std::string hiddenfpre = "";
#else
const std::string hiddenfpre = "";
#endif /* __linux__ || __APPLE__ */

//--------------------- Overhead Structures	 		-----------------------
typedef uint16_t gene_id;

// A particular gene's regulon
typedef std::unordered_map<gene_id, std::unordered_set<gene_id>> hash_network;

// Data structures for hashing to vectors or float values
typedef std::unordered_map<gene_id, std::vector<float>> hash_float_vec; // used explicitly for gexp matrix storage
typedef std::unordered_map<gene_id, std::vector<uint16_t>> hash_short_vec; // used explicitly for ranked gexp matrix storage

typedef std::unordered_map<gene_id, std::unordered_map<gene_id, float>> hash_edge_float; // used for mi/scc edge storage

typedef struct modulator_df {
	const gene_id regulator;
	const gene_id target;
	const gene_id modulator;
	const float reg_tar_mi;
	const float reg_mod_mi;
	const float mod_tar_mi;
	const float reg_tar_scc;
	const float reg_mod_scc;
	const float mod_tar_scc;
	const float cond_mi;
	modulator_df (const gene_id& r, const gene_id& t, const gene_id& m, const float& rtmi, const float& rmmi, const float& mtmi, const float& rtscc, const float& rmscc, const float& mtscc, const float& apcmi) : regulator(r), target(t), modulator(m), reg_tar_mi(rtmi), reg_mod_mi(rmmi), mod_tar_mi(mtmi), reg_tar_scc(rtscc), reg_mod_scc(rmscc), mod_tar_scc(mtscc), cond_mi(apcmi) {};
} modulator_df;
