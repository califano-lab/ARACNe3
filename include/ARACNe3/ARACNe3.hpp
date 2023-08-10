#pragma once

#include <vector>
#include <unordered_set>
#include <string>

#if defined __linux__ || defined __APPLE__
const std::string hiddenfpre = ".";
const char directory_slash = '/';
#elif defined _WIN32 // in windows you must set a hidden file via properties
const char directory_slash = '\\';
const std::string hiddenfpre = "";
#else
const char directory_slash = '/';
const std::string hiddenfpre = "";
#endif /* __linux__ || __APPLE__ */


//--------------------- Overhead Data Structures	 -----------------------
/*
 Switching to this type alias to reduce confusion in other parts of the program
 */
typedef uint16_t gene_id;
typedef std::unordered_set<gene_id> geneset;

// Maps gene to regulon
typedef std::unordered_map<gene_id, geneset>
    gene_to_geneset;

// Used for gexp and ranked gexp matrix storage
typedef std::unordered_map<gene_id, std::vector<float>> gene_to_floats;
typedef std::unordered_map<gene_id, std::vector<uint16_t>> gene_to_shorts;

// used for edge strength
typedef std::unordered_map<gene_id, std::unordered_map<gene_id, float>>
    gene_to_gene_to_float;
typedef std::unordered_map<gene_id, std::unordered_map<gene_id, float>>
    gene_to_gene_to_float;
