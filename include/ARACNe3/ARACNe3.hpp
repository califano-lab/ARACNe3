#pragma once

#include <string>
#include <unordered_set>
#include <unordered_map>
#include <vector>

/*
 Switching to this type alias to reduce confusion in other parts of the program
 */
typedef uint16_t gene_id;
typedef std::unordered_set<gene_id> geneset;

// Maps gene to regulon
typedef std::unordered_map<gene_id, geneset> gene_to_geneset;

// Used for gexp and ranked gexp matrix storage
typedef std::vector<std::vector<float>> gene_to_floats;
typedef std::vector<std::vector<uint16_t>> gene_to_shorts;

// used for network storage
typedef std::unordered_map<gene_id, float> gene_to_float;
typedef std::unordered_map<gene_id, gene_to_float>
    gene_to_gene_to_float;

using vv_float = std::vector<std::vector<float>>;
