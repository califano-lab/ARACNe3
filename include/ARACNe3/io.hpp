#ifndef io_hpp
#define io_hpp

#include "MINDy3.hpp"
#include <string>
#include <vector>
#include <fstream>
#include <random>

std::string makeUnixDirectoryNameUniversal(std::string dir_name);
void makeDir(std::string dir_name);
std::pair<hash_float_vec, hash_short_vec> readExpMatrixAndCopulaTransform(std::string &filename, std::mt19937 &rand);
hash_network readNetworkFile(std::string &filename);
std::set<gene_id> readModList(std::string &filename);
int outputModulatorDf(const std::vector<modulator_df> &output_df, std::ofstream &ofs);

#endif /* io_h */
