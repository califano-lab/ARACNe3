#ifndef apcmi_nullmodel_hpp
#define apcmi_nullmodel_hpp

#include <vector>
#include <string> 
#include <random>
#include <filesystem>
#include <algorithm>
#include <fstream>
#include <iterator>

class APCMINullModel {
private:
	std::vector<float> null_mis;
	float m, b;
	std::string nulls_filename_no_extension, OLS_coefs_filename_no_extension;

public:
	APCMINullModel(const APCMINullModel& copied); //copy ctor
	// rand should be passed from main based on seed for predictable behavior.
	APCMINullModel(const uint16_t sample_size, const uint8_t n_bins, const uint32_t num_null_marginals, std::mt19937 &rand); // ctor for APCMINullModel generation
	APCMINullModel(const uint16_t sample_size, const uint8_t n_bins, const uint32_t num_null_marginals, std::mt19937 &rand, const std::string cached_dir); // ctor for cached APCMINullModel
	~APCMINullModel(); // dtor
	void cacheNullModel(const std::string cached_dir); // cache vec, m, and b
	const float getCMIPVal(const float& cmi, const float& p_precise = 0.001f); // return p value
};

#endif /* apcmi_nullmodel_hpp */
