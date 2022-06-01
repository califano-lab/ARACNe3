/*
 A reg_web is an unordered hash table that maps a regulator identifier to a vector of edges and targets.  An "edge" is the mutual information value associated with the regulator and target.  Therefore, it is convenient to think of a reg_web as a regulator mapped to all its targets.  This file is a series of functions intended to simplify extracting information from that reg_web
 */

#include "ARACNe3.hpp"
namespace regweb {
/*
 Simply returns a vector of the same regID duplicated for every edge it has associated
 */
std::vector<uint16_t> get_regulator_dupe(uint16_t regID, reg_web &regweb) {
	std::vector<edge_tar> &edges = regweb[regID];
	return std::vector<uint16_t>(edges.size(), regID);
}

/*
 Returns a vector of the targetIDs in the regweb
 */
std::vector<uint16_t> get_targets(uint16_t regID, reg_web &regweb) {
	std::vector<edge_tar> &edges = regweb[regID];
	std::vector<uint16_t> tars;
	tars.reserve(edges.size());
	for (auto &et : edges) {
		tars.emplace_back(et.target);
	}
	return tars;
}

/*
 Returns a vector of the MIs in the regweb
 */
std::vector<float> get_MIs(uint16_t regID, reg_web &regweb) {
	std::vector<edge_tar> &edges = regweb[regID];
	std::vector<float> mis;
	mis.reserve(edges.size());
	for (auto &et : edges) {
		mis.emplace_back(et.mi);
	}
	return mis;
}

/*
 Returns a vector of the p-values in the regweb, must be computed as this is not reg_web_p
 */
std::vector<float> get_p_vals(uint16_t regID, reg_web &regweb) {
	checkInitNullMIs();
	const std::vector<float> mis = get_MIs(regID, regweb);
	return getMIPVals(mis);
}
}
