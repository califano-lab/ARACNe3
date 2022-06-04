/*
 A reg_web is an unordered hash table that maps a regulator identifier to a vector of edges and targets.  An "edge" is the mutual information value associated with the regulator and target.  Therefore, it is convenient to think of a reg_web as a regulator mapped to all its targets.  This file is a series of functions intended to simplify extracting information from that reg_web
 */

#include "ARACNe3.hpp"

extern uint16_t tot_num_regulators;
extern uint16_t tot_num_samps;

/*
 Simply returns a vector of the same regID duplicated for every edge it has associated
 */
std::vector<reg_id_t> get_regulator_dupe(reg_id_t regID, reg_web &regweb) {
	std::vector<edge_tar> &edges = regweb[regID];
	return std::vector<reg_id_t>(edges.size(), regID);
}

/*
 Returns a vector of the targetIDs in the regweb
 */
std::vector<reg_id_t> get_targets(reg_id_t regID, reg_web &regweb) {
	std::vector<edge_tar> &edges = regweb[regID];
	std::vector<reg_id_t> tars;
	tars.reserve(edges.size());
	for (auto &et : edges) {
		tars.emplace_back(et.target);
	}
	return tars;
}

/*
 Returns a vector of the MIs in the regweb
 */
std::vector<float> get_MIs(reg_id_t regID, reg_web &regweb) {
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
std::vector<float> get_p_vals(reg_id_t regID, reg_web &regweb) {
	initNullMIs(tot_num_samps);
	const std::vector<float> mis = get_MIs(regID, regweb);
	return getMIPVals(mis);
}

/*
 Sorts all the edge_tar vectors in a reg_web based on the target identifier, from smallest to largest.
 */
reg_web sort_edge_tars(reg_web &regweb) {
	for (reg_id_t i = 0; i < tot_num_regulators; ++i) {
		std::sort(regweb[i].begin(), regweb[i].end());
	}
	return regweb;
}

/*
 Takes a reg_web and returns a hash map of a hash map (for MaxEnt pruning step)
 */
map_map regweb_to_mapmap(reg_web &network) {
	map_map mapmap;
	for (reg_id_t reg = 0; reg < tot_num_regulators; ++reg) {
		for (edge_tar &et : network[reg])
			mapmap[reg][et.target] = et.mi;
	}
	return mapmap;
}
