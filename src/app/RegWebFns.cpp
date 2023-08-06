/*
 A gene_to_edge_tars is an unordered hash table that maps a regulator identifier to a vector of edges and targets.  An "edge" is the mutual information value associated with the regulator and target.  Therefore, it is convenient to think of a gene_to_edge_tars as a regulator mapped to all its targets.  This file is a series of functions intended to simplify extracting information from that gene_to_edge_tars
 */

#include "ARACNe3.hpp"

extern uint16_t tot_num_regulators;
extern uint16_t tot_num_subsample;
extern uint32_t num_null_marginals;

/*
 Simply returns a vector of the same regID duplicated for every edge it has associated
 */
std::vector<gene_id> get_regulator_dupe(gene_id regID, gene_to_edge_tars &regweb) {
	std::vector<edge_tar> &edges = regweb[regID];
	return std::vector<gene_id>(edges.size(), regID);
}

/*
 Returns a vector of the targetIDs in the regweb
 */
std::vector<gene_id> get_targets(gene_id regID, gene_to_edge_tars &regweb) {
	std::vector<edge_tar> &edges = regweb[regID];
	std::vector<gene_id> tars;
	tars.reserve(edges.size());
	for (auto &et : edges)
		tars.emplace_back(et.target);
	return tars;
}

/*
 Returns a vector of the MIs in the regweb
 */
std::vector<float> get_MIs(gene_id regID, gene_to_edge_tars &regweb) {
	std::vector<edge_tar> &edges = regweb[regID];
	std::vector<float> mis;
	mis.reserve(edges.size());
	for (auto &et : edges)
		mis.emplace_back(et.mi);
	return mis;
}

/*
 Returns a vector of the p-values in the regweb, must be computed as this is not gene_to_edge_tars_p
 */
std::vector<float> get_p_vals(gene_id regID, gene_to_edge_tars &regweb) {
	initNullMIs(tot_num_subsample);
	const std::vector<float> mis = get_MIs(regID, regweb);
	return getMIPVals(mis, 100.0/num_null_marginals);
}

/*
 Sorts all the edge_tar vectors in a gene_to_edge_tars based on the target identifier, from smallest to largest.
 */
gene_to_edge_tars sort_edge_tars(gene_to_edge_tars &regweb) {
	for (gene_id i = 0; i < tot_num_regulators; ++i)
		std::sort(regweb[i].begin(), regweb[i].end());
	return regweb;
}

/*
 Takes a gene_to_edge_tars and returns a hash map of a hash map (for MaxEnt pruning step)
 */
gene_to_gene_to_float regweb_to_mapmap(gene_to_edge_tars &network) {
	gene_to_gene_to_float mapmap;
	for (gene_id reg = 0; reg < tot_num_regulators; ++reg)
		for (edge_tar &et : network[reg])
			mapmap[reg][et.target] = et.mi;
	return mapmap;
}
