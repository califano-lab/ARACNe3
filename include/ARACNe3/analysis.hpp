#pragma once

#include <string>

#include "aracne3io.hpp"
#include "logger.hpp"

int ARACNe3Analysis(const ARACNe3IOHandler &io, const std::string &version,
             const std::string &runid, const uint32_t seed,
             const uint8_t threads, const bool verbose, const float alpha,
             const float subsamp_pct, const uint32_t n_nulls,
             const std::string &method, const bool prune_alpha,
             const bool prune_MaxEnt, const bool save_subnets,
             const std::string &subnets_log_dir, const bool adaptive,
             const uint16_t min_subnets, const uint16_t max_subnets,
             const uint16_t min_regulon_occpuancy, const bool consolidate_mode,
             const bool skip_consolidate, Logger *const aracne3_logger,
             const std::string &subnets_dir, const std::string &cached_dir);
