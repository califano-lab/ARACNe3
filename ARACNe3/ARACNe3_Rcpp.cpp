/*
 This file is intended to contain functions callable from R that will allow for
 Rcpp usage and a custom analysis.
 
 This file is not ready for any usage or release.
 */

// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;

// [[Rcpp::export]]
//float APMI(vector<float> &vec_x, vector<float> &vec_y,
