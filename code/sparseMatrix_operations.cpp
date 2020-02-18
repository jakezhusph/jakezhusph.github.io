#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>

using namespace std;
using namespace arma;
using namespace Rcpp;


//------------------------------------
// rowMeans, colMeans of a sparse matrix 
//------------------------------------

// [[Rcpp::export]]
arma::rowvec sp_means(arma::sp_mat sp_data, bool rowMeans = false) {

  arma::sp_mat norm_col_sums;

  arma::mat tmp_mat;

  if (rowMeans) {

    norm_col_sums = arma::mean(sp_data, 1);

    tmp_mat = arma::conv_to< arma::mat >::from(norm_col_sums.col(0));}

  else {

    norm_col_sums = arma::mean(sp_data, 0);

    tmp_mat = arma::conv_to< arma::mat >::from(norm_col_sums.row(0));
  }

  arma::rowvec tmp_vec = arma::conv_to< arma::rowvec >::from(tmp_mat);

  return tmp_vec;
}



//------------------------------------
// rowSums, colSums of a sparse matrix
//------------------------------------


// [[Rcpp::export]]
arma::rowvec sp_sums(arma::sp_mat sp_data, bool rowSums = false) {

  arma::mat tmp_mat;

  arma::sp_mat norm_col_sums;

  if (rowSums) {

    norm_col_sums = arma::sum(sp_data, 1);

    tmp_mat = arma::conv_to< arma::mat >::from(norm_col_sums.col(0));}

  else {

    norm_col_sums = arma::sum(sp_data, 0);

    tmp_mat = arma::conv_to< arma::mat >::from(norm_col_sums.row(0));
  }

  arma::rowvec tmp_vec = arma::conv_to< arma::rowvec >::from(tmp_mat);

  return tmp_vec;
}



//------------------------------------
// Count non-zero elements of sparse matrix
//------------------------------------


// [[Rcpp::export]]
arma::rowvec sp_nz_count(arma::sp_mat sp_data, bool rowSums = false) {

  arma::sp_mat sp_data_binary = spones(sp_data);
  arma::mat tmp_mat;

  arma::sp_mat norm_col_sums;

  if (rowSums) {

    norm_col_sums = arma::sum(sp_data_binary, 1);

    tmp_mat = arma::conv_to< arma::mat >::from(norm_col_sums.col(0));}

  else {

    norm_col_sums = arma::sum(sp_data_binary, 0);

    tmp_mat = arma::conv_to< arma::mat >::from(norm_col_sums.row(0));
  }

  arma::rowvec tmp_vec = arma::conv_to< arma::rowvec >::from(tmp_mat);

  return tmp_vec;
}



//------------------------------------
// rowVariance, colVariance of a sparse matrix
//------------------------------------


// [[Rcpp::export]]
arma::rowvec sp_vars(arma::sp_mat sp_data, bool rowVars = false) {

  arma::sp_mat norm_col_Vars;

  arma::mat tmp_mat;

  if (rowVars) {

    norm_col_Vars = arma::var(sp_data,0, 1);

    tmp_mat = arma::conv_to< arma::mat >::from(norm_col_Vars.col(0));}

  else {

    norm_col_Vars = arma::var(sp_data, 0,0);

    tmp_mat = arma::conv_to< arma::mat >::from(norm_col_Vars.row(0));
  }

  arma::rowvec tmp_vec = arma::conv_to< arma::rowvec >::from(tmp_mat);

  return tmp_vec;
}




//------------------------------------
// For Fast KDC test statistics calculation
//------------------------------------


// [[Rcpp::export]]
double ComputeTestQuantRcpp(arma::vec gene_vec,double deno, arma::mat kmat0, arma::mat loc_mat){

  const int m = gene_vec.n_elem;

  double score = 0;
  arma::vec a = zeros<arma::vec>(m);
  arma::vec b = zeros<arma::vec>(m);

  for(int icell = 1; icell < m; ++icell) {
    a = gene_vec(icell)*gene_vec*deno;
    // Rcpp::Rcout << "error after here " << std::endl;
    b = loc_mat*kmat0.row(icell).t();
    double tmp = dot(a, b);
    score = score + tmp;
  } 
  return score;
}




