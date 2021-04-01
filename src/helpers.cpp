#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
NumericMatrix CoxHessianHelper(const NumericMatrix x, const NumericVector multiplier, const NumericVector eta, const IntegerVector order_b1, const IntegerVector rankmin_b1)
{
    int n = x.nrow();
    int p = x.ncol();
    if(multiplier.size() != n || order_b1.size() != n || rankmin_b1.size() != n || eta.size() != n){
        stop("Input dimensions does not match.");
    }

    NumericMatrix result(n, p);

    // row-wise reverse cumulative sum
    for(int j = 0; j < p; ++j){
        double sum = 0;
        for(int i = n - 1; i >= 0; --i){
            int row_ind = order_b1[i] - 1; // base one
            sum += x(row_ind, j) * eta[i];
            result(row_ind, j) = sum;
        }
    }

    // Adjust for ties
    for(int i = 0; i < n; ++i){
        int row_ind = order_b1[i] - 1; // base one
        int rankmin_ind = order_b1[rankmin_b1[i] - 1] - 1; // base one
        // double mult = multiplier[i];
        for(int j = 0; j < p; ++j){
            result(row_ind, j) = result(rankmin_ind, j);
        }
    }

    // Get the multiplier
    for(int i = 0; i < n; ++i){
        int row_ind = order_b1[i] - 1; // base one
        double mult = multiplier[i];
        for(int j = 0; j < p; ++j){
            result(row_ind, j) *= mult;
        }
    }
    return result;
}

//' @export
// [[Rcpp::export]]
SEXP eigenCrossProd(const Eigen::Map<Eigen::MatrixXd> A)
{
    // Return A^T A 
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(A.cols(), A.cols());
    result.selfadjointView<Eigen::Lower>().rankUpdate(A.transpose());
    result.triangularView<Eigen::Upper>() = result.transpose();
    return Rcpp::wrap(result);
}