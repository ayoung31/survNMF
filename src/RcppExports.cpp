// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// f_grad_rcpp
arma::mat f_grad_rcpp(arma::mat& X, arma::mat& W, arma::mat& H, arma::mat& beta, float alpha, arma::mat& y, arma::mat& delta, float theta, int j);
RcppExport SEXP _survNMF_f_grad_rcpp(SEXP XSEXP, SEXP WSEXP, SEXP HSEXP, SEXP betaSEXP, SEXP alphaSEXP, SEXP ySEXP, SEXP deltaSEXP, SEXP thetaSEXP, SEXP jSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type H(HSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< float >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< float >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    rcpp_result_gen = Rcpp::wrap(f_grad_rcpp(X, W, H, beta, alpha, y, delta, theta, j));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_survNMF_f_grad_rcpp", (DL_FUNC) &_survNMF_f_grad_rcpp, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_survNMF(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
