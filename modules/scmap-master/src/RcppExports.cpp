// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// subdistsmult
arma::mat subdistsmult(const Rcpp::List& subcentroids, const Rcpp::List& query_chunks, const int& M, const int& k, const int& cellnum);
RcppExport SEXP _scmap_subdistsmult(SEXP subcentroidsSEXP, SEXP query_chunksSEXP, SEXP MSEXP, SEXP kSEXP, SEXP cellnumSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type subcentroids(subcentroidsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type query_chunks(query_chunksSEXP);
    Rcpp::traits::input_parameter< const int& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int& >::type cellnum(cellnumSEXP);
    rcpp_result_gen = Rcpp::wrap(subdistsmult(subcentroids, query_chunks, M, k, cellnum));
    return rcpp_result_gen;
END_RCPP
}
// normalise
arma::mat normalise(const arma::mat& dat);
RcppExport SEXP _scmap_normalise(SEXP datSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type dat(datSEXP);
    rcpp_result_gen = Rcpp::wrap(normalise(dat));
    return rcpp_result_gen;
END_RCPP
}
// NN
Rcpp::List NN(const int& w, const int& k, const Rcpp::List& subcentroids, const arma::mat& subclusters, const Rcpp::List& query_chunks, const int& M, const arma::vec& SqNorm);
RcppExport SEXP _scmap_NN(SEXP wSEXP, SEXP kSEXP, SEXP subcentroidsSEXP, SEXP subclustersSEXP, SEXP query_chunksSEXP, SEXP MSEXP, SEXP SqNormSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type subcentroids(subcentroidsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type subclusters(subclustersSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type query_chunks(query_chunksSEXP);
    Rcpp::traits::input_parameter< const int& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type SqNorm(SqNormSEXP);
    rcpp_result_gen = Rcpp::wrap(NN(w, k, subcentroids, subclusters, query_chunks, M, SqNorm));
    return rcpp_result_gen;
END_RCPP
}
// EuclSqNorm
Rcpp::NumericVector EuclSqNorm(const arma::mat& dat);
RcppExport SEXP _scmap_EuclSqNorm(SEXP datSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type dat(datSEXP);
    rcpp_result_gen = Rcpp::wrap(EuclSqNorm(dat));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_scmap_subdistsmult", (DL_FUNC) &_scmap_subdistsmult, 5},
    {"_scmap_normalise", (DL_FUNC) &_scmap_normalise, 1},
    {"_scmap_NN", (DL_FUNC) &_scmap_NN, 7},
    {"_scmap_EuclSqNorm", (DL_FUNC) &_scmap_EuclSqNorm, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_scmap(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
