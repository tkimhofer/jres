// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// convolve3cpp
arma::mat convolve3cpp(arma::mat A, arma::mat B);
RcppExport SEXP _jres_convolve3cpp(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(convolve3cpp(A, B));
    return rcpp_result_gen;
END_RCPP
}
// d2gauss_cpp
Rcpp::NumericVector d2gauss_cpp(NumericVector x, NumericVector y, NumericVector mu, NumericVector sigma);
RcppExport SEXP _jres_d2gauss_cpp(SEXP xSEXP, SEXP ySEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(d2gauss_cpp(x, y, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// expandGrid_rcpp
NumericMatrix expandGrid_rcpp(NumericVector x, NumericVector y);
RcppExport SEXP _jres_expandGrid_rcpp(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(expandGrid_rcpp(x, y));
    return rcpp_result_gen;
END_RCPP
}
// matsymminf2f1
NumericMatrix matsymminf2f1(NumericMatrix A, int idx_cent_row, int idx_cent_col);
RcppExport SEXP _jres_matsymminf2f1(SEXP ASEXP, SEXP idx_cent_rowSEXP, SEXP idx_cent_colSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< int >::type idx_cent_row(idx_cent_rowSEXP);
    Rcpp::traits::input_parameter< int >::type idx_cent_col(idx_cent_colSEXP);
    rcpp_result_gen = Rcpp::wrap(matsymminf2f1(A, idx_cent_row, idx_cent_col));
    return rcpp_result_gen;
END_RCPP
}
// matsymminf2
NumericMatrix matsymminf2(NumericMatrix A, int idx_cent);
RcppExport SEXP _jres_matsymminf2(SEXP ASEXP, SEXP idx_centSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< int >::type idx_cent(idx_centSEXP);
    rcpp_result_gen = Rcpp::wrap(matsymminf2(A, idx_cent));
    return rcpp_result_gen;
END_RCPP
}
// getBoundary
int getBoundary(int i, int add_int, int mat_end, char type);
RcppExport SEXP _jres_getBoundary(SEXP iSEXP, SEXP add_intSEXP, SEXP mat_endSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type add_int(add_intSEXP);
    Rcpp::traits::input_parameter< int >::type mat_end(mat_endSEXP);
    Rcpp::traits::input_parameter< char >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(getBoundary(i, add_int, mat_end, type));
    return rcpp_result_gen;
END_RCPP
}
// matsum
NumericMatrix matsum(NumericMatrix A, NumericMatrix B);
RcppExport SEXP _jres_matsum(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(matsum(A, B));
    return rcpp_result_gen;
END_RCPP
}
// matpow
NumericMatrix matpow(NumericMatrix A, const double exp);
RcppExport SEXP _jres_matpow(SEXP ASEXP, SEXP expSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< const double >::type exp(expSEXP);
    rcpp_result_gen = Rcpp::wrap(matpow(A, exp));
    return rcpp_result_gen;
END_RCPP
}
// matsumup
double matsumup(NumericMatrix A);
RcppExport SEXP _jres_matsumup(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(matsumup(A));
    return rcpp_result_gen;
END_RCPP
}
// minmaxMat
NumericMatrix minmaxMat(NumericMatrix A);
RcppExport SEXP _jres_minmaxMat(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(minmaxMat(A));
    return rcpp_result_gen;
END_RCPP
}
// getBBdim_rcpp
IntegerVector getBBdim_rcpp(NumericVector ppm, int i, double add);
RcppExport SEXP _jres_getBBdim_rcpp(SEXP ppmSEXP, SEXP iSEXP, SEXP addSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type ppm(ppmSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< double >::type add(addSEXP);
    rcpp_result_gen = Rcpp::wrap(getBBdim_rcpp(ppm, i, add));
    return rcpp_result_gen;
END_RCPP
}
// lapOfG
List lapOfG(NumericMatrix sub, NumericVector f1hz, NumericVector f2hz, double cent_f1hz, double cent_f2hz, double sf, NumericVector npix);
RcppExport SEXP _jres_lapOfG(SEXP subSEXP, SEXP f1hzSEXP, SEXP f2hzSEXP, SEXP cent_f1hzSEXP, SEXP cent_f2hzSEXP, SEXP sfSEXP, SEXP npixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type sub(subSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type f1hz(f1hzSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type f2hz(f2hzSEXP);
    Rcpp::traits::input_parameter< double >::type cent_f1hz(cent_f1hzSEXP);
    Rcpp::traits::input_parameter< double >::type cent_f2hz(cent_f2hzSEXP);
    Rcpp::traits::input_parameter< double >::type sf(sfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type npix(npixSEXP);
    rcpp_result_gen = Rcpp::wrap(lapOfG(sub, f1hz, f2hz, cent_f1hz, cent_f2hz, sf, npix));
    return rcpp_result_gen;
END_RCPP
}
// pickPeaks_rcpp
List pickPeaks_rcpp(NumericMatrix jr, NumericVector f1hz, NumericVector f2ppm, double noise, double boundary, double sf);
RcppExport SEXP _jres_pickPeaks_rcpp(SEXP jrSEXP, SEXP f1hzSEXP, SEXP f2ppmSEXP, SEXP noiseSEXP, SEXP boundarySEXP, SEXP sfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type jr(jrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type f1hz(f1hzSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type f2ppm(f2ppmSEXP);
    Rcpp::traits::input_parameter< double >::type noise(noiseSEXP);
    Rcpp::traits::input_parameter< double >::type boundary(boundarySEXP);
    Rcpp::traits::input_parameter< double >::type sf(sfSEXP);
    rcpp_result_gen = Rcpp::wrap(pickPeaks_rcpp(jr, f1hz, f2ppm, noise, boundary, sf));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_jres_convolve3cpp", (DL_FUNC) &_jres_convolve3cpp, 2},
    {"_jres_d2gauss_cpp", (DL_FUNC) &_jres_d2gauss_cpp, 4},
    {"_jres_expandGrid_rcpp", (DL_FUNC) &_jres_expandGrid_rcpp, 2},
    {"_jres_matsymminf2f1", (DL_FUNC) &_jres_matsymminf2f1, 3},
    {"_jres_matsymminf2", (DL_FUNC) &_jres_matsymminf2, 2},
    {"_jres_getBoundary", (DL_FUNC) &_jres_getBoundary, 4},
    {"_jres_matsum", (DL_FUNC) &_jres_matsum, 2},
    {"_jres_matpow", (DL_FUNC) &_jres_matpow, 2},
    {"_jres_matsumup", (DL_FUNC) &_jres_matsumup, 1},
    {"_jres_minmaxMat", (DL_FUNC) &_jres_minmaxMat, 1},
    {"_jres_getBBdim_rcpp", (DL_FUNC) &_jres_getBBdim_rcpp, 3},
    {"_jres_lapOfG", (DL_FUNC) &_jres_lapOfG, 7},
    {"_jres_pickPeaks_rcpp", (DL_FUNC) &_jres_pickPeaks_rcpp, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_jres(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
