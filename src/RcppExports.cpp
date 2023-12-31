// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// SigCoordPiso
NumericMatrix SigCoordPiso(NumericMatrix coord, int resolutionValue, double gamaValue, double increase, double contrastValue, double zoomValue, int numberCoord);
RcppExport SEXP _levi_SigCoordPiso(SEXP coordSEXP, SEXP resolutionValueSEXP, SEXP gamaValueSEXP, SEXP increaseSEXP, SEXP contrastValueSEXP, SEXP zoomValueSEXP, SEXP numberCoordSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< int >::type resolutionValue(resolutionValueSEXP);
    Rcpp::traits::input_parameter< double >::type gamaValue(gamaValueSEXP);
    Rcpp::traits::input_parameter< double >::type increase(increaseSEXP);
    Rcpp::traits::input_parameter< double >::type contrastValue(contrastValueSEXP);
    Rcpp::traits::input_parameter< double >::type zoomValue(zoomValueSEXP);
    Rcpp::traits::input_parameter< int >::type numberCoord(numberCoordSEXP);
    rcpp_result_gen = Rcpp::wrap(SigCoordPiso(coord, resolutionValue, gamaValue, increase, contrastValue, zoomValue, numberCoord));
    return rcpp_result_gen;
END_RCPP
}
// matrix_entrada
List matrix_entrada(NumericMatrix coordPiso, NumericMatrix SignalOut, NumericMatrix signalExp, NumericMatrix signalCtrl, NumericMatrix coord, int resolutionValue, double increase, double zoomValue, int numberCoord);
RcppExport SEXP _levi_matrix_entrada(SEXP coordPisoSEXP, SEXP SignalOutSEXP, SEXP signalExpSEXP, SEXP signalCtrlSEXP, SEXP coordSEXP, SEXP resolutionValueSEXP, SEXP increaseSEXP, SEXP zoomValueSEXP, SEXP numberCoordSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type coordPiso(coordPisoSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type SignalOut(SignalOutSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type signalExp(signalExpSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type signalCtrl(signalCtrlSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< int >::type resolutionValue(resolutionValueSEXP);
    Rcpp::traits::input_parameter< double >::type increase(increaseSEXP);
    Rcpp::traits::input_parameter< double >::type zoomValue(zoomValueSEXP);
    Rcpp::traits::input_parameter< int >::type numberCoord(numberCoordSEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_entrada(coordPiso, SignalOut, signalExp, signalCtrl, coord, resolutionValue, increase, zoomValue, numberCoord));
    return rcpp_result_gen;
END_RCPP
}
// matrix_saida
List matrix_saida(NumericMatrix matrixIn, int resolutionValue, double smoothValue, double gamaValue, double increase, double zoomValue, int h);
RcppExport SEXP _levi_matrix_saida(SEXP matrixInSEXP, SEXP resolutionValueSEXP, SEXP smoothValueSEXP, SEXP gamaValueSEXP, SEXP increaseSEXP, SEXP zoomValueSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type matrixIn(matrixInSEXP);
    Rcpp::traits::input_parameter< int >::type resolutionValue(resolutionValueSEXP);
    Rcpp::traits::input_parameter< double >::type smoothValue(smoothValueSEXP);
    Rcpp::traits::input_parameter< double >::type gamaValue(gamaValueSEXP);
    Rcpp::traits::input_parameter< double >::type increase(increaseSEXP);
    Rcpp::traits::input_parameter< double >::type zoomValue(zoomValueSEXP);
    Rcpp::traits::input_parameter< int >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_saida(matrixIn, resolutionValue, smoothValue, gamaValue, increase, zoomValue, h));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_levi_SigCoordPiso", (DL_FUNC) &_levi_SigCoordPiso, 7},
    {"_levi_matrix_entrada", (DL_FUNC) &_levi_matrix_entrada, 9},
    {"_levi_matrix_saida", (DL_FUNC) &_levi_matrix_saida, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_levi(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
