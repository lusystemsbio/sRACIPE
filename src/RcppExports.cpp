// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// generateThresholds
size_t generateThresholds(const Rcpp::IntegerMatrix geneInteraction, Rcpp::NumericVector thresholdGene, Rcpp::List config);
RcppExport SEXP _sRACIPE_generateThresholds(SEXP geneInteractionSEXP, SEXP thresholdGeneSEXP, SEXP configSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerMatrix >::type geneInteraction(geneInteractionSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type thresholdGene(thresholdGeneSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type config(configSEXP);
    rcpp_result_gen = Rcpp::wrap(generateThresholds(geneInteraction, thresholdGene, config));
    return rcpp_result_gen;
END_RCPP
}
// limitcyclesGRC
int limitcyclesGRC(Rcpp::IntegerMatrix geneInteraction, String outFileLC, String outFileLCIC, Rcpp::List config, Rcpp::LogicalVector& modelConverg, String inFileParams, String inFileGE, Rcpp::NumericVector geneTypes, const int stepper);
RcppExport SEXP _sRACIPE_limitcyclesGRC(SEXP geneInteractionSEXP, SEXP outFileLCSEXP, SEXP outFileLCICSEXP, SEXP configSEXP, SEXP modelConvergSEXP, SEXP inFileParamsSEXP, SEXP inFileGESEXP, SEXP geneTypesSEXP, SEXP stepperSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type geneInteraction(geneInteractionSEXP);
    Rcpp::traits::input_parameter< String >::type outFileLC(outFileLCSEXP);
    Rcpp::traits::input_parameter< String >::type outFileLCIC(outFileLCICSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type config(configSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector& >::type modelConverg(modelConvergSEXP);
    Rcpp::traits::input_parameter< String >::type inFileParams(inFileParamsSEXP);
    Rcpp::traits::input_parameter< String >::type inFileGE(inFileGESEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type geneTypes(geneTypesSEXP);
    Rcpp::traits::input_parameter< const int >::type stepper(stepperSEXP);
    rcpp_result_gen = Rcpp::wrap(limitcyclesGRC(geneInteraction, outFileLC, outFileLCIC, config, modelConverg, inFileParams, inFileGE, geneTypes, stepper));
    return rcpp_result_gen;
END_RCPP
}
// simulateGRCCpp
int simulateGRCCpp(Rcpp::IntegerMatrix geneInteraction, Rcpp::List config, String outFileGE, String outFileParams, String outFileIC, String outFileConverge, Rcpp::NumericVector geneTypes, Rcpp::NumericMatrix signalVals, Rcpp::NumericVector signalingTypes, const int stepper);
RcppExport SEXP _sRACIPE_simulateGRCCpp(SEXP geneInteractionSEXP, SEXP configSEXP, SEXP outFileGESEXP, SEXP outFileParamsSEXP, SEXP outFileICSEXP, SEXP outFileConvergeSEXP, SEXP geneTypesSEXP, SEXP signalValsSEXP, SEXP signalingTypesSEXP, SEXP stepperSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type geneInteraction(geneInteractionSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type config(configSEXP);
    Rcpp::traits::input_parameter< String >::type outFileGE(outFileGESEXP);
    Rcpp::traits::input_parameter< String >::type outFileParams(outFileParamsSEXP);
    Rcpp::traits::input_parameter< String >::type outFileIC(outFileICSEXP);
    Rcpp::traits::input_parameter< String >::type outFileConverge(outFileConvergeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type geneTypes(geneTypesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type signalVals(signalValsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type signalingTypes(signalingTypesSEXP);
    Rcpp::traits::input_parameter< const int >::type stepper(stepperSEXP);
    rcpp_result_gen = Rcpp::wrap(simulateGRCCpp(geneInteraction, config, outFileGE, outFileParams, outFileIC, outFileConverge, geneTypes, signalVals, signalingTypes, stepper));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_sRACIPE_generateThresholds", (DL_FUNC) &_sRACIPE_generateThresholds, 3},
    {"_sRACIPE_limitcyclesGRC", (DL_FUNC) &_sRACIPE_limitcyclesGRC, 9},
    {"_sRACIPE_simulateGRCCpp", (DL_FUNC) &_sRACIPE_simulateGRCCpp, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_sRACIPE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
