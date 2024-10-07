#ifndef HEADER_H
#define HEADER_H

#include <Rcpp.h>
#include <iostream>
#include <random>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <sstream>
#include <string>

//uniformly distributed random number generator in (0,1) range
// Shifted hill function
extern double Hs_Racipe(double A, double AB0, int n_ab, double lambda_ab);
extern std::mt19937_64 u_generator;
extern std::uniform_real_distribution<double> u_distribution;
extern std::mt19937_64 g_generator;
// Gaussian distributed random number generator with mean 0 and 1
//standard deviation
extern std::normal_distribution<double> g_distribution;

// Rcpp::IntegerMatrix readTopology(
//     Rcpp::IntegerMatrix gene_interaction, const Rcpp::String filepath,
//     const Rcpp::String filename, Rcpp::StringVector geneNames);

int generateThresholds(
    const Rcpp::IntegerMatrix gene_interaction,
    Rcpp::NumericVector threshold_gene, const double g_min, const double g_max,
    const double k_min, const double k_max, const int possible_interactions,
    const long model_count_max, const long threshold_max, const double h,
    const double lambda_min,
    const double lambda_max, const int n_min, const int n_max,
    const double standard_deviation_factor);

void stepEM( std::vector <double> &expression_gene,
             std::ofstream &out_GE,
             const double &tot_time,
             const int &number_gene,
             const Rcpp::IntegerMatrix gene_interaction,
             const std::vector<double> &g_gene,
             const std::vector<double> &k_gene,
             const std::vector<std::vector<int> > &n_gene,
             const std::vector<std::vector<double> > &lambda_gene,
             const std::vector<std::vector<double> > &threshold_gene_log,
             const int &possible_interactions,
             const double &standard_deviation_factor,
             const double &D_shot_scaling,
             const std::vector<double> &Darray,
             const int &output_precision,
             const double &print_start, const double &print_interval,
             const double &D,
             const double &h,
             const double &signalRate,
             const Rcpp::NumericVector &geneTypes,
             const bool &isTimeVarying,
             const std::vector<double> &timePoints,
             const Rcpp::NumericMatrix &signalVals,
             const Rcpp::NumericVector &signalingTypes,
             std::unordered_map<int, std::vector<double>> &clamps,
             const int &modelNo);
void stepEM_OU( std::vector <double> &expression_gene,
            std::ofstream &out_GE,
            const double &tot_time,
            const int &number_gene,
            const Rcpp::IntegerMatrix gene_interaction,
            const std::vector<double> &g_gene,
            const std::vector<double> &k_gene,
            const std::vector<std::vector<int> > &n_gene,
            const std::vector<std::vector<double> > &lambda_gene,
            const std::vector<std::vector<double> > &threshold_gene_log,
            const int &possible_interactions,
            const double &standard_deviation_factor,
            const double &D_shot_scaling,
            const std::vector<double> &Darray,
            const int &output_precision,
            const double &print_start, 
            const double &print_interval,
            const double &D,
            const double &h,
            const double &ouNoise_tcorr,
            const double &signalRate,
            const Rcpp::NumericVector &geneTypes,
            const bool &isTimeVarying,
            const std::vector<double> &timePoints,
            const Rcpp::NumericMatrix &signalVals,
            const Rcpp::NumericVector &signalingTypes,
            std::unordered_map<int, std::vector<double>> &clamps,
            const int &modelNo);
void stepRK4( std::vector <double> &expression_gene,
        std::ofstream &out_GE,
        const double &tot_time,
        const int &number_gene,
        const Rcpp::IntegerMatrix gene_interaction,
        const std::vector<double> &g_gene,
        const std::vector<double> &k_gene,
        const std::vector<std::vector<int> > &n_gene,
        const std::vector<std::vector<double> > &lambda_gene,
        const std::vector<std::vector<double> > &threshold_gene_log,
        const int &possible_interactions,
        const double &standard_deviation_factor,
        const int &output_precision,
        const double &print_start, const double &print_interval,
        const double &h,
        const double &signalRate,
        const Rcpp::NumericVector &geneTypes,
        const bool &isTimeVarying,
        const std::vector<double> &timePoints,
        const Rcpp::NumericMatrix &signalVals,
        const Rcpp::NumericVector &signalingTypes,
        std::unordered_map<int, std::vector<double>> &clamps,
        const int &modelNo);

 void stepDP( std::vector <double> &expression_gene,
          std::ofstream &out_GE,
          const double &tot_time,
          const int &number_gene,
          const Rcpp::IntegerMatrix gene_interaction,
          const std::vector<double> &g_gene,
          const std::vector<double> &k_gene,
          const std::vector<std::vector<int> > &n_gene,
          const std::vector<std::vector<double> > &lambda_gene,
          const std::vector<std::vector<double> > &threshold_gene_log,
          const int &possible_interactions,
          const double &standard_deviation_factor,
          const int &output_precision,
          const double &print_start, const double &print_interval,
          double h, const double &rk_tolerance,
          const double &signalRate,
          const Rcpp::NumericVector &geneTypes,
          const bool &isTimeVarying,
          const std::vector<double> &timePoints,
          const Rcpp::NumericMatrix &signalVals,
          const Rcpp::NumericVector &signalingTypes,
          std::unordered_map<int, std::vector<double>> &clamps,
          const int &modelNo);

void stepEMconv( std::vector <double> &expression_gene,
             std::ofstream &out_GE,
             std::ofstream &out_Conv,
             const int &number_gene,
             const Rcpp::IntegerMatrix gene_interaction,
             const std::vector<double> &g_gene,
             const std::vector<double> &k_gene,
             const std::vector<std::vector<int> > &n_gene,
             const std::vector<std::vector<double> > &lambda_gene,
             const std::vector<std::vector<double> > &threshold_gene_log,
             const int &possible_interactions,
             const double &standard_deviation_factor,
             const double &D_shot_scaling,
             const std::vector<double> &Darray,
             const int &output_precision,
             const double &D,
             const double &h,
             const double &signalRate,
             const Rcpp::NumericVector &geneTypes,
             const long double &convergThresh,
             const int &numStepsConverge,
             const int &numConvergenceTests,
             std::unordered_map<int, std::vector<double>> &clamps,
             const int &modelNo);
void stepRK4conv( std::vector <double> &expression_gene,
        std::ofstream &out_GE,
        std::ofstream &out_Conv,
        const int &number_gene,
        const Rcpp::IntegerMatrix gene_interaction,
        const std::vector<double> &g_gene,
        const std::vector<double> &k_gene,
        const std::vector<std::vector<int> > &n_gene,
        const std::vector<std::vector<double> > &lambda_gene,
        const std::vector<std::vector<double> > &threshold_gene_log,
        const int &possible_interactions,
        const double &standard_deviation_factor,
        const int &output_precision,
        const double &h,
        const double &signalRate,
        const Rcpp::NumericVector &geneTypes,
        const long double &convergThresh,
        const int &numStepsConverge,
        const int &numConvergenceTests,
        std::unordered_map<int, std::vector<double>> &clamps,
        const int &modelNo);

 void stepDPconv( std::vector <double> &expression_gene,
          std::ofstream &out_GE,
          std::ofstream &out_Conv,
          const int &number_gene,
          const Rcpp::IntegerMatrix gene_interaction,
          const std::vector<double> &g_gene,
          const std::vector<double> &k_gene,
          const std::vector<std::vector<int> > &n_gene,
          const std::vector<std::vector<double> > &lambda_gene,
          const std::vector<std::vector<double> > &threshold_gene_log,
          const int &possible_interactions,
          const double &standard_deviation_factor,
          const int &output_precision,
          double h, const double &rk_tolerance,
          const double &signalRate,
          const Rcpp::NumericVector &geneTypes,
          const long double &convergThresh,
          const int &numStepsConverge,
          const int &numConvergenceTests,
          const double &testTime,
          std::unordered_map<int, std::vector<double>> &clamps,
          const int &modelNo);


extern size_t convertAdjMatToVector(
     Rcpp::IntegerMatrix gene_interaction, std::vector<size_t>& tgtGene,
     std::vector<std::pair<size_t,size_t> >& intSrcType);

double sum_delta (std::vector<double> &exprxGene, 
                  std::vector<double> &exprxGeneH, int numberGene);

void calcMultiplier(const int& geneCount1, const int& geneCount2,
                    double& growthMultiplier,
                    double& degMultiplier,
                    double& geneValue,
                    const Rcpp::IntegerMatrix geneInteraction,
                    const int& geneN,
                    double& geneLambda,
                    const double& geneThreshold);

void readParameters(Rcpp::IntegerMatrix geneInteraction, 
                        const int &numberGene,
                        std::vector<double> &gGene,
                        std::vector<double> &kGene,
                        std::vector<std::vector<int> > &nGene,
                        std::vector<std::vector<double> > &lambdaGene,
                        std::vector<std::vector<double> > &threshGeneLog,
                        std::ifstream &inParams);

#endif
