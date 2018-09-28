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

// All definitions in threshold_generator.cpp

//uniformly distributed random number generator in (0,1) range
// Shifted hill function
extern double Hs_Racipe(double A, double AB0, int n_ab, double lambda_ab);
extern std::mt19937_64 u_generator;
extern std::uniform_real_distribution<double> u_distribution;
extern std::mt19937_64 g_generator;
extern std::normal_distribution<double> g_distribution;

#endif
