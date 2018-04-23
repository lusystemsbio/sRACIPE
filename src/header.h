
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


// Using system time for seeds of gaussian and uniform random number generators



//uniformly distributed random number generator in (0,1) range
// Shifted hill function
extern double Hs_Racipe(double A, double AB0, int n_ab, double lambda_ab); //defined in multiGeneCircuit_RK_detereministic
extern std::mt19937_64 u_generator; //defined in multiGeneCircuit_RK_deterministic
extern std::uniform_real_distribution<double> u_distribution; //defined in multiGeneCircuit_RK_detereministic
