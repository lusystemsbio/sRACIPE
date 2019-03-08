#include"header.h"


using namespace Rcpp;


unsigned u_seed = 123;// = std::chrono::system_clock::now().time_since_epoch().count();
unsigned g_seed = 56454;//std::chrono::system_clock::now().time_since_epoch().count()*M_PI_4;

extern std::mt19937_64 u_generator (u_seed);
extern std::mt19937_64 g_generator (g_seed);


//uniformly distributed random number generator in (0,1) range
extern std::uniform_real_distribution<double> u_distribution(0.0,1.0);

// Gaussian distributed random number generator with mean 0 and 1 standard deviation
extern std::normal_distribution<double> g_distribution(0.0,1.0);

extern double Hs_Racipe(double A, double AB0, int n_ab, double lambda_ab)
{
  return lambda_ab+(1-lambda_ab)*1/(1+pow((A/AB0),n_ab));
}
