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

using namespace Rcpp;

unsigned u_seed = std::chrono::system_clock::now().time_since_epoch().count();
unsigned g_seed = std::chrono::system_clock::now().time_since_epoch().count()*M_PI_4;


std::mt19937_64 u_generator (u_seed);
std::mt19937_64 g_generator (g_seed);

std::uniform_real_distribution<double> u_distribution(0.0,1.0);
std::normal_distribution<double> g_distribution(0.0,1.0);

double Hs_Racipe(double A, double AB0, int n_ab, double lambda_ab)
{
  return lambda_ab+(1-lambda_ab)*1/(1+pow((A/AB0),n_ab));
}

//////////////////////////////////////////////////////////////
// [[Rcpp::export]]
//////////////////////////////////////////////////////////////
int threshold_calculator_uniform(IntegerMatrix gene_interaction, NumericVector threshold_gene, double g_min, double g_max,
                                           double k_min, double k_max, int possible_interactions,  long model_count_max, long threshold_max, double h, double lambda_min,
                                           double lambda_max, int n_min, int n_max, double median_range, double standard_deviation_factor)
{
  int number_gene=gene_interaction.nrow();
  Rcout<<"generating thresholds for uniform distribution..."<<"\n";
  double gene_isolated_median=(g_min+g_max)/(k_min+k_max);
  for(int i=0;i<number_gene;i++){threshold_gene[i]= gene_isolated_median;}
  //for(int it1=0;it1<number_gene;it1++)for(int it2=0;it2<number_gene;it2++){gene_interaction2[it1][it2]=gene_interaction[it1][it2]; if(it1!=0) gene_interaction2[it1][it2]=0;}
  //for(int i=0;i<number_gene;i++){for(int j=0;j<number_gene;j++)cout<<gene_interaction2[i][j]<<"\t";cout<<endl;}
  //multiGeneCircuit_threshold(gene_interaction2, threshold_gene, rng1);
  //return;

  int interactions_number[number_gene][possible_interactions]; for(int i=0;i<number_gene;i++)for(int j=0;j<possible_interactions;j++)interactions_number[i][j]=0;
  for(int it1=0;it1<number_gene;it1++){for(int it2=0;it2<number_gene;it2++)
  {int g_1=gene_interaction(it1,it2);
    interactions_number[it1][g_1]+=1;}}
  //for(int i=0;i<number_gene;i++){for(int j=0;j<possible_interactions;j++)cout<<interactions_number[i][j]<<"\t";cout<<endl;}
  const long model_count_max2=std::max(model_count_max,threshold_max);
  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
  {
    std::vector<double> Af;
    for(long model_count=0;model_count<model_count_max2;model_count++)
    {
      double g_a=g_min+(g_max-g_min)*u_distribution(u_generator);
      double k_a=k_min+(k_max-k_min)*u_distribution(u_generator);
      Af.push_back(g_a/k_a);
      for(int counter_interaction1=0;counter_interaction1<interactions_number[gene_count1][1];counter_interaction1++)
      {
        double g_b=g_min+(g_max-g_min)*u_distribution(u_generator);
        double k_b=k_min+(k_max-k_min)*u_distribution(u_generator);
        double BA0=gene_isolated_median*(1-standard_deviation_factor*sqrt(3))+(2*sqrt(3)*standard_deviation_factor*gene_isolated_median)*u_distribution(u_generator);

        int n_ba=int(u_distribution(u_generator)*(n_max-n_min))+n_min;
        double gene_lambda=(lambda_max-lambda_min)*u_distribution(u_generator)+lambda_min;

        Af[model_count]=Af[model_count]*Hs_Racipe(g_b/k_b, BA0, n_ba, gene_lambda)/gene_lambda;
      }

      for(int counter_interaction1=0;counter_interaction1<interactions_number[gene_count1][2];counter_interaction1++)
      {
        double g_b=g_min+(g_max-g_min)*u_distribution(u_generator);
        double k_b=k_min+(k_max-k_min)*u_distribution(u_generator);
        //double BA0=exp(BA0_min+(BA0_max)*pcg32_double(rng1));
        double BA0=gene_isolated_median*(1-standard_deviation_factor*sqrt(3))+(2*sqrt(3)*standard_deviation_factor*gene_isolated_median)*u_distribution(u_generator);

        int n_ba=int(u_distribution(u_generator)*(n_max-n_min))+n_min;
        double gene_lambda=1./((lambda_max-lambda_min)*(n_max-n_min)+lambda_min);
        double hill_eval=Hs_Racipe(g_b/k_b, BA0, n_ba, gene_lambda);
        Af[model_count] = Af[model_count]*hill_eval;
      }
    }
    //Find the median of this array and select that as the median
    // size_t size = sizeof(Af) / sizeof(Af[0]);
    std::sort(Af.begin(), Af.end());
    //Rcout<<Af[int(float(model_count_max)*0.5)]<<"\n";
    threshold_gene[gene_count1]=0.5*(Af[int(float(model_count_max)*0.5)]+Af[int(float(model_count_max+2)*0.5)]);
    Af.clear();

  }
  return 0;
  //double A=exp(log(min_gene[gene_count1])+(log(max_gene[gene_count1])-log(min_gene[gene_count1]))*pcg32_double(rng1));
}
