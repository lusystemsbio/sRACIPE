#include"header.h"


using namespace Rcpp;


extern unsigned u_seed;// = std::chrono::system_clock::now().time_since_epoch().count();
//unsigned g_seed1 = std::chrono::system_clock::now().time_since_epoch().count()*M_PI_4;


extern std::mt19937_64 u_generator;
//std::mt19937_64 g_generator1 (g_seed1);

//uniformly distributed random number generator1 in (0,1) range
extern std::uniform_real_distribution<double> u_distribution;

// Gaussian distributed random number generator1 with mean 0 and 1 standard deviation
//std::normal_distribution<double> g_distribution1(0.0,1.0);


extern double Hs_Racipe(double A, double AB0, int n_ab, double lambda_ab);
// [[Rcpp::export]]

int threshold_calculator_uniform(IntegerMatrix gene_interaction, NumericVector threshold_gene, double g_min, double g_max,
                                           double k_min, double k_max, int possible_interactions,  long model_count_max, long threshold_max, double h, double lambda_min,
                                           double lambda_max, int n_min, int n_max, double median_range, double standard_deviation_factor)
{


  int number_gene=gene_interaction.nrow();
  Rcout<<"generating thresholds for uniform distribution1..."<<"\n";
  double gene_isolated_median=(g_min+g_max)/(k_min+k_max);
  for(int i=0;i<number_gene;i++){threshold_gene[i]= gene_isolated_median;}

  int interactions_number[number_gene][possible_interactions]; for(int i=0;i<number_gene;i++)for(int j=0;j<possible_interactions;j++)interactions_number[i][j]=0;
  for(int it1=0;it1<number_gene;it1++){for(int it2=0;it2<number_gene;it2++)
  {int g_1=gene_interaction(it1,it2);
    interactions_number[it1][g_1]+=1;}}

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
        double BA0=gene_isolated_median*(1-standard_deviation_factor*sqrt(3))+(2*sqrt(3)*standard_deviation_factor*gene_isolated_median)*u_distribution(u_generator);

        int n_ba=int(u_distribution(u_generator)*(n_max-n_min))+n_min;
        double gene_lambda=1./((lambda_max-lambda_min)*(n_max-n_min)+lambda_min);
        double hill_eval=Hs_Racipe(g_b/k_b, BA0, n_ba, gene_lambda);
        Af[model_count] = Af[model_count]*hill_eval;
      }
    }

    std::sort(Af.begin(), Af.end());
    threshold_gene[gene_count1]=0.5*(Af[int(float(model_count_max2)*0.5)]+Af[int(float(model_count_max2+2)*0.5)]);
    Af.clear();

  }
  return 0;

}
