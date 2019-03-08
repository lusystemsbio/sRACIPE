#include"header.h"

using namespace Rcpp;


// [[Rcpp::export]]

int generateThresholds(
    const Rcpp::IntegerMatrix gene_interaction,
    Rcpp::NumericVector threshold_gene, const double g_min,
    const double g_max, const double k_min, const double k_max,
    const int possible_interactions,
    const long model_count_max, const long threshold_max, const double h,
    const double lambda_min,
    const double lambda_max, const int n_min, const int n_max,
    const double standard_deviation_factor)
{

// Generate Thresholds
  int number_gene=gene_interaction.nrow();
  Rcout<<"generating thresholds for uniform distribution1..."<<"\n";
  double gene_isolated_median=(g_min+g_max)/(k_min+k_max);
  for(int i=0;i<number_gene;i++){threshold_gene[i]= gene_isolated_median;}

  int interactions_number[number_gene][possible_interactions];
  for(int i=0;i<number_gene;i++)for(int j=0;j<possible_interactions;j++)
    interactions_number[i][j]=0;

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
      for(int counter_interaction1=0;
          counter_interaction1<interactions_number[gene_count1][1];
          counter_interaction1++)
      {
        double g_b=g_min+(g_max-g_min)*u_distribution(u_generator);
        double k_b=k_min+(k_max-k_min)*u_distribution(u_generator);
        double BA0=gene_isolated_median*(1-standard_deviation_factor*sqrt(3)) +
          (2*sqrt(3)*standard_deviation_factor*gene_isolated_median)
          *u_distribution(u_generator);

        int n_ba=int(u_distribution(u_generator)*(n_max-n_min))+n_min;
        double gene_lambda=(lambda_max-lambda_min)*u_distribution(u_generator) +
          lambda_min;

        Af[model_count]=Af[model_count]*Hs_Racipe(g_b/k_b, BA0, n_ba,
                                                  gene_lambda)/gene_lambda;
      }

      for(int counter_interaction1=0;
          counter_interaction1<interactions_number[gene_count1][2];
          counter_interaction1++)
      {
        double g_b=g_min+(g_max-g_min)*u_distribution(u_generator);
        double k_b=k_min+(k_max-k_min)*u_distribution(u_generator);
        double BA0=gene_isolated_median*(1-standard_deviation_factor*sqrt(3)) +
          (2*sqrt(3)*standard_deviation_factor*gene_isolated_median)
          *u_distribution(u_generator);

        int n_ba=int(u_distribution(u_generator)*(n_max-n_min))+n_min;
        double gene_lambda=1./((lambda_max-lambda_min)*(n_max-n_min)+lambda_min);
        double hill_eval=Hs_Racipe(g_b/k_b, BA0, n_ba, gene_lambda);
        Af[model_count] = Af[model_count]*hill_eval;
      }
    }

    std::sort(Af.begin(), Af.end());
    // take the median as threshold
    threshold_gene[gene_count1]=0.5*(Af[int(float(model_count_max2)*0.5)] +
      Af[int(float(model_count_max2+2)*0.5)]);

    Af.clear();
  }
  return 0;

}
