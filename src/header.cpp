#include "header.h"


using namespace Rcpp;


unsigned u_seed = 123;// std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937_64 u_generator (u_seed);

//uniformly distributed random number generator in (0,1) range
std::uniform_real_distribution<double> u_distribution(0.0,1.0);

unsigned g_seed = 123456; //
std::mt19937_64 g_generator (g_seed);
std::normal_distribution<double> g_distribution(0.0,1.0);

extern double Hs_Racipe(double A, double AB0, int n_ab, double lambda_ab)
{
  return lambda_ab+(1-lambda_ab)*1/(1+std::pow((A/AB0),n_ab));
}

// function to convert gene_interaction matrix to vector of vectors
// source_gene vector will contain the sources of interactions
// interactionType will contain the type of interaction
size_t convertAdjMatToVector(
    IntegerMatrix gene_interaction, std::vector<size_t>& tgtGene,
    std::vector<std::pair<size_t,size_t> >& intSrcType){
  size_t nGene = gene_interaction.nrow();
  size_t nInteractions = 0;
  //  Rcout<<nGene<<"\n";
  for(size_t tmpTgtGene = 0; tmpTgtGene < nGene; tmpTgtGene++ ){
    for(size_t tmpSrcGene = 0; tmpSrcGene < nGene; tmpSrcGene++ ){
      if(gene_interaction(tmpTgtGene,tmpSrcGene) > 0){
        nInteractions++;
        tgtGene.push_back(tmpTgtGene);
        intSrcType.push_back(
          std::make_pair(tmpSrcGene,(gene_interaction(tmpTgtGene,tmpSrcGene))));

      }
    }
  }
  return nInteractions;
}
