#include "header.h"
#include <Rcpp.h>
using namespace Rcpp;



void stepEMconv( std::vector <double> &exprxGene,
             std::ofstream &outGE,
             const int &numberGene,
             IntegerMatrix geneInteraction,
             const std::vector<double> &gGene,
             const std::vector<double> &kGene,
             const std::vector<std::vector<int> > &NGene,
             const std::vector<std::vector<double> > &lambda_gene,
             const std::vector<std::vector<double> > &threshold_gene_log,
             const int &possible_interactions,
             const double &standard_deviation_factor,
             const double &D_shot_scaling,
             const std::vector<double> &Darray,
             const int &outputPrecision,
             const double &D,
             const double &h,
             const double &signalRate,
             const NumericVector &geneTypes,
             const long double &convergThresh,
             const int &numStepsConverge,
             const int &numConvergenceTests){

  double exprxGeneH[numberGene]; //array for temp gene expression values
  std::vector<double> prevExprxGene(numberGene);
  for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
  {
    exprxGeneH[geneCount1] = exprxGene[geneCount1];
  }

  double i=0.0;
  
  //For each test, we run the simulation for numStepsConverge iterations
  //and check if the system has changed state in that time
  for(int testIter=0; testIter<numConvergenceTests; testIter++){

    for(int geneCount1=0;geneCount1<numberGene;geneCount1++){
        prevExprxGene[geneCount1]=exprxGene[geneCount1];}
    
    for(int j=0;j<numStepsConverge;j++){
        i+=h;
        for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
        {
          double growthMultiplier=1;
          double degMultiplier=1;

          for(int geneCount2=0;geneCount2<numberGene;geneCount2++)
          {
          double geneValue=exprxGene[geneCount2];
          double geneThreshold=threshold_gene_log[geneCount1][geneCount2];
          int geneN=NGene[geneCount1][geneCount2];
          double geneLambda=lambda_gene[geneCount1][geneCount2];
          calcMultiplier(geneCount1, geneCount2, growthMultiplier, degMultiplier,
                            geneValue, geneInteraction, geneN, geneLambda,
                            geneThreshold);
          }
        if (geneTypes[geneCount1] == 2){
          growthMultiplier = growthMultiplier*signalRate;
          degMultiplier = degMultiplier*signalRate;
        }
        exprxGeneH[geneCount1] = exprxGene[geneCount1] +
          h*(gGene[geneCount1]*growthMultiplier-kGene[geneCount1]*
          exprxGene[geneCount1]*degMultiplier) +
          D*sqrt(h)*g_distribution(g_generator)*Darray[geneCount1]+
          D_shot_scaling*D*sqrt(h)*g_distribution(g_generator)*
          Darray[geneCount1]*exprxGene[geneCount1];
        if(exprxGeneH[geneCount1]<0) exprxGeneH[geneCount1]=0;
      }
        for(int geneCount1=0;geneCount1<numberGene;geneCount1++){
        exprxGene[geneCount1]=exprxGeneH[geneCount1];}
    }

    if(abs(D) < 1e-5){
      double test_delta;
      test_delta = sum_delta(prevExprxGene, exprxGene, numberGene);
      if (test_delta < convergThresh){
        for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
        {
          outGE<<std::setprecision(outputPrecision)
          <<exprxGene[geneCount1]<<"\t";
        }
        break;
      }
    }


  };

  // for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
  // {
  //   outGE<<std::setprecision(outputPrecision)<<exprxGene[geneCount1]<<"\t";
  // }

}

///////////////////////////////////////////////////////////////////////////////

//Runge Kutta fourth order
///////////////////////////////////////////////////////////////////////////////
void stepRK4conv( std::vector <double> &exprxGene,
             std::ofstream &outGE,
             const int &numberGene,
             IntegerMatrix geneInteraction,
             const std::vector<double> &gGene,
             const std::vector<double> &kGene,
             const std::vector<std::vector<int> > &NGene,
             const std::vector<std::vector<double> > &lambda_gene,
             const std::vector<std::vector<double> > &threshold_gene_log,
             const int &possible_interactions,
             const double &standard_deviation_factor,
             const int &outputPrecision,
             const double &h,
             const double &signalRate,
             const NumericVector &geneTypes,
             const long double &convergThresh,
             const int &numStepsConverge,
             const int &numConvergenceTests){

double exprxGeneH1[numberGene]; //array for temp gene expression values
double exprxGeneH2[numberGene]; //array for temp gene expression values
double exprxGeneH3[numberGene]; //array for temp gene expression values
double exprxGeneH4[numberGene]; //array for temp gene expression values

std::vector<double> prevExprxGene(numberGene); //Array for convergence testing

for(int geneCountTmp=0;geneCountTmp<numberGene;geneCountTmp++)
{
  exprxGeneH1[geneCountTmp]=exprxGene[geneCountTmp];
  exprxGeneH2[geneCountTmp]=exprxGene[geneCountTmp];
  exprxGeneH3[geneCountTmp]=exprxGene[geneCountTmp];
  exprxGeneH4[geneCountTmp]=exprxGene[geneCountTmp];
}
double i=0.0;

//For each test, we run the simulation for numStepsConverge iterations
//and check if the system has changed state in that time
for(int testIter=0; testIter<numConvergenceTests; testIter++){

  for(int geneCount1=0;geneCount1<numberGene;geneCount1++){
        prevExprxGene[geneCount1]=exprxGene[geneCount1];}
      
  for(int j=0;j<numStepsConverge;j++){
    i+=h;
      
    for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
    {
      double growthMultiplier=1;
      double degMultiplier=1;

      for(int geneCount2=0;geneCount2<numberGene;geneCount2++)
      {
        double geneValue=exprxGene[geneCount2];
        double geneThreshold=threshold_gene_log[geneCount1][geneCount2];
        int geneN=NGene[geneCount1][geneCount2];
        double geneLambda=lambda_gene[geneCount1][geneCount2];
        calcMultiplier(geneCount1, geneCount2, growthMultiplier, degMultiplier,
                       geneValue, geneInteraction, geneN, geneLambda,
                       geneThreshold);
      }

      if (geneTypes[geneCount1] == 2){
        growthMultiplier = growthMultiplier*signalRate;
        degMultiplier = degMultiplier*signalRate;
      }

      exprxGeneH1[geneCount1]=h*(gGene[geneCount1]*growthMultiplier -
        kGene[geneCount1]*exprxGene[geneCount1]*degMultiplier);
    }

    for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
    {
      double growthMultiplier=1;
      double degMultiplier=1;

      for(int geneCount2=0;geneCount2<numberGene;geneCount2++)
      {
        double geneValue=exprxGene[geneCount2] +
          0.5*exprxGeneH1[geneCount2];
        double geneThreshold=threshold_gene_log[geneCount1][geneCount2];
        int geneN=NGene[geneCount1][geneCount2];
        double geneLambda=lambda_gene[geneCount1][geneCount2];
        calcMultiplier(geneCount1, geneCount2, growthMultiplier, degMultiplier,
                       geneValue, geneInteraction, geneN, geneLambda,
                       geneThreshold);
      }

      if (geneTypes[geneCount1] == 2){
        growthMultiplier = growthMultiplier*signalRate;
        degMultiplier = degMultiplier*signalRate;
      }

      exprxGeneH2[geneCount1]=h*((gGene[geneCount1])*
        growthMultiplier-kGene[geneCount1]*(exprxGene[geneCount1] +
        0.5*exprxGeneH1[geneCount1])*degMultiplier);
    }

    for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
    {
      double growthMultiplier=1;
      double degMultiplier=1;

      for(int geneCount2=0;geneCount2<numberGene;geneCount2++)
      {
        double geneValue=exprxGene[geneCount2] +
          0.5*exprxGeneH2[geneCount2];
        double geneThreshold=threshold_gene_log[geneCount1][geneCount2];
        int geneN=NGene[geneCount1][geneCount2];
        double geneLambda=lambda_gene[geneCount1][geneCount2];
        calcMultiplier(geneCount1, geneCount2, growthMultiplier, degMultiplier,
                       geneValue, geneInteraction, geneN, geneLambda,
                       geneThreshold);
      }

      if (geneTypes[geneCount1] == 2){
        growthMultiplier = growthMultiplier*signalRate;
        degMultiplier = degMultiplier*signalRate;
      }

      exprxGeneH3[geneCount1]=h*((gGene[geneCount1])*growthMultiplier-
        kGene[geneCount1]*(exprxGene[geneCount1] +
        0.5*exprxGeneH2[geneCount1])*degMultiplier);
    }


    for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
    {
      double growthMultiplier=1;
      double degMultiplier=1;

      for(int geneCount2=0;geneCount2<numberGene;geneCount2++)
      {
        double geneValue=exprxGene[geneCount2] +
          exprxGeneH3[geneCount2];
        double geneThreshold=threshold_gene_log[geneCount1][geneCount2];
        int geneN=NGene[geneCount1][geneCount2];
        double geneLambda=lambda_gene[geneCount1][geneCount2];
        calcMultiplier(geneCount1, geneCount2, growthMultiplier, degMultiplier,
                      geneValue, geneInteraction, geneN, geneLambda,
                      geneThreshold);
      }

      if (geneTypes[geneCount1] == 2){
        growthMultiplier = growthMultiplier*signalRate;
        degMultiplier = degMultiplier*signalRate;
      }

      exprxGeneH4[geneCount1]=h*((gGene[geneCount1])*growthMultiplier-
        kGene[geneCount1]*(exprxGene[geneCount1]+
        exprxGeneH3[geneCount1])*degMultiplier);
    }
      /////////////////////////////////////////////////////////////
    for(int geneCount1=0;geneCount1<numberGene;geneCount1++){
      exprxGene[geneCount1] = exprxGene[geneCount1]+
        (exprxGeneH1[geneCount1]+2*exprxGeneH2[geneCount1]+
        2*exprxGeneH3[geneCount1]+exprxGeneH4[geneCount1])/6;
      if(exprxGene[geneCount1]<0) exprxGene[geneCount1]=0;
    }
  }

  
    double test_delta;
    test_delta = sum_delta(prevExprxGene, exprxGene, numberGene);
    if (test_delta < convergThresh){
      for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
        {
        outGE<<std::setprecision(outputPrecision)
        <<exprxGene[geneCount1]<<"\t";
        }
      break;
      }

};
}


