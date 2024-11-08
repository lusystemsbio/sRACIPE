#include "header.h"
#include <Rcpp.h>
using namespace Rcpp;



void stepEMconv( std::vector <double> &exprxGene,
             std::ofstream &outGE,
             std::ofstream &outConv,
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
             const int &numConvergenceIter,
             const bool &noClamps,
             std::unordered_map<int, std::vector<double>> &clamps,
             const int &modelNo){

  double exprxGeneH[numberGene]; //array for temp gene expression values
  std::vector<double> prevExprxGene(numberGene);
  for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
  {
    exprxGeneH[geneCount1] = exprxGene[geneCount1];
  }

  double i=0.0;
  bool isConverged = false;
  bool isNegative = false; //Checking for negative values due to stiffness
  
  //For each test, we run the simulation for numStepsConverge iterations
  //and check if the system has changed state in that time
  for(int testIter=0; testIter<numConvergenceIter; testIter++){

    for(int geneCount1=0;geneCount1<numberGene;geneCount1++){
        prevExprxGene[geneCount1]=exprxGene[geneCount1];}
    
    for(int j=0;j<numStepsConverge;j++){
        i+=h;

        for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
        {
          double growthMultiplier=1;
          double degMultiplier=1;

          if(!noClamps){
            auto it = clamps.find(geneCount1);
            if (it != clamps.end()){
              // If clamped, set to the clamped value
              exprxGeneH[geneCount1] = it->second[modelNo];
            }
          }else{
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
              exprxGene[geneCount1]*degMultiplier);
            if(exprxGeneH[geneCount1]<0) exprxGeneH[geneCount1]=0;
          }
      }
        for(int geneCount1=0;geneCount1<numberGene;geneCount1++){
        exprxGene[geneCount1]=exprxGeneH[geneCount1];}
    }

      double test_delta;
      test_delta = sum_delta(prevExprxGene, exprxGene, numberGene);
      if (test_delta < convergThresh){
        for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
        {
          outGE<<std::setprecision(outputPrecision)
          <<exprxGene[geneCount1]<<"\t";
        }
        isConverged = true;
        outConv<<isConverged<<"\t" << testIter <<"\n";
        break;
      }
      else if(isNegative){
        for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
          {
          outGE<<std::setprecision(outputPrecision)
          <<exprxGene[geneCount1]<<"\t";
          }
        isConverged = true;  //Necessary to stop double reporting
        outConv<<3<<"\t" << testIter <<"\n";
        break;
      }

  };

  if(isConverged == false){
    for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
    {
      outGE<<std::setprecision(outputPrecision)
      <<exprxGene[geneCount1]<<"\t";
    }
    outConv<<isConverged<<"\t" << numConvergenceIter <<"\n";
  }


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
             std::ofstream &outConv,
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
             const int &numConvergenceIter,
             const bool &noClamps,
             std::unordered_map<int, std::vector<double>> &clamps,
             const int &modelNo){

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
bool isConverged = false;
bool isBlowup = false; //Checking for explosions due to stiffness
bool isNegative = false; //Checking for negative values due to stiffness

//For each test, we run the simulation for numStepsConverge iterations
//and check if the system has changed state in that time
for(int testIter=0; testIter<numConvergenceIter; testIter++){

  for(int geneCount1=0;geneCount1<numberGene;geneCount1++){
        prevExprxGene[geneCount1]=exprxGene[geneCount1];}
      
  for(int j=0;j<numStepsConverge;j++){
    i+=h;
      
    for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
    {
      double growthMultiplier=1;
      double degMultiplier=1;

      if(!noClamps){
        auto it = clamps.find(geneCount1);
        if (it != clamps.end()){
          // If clamped, set derivative to 0
          exprxGeneH1[geneCount1] = 0;
        }
      }else{
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
    }

    for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
    {
      double growthMultiplier=1;
      double degMultiplier=1;

      if(!noClamps){
        auto it = clamps.find(geneCount1);
        if (it != clamps.end()){
          // If clamped, derivative to 0
          exprxGeneH2[geneCount1] = 0.0;
      }
      }else{
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
    }

    for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
    {
      double growthMultiplier=1;
      double degMultiplier=1;

      if(!noClamps){
        auto it = clamps.find(geneCount1);
        if (it != clamps.end()){
          // If clamped, set derivative to 0
          exprxGeneH3[geneCount1] = 0.0;
        }
      }else{
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
    }


    for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
    {
      double growthMultiplier=1;
      double degMultiplier=1;

      if(!noClamps){
        auto it = clamps.find(geneCount1);
        if (it != clamps.end()){
          // If clamped, set derivative to 0
          exprxGeneH4[geneCount1] = 0.0;
        }
      }else{
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
    }
      /////////////////////////////////////////////////////////////
    for(int geneCount1=0;geneCount1<numberGene;geneCount1++){
      exprxGene[geneCount1] = exprxGene[geneCount1]+
        (exprxGeneH1[geneCount1]+2*exprxGeneH2[geneCount1]+
        2*exprxGeneH3[geneCount1]+exprxGeneH4[geneCount1])/6;
      if(exprxGene[geneCount1]<0) isNegative = true;
      if(std::isinf(exprxGene[geneCount1])) isBlowup = true;
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
      isConverged = true;  
      outConv<<isConverged<<"\t" << testIter <<"\n";
      break;
      }
    else if(isBlowup || isNegative){
      for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
        {
        outGE<<std::setprecision(outputPrecision)
        <<exprxGene[geneCount1]<<"\t";
        }
      isConverged = true;  //Necessary to stop double reporting
      outConv<<3<<"\t" << testIter <<"\n";
      break;
    }

};
if(isConverged == false){
  for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
  {
    outGE<<std::setprecision(outputPrecision)
    <<exprxGene[geneCount1]<<"\t";
  }
  outConv<<isConverged<<"\t" << numConvergenceIter <<"\n";
}


}

// Dormand-Prince 5th order adaptive step integrator with
// convergence tests

void stepDPconv( std::vector <double> &exprxGene,
              std::ofstream &outGE,
              std::ofstream &outConv,
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
              double h, const double &rkTolerance,
              const double &signalRate,
              const NumericVector &geneTypes,
              const long double &convergThresh,
              const int &numStepsConverge,
              const int &numConvergenceIter,
              const double &testTime,
              const bool &noClamps,
              std::unordered_map<int, std::vector<double>> &clamps,
              const int &modelNo){
  double exprxGeneH[numberGene]; //array for temp gene expression values
  double exprxGeneH1[numberGene]; //array for temp gene expression values
  double exprxGeneH2[numberGene]; //array for temp gene expression values
  double exprxGeneH3[numberGene]; //array for temp gene expression values
  double exprxGeneH4[numberGene]; //array for temp gene expression values
  double exprxGeneH5[numberGene]; //array for temp gene expression values
  double exprxGeneH6[numberGene]; //array for temp gene expression values
  double exprxGeneH7[numberGene]; //array for temp gene expression values

  bool isConverged = false;
  bool isNegative = false;

  std::vector<double> prevExprxGene(numberGene); //Array for convergence testing

  for(int geneCountTmp=0;geneCountTmp<numberGene;geneCountTmp++)
  {
    exprxGeneH[geneCountTmp]=exprxGene[geneCountTmp];
    exprxGeneH1[geneCountTmp]=exprxGene[geneCountTmp];
    exprxGeneH2[geneCountTmp]=exprxGene[geneCountTmp];
    exprxGeneH3[geneCountTmp]=exprxGene[geneCountTmp];
    exprxGeneH4[geneCountTmp]=exprxGene[geneCountTmp];
    exprxGeneH5[geneCountTmp]=exprxGene[geneCountTmp];
    exprxGeneH6[geneCountTmp]=exprxGene[geneCountTmp];
    exprxGeneH7[geneCountTmp]=exprxGene[geneCountTmp];
  }

  double i=0.0;

  //For each test, we run the simulation for numStepsConverge iterations
  //and check if the system has changed state in that time
  for(int testIter=0; testIter<numConvergenceIter; testIter++){
    for(int geneCount1=0;geneCount1<numberGene;geneCount1++){
        prevExprxGene[geneCount1]=exprxGene[geneCount1];}
      
    for(int j=0;j<numStepsConverge;j++){
      for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
      {
        double growthMultiplier=1;
        double degMultiplier=1;

        if(!noClamps){
          auto it = clamps.find(geneCount1);
          if (it != clamps.end()){
            // If clamped, set derivative to 0
            exprxGeneH1[geneCount1] = 0.0;
          }
        }else{
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

          exprxGeneH1[geneCount1]=h*(gGene[geneCount1]*growthMultiplier-
            kGene[geneCount1]*exprxGene[geneCount1]*degMultiplier);
        }
      }

      for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
      {
        double growthMultiplier=1;
        double degMultiplier=1;

        if(!noClamps){
          auto it = clamps.find(geneCount1);
          if (it != clamps.end()){
            // If clamped, set derivative to 0
            exprxGeneH2[geneCount1] = 0.0;
          }
        }else{
          for(int geneCount2=0;geneCount2<numberGene;geneCount2++)
          {
            double geneValue=exprxGene[geneCount2]+
              0.2*exprxGeneH1[geneCount2];
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

          exprxGeneH2[geneCount1]=h*((gGene[geneCount1])*growthMultiplier-
            kGene[geneCount1]*(exprxGene[geneCount1] +
            0.2*exprxGeneH1[geneCount1])*degMultiplier);
        }
      }

      for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
      {
        double growthMultiplier=1;
        double degMultiplier=1;

        if(!noClamps){
          auto it = clamps.find(geneCount1);
          if (it != clamps.end()){
            // If clamped, set derivative to 0
            exprxGeneH3[geneCount1] = 0.0;
          }
        }else{
          for(int geneCount2=0;geneCount2<numberGene;geneCount2++)
          {
            double geneValue=exprxGene[geneCount2]+
              (0.25*exprxGeneH1[geneCount2]+
              0.75*exprxGeneH2[geneCount2])*0.3;
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
            kGene[geneCount1]*(exprxGene[geneCount1]+
            (0.25*exprxGeneH1[geneCount1]+
            0.75*exprxGeneH2[geneCount1])*0.3)*degMultiplier);
        }
      }


      for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
      {
        double growthMultiplier=1;
        double degMultiplier=1;

        if(!noClamps){
          auto it = clamps.find(geneCount1);
          if (it != clamps.end()){
            // If clamped, set derivative to 0
            exprxGeneH4[geneCount1] = 0.0;
          }
        }else{
          for(int geneCount2=0;geneCount2<numberGene;geneCount2++)
          {
            double geneValue=exprxGene[geneCount2]+
              0.8*(((11./9.)*exprxGeneH1[geneCount2]+
              (-14./3.)*exprxGeneH2[geneCount2])+
              (40./9.)*exprxGeneH3[geneCount2]);
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
            0.8*(((11./9.)*exprxGeneH1[geneCount1]+
            (-14./3.)*exprxGeneH2[geneCount1])+
            (40./9.)*exprxGeneH3[geneCount1]))*degMultiplier);
        }
      }

      for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
      {
        double growthMultiplier=1;
        double degMultiplier=1;

        if(!noClamps){
          auto it = clamps.find(geneCount1);
          if (it != clamps.end()){
            // If clamped, set derivative to 0
            exprxGeneH5[geneCount1] = 0.0;
          }
        }else{
          for(int geneCount2=0;geneCount2<numberGene;geneCount2++)
          {
            double geneValue=exprxGene[geneCount2]+(8./9.)*
              ((4843./1458.)*exprxGeneH1[geneCount2]+
              (-3170./243.)*exprxGeneH2[geneCount2]+
              (8056./729.)*exprxGeneH3[geneCount2]+
              (-53./162.)*exprxGeneH4[geneCount2]);

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

          exprxGeneH5[geneCount1]=h*((gGene[geneCount1])*growthMultiplier-
            kGene[geneCount1]*(exprxGene[geneCount1]+(8./9.)*
            ((4843./1458.)*exprxGeneH1[geneCount1]+
            (-3170./243.)*exprxGeneH2[geneCount1]+
            (8056./729.)*exprxGeneH3[geneCount1]+
            (-53./162.)*exprxGeneH4[geneCount1]))*degMultiplier);
        }
      }


      for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
      {
        double growthMultiplier=1;
        double degMultiplier=1;

        if(!noClamps){
          auto it = clamps.find(geneCount1);
          if (it != clamps.end()){
            // If clamped, set derivative to 0
            exprxGeneH6[geneCount1] = 0.0;
          }
        }else{
          for(int geneCount2=0;geneCount2<numberGene;geneCount2++)
          {
            double geneValue=exprxGene[geneCount2]+
              ((9017./3168.)*exprxGeneH1[geneCount2]+
              (-355./33.)*exprxGeneH2[geneCount2]+
              (46732./5247.)*exprxGeneH3[geneCount2]+
              (49./176.)*exprxGeneH4[geneCount2]+
              (-5103./18656.)*exprxGeneH5[geneCount2]);
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

          exprxGeneH6[geneCount1]=h*((gGene[geneCount1])*growthMultiplier-
            kGene[geneCount1]*(exprxGene[geneCount1]+
            (9017./3168.)*exprxGeneH1[geneCount1]+
            (-355./33.)*exprxGeneH2[geneCount1]+
            (46732./5247.)*exprxGeneH3[geneCount1]+
            (49./176.)*exprxGeneH4[geneCount1]+
           (-5103./18656.)*exprxGeneH5[geneCount1])*degMultiplier);
        }
      }



      for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
      {
        double growthMultiplier=1;
        double degMultiplier=1;

        if(!noClamps){
          auto it = clamps.find(geneCount1);
          if (it != clamps.end()){
            // If clamped, set derivative to 0
            exprxGeneH7[geneCount1] = 0.0;
          }
        }else{
          for(int geneCount2=0;geneCount2<numberGene;geneCount2++)
          {
            double geneValue=exprxGene[geneCount2]+
              ((35./384.)*exprxGeneH1[geneCount2]+
              (500./113.)*exprxGeneH3[geneCount2]+
              (125./192.)*exprxGeneH4[geneCount2]+
              (-2187./6784.)*exprxGeneH5[geneCount2]+
                (11./84.)*exprxGeneH6[geneCount2]);
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

          exprxGeneH7[geneCount1]=h*((gGene[geneCount1])*growthMultiplier-
            kGene[geneCount1]*(exprxGene[geneCount1]+
            (35./384.)*exprxGeneH1[geneCount1]+
            (500./113.)*exprxGeneH3[geneCount1]+
            (125./192.)*exprxGeneH4[geneCount1]+
              (-2187./6784.)*exprxGeneH5[geneCount1]+
            (11./84.)*exprxGeneH6[geneCount1])*degMultiplier);
        }
      }
      double max_diff_o4_o5=0;
      /////////////////////////////////////////////////////////////
      for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
      {
        exprxGeneH[geneCount1]=exprxGene[geneCount1]+
          (5179./57600.)*exprxGeneH1[geneCount1]+
          (7571./16695.)*exprxGeneH3[geneCount1]+
          (393./640.)*exprxGeneH4[geneCount1]+
          (-92097./339200.)*exprxGeneH5[geneCount1]+
          (187./2100.)*exprxGeneH6[geneCount1]+
          (1./40.)*exprxGeneH7[geneCount1];
        if(exprxGeneH[geneCount1]<0) exprxGeneH[geneCount1]=0;

        exprxGene[geneCount1]=exprxGene[geneCount1]+
          (35./384.)*exprxGeneH1[geneCount1]+
          (500./1113.)*exprxGeneH3[geneCount1]+
          (125./192.)*exprxGeneH4[geneCount1]+
          (-2187./6784.)*exprxGeneH5[geneCount1]+
          (11./84.)*exprxGeneH6[geneCount1];

        if(exprxGene[geneCount1]<0) isNegative = true;

        double diff_o4_o5=exprxGene[geneCount1] -
          exprxGeneH[geneCount1];

        diff_o4_o5 = diff_o4_o5 >= 0 ? diff_o4_o5 : -diff_o4_o5;
        max_diff_o4_o5 = max_diff_o4_o5 > diff_o4_o5 ? max_diff_o4_o5 :
          diff_o4_o5;

      }
      
      double s_rk = h*rkTolerance/(2*(testTime)*max_diff_o4_o5);

      s_rk=std::pow(s_rk,0.25);
      //Rcout<<s_rk<<"\n";
      if(s_rk<1) {
        if(h>0.00001)
        {h=0.5*h;}
        //Rcout<<"decrease h\n";
      }
      else {
        //Rcout<<"increase h\n";
        i+=h;
        if(s_rk>2) {if(h<0.5) {h=2*h;}}

      }

      for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
      {
        exprxGeneH[geneCount1]=exprxGene[geneCount1];
      }

      //i=i+h;
    }


      double test_delta;
      test_delta = sum_delta(prevExprxGene, exprxGene, numberGene);
      if (test_delta < convergThresh){
        for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
        {
          outGE<<std::setprecision(outputPrecision)
          <<exprxGene[geneCount1]<<"\t";
        }
        isConverged = true;  
        outConv<<isConverged<<"\t" << testIter <<"\n";
          break;
      }
      else if(isNegative){
        for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
          {
          outGE<<std::setprecision(outputPrecision)
          <<exprxGene[geneCount1]<<"\t";
          }
        isConverged = true;  //Necessary to stop double reporting
        outConv<<3<<"\t" << testIter <<"\n";
        break;
      }

  };
  if(isConverged == false){
  for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
  {
    outGE<<std::setprecision(outputPrecision)
    <<exprxGene[geneCount1]<<"\t";
  }
  outConv<<isConverged<<"\t" << numConvergenceIter <<"\n";
}
}