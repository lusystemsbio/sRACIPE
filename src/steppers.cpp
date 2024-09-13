#include "header.h"
#include <Rcpp.h>
using namespace Rcpp;



void calcMultiplier(const int& geneCount1, const int& geneCount2,
                    double& growthMultiplier,
                    double& degMultiplier,
                    double& geneValue,
                    IntegerMatrix geneInteraction,
                    const int& geneN,
                    double& geneLambda,
                    const double& geneThreshold
){

    double geneActMultiplier=1;
    double geneActMultiplier2=1;

    switch(geneInteraction(geneCount1,geneCount2))
    {
    case 0:
      geneActMultiplier=1.0;
      geneActMultiplier2=1.0;
      break;

    case 1:
      geneActMultiplier=(geneLambda+(1.-geneLambda)*
        1./(1.+std::pow((geneValue/geneThreshold),geneN)))/geneLambda;
      geneActMultiplier2=1.0;
      break;

    case 2:
      geneLambda = 1./geneLambda;
      geneActMultiplier = geneLambda+(1.-geneLambda)*
        1./(1.+std::pow((geneValue/geneThreshold),geneN));
        geneActMultiplier2=1.0;
      break;

    case 3:
      geneLambda = 1./geneLambda;
      geneActMultiplier2 = geneLambda+(1.-geneLambda)*
        1./(1.+std::pow((geneValue/geneThreshold),geneN));
        geneActMultiplier=1.0;
      break;

    case 4:
      geneActMultiplier2=(geneLambda+(1.-geneLambda)*
        1./(1.+std::pow((geneValue/geneThreshold),geneN)));
        geneActMultiplier=1.0;
      break;

    case 5:
      geneActMultiplier=(geneLambda+(1.-geneLambda)*
        1./(1.+std::pow((geneValue/geneThreshold),geneN)))/geneLambda;
      geneActMultiplier2=1.0;
      break;

    case 6:
      geneLambda = 1./geneLambda;
      geneActMultiplier = geneLambda+(1.-geneLambda)*
        1./(1.+std::pow((geneValue/geneThreshold),geneN));
        geneActMultiplier2=1.0;
      break;  

    default :
      Rcout << "Invalid Interation code for Gene"<<geneCount1
            <<" and gene"<<geneCount2<<" interaction"<<"\n";
    }

    growthMultiplier=growthMultiplier*geneActMultiplier;
    degMultiplier=degMultiplier*geneActMultiplier2;
}

//This function does step interpolation given a time value and a list of steps
void calcSigValues(const std::vector<double> &timePoints,
                                    const double &t,
                                    const Rcpp::NumericVector &signalVals,
                                    double &currVal){
    
    int nVals = signalVals.size();
    for(int i=0; i<nVals; i++){
      if(t <= timePoints[i+1]){
          currVal = signalVals[i] + ((t-timePoints[i])/(timePoints[i+1]-timePoints[i]))
            *(signalVals[i+1]-signalVals[i]);
          break;
      }
    }
}

void stepEM( std::vector <double> &exprxGene,
             std::ofstream &outGE,
             const double &totTime,
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
             const double &printStart, const double &printInterval,
             const double &D,
             const double &h,
             const double &signalRate,
             const NumericVector &geneTypes,
             const bool &isTimeVarying,
             const std::vector<double> &timePoints,
             const Rcpp::NumericMatrix &signalVals,
             const Rcpp::NumericVector &signalingTypes){

  double exprxGeneH[numberGene]; //array for temp gene expression values
  for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
  {
    exprxGeneH[geneCount1] = exprxGene[geneCount1];
  }

  double i=0.0;
  double printTime = printStart;

  //These vectors are used for time varying param signaling
  std::vector<double> gMults(numberGene, 1);
  std::vector<double> kMults(numberGene, 1);

  do
  {
    i+=h;
    if(isTimeVarying){ //Apply time-varying gene parm changes
      int colIdx = 1;
      for(int geneCount1=0;geneCount1<numberGene;geneCount1++){
          if(signalingTypes[geneCount1] == 1){ //g signaling
            calcSigValues(timePoints, i, 
              signalVals( _ , colIdx), gMults[geneCount1]);
            colIdx++;
          }
          else if(signalingTypes[geneCount1] == 2){ //k signaling
            calcSigValues(timePoints, i, 
              signalVals( _ , colIdx), kMults[geneCount1]);
            colIdx++;
          }
      }
    }

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
        h*(gGene[geneCount1]*growthMultiplier*gMults[geneCount1]
        -kGene[geneCount1]*kMults[geneCount1]*
        exprxGene[geneCount1]*degMultiplier) +
        D*sqrt(h)*g_distribution(g_generator)*Darray[geneCount1]+
        D_shot_scaling*D*sqrt(h)*g_distribution(g_generator)*
        Darray[geneCount1]*exprxGene[geneCount1];
      if(exprxGeneH[geneCount1]<0) exprxGeneH[geneCount1]=0;
    }

    for(int geneCount1=0;geneCount1<numberGene;geneCount1++){
      exprxGene[geneCount1]=exprxGeneH[geneCount1];}

    if((i> printTime) &&
       (i <= (printTime + printInterval)))
    {
      printTime +=printInterval;
      for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
      {
        outGE<<std::setprecision(outputPrecision)
        <<exprxGene[geneCount1]<<"\t";
      }
      //outGE<<"\n";
    }
  }while(i<totTime);

  // for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
  // {
  //   outGE<<std::setprecision(outputPrecision)<<exprxGene[geneCount1]<<"\t";
  // }

}

void stepEM_OU( std::vector <double> &exprxGene,
             std::ofstream &outGE,
             const double &totTime,
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
             const double &printStart, const double &printInterval,
             const double &D,
             const double &h,
             const double &ouNoise_t,
             const double &signalRate,
             const Rcpp::NumericVector &geneTypes,
             const bool &isTimeVarying,
             const std::vector<double> &timePoints,
             const Rcpp::NumericMatrix &signalVals,
             const Rcpp::NumericVector &signalingTypes){

  double exprxGeneH[numberGene]; //array for temp gene expression values
  double currNoise[numberGene]; //array for temp gene expression values
  std::vector <double> prevNoise(numberGene, 0.0); //array for current noise

  for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
  {
    exprxGeneH[geneCount1] = exprxGene[geneCount1];
    prevNoise[geneCount1] = D*Darray[geneCount1] * g_distribution(g_generator);
  }

  double i=0.0;
  double printTime = printStart;
  
  //These vectors are used for time varying param signaling
  std::vector<double> gMults(numberGene, 1);
  std::vector<double> kMults(numberGene, 1);

  do
  {
    i+=h;
    if(isTimeVarying){ //Apply time-varying gene parm changes
      int colIdx = 1;
      for(int geneCount1=0;geneCount1<numberGene;geneCount1++){
          if(signalingTypes[geneCount1] == 1){ //g signaling
            calcSigValues(timePoints, i, 
              signalVals( _ , colIdx), gMults[geneCount1]);
            colIdx++;
          }
          else if(signalingTypes[geneCount1] == 2){ //k signaling
            calcSigValues(timePoints, i, 
              signalVals( _ , colIdx), kMults[geneCount1]);
            colIdx++;
          }
      }
    }

    for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
    {
      double growthMultiplier=1;
      double degMultiplier=1;
      currNoise[geneCount1] = prevNoise[geneCount1] * exp(-h/ouNoise_t) + 
        D*Darray[geneCount1] * sqrt(1-exp(-2*h/ouNoise_t)) * g_distribution(g_generator);

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
        h*(gGene[geneCount1]*growthMultiplier*gMults[geneCount1]
        -kGene[geneCount1]*kMults[geneCount1]*
        exprxGene[geneCount1]*degMultiplier) +
        h*currNoise[geneCount1]+
        D_shot_scaling*D*sqrt(h)*g_distribution(g_generator)*
        Darray[geneCount1]*exprxGene[geneCount1];
      if(exprxGeneH[geneCount1]<0) exprxGeneH[geneCount1]=0;
    }

    for(int geneCount1=0;geneCount1<numberGene;geneCount1++){
      exprxGene[geneCount1]=exprxGeneH[geneCount1];
      prevNoise[geneCount1]=currNoise[geneCount1];}

    if((i> printTime) &&
       (i <= (printTime + printInterval)))
    {
      printTime +=printInterval;
      for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
      {
        outGE<<std::setprecision(outputPrecision)
        <<exprxGene[geneCount1]<<"\t";
      }
      //outGE<<"\n";
    }
  }while(i<totTime);

  // for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
  // {
  //   outGE<<std::setprecision(outputPrecision)<<exprxGene[geneCount1]<<"\t";
  // }

}

///////////////////////////////////////////////////////////////////////////////

//Runge Kutta fourth order
///////////////////////////////////////////////////////////////////////////////
void stepRK4( std::vector <double> &exprxGene,
             std::ofstream &outGE,
             const double &totTime,
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
             const double &printStart, const double &printInterval,
             const double &h,
             const double &signalRate,
             const NumericVector &geneTypes,
             const bool &isTimeVarying,
             const std::vector<double> &timePoints,
             const Rcpp::NumericMatrix &signalVals,
             const Rcpp::NumericVector &signalingTypes){

double exprxGeneH1[numberGene]; //array for temp gene expression values
double exprxGeneH2[numberGene]; //array for temp gene expression values
double exprxGeneH3[numberGene]; //array for temp gene expression values
double exprxGeneH4[numberGene]; //array for temp gene expression values

for(int geneCountTmp=0;geneCountTmp<numberGene;geneCountTmp++)
{
  exprxGeneH1[geneCountTmp]=exprxGene[geneCountTmp];
  exprxGeneH2[geneCountTmp]=exprxGene[geneCountTmp];
  exprxGeneH3[geneCountTmp]=exprxGene[geneCountTmp];
  exprxGeneH4[geneCountTmp]=exprxGene[geneCountTmp];
}
double i=0.0;
double printTime = printStart;

//These vectors are used for time varying param signaling
std::vector<double> gMults(numberGene, 1);
std::vector<double> kMults(numberGene, 1);

do
{
  i+=h;
  if(isTimeVarying){ //Apply time-varying gene parm changes
    int colIdx = 1;
    for(int geneCount1=0;geneCount1<numberGene;geneCount1++){
        if(signalingTypes[geneCount1] == 1){ //g signaling
          calcSigValues(timePoints, i, 
            signalVals( _ , colIdx), gMults[geneCount1]);
          colIdx++;  
        }
        else if(signalingTypes[geneCount1] == 2){ //k signaling
          calcSigValues(timePoints, i, 
            signalVals( _ , colIdx), kMults[geneCount1]);
          colIdx++;
        }
    }
  }
  
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

    exprxGeneH1[geneCount1]=h*(gGene[geneCount1]*growthMultiplier
      *gMults[geneCount1] -kGene[geneCount1]*kMults[geneCount1]*
      exprxGene[geneCount1]*degMultiplier);
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
      growthMultiplier*gMults[geneCount1]-
      kGene[geneCount1]*kMults[geneCount1]*(exprxGene[geneCount1] +
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

    exprxGeneH3[geneCount1]=h*((gGene[geneCount1])*growthMultiplier
      *gMults[geneCount1]-kGene[geneCount1]*
      kMults[geneCount1]*(exprxGene[geneCount1] +
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

    exprxGeneH4[geneCount1]=h*((gGene[geneCount1])*growthMultiplier*
      gMults[geneCount1]-kGene[geneCount1]*
      kMults[geneCount1]*(exprxGene[geneCount1]+
      exprxGeneH3[geneCount1])*degMultiplier);
  }


  /////////////////////////////////////////////////////////////
  for(int geneCount1=0;geneCount1<numberGene;geneCount1++){
    exprxGene[geneCount1] = exprxGene[geneCount1]+
      (exprxGeneH1[geneCount1]+2*exprxGeneH2[geneCount1]+
      2*exprxGeneH3[geneCount1]+exprxGeneH4[geneCount1])/6;
    if(exprxGene[geneCount1]<0) exprxGene[geneCount1]=0;
    }


  if((i> printTime) &&
     (i <= (printTime + printInterval)))
  {
    printTime +=printInterval;
    //std::cout<<i<<"\n";
    for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
    {
      outGE<<std::setprecision(outputPrecision)
            <<exprxGene[geneCount1]<<"\t";
    }
    //outGE<<"\n";
  }
}while(i<totTime);
}



void stepDP( std::vector <double> &exprxGene,
              std::ofstream &outGE,
              const double &totTime,
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
              const double &printStart, const double &printInterval,
              double h, const double &rkTolerance,
              const double &signalRate,
              const NumericVector &geneTypes,
              const bool &isTimeVarying,
              const std::vector<double> &timePoints,
              const Rcpp::NumericMatrix &signalVals,
              const Rcpp::NumericVector &signalingTypes){
  double exprxGeneH[numberGene]; //array for temp gene expression values
  double exprxGeneH1[numberGene]; //array for temp gene expression values
  double exprxGeneH2[numberGene]; //array for temp gene expression values
  double exprxGeneH3[numberGene]; //array for temp gene expression values
  double exprxGeneH4[numberGene]; //array for temp gene expression values
  double exprxGeneH5[numberGene]; //array for temp gene expression values
  double exprxGeneH6[numberGene]; //array for temp gene expression values
  double exprxGeneH7[numberGene]; //array for temp gene expression values

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
  int printCounter = 0;

  
  //These vectors are used for time varying param signaling
  std::vector<double> gMults(numberGene, 1);
  std::vector<double> kMults(numberGene, 1);

  do
  {
    if(isTimeVarying){ //Apply time-varying gene parm changes
      int colIdx = 1;
      for(int geneCount1=0;geneCount1<numberGene;geneCount1++){
          if(signalingTypes[geneCount1] == 1){ //g signaling
            calcSigValues(timePoints, i, 
              signalVals( _ , colIdx), gMults[geneCount1]);
            colIdx++;
          }
          else if(signalingTypes[geneCount1] == 2){ //k signaling
            calcSigValues(timePoints, i, 
              signalVals( _ , colIdx), kMults[geneCount1]);
            colIdx++;
          }
      }
    }

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

      exprxGeneH1[geneCount1]=h*(gGene[geneCount1]*growthMultiplier
        *gMults[geneCount1]-kGene[geneCount1]
        *kMults[geneCount1]*exprxGene[geneCount1]*degMultiplier);
    }

    for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
    {
      double growthMultiplier=1;
      double degMultiplier=1;

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

      exprxGeneH2[geneCount1]=h*((gGene[geneCount1])*growthMultiplier
        *gMults[geneCount1]- kGene[geneCount1]*kMults[geneCount1]
        *(exprxGene[geneCount1] +
        0.2*exprxGeneH1[geneCount1])*degMultiplier);
    }

    for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
    {
      double growthMultiplier=1;
      double degMultiplier=1;

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

      exprxGeneH3[geneCount1]=h*((gGene[geneCount1])*growthMultiplier
        *gMults[geneCount1]-kGene[geneCount1]*kMults[geneCount1]
        *(exprxGene[geneCount1]+(0.25*exprxGeneH1[geneCount1]+
        0.75*exprxGeneH2[geneCount1])*0.3)*degMultiplier);
    }


    for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
    {
      double growthMultiplier=1;
      double degMultiplier=1;

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

      exprxGeneH4[geneCount1]=h*((gGene[geneCount1])*growthMultiplier
        *gMults[geneCount1]-kGene[geneCount1]*kMults[geneCount1]
        *(exprxGene[geneCount1]+
        0.8*(((11./9.)*exprxGeneH1[geneCount1]+
        (-14./3.)*exprxGeneH2[geneCount1])+
        (40./9.)*exprxGeneH3[geneCount1]))*degMultiplier);
    }

    for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
    {
      double growthMultiplier=1;
      double degMultiplier=1;

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

      exprxGeneH5[geneCount1]=h*((gGene[geneCount1])*growthMultiplier
        *gMults[geneCount1]-kGene[geneCount1]*kMults[geneCount1]
        *(exprxGene[geneCount1]+(8./9.)*
        ((4843./1458.)*exprxGeneH1[geneCount1]+
        (-3170./243.)*exprxGeneH2[geneCount1]+
        (8056./729.)*exprxGeneH3[geneCount1]+
        (-53./162.)*exprxGeneH4[geneCount1]))*degMultiplier);
    }


    for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
    {
      double growthMultiplier=1;
      double degMultiplier=1;

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

      exprxGeneH6[geneCount1]=h*((gGene[geneCount1])*growthMultiplier
        *gMults[geneCount1]-kGene[geneCount1]*kMults[geneCount1]
        *(exprxGene[geneCount1]+
        (9017./3168.)*exprxGeneH1[geneCount1]+
        (-355./33.)*exprxGeneH2[geneCount1]+
        (46732./5247.)*exprxGeneH3[geneCount1]+
        (49./176.)*exprxGeneH4[geneCount1]+
        (-5103./18656.)*exprxGeneH5[geneCount1])*degMultiplier);
    }



    for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
    {
      double growthMultiplier=1;
      double degMultiplier=1;

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

      exprxGeneH7[geneCount1]=h*((gGene[geneCount1])*growthMultiplier
        *gMults[geneCount1]-kGene[geneCount1]*kMults[geneCount1]
        *(exprxGene[geneCount1]+
        (35./384.)*exprxGeneH1[geneCount1]+
        (500./113.)*exprxGeneH3[geneCount1]+
        (125./192.)*exprxGeneH4[geneCount1]+
        (-2187./6784.)*exprxGeneH5[geneCount1]+
        (11./84.)*exprxGeneH6[geneCount1])*degMultiplier);
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

      if(exprxGene[geneCount1]<0) exprxGene[geneCount1]=0;

      double diff_o4_o5=exprxGene[geneCount1] -
        exprxGeneH[geneCount1];

      diff_o4_o5 = diff_o4_o5 >= 0 ? diff_o4_o5 : -diff_o4_o5;
      max_diff_o4_o5 = max_diff_o4_o5 > diff_o4_o5 ? max_diff_o4_o5 :
        diff_o4_o5;

    }
    double s_rk = h*rkTolerance/(2*(totTime)*max_diff_o4_o5);

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
    if((i> (printStart + printInterval*printCounter)) &&
       i <= (h+printStart + printInterval*printCounter))
    {

      printCounter++;
      //std::cout<<i<<"\n";
      for(int geneCount1=0;geneCount1<numberGene;geneCount1++)
      {
        outGE<<std::setprecision(outputPrecision)
              <<exprxGene[geneCount1]<<"\t";
      }
      //outGE<<"\n";
    }
    //i=i+h;
  }while(i<totTime);

}