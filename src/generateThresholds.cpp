#include"header.h"
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;


// [[Rcpp::export]]
size_t generateThresholds(
    const Rcpp::IntegerMatrix geneInteraction,
    Rcpp::NumericVector thresholdGene,
    Rcpp::List config)
{
    unsigned int seed =  static_cast<unsigned int>(Rcpp::sample(32000,1,true)(0));
    std::mt19937_64 u_generator (seed);
    
    NumericVector simulationParameters = as<NumericVector>(config["simParams"]);
    NumericVector stochasticParameters = as<NumericVector>(config["stochParams"]);
    NumericVector hyperParameters = as<NumericVector>(config["hyperParams"]);
    LogicalVector options = as<LogicalVector>(config["options"]);

  size_t modelCountMax = static_cast<size_t>(simulationParameters[0]);

  double gMin = hyperParameters[0];
  double gMax = hyperParameters[1];
  double kMin = hyperParameters[2];
  double kMax = hyperParameters[3];
  double lambdaMin = hyperParameters[4];
  double lambdaMax = hyperParameters[5];
  size_t nMin = static_cast<size_t>(hyperParameters[6]);
  size_t nMax = static_cast<size_t>(hyperParameters[7]);
  size_t interactionTypes = static_cast<size_t>(hyperParameters[8]);
  size_t thresholdMax = static_cast<size_t>(hyperParameters[9]);
  double sdFactor = hyperParameters[10];
  double lambdaMinDeg = hyperParameters[12];
  double lambdaMaxDeg = hyperParameters[13];

// Generate Thresholds
  size_t numberGene=geneInteraction.nrow();
  Rcout<<"generating thresholds for uniform distribution1..."<<"\n";
  double geneIsolatedMedian=(gMin+gMax)/(kMin+kMax);
  for(size_t i=0;i<numberGene;i++){thresholdGene[i]= geneIsolatedMedian;}

  size_t interactionsNumber[numberGene][interactionTypes];
  for(size_t i=0;i<numberGene;i++)for(size_t j=0;j<interactionTypes;j++)
    interactionsNumber[i][j]=0;

  for(size_t it1=0;it1<numberGene;it1++){for(size_t it2=0;it2<numberGene;it2++)
  {size_t g_1=geneInteraction(it1,it2);
    interactionsNumber[it1][g_1]+=1;}}

  const long modelCountMax2=std::max(modelCountMax,thresholdMax);
  for(size_t geneCount1=0;geneCount1<numberGene;geneCount1++)
  {
    std::vector<double> Af;
    for(long modelCount=0;modelCount<modelCountMax2;modelCount++)
    {
      double gA=gMin+(gMax-gMin)*u_distribution(u_generator);
      double kA=kMin+(kMax-kMin)*u_distribution(u_generator);
      Af.push_back(gA/kA);
      for(size_t counterInteraction1=0;
          counterInteraction1<interactionsNumber[geneCount1][1];
          counterInteraction1++)
      {
        double gB=gMin+(gMax-gMin)*u_distribution(u_generator);
        double kB=kMin+(kMax-kMin)*u_distribution(u_generator);
        double BA0=geneIsolatedMedian*(1-sdFactor*std::sqrt(3)) +
          (2*std::sqrt(3)*sdFactor*geneIsolatedMedian)
          *u_distribution(u_generator);

        size_t nBA= static_cast<size_t>(u_distribution(u_generator)*
        (nMax-nMin))+nMin;
        double geneLambda=(lambdaMax-lambdaMin)*u_distribution(u_generator) +
          lambdaMin;

        Af[modelCount]=Af[modelCount]*Hs_Racipe(gB/kB, BA0, nBA,
                                                  geneLambda)/geneLambda;
      }

      for(size_t counterInteraction1=0;
          counterInteraction1<interactionsNumber[geneCount1][2];
          counterInteraction1++)
      {
        double gB=gMin+(gMax-gMin)*u_distribution(u_generator);
        double kB=kMin+(kMax-kMin)*u_distribution(u_generator);
        double BA0=geneIsolatedMedian*(1-sdFactor*std::sqrt(3)) +
          (2*std::sqrt(3)*sdFactor*geneIsolatedMedian)
          *u_distribution(u_generator);

        size_t nBA=static_cast<size_t>(u_distribution(u_generator)*
        (nMax-nMin))+nMin;
        double geneLambda=1./((lambdaMax-lambdaMin)*u_distribution(u_generator)+lambdaMin);

        double hillEval=Hs_Racipe(gB/kB, BA0, nBA, geneLambda);
        Af[modelCount] = Af[modelCount]*hillEval;
      }

      for(size_t counterInteraction1=0;
          counterInteraction1<interactionsNumber[geneCount1][3];
          counterInteraction1++)
      {
        double gB=gMin+(gMax-gMin)*u_distribution(u_generator);
        double kB=kMin+(kMax-kMin)*u_distribution(u_generator);
        double BA0=geneIsolatedMedian*(1-sdFactor*std::sqrt(3)) +
          (2*std::sqrt(3)*sdFactor*geneIsolatedMedian)
          *u_distribution(u_generator);

        size_t nBA= static_cast<size_t>(u_distribution(u_generator)*
        (nMax-nMin))+nMin;
        double geneLambda=1./((lambdaMaxDeg-lambdaMinDeg)*u_distribution(u_generator)+lambdaMinDeg);
        Af[modelCount]=Af[modelCount]/Hs_Racipe(gB/kB, BA0, nBA,
                                                  geneLambda);
      }

      for(size_t counterInteraction1=0;
          counterInteraction1<interactionsNumber[geneCount1][4];
          counterInteraction1++)
      {
        double gB=gMin+(gMax-gMin)*u_distribution(u_generator);
        double kB=kMin+(kMax-kMin)*u_distribution(u_generator);
        double BA0=geneIsolatedMedian*(1-sdFactor*std::sqrt(3)) +
          (2*std::sqrt(3)*sdFactor*geneIsolatedMedian)
          *u_distribution(u_generator);

        size_t nBA=static_cast<size_t>(u_distribution(u_generator)*
        (nMax-nMin))+nMin;
        double geneLambda=(lambdaMaxDeg-lambdaMinDeg)*u_distribution(u_generator) +
          lambdaMinDeg;
        double hillEval=Hs_Racipe(gB/kB, BA0, nBA, geneLambda);
        Af[modelCount] = Af[modelCount]/hillEval;
      }

      for(size_t counterInteraction1=0;
          counterInteraction1<interactionsNumber[geneCount1][5];
          counterInteraction1++)
      {
        double gB=gMin+(gMax-gMin)*u_distribution(u_generator);
        double kB=kMin+(kMax-kMin)*u_distribution(u_generator);
        double BA0=geneIsolatedMedian*(1-sdFactor*std::sqrt(3)) +
          (2*std::sqrt(3)*sdFactor*geneIsolatedMedian)
          *u_distribution(u_generator);

        size_t nBA= static_cast<size_t>(u_distribution(u_generator)*
        (nMax-nMin))+nMin;
        double geneLambda=(lambdaMax-lambdaMin)*u_distribution(u_generator) +
          lambdaMin;

        Af[modelCount]=Af[modelCount]*Hs_Racipe(gB/kB, BA0, nBA,
                                                  geneLambda)/geneLambda;
      }

      for(size_t counterInteraction1=0;
          counterInteraction1<interactionsNumber[geneCount1][6];
          counterInteraction1++)
      {
        double gB=gMin+(gMax-gMin)*u_distribution(u_generator);
        double kB=kMin+(kMax-kMin)*u_distribution(u_generator);
        double BA0=geneIsolatedMedian*(1-sdFactor*std::sqrt(3)) +
          (2*std::sqrt(3)*sdFactor*geneIsolatedMedian)
          *u_distribution(u_generator);

        size_t nBA=static_cast<size_t>(u_distribution(u_generator)*
        (nMax-nMin))+nMin;
        double geneLambda=1./((lambdaMax-lambdaMin)*u_distribution(u_generator)+lambdaMin);
        double hillEval=Hs_Racipe(gB/kB, BA0, nBA, geneLambda);
        Af[modelCount] = Af[modelCount]*hillEval;
      }
    }

    std::sort(Af.begin(), Af.end());
    // take the median as threshold
    thresholdGene[geneCount1]=0.5*(Af[static_cast<size_t>((modelCountMax2)*0.5)]
                            + Af[static_cast<size_t>((modelCountMax2+2)*0.5)]);

    Af.clear();
  }
  return 0;

}
