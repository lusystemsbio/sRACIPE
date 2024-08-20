#include "header.h"
#include <Rcpp.h>
using namespace Rcpp;


//This function calculates the squared Euclidean norm of a vector.
double cal_norm(
    std::vector<double> &vector,
    const int &number_genes
){
    double ssq = 0.0;
    for(int i=0;i<number_genes;i++){
        ssq += pow(vector[i], 2);
    }
    return ssq;
}

//This function calculates derivatives using the Euler method. This allows for 
//limit cycle detection without having to use the main steppers or integrating them
//into the limit cycle functionality.
void cal_fX(const int &number_gene,
            const double &signalRate,
            const Rcpp::IntegerMatrix geneInteraction,
            const std::vector<double> &g_gene,
            const std::vector<double> &k_gene,
            const std::vector<std::vector<int> > &n_gene,
            const std::vector<std::vector<double> > &lambda_gene,
            const std::vector<std::vector<double> > &threshold_gene_log,
            const Rcpp::NumericVector &geneTypes,
            std::vector<double> &fX_arr,
            std::vector <double> &exprxGene
){
    double fX2_arr[number_gene]; //fX2_arr saves fX values for degradation part

    //fX=Gx:
    //fX2=kX*X:
    for (int i=0;i<number_gene;++i) {
        fX_arr[i]=g_gene[i]; 
        fX2_arr[i]=k_gene[i]*exprxGene[i]; 
    }

    for(int geneCount1=0;geneCount1<number_gene;geneCount1++){
        double growthMultiplier=1;
        double degMultiplier=1;
        for(int geneCount2=0;geneCount2<number_gene;geneCount2++)
        {
            double geneValue=exprxGene[geneCount2];
            double geneThreshold=threshold_gene_log[geneCount1][geneCount2];
            int geneN=n_gene[geneCount1][geneCount2];
            double geneLambda=lambda_gene[geneCount1][geneCount2];
            calcMultiplier(geneCount1, geneCount2, growthMultiplier, degMultiplier,
                            geneValue, geneInteraction, geneN, geneLambda,
                            geneThreshold);
        }
        if (geneTypes[geneCount1] == 2){
            growthMultiplier = growthMultiplier*signalRate;
            degMultiplier = degMultiplier*signalRate;
        }
        fX_arr[geneCount1] *= growthMultiplier;
        fX2_arr[geneCount1] *= degMultiplier;
    }

    //adjust fX by kX*X:
    for(int i=0; i<number_gene; i++){
        fX_arr[i]-=fX2_arr[i];
    }

}


double cal_limitcycle(int number_gene,
            const double &LCSimStepSize,
            int LCSimSteps,
            const long double &convergThresh,
            const Rcpp::IntegerMatrix geneInteraction,
            const std::vector<double> &g_gene,
            const std::vector<double> &k_gene,
            const std::vector<std::vector<int> > &n_gene,
            const std::vector<std::vector<double> > &lambda_gene,
            const std::vector<std::vector<double> > &threshold_gene_log,
            const double &signalRate,
            const Rcpp::NumericVector &geneTypes,
            std::vector<double> start_exp_arr,
            std::vector<std::vector<double>> LC_exp_arr
){
    //variables
    std::vector<double> fX_arr(number_gene); //fX_arr saves deriv values in each Euler step
    std::vector<double> curr_exp_arr(number_gene);
    std::vector<double> next_exp_arr(number_gene);
    double max_dist = 0.0;

    //copy the start expressions to a local variable:
    for (int i=0;i<number_gene;++i) {
        curr_exp_arr[i]=start_exp_arr[i]; 
        LC_exp_arr[0][i]=curr_exp_arr[i];
    }

    int countStep=0;
    while(countStep<LCSimSteps){
        cal_fX(number_gene,
                signalRate,
                geneInteraction,
                g_gene, k_gene,
                n_gene, lambda_gene,
                threshold_gene_log,
                geneTypes, fX_arr,
                curr_exp_arr);

        //calculate next expressions:
        for (int i=0;i<number_gene;++i) {
             next_exp_arr[i]=curr_exp_arr[i]+fX_arr[i]*LCSimStepSize;
        }        

        double tmp_dist=sum_delta(curr_exp_arr,next_exp_arr,number_gene);
        if(tmp_dist>max_dist) max_dist=tmp_dist;

        countStep+=1;
        //copy next exp to curr exp and save curr exp
        for (int i=0;i<number_gene;++i){
            curr_exp_arr[i]=next_exp_arr[i]>=0?next_exp_arr[i]:0;
            LC_exp_arr[countStep][i]=curr_exp_arr[i];
        }
    }
    return max_dist;
}


int cal_period(const int &number_gene,
            const int &MaxPeriods,
            int countMinExp,
            std::vector<int> &minIdxArr,
            std::vector<std::vector<double> > min_exp_arr,
            std::vector<double> &LC_start_exp_arr,
            const int &NumSampledPeriods,
            const int &AllowedPeriodError,
            const double &SamePointProximity
){
    // Local Variables
    int period=0;
    int countSameMinDist=0;
    int countPeriod=0; //counts number of periods found
    int periodT[MaxPeriods]; //saves the periods
    int idx_lastPeriod = countMinExp - 1; //keeps track of the index of the last period
    double diffExp = 0.0; //saves diff between expr levels at two time points
    bool found = false; //marks whether period is found

    //check in backward direction starting from
    //the last element in min_dist_arr:
    for (int i=countMinExp-2;i>=0;i--){
       diffExp=sum_delta(min_exp_arr[countMinExp-1],
                          min_exp_arr[i],number_gene);
       if(diffExp<=SamePointProximity){
          countSameMinDist++;
          periodT[countPeriod]=minIdxArr[idx_lastPeriod]-minIdxArr[i];
          //save position of the valley where the last period was found:
          idx_lastPeriod=i;
          countPeriod++; //increment count of the periods found
       }
    }
    if(countPeriod>=NumSampledPeriods){
        found = true;
        //check whether each of these periods are 
        //within the allowed error limit:
        for (int i=countPeriod-1;i>=1;i--){
            if(abs(periodT[i]-periodT[i-1])>AllowedPeriodError){
                found = false;
                break;
            }
        }
    }
    if(found){
        //copy the last sampled period to period varialbe for returning:
        period=periodT[countPeriod-1];
        //copy the last minimum expression to the 
        //LC start variable for returning:
        for(int i=0;i<number_gene;i++){
           LC_start_exp_arr[i]=min_exp_arr[countMinExp-1][i];
        }
    }
    return period;

}

int detect_limitcycle(const int &number_gene,
                    const double &LCSimStepSize,
                    int LCSimSteps,
                    const long double &convergThresh,
                    const Rcpp::IntegerMatrix geneInteraction,
                    const std::vector<double> &g_gene,
                    const std::vector<double> &k_gene,
                    const std::vector<std::vector<int> > &n_gene,
                    const std::vector<std::vector<double> > &lambda_gene,
                    const std::vector<std::vector<double> > &threshold_gene_log,
                    const double &signalRate,
                    const Rcpp::NumericVector &geneTypes,
                    const int &LCIter,
                    const int &MaxPeriods,
                    const int &NumSampledPeriods,
                    const int &AllowedPeriodError,
                    const double &SamePointProximity,
                    std::vector<double> start_exp_arr,
                    std::vector<double> LC_start_exp_arr
){

    //variables
    std::vector<double> fX_arr(number_gene); //fX_arr saves deriv values in each Euler step
    std::vector<double> prev_exp_arr(number_gene);
    std::vector<double> curr_exp_arr(number_gene);
    std::vector<double> next_exp_arr(number_gene);

    //storage for expressions at each peak:
    std::vector<std::vector<double>> max_exp_arr(MaxPeriods, std::vector<double>(number_gene));
    //storage for expressions at each valley:
    std::vector<std::vector<double>> min_exp_arr(MaxPeriods, std::vector<double>(number_gene));
    //storage for index at each peak:
    std::vector<int> maxIdxArr(MaxPeriods);
    //storage for index at each valley:
    std::vector<int> minIdxArr(MaxPeriods);
    //store for distances from the start to each valley:
    std::vector<double> min_dist_arr(MaxPeriods);

    //counts for peaks and valleys:
    int countMaxExp, countMinExp;

    double prev_dist, curr_dist; //keeps track of distance from the start
    bool moving_uphill = true; //keeps track of simulation direction

    int count_exp = 0; //global count for simulation steps
    int period = 0; //period for the LC

    //copy the start expressions to a local variable:
    for (int i=0;i<number_gene;++i) {
        curr_exp_arr[i]=start_exp_arr[i]; 
        prev_exp_arr[i]=curr_exp_arr[i]; 
    }

    prev_dist = 0.0;
    curr_dist = prev_dist;
    countMaxExp = 0;
    countMinExp = countMaxExp;



    //calculate ODE-outer loop:
    for(int countIter=0; countIter<LCIter; countIter++){

        //Calculate ODE-inner loop:
        int countStep = 0;
        while(countStep<LCSimSteps){
            cal_fX(number_gene,
                signalRate,
                geneInteraction,
                g_gene, k_gene,
                n_gene, lambda_gene,
                threshold_gene_log,
                geneTypes, fX_arr,
                curr_exp_arr);

            //calculate next expressions:
            for (int i=0;i<number_gene;++i) {
                next_exp_arr[i]=curr_exp_arr[i]+fX_arr[i]*LCSimStepSize;
            }

            curr_dist=sum_delta(start_exp_arr,next_exp_arr,number_gene);
            if(moving_uphill){ //while moving uphill:
                if (curr_dist<prev_dist){
                    //hit the peak already. now, switch to downhill:
                    moving_uphill=false;
                }
            }
            else { //while moving downhill:
                if (curr_dist>prev_dist){
                    //hit the bottom already. now, switch to uphill:
                    moving_uphill=true;
                    //save information about this valley:
                    for (int i=0;i<number_gene;++i){
                        min_exp_arr[countMinExp][i]=curr_exp_arr[i];
                    }
                    minIdxArr[countMinExp]=count_exp;
                    min_dist_arr[countMinExp]=prev_dist;
                    countMinExp++;
                }
            }//end of the block for moving downhill
            if((countMaxExp>=MaxPeriods) || (countMinExp>=MaxPeriods)) break;

            prev_dist=curr_dist; //uphill motion
            //copy next exp to curr exp and save curr exp
            for (int i=0;i<number_gene;++i){
                curr_exp_arr[i]=next_exp_arr[i]>=0?next_exp_arr[i]:0;
            }
            count_exp++;
            countStep+=1;
        }//end of while(countStep<LCSimSteps)

        double fX_norm=cal_norm(fX_arr,number_gene);
        //if a stable state is found or a boundary situation is met,
        //then return -1:
        if(fX_norm<=convergThresh) return -1;

        period = cal_period(
                    number_gene,
                    MaxPeriods,
                    countMinExp,
                    minIdxArr,
                    min_exp_arr,
                    LC_start_exp_arr,
                    NumSampledPeriods,
                    AllowedPeriodError,
                    SamePointProximity);
        //if nonzero period found, then break of the loop:
        if (period!=0) break;
    }

    return period;

}


int find_limitcycles(std::vector<std::vector<double> > &exprxGene,
             std::ofstream &outLC,
             const size_t &number_gene,
             const size_t &nIC,
             const Rcpp::IntegerMatrix geneInteraction,
             const std::vector<double> &g_gene,
             const std::vector<double> &k_gene,
             const std::vector<std::vector<int> > &n_gene,
             const std::vector<std::vector<double> > &lambda_gene,
             const std::vector<std::vector<double> > &threshold_gene_log,
             const double &h,
             const double &signalRate,
             const Rcpp::NumericVector &geneTypes,
             const size_t &modelCount,
             const long double &convergThresh,
             const size_t outputPrecision,
             const double &LCSimTime,
             const double &LCSimStepSize,
             const int &maxLCs,
             const int &LCIter,
             const int &MaxPeriods,
             const int &NumSampledPeriods,
             const int &AllowedPeriodError,
             const double &SamePointProximity,
             Rcpp::LogicalVector &convergBool

){
    int LCSimSteps = (int) (LCSimTime/LCSimStepSize);
    int LC_period_arr[maxLCs];
    std::vector<double> start_exp_arr(number_gene); //start point for calculating limit cycle
    std::vector<double> LC_start_exp_arr(number_gene);

    int period=0; //period of the limit cycle
    int countLC=0;

    //initialize LC_period_arr with zeros:
    for (int i=0;i<maxLCs;i++){
        LC_period_arr[i]=0; 
    }

    for (size_t i=0;i<number_gene;i++){
       LC_start_exp_arr[i]=0.0;
    }

    for(size_t i=0; i<nIC; i++){
        if (convergBool[i] == FALSE) continue;
        period = detect_limitcycle(number_gene,
                LCSimStepSize, LCSimSteps,
                convergThresh, geneInteraction, 
                g_gene, k_gene,
                n_gene, lambda_gene,
                threshold_gene_log,
                signalRate, geneTypes,
                LCIter, MaxPeriods,
                NumSampledPeriods,
                AllowedPeriodError,
                SamePointProximity,
                exprxGene[i],
                LC_start_exp_arr
            );
        
        if(period>0){
            //check for repeated limit cycle:
            //NewSimSteps: simulation steps for calculating expression 
            //levels along limit cycle:
            int NewSimSteps;

            //MAX_DIST: used as an allowed error between the converged expression 
            //levels (based on Initial Conditions) and the expression 
            //level along the limit cycle:
            double MAX_DIST; 
            double tmp_norm=0.0;
            bool is_same;

            is_same=false;
            for(int k=0;k<countLC;k++){
                if(abs(period-LC_period_arr[k])<=AllowedPeriodError){
                    is_same=true;
                    break;
                }
            }
            //if it is the same limit cycle, skip subsequent calculation:
            if(is_same) continue;

            LC_period_arr[countLC]=period;
            countLC++;
            std::vector<std::vector<double>> LC_Exp_Arr(period+2, std::vector<double>(number_gene, 0.0));
            NewSimSteps = period + 1;

            MAX_DIST = cal_limitcycle(number_gene,
                    LCSimStepSize, NewSimSteps,
                    convergThresh, geneInteraction,
                    g_gene, k_gene,
                    n_gene, lambda_gene,
                    threshold_gene_log,
                    signalRate, geneTypes,
                    LC_start_exp_arr,
                    LC_Exp_Arr);

            //relax MAX_DIST by factoring it by 10:
            MAX_DIST *= 10.0;

            //change convergBool status to false for all the positions that 
            //all give to this limit cycle:
            for (size_t j=0;j<nIC;j++){
	            if (convergBool[j]==false) continue;
                for(int k=0;k<period+2;k++){
                    tmp_norm=sum_delta(exprxGene[j],LC_Exp_Arr[k],number_gene);
                    if(tmp_norm<=MAX_DIST) {
	                    convergBool[j]=false;
                        break;
                    }
                }
             }

            if(countLC>=maxLCs) break;

            for(int j=0;j<NewSimSteps +1; j++){
                //Write model number and LC number
                outLC<<modelCount<<"\t"<<countLC<<"\t";
                //Write Period
                outLC<<period<<"\t";
                for(size_t k=0; k<number_gene;k++){
                    outLC<<std::setprecision(outputPrecision)
                    <<LC_Exp_Arr[j][k]<<"\t";
                }
            outLC<<"\n";
            }
        }
    }

    return countLC;

}

// [[Rcpp::export]]

int limitcyclesGRC(Rcpp::IntegerMatrix geneInteraction,
    String outFileLC, Rcpp::List config,
    Rcpp::LogicalVector &modelConverg,
    String inFileParams, String inFileGE,
    Rcpp::NumericVector geneTypes
){
    // Count total limit cycles found
    int totLCs;

    // Initialize the network
    size_t numberGene = geneInteraction.ncol();
    
    Rcpp::NumericVector simulationParameters =
        as<NumericVector>(config["simParams"]);
    Rcpp::NumericVector hyperParameters =
        as<NumericVector>(config["hyperParams"]);
    Rcpp::NumericVector LCParameters = 
        as<NumericVector>(config["LCParams"]);

    size_t numModels = static_cast<size_t>(simulationParameters[0]);
    double h = simulationParameters[2];
    size_t nIC = static_cast<size_t> (simulationParameters[4]);
    size_t outputPrecision = static_cast<size_t> (simulationParameters[5]);
    //double rkTolerance = simulationParameters[6];
    long double convergThresh = simulationParameters[9];
    double signalRate = hyperParameters[11];
    double LCSimTime = LCParameters[0];
    double LCSimStepSize = LCParameters[1];
    int maxLCs = LCParameters[2];
    int LCIter = LCParameters[3];
    int MaxPeriods = LCParameters[4];
    int NumSampledPeriods = LCParameters[5];
    int AllowedPeriodError = LCParameters[6];
    double SamePointProximity = LCParameters[7];

    std::string fileNameGE = inFileGE;
    std::string fileNameParam = inFileParams;
    std::string fileNameLC = outFileLC;

    std::ifstream inParams;
    inParams.open(fileNameParam,
                         std::ifstream::in);
    if(!inParams.is_open()) {     Rcout <<fileNameParam
        << "Cannot open input file for reading parameters.\n";  return -1;
      }

    std::ifstream inGE;
    inGE.open(fileNameGE,
                         std::ifstream::in);
    if(!inGE.is_open()) {     Rcout <<fileNameParam
        << "Cannot open input file for reading expressions.\n";  return -1;
      }

    std::ofstream outLC(fileNameLC, std::ios::out);
    if(!outLC.is_open()) {     Rcout << "Cannot open output file.\n";
      return -1;}
    
    for(size_t modelCount=0; modelCount<numModels; modelCount++){
        //Initialize production rate of genes
        std::vector<double> gGene(numberGene);

        //Initialize degradation rate of genes
        std::vector<double> kGene(numberGene);

        //Initialize hill coefficient for each interaction
        std::vector<std::vector<int> >
          nGene(numberGene, std::vector<int>(numberGene));

        //Initialize fold change for each interaction
        std::vector<std::vector<double> >
          lambdaGene(numberGene, std::vector<double>(numberGene));

        //Initialize threshold for each interaction
        std::vector<std::vector<double> >
          threshGeneLog(numberGene, std::vector<double>(numberGene));

        readParameters( geneInteraction, numberGene, gGene,
                    kGene, nGene,
                    lambdaGene,
                    threshGeneLog, inParams);
        
        std::vector<std::vector<double> > 
            exprxGene(nIC, std::vector<double>(numberGene));
        for(size_t i=0; i<nIC; i++){
            for(size_t j=0; j<numberGene; j++){
                inGE >> exprxGene[i][j];
            }
        }

        Rcpp::LogicalVector convergBool(nIC);
        for(size_t i=0; i<nIC; i++){
            convergBool[i] = modelConverg[i + modelCount];
        }

        int count;
        count = find_limitcycles(exprxGene,
                                outLC, numberGene, nIC,
                                geneInteraction,
                                gGene, kGene, nGene, lambdaGene,
                                threshGeneLog, h, signalRate, geneTypes,
                                modelCount, convergThresh, outputPrecision, LCSimTime, 
                                LCSimStepSize, maxLCs, LCIter,
                                MaxPeriods, NumSampledPeriods, 
                                AllowedPeriodError, SamePointProximity,
                                convergBool);
        totLCs += count;
    }

    inGE.close();
    inParams.close();
    outLC.close();

    return totLCs;
}