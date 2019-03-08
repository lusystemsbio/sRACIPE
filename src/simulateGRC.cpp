#include "header.h"
#include <Rcpp.h>


using namespace Rcpp;



void writeParameters(const int &number_gene, const int &output_precision,
                     const std::vector<double> &g_gene,
                      const std::vector<double> &k_gene,
                      const std::vector<std::vector<int> > &n_gene,
                      const std::vector<std::vector<double> > &lambda_gene,
                      const std::vector<std::vector<double> > &threshold_gene_log,
                      std::fstream &out_param)
{
  // out_param<<"Parameters"<<"\n";
  /////////////////////////////////////////////////////////////////////////////

  //Writing parameters to file
  /////////////////////////////////////////////////////////////////////////////
  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
  {out_param<<std::setprecision(output_precision)<<g_gene[gene_count1]<<"\t";}
  //production rate of each gene

  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
  {out_param<<std::setprecision(output_precision)<<k_gene[gene_count1]<<"\t";}
  // degradation rate of each gene

  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
    for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
  {if(threshold_gene_log[gene_count1][gene_count2]>0)
    out_param<<std::setprecision(output_precision)
    <<threshold_gene_log[gene_count1][gene_count2]<<"\t";}
  // above--thresholds for the inteaction links, thresholds for
  // inward links for genes are written, starting from gene 1
  // (thresholds for all inward links of gene 1 and so on)
  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
    for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
  {if(n_gene[gene_count1][gene_count2]>0)
    out_param<<std::setprecision(1)<<n_gene[gene_count1][gene_count2]<<"\t";}
  // above--n for the inteaction links,
  // n for inward links for genes are written, starting from gene 1
  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
    for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
  {if(lambda_gene[gene_count1][gene_count2]>0)
    out_param<<std::setprecision(output_precision)
    <<lambda_gene[gene_count1][gene_count2]<<"\t";}
  // above--lambda for the inteaction links, lambda for
  // inward links for genes are written, starting from gene 1
  out_param<<"\n";

}

void readParameters(IntegerMatrix gene_interaction, const int &number_gene,
                    std::vector<double> &g_gene,
                     std::vector<double> &k_gene,
                     std::vector<std::vector<int> > &n_gene,
                     std::vector<std::vector<double> > &lambda_gene,
                     std::vector<std::vector<double> > &threshold_gene_log,
                     std::ifstream &in_parameters)
{
  /////////////////////////////////////////////////////////////////////////////

  //Reading parameters from file
  /////////////////////////////////////////////////////////////////////////////
  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
  {    in_parameters >>  g_gene[gene_count1];}

  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
  { in_parameters >>  k_gene[gene_count1]; }

  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
  {
    for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
    {
      if(gene_interaction(gene_count1,gene_count2)!=0)
      {
        in_parameters >>  threshold_gene_log[gene_count1][gene_count2];
        //Rcout<<"Here"<< threshold_gene_log[gene_count1][gene_count2]<<"\n";
      }
    }
  }

  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
    for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
  {double test;
    if(gene_interaction(gene_count1,gene_count2)!=0) {
      in_parameters >>  test;
      n_gene[gene_count1][gene_count2] = std::round(test);}}

  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
    for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
  {if(gene_interaction(gene_count1,gene_count2)!=0)
    in_parameters >>  lambda_gene[gene_count1][gene_count2];}

}

void selectParameters(Rcpp::IntegerMatrix gene_interaction,
                      Rcpp::NumericVector threshold_gene,
                      const double g_min, const double g_max,
                      const double k_min, const double k_max,
                      const int interaction_types,
                       const long model_count_max,const long threshold_max,
                       const double h, const double lambda_min,
                       const double lambda_max, const int n_min, const int n_max,
                       const double tot_time,
                       const double sd_multiplier, int number_gene,
                       const double D_max,  const double D_shot_scaling,
                       const int scaled_noise,
                       const int D_levels, const double D_scaling,
                       const int output_precision, const bool ANNEALING,
                       const int initial_conditions, String filename,
                       std::vector<double> &g_gene,
                       std::vector<double> &k_gene,
                       std::vector<std::vector<int> > &n_gene,
                       std::vector<std::vector<double> > &lambda_gene,
                       std::vector<std::vector<double> > &threshold_gene_log,
                       std::fstream &out_param)
{
  for(int gene_count1=0;gene_count1<number_gene;gene_count1++){
    g_gene[gene_count1]=g_min+(g_max-g_min)*u_distribution(u_generator);
    }

  for(int gene_count1=0;gene_count1<number_gene;gene_count1++){
    k_gene[gene_count1]=k_min+(k_max-k_min)*u_distribution(u_generator);
    }

  for(int gene_count1=0;gene_count1<number_gene;gene_count1++){
    for(int gene_count2=0;gene_count2<number_gene;gene_count2++){
      if(gene_interaction(gene_count1,gene_count2)==0){
        n_gene[gene_count1][gene_count2]=0;
      } else {
        n_gene[gene_count1][gene_count2] =
          int((n_max-n_min)*u_distribution(u_generator))+n_min;}
    }
  }

  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
  {
    for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
    {
      if(gene_interaction(gene_count1,gene_count2)==0)
      {
        lambda_gene[gene_count1][gene_count2]=0;
      } else {
        lambda_gene[gene_count1][gene_count2]=
        (lambda_max-lambda_min)*u_distribution(u_generator)+lambda_min;}
    }
  }

  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
  {
    for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
    {
      if(gene_interaction(gene_count1,gene_count2)==0)
      {
        threshold_gene_log[gene_count1][gene_count2]=0;
      }
      else
      {

        threshold_gene_log[gene_count1][gene_count2] =
          (1-sd_multiplier*sqrt(3))*threshold_gene[gene_count2] +
          (2*sqrt(3)*sd_multiplier*threshold_gene[gene_count2])
        *u_distribution(u_generator);
      }
    }
  }


}


void selectIcRange(const int number_gene, IntegerMatrix gene_interaction,
                    const std::vector<double> &g_gene,
                    const std::vector<double> &k_gene,
                    const std::vector<std::vector<double> > &lambda_gene,
                    std::vector<double> &max_gene,
                    std::vector<double> &min_gene)
{
  for(int gene_count1=0;gene_count1<number_gene;gene_count1++){
    max_gene[gene_count1]=g_gene[gene_count1]/k_gene[gene_count1];}

  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
  {
    double min_gene_multiplier_final=1;

    for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
    {
      double gene_min_multiplier=1;
      double gene_lambda=lambda_gene[gene_count1][gene_count2];
      switch(gene_interaction(gene_count1,gene_count2))
      {
      case 0:
        gene_min_multiplier=1.0;
        break;

      case 2:
        gene_lambda=1./gene_lambda;
        gene_min_multiplier=gene_lambda;
        break;

      case 1:
        gene_min_multiplier=1./gene_lambda;
        break;
      default :
        Rcout << "Invalid Interation code for Gene"<<gene_count1
        <<" and gene"<<gene_count2<<" interaction"<<"\n";
      }

      min_gene_multiplier_final = min_gene_multiplier_final*gene_min_multiplier;
    }

    min_gene[gene_count1] =
      g_gene[gene_count1]/k_gene[gene_count1]*min_gene_multiplier_final;

  }

}


// [[Rcpp::export]]

int simulateGRCCpp(Rcpp::IntegerMatrix gene_interaction,
                Rcpp::NumericVector threshold_gene,
             const double g_min, const double g_max,
             const double k_min, const double k_max,
             const int interaction_types,
             const long model_count_max,const long threshold_max,
             const double h, const double lambda_min,
             const double lambda_max, const int n_min, const int n_max,
             const double tot_time,
             const double sd_multiplier, const int number_gene,
             const double D_max,  const double D_shot_scaling,
             const int scaled_noise,
             const int D_levels, const double D_scaling,
             const int output_precision, const bool ANNEALING,
             const int initial_conditions, const String filename,
             const double print_start,
                 const double print_interval,
                 const bool integrate = true, const bool genParams = true,
                 const bool genIC = true, const int stepper = 1,
                 const double rk_tolerance = 0.001)

{

    std::string file_name= filename;
//    Rcout<<"Running time evolution simulations for "
//    <<std::to_string(number_gene)<<" genes..."<<"\n";
    //  if(initial_conditions>1)
    double D=D_max; // setting noise to maximum noise level
    std::vector<double> Darray(number_gene);
    //array to scale the noise level in each gene
    //Rcout<<"parameter_file"<<parameters_file<<"\n";
    // if(parameters_file) Rcout<<"If true"<<"\n";
    // Scale the noise level in each gene if scaled_noise = 1
    if(scaled_noise==1){

      for(int i=0; i<number_gene;i++)
      {
        Darray[i]=threshold_gene[i];}
 //     Rcout<<"Using a noise level that is proportional to median expression of the gene"<<"\n";

    }
    else{
      for(int i=0; i<number_gene;i++)
      {
        Darray[i]=1.0;
      }
//      Rcout<<"Using same noise level for each gene"<<"\n";

    }

    //Create output files if not there already
    std::fstream out_GE("./tmp/" + file_name +"_geneExpression.txt",
                        std::ios::out);
    if(!out_GE) {     Rcout << "Cannot open output file.\n";  return 1;}

    std::ifstream in_parameters;
    std::fstream out_param;

    if(genParams){
      out_param.open("./tmp/" + file_name + "_parameters.txt",std::ios::out);
    }
    else {
      in_parameters.open("./tmp/" + file_name + "_parameters.txt",
                         std::ifstream::in);
      if(!in_parameters) {     Rcout << "./tmp/" + file_name + "_parameters.txt"
        << "Cannot open input file for reading parameters.\n";  return 1;
      }
    }

    std::ifstream in_ic;
    std::fstream out_ic;
    if(genIC){
      out_ic.open("./tmp/" + file_name + "_IC.txt",std::ios::out);
    }
    else
    {
      in_ic.open("./tmp/" + file_name + "_IC.txt",std::ifstream::in);
      if(!in_ic) {
        Rcout <<"Cannot open input file for reading initial conditions.\n";
        return 1;
        }
    }

          for(long model_count=0;model_count<model_count_max;model_count++)
      {
            if((static_cast<long> (20*model_count) % model_count_max) == 0){
              Rcout<<"====";
            }

            // Check for user interrupt and exit if there is any.
        if (model_count % 100 == 0)
          Rcpp::checkUserInterrupt();

        //Initialize production rate of genes
        std::vector<double> g_gene(number_gene);

        //Initialize degradation rate of genes
        std::vector<double> k_gene(number_gene);

        //Initialize hill coefficient for each interaction
        std::vector<std::vector<int> >
          n_gene(number_gene, std::vector<int>(number_gene));

        //Initialize fold change for each interaction
        std::vector<std::vector<double> >
          lambda_gene(number_gene, std::vector<double>(number_gene));

        //Initialize threshold for each interaction
        std::vector<std::vector<double> >
          threshold_gene_log(number_gene, std::vector<double>(number_gene));

        if(in_parameters.is_open())
        {
          readParameters( gene_interaction, number_gene, g_gene,
                           k_gene, n_gene,
                           lambda_gene,
                           threshold_gene_log, in_parameters);
        }
        else
        {

          selectParameters( gene_interaction,  threshold_gene,
                             g_min,  g_max,
                             k_min,  k_max,  interaction_types,
                             model_count_max, threshold_max,
                             h,  lambda_min,
                             lambda_max,  n_min,  n_max,
                             tot_time,
                             sd_multiplier,  number_gene,
                             D_max,   D_shot_scaling,
                             scaled_noise,
                             D_levels,  D_scaling,
                             output_precision,  ANNEALING,
                             initial_conditions, filename, g_gene,
                             k_gene, n_gene,
                             lambda_gene,
                             threshold_gene_log, out_param);


          writeParameters( number_gene, output_precision, g_gene,
                           k_gene, n_gene,
                           lambda_gene,
                           threshold_gene_log, out_param);
        }


        /////////////////////////////////////////////////////////////////////

        //Initial condition range selection
        /////////////////////////////////////////////////////////////////////
        std::vector<double> max_gene(number_gene);
        std::vector<double> min_gene(number_gene);
        if(genIC)
        {
          selectIcRange( number_gene,  gene_interaction, g_gene, k_gene,
                         lambda_gene, max_gene, min_gene);
        }

        ///////////////////////////////////////////////////////////////////////

        //Initial condition  selection
        ///////////////////////////////////////////////////////////////////////
        for(int ic_count=0;ic_count<initial_conditions;ic_count++)
        {
          std::vector <double> expression_gene(number_gene);
          //array for current gene expression
          std::vector <double> expression_gene0(number_gene);
          //array for initial gene expression
          if(!genIC)
          {
            for(size_t ic_counter=0;ic_counter <number_gene; ic_counter++)
            {
              in_ic >> expression_gene0[ic_counter];
            }
          }

          else
          {
            for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
            {
              expression_gene0[gene_count1]=exp(log(min_gene[gene_count1]) +
                (log(max_gene[gene_count1]) -
                log(min_gene[gene_count1]))*u_distribution(u_generator));

            }

            ///////////////////////////////////////////////////////////////////

            //Writing initial condition to file
            //////////////////////////////////////////////////////////////////

            for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
            {
              out_ic<<std::setprecision(output_precision)
              <<expression_gene0[gene_count1]<<"\t";
            }
            out_ic<<"\n";

          }

          ///////////////////////////////////////////////////////////////////

          //Time Evolution
          ///////////////////////////////////////////////////////////////////
          D=D_max; //Start with maximum noise level for each model

          for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
          {
            expression_gene[gene_count1]=expression_gene0[gene_count1];
          }

          for(int file_count=0;file_count<D_levels;file_count++)
          {
            if(file_count==D_levels-1){D=0;}



            if(ANNEALING) {}
            else {
              for(int gene_count_temp=0;gene_count_temp<number_gene;
              gene_count_temp++) {
                expression_gene[gene_count_temp] =
                  expression_gene0[gene_count_temp];
              }
            }
            if(integrate) {
            switch(stepper){
            case 1:
              // Euler Maruyama method
              stepEM( expression_gene, out_GE, tot_time,
                    number_gene, gene_interaction, g_gene, k_gene, n_gene,
                    lambda_gene, threshold_gene_log, interaction_types,
                    sd_multiplier, D_shot_scaling, Darray,
                    output_precision, print_start, print_interval, D, h);
              break;
            case 4:
              //fourth order Runge-Kutta
//              Rcout<<"RK4";
              stepRK4( expression_gene, out_GE, tot_time, number_gene,
                      gene_interaction, g_gene, k_gene, n_gene, lambda_gene,
                       threshold_gene_log, interaction_types,
                       sd_multiplier,
                       output_precision,
                       print_start,  print_interval, h);
              break;

            case 5:
//              Rcout<<"DP";
              // adaptive Dormand Prince
              stepDP( expression_gene,out_GE,tot_time,number_gene,
                      gene_interaction,g_gene,k_gene,n_gene,lambda_gene,
                      threshold_gene_log,interaction_types,
                      sd_multiplier,
                      output_precision,print_start, print_interval,h,
                      rk_tolerance);
              break;

            default:
              Rcout<< "Error in specifying the stepper.\n";

            }
            }
            //Rcout<<"D="<<D<<"\n";
            // Rcout<< "Noise Level" << file_count<<"\t"<<D<<"\n";
            D=D*D_scaling;

          }

          out_GE<<"\n";
        }
      }


    out_GE.close();
    out_param.close();
    out_ic.close();
    if(in_ic.is_open()) in_ic.close();
    if(in_parameters.is_open()) in_parameters.close();

// Rcout<<"Simulations completed successfully.\n";

  return 0;
}

