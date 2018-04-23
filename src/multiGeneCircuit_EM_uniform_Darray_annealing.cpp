#include "header.h"
#include <Rcpp.h>


using namespace Rcpp;

extern unsigned u_seed;//2 = std::chrono::system_clock::now().time_since_epoch().count();
unsigned g_seed = std::chrono::system_clock::now().time_since_epoch().count()*M_PI_4;

extern std::mt19937_64 u_generator;// (u_seed2);
std::mt19937_64 g_generator (g_seed);

//uniformly distributed random number generator2 in (0,1) range
extern std::uniform_real_distribution<double> u_distribution;

// Gaussian distributed random number generator2 with mean 0 and 1 standard deviation
std::normal_distribution<double> g_distribution(0.0,1.0);


// Shifted hill function
extern double Hs_Racipe(double A, double AB0, int n_ab, double lambda_ab);




// [[Rcpp::export]]

int multiGeneCircuit_EM_uniform_Darray_annealing(IntegerMatrix gene_interaction, NumericVector threshold_gene,
                                                 double g_min, double g_max,
                                                 double k_min, double k_max, int possible_interactions,
                                                 long model_count_max,long threshold_max,
                                                 double h, double lambda_min,
                                                 double lambda_max, int n_min, int n_max,
                                                 double tot_time, double median_range,
                                                 double standard_deviation_factor, int number_gene,
                                                 double D_max,  double D_shot_scaling,
                                                 int GENE_NOISE_SCALING, int file_writing_interval,
                                                 int D_levels, double D_scaling,
                                                 int output_precision, int ANNEALING, int CONSTANT_NOISE, int INITIAL_CONDITIONS, String filename)

{
  std::string file_name= filename;
  Rcout<<"Running time evolution simulations for "<<std::to_string(number_gene)<<" genes..."<<"\n";
  if(INITIAL_CONDITIONS>1) file_writing_interval=1;
  D_max=D_max/sqrt(number_gene); //Scale the maximum noise by the number of genes.
  double D=D_max; // setting noise to maximum noise level
  double Darray[number_gene]; //array to scale the noise level in each gene

  // Scale the noise level in each gene if GENE_NOISE_SCALING=1
  if(GENE_NOISE_SCALING==1){

    for(int i=0; i<number_gene;i++)
    {
      Darray[i]=threshold_gene[i];}
    Rcout<<"Using a noise level that is proportional to median expression of the gene"<<"\n";

  }
  else{
    for(int i=0; i<number_gene;i++)
    {
      Darray[i]=1.0;
    }
    Rcout<<"Using same noise level for each gene"<<"\n";

  }

  //Create output file if not there already


  std::fstream out_0("./results/sRACIPE_EM_"+file_name+"_g"+std::to_string(number_gene)+"_Annealing_"+std::to_string(ANNEALING)+"_parameters.txt",std::ios::out);
  if(!out_0) {     Rcout << "Cannot open output file for writing parameters.\n";  }

  std::fstream out_ic("./results/sRACIPE_EM_"+file_name+"_g"+std::to_string(number_gene)+"_Annealing_"+std::to_string(ANNEALING)+"_IC.txt",std::ios::out);
  if(!out_ic) {     Rcout << "Cannot open output file for writing parameters.\n";  }

  std::fstream out_1("./results/sRACIPE_EM_"+file_name+"_g"+std::to_string(number_gene)+"_Annealing_"+std::to_string(ANNEALING)+"_output.txt",std::ios::out);
  if(!out_1) {     Rcout << "Cannot open output file.\n";  }

  // std::fstream out_D("./results/multiGeneCircuit_EM_g"+std::to_string(number_gene)+"Annealing_"+std::to_string(ANNEALING)+"_noise.txt",std::ios::out);
  //  if(!out_D) {     Rcout << "Cannot open output file.\n";  }


  //Write noise levels at the top of file. This part moved to R.
  // for(int file_count=0;file_count<D_levels;file_count++)
  // {
  //   //for(int gene_count2=0; gene_count2<number_gene; gene_count2++)
  //   {
  //     //out_0<<"supercount="<<super_count<<endl;
  //     if(file_count==D_levels-1) {out_0<<std::setprecision(output_precision)<<0<<"\t";}
  //     else
  //     {out_D<<std::setprecision(output_precision)<<D_max*pow(D_scaling,file_count)<<"\t";}
  //   }
  // }
  // out_D<<"\n";
  // out_D.close();


  //File writing is minimized to speed up the program. Output is written to file after every FILE_WRITING_INTERVAL times.
  int super_count_max=int(model_count_max/file_writing_interval);
  Rcout<<"Super Count="<<super_count_max<<"\n";
  model_count_max=file_writing_interval;
  for(int super_count=0;super_count<super_count_max;super_count++)
  {

    //A temporary array to store gene expression values
    double expression_gene_final[file_writing_interval][D_levels][number_gene];
    for(int fwi_count=0;fwi_count<file_writing_interval;fwi_count++)for(int fn_count=0;fn_count<D_levels;fn_count++)for(int gf=0;gf<number_gene;gf++)expression_gene_final[fwi_count][fn_count][gf]=0;

    for(long model_count=0;model_count<model_count_max;model_count++)
    {
      D=D_max; //Start with maximum noise level for each model

      //Initialize production rate of genes
      double g_gene[number_gene];
      for(int gene_count1=0;gene_count1<number_gene;gene_count1++){g_gene[gene_count1]=g_min+(g_max-g_min)*u_distribution(u_generator);}

      //Initialize degradation rate of genes
      double k_gene[number_gene];//
      for(int gene_count1=0;gene_count1<number_gene;gene_count1++){k_gene[gene_count1]=k_min+(k_max-k_min)*u_distribution(u_generator);}

      //Initialize hill coefficient for each interaction

      int n_gene[number_gene][number_gene];
      for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
      {
        for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
        {
          if(gene_interaction(gene_count1,gene_count2)==0)
          {
            n_gene[gene_count1][gene_count2]=0;
          }
          else n_gene[gene_count1][gene_count2]=int((n_max-n_min)*u_distribution(u_generator))+n_min;
        }
      }

      //Initialize fold change for each interaction

      double lambda_gene[number_gene][number_gene];
      for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
      {
        for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
        {
          if(gene_interaction(gene_count1,gene_count2)==0)
          {
            lambda_gene[gene_count1][gene_count2]=0;
          }
          else lambda_gene[gene_count1][gene_count2]=(lambda_max-lambda_min)*u_distribution(u_generator)+lambda_min;
        }
      }

      //Initialize threshold for each interaction

      double threshold_gene_log[number_gene][number_gene];
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

            threshold_gene_log[gene_count1][gene_count2]=(1-standard_deviation_factor*sqrt(3))*threshold_gene[gene_count2]+(2*sqrt(3)*standard_deviation_factor*threshold_gene[gene_count2])*u_distribution(u_generator);
          }
        }
      }
      ///////////////////////////////////////////////////////////////////////////////////////

      //Writing parameters to file
      ///////////////////////////////////////////////////////////////////////////////////////
      for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
      {out_0<<std::setprecision(output_precision)<<g_gene[gene_count1]<<"\t";} //production rate of each gene

      for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
      {out_0<<std::setprecision(output_precision)<<k_gene[gene_count1]<<"\t";} // degradation rate of each gene

      for(int gene_count1=0;gene_count1<number_gene;gene_count1++)for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
      {if(threshold_gene_log[gene_count1][gene_count2]>0)
        out_0<<std::setprecision(output_precision)<<threshold_gene_log[gene_count1][gene_count2]<<"\t";}
      // above--thresholds for the inteaction links, thresholds for inward links for genes are written, starting from gene 1 (thresholds for all inward links of gene 1 and so on)
      for(int gene_count1=0;gene_count1<number_gene;gene_count1++)for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
      {if(n_gene[gene_count1][gene_count2]>0)
        out_0<<std::setprecision(1)<<n_gene[gene_count1][gene_count2]<<"\t";}
      // above--n for the inteaction links, n for inward links for genes are written, starting from gene 1
      for(int gene_count1=0;gene_count1<number_gene;gene_count1++)for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
      {if(lambda_gene[gene_count1][gene_count2]>0)
        out_0<<std::setprecision(output_precision)<<lambda_gene[gene_count1][gene_count2]<<"\t";}
      // above--lambda for the inteaction links, lambda for inward links for genes are written, starting from gene 1
      out_0<<"\n";
      ///////////////////////////////////////////////////////////////////////////////////////

      //Initial condition range selection
      ///////////////////////////////////////////////////////////////////////////////////////
      double max_gene[number_gene];for(int gene_count1=0;gene_count1<number_gene;gene_count1++){max_gene[gene_count1]=g_gene[gene_count1]/k_gene[gene_count1];}

      double min_gene[number_gene]; //=lambda_ba/lambda_aa;

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
            Rcout << "Invalid Interation code for Gene"<<gene_count1<<" and gene"<<gene_count2<<" interaction"<<"\n";
          }

          min_gene_multiplier_final=min_gene_multiplier_final*gene_min_multiplier;
        }

        min_gene[gene_count1]=g_gene[gene_count1]/k_gene[gene_count1]*min_gene_multiplier_final;

      }


      ///////////////////////////////////////////////////////////////////////////////////////

      //Initial condition  selection
      ///////////////////////////////////////////////////////////////////////////////////////
      for(int ic_count=0;ic_count<INITIAL_CONDITIONS;ic_count++)
      {
        double expression_gene[number_gene]; //array for current gene expression
        double expression_gene0[number_gene]; //array for initial gene expression

        for(int gene_count1=0;gene_count1<number_gene;gene_count1++){expression_gene0[gene_count1]=exp(log(min_gene[gene_count1])+(log(max_gene[gene_count1])-log(min_gene[gene_count1]))*u_distribution(u_generator));
          expression_gene[gene_count1]=expression_gene0[gene_count1];}


        ///////////////////////////////////////////////////////////////////////////////////////

        //Writing initial condition to file
        ///////////////////////////////////////////////////////////////////////////////////////


        //Rcout<<"Written"<<"\n";

       // for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
        //{out_ic<<std::setprecision(output_precision)<<expression_gene0[gene_count1]<<"\t";} //initial condition of each gene




        ///////////////////////////////////////////////////////////////////////////////////////

        //Time Evolution
        ///////////////////////////////////////////////////////////////////////////////////////

        for(int file_count=0;file_count<D_levels;file_count++)
        {

          double expression_gene_h[number_gene]; //array for temp gene expression values
          if(ANNEALING==1) {}
          else {
            for(int gene_count_temp=0;gene_count_temp<number_gene;gene_count_temp++)
            { expression_gene[gene_count_temp]=expression_gene0[gene_count_temp];

              expression_gene_h[gene_count_temp]=expression_gene0[gene_count_temp];
             // out_ic<<std::setprecision(output_precision)<<expression_gene0[gene_count_temp]<<"\t";
              }
          }
          for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
            {out_ic<<std::setprecision(output_precision)<<expression_gene0[gene_count1]<<"\t";} //initial condition of each gene
          //Rcout<<"D="<<D<<"\n";
          double i=0.0;
          do
          {
            for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
            {
              double final_multiplier=1;

              for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
              {
                double gene_value=expression_gene[gene_count2];
                double gene_threshold=threshold_gene_log[gene_count1][gene_count2];
                int gene_n=n_gene[gene_count1][gene_count2];
                double gene_lambda=lambda_gene[gene_count1][gene_count2];
                double gene_activation_multiplier=1;

                switch(gene_interaction(gene_count1,gene_count2))
                {
                case 0:
                  gene_activation_multiplier=1.0;
                  break;

                case 2:
                  gene_lambda=1./gene_lambda;
                  gene_activation_multiplier=gene_lambda+(1.-gene_lambda)*1./(1.+pow((gene_value/gene_threshold),gene_n));
                  break;

                case 1:
                  gene_activation_multiplier=(gene_lambda+(1.-gene_lambda)*1./(1.+pow((gene_value/gene_threshold),gene_n)))/gene_lambda;
                  break;

                default :
                  Rcout << "Invalid Interation code for Gene"<<gene_count1<<" and gene"<<gene_count2<<" interaction"<<"\n";
                }

                final_multiplier=final_multiplier*gene_activation_multiplier;
              }


              expression_gene_h[gene_count1]=expression_gene[gene_count1]+h*(g_gene[gene_count1]*final_multiplier-k_gene[gene_count1]*expression_gene[gene_count1])+
                D*sqrt(h)*g_distribution(g_generator)*Darray[gene_count1]+
                D_shot_scaling*D*sqrt(h)*g_distribution(g_generator)*Darray[gene_count1]*expression_gene[gene_count1];
              if(expression_gene_h[gene_count1]<0) expression_gene_h[gene_count1]=0;
            }

            for(int gene_count1=0;gene_count1<number_gene;gene_count1++){
              expression_gene[gene_count1]=expression_gene_h[gene_count1];}



            i+=h;


          }while(i<tot_time);
          if(file_count==D_levels-2){D=0;}
          else {D=D*D_scaling;}
          if(file_writing_interval==1){
            //Rcout<<"IC";
            for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
            {
              out_1<<std::setprecision(output_precision)<<expression_gene[gene_count1]<<"\t";
            }
            out_1<<"\n";
          }
          else {

            for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
            {
              expression_gene_final[model_count][file_count][gene_count1]=expression_gene[gene_count1];
            }
          }

        }
        out_ic<<"\n";

      }
    }
    ///////////////////////////////////////////////////////////////////////////////////////

    //Writing to file
    ///////////////////////////////////////////////////////////////////////////////////////

    if(file_writing_interval>1){
      for(int fwi_count=0;fwi_count<file_writing_interval;fwi_count++)
      {
        for(int file_count=0;file_count<D_levels;file_count++)
        {

          for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
          {
            out_1<<std::setprecision(output_precision)<<expression_gene_final[fwi_count][file_count][gene_count1]<<"\t";

          }
        }
        out_1<<"\n";

      }
    }


  }
  out_1.close();
  out_0.close();
  Rcout<<"Simulations completed successfully. Data files are in results folder.\n";
  return 0;
}


