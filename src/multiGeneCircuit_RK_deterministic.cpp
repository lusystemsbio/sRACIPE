#include"header.h"
using namespace Rcpp;


// extern unsigned u_seed = 3654734;//std::chrono::system_clock::now().time_since_epoch().count();
//extern unsigned g_seed = std::chrono::system_clock::now().time_since_epoch().count()*M_PI_4;

//extern std::mt19937_64 u_generator;// (u_seed);
// std::mt19937_64 g_generator (g_seed);
// Shifted hill function

//extern double Hs_Racipe(double A, double AB0, int n_ab, double lambda_ab)
//{
//  return lambda_ab+(1-lambda_ab)*1/(1+pow((A/AB0),n_ab));
//}


//uniformly distributed random number generator in (0,1) range
//extern std::uniform_real_distribution<double> u_distribution(0.0,1.0);

// Gaussian distributed random number generator with mean 0 and 1 standard deviation
 //std::normal_distribution<double> g_distribution(0.0,1.0);


// [[Rcpp::export]]

int multiGeneCircuit_RK_deterministic(IntegerMatrix gene_interaction, NumericVector threshold_gene,
                                      double g_min, double g_max,
                                      double k_min, double k_max, int possible_interactions,
                                      long model_count_max,long threshold_max,
                                      double h, double lambda_min,
                                      double lambda_max, int n_min, int n_max,
                                      double tot_time, double median_range,
                                      double standard_deviation_factor, int number_gene,
                                      int output_precision, int INITIAL_CONDITIONS, String filename)

{
  std::string file_name= filename;

    Rcout<<"Running time evolution simulations for "<<file_name<<" consisting of"<<std::to_string(number_gene)<<" genes for the deterministic case using Runge-Kutta integration method."<<"\n";


    std::fstream out_0("./results/sRACIPE_RK_"+file_name+"_g"+std::to_string(number_gene)+"_parameters.txt",std::ios::out);
    if(!out_0) {     Rcout << "Cannot open output file for writing parameters.\n";  }

    std::fstream out_ic("./results/sRACIPE_RK_"+file_name+"_g"+std::to_string(number_gene)+"_IC.txt",std::ios::out);
    if(!out_ic) {     Rcout << "Cannot open output file for writing parameters.\n";  }

    std::fstream out_1("./results/sRACIPE_RK_"+file_name+"_g"+std::to_string(number_gene)+"_output.txt",std::ios::out);
    if(!out_1) {     Rcout << "Cannot open output file.\n";  }




    for(long model_count=0;model_count<model_count_max;model_count++)
    {

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

            for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
            {out_ic<<std::setprecision(output_precision)<<expression_gene0[gene_count1]<<"\t";} //initial condition of each gene

            out_ic<<"\n";


            ///////////////////////////////////////////////////////////////////////////////////////

            //Time Evolution
            ///////////////////////////////////////////////////////////////////////////////////////


            double expression_gene_h[number_gene]; //array for temp gene expression values
            double expression_gene_h1[number_gene]; //array for temp gene expression values
            double expression_gene_h2[number_gene]; //array for temp gene expression values
            double expression_gene_h3[number_gene]; //array for temp gene expression values
            double expression_gene_h4[number_gene]; //array for temp gene expression values

            for(int gene_count_temp=0;gene_count_temp<number_gene;gene_count_temp++)
            {
                expression_gene[gene_count_temp]=expression_gene0[gene_count_temp];
                expression_gene_h[gene_count_temp]=expression_gene0[gene_count_temp];
                expression_gene_h1[gene_count_temp]=expression_gene0[gene_count_temp];
                expression_gene_h2[gene_count_temp]=expression_gene0[gene_count_temp];
                expression_gene_h3[gene_count_temp]=expression_gene0[gene_count_temp];
                expression_gene_h4[gene_count_temp]=expression_gene0[gene_count_temp];
            }

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


                    expression_gene_h1[gene_count1]=h*(g_gene[gene_count1]*final_multiplier-k_gene[gene_count1]*expression_gene[gene_count1]);
                }

                for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
                {
                    double final_multiplier=1;

                    for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
                    {
                        double gene_value=expression_gene[gene_count2]+0.5*expression_gene_h1[gene_count2];
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


                    expression_gene_h2[gene_count1]=h*((g_gene[gene_count1])*final_multiplier-k_gene[gene_count1]*(expression_gene[gene_count1]+0.5*expression_gene_h1[gene_count1]));
                }

                for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
                {
                    double final_multiplier=1;

                    for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
                    {
                        double gene_value=expression_gene[gene_count2]+0.5*expression_gene_h2[gene_count2];
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


                    expression_gene_h3[gene_count1]=h*((g_gene[gene_count1])*final_multiplier-k_gene[gene_count1]*(expression_gene[gene_count1]+0.5*expression_gene_h2[gene_count1]));
                }


                for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
                {
                    double final_multiplier=1;

                    for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
                    {
                        double gene_value=expression_gene[gene_count2]+expression_gene_h3[gene_count2];
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


                    expression_gene_h4[gene_count1]=h*((g_gene[gene_count1])*final_multiplier-k_gene[gene_count1]*(expression_gene[gene_count1]+expression_gene_h3[gene_count1]));
                }


                /////////////////////////////////////////////////////////////
                for(int gene_count1=0;gene_count1<number_gene;gene_count1++){
                    expression_gene_h[gene_count1]=expression_gene[gene_count1]+(expression_gene_h1[gene_count1]+2*expression_gene_h2[gene_count1]+2*expression_gene_h3[gene_count1]+expression_gene_h4[gene_count1])/6;

                    if(expression_gene_h[gene_count1]<0) expression_gene_h[gene_count1]=0;
                    expression_gene[gene_count1]=expression_gene_h[gene_count1];}



                i+=h;


            }while(i<tot_time);


            for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
            {
                out_1<<std::setprecision(output_precision)<<expression_gene[gene_count1]<<"\t";
            }
            out_1<<"\n";

        }

    }
    out_ic<<"\n";

    out_1.close();
    out_0.close();
    Rcout<<"Deterministic simulations using Runge-Kutta completed successfully. Data files are in results folder.\n";
    return 0;
}


