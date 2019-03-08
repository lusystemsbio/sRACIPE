#include"header.h"
#include <Rcpp.h>
using namespace Rcpp;

void stepEM( std::vector <double> &expression_gene,
             std::fstream &out_GE,
             const double &tot_time,
             const int &number_gene,
             IntegerMatrix gene_interaction,
             const std::vector<double> &g_gene,
             const std::vector<double> &k_gene,
             const std::vector<std::vector<int> > &n_gene,
             const std::vector<std::vector<double> > &lambda_gene,
             const std::vector<std::vector<double> > &threshold_gene_log,
             const int &possible_interactions,
             const double &standard_deviation_factor,
             const double &D_shot_scaling,
             const std::vector<double> &Darray,
             const int &output_precision,
             const double &print_start, const double &print_interval,
             const double &D,
             const double &h){

  double expression_gene_h[number_gene]; //array for temp gene expression values

  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
  {
    expression_gene_h[gene_count1] = expression_gene[gene_count1];
  }

  double i=0.0;
  int print_counter = 0;
  do
  {
    i+=h;
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
          gene_activation_multiplier = gene_lambda+(1.-gene_lambda)*
            1./(1.+pow((gene_value/gene_threshold),gene_n));
          break;

        case 1:
          gene_activation_multiplier=(gene_lambda+(1.-gene_lambda)*
            1./(1.+pow((gene_value/gene_threshold),gene_n)))/gene_lambda;
          break;

        default :
          Rcout << "Invalid Interation code for Gene"<<gene_count1
          <<" and gene"<<gene_count2<<" interaction"<<"\n";
        }

        final_multiplier=final_multiplier*gene_activation_multiplier;
      }
      expression_gene_h[gene_count1] = expression_gene[gene_count1] +
        h*(g_gene[gene_count1]*final_multiplier-k_gene[gene_count1]*
        expression_gene[gene_count1]) +
        D*sqrt(h)*g_distribution(g_generator)*Darray[gene_count1]+
        D_shot_scaling*D*sqrt(h)*g_distribution(g_generator)*
        Darray[gene_count1]*expression_gene[gene_count1];
      if(expression_gene_h[gene_count1]<0) expression_gene_h[gene_count1]=0;
    }

    for(int gene_count1=0;gene_count1<number_gene;gene_count1++){
      expression_gene[gene_count1]=expression_gene_h[gene_count1];}

    if((i> (print_start + print_interval*print_counter)) &&
       i <= (h+print_start + print_interval*print_counter))
    {
      print_counter++;
      //std::cout<<i<<"\n";
      for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
      {
        out_GE<<std::setprecision(output_precision)
        <<expression_gene[gene_count1]<<"\t";
      }
      //out_GE<<"\n";
    }
  }while(i<tot_time);

  // for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
  // {
  //   out_GE<<std::setprecision(output_precision)<<expression_gene[gene_count1]<<"\t";
  // }

}

///////////////////////////////////////////////////////////////////////////////////////

//Runge Kutta fourth order
///////////////////////////////////////////////////////////////////////////////////////
void stepRK4( std::vector <double> &expression_gene,
             std::fstream &out_GE,
             const double &tot_time,
             const int &number_gene,
             IntegerMatrix gene_interaction,
             const std::vector<double> &g_gene,
             const std::vector<double> &k_gene,
             const std::vector<std::vector<int> > &n_gene,
             const std::vector<std::vector<double> > &lambda_gene,
             const std::vector<std::vector<double> > &threshold_gene_log,
             const int &possible_interactions,
             const double &standard_deviation_factor,
             const int &output_precision,
             const double &print_start, const double &print_interval,
             const double &h){

double expression_gene_h1[number_gene]; //array for temp gene expression values
double expression_gene_h2[number_gene]; //array for temp gene expression values
double expression_gene_h3[number_gene]; //array for temp gene expression values
double expression_gene_h4[number_gene]; //array for temp gene expression values

for(int gene_count_temp=0;gene_count_temp<number_gene;gene_count_temp++)
{
  expression_gene_h1[gene_count_temp]=expression_gene[gene_count_temp];
  expression_gene_h2[gene_count_temp]=expression_gene[gene_count_temp];
  expression_gene_h3[gene_count_temp]=expression_gene[gene_count_temp];
  expression_gene_h4[gene_count_temp]=expression_gene[gene_count_temp];
}
double i=0.0;
int print_counter = 0;
do
{
  i+=h;
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
        gene_activation_multiplier=gene_lambda+(1.-gene_lambda)*
          1./(1.+pow((gene_value/gene_threshold),gene_n));
        break;

      case 1:
        gene_activation_multiplier=(gene_lambda+(1.-gene_lambda)*
          1./(1.+pow((gene_value/gene_threshold),gene_n)))/gene_lambda;
        break;

      default :
        Rcout << "Invalid Interation code for Gene"<<gene_count1
        <<" and gene"<<gene_count2<<" interaction"<<"\n";
      }

      final_multiplier=final_multiplier*gene_activation_multiplier;
    }


    expression_gene_h1[gene_count1]=h*(g_gene[gene_count1]*final_multiplier -
      k_gene[gene_count1]*expression_gene[gene_count1]);
  }

  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
  {
    double final_multiplier=1;

    for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
    {
      double gene_value=expression_gene[gene_count2] +
        0.5*expression_gene_h1[gene_count2];
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
        gene_activation_multiplier=gene_lambda+(1.-gene_lambda)*
          1./(1.+pow((gene_value/gene_threshold),gene_n));
        break;

      case 1:
        gene_activation_multiplier=(gene_lambda+(1.-gene_lambda)*
          1./(1.+pow((gene_value/gene_threshold),gene_n)))/gene_lambda;
        break;

      default :
        Rcout << "Invalid Interation code for Gene"<<gene_count1
        <<" and gene"<<gene_count2<<" interaction"<<"\n";
      }

      final_multiplier=final_multiplier*gene_activation_multiplier;
    }


    expression_gene_h2[gene_count1]=h*((g_gene[gene_count1])*
      final_multiplier-k_gene[gene_count1]*(expression_gene[gene_count1] +
      0.5*expression_gene_h1[gene_count1]));
  }

  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
  {
    double final_multiplier=1;

    for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
    {
      double gene_value=expression_gene[gene_count2] +
        0.5*expression_gene_h2[gene_count2];
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
        gene_activation_multiplier=gene_lambda+(1.-gene_lambda)*
          1./(1.+pow((gene_value/gene_threshold),gene_n));
        break;

      case 1:
        gene_activation_multiplier=(gene_lambda+(1.-gene_lambda)*
          1./(1.+pow((gene_value/gene_threshold),gene_n)))/gene_lambda;
        break;

      default :
        Rcout << "Invalid Interation code for Gene"<<gene_count1
        <<" and gene"<<gene_count2<<" interaction"<<"\n";
      }

      final_multiplier=final_multiplier*gene_activation_multiplier;
    }


    expression_gene_h3[gene_count1]=h*((g_gene[gene_count1])*final_multiplier-
      k_gene[gene_count1]*(expression_gene[gene_count1] +
      0.5*expression_gene_h2[gene_count1]));
  }


  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
  {
    double final_multiplier=1;

    for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
    {
      double gene_value=expression_gene[gene_count2] +
        expression_gene_h3[gene_count2];
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
        gene_activation_multiplier=gene_lambda+(1.-gene_lambda)*
          1./(1.+pow((gene_value/gene_threshold),gene_n));
        break;

      case 1:
        gene_activation_multiplier=(gene_lambda+(1.-gene_lambda)*
          1./(1.+pow((gene_value/gene_threshold),gene_n)))/gene_lambda;
        break;

      default :
        Rcout << "Invalid Interation code for Gene"<<gene_count1
        <<" and gene"<<gene_count2<<" interaction"<<"\n";
      }

      final_multiplier=final_multiplier*gene_activation_multiplier;
    }


    expression_gene_h4[gene_count1]=h*((g_gene[gene_count1])*final_multiplier-
      k_gene[gene_count1]*(expression_gene[gene_count1]+
      expression_gene_h3[gene_count1]));
  }


  /////////////////////////////////////////////////////////////
  for(int gene_count1=0;gene_count1<number_gene;gene_count1++){
    expression_gene[gene_count1] = expression_gene[gene_count1]+
      (expression_gene_h1[gene_count1]+2*expression_gene_h2[gene_count1]+
      2*expression_gene_h3[gene_count1]+expression_gene_h4[gene_count1])/6;
    if(expression_gene[gene_count1]<0) expression_gene[gene_count1]=0;
    }


  if((i> (print_start + print_interval*print_counter)) &&
     i <= (h+print_start + print_interval*print_counter))
  {
    print_counter++;
    //std::cout<<i<<"\n";
    for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
    {
      out_GE<<std::setprecision(output_precision)
            <<expression_gene[gene_count1]<<"\t";
    }
    //out_GE<<"\n";
  }
}while(i<tot_time);
}



void stepDP( std::vector <double> &expression_gene,
              std::fstream &out_GE,
              const double &tot_time,
              const int &number_gene,
              IntegerMatrix gene_interaction,
              const std::vector<double> &g_gene,
              const std::vector<double> &k_gene,
              const std::vector<std::vector<int> > &n_gene,
              const std::vector<std::vector<double> > &lambda_gene,
              const std::vector<std::vector<double> > &threshold_gene_log,
              const int &possible_interactions,
              const double &standard_deviation_factor,
              const int &output_precision,
              const double &print_start, const double &print_interval,
              double h, const double &rk_tolerance){
  double expression_gene_h[number_gene]; //array for temp gene expression values
  double expression_gene_h1[number_gene]; //array for temp gene expression values
  double expression_gene_h2[number_gene]; //array for temp gene expression values
  double expression_gene_h3[number_gene]; //array for temp gene expression values
  double expression_gene_h4[number_gene]; //array for temp gene expression values
  double expression_gene_h5[number_gene]; //array for temp gene expression values
  double expression_gene_h6[number_gene]; //array for temp gene expression values
  double expression_gene_h7[number_gene]; //array for temp gene expression values

  for(int gene_count_temp=0;gene_count_temp<number_gene;gene_count_temp++)
  {
    expression_gene_h[gene_count_temp]=expression_gene[gene_count_temp];
    expression_gene_h1[gene_count_temp]=expression_gene[gene_count_temp];
    expression_gene_h2[gene_count_temp]=expression_gene[gene_count_temp];
    expression_gene_h3[gene_count_temp]=expression_gene[gene_count_temp];
    expression_gene_h4[gene_count_temp]=expression_gene[gene_count_temp];
    expression_gene_h5[gene_count_temp]=expression_gene[gene_count_temp];
    expression_gene_h6[gene_count_temp]=expression_gene[gene_count_temp];
    expression_gene_h7[gene_count_temp]=expression_gene[gene_count_temp];
  }

  double i=0.0;
  int print_counter = 0;

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
          gene_activation_multiplier=gene_lambda+(1.-gene_lambda)*
            1./(1.+pow((gene_value/gene_threshold),gene_n));
          break;

        case 1:
          gene_activation_multiplier=(gene_lambda+(1.-gene_lambda)*
            1./(1.+pow((gene_value/gene_threshold),gene_n)))/gene_lambda;
          break;

        default :
          Rcout << "Invalid Interation code for Gene"<<gene_count1
          <<" and gene"<<gene_count2<<" interaction"<<"\n";
        }

        final_multiplier=final_multiplier*gene_activation_multiplier;
      }


      expression_gene_h1[gene_count1]=h*(g_gene[gene_count1]*final_multiplier-
        k_gene[gene_count1]*expression_gene[gene_count1]);
    }

    for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
    {
      double final_multiplier=1;

      for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
      {
        double gene_value=expression_gene[gene_count2]+
          0.2*expression_gene_h1[gene_count2];
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
          gene_activation_multiplier=gene_lambda+(1.-gene_lambda)*
            1./(1.+pow((gene_value/gene_threshold),gene_n));
          break;

        case 1:
          gene_activation_multiplier=(gene_lambda+(1.-gene_lambda)*
            1./(1.+pow((gene_value/gene_threshold),gene_n)))/gene_lambda;
          break;

        default :
          Rcout << "Invalid Interation code for Gene"<<gene_count1
          <<" and gene"<<gene_count2<<" interaction"<<"\n";
        }

        final_multiplier=final_multiplier*gene_activation_multiplier;
      }


      expression_gene_h2[gene_count1]=h*((g_gene[gene_count1])*final_multiplier-
        k_gene[gene_count1]*(expression_gene[gene_count1] +
        0.2*expression_gene_h1[gene_count1]));
    }

    for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
    {
      double final_multiplier=1;

      for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
      {
        double gene_value=expression_gene[gene_count2]+
          (0.25*expression_gene_h1[gene_count2]+
          0.75*expression_gene_h2[gene_count2])*0.3;
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
          gene_activation_multiplier=gene_lambda+(1.-gene_lambda)*
            1./(1.+pow((gene_value/gene_threshold),gene_n));
          break;

        case 1:
          gene_activation_multiplier=(gene_lambda+(1.-gene_lambda)*
            1./(1.+pow((gene_value/gene_threshold),gene_n)))/gene_lambda;
          break;

        default :
          Rcout << "Invalid Interation code for Gene"<<gene_count1
          <<" and gene"<<gene_count2<<" interaction"<<"\n";
        }

        final_multiplier=final_multiplier*gene_activation_multiplier;
      }


      expression_gene_h3[gene_count1]=h*((g_gene[gene_count1])*final_multiplier-
        k_gene[gene_count1]*(expression_gene[gene_count1]+
        (0.25*expression_gene_h1[gene_count1]+
        0.75*expression_gene_h2[gene_count1])*0.3));
    }


    for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
    {
      double final_multiplier=1;

      for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
      {
        double gene_value=expression_gene[gene_count2]+
          0.8*(((11./9.)*expression_gene_h1[gene_count2]+
          (-14./3.)*expression_gene_h2[gene_count2])+
          (40./9.)*expression_gene_h3[gene_count2]);
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
          gene_activation_multiplier=gene_lambda+(1.-gene_lambda)*
            1./(1.+pow((gene_value/gene_threshold),gene_n));
          break;

        case 1:
          gene_activation_multiplier=(gene_lambda+(1.-gene_lambda)*
            1./(1.+pow((gene_value/gene_threshold),gene_n)))/gene_lambda;
          break;

        default :
          Rcout << "Invalid Interation code for Gene"<<gene_count1
          <<" and gene"<<gene_count2<<" interaction"<<"\n";
        }

        final_multiplier=final_multiplier*gene_activation_multiplier;
      }


      expression_gene_h4[gene_count1]=h*((g_gene[gene_count1])*final_multiplier-
        k_gene[gene_count1]*(expression_gene[gene_count1]+
        0.8*(((11./9.)*expression_gene_h1[gene_count1]+
        (-14./3.)*expression_gene_h2[gene_count1])+
        (40./9.)*expression_gene_h3[gene_count1])));
    }

    for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
    {
      double final_multiplier=1;

      for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
      {
        double gene_value=expression_gene[gene_count2]+(8./9.)*
          ((4843./1458.)*expression_gene_h1[gene_count2]+
          (-3170./243.)*expression_gene_h2[gene_count2]+
          (8056./729.)*expression_gene_h3[gene_count2]+
          (-53./162.)*expression_gene_h4[gene_count2]);

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
          gene_activation_multiplier=gene_lambda+(1.-gene_lambda)*
            1./(1.+pow((gene_value/gene_threshold),gene_n));
          break;

        case 1:
          gene_activation_multiplier=(gene_lambda+(1.-gene_lambda)*
            1./(1.+pow((gene_value/gene_threshold),gene_n)))/gene_lambda;
          break;

        default :
          Rcout << "Invalid Interation code for Gene"<<gene_count1
          <<" and gene"<<gene_count2<<" interaction"<<"\n";
        }

        final_multiplier=final_multiplier*gene_activation_multiplier;
      }


      expression_gene_h5[gene_count1]=h*((g_gene[gene_count1])*final_multiplier-
        k_gene[gene_count1]*(expression_gene[gene_count1]+(8./9.)*
        ((4843./1458.)*expression_gene_h1[gene_count1]+
        (-3170./243.)*expression_gene_h2[gene_count1]+
        (8056./729.)*expression_gene_h3[gene_count1]+
        (-53./162.)*expression_gene_h4[gene_count1])));
    }


    for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
    {
      double final_multiplier=1;

      for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
      {
        double gene_value=expression_gene[gene_count2]+
          ((9017./3168.)*expression_gene_h1[gene_count2]+
          (-355./33.)*expression_gene_h2[gene_count2]+
          (46732./5247.)*expression_gene_h3[gene_count2]+
          (49./176.)*expression_gene_h4[gene_count2]+
          (-5103./18656.)*expression_gene_h5[gene_count2]);
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
          gene_activation_multiplier=gene_lambda+(1.-gene_lambda)*
            1./(1.+pow((gene_value/gene_threshold),gene_n));
          break;

        case 1:
          gene_activation_multiplier=(gene_lambda+(1.-gene_lambda)*
            1./(1.+pow((gene_value/gene_threshold),gene_n)))/gene_lambda;
          break;

        default :
          Rcout << "Invalid Interation code for Gene"<<gene_count1
          <<" and gene"<<gene_count2<<" interaction"<<"\n";
        }

        final_multiplier=final_multiplier*gene_activation_multiplier;
      }


      expression_gene_h6[gene_count1]=h*((g_gene[gene_count1])*final_multiplier-
        k_gene[gene_count1]*(expression_gene[gene_count1]+
        (9017./3168.)*expression_gene_h1[gene_count1]+
        (-355./33.)*expression_gene_h2[gene_count1]+
        (46732./5247.)*expression_gene_h3[gene_count1]+
        (49./176.)*expression_gene_h4[gene_count1]+
        (-5103./18656.)*expression_gene_h5[gene_count1]));
    }



    for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
    {
      double final_multiplier=1;

      for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
      {
        double gene_value=expression_gene[gene_count2]+
          ((35./384.)*expression_gene_h1[gene_count2]+
          (500./113.)*expression_gene_h3[gene_count2]+
          (125./192.)*expression_gene_h4[gene_count2]+
          (-2187./6784.)*expression_gene_h5[gene_count2]+
          (11./84.)*expression_gene_h6[gene_count2]);
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
          gene_activation_multiplier=gene_lambda+(1.-gene_lambda)*
            1./(1.+pow((gene_value/gene_threshold),gene_n));
          break;

        case 1:
          gene_activation_multiplier=(gene_lambda+(1.-gene_lambda)*
            1./(1.+pow((gene_value/gene_threshold),gene_n)))/gene_lambda;
          break;

        default :
          Rcout << "Invalid Interation code for Gene"<<gene_count1
          <<" and gene"<<gene_count2<<" interaction"<<"\n";
        }

        final_multiplier=final_multiplier*gene_activation_multiplier;
      }


      expression_gene_h7[gene_count1]=h*((g_gene[gene_count1])*final_multiplier-
        k_gene[gene_count1]*(expression_gene[gene_count1]+
        (35./384.)*expression_gene_h1[gene_count1]+
        (500./113.)*expression_gene_h3[gene_count1]+
        (125./192.)*expression_gene_h4[gene_count1]+
        (-2187./6784.)*expression_gene_h5[gene_count1]+
        (11./84.)*expression_gene_h6[gene_count1]));
    }
    double max_diff_o4_o5=0;
    /////////////////////////////////////////////////////////////
    for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
    {
      expression_gene_h[gene_count1]=expression_gene[gene_count1]+
        (5179./57600.)*expression_gene_h1[gene_count1]+
        (7571./16695.)*expression_gene_h3[gene_count1]+
        (393./640.)*expression_gene_h4[gene_count1]+
        (-92097./339200.)*expression_gene_h5[gene_count1]+
        (187./2100.)*expression_gene_h6[gene_count1]+
        (1./40.)*expression_gene_h7[gene_count1];
      if(expression_gene_h[gene_count1]<0) expression_gene_h[gene_count1]=0;

      expression_gene[gene_count1]=expression_gene[gene_count1]+
        (35./384.)*expression_gene_h1[gene_count1]+
        (500./1113.)*expression_gene_h3[gene_count1]+
        (125./192.)*expression_gene_h4[gene_count1]+
        (-2187./6784.)*expression_gene_h5[gene_count1]+
        (11./84.)*expression_gene_h6[gene_count1];

      if(expression_gene[gene_count1]<0) expression_gene[gene_count1]=0;

      double diff_o4_o5=expression_gene[gene_count1] -
        expression_gene_h[gene_count1];

      diff_o4_o5 = diff_o4_o5 >= 0 ? diff_o4_o5 : -diff_o4_o5;
      max_diff_o4_o5 = max_diff_o4_o5 > diff_o4_o5 ? max_diff_o4_o5 : diff_o4_o5;

    }
    double s_rk = h*rk_tolerance/(2*(tot_time)*max_diff_o4_o5);

    s_rk=pow(s_rk,0.25);
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

    for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
    {
      expression_gene_h[gene_count1]=expression_gene[gene_count1];


    }
    if((i> (print_start + print_interval*print_counter)) &&
       i <= (h+print_start + print_interval*print_counter))
    {
      print_counter++;
      //std::cout<<i<<"\n";
      for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
      {
        out_GE<<std::setprecision(output_precision)
              <<expression_gene[gene_count1]<<"\t";
      }
      //out_GE<<"\n";
    }
    //i=i+h;
  }while(i<tot_time);

}
