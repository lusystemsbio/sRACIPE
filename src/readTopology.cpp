#include "header.h"
#include <Rcpp.h>

using namespace Rcpp;


// [[Rcpp::export]]

Rcpp::IntegerMatrix readTopology(
    Rcpp::IntegerMatrix gene_interaction, const Rcpp::String filepath,
    const Rcpp::String filename, Rcpp::StringVector geneNames)
  {

  std::string file_name = filename;
//  std::fstream out_2("./tmp/gene_interaction_topology_" + file_name + ".txt",
//                     std::fstream::out);
  file_name = filepath;
  std::ifstream infile ( file_name, std::ifstream::in);
  if(!infile) {
    Rcout << "Cannot open input file.\n";
    return 1;
  }
  std::vector <std::string> gene_names;
  std::string word;
  int i=0; int first_gene=0,second_gene=0;
  while (infile >> word)
  {
    if(i<3){}
    else if(i%3!=2)
    {
      if (std::find(gene_names.begin(), gene_names.end(),word)!=gene_names.end())
      {
        if(i%3==0)
        {first_gene = distance(gene_names.begin(), find(gene_names.begin(),
                                                gene_names.end(), word));}
        else second_gene = distance(gene_names.begin(), find(gene_names.begin(),
                                                     gene_names.end(), word));
      }
      else
      {
        gene_names.push_back(word);
        if((i%3==0))first_gene=gene_names.size()-1;
        else second_gene=gene_names.size()-1;
      }

    }
    else gene_interaction(second_gene,first_gene)=std::stoi(word);
    i++;
  }

//  out_2<<"gene"<<"\n";
  for (int i=0; i<gene_names.size(); i++)
  {
//    out_2 << gene_names[i] << "\n";
    geneNames(i) = gene_names[i];
  }

  return gene_interaction;

}
