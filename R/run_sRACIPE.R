

print("Thank you for using sRACIPE!")
working_directory<-getwd()
input_directory<-file.path(working_directory, "inputs")
setwd(input_directory)
print("Reading configuration file from inputs folder")

f_inputs <- read.table(dir(pattern = ".cfg"),
                       header = TRUE)

f_inputs

print("Reading topology file from inputs folder")

tpo_filename<-dir(pattern = ".tpo")
f_tpo <- read.table(tpo_filename, header = TRUE)
f_tpo
setwd(working_directory)

output_directory<-file.path(working_directory, "results")
#setwd(output_directory)

g_min<-f_inputs[f_inputs[,1]=="MPR_MIN",2]
g_max<-f_inputs[f_inputs[,1]=="MPR_MAX",2]
k_min<-f_inputs[f_inputs[,1]=="DNR_MIN",2]
k_max<-f_inputs[f_inputs[,1]=="DNR_MAX",2]
lambda_min<-f_inputs[f_inputs[,1]=="FCH_MIN",2]
lambda_max<-f_inputs[f_inputs[,1]=="FCH_MAX",2]
n_min<-f_inputs[f_inputs[,1]=="HCO_MIN",2]
n_max<-f_inputs[f_inputs[,1]=="HCO_MAX",2]
h<-f_inputs[f_inputs[,1]=="EULER_STEP_SIZE",2]
tot_time <-f_inputs[f_inputs[,1]=="EULER_SIM_TIME",2]

number_gene<-length(levels(f_tpo$Source))
#possible_interactions=1+(unique(f_tpo$Type))
possible_interactions=3;
median_range <-f_inputs[f_inputs[,1]=="MEDIAN_RANGE",2]
model_count_max <-f_inputs[f_inputs[,1]=="NUM_MODELS",2]
threshold_max <-f_inputs[f_inputs[,1]=="THRESHOLD_MODELS",2]
threshold_max=max(model_count_max,threshold_max)
D_levels<-f_inputs[f_inputs[,1]=="NOISE_LEVELS",2]
D_max<-f_inputs[f_inputs[,1]=="MAX_NOISE",2]/sqrt(number_gene)
D_scaling<-f_inputs[f_inputs[,1]=="NOISE_SCALING_FACTOR",2]
standard_deviation_factor=0.5658;
D_shot_scaling<-f_inputs[f_inputs[,1]=="SHOT_NOISE_SCALING",2]
gene_noise_scaling<-f_inputs[f_inputs[,1]=="GENE_NOISE_SCALING",2]
file_writing_interval<-f_inputs[f_inputs[,1]=="FILE_WRITING_INTERVAL",2]
output_precision<-f_inputs[f_inputs[,1]=="OUTPUT_PRECISION",2]
ANNEALING<-f_inputs[f_inputs[,1]=="ANNEALING",2]
#CONSTANT_NOISE<-f_inputs[f_inputs[,1]=="CONSTANT_NOISE",2]
CONSTANT_NOISE<-1
GENE_NOISE_SCALING<-f_inputs[f_inputs[,1]=="GENE_NOISE_SCALING",2]
n_IC<-f_inputs[f_inputs[,1]=="INITIAL_CONDITIONS",2]


gene_interaction <- matrix(0, nrow = number_gene, ncol = number_gene)
gene_interaction
threshold_gene<- rep(0, number_gene)



Rcpp::sourceCpp("interaction_reader.cpp")


gene_interaction<-interaction_reader(gene_interaction, tpo_filename, number_gene)

print("Generating gene thresholds")

#
# g_min=1.
# g_max=100.
# k_min=0.1
# k_max=1.
# lambda_min=1.
# lambda_max=100.;

# n_min=1
# n_max=6;
# h = .02;
# tot_time = 50.1	#//number of iterations to be carried out
# number_gene=2;
#
# median_range=100
# model_count_max=1000
# threshold_max=1000
# D_counter=1
# D=1.0
# gene_interaction <- matrix(c(1, 2, 2, 0), nrow = 2)
# threshold_gene<- c(0,0)

Rcpp::sourceCpp("threshold_generator.cpp")

threshold_test<- threshold_calculator_uniform( gene_interaction,  threshold_gene,  g_min,  g_max,
                              k_min,  k_max,  possible_interactions,   model_count_max,  threshold_max,  h,  lambda_min,
                              lambda_max,  n_min,  n_max,  median_range,  standard_deviation_factor)
if(threshold_test!=0)
  print("Error in time evolution function")

# Rcpp::sourceCpp("threshold_test.cpp")
#
# threshold_gene2<-threshold_calculator_uniform_test( gene_interaction,  threshold_gene,  g_min,  g_max,
#                                    k_min,  k_max,  possible_interactions,
#                                    model_count_max, threshold_max,  h,  lambda_min,
#                                    lambda_max,  n_min,  n_max,  tot_time,  median_range,    standard_deviation_factor,
#                                    number_gene,  D,  D_counter,  D_shot_scaling,
#                                    GENE_NOISE_SCALING,  file_writing_interval,  D_levels,  D_scaling,  output_precision)

print("Running the simulations for models")

Rcpp::sourceCpp("multiGeneCircuit_EM_uniform_Darray_annealing.cpp")

Time_evolution_test<- multiGeneCircuit_EM_uniform_Darray_annealing( gene_interaction,  threshold_gene,  g_min,  g_max,
                                                   k_min,  k_max,  possible_interactions,
                                                   model_count_max, threshold_max,  h,  lambda_min,
                                                   lambda_max,  n_min,  n_max,  tot_time,  median_range,    standard_deviation_factor,
                                                   number_gene,  D_max,  D_shot_scaling,
                                                   GENE_NOISE_SCALING,  file_writing_interval,  D_levels,  D_scaling,  output_precision, ANNEALING, CONSTANT_NOISE)

if(Time_evolution_test!=0)
  print("Error in time evolution function")
