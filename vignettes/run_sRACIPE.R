rm(list = ls())
library(sRACIPEv03)
working_directory <- getwd()
config_file = "inputs/sRACIPE.cfg"
configuration <- sRACIPE_load_configuration(config_file = config_file)

topology_file <- "inputs/test.tpo"
topology <- sRACIPE_load_topology(topology_file = topology_file)

plot_network(topology)

output_file <- sRACIPE_RK_deterministic(topology_file = topology_file)

plot_data_deterministic(output_file, plot_filename = topology$filename)

# Works for data generated using sRACIPE_RK_deterministic, sRACIPE_RK_adaptive_deterministic

# To use adaptive time step methjod
#output_file <- sRACIPE_RK_adaptive_deterministic(topology_file = topology_file)

# For stochastic simulations
#output_file <- sRACIPE_stochastic(topology_file = topology_file)

# For stochastic simulations and evaluation of whether the steady state soultions are fixed point or not
#output_file <- sRACIPE_stochastic_multiprint(topology_file = topology_file, NUM_MODELS = 10,ANNEAL = T, NOISE_LEVELS = 4, PRINT_START = 40, PRINT_INTERVAL = 4)


#plot_data_stochastic(output_file, plot_filename=topology$filename, topology_df=topology, config = configuration, bin_count=40)

## Works for data generated using sRACIPE_stochastic, sRACIPE_stochastic_multiprint

#plot_interactive_heatmap(output_file, plot_filename = topology$filename)

# ###########################################################################
# #Knockout Simulations
#
# topology_file <- "/Users/koharv/Documents/Work/ScenicMay/IPSC/output/IPSC_network_w001_corr002_Target.txt"
# topology <- sRACIPE_load_topology(topology_file)
# sRACIPEv03::plot_network(topology)
# output_file <- sRACIPEv03::sRACIPE_RK_deterministic(topology_file =  topology_file, NUM_MODELS = 5000)
#
# KO = "MAZ"
# output_file_knockout <- sRACIPEv03::sRACIPE_RK_deterministic_knockout(topology_file =  topology_file, NUM_MODELS = 5000, KNOCKOUT = KO)
# sRACIPEv03::plot_data_knockout_single(output_file, output_file_knockout, KNOCKOUT = KO,plot_filename = topology$filename, topology_df = topology)
# #sRACIPEv03::plot_interactive_heatmap(output_file,plot_filename = topology$filename, topology_df = topology)
#
