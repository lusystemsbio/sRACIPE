rm(list = ls())
library(sRACIPE)
working_directory <- getwd()
config_file = "inputs/sRACIPE.cfg"
configuration <- sRACIPE_load_configuration(config_file = config_file)

topology_file <- "inputs/test.tpo"
topology <- sRACIPE_load_topology(topology_file = topology_file)

output_file <- sRACIPE_RK_deterministic(topology_file = topology_file)

plot_data_deterministic(output_file, plot_filename = topology$filename)

# Works for data generated using sRACIPE_RK_deterministic, sRACIPE_RK_adaptive_deterministic

# To use adaptive time step methjod
#output_file <- sRACIPE_RK_adaptive_deterministic(topology_file = topology_file)

# For stochastic simulations
#output_file <- sRACIPE_stochastic(topology_file = topology_file)

#plot_data_stochastic(output_file, plot_filename=topology$filename, topology_df=topology, config = configuration, bin_count=40)

## Works for data generated using sRACIPE_stochastic, sRACIPE_stochastic_multiprint

#plot_network(topology)
