rm(list = ls())
working_directory <- getwd()
topology_file = "inputs/test.tpo"

config_file = "inputs/sRACIPE.cfg"
configuration <- sRACIPE_load_configuration(config_file = config_file)

topology <- sRACIPE_load_topology(topology_file = topology_file)

results_directory <- ifelse(!dir.exists(file.path(working_directory, "results")), dir.create(file.path(working_directory, "results")), file.path(working_directory, "results"))
plot_network(topology)

output_file <- sRACIPE_RK_deterministic(topology_file = topology_file)

#output_file <- sRACIPE_RK_adaptive_deterministic(topology_file = topology_file)

#output_file <- sRACIPE_stochastic(topology_file = topology_file)

plot_data_deterministic(output_file, plot_filename = topology$filename)

#plot_data_stochastic(output_file, plot_filename=topology$filename, topology_df=topology, config = configuration, bin_count=40)

plot_interactive_heatmap(output_file, plot_filename = topology$filename)

