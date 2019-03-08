# Load the package
library(sRACIPE)

# Call the simulate function. It will simulate the circuit and plot the
# simulated data. The results will be saved in the results folder in the
# working directory. Use different arguments to modify the simulation parameters.

rSet <- sRACIPE::simulateGRC(circuit = "inputs/test.tpo")

# Perform in-silico knockdown analysis and plot the results.
sRACIPE::knockdownAnalysis(rSet = rSet)


# Perform in-silico knockdown analysis and plot the results.
sRACIPE::overExprAnalysis(rSet = rSet)
