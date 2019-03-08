
#' @export
#' @title Calculates the similarity between two gene expression data.
#' @description Comparison is done across columns, i.e., how similar are the columns in the two dataset.
#' For gene expression data, format data so that gene names are in rows and samples in columns.
#' @param data.reference Matrix. The reference data matrix, for example, the experimental gene expression values
#' @param data.simulation Matrix. The data matrix to be compared.
#' @param n.clusters (optional) Integer. The number of clusters in which the
#' reference data should be clustered for comparison.
#' Not needed if cluster.cut is provided.
#' @param p.value (optional) Numeric. p-value to consider two gene expression
#'  sets as belonging to same cluster.
#' Ward's method with spearman correlation is used to determine if a
#' model belongs to a specific cluster.
#' @param permuted.var (optional) Similarity scores computed after permutations.
#' @param cluster.cut (optional) Integer vector. Clsuter numbers assigned to reference data.
#' If cluster.cut is missing, hierarchical clustering using /code{ward.D2}
#' and /code{distance  = (1-cor(x, method = "spear"))/2} will be used to cluster the reference data.
#' @param cluster.method (optional) Character - default \code{ward.D2}, other
#' options include \code{complete}. Clustering method to be used to cluster the experimental data.
#' \code{\link[stats]{hclust}} for other options.
#' @param cor.method (optional) Correlation method. Default method is "spearman". For single cell data, use "kendall"
#' @param permutation.method "sample" or "reference"
#' @return A list containing the KL distance of new cluster distribution from reference data and
#' the probability of each cluster in the reference and simulated data.
#'
#' @section Related Functions:
#'
#' \code{\link{simulateGRC}},  \code{\link{knockdownAnalysis}},
#' \code{\link{overExprAnalysis}},  \code{\link{plotData}},
#' \code{\link{calcHeatmapSimilarity}}


calcHeatmapSimilarity = function(data.reference, data.simulation, cluster.cut = NULL, n.clusters = 3,
                             p.value=0.05, permuted.var, permutations = 1000,
                             cor.method = "spearman", cluster.method = "ward.D2", method = "pvalue", buffer = 0.0001, permutation.method = "simulation", return.data = FALSE) {

  #'


  ########
  # _7: Each model is compared with permutations of random data
  # Two new functions ClustFunction and ModelPValue
  # Assigns clusters based on p values
  # should work with all mean, min, z score type assignments as p values are related
  # A separate funtion to change the assignment type. Change ClustFunction

  message("Calculating the similarity index")
  #  n.clusters = 3
  n.models <- dim(data.reference)[2]
  n.models.KO <- dim(data.simulation)[2]

  if (missing(permutations)) {
    permutations = 1000
  }

  if (missing(cor.method)) {
    cor.method <- "spearman"
  }

  ref.cor <- cor((data.reference), method = cor.method)

  if (missing(cluster.cut)) {
    if(missing(n.clusters)){
      stop("Please specify the number of clusters using n.clusters or
           cluster assignments using cluster.cut")
    }

    # cluster the reference data if the clutering assignments has not been provided.
    distance <- as.dist((1-ref.cor)/2)
    clusters <- hclust(distance, method = cluster.method)
    #plot(clusters)
    cluster.cut <- cutree(clusters, n.clusters)

    } else {
      if(!missing(n.clusters)){
        warnings("Neglecting n.clusters. The number of clusters will be determined from cluster.cut.")
      }
      n.clusters <- length(unique(cluster.cut))
    }

  # find the variance within each cluster
  #TO DO Will standard deviation be better? shouldn't be with ward method.

  ref.cluster.var <- c(rep(0,n.clusters))
  for(j in 1:n.clusters)
  {
  #  print(j)
    temp.cluster.var <- (((1 - ref.cor[which(cluster.cut==j), which(cluster.cut==j)])/2)^2)
    ref.cluster.var[j] <- ClustFunction(temp.cluster.var[upper.tri(temp.cluster.var, diag = FALSE)])
    temp.cluster.var <- NULL
  }


  #  cluster.cut <- cluster.cut[1:10]
  #  data.reference <- data.reference[,1:10]
  simulated.ref.cor <- t(cor(data.reference, data.simulation, method = cor.method))

  #clusterFreq <- table(CLUSTERCUT)/n_models

  if (sum(is.na(simulated.ref.cor)) > 0) {
    message("Error in correlation. Please verify the data")
  }

  simulated.cluster.var <- matrix(0, nrow=n.models.KO, ncol = n.clusters)

  for(i in 1:n.models.KO){
    for(j in 1:n.clusters)
    {
      temp.cluster.var <- ((1 - simulated.ref.cor[i, which(cluster.cut==j)])/2)^2
      simulated.cluster.var[i,j] <- ClustFunction(temp.cluster.var )
      temp.cluster.var <- NULL
    }
  }

  # testing clustering robustness
  # ref.cluster.var <- matrix(0, nrow = n.models, ncol =  n.clusters)
  # for(i in 1:n.models)
  # for(j in 1:n.clusters)
  # {
  #   ref.cluster.var[i,j] <- mean(((1 - ref.cor[i, which(cluster.cut==j)])/2)^2)
  # }

  if (method == "variance") {
    simulated.cluster <- matrix(0, nrow =  n.models.KO, ncol = 2)
    simulated.cluster[, 2] <- apply(simulated.cluster.var,1,min)
    # simulated.cluster.allowed <- simulated.cluster.var < ref.cluster.var
    simulated.cluster[, 1] <- apply(simulated.cluster.var,1,which.min)
    simulated.cluster[which(3*ref.cluster.var[simulated.cluster[,1]] < simulated.cluster[, 2]), 1] <- 0
    simulated.cluster <- simulated.cluster[,-2]

  }
  #  permutations = 1000
  if(missing(method)) {
    method = "pvalue"
  }
  if (method == "pvalue" ) {
    message("pvalue method")
    if(missing(permuted.var )) {
      if(permutation.method == "reference"){
      permuted.var <- PermutedVar(simulated.ref.cor, cluster.cut, permutations, ref.cluster.var)
      simulated.var.P.value <- SimulatedVarPValue(permuted.var, p.value)
      #rowSums(simulated.cluster.allowed)
      #simulated.cluster.var.sorted <- sort(simulated.cluster.var, index.return = TRUE )
      # simulated.cluster.allowed <- simulated.cluster.var < simulated.var.P.value
      simulated.cluster <- matrix(0, nrow =  n.models.KO, ncol = 2)
      simulated.cluster[, 2] <- apply(simulated.cluster.var,1,min)
      simulated.cluster[, 1] <- apply(simulated.cluster.var,1,which.min)
      simulated.cluster[which(simulated.var.P.value[simulated.cluster[,1]] < simulated.cluster[, 2]), 1] <- 0
      simulated.cluster <- simulated.cluster[,-2]

      } else {
        message("simulation permutation")

        p.value.mat <- ModelPvalue(data.simulation, data.reference, cluster.cut, permutations,
                                   ref.cluster.var, cor.method, simulated.cluster.var)
        simulated.cluster <- matrix(0, nrow =  n.models.KO, ncol = 2)
        simulated.cluster[, 2] <- apply(p.value.mat,1,min)
        simulated.cluster[, 1] <- apply(p.value.mat,1,which.min)
        simulated.cluster[which(simulated.cluster[,2] > p.value), 1] <- 0
        simulated.cluster <- simulated.cluster[,-2]

        }
    }


  }

  similarity <- list()
  cluster.names <- unique(cluster.cut)
  #print(c("Original Clusters", cluster.names))
  cluster.names <- c(0, cluster.names) #test

  bufferEnteriesPerCluster <- max(1,as.integer(buffer*n.models))
  cluster.cut.adjusted <- c(cluster.cut, rep(0,bufferEnteriesPerCluster))

  simulated.cluster.names <- unique(simulated.cluster)
 # print(c("Simulated Clusters", simulated.cluster.names))
  missing.ref.clusters <- setdiff(cluster.names, simulated.cluster.names)
  #print(c("Missing Clusters", missing.ref.clusters))
  bufferEnteriesPerCluster <- max(1,as.integer(buffer*n.models.KO))
  missing.ref.clusters.add <- numeric() #c(rep(0,bufferEnteriesPerCluster*length(missing.ref.clusters)))
  if (length(missing.ref.clusters) > 0) {
  for(i in 1:length(missing.ref.clusters))
  {
    missing.ref.clusters.add <- c(missing.ref.clusters.add, rep(missing.ref.clusters[i],bufferEnteriesPerCluster))
  }
  }
  simulated.cluster.adjusted <- c(simulated.cluster, missing.ref.clusters.add)




  ref.cluster.freq <- table(cluster.cut.adjusted)/(length(cluster.cut.adjusted))
  # similarity$ref.cluster.freq <- table(cluster.cut)/n.models
  similarity$ref.cluster.freq <- ref.cluster.freq

  simulated.cluster.freq <- table(simulated.cluster.adjusted)/length(simulated.cluster.adjusted)

  #similarity$simulated.cluster.freq <- table(simulated.cluster)/n.models.KO
  similarity$simulated.cluster.freq <- simulated.cluster.freq

  similarity$cluster.similarity <- simulated.cluster.freq*log(simulated.cluster.freq/ref.cluster.freq)
  similarity$KL <- sum(similarity$cluster.similarity )

  if(return.data){
  similarity$data.reference <- data.reference
  colnames(similarity$data.reference) <- cluster.cut
  similarity$data.reference <- similarity$data.reference[,order(colnames(similarity$data.reference))]


  similarity$data.simulation <- data.simulation[,which(simulated.cluster>0)]
  colnames(similarity$data.simulation) <- simulated.cluster[which(simulated.cluster>0)]
  similarity$data.simulation <- similarity$data.simulation[,order(colnames(similarity$data.simulation))]
  ref.sim.cor <- numeric()
  previous.cluster.size <- 0
  ref.sim.cor.ref <- numeric()
  previous.cluster.size.ref <- 0

  for(i in 1:(length(unique(colnames(similarity$data.simulation)))))
  {
    temp.ref <- similarity$data.reference[,which(colnames(similarity$data.reference)==i)]
    temp.sim <- similarity$data.simulation[,which(colnames(similarity$data.simulation)==i)]

    temp.ref.sim.cor <- cor(temp.ref,temp.sim, method = cor.method)
    ref.sim.cor <- c(ref.sim.cor,previous.cluster.size +
                       sort(colMeans(temp.ref.sim.cor), decreasing = T, index.return = T)$ix)
    previous.cluster.size <- previous.cluster.size + dim(temp.sim)[2]

    ref.sim.cor.ref <- c(ref.sim.cor.ref, previous.cluster.size.ref +
                       sort(rowMeans(temp.ref.sim.cor), decreasing = T, index.return = T)$ix)
    previous.cluster.size.ref <- previous.cluster.size.ref + dim(temp.ref)[2]


  }

  similarity$data.simulation <- similarity$data.simulation[,ref.sim.cor]
  tmp <- data.simulation[,which(simulated.cluster == 0)]
  colnames(tmp) <- rep(0, dim(tmp)[2])

  similarity$data.simulation <- cbind(similarity$data.simulation[,ref.sim.cor], tmp)

  similarity$data.reference <- similarity$data.reference[,ref.sim.cor.ref]

  #TO DO : This invovlves repeat calculation of cor--can be optimized
  similarity$simulated.ref.cor <- t(cor(similarity$data.reference, similarity$data.simulation, method = cor.method))
}
  #image(similarity$simulated.ref.cor, col = plot_color)
  return(similarity)
}

#########################################################
# Helper functions
#########################################################
#' @title Find nth minimum value from a vector
#' @description A utility function to find the nth minimum
#'
#' @param x the given unsorted vector
#' @param index N.
#' @return the nth minimum element of the vector
#'
NthMin <- function(x,index) {

  return (sort(x, decreasing = FALSE, partial = index)[index])

}

#############################################

ClustFunction <- function(x){
  #return (mean(x))
  return (min(x))
}

#############################################
#' @title Find variance of permutations
#' @description A utility function to generate permutations
#'
#' @param simulated.ref.cor Correlation matrix of simulated and reference data
#' @param cluster.cut The original cluster assignments
#' @param permutations The number of permutations
#' @return An array of dimension n.models by n.clusters by permutations
#'
PermutedVar <- function(simulated.ref.cor, cluster.cut, permutations, ref.cluster.var){

  n.clusters <- length(unique(cluster.cut))
  n.models.KO <- dim(simulated.ref.cor)[1]
  permuted.var <- array(0, c(n.models.KO, n.clusters, permutations))
  for(k in 1:permutations){
    cluster.cut.permuted <- sample(cluster.cut)
    for(j in 1:n.clusters)
    {
      cor.mat <- simulated.ref.cor[,which(cluster.cut.permuted == j)]
      for(i in 1:n.models.KO){
        temp.cluster.var <- (((1-cor.mat[i, ])/2)^2)
        permuted.var[i,j,k] <- (mean(temp.cluster.var)/ (ref.cluster.var[j]))
      }
    }
  }
  return(permuted.var)
}
#############################################
#' @title Returns the variance array after permutations.
#' @description A utility function.
#' @param data.simulation Simulation data to be compared.
#' @param data.reference Reference data with which to compare.
#' @param cluster.cut The original cluster assignments.
#' @param permutations The number of permutations.
#' @param ref.cluster.var SD of the clusters.
#' @return An array of dimension n.models by n.clusters by permutations.
#'
ModelPvalue <- function(data.simulation, data.reference, cluster.cut, permutations,
                        ref.cluster.var, cor.method, simulated.cluster.var){

  n.clusters <- length(unique(cluster.cut))
  n.models.KO <- dim(data.simulation)[2]
  n.gene <- dim(data.simulation)[1]
  p.value.mat <- matrix(0, nrow = n.models.KO, ncol = n.clusters)
  random.models <-  matrix(rep(seq(1:n.gene),permutations), n.gene, permutations)
  random.models <- apply(random.models,2,sample)
  #random.models <- matrix(0, nrow = n.gene, ncol = permutations)
  #random.models <- t(apply(data.simulation,1,function(x)  sample(x, replace = TRUE, size = permutations)))
  permuted.ref.cor <- matrix(0,nrow = permutations, ncol = dim(data.reference)[2])
  permuted.ref.cor <- cor(random.models, data.reference, method = cor.method)
  for(j in 1:n.clusters)
  {
    dist.mat <- ((1 - permuted.ref.cor[,which(cluster.cut == j)])/2)^2
    temp.vector <- sort(apply(dist.mat,1,ClustFunction))
    for (i in 1:n.models.KO) {
      p.value.mat[i,j] <- (which(abs(temp.vector - simulated.cluster.var[i,j])==min(abs(temp.vector - simulated.cluster.var[i,j])))[1] - 1)/permutations
      #[1] as sometimes which() might satisfy for multiple values
    }

  }

  return(p.value.mat)
}

#############################################
#' @title Finds the variance corresponding to a given value.
#' @description  A utility function for calculating the p values.
#'
#' @param permuted.var An array containing the distance of clusters for each model for every permutation.
#' @param p.value Cut off p vlaue.
#' @return p-values for each model.
#'
SimulatedVarPValue <- function(permuted.var, p.value){

  permutations <- dim(permuted.var)[3]
  n.models.KO <- dim(permuted.var)[1]
  n.clusters <- dim(permuted.var)[2]
  selected.index = as.integer(permutations*p.value)

  if(selected.index==0) {stop("Number of permutations is not sufficient to achieve the required p.value.
                              Please increase the permutations")}
  #print(selected.index)
  simulated.var.P.value <- matrix(0, nrow=n.models.KO, ncol = n.clusters)

  for (i in 1:n.models.KO) {
    for(j in 1:n.clusters)    {
      simulated.var.P.value[i,j] <- NthMin(permuted.var[i,j,],selected.index)
    }
  }

  return(simulated.var.P.value)
  }
#############################################

#' @title Finds the variance corresponding to a given value.
#' @description A utility function to calculate the absolute p values.
#' @param permuted.var An array containing the distance of clusters for each model for every permutation.
#' @param p.value Cut off p vlaue.
#' @return p-values for each model.
#'
SimulatedPValueAbs <- function(permuted.var, simulated.cluster.var){

  permutations <- dim(permuted.var)[3]
  n.models.KO <- dim(permuted.var)[1]
  n.clusters <- dim(permuted.var)[2]

  simulated.var.P.value <- matrix(0, nrow=n.models.KO, ncol = n.clusters)

  for (i in 1:n.models.KO) {
    for(j in 1:n.clusters)    {
      temp.vector <- sort(permuted.var[i,j,])

      simulated.var.P.value[i,j] <- which(abs(temp.vector - simulated.cluster.var[i,j])==min(abs(temp.vector - simulated.cluster.var[i,j])))/permutations
    }
  }

  return(simulated.var.P.value)
}



