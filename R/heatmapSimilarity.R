
#' @export
#' @title Calculates the similarity between two gene expression data.
#' @description Comparison is done across columns, i.e., 
#' how similar are the columns in the two dataset.
#' For gene expression data, format data so that gene names are in rows and 
#' samples in columns.
#' @param dataReference Matrix. The reference data matrix, for example, 
#' the experimental gene expression values
#' @param dataSimulation Matrix. The data matrix to be compared.
#' @param nClusters (optional) Integer. The number of clusters in which the
#' reference data should be clustered for comparison.
#' Not needed if clusterCut is provided.
#' @param pValue (optional) Numeric. p-value to consider two gene expression
#'  sets as belonging to same cluster.
#' Ward's method with spearman correlation is used to determine if a
#' model belongs to a specific cluster.
#' @param permutedVar (optional) Similarity scores computed after permutations.
#' @param clusterCut (optional) Integer vector. Clsuter numbers assigned 
#' to reference data.
#' If clusterCut is missing, hierarchical clustering using /code{ward.D2}
#' and /code{distance  = (1-cor(x, method = "spear"))/2} will be used to 
#' cluster the reference data.
#' @param clusterMethod (optional) Character - default \code{ward.D2}, other
#' options include \code{complete}. Clustering method to be used to cluster the
#'  experimental data. \code{\link[stats]{hclust}} for other options.
#' @param corMethod (optional) Correlation method. Default method is "spearman".
#'  For single cell data, use "kendall"
#' @param permutMethod "sample" or "reference"
#' @param permutations (optional) Integer. Default \code{1000}.
#'  Number of gene permutations to generate the null distibution. 
#' @param returnData (optional) Logical. Default \code{FALSE}. Whether to
#' return the sorted and clustered data.
#' @param buffer (optional) Numeric. Default \code{0.001}. The fraction of
#' models to be assigned to clusters to which no samples could be assigned.
#' For example, a minimum of 1 ghost sample in reference is assigned to
#'  NULL cluster.
#' @param method (optional) character. Method to compare the gene expressions.
#'  Default \code{pvalue}. One can use \code{variance} as well which assigns
#'  clusters based on the cluster whose samples have minimum variance with
#'  the simulated sample.
#' @return A list containing the KL distance of new cluster distribution from 
#' reference data and
#' the probability of each cluster in the reference and simulated data.
#'
#' @section Related Functions:
#'
#' \code{\link{sracipeSimulate}},  \code{\link{sracipeKnockDown}},
#' \code{\link{sracipeOverExp}},  \code{\link{sracipePlotData}},
#' \code{\link{sracipeHeatmapSimilarity}}


sracipeHeatmapSimilarity = function(
  dataReference, dataSimulation, clusterCut = NULL, nClusters = 3, pValue=0.05, 
  permutedVar, permutations = 1000, corMethod = "spearman", 
  clusterMethod = "ward.D2", method = "pvalue", buffer = 0.001, 
  permutMethod = "simulation", returnData = FALSE) {
  
  #'
  
  
  ########
  # _7: Each model is compared with permutations of random data
  # Two new functions ClustFunction and ModelPValue
  # Assigns clusters based on p values
  # should work with all mean, min, z score type assignments 
  # as p values are related
  # A separate funtion to change the assignment type. Change ClustFunction
  
  commonGenes <- intersect(rownames(dataSimulation), rownames(dataReference))
  if(length(commonGenes) == 0) {
    message(" No Common genes found between the simulated and reference data.")
    return()
  }
 if(is.null(colnames(dataReference))){
   colnames(dataReference) <- seq_len(ncol(dataReference))
 }
  if(is.null(colnames(dataSimulation))){
    colnames(dataSimulation) <- seq_len(ncol(dataSimulation))
  }
  
  message("Calculating the similarity index")
  #  nClusters = 3
  n.models <- dim(dataReference)[2]
  nModelsKO <- dim(dataSimulation)[2]
  
  if (missing(permutations)) {
    permutations = 1000
  }
  
  if (missing(corMethod)) {
    corMethod <- "spearman"
  }
  
  refCor <- cor((dataReference), method = corMethod)
  
  if (missing(clusterCut)) {
    if(missing(nClusters)){
      stop("Please specify the number of clusters using nClusters or
           cluster assignments using clusterCut")
    }
    
    # cluster the reference data if the clutering assignments has not been 
    # provided.
    distance <- as.dist((1-refCor)/2)
    clusters <- hclust(distance, method = clusterMethod)
    #plot(clusters)
    clusterCut <- cutree(clusters, nClusters)
    
    } else {
      if(!missing(nClusters)){
        warnings("Neglecting nClusters. The number of clusters will be 
                 determined from clusterCut.")
      }
      nClusters <- length(unique(clusterCut))
    }
  # Use only selected genes for comparison.
  # Clustering is done using all genes. This can lead to differences in
  # how the clusters look.
  dataReference <- dataReference[commonGenes, ]
  dataSimulation <- dataSimulation[commonGenes, ]
  # find the variance within each cluster
  #TO DO Will standard deviation be better? shouldn't be with ward method.
  
  refClusterVar <- c(rep(0,nClusters))
  for(j in seq_len(nClusters))
  {
    #  print(j)
    temp.cluster.var <- (((1 - refCor[which(clusterCut==j), 
                                      which(clusterCut==j)])/2)^2)
    refClusterVar[j] <- .ClustFunction(
      temp.cluster.var[upper.tri(temp.cluster.var, diag = FALSE)])
    temp.cluster.var <- NULL
  }
  
  
  #  clusterCut <- clusterCut[1:10]
  #  dataReference <- dataReference[,1:10]
  simulated.refCor <- t(cor(dataReference, dataSimulation, method = corMethod))
  
  #clusterFreq <- table(CLUSTERCUT)/n_models
  
  if (sum(is.na(simulated.refCor)) > 0) {
    message("Error in correlation. Please verify the data")
  }
  
  simulatedClusterVar <- matrix(0, nrow=nModelsKO, ncol = nClusters)
  
  for(i in seq_len(nModelsKO)){
    for(j in seq_len(nClusters))
    {
      temp.cluster.var <- ((1 - simulated.refCor[i, which(clusterCut==j)])/2)^2
      simulatedClusterVar[i,j] <- .ClustFunction(temp.cluster.var )
      temp.cluster.var <- NULL
    }
  }
  

  if (method == "variance") {
    simulated.cluster <- matrix(0, nrow =  nModelsKO, ncol = 2)
    simulated.cluster[, 2] <- apply(simulatedClusterVar,1,min)
    # simulated.cluster.allowed <- simulatedClusterVar < refClusterVar
    simulated.cluster[, 1] <- apply(simulatedClusterVar,1,which.min)
    simulated.cluster[which(3*refClusterVar[simulated.cluster[,1]] < 
                              simulated.cluster[, 2]), 1] <- 0
    simulated.cluster <- simulated.cluster[,-2]
    
  }
  #  permutations = 1000
  if(missing(method)) {
    method = "pvalue"
  }
  if (method == "pvalue" ) {
    message("pvalue method")
    if(missing(permutedVar )) {
      if(permutMethod == "reference"){
        permutedVar <- .PermutedVar(simulated.refCor, clusterCut, permutations, 
                                   refClusterVar)
        simulatedVarPValue <- .SimulatedVarPValue(permutedVar, pValue)
        # rowSums(simulated.cluster.allowed)
        # simulatedClusterVar.sorted <- sort(simulatedClusterVar, 
        # index.return = TRUE )
        # simulated.cluster.allowed <- simulatedClusterVar < simulatedVarPValue
        simulated.cluster <- matrix(0, nrow =  nModelsKO, ncol = 2)
        simulated.cluster[, 2] <- apply(simulatedClusterVar,1,min)
        simulated.cluster[, 1] <- apply(simulatedClusterVar,1,which.min)
        simulated.cluster[which(simulatedVarPValue[simulated.cluster[,1]] < 
                                  simulated.cluster[, 2]), 1] <- 0
        simulated.cluster <- simulated.cluster[,-2]
        
      } else {
        message("simulation permutation")
        
        pValueMat <- .ModelPvalue(
          dataSimulation, dataReference, clusterCut, permutations,
          refClusterVar, corMethod, simulatedClusterVar)
        simulated.cluster <- matrix(0, nrow =  nModelsKO, ncol = 2)
        simulated.cluster[, 2] <- apply(pValueMat,1,min)
        simulated.cluster[, 1] <- apply(pValueMat,1,which.min)
        simulated.cluster[which(simulated.cluster[,2] > pValue), 1] <- 0
        simulated.cluster <- simulated.cluster[,-2]
        
      }
    }
    
    
  }
  
  similarity <- list()
  similarity$simClusters <- sort(simulated.cluster)
  similarity$simClusters <- c(similarity$simClusters[similarity$simClusters>0],similarity$simClusters[similarity$simClusters==0])
  cluster.names <- unique(clusterCut)
  #print(c("Original Clusters", cluster.names))
  cluster.names <- c(0, cluster.names) #test
  
  bufferEnteriesPerCluster <- max(1,as.integer(buffer*n.models))
  clusterCut.adjusted <- c(clusterCut, rep(0,bufferEnteriesPerCluster))
  
  simulated.cluster.names <- unique(simulated.cluster)
  # print(c("Simulated Clusters", simulated.cluster.names))
  missing.ref.clusters <- setdiff(cluster.names, simulated.cluster.names)
  #print(c("Missing Clusters", missing.ref.clusters))
  bufferEnteriesPerCluster <- max(1,as.integer(buffer*nModelsKO))
  missing.ref.clusters.add <- numeric() 
  #c(rep(0,bufferEnteriesPerCluster*length(missing.ref.clusters)))
  if (length(missing.ref.clusters) > 0) {
    for(i in seq_along(missing.ref.clusters))
    {
      missing.ref.clusters.add <- c(missing.ref.clusters.add, 
                                    rep(missing.ref.clusters[i],
                                        bufferEnteriesPerCluster))
    }
  }
  simulated.cluster.adjusted <- c(simulated.cluster, missing.ref.clusters.add)
  
  
  
  
  ref.cluster.freq <- table(clusterCut.adjusted)/(length(clusterCut.adjusted))
  # similarity$ref.cluster.freq <- table(clusterCut)/n.models
  similarity$ref.cluster.freq <- ref.cluster.freq
  
  simulated.cluster.freq <- 
    table(simulated.cluster.adjusted)/length(simulated.cluster.adjusted)
  
  #similarity$simulated.cluster.freq <- table(simulated.cluster)/nModelsKO
  similarity$simulated.cluster.freq <- simulated.cluster.freq
  
  similarity$cluster.similarity <- 
    simulated.cluster.freq*log(simulated.cluster.freq/ref.cluster.freq)
  similarity$KL <- sum(similarity$cluster.similarity )
  
  if(returnData){
     # similarity$dataReference <- dataReference
     dataRefSamples <- colnames(dataReference)
 #    print(dataRefSamples)
     dataRefSamples <- dataRefSamples[order(clusterCut)]
#     print(dataRefSamples)
     clusterCut <- clusterCut[order(clusterCut)]
     similarity$refClusters <- clusterCut
     
     #colnames(similarity$dataReference) <- clusterCut
     similarity$dataReference <-
      dataReference[,dataRefSamples]
     #print(dim(dataReference))
     # print(dim(similarity$dataReference))


     similarity$dataSimulation <- dataSimulation[,which(simulated.cluster>0)]
    colnames(similarity$dataSimulation) <-
      simulated.cluster[which(simulated.cluster>0)]
     similarity$dataSimulation <-
      similarity$dataSimulation[,order(colnames(similarity$dataSimulation))]
    refSimCor <- numeric()
    previous.cluster.size <- 0
    refSimCor.ref <- numeric()
    previous.cluster.size.ref <- 0
  #  print(colnames(similarity$dataSimulation))
    # print(clusterCut)
    for(i in seq_len((nClusters+1) ))#(length(unique(colnames(similarity$dataSimulation)))))
    {
   #   print(i)
      temp.ref <- similarity$dataReference[,which(
        clusterCut==i)]
      temp.sim <- similarity$dataSimulation[,which(
        colnames(similarity$dataSimulation)==i)]
      #similarity$simCluster <- colnames(similarity$dataSimulation)
      
      temp.refSimCor <- cor(temp.ref,temp.sim, method = corMethod)
      refSimCor <- c(refSimCor,previous.cluster.size +
                         sort(colMeans(temp.refSimCor), 
                              decreasing = TRUE, index.return = TRUE)$ix)
      previous.cluster.size <- previous.cluster.size + dim(temp.sim)[2]
      
      refSimCor.ref <- c(refSimCor.ref, previous.cluster.size.ref +
                             sort(rowMeans(temp.refSimCor), decreasing = TRUE, 
                                  index.return = TRUE)$ix)
      previous.cluster.size.ref <- previous.cluster.size.ref + dim(temp.ref)[2]
      
      
    }
    
    similarity$dataSimulation <- similarity$dataSimulation[,refSimCor]
    tmp <- dataSimulation[,which(simulated.cluster == 0)]
    colnames(tmp) <- rep(0, dim(tmp)[2])
    
    similarity$dataSimulation <- cbind(similarity$dataSimulation[,refSimCor],
                                       tmp)
    colnames(similarity$dataSimulation) <- seq_len(ncol(similarity$dataSimulation))
    #print(dim(similarity$dataReference))
    # similarity$dataReference <- similarity$dataReference[,refSimCor.ref]
    #print(dim(similarity$dataReference))
    #TO DO : This invovlves repeat calculation of cor--can be optimized
 #   print(similarity$dataReference)
 #   print(similarity$dataSimulation)
    similarity$simulated.refCor <- t(cor(similarity$dataReference, 
                                         similarity$dataSimulation, 
                                         method = corMethod))
  }
  #image(similarity$simulated.refCor, col = plot_color)
  return(similarity)
}

#########################################################
# Helper functions
#########################################################
#' @title Find nth minimum value from a vector
#' @description A utility function to find the nth minimum
#' @param x the given unsorted vector
#' @param index N.
#' @return the nth minimum element of the vector
#'
.NthMin <- function(x,index) {
  
  return (sort(x, decreasing = FALSE, partial = index)[index])
  
}

#############################################

.ClustFunction <- function(x){
  #return (mean(x))
  return (min(x))
}


#' @title Find variance of permutations
#' @description A utility function to generate permutations 
#' @param simulated.refCor Correlation matrix of simulated and reference data
#' @param clusterCut The original cluster assignments
#' @param permutations The number of permutations
#' @param refClusterVar Reference Cluster Variance
#' @return An array of dimension n.models by nClusters by permutations
#'
.PermutedVar <- function(simulated.refCor, clusterCut, permutations, 
                        refClusterVar){
  
  nClusters <- length(unique(clusterCut))
  nModelsKO <- dim(simulated.refCor)[1]
  permutedVar <- array(0, c(nModelsKO, nClusters, permutations))
  for(k in seq_len(permutations)){
    clusterCut.permuted <- sample(clusterCut)
    for(j in seq_len(nClusters))
    {
      cor.mat <- simulated.refCor[,which(clusterCut.permuted == j)]
      for(i in seq_len(nModelsKO)){
        temp.cluster.var <- (((1-cor.mat[i, ])/2)^2)
        permutedVar[i,j,k] <- (mean(temp.cluster.var)/ (refClusterVar[j]))
      }
    }
  }
  return(permutedVar)
}
#############################################
#' @title Returns the variance array after permutations.
#' @description A utility function.
#' @param dataSimulation Simulation data to be compared.
#' @param dataReference Reference data with which to compare.
#' @param clusterCut The original cluster assignments.
#' @param permutations The number of permutations.
#' @param refClusterVar SD of the clusters.
#' @param corMethod Correlation method to be used.
#' @param simulatedClusterVar Variance of simulated clusters
#' @return An array of dimension n.models by nClusters by permutations.
#'
.ModelPvalue <- function(dataSimulation, dataReference, clusterCut, permutations,
                        refClusterVar, corMethod, simulatedClusterVar){
  
  nClusters <- length(unique(clusterCut))
  nModelsKO <- dim(dataSimulation)[2]
  n.gene <- dim(dataSimulation)[1]
  pValueMat <- matrix(0, nrow = nModelsKO, ncol = nClusters)
  randomModels <-  matrix(rep(seq_len(n.gene),permutations), n.gene, permutations)
  randomModels <- apply(randomModels,2,sample)
  #randomModels <- matrix(0, nrow = n.gene, ncol = permutations)
  #randomModels <- t(apply(dataSimulation,1,function(x)  sample(
  #x, replace = TRUE, size = permutations)))
  permutedRefCor <- matrix(0,nrow = permutations, ncol = dim(dataReference)[2])
  permutedRefCor <- cor(randomModels, dataReference, method = corMethod)
  for(j in seq_len(nClusters))
  {
    dist.mat <- ((1 - permutedRefCor[,which(clusterCut == j)])/2)^2
    tempVector <- sort(apply(dist.mat,1,.ClustFunction))
    for (i in seq_len(nModelsKO)) {
      pValueMat[i,j] <- (which(abs(
        tempVector - simulatedClusterVar[i,j]) ==min(abs(
          tempVector - simulatedClusterVar[i,j])))[1] - 1)/permutations
      #[1] as sometimes which() might satisfy for multiple values
    }
    
  }
  
  return(pValueMat)
}

#############################################
#' @title Finds the variance corresponding to a given value.
#' @description  A utility function for calculating the p values.
#' @param permutedVar An array containing the distance of clusters for 
#' each model for every permutation.
#' @param pValue Cut off p vlaue.
#' @return p-values for each model.
#'

.SimulatedVarPValue <- function(permutedVar, pValue){
  
  permutations <- dim(permutedVar)[3]
  nModelsKO <- dim(permutedVar)[1]
  nClusters <- dim(permutedVar)[2]
  selectedIndex = as.integer(permutations*pValue)
  
  if(selectedIndex==0) {stop("Number of permutations is not sufficient 
to achieve the required pValue.
                              Please increase the permutations")}
  #print(selectedIndex)
  simulatedVarPValue <- matrix(0, nrow=nModelsKO, ncol = nClusters)
  
  for (i in seq_len(nModelsKO)) {
    for(j in seq_len(nClusters))    {
      simulatedVarPValue[i,j] <- .NthMin(permutedVar[i,j,],selectedIndex)
    }
  }
  
  return(simulatedVarPValue)
  }
#############################################

#' @title Finds the variance corresponding to a given value.
#' @description A utility function to calculate the absolute p values.
#' @param permutedVar An array containing the distance of clusters for 
#' each model for every permutation.
#' @param simulatedClusterVar Variance of simulated clusters
#' @return p-values for each model.
#'
.SimulatedPValueAbs <- function(permutedVar, simulatedClusterVar){
  
  permutations <- dim(permutedVar)[3]
  nModelsKO <- dim(permutedVar)[1]
  nClusters <- dim(permutedVar)[2]
  
  simulatedVarPValue <- matrix(0, nrow=nModelsKO, ncol = nClusters)
  
  for (i in seq_len(nModelsKO)) {
    for(j in seq_len(nClusters) )   {
      tempVector <- sort(permutedVar[i,j,])
      
      simulatedVarPValue[i,j] <- which(abs(
        tempVector - simulatedClusterVar[i,j])==min(abs(
          tempVector - simulatedClusterVar[i,j])))/permutations
    }
  }
  
  return(simulatedVarPValue)
}



