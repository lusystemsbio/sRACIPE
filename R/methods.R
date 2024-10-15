#' @export
#' @rdname sracipeCircuit
#' @aliases sracipeCircuit
setMethod(f="sracipeCircuit",
          signature="RacipeSE",
          definition=function(.object)
          {
            interactions <- as.matrix(rowData(.object))
            circuitInteractions <- getInteractions(interactions)
            return(circuitInteractions)
          }
)
#' @rdname sracipeCircuit-set
#' @aliases sracipeCircuit-set
setMethod("sracipeCircuit<-", "RacipeSE",
          function(.object, value) {
            circuitTable <- data.frame()
            filename <- character()
            if(is(value, "character") | is(value, "data.frame")){
            if(is(value, "character")){
              if(file.exists(value)){
                circuitTable <- utils::read.table(value, header = TRUE,
                                           stringsAsFactors = FALSE)
                colnames(circuitTable) <- c("Source","Target","Type")
                filename <- basename(value)

                filename <- tools::file_path_sans_ext(filename)

              }
              else {
                message("File not found!")
              }
            }
            if(is(value, "data.frame")){
              filename <- deparse(substitute(value))
              if(dim(value)[2]!=3)
              {
                message("Incorrect number of columns in circuit")
                return()
                }
              storage.mode(value[,3]) <- "integer"
              if(sum(!(value[,3] %in% c(1,2,3,4,5,6)))>0){
                message("Incorrect interactions (only 1,2,3,4,5,6 allowed)")
                return()
              }
              circuitTable <- value
              colnames(circuitTable) <- c("Source","Target","Type")
              }
            }
            else{
              stop("Incorrect circuit!
                   The circuit should either be a dataframe or filename.")
            }
            ## Check for special characters in gene names
            genes <- circuitTable$Source
            genes <- gsub("[[:punct:]]", "", genes)
            circuitTable$Source <- genes
            genes <- circuitTable$Target
            genes <- gsub("[[:punct:]]", "", genes)
            circuitTable$Target <- genes
            circuitGenes <-  unique(c(circuitTable$Source,circuitTable$Target))

            nGenes <- length(circuitGenes)
            nInteraction <- length(circuitTable$Source)
            geneTypes = rep(1, nGenes)
            circuitAdjMat <- matrix(data = 0,nrow = nGenes,
                                               ncol = nGenes)
            storage.mode(circuitAdjMat) <- "integer"
            rownames(circuitAdjMat) <- circuitGenes
            colnames(circuitAdjMat) <- circuitGenes

            for(i in seq_len(dim(circuitTable)[1])){
              circuitAdjMat[circuitTable[i,2], circuitTable[i,1]] <-
                circuitTable[i,3]
            }

            for(i in seq_len(nGenes)){
              if (all(circuitAdjMat[i,] %in% c(0,5,6))){
                if (sum(circuitAdjMat[i,]) > 0){
                  geneTypes[i] <- 2
                }
              }
            }

            configData <- NULL
            data("configData",envir = environment(), package = "sRACIPE")


            .object <- RacipeSE(
              assays = SimpleList(matrix(NA, nrow = nGenes,ncol = configData$simParams["numModels"])),
              rowData = DataFrame(circuitAdjMat),
              colData = DataFrame(matrix(NA,nrow = configData$simParams["numModels"],ncol=0)),
                                 metadata = list(
                                   annotation = filename,
                                   nInteractions = nInteraction,
                                   config = configData,
                                   geneTypes = geneTypes)
                                 )
            message("circuit file successfully loaded")
            return(.object)
          }
)

#' @export
#' @rdname sracipeGetTS
#' @aliases sracipeGetTS
setMethod(f="sracipeGetTS",
          signature="RacipeSE",
          definition=function(.object)
          {
            return(metadata(.object)$timeSeries)
          }
)

#' @rdname sracipeConfig
#' @aliases sracipeConfig
setMethod(f="sracipeConfig",
          signature="RacipeSE",
          definition=function(.object)
          {
            return(metadata(.object)$config)
          }
)

#' @rdname sracipeConfig-set
#' @aliases sracipeConfig-set
setMethod("sracipeConfig<-", "RacipeSE",
          function(.object, value) {
            metadata(.object)$config <- value
            return(.object)

          }
)

#' @rdname sracipeParams
#' @aliases sracipeParams
setMethod("sracipeParams", signature("RacipeSE"), function(.object) {
    return(as(colData(.object)[,seq_len(
        2*length(names(.object)) + 3*metadata(.object)$nInteractions),drop=FALSE], "data.frame"))
})

#' @rdname sracipeParams-set
#' @aliases sracipeParams-set
setMethod(f="sracipeParams<-",
          signature="RacipeSE",
          definition=function(.object, value)
          {
            colData(.object)[,seq_len(2*length(names(.object)) +
                     3*metadata(.object)$nInteractions)] <-
                S4Vectors::DataFrame(value)
            return(.object)
          }

)
#' @rdname sracipeIC
#' @aliases sracipeIC
setMethod(f="sracipeIC",
          signature="RacipeSE",
          definition=function(.object)
          {
            return(t(as.data.frame(colData(.object)[,(2*length(names(.object)) +
              3*metadata(.object)$nInteractions+1):(dim(colData(.object))[2])])))
          }
)
#' @rdname sracipeIC-set
#' @aliases sracipeIC-set
setMethod(f="sracipeIC<-",
          signature="RacipeSE",
          definition=function(.object, value)
          {
            value <- t(value)
            colData(.object)[,(2*length(names(.object)) +
        3*metadata(.object)$nInteractions +1):
            (dim(colData(.object))[2])] <- S4Vectors::DataFrame(value)
            return(.object)
          }
)
#' @export
#' @import SummarizedExperiment
#' @importMethodsFrom BiocGenerics annotation
#' @title  A method to get the annotation
#' @param object RacipeSE object
#' @param ... Additional arguments, for use in specific methods.
#' @examples
#'
#' data("demoCircuit")
#' rSet <- sRACIPE::sracipeSimulate(circuit = demoCircuit, numModels = 100)
#' ann <- annotation(rSet)
#' annotation(rSet) <- ann
#'
#' @return character
#'
setMethod("annotation", signature(object = "RacipeSE"), function(object,...) {
    return(metadata(object)$annotation)
    })


#' @export
#' @import SummarizedExperiment
#' @importMethodsFrom BiocGenerics annotation<-
#' @title  A method to set the circuit name or annotation
#' @param object RacipeSE object
#' @param value annotation character
#' @param ... Additional arguments, for use in specific methods.
#' @examples
#' data("demoCircuit")
#' rSet <- sRACIPE::sracipeSimulate(circuit = demoCircuit, numModels = 100)
#' ann <- annotation(rSet)
#' annotation(rSet) <- ann
#' @return A RacipeSE object
setMethod("annotation<-", signature(object = "RacipeSE"),
          function(object,...) {
              metadata(object)$annotation <- value
              return(object)
              })

#' @rdname sracipeNormalize
#' @aliases sracipeNormalize
setMethod(f="sracipeNormalize",
          signature="RacipeSE",
          definition=function(.object)
          {
            metadataTmp <- metadata(.object)
            assayDataTmp <- assays(.object)
  geneExpression <- assayDataTmp[[1]]
  geneExpression <- log2(1+geneExpression)
  means <- rowMeans(geneExpression)
  sds <-  apply(geneExpression, 1, sd)
  geneExpression <- sweep(geneExpression, 1, means, FUN = "-")
  geneExpression <- sweep(geneExpression, 1, sds, FUN = "/")
  assayDataTmp2 <- SimpleList(deterministic = geneExpression)
#  geneExpression <- data.frame(geneExpression)
#  assayDataTmp2 <- list(deterministic = geneExpression)

  tsSims <- 0
  if(!is.null(metadata(.object)$tsSimulations)){
    tsSims <- length(metadata(.object)$tsSimulations)
    tsSimulations <- assayDataTmp[2:(tsSims + 1)]
    tsSimulations <- lapply(tsSimulations,function(x) (1+x))
    tsSimulations <- lapply(tsSimulations,log2)
    #   stochasticSimulations <-
    #      lapply(stochasticSimulations, function(x) x[,is.finite(colMeans(x))])

    tsSimulations <- lapply(tsSimulations,
                            function(x) sweep(x, 1, means, FUN = "-"))
    tsSimulations <- lapply(tsSimulations,
                            function(x) sweep(x, 1, sds, FUN = "/"))
    assayDataTmp2 <- c(assayDataTmp2, tsSimulations)
  }

  stochSims <- 0
  if(!is.null(metadata(.object)$stochasticSimulations)){
    stochSims <- length(metadata(.object)$stochasticSimulations)
    stochasticSimulations <- assayDataTmp[2:(stochSims + 1)]
    stochasticSimulations <- lapply(stochasticSimulations,function(x) (1+x))
    stochasticSimulations <- lapply(stochasticSimulations,log2)
 #   stochasticSimulations <-
 #      lapply(stochasticSimulations, function(x) x[,is.finite(colMeans(x))])

    stochasticSimulations <- lapply(stochasticSimulations,
                                    function(x) sweep(x, 1, means, FUN = "-"))
    stochasticSimulations <- lapply(stochasticSimulations,
                                    function(x) sweep(x, 1, sds, FUN = "/"))
    assayDataTmp2 <- c(assayDataTmp2, stochasticSimulations)
  }



  if(!is.null(metadata(.object)$knockOutSimulations)){
    knockOutSimulations <- .object$knockOutSimulations
    koSims <- length(metadata(.object)$knockOutSimulations)
    knockOutSimulations <- assayDataTmp[(2+stochSims):(koSims +1 +stochSims)]
    for(ko in seq_len(length(knockOutSimulations))){
      simData <- knockOutSimulations[[ko]]
      tmpGene <- names(knockOutSimulations[ko])
      tmpMeans <- means
      tmpMeans[which(names(simData) == tmpGene)] <- 0
      tmpSds <- sds
      tmpSds[which(names(simData) == tmpGene)] <- 1

      simData <- log2(1+simData)
      simData[,which(names(simData) == tmpGene)] <- 0
      simData <- sweep(simData, 1, tmpMeans, FUN = "-")
      simData <- sweep(simData, 1, tmpSds, FUN = "/")
      knockOutSimulations[[ko]] <- simData
    }
    assayDataTmp2 <- c(assayDataTmp2, knockOutSimulations)

  }
  metadataTmp$normalized <- TRUE
  assays(.object) <- assayDataTmp2
  metadata(.object) <- metadataTmp
  return(.object)
}
)


#' @rdname sracipePlotCircuit
#' @aliases sracipePlotCircuit
setMethod(f="sracipePlotCircuit",
          signature="RacipeSE",
          definition=function(.object, plotToFile = TRUE)
          {
  topology <- sracipeCircuit(.object)

  if(plotToFile){


    net_file <- paste(getwd(),
                      "/network_",
                      annotation(.object),
                      ".html",
                      sep = "")
    # setwd(net_file)
  }
  node_list <-
    unique(c(topology[, 1], topology[, 2]))

  nodes <-
    data.frame(
      id = node_list,
      label = node_list,
      font.size = 50,
      value = c(rep(1, length(node_list)))
    )
  edge_col <- data.frame(c(1, 2, 3, 4, 5, 6), c("blue", "darkred", "cyan", "deeppink", "blueviolet", "darkorange"))
  arrow_type <- data.frame(c(1, 2, 3, 4, 5, 6), c("arrow", "circle", "arrow", "circle", "arrow", "circle"))
  colnames(arrow_type) <- c("type", "color")
  colnames(edge_col) <- c("type", "color")
  edges <-
    data.frame(
      from = c(topology[, 1]),
      to = c(topology[, 2]),
      arrows.to.type	= arrow_type$color[c(as.numeric(topology[, 3]))],
      width = 3,
      color = edge_col$color[c(as.numeric(topology[, 3]))]
    )

  #visNetwork(nodes, edges, height = "500px", width = "100%")  %>%
  #visEdges(arrows = "to") %>%
  #visOptions(manipulation = TRUE) %>%
  #visLayout(randomSeed = 123) %>%
  #visPhysics(solver = "forceAtlas2Based")

  network <-
    visNetwork::visNetwork(nodes, edges, height = "1000px", width = "100%") %>%
    visEdges(arrows = "to") %>%
    visOptions(manipulation = FALSE) %>%
    visLayout(randomSeed = 123) %>%
    #visNodes(scaling = list(label = list(enabled = T))) %>%
    visPhysics(solver = "forceAtlas2Based", stabilization = FALSE)
  if(plotToFile){
    visNetwork::visSave(network, file = net_file, selfcontained = FALSE)
  } else {network}

}

)


#' @rdname sracipePlotData
#' @aliases sracipePlotData
setMethod(f="sracipePlotData",
          signature="RacipeSE",
          definition=function(.object, plotToFile = TRUE, nClusters = 2,
                              heatmapPlot = TRUE,
                              pcaPlot = TRUE, umapPlot = TRUE,
                              networkPlot = TRUE,
                              clustMethod = "ward.D2", col = col,
                              distType = "euclidean",
                              assignedClusters = NULL,
                              corMethod = "spearman", ...)
          {


  if(missing(col)){
    col <-  grDevices::colorRampPalette(rev(
      RColorBrewer::brewer.pal(11, 'Spectral')))
  }
  col2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
            "#D55E00", "#CC79A7")
  p <- list()
  ts <- list()
  tsCounter <- 1

  stoch <- list()

  i=1;
  koPlot <- list()
  koPlotCounter = 1
  koHeatmap <- list()
  if(!metadata(.object)$normalized) {.object <- sracipeNormalize(.object)}

  metadataTmp <- metadata(.object)
  assayDataTmp <- assays(.object)

  if(missing(distType)){
    corCol <- stats::cor((assayDataTmp[[1]]), method = corMethod)
    distanceCol <- stats::as.dist((1 - corCol) / 2)
    corRow <- stats::cor(t(assayDataTmp[[1]]), method = corMethod)
    distanceRow <- stats::as.dist((1 - corRow) / 2)
  }
  else{
    # distType = "manhattan"
    distanceCol <- stats::dist(t(assayDataTmp[[1]]), method = distType)
    distanceRow <- stats::dist((assayDataTmp[[1]]), method = distType)
  }
  clustersCol <- stats::hclust(distanceCol, method = clustMethod)
  ddCol <- as.dendrogram(clustersCol)


  clustersRow <- stats::hclust(distanceRow, method = clustMethod)
  ddRow <- stats::as.dendrogram(clustersRow)
  if(is.null(assignedClusters)){
    clustCut <- stats::cutree(clustersCol, nClusters)
    clustColors <- col2[clustCut]
    assignedClusters <- clustCut
  }

  if(!missing(assignedClusters)){
    clustNames <- unique(assignedClusters)
    nClusters <- length(clustNames)
    clustColors <- numeric(length(assignedClusters))
    for(tmp1 in seq_len(length(clustColors))){
      clustColors[tmp1] <- which(clustNames == assignedClusters[tmp1] )
    }
    clustColors <- col2[clustColors]
  }
  names(clustColors) <- assignedClusters

  if(plotToFile){
    fileName <- paste0(annotation(.object),"_heatmap.pdf")

    pdf(fileName, onefile = TRUE)
  }
  if(heatmapPlot) {
    gplots::heatmap.2((assayDataTmp[[1]]),
                    Colv = ddCol,
                    Rowv = ddRow,
                    trace = "none",
                    col = col,
                    ColSideColors = clustColors
    )
  }


    if(plotToFile){

    dev.off()
  }
  V1 <- NULL
  V2 <- NULL
  PC1 <- NULL
    PC2 <- NULL
  if(umapPlot){
    umapGE <- umap::umap(t(assayDataTmp[[1]]))
    p[[i]] <-
      ggplot2::ggplot(data = as.data.frame(umapGE$layout)) +
      geom_point(aes(x = V1, y=V2), color = clustColors, shape = 1) +
      labs(x = "Umap1", y="Umap2") +
      theme(text = element_text(size=30),
            panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.major = element_line(color="gray", size=0.25))
    #  panel.border = element_rect(color = "black"))
    i <- i+1
  }

  if(pcaPlot){

    pca1 = summary(prcomp(t(assayDataTmp[[1]]), scale. = FALSE))
    p[[i]] <-
      ggplot2::ggplot(data = as.data.frame(pca1$x)) +
      geom_point(aes(x = PC1, y=PC2), color = clustColors, shape = 1) +
      labs(x = paste0("PC1(",100*pca1$importance[2,1],"%)"),
           y=paste0("PC2(",100*pca1$importance[2,2],"%)")) +
      theme(text = element_text(size=30),
            panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.major = element_line(color="gray", size=0.25))

    if(!is.null(metadataTmp$tsSimulations)){
      tsPca <- assayDataTmp[
        2:(1+length(metadataTmp$tsSimulations))]
      for(j in seq_len(length(tsPca))){
        tsPca[[j]] <-
          t(scale(t(tsPca[[j]]), pca1$center, pca1$scale) %*%
              pca1$rotation)

        ts[[j]] <-
          ggplot2::ggplot(data = as.data.frame(t(tsPca[[j]]))) +
          geom_point(aes(x = PC1, y=PC2), shape = 1) +
          labs(x = paste0("PC1(",100*pca1$importance[2,1],"%)"),
               y=paste0("PC2(",100*pca1$importance[2,2],"%)"),
               title = names(metadataTmp$tsSimulations)[j]) +
          theme(text = element_text(size=30),
                panel.background = element_rect(fill = "white",
                                                color = "black"),
                panel.grid.major = element_line(color="gray", size=0.25))
      }

    }

    if(!is.null(metadataTmp$stochasticSimulations)){
      stochasticPca <- assayDataTmp[
        2:(1+length(metadataTmp$stochasticSimulations))]
      for(j in seq_len(length(stochasticPca))){
        stochasticPca[[j]] <-
          t(scale(t(stochasticPca[[j]]), pca1$center, pca1$scale) %*%
              pca1$rotation)

        stoch[[j]] <-
          ggplot2::ggplot(data = as.data.frame(t(stochasticPca[[j]]))) +
          geom_point(aes(x = PC1, y=PC2), shape = 1) +
          labs(x = paste0("PC1(",100*pca1$importance[2,1],"%)"),
               y=paste0("PC2(",100*pca1$importance[2,2],"%)"),
               title = names(metadataTmp$stochasticSimulations)[j]) +
          theme(text = element_text(size=30),
                panel.background = element_rect(fill = "white",
                                                color = "black"),
                panel.grid.major = element_line(color="gray", size=0.25))
      }

    }
  }
    if(!is.null(metadataTmp$knockOutSimulations)){
      knockOutSimulations <- assayDataTmp[
        (2+length(metadataTmp$stochasticSimulations)):
          (1+length(metadataTmp$stochasticSimulations) +
             length(metadataTmp$knockOutSimulations))]
      if(plotToFile){
        fileName <- paste0(annotation(.object),"_KO_heatmaps.pdf")
        pdf(fileName, onefile = TRUE)
      }
      for(ko in seq_len(length(knockOutSimulations))){

        simData <- knockOutSimulations[[ko]]
        geneExpression <- assayDataTmp[[1]]
        tmpGene <- names(knockOutSimulations[ko])
        simData <- simData[-which(rownames(simData) == tmpGene),]
        geneExpression <- geneExpression[-which(rownames(geneExpression) ==
                                                  tmpGene),]
        if(pcaPlot){
        pcaKo <- summary(prcomp(t(geneExpression), scale. = FALSE))
        koPlot[[koPlotCounter]] <-
          ggplot2::ggplot(data = as.data.frame(pcaKo$x)) +
          geom_point(aes(x = PC1, y=PC2), color = clustColors, shape = 1) +
          labs(x = paste0("PC1(",100*pcaKo$importance[2,1],"%)"),
               y=paste0("PC2(",100*pcaKo$importance[2,2],"%)")) +
          theme(text = element_text(size=30),
                panel.background = element_rect(fill = "white",
                                                color = "black"),
                panel.grid.major = element_line(color="gray", size=0.25))
        koPlotCounter <- koPlotCounter + 1

        simDataPca <-
          t(scale(t(simData), pcaKo$center, pcaKo$scale) %*% pcaKo$rotation)
        koPlot[[koPlotCounter]] <-
          ggplot2::ggplot(data = as.data.frame(t(simDataPca))) +
          geom_point(aes(x = PC1, y=PC2), shape = 1) +
          labs(x = paste0("PC1(",100*pcaKo$importance[2,1],"%)"),
               y=paste0("PC2(",100*pcaKo$importance[2,2],"%)"),
               title = names(knockOutSimulations[ko])) +
          theme(text = element_text(size=30),
                panel.background = element_rect(fill = "white",
                                                color = "black"),
                panel.grid.major = element_line(color="gray", size=0.25))
        koPlotCounter <- koPlotCounter + 1
        }
        if(heatmapPlot){

          gplots::heatmap.2(
            geneExpression, trace = "none", col = col, main = "WT",
            hclustfun = function(x) hclust(x, method = 'ward.D2'),
            distfun=function(x) as.dist((1-cor(t(x), method = "spear"))/2)
          )
          gplots::heatmap.2(
            simData, trace = "none", col = col,
            main = names(knockOutSimulations[ko]),
            hclustfun = function(x) hclust(x, method = 'ward.D2'),
            distfun=function(x) as.dist((1-cor(t(x), method = "spear"))/2)
          )
        }
      }
      if(heatmapPlot) dev.off()
    }

  if(plotToFile){
    fileName <- paste0(annotation(.object),"_plots.pdf")
    pdf(fileName, onefile = TRUE)
  }

  for (i in seq(length(p))) {
     gridExtra::grid.arrange(p[[i]])
    # do.call("grid.arrange", p[[i]])
  }

  if(plotToFile){
    message("Plots saved as pdf files in the working directory.")
    dev.off()
  }
    if(!is.null(metadataTmp$tsSimulations)){
      if(plotToFile){
        fileName <- paste0(annotation(.object),"_tsPlots.pdf")
        pdf(fileName, onefile = TRUE)
      }
      for (i in seq_along(ts)) {
        gridExtra::grid.arrange(ts[[i]])
        #   do.call("grid.arrange", p[[i]])
      }
      if(plotToFile){
        dev.off()
      }
    }

  if(!is.null(metadataTmp$stochasticSimulations)){
    if(plotToFile){
      fileName <- paste0(annotation(.object),"_stochasticPlots.pdf")
      pdf(fileName, onefile = TRUE)
    }
    for (i in seq_along(stoch)) {
      gridExtra::grid.arrange(stoch[[i]])
      #   do.call("grid.arrange", p[[i]])
    }
    if(plotToFile){
      dev.off()
    }
  }

  if(!is.null(metadataTmp$knockOutSimulations)){
    if(plotToFile){
      fileName <- paste0(annotation(.object),"_KO_Plots.pdf")
      pdf(fileName, onefile = TRUE)
    }
    for (i in seq(length(koPlot)/2)) {
      gridExtra::grid.arrange(koPlot[[2*i-1]],koPlot[[2*i]], nrow = 2)
      #   do.call("grid.arrange", p[[i]])
    }
    if(plotToFile){
      dev.off()
    }
  }


  if(networkPlot){
    sracipePlotCircuit(.object, plotToFile = plotToFile)
  }
  if(umapPlot)metadataTmp$umap <- umapGE
  if(pcaPlot) metadataTmp$pca <- pca1
  metadataTmp$assignedClusters <- assignedClusters
  metadata(.object) <- metadataTmp
  return(.object)
}

)

#' @rdname sracipePlotParamBifur
#' @aliases sracipePlotParamBifur
setMethod(f="sracipePlotParamBifur",
          signature="RacipeSE",
          definition=function(.object, paramName, data = NULL,
                              paramValue = NULL, assignedClusters = NULL,
                              plotToFile = FALSE){

  if(missing(paramValue)){
    paramValue <- as.matrix(sracipeParams(.object))
    paramValue <- paramValue[,paramName]
  }
  if(missing(data)){
    data = t(assays(.object)[[1]])
  }
  paramValue <- rep(paramValue, ncol(data))
  data <- reshape2::melt(data)[,2:3]
  Expression <- NULL
  Gene <- NULL
  colnames(data) <- c("Gene", "Expression")
  if(missing(assignedClusters)){
    if(is.null(metadata(.object)$assignedClusters)){
      assignedClusters <- data$Gene
    } else {
      assignedClusters <- metadata(.object)$assignedClusters
      assignedClusters <- as.factor(rep(assignedClusters, ncol(data)))

    }
  } else {
    assignedClusters <- as.factor(rep(assignedClusters, ncol(data)))
  }
  if(plotToFile){
    fileName <- paste0(annotation(.object),"_",
                       paramName,"_BifurPlot.pdf")
    pdf(fileName, onefile = TRUE)
  }
  Cluster <- assignedClusters
  p <- ggplot(data, aes(x = paramValue, y=Expression, color = Gene,
                        shape = Cluster) ) +
    geom_point() +
    theme_bw() +
    labs(x=paramName) +
    theme(text = element_text(size=20))
  gridExtra::grid.arrange(p)

  if(plotToFile){
    dev.off()
  }
  return()
          }

)

#' @rdname sracipeOverExp
#' @aliases sracipeOverExp
setMethod(f="sracipeOverExp",
          signature="RacipeSE",
          definition=function(.object, overProduction = 10,
                              nClusters = 2,
                              clusterOfInterest = 2,
                              plotFilename = NULL,
                              plotHeatmap = TRUE,
                              plotBarPlot = TRUE,
                              clusterCut = NULL,
                              plotToFile = FALSE){

            if(!metadata(.object)$normalized) {.object <- sracipeNormalize(.object)}

              dataSimulation <- assays(.object)[[1]]
              geneNames <- names(.object)
              params <- sracipeParams(.object)
              params <- params[, seq_len(nrow(dataSimulation))]

            if (missing(plotFilename))
              plotFilename <- annotation(.object)
            filename <-
              (paste(plotFilename, "_overExpr.pdf", sep = ""))

            #library(htmlwidgets)
            #library(d3heatmap)
            rf <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
            plot_color <- rf(32)
            plotColor2 <- c("#000000", "#E69F00", "#56B4E9",
                            "#009E73", "#F0E442", "#0072B2",
                             "#D55E00", "#CC79A7")
            dataSimulation <- t(dataSimulation)

            # Find cluster assignments
            if(missing(clusterCut)){
              ref_cor <- cor((dataSimulation), method = "spearman")
              distance <- as.dist((1 - ref_cor) / 2)
              clusters <- hclust(distance, method = "ward.D2")
              clusterCut <- cutree(clusters, nClusters)
            }

            # List to hold the cluster percentages for
            # each gene knockdown and REF-reference or no knockdown
            knockdownsName <- c("WT", geneNames)
            perturbationKd <- as.list(knockdownsName)
            names(perturbationKd) <- knockdownsName
            # No knockdown cluster ratios
            perturbationKd[[1]] <-
              table((clusterCut)) / sum(table((clusterCut)))
            clusterNames <- as.character(seq_len(nClusters))

            #knockdown cluster ratios
            for (i in seq_len(length(knockdownsName) - 1)) {
              perturbationKd[[i + 1]] <-
                table(clusterCut[params[, i] > (
                  max(params[, i]) - overProduction * 0.01 *
                    (max(params[, i]) - min(params[, i])))]) /
                sum(table(clusterCut[
                  params[, i] > (max(params[, i]) - overProduction * 0.01 *
                                   (max(params[, i])-min(params[, i])))]))

              # check if any cluster is missing after knockdown
              if (nrow(perturbationKd[[i + 1]]) < nClusters) {
                missingClusters <-
                  clusterNames[!((clusterNames %in%
                                    names(perturbationKd[[i + 1]])))]
                # add 1 model of each missing cluster
                table(c(clusterCut[params[, i] < (
                  min(params[, i]) + overProduction * 0.01 * (
                    max(params[, i]) -min(params[, i])))], missingClusters)) /
                  sum(table(c(clusterCut[params[, i] < (
                    min(params[, i]) + overProduction * 0.01 * (
                      max(params[, i])-min(params[, i])))], missingClusters)))
              }
            }

            if (plotBarPlot) {
              # Restructure dataframe for plotting
              meltPKd <- reshape2::melt(perturbationKd)
              meltPKd$Var1 <- factor(meltPKd$Var1)
              L1 <- NULL
              value <- NULL
              Cluster <- NULL
              names(meltPKd) <- c("Cluster",  "value", "L1")
              # convert ratio to percentages
              meltPKd$value <- 100 * meltPKd$value
              # sort transcription factors based on increasing ratio of
              # last cluster
              if(missing(clusterOfInterest)) clusterOfInterest <- nClusters
              meltPKd$L1 <-
                factor(meltPKd$L1, levels = meltPKd$L1[(
                  order(meltPKd$value[meltPKd$Cluster ==
                                        as.character(clusterOfInterest)])) *
                                                         (nClusters)])
              # plot the histogram
              if(plotToFile){
                pdf(paste(plotFilename, "_overExprBarplot.pdf",
                          sep = ""))
              }
              p <- ggplot2::ggplot(data = meltPKd, aes(x = L1, y = value,
                                                       fill = Cluster)) +
                geom_bar(stat = "identity") +
                coord_flip() +
                labs(title = "TF Over Expression Analysis") +
                xlab("Transcription Factor") +
                ylab("Cluster Percentage") +
                scale_fill_manual(values = plotColor2) +
                theme(axis.text.x = element_text(angle = 0, hjust = 1),
                      text = element_text(size = 12))
              gridExtra::grid.arrange(p)
              if(plotToFile){
                dev.off()
              }
            }

            if (plotHeatmap) {
              if(plotToFile){
                pdf(paste(plotFilename, "_heatmap.pdf", sep = ""))
              }
              gplots::heatmap.2((dataSimulation),
                                col = plot_color,
                                hclustfun = function(x)
                                  hclust(x, method = 'ward.D2'),
                                distfun = function(x)
                                  as.dist((1 - cor(t(
                                    x
                                  ), method = "spear")) / 2),
                                trace = "none",
                                ColSideColors = plotColor2[clusterCut]
              )
              if(plotToFile){
                dev.off()
              }
            }
            return(perturbationKd)
          }
)
#' @export
#' @rdname sracipeKnockDown
#' @aliases sracipeKnockDown
setMethod(f="sracipeKnockDown",
          signature="RacipeSE",
          definition=function(.object, reduceProduction = 10,
                              nClusters = 2,
                              clusterOfInterest = 2,
                              plotFilename = NULL,
                              plotHeatmap = TRUE,
                              plotBarPlot = TRUE,
                              clusterCut = NULL,
                              plotToFile = FALSE){
            if(!metadata(.object)$normalized) {.object <- sracipeNormalize(.object)}
            dataSimulation <- assays(.object)[[1]]
            geneNames <- names(.object)
            params <- sracipeParams(.object)
            params <- params[, seq_len(nrow(dataSimulation))]

            if (missing(plotFilename))
              plotFilename <- annotation(.object)

            filename <-
              (paste(plotFilename, "_knockdown.pdf", sep = ""))

            #library(htmlwidgets)
            #library(d3heatmap)
            rf <- grDevices::colorRampPalette(rev(
              RColorBrewer::brewer.pal(11, 'Spectral')))
            plot_color <- rf(32)
            plot_color2 <- c("#000000", "#E69F00", "#56B4E9",
                             "#009E73", "#F0E442", "#0072B2",
                             "#D55E00", "#CC79A7")


            # Find cluster assignments
            if(missing(clusterCut)){
              ref_cor <- cor((dataSimulation), method = "spearman")
              distance <- as.dist((1 - ref_cor) / 2)
              clusters <- hclust(distance, method = "ward.D2")
              clusterCut <- cutree(clusters, nClusters)
            }

            # List to hold the cluster percentages for each gene
            #knockdown and REF-reference or no knockdown
            knockdownsName <- c("WT", geneNames)
            perturbationKd <- as.list(knockdownsName)
            names(perturbationKd) <- knockdownsName
            # No knockdown cluster ratios
            perturbationKd[[1]] <-
              table((clusterCut)) / sum(table((clusterCut)))
            clusterNames <- as.character(seq_len(nClusters))

            #knockdown cluster ratios
            for (i in seq_len(length(knockdownsName) - 1)) {
              perturbationKd[[i + 1]] <-
                table(clusterCut[params[, i] < (
                  min(params[, i]) + reduceProduction * 0.01 * (
                    max(params[, i]) - min(params[, i])))]) /
                sum(table(clusterCut[params[, i] < (
                  min(params[, i]) + reduceProduction * 0.01 * (
                    max(params[, i])-min(params[, i])))]))

              # check if any cluster is missing after knockdown
              if (nrow(perturbationKd[[i + 1]]) < nClusters) {
                missingClusters <-
                  clusterNames[!((clusterNames %in% names(
                    perturbationKd[[i + 1]])))]
                # add 1 model of each missing cluster
                table(c(clusterCut[params[, i] < (
                  min(params[, i]) + reduceProduction * 0.01 * (
                    max(params[, i])-min(params[, i])))], missingClusters)) /
                  sum(table(c(clusterCut[params[, i] < (
                    min(params[, i]) + reduceProduction * 0.01 * (
                      max(params[, i])-min(params[, i])))], missingClusters)))
              }
            }

            if (plotBarPlot) {
              # Restructure dataframe for plotting
              meltPKd <- reshape2::melt(perturbationKd)
              meltPKd$Var1 <- factor(meltPKd$Var1)
              L1 <- NULL
              value <- NULL
              Cluster <- NULL
              names(meltPKd) <- c("Cluster",  "value", "L1")
              # convert ratio to percentages
              meltPKd$value <- 100 * meltPKd$value
              # sort transcription factors based on
              # increasing ratio of last cluster
              if(missing(clusterOfInterest)) clusterOfInterest <- nClusters
              meltPKd$L1 <-
                factor(meltPKd$L1, levels = meltPKd$L1[(
                  order(meltPKd$value[meltPKd$Cluster ==
                                        as.character(clusterOfInterest)])) *
                                                         (nClusters)])
              # plot the histogram

              p <- ggplot2::ggplot(data = meltPKd, aes(x = L1, y = value,
                                                       fill = Cluster)) +
                geom_bar(stat = "identity") +
                coord_flip() +
                labs(title = "TF Knockout Analysis") +
                xlab("Transcription Factor") +
                ylab("Cluster Percentage") +
                scale_fill_manual(values = plot_color2) +
                theme(axis.text.x = element_text(angle = 0, hjust = 1),
                      text = element_text(size = 12))
              if(plotToFile){
                fileName <- paste(plotFilename, "_kdBarplot.pdf",
                                  sep = "")
                pdf(fileName, onefile = TRUE)
              }
              gridExtra::grid.arrange(p)
              if(plotToFile){
                dev.off()
              }


            }

            if (plotHeatmap) {
              if(plotToFile){
                pdf(paste(plotFilename, "_heatmap.pdf", sep = ""))
              }
              gplots::heatmap.2((dataSimulation),
                                col = plot_color,
                                hclustfun = function(x)
                                  hclust(x, method = 'ward.D2'),
                                distfun = function(x)
                                  as.dist((1 - cor(t(
                                    x
                                  ), method = "spear")) / 2),
                                trace = "none",
                                ColSideColors = plot_color2[clusterCut]
                )
              if(plotToFile){
                dev.off()
              }
            }
            return(perturbationKd)
          }


)

#' @export
#' @rdname sracipeConvergeDist
#' @aliases sracipeConvergeDist
setMethod(f="sracipeConvergeDist",
          signature="RacipeSE",
          definition=function(.object, plotToFile = FALSE)
          {
            metadataTmp <- metadata(.object)
            #Checks if modelConvergence data exists in the object metadata
            dataExists <- "modelConvergence" %in% names(metadataTmp)
            if(!dataExists){
              message("Cannot plot without convergence data")
              return(.object)
            }
            convergenceData <- metadataTmp$modelConvergence
            configuration <- metadataTmp$config
            numModels <- configuration$simParams["numModels"]
            nIC <- configuration$simParams["nIC"]
            numConvergenceTests <- configuration$simParams["numConvergenceTests"]
            numExprx <- numModels #Check RacipeSE() constructor for why this is true
            lc <- configuration$options["limitcycles"]

            #Initialize proportions
            convergedProportions <- numeric(numConvergenceTests)

            if(plotToFile){
              fileName <- paste0(annotation(.object),"_ConvergDist.pdf")
              pdf(fileName) #Opens graphics object to store file in
            }

            #Getting rid of non-converged models
            convergedICs <- convergenceData[convergenceData[, 1] == 1, ]
            testScores <- convergedICs[,2]

            #Removing limit cycles from consideration
            if(lc){
              numLCICs <- length(which(convergenceData[,1] == 2))
              numExprx <- numExprx - numLCICs
            }

            #Removing models with NaN values
            if(!all(convergenceData[, 1] != 3)){
              numNaN <- sum(convergenceData[, 1] == 3)
              numExprx <- numExprx - numNaN
            }

            for (i in 1:numConvergenceTests){
              convergedProportions[i] <- sum(testScores <= i) / numExprx
            }

            title = paste0("Ratio of Stable ", annotation(.object), " Models over number of convergence tests")
            plot(seq(1,numConvergenceTests), convergedProportions, type="l", col="blue",
                 xlab="# Convergence Tests", ylab = "% Converged Models",
                 main = title)

            if(plotToFile){
              message("Plot saved as pdf files in the working directory.")
              dev.off() #closes graphics object and send it to working directory
            }

            metadataTmp$stableProportion <- convergedProportions[numConvergenceTests]
            if(convergedProportions[numConvergenceTests] > 0.99){
              ninetyNineIndex <- which(convergedProportions > 0.99)[1]
              metadataTmp$ninetyNineConvergedNum <- ninetyNineIndex
            }
            metadata(.object) <- metadataTmp

            return(.object)
          }

)

#' @export
#' @rdname sracipeCombineRacipeSE
#' @aliases sracipeCombineRacipeSE
setMethod(f="sracipeCombineRacipeSE",
          signature = "list",
          definition = function(.object){
            #Validate sameness of objects using first element
            validationConfig <- sracipeConfig(.object[[1]])
            validCircuit <- sracipeCircuit(.object[[1]])

            #Can't use full config for validation b/c threshold vals are expected
            #to be different, so use param vectors instead
            validSim <- validationConfig$simParams
            validStoch <- validationConfig$stochParams
            validHyper <- validationConfig$hyperParams
            validOptions <- validationConfig$options
            validLC <- validationConfig$LCParams

            nIC <- validSim["nIC"]
            numModels <- validSim["numModels"] / nIC

            #Create lists for storing values
            statesList <- list()
            paramsList <- list()
            icList <- list()
            if(validOptions["convergTesting"]){
              convergList <- list()
              if(nIC > 1){uniqueCountList <- list()}
              if(validOptions["limitcycles"]){
                LCList <- list()
                LCNum <- 0
              }
            }

            for(racipeObj in .object){
              #validation
              objConfig <- sracipeConfig(racipeObj)
              if(!all(sracipeCircuit(racipeObj) == validCircuit,
                      objConfig$simParams == validSim, objConfig$stochParams == validStoch,
                      objConfig$hyperParams == validHyper, objConfig$options == validOptions,
                      objConfig$LCParams == validLC)){
                message("One of the provided RacipeSE objects has different params or topo
                        from the first object")
                return()
              }

              statesList <- c(statesList, assay(racipeObj))
              paramsList <- c(paramsList, sracipeParams(racipeObj))
              icList <- c(icList, sracipeIC(racipeObj))
              if(validOptions["convergTesting"]){
                objMetadata <- metadata(racipeObj)
                convergList <- c(convergList, objMetadata$modelConvergence)
                if(nIC > 1){
                  uniqueCountList <- c(uniqueCountList, objMetadata$uniqueStateCounts)
                }
                if(validOptions["limitcycles"]){
                  if(!("LCData" %in% names(objMetadata))){
                    next
                  }
                  objLC <- objMetadata$LCData

                  #Fixing modelCount info in LCData for each object
                  objIdx <- which(.object == racipeObj)
                  if(objIdx > 1){
                    objLC[,1] <- objLC[,1]*(objIdx-1)*numModels
                  }
                  LCList <- c(LCList, objLC)
                  LCNum <- LCNum + objMetadata["totalNumofLCs"]
                }
              }
            }

            #Gluing things together
            print("here1")
            combinedStates <- do.call(cbind, statesList)
            combinedParams <- do.call(rbind, paramsList)
            print("here2")
            combinedICs <- do.call(cbind, icList)

            print("here3")
            col <- cbind(combinedParams, t(combinedICs))

            metadataTmp <- metadata(.object[[1]])
            if(validOptions["convergTesting"]){
              metadataTmp$modelConvergence <- do.call(rbind, convergList)
              if(nIC > 1){
                combinedUniqueCounts <- do.call(rbind, uniqueCountList)
                combinedUniqueCounts[,1] <- seq_len(numModels*(length(.object)))
                metadataTmp$uniqueStateCounts <- combinedUniqueCounts
              }
              if(validOptions["limitcycles"]){
                metadataTmp$LCData <- do.call(rbind, LCList)
              }
            }


            rSet <- RacipeSE(rowData = rowData(.object[[1]]), colData = col,
                             assays = combinedStates, metadata = metadataTmp)

            return(rSet)

          }
)
