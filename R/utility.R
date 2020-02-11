getInteractions <- function(value) {
  interaction = 1
  nInteractions <- sum(value > 0)
  circuitInteractions <-
    data.frame(matrix(nrow = nInteractions, ncol = 3))
  geneNames <- rownames(value)
  for (i in seq_len(nrow(value))) {
    for (j in seq_len(nrow(value))) {
      if (as.integer(value[i, j]) > 0) {
        circuitInteractions[interaction, 2] <- geneNames[i]
        circuitInteractions[interaction, 1] <- geneNames[j]
        circuitInteractions[interaction, 3] <-
          as.integer(value[i, j])
        interaction <- interaction + 1
      }
    }
  }
  names(circuitInteractions) <- c("Source", "Target", "Type")
  return(circuitInteractions)
}

#' @export
#' @title Generate parameter names for a circuit
#' @param circuit RacipeSE object or topology as data.frame or filename
#' @examples
#' rSet <- RacipeSE()
#' data("demoCircuit")
#' sracipeCircuit(rSet) <- demoCircuit
#' paramNames <- sRACIPE::sracipeGenParamNames(rSet)
#'
#' @return list
sracipeGenParamNames <- function(circuit = "inputs/test.tpo") {
  if (methods::is(circuit, "RacipeSE")) {
    rSet <- circuit
  }  else {
    rSet <- RacipeSE()
    sracipeCircuit(rSet) <- circuit
  }
  geneInteraction <- as.matrix(rowData(rSet))
  geneNames <- names(rSet)
  
  paramList <- list()
  tmp <- lapply(geneNames,  function(x)
    paste("G_", x, sep = ""))
  paramList <- append(paramList, tmp)
  tmp <- lapply(geneNames,  function(x)
    paste("K_", x, sep = ""))
  paramList <- append(paramList, tmp)
  tmp <- list()
  tmp2 <- list()
  tmp3 <- list()
  for (gene1 in seq_along(geneNames))
  {
    for (gene2 in seq_along(geneNames))
    {
      if (geneInteraction[gene1, gene2] > 0) {
        tmp <- append(tmp,
                      paste("TH_", geneNames[[gene2]], "_",
                            geneNames[[gene1]], sep = ""))
        tmp2 <- append(tmp2,
                       paste("N_", geneNames[[gene2]], "_",
                             geneNames[[gene1]], sep = ""))
        tmp3 <- append(tmp3,
                       paste("FC_", geneNames[[gene2]], "_",
                             geneNames[[gene1]], sep = ""))
      }
    }
  }
  paramList <- do.call(c, list(paramList, tmp, tmp2, tmp3))
  
  return(as.character(paramList))
}

.loadNetworkFile <- function( networkFile = "inputs/test.net") {
  #' Loads the network/topology file.
  #'
  #' The network file should contain three columns with headers,
  #' "Source" "Target" "Type"
  #' Here "Source" and "Target" are the names of the genes and "Type" refers to
  #' the regulation, "1" if source activates target and "2" if source inhibits
  #' target.
  #' @param networkFile Network file name
  #' @return Network as a dataframe
  #'
  #'
  
  if(missing(networkFile)){
    stop("Please specify the network file!")
  }
  
  if(file.exists(networkFile)){
    networkTable <- read.table(networkFile, header = TRUE,
                               stringsAsFactors = FALSE)
    colnames(networkTable) <- c("Source","Target","Type")
    return(networkTable)
  } else {
    stop("Network file not found!")
  }
  
}


#' @export
#' @title Density Plot
#' @description Plot the density of points as an image alongwith histograms on
#' the sides.
#' @param plotData Dataframe containing the data.
#' @param binCount (optional) Integer. Default 40. The number of bins to be used for
#' dividing the data along an axis.
#' @param plotColor (optional) The color palette.
#' @param ... any additional arguments
#' @return plot
#' @section Related Functions:
#'
#' \code{\link{sracipeSimulate}},  \code{\link{sracipeKnockDown}},
#' \code{\link{sracipeOverExp}},  \code{\link{sracipePlotData}},
#' \code{\link{sracipeHeatmapSimilarity}}

densityPlot = function(plotData, binCount=40, plotColor=NULL, ...) {
  colnames(plotData) <- c("x", "y")
  if(is.null(plotColor)){
    rf <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'Spectral')))
    plotColor <- rf(32)
  }
  h1 <- hist(plotData$x, breaks = binCount, plot = FALSE)
  h2 <- hist(plotData$y, breaks = binCount, plot = FALSE)
  top <- max(h1$counts, h2$counts)
  k <- MASS::kde2d(plotData$x, plotData$y, n = binCount)
  
  # margins
  oldpar <- par()
  par(mar = c(3, 3, 1, 1))
  layout(matrix(c(2, 0, 1, 3), 2, 2, byrow = TRUE), c(3, 1), c(1, 3))
  image(k, col = plotColor) #plot the image
  par(mar = c(0, 2, 1, 0))
  barplot(
    h1$counts,
    axes = FALSE,
    ylim = c(0, top),
    space = 0,
    col = 'red'
  )
  par(mar = c(2, 0, 0.5, 1))
  barplot(
    h2$counts,
    axes = FALSE,
    xlim = c(0, top),
    space = 0,
    col = 'red',
    horiz = TRUE
  )
}
