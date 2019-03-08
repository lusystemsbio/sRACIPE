# This file contains deprecated functions which will be removed in a future
# version

#' @export
#' @title Load Data (Deprecated).
#' @description Use \code{\link{simulateGRC}}
#'
#' \code{\link{simulateGRC}},  \code{\link{knockdownAnalysis}},
#' \code{\link{overExprAnalysis}},  \code{\link{plotData}},
#' \code{\link{calcHeatmapSimilarity}}
#'
load_data = function(data_simulation = data_simulation,
                     topology_df = topology) {
  working_directory <- getwd()
  name_genes <-
    read.table(
      paste(
        working_directory,
        "/results/gene_interaction_topology_",
        topology$filename,
        ".txt",
        sep = ""
      ),
      header = TRUE,
      stringsAsFactors = FALSE
    )
  name_genes <- t(as.matrix(name_genes))
  data_simulation <- log2(data_simulation)
  data_simulation <- scale(data_simulation)
  name_models <- seq(1:nrow(data_simulation))
  row.names(data_simulation) <- name_models
  colnames(data_simulation) <- name_genes
  return(data_simulation)
}



#' @export
#' @title Plot Stochastic Data (Deprecated)
#' @description Use \code{\link{plotData}} instead.
#'
#' \code{\link{simulateGRC}},  \code{\link{knockdownAnalysis}},
#' \code{\link{overExprAnalysis}},  \code{\link{plotData}},
#' \code{\link{calcHeatmapSimilarity}}
plot_data_stochastic = function(output_file,
                                plot_filename = filename,
                                topology_df = topology,
                                config = configuration,
                                bin_count = 40) {
  message("Plotting the stochastic results")
  rf <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
  plot_color <- rf(32)

  if (missing(bin_count))
    bin_count <- 40
  if (missing(plot_filename))
    output_filename <- "plots"

  working_directory <- getwd()
  name_genes <-
    read.table(
      paste(
        working_directory,
        "/results/gene_interaction_topology_",
        topology$filename,
        ".txt",
        sep = ""
      ),
      header = T,
      stringsAsFactors = F
    )

  data_simulation_all <-
    as.data.frame(read.table(output_file, header = F))
  col_start <- topology$number_gene * (configuration$NOISE_LEVELS - 1) +
    1
  col_end <- topology$number_gene * (configuration$NOISE_LEVELS)
  data_simulation <-
    as.data.frame(data_simulation_all[, col_start:col_end])
  data_simulation <- log2(data_simulation)
  data_simulation <-
    data_simulation[is.finite(rowMeans(data_simulation)),]

  means <- colMeans(data_simulation)
  sds <- apply(data_simulation, 2, sd)
  data_simulation <- sweep(data_simulation, 2, means, FUN = "-")
  data_simulation <- sweep(data_simulation, 2, sds, FUN = "/")

  name_models <- seq(1:nrow(data_simulation))
  name_genes <- t(as.matrix(name_genes))
  row.names(data_simulation) <- name_models
  colnames(data_simulation) <- name_genes


  pca = prcomp(data_simulation, scale. = FALSE, center = FALSE)
  #pca_data <- data.frame(x=pca$x[,1],y=pca$x[,2])

  pdf(paste("results/", plot_filename, "_stochastic_pca.pdf", sep = ""))

  for (i in 1:configuration$NOISE_LEVELS)
  {
    col_start <- (topology$number_gene) * (configuration$NOISE_LEVELS - i) +
      1
    col_end <-
      (topology$number_gene) * (configuration$NOISE_LEVELS - i + 1)
    data_simulation <-
      as.data.frame(data_simulation_all[, col_start:col_end])
    data_simulation <- log2(data_simulation)
    data_simulation <-
      data_simulation[is.finite(rowMeans(data_simulation)),]

    data_simulation <- sweep(data_simulation, 2, means, FUN = "-")
    data_simulation <- sweep(data_simulation, 2, sds, FUN = "/")

    name_models <- seq(1:nrow(data_simulation))
    row.names(data_simulation) <- name_models
    colnames(data_simulation) <- name_genes

    pca_data <-
      scale(data_simulation, pca$center, pca$scale) %*% pca$rotation
    pca_data <- data.frame(x = pca_data[, 1], y = pca_data[, 2])
    density_plot(pca_data, bin_count, plot_color)
  }


  dev.off()

  message("Plots in the pdf file in the results folder.")

}

#' @export
#' @title Plot knockout Data (Deprecated)
#' @description Use \code{\link{plotData}} instead.
#'
#' \code{\link{simulateGRC}},  \code{\link{knockdownAnalysis}},
#' \code{\link{overExprAnalysis}},  \code{\link{plotData}},
#' \code{\link{calcHeatmapSimilarity}}
#'
plot_data_knockout_single = function(output_file,
                                     output_file_knockout,
                                     plot_filename = filename,
                                     topology_df = topology,
                                     KNOCKOUT = NA_character_,
                                     config = configuration,
                                     bin_count = 40) {
  message("Plotting the single knockout results")
  rf <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
  plot_color <- rf(32)

  if (missing(bin_count))
    bin_count <- 40
  if (missing(plot_filename))
    output_filename <- "plots"

  working_directory <- getwd()
  name_genes <-
    read.table(
      paste(
        working_directory,
        "/results/gene_interaction_topology_",
        topology$filename,
        ".txt",
        sep = ""
      ),
      header = T,
      stringsAsFactors = F
    )
  name_genes <- t(as.matrix(name_genes))
  data_simulation <- read.table(output_file, header = F)
  name_models <- seq(1:nrow(data_simulation))
  row.names(data_simulation) <- name_models
  colnames(data_simulation) <- name_genes
  knockout_number <- as.integer(which(name_genes == KNOCKOUT))

  data_simulation <- log2(data_simulation)


  data_simulation <- data_simulation[, -knockout_number]

  col.means <- colMeans(data_simulation)
  col.sds <- apply(data_simulation, 2, sd)
  data_simulation <- sweep(data_simulation, 2, col.means, FUN = "-")
  data_simulation <- sweep(data_simulation, 2, col.sds, FUN = "/")


  data_simulation_knockout <-
    read.table(output_file_knockout, header = F)
  name_models <- seq(1:nrow(data_simulation_knockout))
  row.names(data_simulation_knockout) <- name_models
  colnames(data_simulation_knockout) <- name_genes

  data_simulation_knockout <-
    data_simulation_knockout[, -knockout_number]

  data_simulation_knockout <- log2(data_simulation_knockout)
  data_simulation_knockout <-
    sweep(data_simulation_knockout, 2, col.means, FUN = "-")
  data_simulation_knockout <-
    sweep(data_simulation_knockout, 2, col.sds, FUN = "/")


  pca = prcomp(data_simulation, scale. = FALSE, center = FALSE)
  pca_data <- data.frame(x = pca$x[, 1], y = pca$x[, 2])

  pdf(paste("results/", plot_filename, KNOCKOUT, "_knockout.pdf", sep = ""))
  density_plot(pca_data, bin_count, plot_color)

  pca_data_knockout <-
    scale(data_simulation_knockout, pca$center, pca$scale) %*% pca$rotation
  pca_data_knockout <-
    data.frame(x = pca_data_knockout[, 1], y = pca_data_knockout[, 2])
  density_plot(pca_data_knockout, bin_count, plot_color)
  #col_count <- as.integer(65536/dim(data_simulation)[2])
  #heatmap_data <- t(data_simulation[1:col_count,])
  #dist <- as.dist((1-cor(t(heatmap_data), method = "spear"))/2)
  #cluster <- hclust(dist,method = 'ward.D2')
  #dendrogram <- as.dendrogram(cluster)
  full <-
    heatmap.2(
      t(data_simulation),
      col = plot_color,
      hclustfun = function(x)
        hclust(x, method = 'ward.D2'),
      distfun = function(x)
        as.dist((1 - cor(t(
          x
        ), method = "spear")) / 2),
      trace = "none"
    )

  heatmap.2(
    t(data_simulation_knockout[, rev(full$rowInd)]),
    Rowv = NA,
    dendrogram = "column",
    breaks = full$breaks,
    col = plot_color,
    hclustfun = function(x)
      hclust(x, method = 'ward.D2'),
    distfun = function(x)
      as.dist((1 - cor(t(
        x
      ), method = "spear")) / 2),
    trace = "none"
  )

  # heatmap.2(t(data_simulation_knockout), col=plot_color, hclustfun = function(x) hclust(x,method = 'ward.D2'), distfun=function(x) as.dist((1-cor(t(x), method = "spear"))/2), trace = "none")

  dev.off()

  message("Plots in the pdf file in the results folder.")

}
#' @export
#' @title Plot knockout data for all genes (Deprecated)
#' @description Use \code{\link{plotData}} instead.
#'
#' \code{\link{simulateGRC}},  \code{\link{knockdownAnalysis}},
#' \code{\link{overExprAnalysis}},  \code{\link{plotData}},
#' \code{\link{calcHeatmapSimilarity}}
plot_data_knockout_all = function(output_file,
                                  output_file_KOpre,
                                  output_file_KOsuf,
                                  tfs_knockout,
                                  plot_filename = filename,
                                  topology_df = topology,
                                  config = configuration,
                                  bin_count = 40) {
  message("Plotting the single knockout results")
  rf <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
  plot_color <- rf(32)

  if (missing(bin_count))
    bin_count <- 40
  if (missing(plot_filename))
    output_filename <- "plots"

  working_directory <- getwd()
  name_genes <-
    read.table(
      paste(
        working_directory,
        "/results/gene_interaction_topology_",
        topology$filename,
        ".txt",
        sep = ""
      ),
      header = T,
      stringsAsFactors = F
    )
  name_genes <- t(as.matrix(name_genes))
  data_simulation <- read.table(output_file, header = F)
  name_models <- seq(1:nrow(data_simulation))
  row.names(data_simulation) <- name_models
  colnames(data_simulation) <- name_genes
  data_simulation <- log2(data_simulation)
  data_simulation1 <- data_simulation
  pdf(paste("results/", plot_filename, "_knockout.pdf", sep = ""))

  for (tf.current in 1:length(tfs_knockout)) {
    KNOCKOUT <- tfs_knockout[tf.current]
    knockout_number <- as.integer(which(name_genes == KNOCKOUT))
    output_file_knockout <-
      paste0(output_file_KOpre, KNOCKOUT, output_file_KOsuf)
    data_simulation <- data_simulation1[, -knockout_number]

    col.means <- colMeans(data_simulation)
    col.sds <- apply(data_simulation, 2, sd)
    data_simulation <- sweep(data_simulation, 2, col.means, FUN = "-")
    data_simulation <- sweep(data_simulation, 2, col.sds, FUN = "/")

    data_simulation_knockout <-
      read.table(output_file_knockout, header = F)
    name_models <- seq(1:nrow(data_simulation_knockout))
    row.names(data_simulation_knockout) <- name_models
    colnames(data_simulation_knockout) <- name_genes

    data_simulation_knockout <-
      data_simulation_knockout[, -knockout_number]

    data_simulation_knockout <- log2(data_simulation_knockout)
    data_simulation_knockout <-
      sweep(data_simulation_knockout, 2, col.means,
            FUN = "-")
    data_simulation_knockout <-
      sweep(data_simulation_knockout, 2, col.sds,
            FUN = "/")


    pca = prcomp(data_simulation, scale. = FALSE, center = FALSE)
    pca_data <- data.frame(x = pca$x[, 1], y = pca$x[, 2])

    par(mfrow = c(2, 1))
    #par(mfrow=c(1,2), pty="s")

    #layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
    h1 <- hist(pca_data$x, breaks = bin_count, plot = F)
    h2 <- hist(pca_data$y, breaks = bin_count, plot = F)
    k <- kde2d(pca_data$x, pca_data$y, n = bin_count)
    image(k, col = plot_color)
    title(main = paste(KNOCKOUT, "_WT", sep = ""))
    #density_plot(pca_data,bin_count,plot_color)

    pca_data_knockout <- scale(data_simulation_knockout, pca$center,
                               pca$scale) %*% pca$rotation
    pca_data_knockout <- data.frame(x = pca_data_knockout[, 1],
                                    y = pca_data_knockout[, 2])
    h1 <- hist(pca_data_knockout$x, breaks = bin_count, plot = F)
    h2 <- hist(pca_data_knockout$y, breaks = bin_count, plot = F)
    k <- MASS::kde2d(pca_data_knockout$x, pca_data_knockout$y, n = bin_count)

    image(k, col = plot_color)
    title(main = paste(KNOCKOUT, "_KO", sep = ""))
    #density_plot(pca_data_knockout,bin_count,plot_color)
    #col_count <- as.integer(65536/dim(data_simulation)[2])
    #heatmap_data <- t(data_simulation[1:col_count,])
    #dist <- as.dist((1-cor(t(heatmap_data), method = "spear"))/2)
    #cluster <- hclust(dist,method = 'ward.D2')
    #dendrogram <- as.dendrogram(cluster)

    full <- heatmap.2(
      t(data_simulation),
      col = plot_color,
      hclustfun = function(x)
        hclust(x, method = 'ward.D2'),
      distfun = function(x)
        as.dist((1 - cor(t(
          x
        ),
        method = "spear")) /
          2),
      trace = "none"
    )
    title(main = paste(KNOCKOUT, "_WT", sep = ""))

    heatmap.2(
      t(data_simulation_knockout[, rev(full$rowInd)]),
      Rowv = NA,
      dendrogram = "column",
      breaks = full$breaks,
      col = plot_color,
      hclustfun = function(x)
        hclust(x, method = 'ward.D2'),
      distfun = function(x)
        as.dist((1 - cor(t(
          x
        ), method = "spear")) / 2),
      trace = "none"
    )
    title(main = paste(KNOCKOUT, "_KO", sep = ""))
    # heatmap.2(t(data_simulation_knockout), col=plot_color, hclustfun = function(x) hclust(x,method = 'ward.D2'), distfun=function(x) as.dist((1-cor(t(x), method = "spear"))/2), trace = "none")
  }
  dev.off()

  message("Plots in the pdf file in the results folder.")

}



#' @export
#' @title Plot the network (Deprecated)
#' @description Use \code{\link{plotCircuit}}.
#'
#' \code{\link{simulateGRC}},  \code{\link{knockdownAnalysis}},
#' \code{\link{overExprAnalysis}},  \code{\link{plotData}},
#' \code{\link{calcHeatmapSimilarity}}
plot_network = function(topology = topology, plotToFile = TRUE) {

  if(plotToFile){
    if (!dir.exists(file.path(getwd(), "results")))
      dir.create(file.path(getwd(), "results"))

    net_file <- paste(getwd(),
                      "/results/network_",
                      topology$filename,
                      ".html",
                      sep = "")
    # setwd(net_file)
  }
  node_list <-
    unique(c(topology$topology[, 1], topology$topology[, 2]))

  nodes <-
    data.frame(
      id = node_list,
      label = node_list,
      font.size = 50,
      value = c(rep(1, length(node_list)))
    )
  edge_col <- data.frame(c(1, 2), c("blue", "darkred"))
  arrow_type <- data.frame(c(1, 2), c("arrow", "circle"))
  colnames(arrow_type) <- c("type", "color")
  colnames(edge_col) <- c("type", "color")
  edges <-
    data.frame(
      from = c(topology$topology[, 1]),
      to = c(topology$topology[, 2])
      #   , arrows = c(c(topology$topology$Target), c(topology$topology$Target))
      #, arrows = "to"
      ,
      arrows.to.type	= arrow_type$color[c(as.numeric(topology$topology[, 3]))]
      ,
      width = 3
      ,
      color = edge_col$color[c(as.numeric(topology$topology[, 3]))]
    )

  #visNetwork(nodes, edges, height = "500px", width = "100%")  %>%
  #visEdges(arrows = "to") %>%
  #visOptions(manipulation = TRUE) %>%
  #visLayout(randomSeed = 123) %>%
  #visPhysics(solver = "forceAtlas2Based")

  network <-
    visNetwork::visNetwork(nodes, edges, height = "1000px", width = "100%") %>%
    #visEdges(arrows = "to") %>%
    visOptions(manipulation = TRUE) %>%
    visLayout(randomSeed = 123) %>%
    #visNodes(scaling = list(label = list(enabled = T))) %>%
    visPhysics(solver = "forceAtlas2Based", stabilization = FALSE)
  if(plotToFile){
    visNetwork::visSave(network, file = net_file, selfcontained = F)
  } else {network}
  #visNetwork(gephi = "test")
}

#' @export
#' @title Plot the density
#' @description Use \code{\link{plotDensity}}.
density_plot <- function(plot_data, bin_count=40, plot_color=NULL){
  densityPlot(plot_data, bin_count, plot_color)
}

#' @export
#' @title Plot the simulated data
#' @description  Use \code{\link{plotData}}.
#'
#' \code{\link{simulateGRC}},  \code{\link{knockdownAnalysis}},
#' \code{\link{overExprAnalysis}},  \code{\link{plotData}},
#' \code{\link{calcHeatmapSimilarity}}
plot_data_deterministic = function(data_file = output_file,
                                   plot_filename = plot_filename,
                                   topology = topology,
                                   bin_count = 40,
                                   plotTsne = FALSE) {
  working_directory <- getwd()
  message("Plotting the results")
  data_simulation <- read.table(output_file, header = FALSE)
  data_simulation <- load_data(data_simulation, topology)

  if (missing(bin_count))
    bin_count <- 40
  if (missing(plot_filename))
    plot_filename <- "plots"

  library(reshape2)
  library(ggplot2)

  library(gplots)
  library(MASS)
  library(RColorBrewer)
  #library(htmlwidgets)
  #library(d3heatmap)
  rf <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
  plot_color <- rf(32)


  pca1 = prcomp(data_simulation, scale. = FALSE)
  pca_data <- data.frame(x = pca1$x[, 1], y = pca1$x[, 2])

  pdf(paste("results/", plot_filename, "_deterministic.pdf", sep = ""))
  density_plot(pca_data, bin_count, plot_color)
  #heatmap(t(data_simulation), col=plot_color, hclustfun = function(x) hclust(x,method = 'complete'))
  heatmap.2(
    t(data_simulation),
    col = plot_color,
    hclustfun = function(x)
      hclust(x, method = 'ward.D2'),
    distfun = function(x)
      as.dist((1 - cor(t(
        x
      ), method = "spear")) / 2),
    trace = "none"
  )
  if (plotTsne) {
    library(Rtsne)
    tsne_out <-
      Rtsne(as.matrix(data_simulation), check_duplicates = FALSE) # Run TSNE
    tsne_data <- data.frame(x = tsne_out$Y[, 1], y = tsne_out$Y[, 2])
    density_plot(tsne_data, bin_count, plot_color)
  }
  dev.off()
  #results_directory <- ifelse(!dir.exists(file.path(working_directory, "results")), dir.create(file.path(working_directory, "results")), file.path(working_directory, "results"))
  #setwd(results_directory)
  #hmap_file <- paste(plot_filename,"_heatmap.html", sep = "")
  #hmap <- d3heatmap(t(data_simulation), col=plot_color)
  #saveWidget(hmap,hmap_file)
  #setwd(working_directory)
  #  data <- list(data_simulation,pca1,tsne_out,pca_data,tsne_data)
  message("Plots in the pdf file in the results folder.")

}

