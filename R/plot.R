load_data = function(data_simulation = data_simulation, topology_df=topology){
  working_directory <- getwd()
  name_genes <- read.table(paste(working_directory,"/results/gene_interaction_topology_",topology$filename,".txt",sep=""), header = T, stringsAsFactors = F)
  name_genes <- t(as.matrix(name_genes))
  data_simulation <- log2(data_simulation)
  data_simulation <- scale(data_simulation)
  name_models <- seq(1:nrow(data_simulation))
  row.names(data_simulation) <- name_models
  colnames(data_simulation) <- name_genes
  return(data_simulation)
}

plot_interactive_heatmap = function(data_file=output_file, plot_filename=plot_filename,  topology_df=topology){
  working_directory <- getwd()
  message("Plotting the interactive heatmap")
  data_simulation <- read.table(output_file, header = F)
  data_simulation <- load_data(data_simulation, topology)

  if(missing(plot_filename)) plot_filename <- "heatmap"

  library(reshape2)
  library(ggplot2)
  library(Rtsne)
  library(gplots)
  library(MASS)
  library(RColorBrewer)
  library(htmlwidgets)
  library(d3heatmap)

  rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
  plot_color <- rf(32)

  results_directory <- ifelse(!dir.exists(file.path(working_directory, "results")), dir.create(file.path(working_directory, "results")), file.path(working_directory, "results"))
  setwd(results_directory)
  hmap_file <- paste(plot_filename,"_heatmap.html", sep = "")
  hmap <- d3heatmap(t(data_simulation), col=plot_color)
  saveWidget(hmap,hmap_file)
  setwd(working_directory)
  message("Heatmaps html file is in the results folder.")
}
plot_data_deterministic = function(data_file=output_file, plot_filename=plot_filename,  topology_df=topology, bin_count=40){
  working_directory <- getwd()
  message("Plotting the results")
  data_simulation <- read.table(output_file, header = F)
  data_simulation <- load_data(data_simulation, topology)

  if(missing(bin_count)) bin_count <- 40
  if(missing(plot_filename)) plot_filename <- "plots"

  library(reshape2)
  library(ggplot2)
  library(Rtsne)
  library(gplots)
  library(MASS)
  library(RColorBrewer)
  #library(htmlwidgets)
  #library(d3heatmap)
  rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
  plot_color <- rf(32)

  tsne_out <- Rtsne(as.matrix(data_simulation), check_duplicates = FALSE) # Run TSNE
  tsne_data <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2])

  pca1 = prcomp(data_simulation, scale. = FALSE)
  pca_data <- data.frame(x=pca1$x[,1],y=pca1$x[,2])

  pdf(paste("results/",plot_filename,"_deterministic.pdf",sep = ""))
  density_plot(pca_data,bin_count,plot_color)
  heatmap(t(data_simulation), col=plot_color, hclustfun = function(x) hclust(x,method = 'complete'))
  density_plot(tsne_data,bin_count, plot_color)
  dev.off()
  #results_directory <- ifelse(!dir.exists(file.path(working_directory, "results")), dir.create(file.path(working_directory, "results")), file.path(working_directory, "results"))
  #setwd(results_directory)
  #hmap_file <- paste(plot_filename,"_heatmap.html", sep = "")
  #hmap <- d3heatmap(t(data_simulation), col=plot_color)
  #saveWidget(hmap,hmap_file)
  #setwd(working_directory)
  data <- list(data_simulation,pca1,tsne_out,pca_data,tsne_data)
  message("Plots in the pdf file in the results folder.")

}

plot_data_stochastic = function(output_file, plot_filename=filename,  topology_df=topology, config = configuration, bin_count=40){
  message("Plotting the stochastic results")
  library(reshape2)
  library(ggplot2)
  library(Rtsne)
  library(gplots)
  library(MASS)
  library(RColorBrewer)

  rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
  plot_color <- rf(32)

  if(missing(bin_count)) bin_count <- 40
  if(missing(plot_filename)) output_filename <- "plots"

  working_directory <- getwd()
  name_genes <- read.table(paste(working_directory,"/results/gene_interaction_topology_",topology$filename,".txt",sep=""), header = T, stringsAsFactors = F)

 data_simulation_all <- as.data.frame(read.table(output_file, header = F))
 col_start <- topology$number_gene*(configuration$NOISE_LEVELS-1)+1
 col_end <- topology$number_gene*(configuration$NOISE_LEVELS)
 data_simulation <- as.data.frame(data_simulation_all[,col_start:col_end])
 data_simulation <- log2(data_simulation)
 data_simulation <- data_simulation[is.finite(rowMeans(data_simulation)), ]

 means <- colMeans(data_simulation)
 sds <- apply(data_simulation, 2, sd)
 data_simulation <- sweep(data_simulation,2,means,FUN = "-")
 data_simulation <- sweep(data_simulation,2,sds,FUN = "/")

 name_models <- seq(1:nrow(data_simulation))
 name_genes <- t(as.matrix(name_genes))
 row.names(data_simulation) <- name_models
 colnames(data_simulation) <- name_genes


 pca = prcomp(data_simulation, scale. = FALSE, center = FALSE)
 #pca_data <- data.frame(x=pca$x[,1],y=pca$x[,2])

 pdf(paste("results/",plot_filename,"_stochastic_pca.pdf",sep = ""))

 for(i in 1:configuration$NOISE_LEVELS)
 {
   col_start <- (topology$number_gene)*(configuration$NOISE_LEVELS-i)+1
   col_end <- (topology$number_gene)*(configuration$NOISE_LEVELS-i+1)
   data_simulation <- as.data.frame(data_simulation_all[,col_start:col_end])
   data_simulation <- log2(data_simulation)
   data_simulation <- data_simulation[is.finite(rowMeans(data_simulation)), ]

   data_simulation <- sweep(data_simulation,2,means,FUN = "-")
   data_simulation <- sweep(data_simulation,2,sds,FUN = "/")

   name_models <- seq(1:nrow(data_simulation))
   row.names(data_simulation) <- name_models
   colnames(data_simulation) <- name_genes

   pca_data <- scale(data_simulation, pca$center, pca$scale) %*% pca$rotation
   pca_data <- data.frame(x=pca_data[,1],y=pca_data[,2])
   density_plot(pca_data,bin_count,plot_color)
 }


  dev.off()

  message("Plots in the pdf file in the results folder.")

}

density_plot = function(plot_data,bin_count, plot_color){
  colnames(plot_data) <- c("x", "y")

  h1 <- hist(plot_data$x, breaks=bin_count, plot=F)
  h2 <- hist(plot_data$y, breaks=bin_count, plot=F)
  top <- max(h1$counts, h2$counts)
  k <- kde2d(plot_data$x, plot_data$y, n=bin_count)

  # margins
  oldpar <- par()
  par(mar=c(3,3,1,1))
  layout(matrix(c(2,0,1,3),2,2,byrow=T),c(3,1), c(1,3))
  image(k, col=plot_color) #plot the image
  par(mar=c(0,2,1,0))
  barplot(h1$counts, axes=F, ylim=c(0, top), space=0, col='red')
  par(mar=c(2,0,0.5,1))
  barplot(h2$counts, axes=F, xlim=c(0, top), space=0, col='red', horiz=T)
}

plot_network = function(topology = topology){
  library(visNetwork)

  working_directory <- getwd()
  #  setwd(paste(working_directory,"/results",sep = ""))


  net_file <- paste(working_directory,"/results/network_",topology$filename,".html",sep="")
  # setwd(net_file)

  node_list <- unique(c(topology$topology[,1], topology$topology[,2]))

  nodes <- data.frame(id = node_list, label = node_list, font.size =50, value=c(rep(1,length(node_list))))
  edge_col <- data.frame(c(1,2),c("blue","darkred"))
  colnames(edge_col) <- c("type", "color")
  edges <- data.frame(from = c(topology$topology[,1]), to = c(topology$topology[,2])
                      #   , arrows = c(c(topology$topology$Target), c(topology$topology$Target))
                      , arrows = "to"
                      , width = 2
                      , color = edge_col$color[c(as.numeric(topology$topology[,3]))]
  )

  #visNetwork(nodes, edges, height = "500px", width = "100%")  %>%
  #visEdges(arrows = "to") %>%
  #visOptions(manipulation = TRUE) %>%
  #visLayout(randomSeed = 123) %>%
  #visPhysics(solver = "forceAtlas2Based")

  network <- visNetwork(nodes, edges, height = "1000px", width = "100%") %>%
    #visEdges(arrows = "to") %>%
    visOptions(manipulation = TRUE) %>%
    visLayout(randomSeed = 123) %>%
    #visNodes(scaling = list(label = list(enabled = T))) %>%
    visPhysics(solver = "forceAtlas2Based", stabilization = FALSE)

  visSave(network, file = net_file, selfcontained = F)
  #visNetwork(gephi = "test")
}

