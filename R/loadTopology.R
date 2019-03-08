
#' @export
#' @title Load the topology from a file.
#' @description Read an input file and return the circuit information
#'  as a \code{list}.
#' A typical topology file looks like
#' \tabular{lcr}{
#'   Source \tab Target \tab Type \cr
#'   geneA \tab geneB \tab 2 \cr
#'   geneB \tab geneC \tab 1 \cr
#'   geneB \tab geneA \tab 2
#' }
#' Here the regulation type is specified by number - activation: \code{1},
#'  inhibition: \code{2}
#' @param circuit \code{Character or data.frame}. Topology file address or as a
#' data.frame.
#' @param header (optional) \code{Logical}. Whether the input file has a header
#' or not.
#' \code{TRUE} by default.
#' @return \code{List}. Contains the number of genes, circuit name,
#' topology (Source, Target, Type), and directory path.
#' @examples
#' topologyInfo <- loadCircuit("inputs/test.tpo")
#'
#' @section Related Functions:
#'
#' \code{\link{simulateGRC}},  \code{\link{knockdownAnalysis}},
#' \code{\link{overExprAnalysis}},  \code{\link{plotData}},
#' \code{\link{calcHeatmapSimilarity}}


loadCircuit <-function( circuit="inputs/test.tpo", header = TRUE){
  if(missing(circuit)){
    print("Please specify the topology file")
    return()
  }
  if(class(circuit) == "data.frame"){

      f_tpo <- circuit
      colnames(f_tpo) <- c("Source","Target","Type")
      filename <- deparse(substitute(circuit))
      inputDir <- ifelse(!dir.exists(file.path(getwd(), "inputs")),
                                             dir.create(file.path(getwd(), "inputs")), TRUE)
      write.table(f_tpo, file = paste0("inputs/",filename,".txt"),
                  quote = FALSE,row.names = FALSE)
      circuit <- paste0("inputs/",filename,".txt")
  }

  if(file.exists(circuit)){
  f_tpo <- read.table(circuit, header = header, stringsAsFactors = FALSE)
  colnames(f_tpo) <- c("Source","Target","Type")
  tpo_filename <- basename(circuit)
  filename <- tools::file_path_sans_ext(tpo_filename)
  number_gene <-  length(unique(c(f_tpo$Source,f_tpo$Target)))
  topologyInfo <- list(number_gene = number_gene, filename=filename,
                       topology=f_tpo, topology_filepath = circuit)
  }

  else
  {
    print("Please check the input filename. Make sure that full path is
          specified. For example, if it is in inputs folder, use inputs/filename
          as argument to the function. ")
  }

#  message("Topology file successfully loaded")
  return(topologyInfo)
}


#' @export
#' @title (deprecated) Load the topology from a file.
#' @description Read an input file and return the circuit information as a list.
#' This function will be removed in future versions.
#' Use the function \code{\link{loadTopology}} instead.
#' @param topology_file Character. Topology file address.
#' @return List. Contains the number of genes, circuit name,
#' topology (Source, Target, Type), and directory path.
#'
sRACIPE_load_topology<-function( topology_file="inputs/test.tpo"){
  loadCircuit(topology_file)
}

