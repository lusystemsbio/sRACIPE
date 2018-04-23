sRACIPE_load_topology<-function( topology_file="inputs/pdx1.tpo"){
  working_directory<-getwd()
  if(missing(topology_file)){
    if (dir.exists(file.path(working_directory, "inputs"))){
      print("Reading topology file from inputs folder")
      input_directory<-file.path(working_directory, "inputs")
      setwd(input_directory)
      tpo_filename<-dir(pattern = ".tpo")
      if(length(tpo_filename)>1){
        print("There are multiple tpo files in the inputs folder. Please specify the filename to be used as input topology file using its filepath.")
        setwd(working_directory)
        return()
      }
      else {topology_file <- file.path(working_directory,"inputs",tpo_filename)}
    }
    else {
      print("Missing inputs folder! Please make sure that you are in right working directory and the directory contains an input folder with a .tpo file. Using topology of a toggle switch with one self activating gene. You can directly modify f_tpo directly as well.")
      setwd(working_directory)
      return()
    }
  }
    else{
      tpo_filename <- basename(topology_file)
    }
  if(file.exists(topology_file)){
  f_tpo <- read.table(topology_file, header = TRUE,stringsAsFactors = FALSE)
  colnames(f_tpo) <- c("Source","Target","Type")
  filename <- tools::file_path_sans_ext(tpo_filename)
  length(unique(c(f_tpo$Source,f_tpo$Target)))
  number_gene <-  length(unique(c(f_tpo$Source,f_tpo$Target)))
  topology <- list(number_gene = number_gene, filename=filename, topology=f_tpo, topology_filepath = topology_file)}
  else
  {
    print("Please check the input filename. Make sure that full path is specified. For example, if it is in inputs folder, use inputs/filename as argument to the function. ")
  }
  setwd(working_directory)
  message("Topology file successfully loaded")
  return(topology)
}
