sRACIPE_load_topology<-function(){
  working_directory<-getwd()
  if (dir.exists(file.path(working_directory, "inputs"))){
    print("Reading topology file from inputs folder")
    input_directory<-file.path(working_directory, "inputs")
    setwd(input_directory)

  tpo_filename<<-dir(pattern = ".tpo")
  f_tpo <<- read.table(tpo_filename, header = TRUE)
  f_tpo

  setwd(working_directory)
  number_gene <<-length(levels(f_tpo$Source))

  } else {
    print("Missing inputs folder! Please make sure that you are in right working directory and the directory contains an input folder with a .tpo file. Using topology of a toggle switch with one self activating gene. You can directly modify f_tpo directly as well.")}



}
