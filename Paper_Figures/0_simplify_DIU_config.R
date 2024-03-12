#!/usr/bin/env Rscript

library("dplyr")

root_dir <- "/lustre/projects/Research_Project-MRC148213/lsl693/"
root_sfari <- "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/"
root_rb_dir <- "/lustre/projects/Research_Project-MRC148213/Rosie/SFARIdevelopmentalgenomics/"
recovered_dir <- "/lustre/recovered/Research_Project-MRC148213/sl693/RBFetal/"
dirnames <- list(
  
  wholetarg_SQ = paste0(root_rb_dir,"6_sqanti3/"),
  
  output = paste0(root_sfari,"/0_output/"),
  utils = paste0(root_sfari,"/0_utils/"),
  protein = paste0(root_sfari, "/8_longReadProteogenomics/longReadProteogenomics"),
  
  DGE = paste0(root_sfari, "10_deseq//1_DGE/"), 
  DTE = paste0(root_sfari, "10_deseq//2_DTE/"),
  DIU = paste0(root_sfari, "10_deseq//3_DIU/"),
  
  # Leung et al. 2021 PacBio HumanCTX dataset
  humanPacBio = "/lustre/projects/Research_Project-MRC148213/lsl693/PacBioPaper/SQANTI2/HumanCTX"
)


read_DIU <- function(inputPath){
  DIU_targeted <- list.files(path=inputPath,full.names = T, pattern = "resultDIU")
  if(length(DIU_targeted) > 0){
    DIU_targeted <- lapply(DIU_targeted,function(x) read.table(x)[-1,])
    DIU_targeted <- do.call(rbind, DIU_targeted)
    colnames(DIU_targeted) <- c("Gene","p.value","FDR","podiumChange","totalChange")
    DIU_targeted <- DIU_targeted %>% mutate(FDR = as.numeric(as.character(FDR)))
  }else{
    return(NULL)
  }
}


DIU <- list(
  #targetedAge = read_DIU(paste0(dirnames$DIU,"targeted/group")),
  #targetedSex = read_DIU(paste0(dirnames$DIU,"targeted/sex")),
  #wholeAge = read_DIU(paste0(dirnames$DIU,"whole/group")),
  #wholeSex = read_DIU(paste0(dirnames$DIU,"whole/sex")),
  wholeAllAge = read_DIU(paste0(dirnames$DIU,"whole/allGroup")),
  wholeAllSex = read_DIU(paste0(dirnames$DIU,"whole/allSex"))
)
DIUSig <- lapply(DIU, function(x) x[x$FDR <= 0.05, ])

save(DIU, file = paste0(dirnames$output,"DIU.RData"))
save(DIUSig, file = paste0(dirnames$output,"DIUSig.RData"))
