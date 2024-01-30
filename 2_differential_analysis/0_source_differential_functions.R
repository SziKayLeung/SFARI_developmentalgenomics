## ---------- Script -----------------
##
## Script name:
##
## Purpose of script: sources functions for generating downstream plots for SFARI dataset
##
## Author: Szi Kay Leung
##
## Email: S.K.Leung@exeter.ac.uk
##
## ---------- Notes -----------------
##
##
##
##
##

## ---------- Packages -----------------

suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(rjson)) # json files
suppressMessages(library(plyr)) # revalue
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(reshape))
suppressMessages(library(gridExtra))
suppressMessages(library(grid))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(viridis))
suppressMessages(library(wesanderson))
suppressMessages(library(extrafont))
suppressMessages(library(tidyr))
suppressMessages(library(purrr))
suppressMessages(library(tibble))
suppressMessages(library(VennDiagram))
suppressMessages(library(directlabels))
suppressMessages(library(cowplot))
suppressMessages(library(readxl))
suppressMessages(library(ggdendro))
suppressMessages(library(pheatmap))
suppressMessages(library(extrafont))
suppressMessages(loadfonts())

## ----------Functions-----------------

# load all the functions
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/General/5_TappAS_Differential/plot_aesthetics.R")
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/General/5_TappAS_Differential/plot_tappas_analysis.R")
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/General/5_TappAS_Differential/sqanti_general.R")


## ----------Plot colours-----------------

# plot label colour
label_colour <- function(genotype){
  if(genotype == "Braak0"){colour = wes_palette("Royal1")[1]}else{
    if(genotype == "Braak1"){colour = wes_palette("Royal2")[3]}else{
      if(genotype == "Braak6"){colour = wes_palette("Royal1")[2]}else{
        if(genotype == "novel"){colour = wes_palette("Darjeeling1")[4]}else{
          if(genotype == "known"){colour = wes_palette("Darjeeling1")[5]}else{
            if(genotype == "targeted"){colour = wes_palette("Darjeeling1")[2]}else{
              if(genotype == "whole"){colour = wes_palette("Darjeeling1")[1]}else{
                if(genotype == "whole+targeted"){colour = wes_palette("Darjeeling2")[1]}else{
                  if(genotype == "Control"){colour = wes_palette("Royal1")[1]}else{
                    if(genotype == "Case"){colour = wes_palette("Royal1")[2]}else{
                      if(genotype == c("Yes")){colour = alpha(wes_palette("Cavalcanti1")[4],0.8)}else{
                        if(genotype == c("No")){colour = alpha(wes_palette("Cavalcanti1")[5],0.5)}else{
                        }}}}}}}}}}}}
  return(colour)
}



label_group <- function(genotype){
  if(genotype == "Case"){group = "Adult"}else{
    if(genotype == c("Control")){group = "Fetal"}else{
      if(genotype == c("Case_timeanalysis")){group = "Male"}else{
        if(genotype == c("Control_timeanalysis")){group = "Female"}else{
          if(genotype == c("0_timeanalysis")){group = "Fetal"}else{
            if(genotype == c("1_timeanalysis")){group = "Adult"}else{
}}}}}}
  return(group)
}
