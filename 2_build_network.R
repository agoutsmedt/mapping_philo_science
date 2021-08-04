#' ---
#' title: "Script for building the different types of networks"
#' author: "Aurélien Goutsmedt"
#' date: "/ Last compiled on `r format(Sys.Date())`"
#' output: 
#'   github_document:
#'     toc: true
#'     number_sections: true
#' ---

#+ r setup, include = FALSE
knitr::opts_chunk$set(eval = FALSE)

#' # What is this script for?
#' 
#' In this script, we build different networks (cocitation, coupling and semantic similarity)
#' for different subperiods. We use the Leiden algorithm to identify communities,
#' and calculate spatial coordinates using Force Atlas 2 algorithm.
#' 
#' > WARNING: This script represents a first step of the project, and some processes have been
#' improved (notably by the creation of new functions).
#' 
#' # Loading packages, paths and data
#' 
#' ## External scripts

source("functions_for_network_analysis.R")
source("packages_and_paths.R")

#' ## Loading Data

nodes <- read_excel(paste0(data_path, "Philo_des_Sciences_2010-2020_Gender_VL.xlsx")) %>% 
  as.data.table()
direct_citations <- readRDS(paste0(data_path, "cleaned_direct_citations.rds"))
