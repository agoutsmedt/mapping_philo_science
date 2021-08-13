# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
############################## LOADING PACKAGES, PATHS AND OBJECTS ####################################--------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

##################### Packages ############################################--------------

cran_list <- c("data.table",
               "tidyverse",
               "readxl",
               "igraph",
               "tidygraph",
               "leidenAlg",
               "ggplot2",
               "ggraph",
               "ggrepel",
               "ggnewscale",
               "tidytext",
               "textstem",
               "tm",
               "quanteda",
               "scico",
               "sigmajs",
               "ragg",
               "htmlwidgets")
for (p in cran_list) {
  if (p %in% installed.packages() == FALSE) {
    install.packages(p, dependencies = TRUE)
  }
  library(p, character.only = TRUE)
}

github_list <- c("ParkerICI/vite", "agoutsmedt/biblionetwork", "agoutsmedt/networkflow")
for (p in github_list) {
  if (gsub(".*/", "", p) %in% installed.packages() == FALSE) {
    devtools::install_github(p)
  }
  library(gsub(".*/", "", p), character.only = TRUE)
}


######################### Paths and data ##########################################------------

data_path <- "C:/Users/agout/Documents/MEGA/Research/R/data/mapping_philo_science/"
picture_path <- "C:/Users/agout/Documents/MEGA/Research/R/projets/mapping_philo_science/pictures/"
