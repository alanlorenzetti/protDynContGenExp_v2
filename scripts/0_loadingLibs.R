# alorenzetti 20200801

# description #####
# this script will load all required libraries
# and set the environment

############ loading packages
# biocmanager is required for loading and installing packages
if(!require("BiocManager")){install.packages("BiocManager"); library("BiocManager")}

# pacman is a nice package manager; make it easier to load and install packages
if(!require("pacman")){install.packages("pacman"); library("pacman")}

# this one isnt stored neither in cran nor in bioconductor
# https://github.com/drostlab/orthologr
library(orthologr) 

# the following is not available at CRAN
# library(devtools) ; install_github("allydunham/tblhelpr")
library(tblhelpr)

# libcurl is required
# list of packages
packs = c("tidyverse",
          "rtracklayer",
          "GenomicRanges",
          "Biostrings",
          "BSgenome",
          "readxl",
          "svglite",
          "ggtree",
          "ComplexHeatmap",
          "circlize",
          "viridis",
          "broom",
          "eulerr",
          "DESeq2",
          "preprocessCore",
          "ggpubr",
          "extrafont",
          "ggtext")
p_load(char = packs)

# loading arial font
font_import(prompt = F)
loadfonts(quiet = T)

# setting ggplot2 theme
mytheme = theme_bw() +
  theme(
    axis.title.x = element_markdown(family = "Arial", size = 10),
    axis.title.y = element_markdown(family = "Arial", size = 10), 
    text = element_text(family = "Arial", size = 10),
    axis.text = element_text(family = "Arial", size = 10),
    strip.text = element_text(family = "Arial", size = 10),
    legend.title = element_text(family = "Arial", size = 10),
    legend.text = element_text(family = "Arial", size = 10),
    title = element_text(family = "Arial", size = 10)
  )

theme_set(mytheme)

# setting working directory
setwd("~/gdrive/documentos/doutorado/isb/protDynContGenExp_v2/")
