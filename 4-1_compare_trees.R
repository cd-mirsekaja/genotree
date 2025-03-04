#!/usr/bin/env Rscript

# remove all variables from environment
rm(list=ls())

# load required libraries
#library(data.tree)
library(ape)
library(tidyverse)
library(ggtree)
library(treeio)
#library(randomcoloR)
library(DBI) # SQLite-Library
library(glue) # for making nice text
library(dendextend)
library(phylogram)
# import custom functions
source("/Users/privatstudium/Documents/Programming/Bachelorarbeit/0_genotree_repository/4-2_load_trees.R")
source("/Users/privatstudium/Documents/Programming/Bachelorarbeit/0_genotree_repository/4-3_functions.R")


# set working directory
setwd("~/Documents/Programming/Bachelorarbeit/tree_recoloring/main")

# set outgroup for rerooting. 298 is Petromyzon marinus, 359 is Asterias rubens, 362 is Pecten maximus
outgroup <- "359"
# set threshold modifier
threshold <- "thr0_35"
# set layout for export tree. Can be 'rectangular', 'roundrect', 'circular'.
tree_layout <- "rectangular"
# set aster version modifier
aster_ver_1 <- "astralpro3"
aster_ver_2 <- "astral4"

treefile1 <- paste("data/treefiles-",threshold,"/all-loci_",aster_ver_1,".treefile", sep = '', collapse = NULL)

tree_ap3_data <- load_supertree(outgroup,threshold,aster_ver_1,tree_layout)
tree_a4_data <- load_supertree(outgroup,threshold,aster_ver_2,tree_layout)

tree_ap3 <- tree_ap3_data$renamed_tree
tree_a4 <- tree_a4_data$renamed_tree

dend_ap3 <- as.dendrogram.phylo(tree_ap3_data$renamed_tree)
dend_a4 <- as.dendrogram.phylo(tree_a4_data$renamed_tree)
clus10 = cutree(dend_ap3,k=4)

plot(as.phylo(dend_ap3), type = "fan", tip.color = colors[clus10],
     label.offset = 1, cex = 0.7)

tang <- tanglegram(dend_ap3, dend_a4, 
                   main_left = "ASTRAL-Pro3", main_right = "ASTRAL-IV",
                   main = "Tanglegram of the two main supertrees",
                   highlight_branches_col = TRUE,
                   margin_top = 6)
