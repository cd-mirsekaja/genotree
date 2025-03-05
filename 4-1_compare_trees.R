#!/usr/bin/env Rscript

# remove all variables from environment
rm(list=ls())

# load required libraries
library(ape)
library(tidyverse)
library(ggtree)
library(treeio)
library(DBI) # SQLite-Library
library(glue)
library(randomcoloR)
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
threshold_1 <- "thr0_35"
threshold_2 <- "thr0_35-realigned"
# set layout for export tree. Can be 'rectangular', 'roundrect', 'circular'.
tree_layout <- "rectangular"
# set aster version modifier
aster_ver_1 <- "astralpro3"
aster_ver_2 <- "astral4"
# set save modifier
savemod <- paste0("r", outgroup, "-",threshold_1)
# set output directory and make it if neccessary
out_dir <- paste0("output/",savemod,"-tanglegrams/")
if (!file.exists(out_dir))
{ dir.create(out_dir) }

# use custom function to import tree objects
tree_ap3_data <- load_supertree(outgroup,threshold_1,aster_ver_1,tree_layout)
tree_ap3_real_data <- load_supertree(outgroup,threshold_2,aster_ver_1,tree_layout)
tree_a4_data <- load_supertree(outgroup,threshold_1,aster_ver_2,tree_layout)
tree_a4_real_data <- load_supertree(outgroup,threshold_2,aster_ver_2,tree_layout)

# get the two renamed trees from imported data
tree_ap3 <- tree_ap3_data$renamed_tree
tree_ap3_real <- tree_ap3_real_data$renamed_tree
tree_a4 <- tree_a4_data$renamed_tree
tree_a4_real <- tree_a4_real_data$renamed_tree

# turn phylo objects into dendrograms
dend_ap3 <- as.dendrogram.phylo(tree_ap3)
dend_ap3_real <- as.dendrogram.phylo(tree_ap3_real)
dend_a4 <- as.dendrogram.phylo(tree_a4)
dend_a4_real <- as.dendrogram.phylo(tree_a4_real)

# construct tanglegram for the two supertrees
# warning: graphics panel needs certain size
supertree_tang_ap3-a4 <- tanglegram(dend1 = dend_ap3, dend2 = dend_a4)
supertree_tang_ap3_real-a4_real <- tanglegram(dend1 = dend_ap3_real, dend2 = dend_a4_real)
supertree_tang_ap3-ap3_real <- tanglegram(dend1 = dend_ap3, dend2 = dend_ap3_real)
supertree_tang_a4-a4_real <- tanglegram(dend1 = dend_a4, dend2 = dend_a4_real)

# save the tanglegram as pdf
save_supertree_tang <- paste0(out_dir,savemod,"-Tanglegram_Compare_AP3_real_A4_real.pdf")
pdf(save_supertree_tang,width=40,height=45)
# plot tanglegram with modifiers
plot(supertree_tang_ap3_real-a4_real,
     main_left = "ASTRAL-Pro3 (re-aligned)", main_right = "ASTRAL-IV (re-aligned)",
     main = "Tanglegram of the two main supertrees",
     highlight_branches_col = TRUE,
     #common_subtrees_color_branches = TRUE,
     margin_top = 10,
     columns_width = c(6,4,6),
     sort = FALSE)
dev.off()

tangle_plots <- function(dend_a, dend_b, out_dir,savemod1,savemod2,treename_a,treename_b)
{
  tangle <- tanglegram(dend1 = dend_a, dend2 = dend_b)
  savepath <- paste0(out_dir,savemod1,"-Tanglegram_Compare-",savemod2)
  pdf(savepath,width=40,heigth=45)
  
  plot(tangle,
       main_left = treename_a, main_right = treename_b,
       main = "Tanglegram of the two main supertrees",
       highlight_branches_col = TRUE,
       #common_subtrees_color_branches = TRUE,
       margin_top = 10,
       columns_width = c(6,4,6),
       sort = FALSE)
  dev.off()
  print(paste0("Tanglegram for ",treename_a," and ",treename_b," saved as\n",savepath))
}


#col_rand <- randomcoloR::distinctColorPalette(k=40)

#unbranch_dend_ap3 <- unbranch(dend_ap3,k=length(col_rand))
#clus = dendextend::cutree(unbranch_dend_ap3, length(col_rand))
#plot(as.phylo(dend_ap3), type = "phylogram", tip.color = col_rand[clus],
#     label.offset = 0, cex = 0.7)
