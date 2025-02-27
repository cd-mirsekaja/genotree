#!/usr/bin/env Rscript

# remove all variables from environment
rm(list=ls())
# set working directory
setwd("~/Documents/Programming/Bachelorarbeit/tree_recoloring/main")

# load required libraries
library(data.tree)
library(ape)
library(tidyverse)
library(ggtree)
library(treeio)
library(phytools)
library(geiger)
library(randomcoloR)
library(DBI) # SQLite-Library
library(glue) # for making nice text
library(svglite) # for exporting trees as svg for manual annotation
# import custom functions
source("4-1_functions.R")

# set outgroup for rerooting. 298 is Petromyzon marinus, 359 is Asterias rubens, 362 is Pecten maximus
outgroup <- "359"
# set threshold modifier
threshold <- "thr0_35"
# set aster version modifier
aster_ver <- "astralpro3"
# set layout for export tree. Can be 'rectangular', 'roundrect', 'radial', 'circular', ' ellipse'.
tree_layout <- "rectangular"


# set import tree and modificator for saving pdfs
treefile <- paste("data/treefiles-",threshold,"/all-loci_",aster_ver,".treefile", sep = '', collapse = NULL)
savemod <- paste("r", outgroup, "-",threshold,"-",aster_ver,"-",tree_layout, sep = "", collapse = NULL)
out_dir <- paste("output/",savemod,"/",sep="")

# make output subdirectory if it doesn't exist
if (!file.exists(out_dir))
{ dir.create(out_dir) }

# open the connection to the log file
log_file <- file(paste(out_dir,savemod,"_README.log",sep=""),open="wt",blocking=FALSE)

# import data matrix from SQLite database
db_conn <- dbConnect(RSQLite::SQLite(), "data/genotree_master_library_corrected.db")
matrix_query <- "SELECT ids.IDX, ids.AccessionNumber, taxonomy.ScientificName, taxonomy.Authority, taxonomy.Phylum, taxonomy.Gigaclass, taxonomy.Class, taxonomy.'Order', taxonomy.Family, taxonomy.taxGroup, traits.isMarine, traits.isBrackish, traits.isFresh, traits.isAllWater, traits.isMarineFresh, traits.isTerrestrial, traits.isMigratory
FROM ids 
INNER JOIN taxonomy 
ON ids.IDX = taxonomy.IDX
INNER JOIN traits
ON ids.IDX = traits.IDX"
data_matrix <- dbGetQuery(db_conn, matrix_query)
outgroup_name <- dbGetQuery(db_conn, glue_sql(.con=db_conn,"SELECT ScientificName FROM taxonomy WHERE IDX={outgroup}"))$ScientificName
dbDisconnect(db_conn)

headline <- glue("Rooted on {outgroup_name} (IDX {outgroup}).\nThreshold value {threshold}, ASTER version of original consensus tree was {aster_ver}.")
writeLines(headline,con=log_file,sep="\n")

# set tip names for all output trees
data_matrix$refTipName <- paste("(",data_matrix$IDX,")",sep="")
data_matrix$mainTipName <- paste(data_matrix$ScientificName," (",data_matrix$IDX,")",sep="")

# import the specified treefile and reroot it at the specified outgroup
rerooted_tree <- import_tree(treefile, outgroup)

# rename taxa in rerooted tree
renamed_tree <- rename_taxa(rerooted_tree, data_matrix, key = 1, value = mainTipName)

# add phylotree with bootstrap values
bs_tibble <- tibble(
  node=1:Nnode(renamed_tree)+Ntip(renamed_tree),
  bootstrap=renamed_tree$node.label
)

bootstrap_plot <- ggtree(renamed_tree, layout=tree_layout) %<+% bs_tibble
bootstrap_plot <- bootstrap_plot + theme_tree2()+geom_rootpoint()+geom_aline()+geom_text(aes(label=bootstrap),hjust=-0.1,size=3)+geom_tiplab(aes(label=label),align = TRUE)+theme_tree()

# make tree visualization basic plot and set basic parameters
base_plot <- ggtree(renamed_tree,layout = tree_layout)
base_plot <- base_plot+theme_tree2()+geom_rootpoint()+geom_tiplab(offset=0)+theme_tree()


# make reference plot without renamed tips
make_reference(rerooted_tree, paste(out_dir, savemod, "_ReferencePlot.pdf", sep = ""))


### Annotated Trees ###
# set output path for colored trees and set vectors with annotation colors
path_colored <- paste(out_dir,savemod,"_ColoredPlot_",sep="")
anno_colors_class <- c("green3","blue","brown","purple4","blue","violet","green2","grey","blue","black","blue","green4","pink4","red4")
anno_colors_taxgroup <- c("green3","blue","brown","purple4","violet","green2","grey","black","green4","pink4","red4","grey","darkgrey")

# make tree annotated with class
annotate_by_taxgroup("Class",base_plot,path_colored,col_vec = anno_colors_class)
annotate_by_taxgroup("Class",base_plot,path_colored,col_vec = anno_colors_class, include = c("Teleostei","Holostei","Chondrostei","Cladistii","Elasmobranchii"))

# make tree annotated with taxGroup
annotate_by_taxgroup("taxGroup", base_plot, path_colored, col_vec = anno_colors_taxgroup)
annotate_by_taxgroup("taxGroup", base_plot, path_colored, col_vec = anno_colors_taxgroup, include=c("Fish"))

# make assorted annotated trees
annotate_by_taxgroup("Phylum",base_plot,path_colored,colorscheme = "blue")
annotate_by_taxgroup("Order",base_plot,path_colored,colorscheme = "blue")
annotate_by_taxgroup("Family", base_plot,path_colored,colorscheme = "blue")

# save uncolored tree as pdf
pdf(file=paste(out_dir,savemod,"_MainPlot.pdf", sep=""), width=35, height=75)
base_plot
dev.off()

# set output path for bootstrapped trees
path_bootstrap <- paste(out_dir,savemod,"_BootstrapPlot_", sep="")

# make annotated bootstrap trees
annotate_by_taxgroup("Class",bootstrap_plot,path_bootstrap,col_vec = anno_colors_class)
annotate_by_taxgroup("Order",bootstrap_plot,path_bootstrap,colorscheme = "blue")

# save uncolored bootstrap tree as pdf
pdf(file=paste(out_dir,savemod,"_BootstrapPlot.pdf", sep=""), width=35, height=75)
bootstrap_plot
dev.off()


### Subtrees ###
path_subtree <- paste(out_dir,savemod,"_SubtreePlot_",sep="")

# check ASTER version and use different nodes for subtrees
if (aster_ver=="astral4") {
  # block for ASTRAL-IV
  make_subtree(407,"Birds",path_subtree,10,10)
  make_subtree(405,"Sauropsida",path_subtree,10,15)
  make_subtree(457,"Mammals",path_subtree,10,10)
  make_subtree(400,"Fish",path_subtree,20,50)
} else if (aster_ver=="astral-pro3") {
  # block for ASTRAL-Pro3
  make_subtree(419,"Birds",path_subtree,10,10)
  make_subtree(417,"Sauropsida",path_subtree,10,15)
  make_subtree(507,"Amphibians",path_subtree,10,10)
  make_subtree(469,"Mammals",path_subtree,10,10)
  make_subtree(412,"Fish",path_subtree,20,50)
  make_subtree(569,"Salmoniformes",path_subtree,10,5)
  make_subtree(527,"Anguilliformes",path_subtree,10,5)
}


### Trait-Trees ###
# make a ggtree phylo plot for annotating with the traits
main_plot_traits <- annotate_by_taxgroup("Order",base_plot,fill="none",export=FALSE)

# convert Index column in data_matrix to character (ChatGPT)
data_matrix$IDX <- as.character(data_matrix$IDX)
# convert tree to a tibble for easier manipulation (ChatGPT)
tree_data <- as_tibble(rerooted_tree)
# join the data matrix with the tree data based on the tip labels (ChatGPT)
tree_data <- tree_data %>%
  left_join(data_matrix, by = c("label" = "IDX"))

# map traits onto its own tree and save the trees as pdfs
trait_tree(paste(out_dir,savemod,"_TraitPlot_Marine.pdf",sep=''),main_plot_traits,"isMarine","blue")
trait_tree(paste(out_dir,savemod,"_TraitPlot_Fresh.pdf",sep=''),main_plot_traits,"isFresh","lightblue")
trait_tree(paste(out_dir,savemod,"_TraitPlot_Brackish.pdf",sep=''),main_plot_traits,"isBrackish","purple")
trait_tree(paste(out_dir,savemod,"_TraitPlot_AllWater.pdf",sep=''),main_plot_traits,"isAllWater","blue")
trait_tree(paste(out_dir,savemod,"_TraitPlot_Terrestrial.pdf",sep=''),main_plot_traits,"isTerrestrial","orange")
trait_tree(paste(out_dir,savemod,"_TraitPlot_Migration.pdf",sep=''),main_plot_traits,"isMigratory","orange")

# close the connection to the log file
close(log_file)
