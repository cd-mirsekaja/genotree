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
matrix_query <- "SELECT ids.IDX, ids.AccessionNumber, taxonomy.ScientificName, taxonomy.Authority, taxonomy.Phylum, taxonomy.Gigaclass, taxonomy.Class, taxonomy.'Order', taxonomy.Family, taxonomy.Genus, taxonomy.taxGroup, traits.isMarine, traits.isBrackish, traits.isFresh, traits.isAllWater, traits.isMarineFresh, traits.isTerrestrial, traits.isMigratory
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
# close the connection to the log file
close(log_file)

# set tip names for all output trees
data_matrix$refTipName <- paste("(",data_matrix$IDX,")",sep="")
data_matrix$mainTipName <- paste(data_matrix$ScientificName," (",data_matrix$IDX,")",sep="")

# import the specified treefile and reroot it at the specified outgroup
rerooted_tree <- import_tree(treefile, outgroup)

# rename taxa in rerooted tree
renamed_tree <- rename_taxa(rerooted_tree, data_matrix, key = 1, value = mainTipName)

# make tree visualization basic plot and set aesthetic parameters
base_plot <- ggtree(renamed_tree,layout = tree_layout)
if (tree_layout=="rectangular" || tree_layout=="roundrect") {
  base_plot <- base_plot+theme_tree2()+geom_rootpoint()+geom_tiplab(offset=0)+theme_tree()
} else if (tree_layout=="circular") {
  base_plot <- base_plot+theme_tree2()+geom_rootpoint()+geom_tiplab(offset=0.5)+geom_treescale()+theme_tree()
}

bootstrap_plot <- get_bootstrap_plot(renamed_tree)




