
# function for rerooting a phylotree object
reroot_tree <- function(phylotree, outgroup = "none")
{
  # if no outgroup is specified, return unrooted tree
  if (outgroup == "none")
  { 
    print("No outgroup supplied, returning unrooted tree...")
    return(phylotree)
  } else {
    # get the node for the outgroup for rerooting.
    reroot_node <- ggtree::MRCA(phylotree, .node1 = outgroup)
    # reroot at specified node in the tree
    rerooted_tree <- ape::root(phylotree, reroot_node)
    print("Reroot successful, returning new tree...")
    return(rerooted_tree)
  }
}

# function for importing the specified tree and rerooting it
import_tree <- function(path, outgroup = "none")
{
  # import treefile
  tree <- read.tree(path)
  # if no outgroup is specified, return unrooted tree
  out_tree <- reroot_tree(tree,outgroup)
  return(out_tree)
}

# function for making a reference plot and saving it as a file
make_reference <- function(phylotree, savepath)
{
  renamed_reftree <- rename_taxa(phylotree, data_matrix, key = 1, value = mainTipName)
  reference_plot = ggtree(renamed_reftree, layout = "rectangular")
  reference_plot = reference_plot + theme_tree2() + geom_tiplab(offset=0.03) + geom_text(aes(label=node),hjust=-.3) + geom_rootpoint()
  # save generated tree as pdf
  ggsave(savepath,reference_plot,device="pdf",dpi=300,width=30,height=75,limitsize=FALSE)
}

# A helper function that gets all descendant tips of a given node.
# This version uses extract.clade to get tip labels.
get_tips <- function(tree, node) {
  clade <- extract.clade(tree, node)
  return(clade$tip.label)
}

# Function to get the highest (maximal) pure clades for a given taxon group.
get_pure_block_nodes <- function(taxgroup, tree, taxgroup_type) {
  # Get the taxa (e.g., species IDs) that should be in this taxonomic group.
  taxon_tips <- intersect(data_matrix$IDX[data_matrix[[taxgroup_type]] == taxgroup], tree$tip.label)
  n_tips <- length(tree$tip.label)
  pure_blocks <- c()
  
  # Loop over all internal nodes
  for (node in (n_tips + 1):(n_tips + tree$Nnode)) {
    node_tips <- get_tips(tree, node)
    # Check if the node is pure (all its descendant tips are from taxon_tips)
    if (all(node_tips %in% taxon_tips)) {
      # Identify the parent node from the tree edge matrix.
      parent <- tree$edge[tree$edge[, 2] == node, 1]
      parent_is_pure <- FALSE
      if (length(parent) > 0) {
        parent_tips <- get_tips(tree, parent)
        if (all(parent_tips %in% taxon_tips)) {
          parent_is_pure <- TRUE
        }
      }
      # Only record nodes that are not nested within a larger pure block.
      if (!parent_is_pure) {
        pure_blocks <- c(pure_blocks, node)
      }
    }
  }
  
  # Return the set of maximal nodes that define pure blocks.
  return(pure_blocks)
}

# Revised main function that uses the pure block approach for paraphyletic groups.
get_mrca <- function(taxgroup, tree, taxgroup_type) {
  # Get the list of taxa IDs corresponding to the taxon group.
  taxon_tips <- intersect(data_matrix$IDX[data_matrix[[taxgroup_type]] == taxgroup], tree$tip.label)
  
  if (length(taxon_tips) > 1) {
    # If the whole group is monophyletic, we get one block.
    if (is.monophyletic(tree, taxon_tips)) {
      return(getMRCA(tree, tip = taxon_tips))
    } else {
      # For paraphyletic groups, return the maximal pure blocks.
      return(get_pure_block_nodes(taxgroup, tree, taxgroup_type))
    }
  } else if (length(taxon_tips) == 1) {
    # For a single tip, simply return that tip's node.
    return(MRCA(tree, .node1 = taxon_tips))
  } else {
    return("not in tree")
  }
}

# function for annotating the plot with colours
apply_annotation <- function(anno_plot, phylotree, lb, cl, fill = "none", cl_type="taxGroup")
{
  # search for most recent common ancestors of given taxon group
  node_vector <- get_mrca(lb, phylotree, cl_type)
  # if the taxon group is not in the tree, return unmodified plot
  if (node_vector[1] == "not in tree" || length(node_vector) == 0)
  {return(anno_plot)}
  
  # apply annotations for all clades in the tree
  for (nd in node_vector)
  {
    # set the label for the given clade
    anno_plot <- anno_plot + geom_cladelabel(node = nd, label = lb, color = cl, offset = .1,align = TRUE)
    
    # highlight nodes with specified parameters
    if (fill == "full")
    {anno_plot <- anno_plot + geom_hilight(node = nd,fill = cl,color = "black")}
    else if (fill == "rim")
    {anno_plot <- anno_plot + geom_hilight(node = nd, fill = "transparent", color = cl, linewidth = 0.75)}
    else if (fill == "none")
    {}
  }
  
  # return the modified plot
  return(anno_plot)
}

# function for annotating a ggtree phylo plot with colors and taxon group names
annotate_by_taxgroup <- function(taxgroup, phyloplot, phylotree, path="", fill="full", colorscheme="monochrome", col_vec=c(), include=c(), wd=30, ht=75, export=TRUE, output_format="pdf")
{
  tax_anno_plot <- phyloplot
  color_count <- length(unique(data_matrix[[taxgroup]]))
  # if the amount of colors given is equal to the amount of taxa
  if (length(col_vec)!=0 && length(col_vec)==color_count)
  { anno_colors <- col_vec }
  # if there are more taxa than colors given
  else if (length(col_vec)!=0 && length(col_vec)<color_count)
  {  
    col_diff <- color_count-length(col_vec)
    col_rand <- randomColor(count=col_diff,hue=colorscheme, luminosity = "dark")
    anno_colors <- c(col_vec,col_rand)
  }
  # if there are less taxa than colors given
  else if (length(col_vec)!=0 && length(col_vec)>color_count)
  { anno_colors <- col_vec[1:color_count] }
  # if no colors were given
  else
  { anno_colors <- randomColor(count=color_count, hue=colorscheme, luminosity = "dark") }
  
  print(paste("=== Generating tree for ",taxgroup," ===",sep = ""))
  index <- 1
  for (entry in unique(data_matrix[[taxgroup]])){
    print(paste(entry,anno_colors[index],sep=": "))
    
    if (length(include)==0)
    { tax_anno_plot <- apply_annotation(tax_anno_plot,phylotree,entry,anno_colors[index],fill=fill,cl_type=taxgroup) }
    else if (length(include)>0 && entry %in% include)
    { tax_anno_plot <- apply_annotation(tax_anno_plot,phylotree,entry,anno_colors[index],fill=fill,cl_type=taxgroup) }
    else
    { tax_anno_plot <- apply_annotation(tax_anno_plot,phylotree,entry,anno_colors[index],fill="none",cl_type=taxgroup) }
    index <- index+1
  }
  
  if (export==TRUE)
  {
    if (length(include)==0)
    { save_path <- paste0(path,taxgroup,".",output_format) }
    else
    { save_path <- paste0(path,taxgroup,"_only-",paste(include,collapse="-"),".",output_format)  }
    # save the created plot as a pdf
    ggsave(save_path,tax_anno_plot,device=output_format,width=wd,height=ht,limitsize=FALSE)
    print(paste("Colored tree saved as",save_path))
  }
  else
  { return(tax_anno_plot) }
  
}

# function for making a bootstrap phyloplot
get_bootstrap_plot <- function(phylotree,filter_value=0.95)
{
  # add phylotree with bootstrap values, by default shows values under 0.95
  bootstrap_tibble <- tibble(
    node=1:Nnode(phylotree)+Ntip(phylotree),
    bootstrap = ifelse(as.numeric(phylotree$node.label) > filter_value, "", phylotree$node.label))
  
  bootstrap_plot <- ggtree(phylotree, layout="rectangular") %<+% bootstrap_tibble
  bootstrap_plot <- bootstrap_plot + theme_tree2()+geom_rootpoint()+geom_text(aes(label=bootstrap),hjust=-0.1,size=3)+geom_tiplab(aes(label=label),offset=0.02)+theme_tree()
  
  return(bootstrap_plot)
}

make_comparison_subtree <- function(input_node,clade,phylotree,path,wd=10,ht=10,output_format="jpeg")
{
  print(paste0("=== Generating Comparison Subtree for ",clade," ==="))
  # create pdf device for saving
  save_path <- paste0(path,"Comp_",clade,".",output_format)
  # create device for saving
  if (output_format=="pdf")
    { pdf(file=save_path,width=wd,height=ht,units="in") }
  else if (output_format=="jpeg")
    { jpeg(file=save_path,width=wd,height=ht,units="in",res=300) }
  
  # get all descendants of input node
  tips <- c(getDescendants(phylotree,input_node))
  # get subtree of input node and title it
  subtree <- zoom(phylotree,tips,subtree = FALSE,align.tip.label=FALSE)
  
  # print the tree out
  print(subtree)
  dev.off()
  
  print(paste("Subtree tree saved as",save_path))
}

make_annotation_subtree <- function(input_node,clade,phylotree,path,anno_group="none",tree_layout="roundrect",reroot_node="none",node_labels=FALSE,anno_color="blue",wd=10,ht=10,export=TRUE,output_format="pdf")
{
  print(paste("=== Generating Annotated Subtree for ",clade," ===",sep=""))
  all_subtrees <- subtrees(phylotree)
  for (entry in all_subtrees) {
    if (entry$name==input_node){
      sub_tree <- entry
    }
  }
  if (!exists("sub_tree")){
    print("Node not found in tree.")
    exit
  } else {
    if (reroot_node!="none")
    { sub_tree <- reroot_tree(sub_tree,reroot_node) }
    
    renamed_subtree <- rename_taxa(sub_tree, data_matrix, key = 1, value = mainTipName)
    sub_plot <- ggtree(renamed_subtree,layout=tree_layout)
    if (node_labels==TRUE) {
    sub_plot <- sub_plot +
      theme_tree() +
      geom_rootpoint() +
      geom_tiplab(offset = 0.0005) +
      geom_text(aes(label=node),hjust=.1) +
      coord_cartesian(clip = "off") +
      hexpand(0.3) + 
      theme(plot.margin = margin(10, 10, 10, 10))
    } else {
      sub_plot <- sub_plot +
        theme_tree() +
        geom_rootpoint() +
        geom_tiplab(offset = 0.0005) +
        coord_cartesian(clip = "off") +
        hexpand(0.3) + 
        theme(plot.margin = margin(10, 10, 10, 10))
    }
  }
  save_path <- paste0(path,"Anno_",clade)
  if (anno_group!="none")
  {
    out_plot <- annotate_by_taxgroup(anno_group,sub_plot,sub_tree,path=save_path,colorscheme=anno_color,wd=wd,ht=ht,export=export,output_format=output_format)
  } else {
    out_plot <- sub_plot
    save_path <- paste0(save_path,".",output_format)
    ggsave(save_path,out_plot,device=output_format,width=wd,height=ht,limitsize=FALSE,units="in")
    print(paste0("Subtree saved as ",save_path))
  }
  return(out_plot)
}


# Helper function: Compute ancestral states in a bottom-up (postorder) pass.
fast_get_ancestral_state <- function(tree, tip_states) {
  Ntip <- length(tree$tip.label)
  Nnode <- tree$Nnode
  total_nodes <- Ntip + Nnode
  # Preallocate a vector for all node states; initializing with NA.
  states <- rep(NA, total_nodes)
  
  # Fill in states for the tip nodes.
  # Ensure tip_states is named with tip labels.
  states[1:Ntip] <- tip_states[tree$tip.label]
  
  # Process internal nodes in numerical order (internal nodes are often > Ntip).
  # Assuming that for a rooted tree, for every internal node the children have lower numbers.
  for (node in (Ntip + 1):total_nodes) {
    # Get children nodes for this internal node from the edge matrix.
    children <- tree$edge[tree$edge[, 1] == node, 2]
    child_states <- states[children]
    # If any child has an NA state or the children disagree, mark as NA.
    if (any(is.na(child_states)) || length(unique(child_states)) > 1) {
      states[node] <- NA
    } else {
      states[node] <- unique(child_states)
    }
  }
  return(states)
}

# Revised version of color_by_habitat using the fast ancestral state calculation.
color_by_trait <- function(col_plot, char_col_name, cl) {
  # Assume that rerooted_tree is your tree and tree_data is a tibble version obtained earlier.
  Ntip <- length(rerooted_tree$tip.label)
  
  # Extract the trait states for the tips from tree_data.
  # This assumes that tree_data is ordered such that the first Ntip rows correspond to tips.
  tip_states <- tree_data[[char_col_name]][1:Ntip]
  # Name the tip states with corresponding tip labels.
  names(tip_states) <- rerooted_tree$tip.label
  # Compute ancestral states for all nodes using the fast bottom-up approach.
  ancestral_states <- fast_get_ancestral_state(rerooted_tree, tip_states)
  # Add the ancestral state information to tree_data.
  tree_data$char_in <- ancestral_states
  
  # Build the plot using ggtree, mapping the computed states onto branches and tips.
  col_plot <- col_plot %<+% tree_data +
    geom_tree(aes(shape = factor(char_in), color = factor(char_in))) +
    geom_tippoint(aes(color = factor(char_in)), size = 3) +
    scale_color_manual(
      values = c("0" = "black", "1" = cl, "NA" = "gray"),
      labels = c("not present", "present", "unknown")
    ) +
    labs(color = char_col_name, shape = char_col_name) +
    theme(
      legend.position = "top",
      legend.title.position = "left",
      legend.justification = c("left", "top"),
      legend.title = element_text(size = 25),
      legend.text = element_text(size = 20),
      legend.background = element_rect(linewidth = 0.35, color = "grey"),
      legend.margin = margin(6, 6, 6, 6),
      legend.box = "vertical"
    )
  
  return(col_plot)
}

# function for generating a tree that has the specified trait mapped onto it
trait_tree <- function(save_path, treeplot, trait, color)
{
  print(paste("=== Generating Trait Tree for ",trait," colored in ",color," ===",sep=""))
  main_plot_trait <- color_by_trait(treeplot, trait, color)
  
  ggsave(save_path,main_plot_trait,device="pdf",width=35,height=75,limitsize=FALSE)
  print(paste("Trait tree saved as",save_path))
}
