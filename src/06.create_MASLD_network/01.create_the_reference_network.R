library(igraph) # 2.1.1
library(OmnipathR) # 3.8.0
library(tidyr) # 1.3.1
library(dplyr) # 1.1.4
library(stringr) # 1.6.0
source("../library.R")

################################################################################
# Description
################################################################################
# This script is used to create the the reference interaction network. It is 
# constructed by integrating curated signalling and metabolic interactions from 
# Reactome, KEGG, HumanCyc, ReconX, provided by Pathway Commons(v12), and OmniPath, 
# accessed using the OmniPathR library (v3.8.0), using the union function of 
# the igraph package. Weakly supported OmniPath interactions were removed, 
# keeping only those where the number of resources and curation efforts were 
# above the 0.5 and 0.75 quantiles of their respective distributions. Dead-end 
# metabolites were also excluded, protein complexes were expanded into fully 
# connected subgraphs, and inhibitory edges were discarded.
# Outputs:
# - The reference network in tsv and RData formats.
################################################################################


################################################################################
# A function to split the complexes and create interactions among its members.
################################################################################
split_complexes <- function(complexes_df) {
  
  final_df <- c() # here the new interactions will be stored
  
  for (i in 1:dim(complexes_df)[1]) {
    row <- complexes_df[i,]
    tmp_interactors <- list()
    tmp_interactions_df <- c()
    
    # check the source if it is a complex 
    if (grepl(pattern = '^COMPLEX', x = row$source)) {
      tmp <- strsplit(x=row$source_genesymbol, split='_')[[1]]
      tmp <- tmp[tmp != ""]
      tmp <- tmp[str_length(tmp) > 1]
      tmp_interactors[['source']] <- tmp
      if (length(tmp_interactors[["source"]]) > 1) {
        comb_df <- as.data.frame(t(combn(x = tmp_interactors[["source"]], m=2)))
        reverse_comb_df <- comb_df[, c('V2', 'V1')]
        names(reverse_comb_df) <- c('V1', 'V2')
        comb_df <- rbind(comb_df, reverse_comb_df)
        colnames(comb_df) <- c('source_genesymbol','target_genesymbol')
        annotation_df <- t(sapply(1:dim(comb_df)[1], function(i){row}))
        annotation_df[,'source_genesymbol'] <- comb_df$source_genesymbol
        annotation_df[,'target_genesymbol'] <- comb_df$target_genesymbol
        annotation_df[,'is_directed'] <- 'complex'
        annotation_df[,'target'] <- annotation_df[,'source']
        tmp_interactions_df <- rbind(tmp_interactions_df, annotation_df)
      }
    } else {
      tmp_interactors[['source']] <- row$source_genesymbol
    }
    
    # check the target if it is a complex 
    if (grepl(pattern = '^COMPLEX', x = row$target)) {
      tmp <- strsplit(x=row$target_genesymbol, split='_')[[1]]
      tmp <- tmp[tmp != ""]
      tmp <- tmp[str_length(tmp) > 1]
      tmp_interactors[['target']] <- tmp
      if (length(tmp_interactors[['target']]) > 1) {
        comb_df <- as.data.frame(t(combn(x = tmp_interactors[["target"]], m=2)))
        reverse_comb_df <- comb_df[, c('V2', 'V1')]
        names(reverse_comb_df) <- c('V1', 'V2')
        comb_df <- rbind(comb_df, reverse_comb_df)
        colnames(comb_df) <- c('source_genesymbol','target_genesymbol')
        annotation_df <- t(sapply(1:dim(comb_df)[1], function(i){row}))
        annotation_df[,'source_genesymbol'] <- comb_df$source_genesymbol
        annotation_df[,'target_genesymbol'] <- comb_df$target_genesymbol
        annotation_df[,'is_directed'] <- 'complex'
        annotation_df[,'source'] <- annotation_df[,'target']
        tmp_interactions_df <- rbind(tmp_interactions_df, annotation_df)
      }
    } else {
      tmp_interactors[['target']] <- row$target_genesymbol
    }
    
    # Create all the combinations between target and source 
    comb_df <- expand.grid(tmp_interactors[["source"]], tmp_interactors[["target"]], 
                           stringsAsFactors = FALSE)
    colnames(comb_df) <- c('source_genesymbol','target_genesymbol')
    annotation_df <- t(sapply(1:dim(comb_df)[1], function(i){row}))
    annotation_df[,'source_genesymbol'] <- comb_df$source_genesymbol
    annotation_df[,'target_genesymbol'] <- comb_df$target_genesymbol
    #annotation_df[,'is_directed'] <- 'interaction_contains_complex'
    tmp_interactions_df <- rbind(tmp_interactions_df, annotation_df)
    
    # Add the created interactions in the final data frame
    final_df <- rbind(final_df, tmp_interactions_df)
  }
  duplicates <- which(duplicated(final_df))
  final_df <- final_df[-duplicates,]
  return(final_df)
}


################################################################################
# A function to modify some fields in the data frame of OmniPath.
################################################################################
omnipath_interaction_information <- function(df) {
  output <- c()
  for (i in 1:dim(df)[1]) {
    row <- df[i,]
    s = ''
    if (row$is_directed == 'complex') {
      s <- as.character(row$is_directed)
    } else {
      if (row$is_directed == 1) {
        s <- 'directed'
      }
      else if (row$is_directed == 0) {
        s <- 'undirected'
      }
      if (row$is_stimulation == 1) {
        s <- paste(s, '_stimulation', sep='')
      }
      if (row$is_inhibition == 1) {
        s <- paste(s, '_inhibition', sep='')
      }
    } # else
    output <- c(output, s)
  } # for
  return(output)
}


################################################################################
# A function to read the tab-delimited files of interaction databases. These 
# files have been extract from PathwayCommons, so they have the same format.
################################################################################
read_sif_file <- function(file, database_name) {
  df <- read.csv(file, sep='\t', header = FALSE)
  colnames(df) <- c('source', database_name, 'target')
  df <- df[,c('source', 'target', database_name)]
  network <- igraph::graph_from_data_frame(df)
  return(network)
}


################################################################################
# A function to remove metabolites which do not function as bridges on the 
# network (so their out degree is equal to zero).
################################################################################
remove_chebi_nodes <- function(network) {
  s = 10
  while (s > 0) {
    degrees_df <- data.frame(degree = degree(network, mode = 'out'))
    indexes <- which(grepl(pattern='CHEBI', row.names(degrees_df)))
    chebi_degrees_df = degrees_df[indexes,,drop=FALSE]
    chebi_nodes <- rownames(chebi_degrees_df[chebi_degrees_df$degree == 0,,drop=FALSE])
    network <- delete_vertices(network, chebi_nodes)
    s = length(chebi_nodes)
  }
  return(network)
}


################################################################################
# A function to remove inhibitory interactions, as their effect cannot be modeled
# with network propagation and maxflow.
################################################################################
remove_inhibitions <- function(network){
  network <- reference_network
  network_df <- igraph::as_data_frame(reference_network)
  network_df$edge_id <- rownames(network_df)
  tmp_df <- apply(network_df, c(1,2), function(x) {grepl('inhibit', x)})
  inhibition_interactions <- which(rowSums(tmp_df) > 0)
  to_remove <- as.integer(network_df[inhibition_interactions,'edge_id'])
  network <- delete_edges(network, edges = to_remove)
  return(network)
}


################################################################################
# Inputs
################################################################################
output_dir = '../../results/ucam_sanyal/networks/'
files <- list(
  'KEGG' = '../../data/interaction_databases/kegg.tsv',
  'ReconX' = '../../data/interaction_databases/reconx.tsv',
  'HumanCyc' = '../../data/interaction_databases/humancyc.tsv',
  'Reactome' = '../../data/interaction_databases/reactome.tsv'
)

################################################################################
# 1. Loading of the interaction databases (sif files derived from PathwayCommons) 
# The goal is to create the respective networks (igraph objects) and save them 
# in a list ('networks_list').
################################################################################
networks_list <- list()
for (name in names(files)) {
  file <- files[[name]]
  network <- read_sif_file(file, name)
  networks_list[[name]] <- network
}
# > lengths(networks_list)
# KEGG   ReconX HumanCyc Reactome 
# 2282     2708     3309    11069 

################################################################################
# 2. Loading of the Omnipath database + filtering + modification of complexes.
# Finally the generated igraph object will be save in the 'networks_list', like
# the previously defined interaction networks.
################################################################################
# A. Loading
interactions <- omnipath_interactions() %>% as_tibble()

# B. Filtering of interaction based on n_resources and curation_effort
#interactions <- interactions[(interactions$n_resources >= 1) & (interactions$curation_effort >= 1),]
q1 <- quantile(interactions$n_resources, probs = c(0.5))[1]
q2 <- quantile(interactions$curation_effort, probs = c(0.75))[1]
interactions <- interactions[(interactions$n_resources >= q1) & (interactions$curation_effort >= q2),]

# C. Finding of interactions which contain complexes
c1 <- which(grepl(pattern = '^COMPLEX', x = interactions$source))
c2 <- which(grepl(pattern = '^COMPLEX', x = interactions$target))
complexes <- unique(c(c2, c1))

# D. Split of complexes
complex_interactions <- split_complexes(data.frame(interactions[complexes,]))

# E. Further filtering of  network interactions
interactions <- data.frame(interactions)[-complexes,]
interactions <- rbind(interactions, complex_interactions)
rownames(interactions) <- NULL

# F. Modification of the interactions information
s <- omnipath_interaction_information(interactions)
interactions <- interactions[,c('source_genesymbol', 'target_genesymbol')]
colnames(interactions) <- c('source', 'target')
interactions[,'OmniPath'] <- s
duplicates <- which(duplicated(interactions))
interactions <- interactions[-duplicates,]
rownames(interactions) <- NULL

# G. Data frame to igraph object (network generation)
network <- igraph::graph_from_data_frame(interactions)
networks_list[['OmniPath']] <- network
# > lengths(networks_list)
# KEGG   ReconX HumanCyc Reactome OmniPath 
# 2282     2708     3309    11069     3600 



################################################################################
# 3. Creation of the reference super-network
################################################################################
# A. Unification of the five networks
reference_network <- igraph::union(networks_list[["KEGG"]], 
                                   networks_list[["ReconX"]],
                                   networks_list[["HumanCyc"]],
                                   networks_list[['Reactome']],
                                   networks_list[['OmniPath']])

# B. Removal of duplicated edges (annotation in more than one database)
reference_network_edges <- igraph::as_data_frame(reference_network)
reference_network_edges <- data.frame(
  reference_network_edges %>% 
    group_by(from, to) %>%
    mutate(KEGG=paste0(na.omit(KEGG), collapse=' & '),
           ReconX=paste0(na.omit(ReconX), collapse=' & '),
           HumanCyc=paste0(na.omit(HumanCyc), collapse=' & '),
           Reactome=paste0(na.omit(Reactome), collapse=' & '),
           OmniPath=paste0(na.omit(OmniPath), collapse=' & ')) %>%
    distinct(from, to, .keep_all = TRUE)
)
reference_network_edges[reference_network_edges == ""] = NA
reference_network <- igraph::graph_from_data_frame(reference_network_edges)
paste('Ref network - integration:', length(V(reference_network)), 'nodes', sep=' ')
paste('Ref network - integration:', length(E(reference_network)), 'edges', sep=' ')
      
# C. Removal of CHEBI nodes, with out degree equal to zero
reference_network <- remove_chebi_nodes(reference_network)

# D. Removal of inhibitions (they cannot be modeled in network propagation & maxflow)
reference_network <- remove_inhibitions(reference_network)

# E. Detection of network components in order to keep the biggest one. We expect
# that the vast majority of nodes are connected in one dominant component.
detected_components <- components(reference_network)
biggest_component <- which.max(detected_components$csize)
nodes <- V(reference_network)[detected_components$membership == biggest_component]
reference_network <- igraph::induced_subgraph(reference_network, nodes)
paste('Ref network - chebi/inhibitions removal:', length(V(reference_network)), 'nodes', sep=' ')
paste('Ref network - chebi/inhibitions:', length(E(reference_network)), 'edges', sep=' ')

# F. Saving of the reference network in a tsv file
write.table(x = igraph::as_data_frame(reference_network), row.names = TRUE, sep='\t',
            file = paste(output_dir, 'reference_network.tsv', sep=''))
reference_network <- reference_network
save(reference_network, file=paste(output_dir, "reference_network.RData", sep=''))


length(V(reference_network))
#12549
length(E(reference_network))
#205706
is_connected(reference_network)
#TRUE
