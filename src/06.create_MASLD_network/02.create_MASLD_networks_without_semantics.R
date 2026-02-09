suppressMessages(library(dplyr)) # 1.1.4
suppressMessages(library(tidyr)) # 1.3.1
suppressMessages(library(OmnipathR)) # 3.8.0
suppressMessages(library(PCSF)) # 0.99.1
suppressMessages(library(igraph)) # 2.2.1
suppressMessages(library(GOSemSim)) # 2.29.2

################################################################################
# Description
################################################################################
# This long script generates components of the MASLD network, using as framework 
# the reference one. The results from WGCNA, differential expression and 
# TF-activity analysis, TF and Reactome enrichment analysis are combined in order 
# to filter the reference network and keep those parts which are highly 
# associated with MASLD variables (Steatosis, Balloning, Fibrosis, NAS score) and
# the molecular differentiations along the trajectory of sliding windows.
# Output files:
#
# A new folder will be created to store the first draft version of MASLD upstream
# and downstream networks (MASLD_network_construction). The following files
# will be included there:
# 
# - downstream_gene_sets_for_MASLD_networks.RData
# - upstream_gene_sets_for_MASLD_networks.RData
# - MASLD_directed_downstream_networks_v0.RData
# - MASLD_directed_upstream_networks_v0.RData
#
# A tab-delimited file will be created to include all the genes (nodes) in 
# these networks: MASLD_nodes_without_semantics.csv
#
# Also, another file which contains the valid (annotated) entrez gene ids
# (regardless their inclusion in MASLD networks) will be created 
# (valid_entrez_ids.csv). This file, as well as that of MASLD nodes, will be 
# used in next Python script to find the annotated gene symbols and their 
# respective entrez ids for which the semantic similarity scores will be 
# calculated.
#
################################################################################


################################################################################
# Inputs
################################################################################
modules_dir = '../../results/ucam_sanyal/wgcna_and_linear_modelling/grid_params/4_180/'
tf_analysis_dir = '../../results/ucam_sanyal/de_and_tf_analysis/'
enrichment_analysis_dir = '../../results/ucam_sanyal/enrichment_analysis/'
output_dir = '../../results/ucam_sanyal/networks/'
tmp_output_dir = '../../results/ucam_sanyal/networks/MASLD_network_construction/'
dir.create(tmp_output_dir, showWarnings = FALSE, recursive = TRUE)
load('../../data/reference_network.RData')
Regulon_file<- read.csv("../../data/collectTRI_network.tsv", header=T, sep='\t')
load("../../data/annotation_databases/annotation_lists.RData")
reactome_terms_df <- read.csv2('../../data/reactome_graph.tsv', sep='\t', header=TRUE)
signalling_terms_df <- read.csv2('../../data/reactome_signalling_terms.tsv', sep='\t', header=TRUE)


################################################################################
# 1. Creation of a list for the module gene sets of MASLD variables. Initially,
# the gene co-expression (GE) module gene sets and their coefficients 
# from the linear models are loaded from the respective tab-delimited files. 
# The most significant GE modules (FDR < 0.05) for each MASLD variable are 
# selected and their gene sets are saved into the "module_sets_per_var" list.
################################################################################

# A. Loading of the GE module gene sets.
modules_df <- read.table(paste(modules_dir, 'filtered_modules.tsv',sep=''), 
                         sep='\t', header=TRUE)
modules_gene_sets <- split(x=modules_df$gene_symbol, f=modules_df$module_color)

# B. Loading of the coefficients from lm(variable ~ modules) models. Statistical 
# significance at 0.05 (removal of age & inflammation from the independent variables).
filename <- paste(modules_dir, 'module_variable_coefficients_FINAL.tsv', sep='')
modules_per_var_df <- read.table(filename, sep='\t', header=TRUE)
modules_per_var_df <- modules_per_var_df[modules_per_var_df$p.value_adjust < 0.05,]
modules_per_var_df <- modules_per_var_df[!modules_per_var_df$term %in% c('age', 'inflammation'),] 
modules_per_var <- sapply(unique(modules_per_var_df$variable), function(var) {
  modules_per_var_df[modules_per_var_df$variable == var, 'term']
}, USE.NAMES = TRUE)

# C. Creation of the "module_sets_per_var" list to save the GE module gene sets.
module_sets_per_var <- sapply(unique(modules_per_var_df$variable), function(var) {
  var_modules <- modules_per_var[[var]]
  gene_set <- c()
  for (module in var_modules) {
    gene_set <- c(gene_set, modules_gene_sets[[module]])
  }
  unique(gene_set)
}, USE.NAMES = TRUE)
print('length of filtered module gene sets:')
print(lengths(module_sets_per_var))
#Fibrosis  Steatosis Ballooning        NAS 
#3138       2540       3820       3377 


################################################################################
# 2. Similarly as above, a gene set related to the differentially activated TFs 
# (FDR<0.05) is created and saved in the 'module_set_for_diff_activated_tfs' 
# variable. Initially the TFs are loaded from the results of msViper and their
# regulons are retrieved from the collectTri database. They are filtered to keep
# gene that are in the GE modules, and finally all of them are unified into a 
# vector of genes.
################################################################################

# A. A unified vector for all genes in the significant modules.
module_genes <- unique(unlist(module_sets_per_var, use.names = FALSE))
length(module_genes)
# 6336

# B. Loading of the results of TF activity analysis (collectTri & msviper) and
# filtering with FDR < 0.05.
filename <- paste(tf_analysis_dir, 'tf_results.tsv', sep='')
diff_tf_df <- read.table(filename, sep='\t', header=TRUE, row.names = 1)
cols <- colnames(diff_tf_df)[grepl(pattern = 'fdr', colnames(diff_tf_df))]
fdr_tf_df <- diff_tf_df[,cols]
scores <- rowSums(fdr_tf_df < 0.05)
diff_activated_tfs <- names(scores[scores > 0]) # FDR < 0.05 at least in one SW
print('length of diff activated tfs:')
print(length(diff_activated_tfs))
# 141

# C. Loading of the collectTri regulons to select those of differentially 
# activated TFs. This gene set is filtered with the 'module_genes' to keep only 
# those genes which are both in the regulons of diff TFs and the significant 
# GE modules.
Regulon_file <- Regulon_file[Regulon_file$tf %in% diff_activated_tfs,]
Regulon_file <- Regulon_file[Regulon_file$target %in% module_genes,]
targets <- unique(Regulon_file$target)
module_set_for_diff_activated_tfs <- targets
print('length of filtered diff_exp gene sets:')
print(length(module_set_for_diff_activated_tfs))
# 1905


################################################################################
# 3. In this step two lists are created to store the TFs and genes in related 
# signalling pathways for each MASLD variable. Initially, the 'pathways_per_var'
# is create to save the results of pathway analysis per MASLD variable. Then, 
# these pathways are filtered to keep only the semantically unique ones (i.e. 
# pathways which do not have any descendant pathway in the results). For these 
# pathways, the related TFs are retrieved from the 'genes' columns and are 
# saved in the "tf_sets_per_var" list. Then the pathways are further 
# filtered to keep only those which are relevant to signalling processes. 
# Finally, the gene sets of the selected pathways (no only the TFs) are retrieved 
# from the global annotation of Reactome to create the "reactome_sets_per_var" 
# list. 
################################################################################

# A. Loading of data related to Reactome terms:
# - Their genomic annotation
# - Their parent-child relations
# - The list of pathways which are related to signalling processes (this is 
# a customized list)
reactome_gene_sets <- annotation_lists[["Reactome_2024"]]


# B. Loading of the Reactome enrichment analysis results for the enriched TFs
# and separation based on the MASLD variables.
filename <- paste(enrichment_analysis_dir, 'reactome_enrichment_for_enriched_tfs.tsv', sep='')
enrichment_df <- read.table(filename, sep='\t', header=TRUE)
pathways_per_var <- sapply(unique(enrichment_df$variable), function(var) {
    enrichment_df[enrichment_df$variable == var,]
}, USE.NAMES = TRUE, simplify = FALSE)

# C. Filtering of the enriched pathways based on their semantic specificity.
unique_pathways_per_var <- sapply(names(pathways_per_var), function(var) {
    terms_df <- pathways_per_var[[var]]
    tmp <- c()
    for (term_id in terms_df$term_id) {
        definition <- terms_df[terms_df$term_id == term_id, 'term']
        descendants <- reactome_terms_df[reactome_terms_df$term_id == term_id, "descendants"]
        descendants <- unlist(strsplit(descendants, ','))
        if (length(descendants) == 0) {
          tmp <- c(tmp, definition)
        } else {
          descendants <- reactome_terms_df[reactome_terms_df$id %in% descendants, 'term_id']
          if (length(intersect(terms_df$term_id, descendants)) == 0) {
            tmp <- c(tmp, definition)
          } else {
            # pass
          }
        } 
    }
    tmp
}, USE.NAMES = TRUE)
print('reactome pathways per var (semantic filtering):')
print(lengths(unique_pathways_per_var))
# Fibrosis Steatosis Ballooning  NAS 
# 129         9         132         46 

# D. Selection of TFs which are associated with the semantically unique 
# pathways. Then these pathways are filtered to keep only those which are 
# related to signalling and their gene sets are retieved from the global 
# annotation of Reactome.
tf_sets_per_var <- list()
reactome_sets_per_var <- list()
for (var in names(unique_pathways_per_var)) {
  # Get the semantically unique terms
  terms <- unique_pathways_per_var[[var]]
  # Find and save the set of related TFs
  tmp_tfs <- c()
  for (pathway in terms){
    row <- enrichment_df[(enrichment_df$term == pathway & enrichment_df$variable == var),]
    tmp_tfs <- c(tmp_tfs, strsplit(split=';', x=row$genes)[[1]])
  }
  tf_sets_per_var[[var]] <- unique(tmp_tfs)
  # Find the semantically unique terms which are related to signalling
  var_enrichment_df <- pathways_per_var[[var]]
  var_enrichment_df <- var_enrichment_df[var_enrichment_df$term %in% terms,]
  var_enrichment_df <- var_enrichment_df[var_enrichment_df$term_id %in% signalling_terms_df$term_id,]
  # Get their gene sets and save them
  tmp <- c()
  for (pathway in var_enrichment_df$term) {
    tmp <- c(tmp, reactome_gene_sets[[pathway]])
  }
  reactome_sets_per_var[[var]] <- unique(tmp)
}
print('length of tfs sets:')
print(lengths(tf_sets_per_var))
# Fibrosis Steatosis Ballooning  NAS 
# 114        24         127         48 

print('length of reactome gene sets (signalling):')
print(lengths(reactome_sets_per_var))
#Fibrosis Steatosis Ballooning  NAS 
# 1074        237        1115   400


################################################################################
# 4. Similarly to the previous step, a vector of genes, related to the 
# signalling pathways which are enriched in the differentially activated TFs 
# (collectTRI-msViper), is created. The results from the enrichment analysis are
# loaded and filtered to keep only the pathways related to signalling processes. 
# Then, the semantic filtering is performed and finally the vector 
# "reactome_set_for_diff_activated_tfs" is created, including the genes 
# of semantically unique and siganlling-related pathways.
################################################################################

# A. Loading of the results from Reactome enrichment analysis for the differentially
# activated TFs to get the terms, which are related to signalling and they are 
# semantically unique ('diff_activated_tfs_pathways').
filename <- paste(enrichment_analysis_dir, 'reactome_enrichment_for_diff_activated_tfs.tsv', sep='')
enrichment_df <- read.table(filename, sep='\t', header=TRUE)
enrichment_df <- enrichment_df[enrichment_df$term_id %in% signalling_terms_df$term_id,]
diff_activated_tfs_pathways <- c()
for (term_id in enrichment_df$term_id) {
  definition <- enrichment_df[enrichment_df$term_id == term_id, 'term']
  descendants <- reactome_terms_df[reactome_terms_df$term_id == term_id, "descendants"]
  descendants <- unlist(strsplit(descendants, ','))
  if (length(descendants) == 0) {
    diff_activated_tfs_pathways <- c(diff_activated_tfs_pathways, definition)
  } else {
    descendants <- reactome_terms_df[reactome_terms_df$id %in% descendants, 'term_id']
    if (length(intersect(enrichment_df$term_id, descendants)) == 0) {
      diff_activated_tfs_pathways <- c(diff_activated_tfs_pathways, definition)
    } else {
      # pass
    }
  } 
}

# B. Definition of the global gene set for the above pathways
# ('reactome_set_for_diff_activated_tfs').
tmp <- c()
for (pathway in diff_activated_tfs_pathways) {
  tmp <- c(tmp, reactome_gene_sets[[pathway]])
}
reactome_set_for_diff_activated_tfs <- unique(tmp)
print('length of reactome gene sets (signalling) for the diff activated TFs:')
print(length(reactome_set_for_diff_activated_tfs))
# 909


################################################################################
# 5. Loading and transformation of reference PPI network into the appropriate 
# format, using the 'construct_interactome' function of the PCSF package. This 
# object will be used as the template for the implementation of 'PCSF_rand'.
################################################################################
base_network_df <- igraph::as_data_frame(reference_network) # unified_network is directed
base_network_df <- base_network_df[,c(1,2)] 
base_network_df[,'cost'] <- 0.001
base_ppi <- construct_interactome(base_network_df) # it is undirected
# Also the 'construct_interactome' function removes duplicated edges in the 
# undirected graph, which have been produced due to the bi-directionality of 
# some interactions, like those of complexes.


################################################################################
# 6. Creation of an upstream network for each MASLD variable, based on 
# the reference one and the PCSF algorithm. Before the construction of these 
# networks, the list 'upstream_gene_sets_per_var' is created to store the seed 
# nodes for each MASLD variable. Specifically, the variable-related TFs and 
# reactome gene sets (defined in the previous steps) are used. These seed sets 
# are used as terminal nodes in the PCSF algorithm to define the 
# upstream networks. A similar network is created for the diff_exp condition.
################################################################################

# A. Creation of the 'upstream_gene_sets_per_var' which contains the TFs and 
# Reactome (signalling-related) gene sets for each MASLD variable.
upstream_gene_sets_per_var <- sapply(names(reactome_sets_per_var), function(var){
  g1 <- reactome_sets_per_var[[var]]
  g2 <- tf_sets_per_var[[var]]
  #g3 <- module_sets_per_var[[var]]
  unique(c(g1, g2))
}, USE.NAMES = TRUE)
upstream_gene_sets_per_var[['diff_exp']] <- unique(c(diff_activated_tfs,
                                                     reactome_set_for_diff_activated_tfs))
print('length of upstream network gene sets:')
lengths(upstream_gene_sets_per_var)
# Fibrosis  Steatosis Ballooning        NAS   diff_exp 
# 1109      250       1160              420   986 
print('length of the unified upstream network gene set:')
length(unique(unlist(upstream_gene_sets_per_var)))
# 1237

# B. Saving of the upstream gene sets (before the creation of the networks).
filename <- paste(tmp_output_dir, "upstream_gene_sets_for_MASLD_networks.RData", sep='')
save(upstream_gene_sets_per_var, file=filename)

# C. Construction of the upstream networks with RCSF_rand. We used the
# randomization technique to extract a more robust solution. However, we use only
# n=10, as there is no difference in nodes scores and edge weights, which
# means that the different solutions will not have significant differences. We
# want to assure that all the important nodes (seeds) will be in the final 
# solution, that's why we run the algorithm 10 times (the final network is the
# union of these 10 solutions).
upstream_networks_per_var <- sapply(names(upstream_gene_sets_per_var), function(var) {
  seeds_vector <- as.vector(upstream_gene_sets_per_var[[var]])
  seeds_scores <- rep(1, length(seeds_vector))
  names(seeds_scores) <- seeds_vector
  # PCSF_rand to get the network - it is undirected.
  network <- PCSF_rand(base_ppi, terminals=seeds_scores, mu=0.01, n=10, r=0.1)
  #network <- PCSF(base_ppi, terminals=seeds_scores, mu=5e-04)
  # Finding of chebi nodes and remove them - we do not need them in the 
  # upstream networks.
  indexes <- which(grepl(pattern='CHEBI', as_ids(V(network))))
  chebi_nodes <- as_ids(V(network))[indexes]
  network <- delete_vertices(network, chebi_nodes)
  # Detection of network components to keep the biggest one (probably all the 
  # other are pairs or small groups so they should be be removed).
  network_components <- components(network)
  main_component <- which.max(network_components$csize)
  to_remove <- network_components$membership[network_components$membership != main_component]
  network <- delete_vertices(network, to_remove)
  # Definition of the part of reference network with the selected nodes. This 
  # network is directed.
  subnetwork <- induced_subgraph(reference_network, vids=as_ids(V(network)), 
                                 impl='create_from_scratch')
  print(paste(var, length(V(subnetwork))))
  subnetwork
}, USE.NAMES = TRUE)

# D. Definition of some vertex attributes which indicate the network where 
# each node is included.
for (name in names(upstream_networks_per_var)) {
  tmp <- upstream_networks_per_var[[name]]
  tmp <- set_vertex_attr(tmp, name=paste('upstream_',name,sep=''), value=TRUE)
  upstream_networks_per_var[[name]] <- tmp
}

# E. Saving of the upstream networks.
filename <- paste(tmp_output_dir, "MASLD_directed_upstream_networks_v0.RData", sep='')
save(upstream_networks_per_var, file=filename)


################################################################################
# 7. Creation of a downstream network for each MASLD variable, based on 
# the reference one and the PCSF algorithm. Before the construction of these 
# networks, the list "downstream_gene_sets_per_var" is created to store the seed 
# nodes for each MASLD variable. Specifically, the variable-related TFs and 
# module gene sets are used for that task. These sets are used as terminal nodes 
# in the PCSF algorithm to define the respective downstream networks.
# A similar network is created for the diff_exp condition.
################################################################################

# A. Creation of the 'downstream_gene_sets_per_var' which contains the TFs and 
# module gene sets for each MASLD variable.
downstream_gene_sets_per_var <- sapply(names(module_sets_per_var), function(var){
  g2 <- tf_sets_per_var[[var]]
  g3 <- module_sets_per_var[[var]]
  unique(c(g2, g3))
}, USE.NAMES = TRUE)
downstream_gene_sets_per_var[['diff_exp']] <- unique(c(diff_activated_tfs,
                                                       module_set_for_diff_activated_tfs))
print('length of downstream network gene sets:')
lengths(downstream_gene_sets_per_var)
# Fibrosis  Steatosis Ballooning        NAS   diff_exp 
# 3224      2561      3912              3413  1998 
print('length of the unified downstream network gene set:')
length(unique(unlist(downstream_gene_sets_per_var)))
# 6464

# B. Saving of the downstream gene sets (before the creation of the networks).
filename <- paste(tmp_output_dir, "downstream_gene_sets_for_MASLD_networks.RData", sep='')
save(downstream_gene_sets_per_var, file=filename)

# C. Construction of the downstream networks with RCSF_rand. We used the
# randomization technique to extract a more robust solution. However, we use only
# n=10, as there is not difference in nodes scores and edge weights, which
# means that the different solutions will not have significant differences. We
# want to assure that all the important nodes (seeds) will be in the final 
# solution, that's why we run the algorithm 10 times (the final network is the
# union of these 10 solutions).
downstream_networks_per_var <- sapply(names(downstream_gene_sets_per_var), function(var) {
  seeds_vector <- as.vector(downstream_gene_sets_per_var[[var]])
  seeds_scores <- rep(1, length(seeds_vector))
  names(seeds_scores) <- seeds_vector
  # Run PCSF_rand to get the network - it is undirected
  network <- PCSF_rand(base_ppi, terminals=seeds_scores, mu=0.01, n=10, r=0.1)
  #network <- PCSF(base_ppi, terminals=seeds_scores, mu=5e-04)
  # Finding of chebi nodes and remove them - we do not need them in the 
  # upstream networks.
  #indexes <- which(grepl(pattern='CHEBI', as_ids(V(network))))
  #chebi_nodes <- as_ids(V(network))[indexes]
  #network <- delete_vertices(network, chebi_nodes)
  # Detection of network components to keep the biggest one (probably all the 
  # other are pairs or small groups so they should be be removed).
  network_components <- components(network)
  main_component <- which.max(network_components$csize)
  to_remove <- network_components$membership[network_components$membership != main_component]
  network <- delete_vertices(network, names(to_remove))
  # Definition of the part of reference network which includes the selected nodes. 
  # This network is directed.
  subnetwork <- induced_subgraph(reference_network, vids=as_ids(V(network)), 
                                 impl='create_from_scratch')
  print(paste(var, length(V(subnetwork))))
  subnetwork
}, USE.NAMES = TRUE)

# D. Definition of some vertex attributes which indicate the network where each 
# node is included.
for (name in names(downstream_networks_per_var)) {
  tmp <- downstream_networks_per_var[[name]]
  tmp <- set_vertex_attr(tmp, name=paste('downstream_',name,sep=''), value=TRUE)
  downstream_networks_per_var[[name]] <- tmp
}

# E. Saving of the downstream networks.
filename <- paste(tmp_output_dir, "MASLD_directed_downstream_networks_v0.RData", sep='')
save(downstream_networks_per_var, file=filename)


################################################################################
# 8. This final step is used to gather together all the genes in upstream and
# downstream networks. This vector will be used to calculate the pairwise
# semantic similarities later.
################################################################################

# A. Gene sets of all upstream networks.
upstream_nodes <- c()
for (var in names(upstream_networks_per_var)) {
  nodes <- as_ids(V(upstream_networks_per_var[[var]]))
  upstream_nodes <- c(upstream_nodes, nodes)
}
upstream_nodes <- unique(upstream_nodes)

# B. Gene sets of all downstream networks.
downstream_nodes <- c()
for (var in names(downstream_networks_per_var)) {
  nodes <- as_ids(V(downstream_networks_per_var[[var]]))
  downstream_nodes <- c(downstream_nodes, nodes)
}
downstream_nodes <- unique(downstream_nodes)

# C. Unified vector of both upstream and downstream genes.
nodes <- unique(c(upstream_nodes, downstream_nodes))
print('length of the unified gene set:')
length(nodes)
# 4294
nodes_df = data.frame(gene_symbol = nodes)
filename <- paste(output_dir, "MASLD_nodes_without_semantics.csv", sep='')
write.table(x = nodes_df, file = filename, sep=',')

# D. Finding and saving of genes  which are annotated with GO BP terms. Only the
# entrez ids are saved here. In the next script, gene names and entrez ids will
# be combined to create the 'MASLD_nodes_without_semantics_filtered.csv' file,
# which will be used for the calculation of semantic similarity scores.
hsGO <- godata('org.Hs.eg.db', ont='BP')
entrez_ids <- hsGO@geneAnno$ENTREZID
entrez_ids <- unique(entrez_ids)
entrez_df = data.frame(entrez_id = entrez_ids)
write.table(x = entrez_df, file = paste(output_dir, "valid_entrez_ids.csv", sep=''),
            sep=',')

