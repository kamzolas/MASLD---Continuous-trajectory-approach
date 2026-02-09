library(igraph) # 2.2.1
library(tidyr) # 1.3.1
library(dplyr) # 1.1.4


################################################################################
# Description
################################################################################
# This script integrates the interactions from upstream and downstream networks
# into a unique super-network (in both directed and undirected formats). The 
# edge and node attributes keep the information of the initial subnetworks from
# which they have been retrieved. For example one interaction could have been 
# retrieved from KEGG and ReconX, and it was included in the downstream network 
# of Fibrosis. This kind of information is still available in the unified 
# networks to keep track for all the elements.
# Output files:
# - There are three outputs (two networks in RData and a csv with network nodes):
#   - MASLD_unified_directed_network_with_semantics.RData
#   - MASLD_unified_undirected_network_with_semantics.RData
#   - MASLD_nodes_with_semantics.csv
################################################################################


################################################################################
# Inputs
################################################################################
main_dir = '../../results/ucam_sanyal/networks/'


################################################################################
# 1. Unification of downstream MASLD networks
################################################################################

# A. Loading of networks and unification with igraph::union
filename <- paste(main_dir, 'MASLD_network_construction/MASLD_directed_downstream_networks_v1.RData', sep='')
load(filename)
downstream_unified_network <- igraph::union(downstream_networks_per_var[["Fibrosis"]],
                                            downstream_networks_per_var[["NAS"]],
                                            downstream_networks_per_var[["Steatosis"]],
                                            downstream_networks_per_var[["Ballooning"]],
                                            downstream_networks_per_var[["diff_exp"]])

# B. Definition of node attributes to keep as information their initial networks.
tmp <- downstream_unified_network
tmp <- set_vertex_attr(tmp, name='downstream', value=TRUE)
attributes <- get.vertex.attribute(tmp)
tmp <- set_vertex_attr(tmp, name='downstream_Fibrosis', value=FALSE, 
                       index=which(is.na(attributes$downstream_Fibrosis)))
tmp <- set_vertex_attr(tmp, name='downstream_NAS', value=FALSE, 
                       index=which(is.na(attributes$downstream_NAS)))
tmp <- set_vertex_attr(tmp, name='downstream_Steatosis', value=FALSE, 
                       index=which(is.na(attributes$downstream_Steatosis)))
tmp <- set_vertex_attr(tmp, name='downstream_Ballooning', value=FALSE, 
                       index=which(is.na(attributes$downstream_Ballooning)))
tmp <- set_vertex_attr(tmp, name='downstream_diff_exp', value=FALSE, 
                       index=which(is.na(attributes$downstream_diff_exp)))
downstream_unified_network <- tmp

# C. Removal of duplicated edge attributes and creation of new ones, which 
# indicate the database from which an interaction was retrieved, as well as the 
# type of that interaction.
attrs <- edge_attr_names(downstream_unified_network)
cols <- c('KEGG', 'ReconX', 'HumanCyc', 'Reactome', 'OmniPath', 'weight')
for (col in cols) {
  col_attrs <- attrs[grep(col, attrs)]
  tmp_df <- c()
  for (attr in col_attrs) {
    tmp_df <- cbind(tmp_df, get.edge.attribute(downstream_unified_network, attr))
  }
  tmp_df <- as.data.frame(tmp_df)
  tmp_df2 <- data.frame(id = seq(1:dim(tmp_df)[1]))
  tmp_df2$count <- rowSums(is.na(tmp_df))
  new_col <- c()
  for (i in 1:dim(tmp_df2)[1]) {
    if (tmp_df2[i, 'count'] == 5) {
      s <- NA
    } else if (col == 'weight') {
      s <- unique(na.exclude(as.numeric(as.vector(tmp_df[i,]))))
    } else {
      s <- unique(na.exclude(as.character(as.vector(tmp_df[i,]))))
    }
    new_col <- c(new_col, s)
  }
  tmp_df2[,col] <- new_col
  for (attr in col_attrs) {
    downstream_unified_network <- delete_edge_attr(downstream_unified_network, attr)
  }
  downstream_unified_network <- set_edge_attr(downstream_unified_network, name=col,
                                              value=tmp_df2[,col])
}


################################################################################
# 2: Unification of upstream MASLD networks
################################################################################

# A. Loading of networks and unification with igraph::union
filename <- paste(main_dir, 'MASLD_network_construction/MASLD_directed_upstream_networks_v1.RData', sep='')
load(filename)
upstream_unified_network <- igraph::union(upstream_networks_per_var[["Fibrosis"]],
                                          upstream_networks_per_var[["NAS"]],
                                          upstream_networks_per_var[["Steatosis"]],
                                          upstream_networks_per_var[["Ballooning"]],
                                          upstream_networks_per_var[["diff_exp"]])

# B. Definition of node attributes to keep as information their initial networks.
tmp <- upstream_unified_network
tmp <- set_vertex_attr(tmp, name='upstream', value=TRUE)
attributes <- get.vertex.attribute(tmp)
tmp <- set_vertex_attr(tmp, name='upstream_Fibrosis', value=FALSE, 
                       index=which(is.na(attributes$upstream_Fibrosis)))
tmp <- set_vertex_attr(tmp, name='upstream_NAS', value=FALSE, 
                       index=which(is.na(attributes$upstream_NAS)))
tmp <- set_vertex_attr(tmp, name='upstream_Steatosis', value=FALSE, 
                       index=which(is.na(attributes$upstream_Steatosis)))
tmp <- set_vertex_attr(tmp, name='upstream_Ballooning', value=FALSE, 
                       index=which(is.na(attributes$upstream_Ballooning)))
tmp <- set_vertex_attr(tmp, name='upstream_diff_exp', value=FALSE, 
                       index=which(is.na(attributes$upstream_diff_exp)))
upstream_unified_network <- tmp

# C. Removal of duplicated edge attributes and creation of new ones, which 
# indicate the database from which an interaction was retrieved, as well as the 
# type of that interaction.
attrs <- edge_attr_names(upstream_unified_network)
cols <- c('KEGG', 'ReconX', 'HumanCyc', 'Reactome', 'OmniPath', 'weight')
for (col in cols) {
  col_attrs <- attrs[grep(col, attrs)]
  tmp_df <- c()
  for (attr in col_attrs) {
    tmp_df <- cbind(tmp_df, get.edge.attribute(upstream_unified_network, attr))
  }
  tmp_df <- as.data.frame(tmp_df)
  tmp_df2 <- data.frame(id = seq(1:dim(tmp_df)[1]))
  tmp_df2$count <- rowSums(is.na(tmp_df))
  new_col <- c()
  for (i in 1:dim(tmp_df2)[1]) {
    if (tmp_df2[i, 'count'] == 5) {
      s <- NA
    } else if (col == 'weight') {
      s <- unique(na.exclude(as.numeric(as.vector(tmp_df[i,]))))
    } else {
      s <- unique(na.exclude(as.character(as.vector(tmp_df[i,]))))
    }
    new_col <- c(new_col, s)
  }
  tmp_df2[,col] <- new_col
  for (attr in col_attrs) {
    upstream_unified_network <- delete_edge_attr(upstream_unified_network, attr)
  }
  upstream_unified_network <- set_edge_attr(upstream_unified_network, name=col,
                                            value=tmp_df2[,col])
}

print('length of the unified downstream network:')
print(length(V(downstream_unified_network)))
# 3201
print('length of the unified upstream network:')
print(length(V(upstream_unified_network)))
# 972


################################################################################
# 3. Construction of the unified network.
################################################################################

# A. Creation of  directed MASLD network.
unified_directed_network <- igraph::union(upstream_unified_network, 
                                          downstream_unified_network)

# B. Modification of edge attributes (as above).
attrs <- edge_attr_names(unified_directed_network)
cols <- c('KEGG', 'ReconX', 'HumanCyc', 'Reactome', 'OmniPath', 'weight')
for (col in cols) {
  col_attrs <- attrs[grep(col, attrs)]
  tmp_df <- c()
  for (attr in col_attrs) {
    tmp_df <- cbind(tmp_df, get.edge.attribute(unified_directed_network, attr))
  }
  tmp_df <- as.data.frame(tmp_df)
  tmp_df2 <- data.frame(id = seq(1:dim(tmp_df)[1]))
  tmp_df2$count <- rowSums(is.na(tmp_df))
  new_col <- c()
  for (i in 1:dim(tmp_df2)[1]) {
    if (tmp_df2[i, 'count'] == 2) {
      s <- NA
    } else if (col == 'weight') {
      s <- unique(na.exclude(as.numeric(as.vector(tmp_df[i,]))))
    } else {
      s <- unique(na.exclude(as.character(as.vector(tmp_df[i,]))))
    }
    new_col <- c(new_col, s)
  }
  tmp_df2[,col] <- new_col
  for (attr in col_attrs) {
    unified_directed_network <- delete_edge_attr(unified_directed_network, attr)
  }
  unified_directed_network <- set_edge_attr(unified_directed_network, name=col,
                                            value=tmp_df2[,col])
}
print('length of the unified network:')
print(length(V(unified_directed_network)))
# 3701

# C. Modification of vertex attributes.
attrs <- vertex_attr_names(unified_directed_network)
for (attr in attrs) {
  values <- get.vertex.attribute(unified_directed_network, name=attr)
  values[is.na(values)] <- FALSE
  unified_directed_network <- set_vertex_attr(unified_directed_network, 
                                              name=attr, value=values)
}

# D. Creation of the undirected version.
unified_undirected_network <- as.undirected(unified_directed_network,
                                            mode="collapse",
                                            edge.attr.comb = "first")
print('length of the unified network:')
print(length(V(unified_undirected_network)))
sum(duplicated(igraph::as_data_frame(unified_undirected_network))) # should be 0

# E. Saving of the unified networks
MASLD_unified_directed_network <- unified_directed_network
MASLD_unified_undirected_network <- unified_undirected_network

filename <- paste(main_dir, "MASLD_unified_directed_network_with_semantics.RData", sep='')
save(MASLD_unified_directed_network, file=filename)
filename <- paste(main_dir, "MASLD_unified_undirected_network_with_semantics.RData", sep='')
save(MASLD_unified_undirected_network, file=filename)

# F. Saving of the nodes of the unified networks in a simple data frame
df <- data.frame(gene_symbol=as_ids(V(MASLD_unified_undirected_network)))
filename <- paste(main_dir, "MASLD_nodes_with_semantics.csv", sep='')
write.table(df, file=filename, quote = FALSE, sep=',')


