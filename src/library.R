library(data.table)
library(dplyr)
library(igraph)

execute_fisher_exact_test <- function(input_gene_set, annotation_list, 
                                      adjust_annotation = FALSE,
                                      gene_set_for_adjustment = NULL) {
  
  # background definition
  background <- unique(unname(unlist(annotation_list)))
  #print(paste('length of', database, 'gene content:', length(annotation_background)))
  if (adjust_annotation == TRUE) {
    background <- intersect(background, gene_set_for_adjustment)
  } else {
    #pass
  }
  #print(paste('length of used background:', length(background)))
  
  # input_gene_set definition
  #print(paste('length of input gene set:', length(input_gene_set)))
  if (adjust_annotation == TRUE) {
    input_gene_set <- intersect(background, input_gene_set)
  } else {
    #pass
  }
  #print(paste('length of analyzed gene set:', length(input_gene_set)))
  
  # Apply Fisher Exact to the annotation list 
  fisher_exact_list <- lapply(annotation_list, function(term_gene_set) {
    if (adjust_annotation == TRUE) {
      term_gene_set <- intersect(background, term_gene_set)
    } else {
      # pass
    }
    a11 <- length(intersect(term_gene_set, input_gene_set)) # genes both term set and sample
    a12 <- length(setdiff(term_gene_set, input_gene_set)) # genes in term set but not in sample
    a21 <- length(intersect(setdiff(background, term_gene_set), input_gene_set)) # genes not in term set but in sample
    a22 <- length(setdiff(setdiff(background, term_gene_set), input_gene_set)) # all the other genes
    contigency_table <- matrix(c(a11, a21, a12, a22), nrow = 2)
    test <- fisher.test(x=contigency_table, alternative='greater', or = 1)
    overlap = paste(as.character(a11), '/', as.character(a11+a12), sep = '') # enrichment
    c(test$p.value, test$estimate, overlap, 
      paste(intersect(term_gene_set, input_gene_set), collapse=';'))
  })
  fisher_exact_df <- transpose(data.frame(fisher_exact_list))
  colnames(fisher_exact_df) <- c('p.value', 'odds_ratio', 'overlap', 'genes')
  fisher_exact_df$p.value <- as.numeric(fisher_exact_df$p.value)
  fisher_exact_df$odds_ratio <- as.numeric(fisher_exact_df$odds_ratio)
  fisher_exact_df$adj_p.value <- p.adjust(fisher_exact_df$p.value, method="BH")
  row.names(fisher_exact_df) <- names(fisher_exact_list)
  fisher_exact_df <- fisher_exact_df[,c('p.value', 'odds_ratio', 'overlap', 'adj_p.value', 'genes')]
  return(fisher_exact_df)
}



laplacian_normalization <- function(network) {
  edges_df <- igraph::as_data_frame(network)
  edges_df$edge_id <- as.integer(row.names(edges_df))
  if (is_directed(network)) {
    # Get the outer degrees for the sources
    degrees_df <- data.frame(degree = degree(network, mode = 'out')) # in out all
    edges_df <- merge(edges_df, degrees_df, by.x='from', by.y=0, all.y=FALSE, sort=FALSE)
    edges_df$degree[edges_df$degree == 0] = 1
    colnames(edges_df)[colnames(edges_df) == 'degree'] = 'from_degree'
    # Get the inner degrees for the sinks
    degrees_df <- data.frame(degree = degree(network, mode = 'in')) # in out all
    edges_df <- merge(edges_df, degrees_df, by.x='to', by.y=0, all.y=FALSE, sort=FALSE)
    edges_df$degree[edges_df$degree == 0] = 1
    colnames(edges_df)[colnames(edges_df) == 'degree'] = 'to_degree'
    # Calculate the new attributes
    edges_df <- edges_df[order(edges_df$edge_id),]
    edges_df$lapl_norm <- 1/(sqrt(edges_df$from_degree * edges_df$to_degree))
    edges_df$lapl_norm_weight <- edges_df$weight/(sqrt(edges_df$from_degree * edges_df$to_degree))
  } else {
    # Get the degrees (mode=all)
    degrees_df <- data.frame(degree = degree(network, mode = 'all')) # in out all
    edges_df <- merge(edges_df, degrees_df, by.x='from', by.y=0, all.y=FALSE, sort=FALSE)
    edges_df$degree[edges_df$degree == 0] = 1
    colnames(edges_df)[colnames(edges_df) == 'degree'] = 'from_degree'
    #
    edges_df <- merge(edges_df, degrees_df, by.x='to', by.y=0, all.y=FALSE, sort=FALSE)
    edges_df$degree[edges_df$degree == 0] = 1
    colnames(edges_df)[colnames(edges_df) == 'degree'] = 'to_degree'
    # Calculate the new attributes
    edges_df <- edges_df[order(edges_df$edge_id),]
    edges_df$lapl_norm <- 1/(sqrt(edges_df$from_degree * edges_df$to_degree))
    edges_df$lapl_norm_weight <- edges_df$weight/(sqrt(edges_df$from_degree * edges_df$to_degree))
  }
  output_network <- set_edge_attr(network, name='lapl_norm',
                                  index=edges_df$edge_id,
                                  value=edges_df$lapl_norm)
  output_network <- set_edge_attr(output_network, name='lapl_norm_weight',
                                  index=edges_df$edge_id,
                                  value=edges_df$lapl_norm_weight)
  return(output_network)
}








