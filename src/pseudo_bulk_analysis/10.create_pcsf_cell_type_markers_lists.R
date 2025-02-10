suppressMessages(library(igraph)) # 1.4.3
suppressMessages(library(PCSF)) # 0.99.1

load('../../results/networks/NAFLD_unified_undirected_network_with_semantics.RData')

base_network_df <- igraph::as_data_frame(NAFLD_unified_undirected_network)
base_network_df <- base_network_df[,c('from','to', 'weight')] 
base_network_df[,'cost'] <- 0.01 * 1/base_network_df[,'weight'] 
base_network_df <- base_network_df[,c('from','to', 'cost')] 
base_ppi <- construct_interactome(base_network_df) # it is undirected

files <- list.files('../../data/sc_liver_cell_atlas/filtered_cell_type_markers/')

output_dir <- '../../data/sc_liver_cell_atlas/pscf_cell_type_markers_markers/'
dir.create(output_dir, showWarnings = FALSE)

for (f in files) {
  filename <- paste('../../data/sc_liver_cell_atlas/filtered_cell_type_markers/', f, sep='')
  df <- read.csv(filename, row.names = 1, sep='\t')
  seeds_scores <- df$avg_log2FC
  names(seeds_scores) <- df$gene
  network <- PCSF_rand(base_ppi, terminals=seeds_scores, mu=5e-04, n=10, r=0.1)
  selected_genes <- as_ids(V(network))
  to_remove <- c(selected_genes[grep('^RPS', selected_genes)],
                 selected_genes[grep('^RPL', selected_genes)],
                 selected_genes[grep('^RPN', selected_genes)],
                 selected_genes[grep('^MT', selected_genes)],
                 selected_genes[grep('^MRPS', selected_genes)],
                 selected_genes[grep('^MRPL', selected_genes)])
  selected_genes <- selected_genes[!selected_genes %in% to_remove]
  output_df <- data.frame(gene_symbol=selected_genes)
  write.table(output_df, file = paste(output_dir, f, sep=''), sep='\t')
}





