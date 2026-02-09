library(igraph)
library(rlang)
library(jsonlite)

################################################################################
# Description
################################################################################
# Extraction of patients stratification for the optimal sw_graph path. The user
# should give as input the id of the optimal path in order to get the 
# stratification of patients into sliding windows.
# Output files:
# - A siple csv file which contains the mapping of patients to sliding windows.
################################################################################


################################################################################
# Inputs
################################################################################
main_dir <- '../../results/ucam_sanyal/finding_of_optimal_SWs_sequence/'
output_dir <- '../../data/ucam_sanyal/'
optimal_path_key <- 'sw_graph_8_24_03_2441'
key = optimal_path_key
output_filename = 'sw_samples.csv'

sorted_samples_filename <- "../../data/ucam_sanyal/PC1_sorted_samples.csv"
sorted_samples_df = read.table(sorted_samples_filename, sep=',', header=TRUE)
rownames(sorted_samples_df) = sorted_samples_df[,1] 
sorted_samples_df[,1] = NULL
sorted_samples_df[,2] = rownames(sorted_samples_df)
colnames(sorted_samples_df) <- c('Sorting_axis', 'Sample_name')
load(paste(main_dir, 'final_paths_in_sw_graph.RData', sep=''))


################################################################################
# Extract the path and create the sw_samples data frame
################################################################################
path <- paths_in_sw_graph[[optimal_path_key]]
path_nodes <- as_ids(V(path))#[2:(length(path)-1)] # remove start and end nodes
sw_data <- list()
sw = 1
for (node in path_nodes) {
  ids_str = strsplit(node, split='_')[[1]][2]
  ids = as.numeric(strsplit(ids_str, split='-')[[1]])
  samples <- rownames(sorted_samples_df[ids,])
  sw_data[[as.character(sw)]] <- samples
  sw = sw + 1
}

sw_df <- data.frame()
for (sw in names(sw_data)) {
  tmp_vector <- paste(sw_data[[sw]], collapse=';')
  sw <- paste('SW', sw, sep = '_')
  tmp_vector <- c(sw, tmp_vector)
  sw_df <- rbind(sw_df, tmp_vector)
}
colnames(sw_df) <- c('sw', 'samples')
write.table(sw_df, file=paste(output_dir, output_filename, sep=''), sep=',', 
            row.names = FALSE, quote = FALSE)

