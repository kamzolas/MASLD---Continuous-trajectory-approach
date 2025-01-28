library(DESeq2) # 1.43.5
library(jsonlite) # 1.8.8
library(igraph) # 1.4.3


adjust_batch <- function(template) {
  b1 <- length(unique(template$Dataset))
  b2 <- length(unique(template$Sex))
  if (b1 == 1 && b2 == 1) {
    for (i in 1:dim(template)[1]) {
      template[i, 'Sex'] = as.character(i)
    }
  }
  return(template)
}


ref_pc_sorting_df = read.table("../../data/PC1_sorted_samples.csv", sep=',', 
                               header=TRUE)
rownames(ref_pc_sorting_df) = ref_pc_sorting_df[,1] 
ref_pc_sorting_df[,1] = NULL
ref_pc_sorting_df[,2] = rownames(ref_pc_sorting_df)
colnames(ref_pc_sorting_df) <- c('PC1_value', 'Sample.name')

# Samples metadata
template_df = read.csv(file = "../../data/metadata.csv")[, -c(1,3)]
rownames(template_df) = template_df$Sample.name

# Counts matrix
cts_ucamsanyal <- read.csv(file = "../../data/merged_counts.csv", sep = ",")
rownames(cts_ucamsanyal) = as.character(cts_ucamsanyal$X)
cts_ucamsanyal = cts_ucamsanyal[,-c(1)]
colnames(cts_ucamsanyal) = template_df$Sample.name
cts_ucamsanyal <- cts_ucamsanyal[(rowSums(cts_ucamsanyal)>dim(cts_ucamsanyal)[2]),]
unique(colnames(cts_ucamsanyal) == template_df$Sample.name)
cts = data.matrix(cts_ucamsanyal)
colnames(cts) == row.names(template_df)



# Get the optimal path and create the SWs
load('../../results/finding_of_optimal_SWs_sequence/final_paths_in_sw_graph.RData')
optimal_sw_graph <- 'sw_graphs_8_24_03_1267'
path <- paths_in_sw_graph[[optimal_sw_graph]]
path_nodes <- as_ids(V(path))
sw_list <- list()
sw = 1
for (node in path_nodes) {
  ids_str = strsplit(node, split='_')[[1]][2]
  ids = as.numeric(strsplit(ids_str, split='-')[[1]])
  samples <- rownames(ref_pc_sorting_df[ids,])
  sw_list[[as.character(sw)]] <- samples
  sw = sw + 1
}


# DESeq2
sws = as.integer(names(sw_list))
sws = sws[order(sws, decreasing = FALSE)]
deseq_results <- list()
for (sw in sws[2:length(sws)]) {

  sw1_samples <- sw_list[[as.character(sw-1)]]
  sw1_samples <- unlist(sw1_samples)
  sw2_samples <- sw_list[[as.character(sw)]]
  sw2_samples <- unlist(sw2_samples)
  
  sw1_samples <- setdiff(sw1_samples, sw2_samples)

  # Subset the counts matrix based on the input samples
  sw1_cts <- cts[,sw1_samples,drop=FALSE]
  sw2_cts <- cts[,sw2_samples,drop=FALSE]
  
  # Customize the template and cts for sw1_samples
  sw1_template_df <- template_df[sw1_samples, c('Dataset', 'Sex'), drop=FALSE]
  sw1_key <- 'SW1'
  sw1_template_df$sample_id <- paste(sw1_key, seq(1,length(sw1_samples)), sep='_')
  sw1_template_df$condition <- rep(sw1_key, length(sw1_samples))
  sw1_template_df <- adjust_batch(sw1_template_df)
  colnames(sw1_cts) = sw1_template_df$sample_id
  
  # Customize the template and cts for sw2_samples
  sw2_template_df <- template_df[sw2_samples, c('Dataset', 'Sex'), drop=FALSE]
  sw2_key <- 'SW2'
  sw2_template_df$sample_id <- paste(sw2_key, seq(1,length(sw2_samples)), sep='_')
  sw2_template_df$condition <- rep(sw2_key, length(sw2_samples))
  sw2_template_df <- adjust_batch(sw2_template_df)
  colnames(sw2_cts) = sw2_template_df$sample_id
  
  # Concatenate the count matrices and templates
  counts_matrix = cbind(sw1_cts, sw2_cts)
  coldata = rbind(sw1_template_df, sw2_template_df)
  rownames(coldata) <- NULL
  coldata$condition <- as.factor(coldata$condition)
  
  # Factorization of variables
  coldata$batch_id <- paste0(coldata$Dataset, coldata$Sex)
  coldata$batch_id <- as.factor(coldata$batch_id)
  
  # Run DESeq and get the results
  dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                colData = coldata,
                                design = ~ batch_id + condition)
  dds$condition <- relevel(dds$condition, ref = sw1_key)
  dds <- DESeq(dds, parallel=TRUE)
  res <- results(dds)
  res <- as.data.frame(res)
  res <- res[!is.na(res$padj),]
  res <- res[,c('log2FoldChange', 'pvalue', 'padj')]
  #degs <- dim(res[res$padj < 0.1,])[1]
  key <- paste('SW', as.character(sw-1), '_vs_', 'SW', as.character(sw), sep='')
  #deseq_results[[sw]] <- degs
  colnames(res) <- paste(key, colnames(res), sep='_')
  deseq_results[[key]] <- res
}



final_df <- data.frame(genes=rownames(deseq_results[[1]]))
rownames(final_df) <- final_df$genes
for (df in deseq_results) {
  final_df <- merge(final_df, df, by=0, all=TRUE)
  row.names(final_df) <- final_df$Row.names
  final_df$Row.names <- NULL
}
final_df$genes <- NULL

idx <- grep('padj',colnames(final_df))
final_df[,idx][is.na(final_df[,idx])] <- 1
idx <- grep('pvalue',colnames(final_df))
final_df[,idx][is.na(final_df[,idx])] <- 1
idx <- grep('log',colnames(final_df))
final_df[,idx][is.na(final_df[,idx])] <- 0

ensembl_mapping <- read.table('../../data/ensembl_mapping.tsv', sep='\t',
                              header=TRUE, row.names = 1)

final_df <- merge(final_df, ensembl_mapping, by.x=0, by.y='ensembl_gene_id',
                  all.x=TRUE, all.y=FALSE)
colnames(final_df)[1] <- 'ensembl_gene_id'
write.table(final_df, file = '../../data/deseq_results_SW_vs_previous.csv', sep=',')
