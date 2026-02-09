library(DESeq2) # 1.46.0
library(jsonlite) # 2.0.0
library(igraph) # 2.2.1
library(assertr) # 3.0.1

################################################################################
# Description
################################################################################
# Application of DESeq2 for the sliding window-based groups of patients. Each SW
# is compared with the previous one.
# Outputs:
# - A csv file which contains the results for all the comparisons along the SW-based
# trajectory
################################################################################


################################################################################
# Inputs
################################################################################
main_dir <- "../../data/ucam_sanyal/"
sorted_samples_filename <- paste(main_dir, "PC1_sorted_samples.csv", sep='')
template_filename <- paste(main_dir, 'metadata.csv', sep='')
counts_filename <- paste(main_dir, "counts_matrix.csv", sep='')
sw_filename <- paste(main_dir, 'sw_samples.csv', sep='')
deseq_results_filename <- paste(main_dir, 'deseq_results.csv', sep='')
ensembl_mapping <- read.table('../../data/ensembl_mapping.tsv', sep='\t',
                              header=TRUE, row.names = 1)

# Samples ranking
sorted_samples_df = read.table(sorted_samples_filename, sep=',', header=TRUE)
rownames(sorted_samples_df) = sorted_samples_df[,1] 
sorted_samples_df[,1] = NULL
sorted_samples_df[,2] = rownames(sorted_samples_df)
colnames(sorted_samples_df) <- c('Sorting_axis', 'Sample_name')

# Samples metadata & counts matrix
sample_names_col <- 'Sample.name'
batches <- c("Dataset", "Sex")
template_df = read.csv(file = template_filename, sep=',')
rownames(template_df) = template_df[,sample_names_col]
template_df <- template_df[rownames(sorted_samples_df),]
counts_matrix_df <- read.csv(file = counts_filename, sep = ",", check.names=FALSE,
                             row.names=1)
counts_matrix_df <- counts_matrix_df[, rownames(sorted_samples_df)]
counts_matrix_df <- counts_matrix_df[(rowSums(counts_matrix_df)>dim(counts_matrix_df)[2]),] #Exclude low expressed counts
counts_matrix = data.matrix(counts_matrix_df)
colnames(counts_matrix) == row.names(template_df)

# SW stratification
sw_df <- read.csv(sw_filename, sep=',')
sw_list <- list()
for (sw in rownames(sw_df)) {
  samples <- strsplit(sw_df[sw, 'samples'], ';')[[1]]
  sw_list[[as.character(sw)]] <- samples
}


################################################################################
# DESeq2 analysis
################################################################################
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
  sw1_counts_matrix <- counts_matrix[,sw1_samples,drop=FALSE]
  sw2_counts_matrix <- counts_matrix[,sw2_samples,drop=FALSE]
  
  # Customize the template and counts_matrix for sw1_samples
  sw1_template_df <- template_df[sw1_samples, batches, drop=FALSE]
  sw1_key <- 'SW1'
  sw1_template_df$sample_id <- paste(sw1_key, seq(1,length(sw1_samples)), sep='_')
  sw1_template_df$condition <- rep(sw1_key, length(sw1_samples))
  colnames(sw1_counts_matrix) = sw1_template_df$sample_id
  
  # Customize the template and counts_matrix for sw2_samples
  sw2_template_df <- template_df[sw2_samples, batches, drop=FALSE]
  sw2_key <- 'SW2'
  sw2_template_df$sample_id <- paste(sw2_key, seq(1,length(sw2_samples)), sep='_')
  sw2_template_df$condition <- rep(sw2_key, length(sw2_samples))
  colnames(sw2_counts_matrix) = sw2_template_df$sample_id
  
  # Concatenate the count matrices and templates
  sw_counts_matrix = cbind(sw1_counts_matrix, sw2_counts_matrix)
  coldata = rbind(sw1_template_df, sw2_template_df)
  coldata$sample_id == colnames(sw_counts_matrix)
  rownames(coldata) <- NULL
  coldata$condition <- as.factor(coldata$condition)
  
  # Factorization of variables
  coldata$batch_id <- as.character(col_concat(coldata[,batches]))
  coldata$batch_id <- as.factor(coldata$batch_id)
  
  # Run DESeq and get the results
  dds <- DESeqDataSetFromMatrix(countData = sw_counts_matrix,
                                colData = coldata,
                                design = ~ batch_id + condition)
  dds$condition <- relevel(dds$condition, ref = sw1_key)
  dds <- DESeq(dds, parallel=TRUE)
  res <- results(dds)
  res <- as.data.frame(res)
  res <- res[!is.na(res$padj),]
  res <- res[,c('log2FoldChange', 'pvalue', 'padj')]
  key <- paste('SW', as.character(sw-1), '_vs_', 'SW', as.character(sw), sep='')
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
final_df <- merge(final_df, ensembl_mapping, by.x=0, by.y='ensembl_gene_id',
                  all.x=TRUE, all.y=FALSE)
colnames(final_df)[1] <- 'ensembl_gene_id'
write.table(final_df, file = deseq_results_filename, sep=',')
