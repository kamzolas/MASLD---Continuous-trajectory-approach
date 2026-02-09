library(readr) # 2.1.6 (read_csv)
library(dplyr) # 1.1.4

################################################################################ 
# Description
################################################################################
# Filtering of WGCNA module gene sets based on the correlation of gene expression 
# profiles with those of eigen-genes. If the module contains more than 1000 genes
# then the 20th quantile of correlation values is used as a cut-off. If the module
# is smaller, then genes with positive correlation values are selected.
# Outputs:
# - A tsv file containing the filtered module gene sets.
################################################################################


################################################################################
# Inputs
################################################################################
args = commandArgs(trailingOnly=TRUE)
deep_split = args[1] # for example 2
min_size = args[2] # for example 30
key <- paste(as.character(deep_split), as.character(min_size), sep='_')
main_dir <- paste('../../results/ucam_sanyal/wgcna_and_linear_modelling/grid_params/', 
                  key, '/', sep='')
counts_filename <- '../../data/ucam_sanyal/batch_corrected_counts_matrix.csv'
ensembl_mapping <- read.table('../../data/ensembl_mapping.tsv', sep='\t',
                              row.names=1, header = TRUE)
to_remove_genes <- ensembl_mapping[ensembl_mapping$external_gene_name == "",'ensembl_gene_id']


################################################################################
# 1. Loading and transformation of the expression dataset.
################################################################################
GeneXData <- read_csv(counts_filename)
colnames(GeneXData)[1] <- "GeneID"
gene.IDs <- GeneXData$GeneID
GeneXData <- as.data.frame(GeneXData)
row.names(GeneXData) <- gene.IDs
GeneXData$GeneName <- NULL
GeneXData$GeneID <- NULL
GeneXData <- as.data.frame(t(GeneXData))
colnames(GeneXData) <- gene.IDs


################################################################################
# 2. Removal of to_remove_genes set from the expression dataset.
################################################################################
cols <- setdiff(colnames(GeneXData), to_remove_genes)
GeneXData <- GeneXData[, cols]


################################################################################
# 3. Loading of the eigen-genes of gene co-expression nodules and their 
# respective gene sets
################################################################################
load(paste(main_dir, 'eigengenes.RData', sep=''))
modules_df <- read.table(paste(main_dir, 'modules.tsv',sep=''), sep='\t', header=TRUE)
modules_gene_sets <- split(x=modules_df$ensembl_gene_id, f=modules_df$module_color)
#lengths(modules_gene_sets)
#length(unlist(modules_gene_sets)) #17020


################################################################################
# 4. Filtering of module gene sets.
################################################################################
selected_genes <- c()
for (module in names(modules_gene_sets)) {
  gene_set <- modules_gene_sets[[module]]
  gene_set_expression <- GeneXData[,gene_set]
  eigen_gene <- MEList$eigengenes[,module,drop=FALSE]
  samples_order <- rownames(eigen_gene)
  gene_set_expression <- gene_set_expression[samples_order,]
  correlations <- t(as.data.frame(cor(eigen_gene, gene_set_expression)))
  correlations <- data.frame(correlations)
  colnames(correlations) <- c('cor_value')
  correlations <- correlations[order(correlations$cor_value, decreasing=TRUE),,drop=FALSE]
  #most_correlated_genes <- rownames(correlations[correlations$cor_value >= 0.4,,drop=FALSE])
  if (dim(correlations)[1] < 1000) {
    most_correlated_genes <- rownames(correlations[correlations$cor_value >= 0,,drop=FALSE])
  } else {
    most_correlated_genes <- rownames(correlations[correlations$cor_value >= quantile(correlations$cor_value, 0.2),,drop=FALSE])
  }
  print(paste(module, length(gene_set), length(most_correlated_genes)))
  selected_genes <- c(selected_genes, most_correlated_genes)
}


################################################################################
# 5. Saving of the new gene sets in a tab-delimited file.
################################################################################
final_df <- modules_df[modules_df$ensembl_gene_id %in% selected_genes,]
rownames(final_df) <- NULL
filename <- paste(main_dir, 'filtered_modules.tsv', sep='')
write.table(x=final_df, file = filename, sep='\t')
lengths(split(final_df$gene_symbol, final_df$module_color))
