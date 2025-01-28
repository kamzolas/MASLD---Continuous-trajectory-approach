library(readr) # 2.1.4 (read_csv)
library(dplyr) # 1.1.4


################################################################################
#
# 1. input arguments
#
################################################################################
args = commandArgs(trailingOnly=TRUE)
deep_split = args[1] #4
min_size = args[2] #30

key <- paste(as.character(deep_split), as.character(min_size), sep='_')
output_dir <- paste('../../results/wgcna_and_linear_modelling/grid_params/', key, '/', sep='')

################################################################################
#
# 2. Loading of the Ensembl mapping to create the 'to_remove_genes' vector.
# Ensembl ids without corresponding external_gene_name will be filtered out
# in the following steps.
#
################################################################################
#remove_genes <- FALSE
ensembl_mapping <- read.table('../../data/ensembl_mapping.tsv', sep='\t',
                              row.names=1, header = TRUE)
to_remove_genes <- ensembl_mapping[ensembl_mapping$external_gene_name == "",'ensembl_gene_id']




################################################################################
#
# 3. Loading and transformation of the expression dataset.
#
################################################################################
GeneXData <- read_csv("../../data/batch_corrected_counts_(dataset+gender).csv")
colnames(GeneXData)[1] <- "GeneID"
gene.IDs <- GeneXData$GeneID
GeneXData <- as.data.frame(GeneXData)
row.names(GeneXData) <- gene.IDs
GeneXData$GeneName <- NULL
GeneXData$GeneID <- NULL
GeneXData <- as.data.frame(t(GeneXData))
colnames(GeneXData) <- gene.IDs


################################################################################
#
# 4. Removal of to_remove_genes set from the expression dataset.
#
################################################################################
cols <- setdiff(colnames(GeneXData), to_remove_genes)
GeneXData <- GeneXData[, cols]


################################################################################
#
# 5. Loading of the eigen-genes of gene co-expression nodules and their 
# respective gene sets
#
################################################################################
load(paste(output_dir, 'eigengenes.RData', sep=''))
modules_df <- read.table(paste(output_dir, 'modules.tsv',sep=''), sep='\t', header=TRUE)
modules_gene_sets <- split(x=modules_df$ensembl_gene_id, f=modules_df$module_color)



################################################################################
#
# 6. Selection of genes with correlation equal/greater than 0.4 with the 
# eigen-gene of their co-expression module.
#
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
  most_correlated_genes <- rownames(correlations[correlations$cor_value >= 0.4,,drop=FALSE])
  print(paste(module, length(gene_set), length(most_correlated_genes)))
  selected_genes <- c(selected_genes, most_correlated_genes)
}

# Module  WholeGeneSet  GeneSetCorr04
# [1] "MEblack 664 552"
# [1] "MEblue 3641 1077"
# [1] "MEbrown 2720 1614"
# [1] "MEgreen 1033 353"
# [1] "MEgreenyellow 233 111"
# [1] "MEmagenta 371 349"
# [1] "MEpink 510 284"
# [1] "MEpurple 273 151"
# [1] "MEred 871 522"
# [1] "MEsalmon 192 157"
# [1] "MEtan 197 182"
# [1] "MEturquoise 3686 1155"
# [1] "MEyellow 2629 845"


################################################################################
#
# 7. Saving of the new gene sets in a tab-delimited file.
#
################################################################################
final_df <- modules_df[modules_df$ensembl_gene_id %in% selected_genes,]
rownames(final_df) <- NULL
filename <- paste(output_dir, 'filtered_modules.tsv', sep='')
write.table(x=final_df, file = filename, sep='\t')
lengths(split(final_df$gene_symbol, final_df$module_color))

#MEblack        MEblue       MEbrown       MEgreen MEgreenyellow     MEmagenta        MEpink      MEpurple         MEred 
#552          1077          1614           353           111           349           284           151           522 
#MEsalmon         MEtan   MEturquoise      MEyellow 
#157           182          1155           845 
