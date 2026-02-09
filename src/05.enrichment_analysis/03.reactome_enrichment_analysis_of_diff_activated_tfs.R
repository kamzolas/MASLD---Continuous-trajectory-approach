suppressMessages(library(dplyr)) # 1.1.2
source("../library.R")

################################################################################
# Description:
################################################################################
# This script performs the same analysis as the previous one, but here the set of
# differentially activated TFs is used as input instead of the set of TFs which
# has been linked to a specific MASLD variable. This group of TFs is defined by the
# results of differential analysis along the SW-based disease trajectory.
# - Outputs:
# - A tsv file which contains all the statistically significant results
################################################################################


################################################################################
# Inputs
################################################################################
input_dir = '../../results/ucam_sanyal/de_and_tf_analysis/'
output_dir = '../../results/ucam_sanyal/enrichment_analysis/'
load('../../data/annotation_databases/annotation_lists.RData')
annotation_list <- annotation_lists[["Reactome_2024"]]


################################################################################
# 1. Loading of the TF activities results and selection of those TFs with at
# least one FDR value lower than 0.05 (at least in one sw).
################################################################################
filename <- paste(input_dir, 'tf_results.tsv', sep='')
diff_tf_df <- read.table(filename, sep='\t', header=TRUE, row.names = 1)
cols <- colnames(diff_tf_df)[grepl(pattern = 'fdr', colnames(diff_tf_df))]
fdr_tf_df <- diff_tf_df[,cols]
scores <- rowSums(fdr_tf_df < 0.05)
diff_activated_tfs <- names(scores[scores > 0])
print(paste('Diff activated TFs:', length(diff_activated_tfs)))


################################################################################
# 2. Reactome enrichment analysis (FDR < 0.05).
################################################################################
res = execute_fisher_exact_test(diff_activated_tfs, annotation_list, FALSE)
res = res[res$adj_p.value < 0.05,]
res$Reactome <- row.names(res)
row.names(res) <- NULL
pathways_enrichment_df <- res


################################################################################
# 3. Modification and saving of the results.
################################################################################
colnames(pathways_enrichment_df) <- c('p.value', 'odds_ratio', 'overlap', 
                                      'adj.p.value', 'genes', 'term')
pathways_enrichment_df$log10pvalue <- -log10(pathways_enrichment_df$adj.p.value)
pathways_enrichment_df$term_id <- sapply(pathways_enrichment_df$term, function(term){
  fields <- unlist(strsplit(term, ' '))
  fields[length(fields)]
})
pathways_enrichment_df <- pathways_enrichment_df[,c('term', 'term_id', 'odds_ratio', 
                                                    'overlap', 'p.value', 'adj.p.value',
                                                    'log10pvalue', 'genes')]
filename <- paste(output_dir, 'reactome_enrichment_for_diff_activated_tfs.tsv', sep='')
write.table(x = pathways_enrichment_df, file = filename, sep='\t', 
            quote = FALSE, row.names=FALSE)
