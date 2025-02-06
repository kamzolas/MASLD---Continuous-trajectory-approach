suppressMessages(library(dplyr)) # 1.1.2
source("../library.R")


################################################################################
# 
################################################################################



################################################################################
#
# 1. input arguments
#
################################################################################
dir = '../../results/networks/'



################################################################################
#
# 2. Loading of the TF activities results and selection of those TFs with at
# least one FDR value lower than 0.05 (at least in one sw).
#
################################################################################
filename <- paste(dir, 'de_and_msviper/msviper_results.tsv', sep='')
diff_tf_df <- read.table(filename, sep='\t', header=TRUE, row.names = 1)
cols <- colnames(diff_tf_df)[grepl(pattern = 'fdr', colnames(diff_tf_df))]
fdr_tf_df <- diff_tf_df[,cols]
scores <- rowSums(fdr_tf_df < 0.05)
diff_activated_tfs <- names(scores[scores > 0])
print(paste('Diff activated TFs:', length(diff_activated_tfs)))


################################################################################
#
# 3. Reactome enrichment analysis (FDR < 0.05).
#
################################################################################
load('../../data/annotation_databases/annotation_lists.RData')
annotation_list <- annotation_lists[["Reactome_2024"]]

res = execute_fisher_exact_test(diff_activated_tfs, annotation_list, FALSE)
res = res[res$adj_p.value < 0.05,]
res$Reactome <- row.names(res)
row.names(res) <- NULL
pathways_enrichment_df <- res



################################################################################
#
# 4. Modification of the results and saving as data frame.
#
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
filename <- paste(dir, 'reactome_enrichment_for_diff_activated_tfs.tsv', sep='')
write.table(x = pathways_enrichment_df, file = filename, sep='\t', 
            quote = FALSE, row.names=FALSE)
