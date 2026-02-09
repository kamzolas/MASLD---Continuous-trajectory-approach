suppressMessages(library(dplyr)) # 1.1.4
source("../library.R")

################################################################################
# Description:
################################################################################
# This script applies Reactome enrichment analysis in order to associate the 
# enriched TFs (from the previous step) with pathways. The analysis is performed
# for each MASLD variable. Firstly, the significant modules for each variable are
# selected from the generated linear models. Secondly, the TFs linked to these 
# modules are used as input gene set to run the enrichment analysis.
# - Outputs:
# - A tsv file which contains all the statistically significant results for 
# MASLD variable
################################################################################


################################################################################
# Inputs
################################################################################
main_dir = '../../results/ucam_sanyal/wgcna_and_linear_modelling/grid_params/4_180/'
output_dir <- '../../results/ucam_sanyal/enrichment_analysis/'
load('../../data/annotation_databases/annotation_lists.RData')
annotation_list <- annotation_lists[["Reactome_2024"]]


################################################################################
# 1. Loading of the TF enrichment analysis results (FDR<0.05).
################################################################################
tfs_enrichment_df <- read.table(paste(output_dir, 'tfs_enrichment.tsv', sep=''), sep='\t', 
                                header=TRUE)
selected_tfs <- tfs_enrichment_df[tfs_enrichment_df$adj_p.value < 0.05,]


################################################################################
# 2. Loading of the linear regression results for gene co-expression modules.
# Filtering with FDR<0.05.
################################################################################
coefs_df <- read.table(paste(main_dir, 'module_variable_coefficients_FINAL.tsv', sep=''), 
                       sep='\t', header=TRUE)
coefs_df <- coefs_df[coefs_df$p.value_adjust < 0.05,]
modules_per_variable <- sapply(unique(coefs_df$variable), function(variable) {
  terms <- coefs_df[coefs_df$variable == variable, 'term']
  terms <- setdiff(terms, c('inflammation', 'age'))
})


################################################################################
# 3. Reactome enrichment analysis (FDR < 0.05).
################################################################################
pathways_enrichment_df <- data.frame()
for (variable in unique(coefs_df$variable)){
  modules <- modules_per_variable[[variable]]
  tfs <- unique(selected_tfs[selected_tfs$module %in% modules, 'tf'])
  gene_set = tfs
  print(paste(variable, length(gene_set)))
  print(gene_set)
  next()
  res = execute_fisher_exact_test(gene_set, annotation_list, FALSE)
  res$variable = rep(variable, dim(res)[1])
  res = res[res$adj_p.value < 0.05,]
  res$Reactome <- row.names(res)
  row.names(res) <- NULL
  if (dim(res)[1] > 0) {
    pathways_enrichment_df <- rbind(pathways_enrichment_df, res)
  } else {
    #pass without recording
  }
}


################################################################################
# 4. Modification and saving of the results
################################################################################
colnames(pathways_enrichment_df) <- c('p.value', 'odds_ratio', 'overlap', 'adj.p.value', 
                                      'genes', 'variable', 'term')
pathways_enrichment_df$log10pvalue <- -log10(pathways_enrichment_df$adj.p.value)
pathways_enrichment_df$term_id <- sapply(pathways_enrichment_df$term, function(term){
  fields <- unlist(strsplit(term, ' '))
  fields[length(fields)]
})
pathways_enrichment_df <- pathways_enrichment_df[,c('term', 'term_id', 'variable',
                                                    'odds_ratio', 'overlap',
                                                    'p.value', 'adj.p.value',
                                                    'log10pvalue', 'genes')]
filename <- paste(output_dir, 'reactome_enrichment_for_enriched_tfs.tsv', sep='')
write.table(x = pathways_enrichment_df, file = filename, sep='\t', quote = FALSE, 
            row.names=FALSE)
