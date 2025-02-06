suppressMessages(library(dplyr)) # 1.1.2
suppressMessages(library(data.table)) # 1.14.18 rbindlist
suppressMessages(library(rlang)) # 1.1.3 duplicate
source("../library.R")


################################################################################
#
################################################################################



################################################################################
#
# 1. Input arguments
#
################################################################################
#args = commandArgs(trailingOnly=TRUE)
deep_split = 4
min_size = 180
main_dir = paste('../../results/wgcna_and_linear_modelling/grid_params/', 
                 as.character(deep_split), '_', as.character(min_size), '/', sep='')



################################################################################
#
# 2. Loading of linear regression results for gene co-expression modules.
# Modules with FDR < 0.05, at least in one model, will be kept as significant.
#
# Additionally their filtered gene content is retrieved, in order to define the
# background gene set for the enrichment analysis.
#
################################################################################
coefs_df <- read.table(paste(main_dir, 'module_variable_coefficients_FINAL.tsv', sep=''), 
                       sep='\t', header=TRUE)
coefs_df <- coefs_df[coefs_df$p.value_adjust<0.05,] # STATISTICAL THRESHOLD 0.05 FOR THE MODULES
coefs_df <- coefs_df[!coefs_df$term %in% c('age', 'inflammation'),] 
modules <- unique(coefs_df$term)

modules_df <- read.table(paste(main_dir, 'filtered_modules.tsv',sep=''), sep='\t', header=TRUE)
modules_gene_sets <- split(x=modules_df$gene_symbol, f=modules_df$module_color)
modules_background <- unique(unname(unlist(modules_gene_sets)))



################################################################################
#
# 3. Loading of the TF-regulon databases and execution of enrichment analysis
# using the fisher exact function (filtering with FDR < 0.05).
#
################################################################################
load('../../data/annotation_databases/annotation_lists.RData')
databases <- c(
  "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
  "Transcription_Factor_PPIs",
  "TRANSFAC_and_JASPAR_PWMs",
  "COLLECTRI"
)

tfs_enrichment <- data.frame()
for (database in databases) {
    annotation_list <- annotation_lists[[database]]
    for (module in modules) {
        print(paste(database, module))
        gene_set = modules_gene_sets[[module]]
        res = execute_fisher_exact_test(gene_set, annotation_list, TRUE, 
                                        modules_background)
        res$module = rep(module, dim(res)[1])
        res$database = rep(database, dim(res)[1])
        res$tf <- rownames(res)
        rownames(res) <- NULL
        # STATISTICAL THRESHOLD 0.05 TO DEFINE THE ENRICHED TFS
        res = res[res$adj_p.value < 0.05,] 
        if (dim(res)[1] > 0) {
            tfs_enrichment <- rbind(tfs_enrichment, res)
        } else {
          #pass without recording
        }
    }
}



################################################################################
#
# 4. Modification of the enrichment analysis results and saving them as a data 
# frame.
#
################################################################################
output_df <- duplicate(tfs_enrichment)
output_df <- output_df[!grepl('mouse', output_df$tf),]
output_df$tf <- unlist(lapply(strsplit(output_df$tf, " "), head, n = 1L))
rownames(output_df) <- NULL
colnames <- c('tf', 'module', 'database', 'p.value', 'adj_p.value', 'overlap', "odds_ratio")
output_df <- output_df[, colnames]

output_dir <- '../../results/networks/'
write.table(output_df, paste(output_dir, 'tfs_enrichment.tsv', sep=''), sep='\t', 
            quote = FALSE, row.names = FALSE)
