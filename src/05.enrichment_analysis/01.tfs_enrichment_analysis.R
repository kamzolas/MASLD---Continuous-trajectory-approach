suppressMessages(library(dplyr)) # 1.1.4
suppressMessages(library(data.table)) # 1.17.8 rbindlist
suppressMessages(library(rlang)) # 1.1.6 duplicate
source("../library.R")

################################################################################
# Description:
################################################################################
# This script applies TF enrichment analysis in order to associate the gene 
# sets of co-expression modules with TFs. Only the modules which have been linked
# with MASLD variables through linear modelling, are used. The analysis is 
# performed using CollecTRI and other TF-regulon databases from EnrichR 
# (ENCODE and ChEA Consensus TFs from ChIP-X, DNA binding preferences from JASPAR, 
# TF protein-protein interactions).
# - Outputs:
# - A tsv file which contains all the statistically significant results for each 
# gene set module.
################################################################################


################################################################################
# Inputs
################################################################################
main_dir = '../../results/ucam_sanyal/wgcna_and_linear_modelling/grid_params/4_180/'
output_dir <- '../../results/ucam_sanyal/enrichment_analysis/'
load('../../data/annotation_databases/annotation_lists.RData')
databases <- c(
  "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
  "Transcription_Factor_PPIs",
  "TRANSFAC_and_JASPAR_PWMs",
  "COLLECTRI"
)


################################################################################
# 1. Loading of linear regression results for gene co-expression modules.
# Modules with FDR < 0.05, at least in one model, will be kept as significant.
#
# Additionally their filtered gene content is selected, in order to define the
# background gene set for the enrichment analysis.
################################################################################
coefs_df <- read.table(paste(main_dir, 'module_variable_coefficients_FINAL.tsv', sep=''), 
                       sep='\t', header=TRUE)
coefs_df <- coefs_df[coefs_df$p.value_adjust<0.05,] # STATISTICAL THRESHOLD 0.05 FOR THE MODULES
coefs_df <- coefs_df[!coefs_df$term %in% c('age', 'inflammation'),] 
modules <- unique(coefs_df$term)

modules_df <- read.table(paste(main_dir, 'filtered_modules.tsv',sep=''), sep='\t', header=TRUE)
modules_df <- modules_df[modules_df$module_color %in% modules,]
modules_gene_sets <- split(x=modules_df$gene_symbol, f=modules_df$module_color)
modules_background <- unique(unname(unlist(modules_gene_sets)))


################################################################################
# 2. Enrichment analysis using the fisher exact function (filtering with FDR < 0.05).
################################################################################
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
# 3. Modification and saving of enrichment analysis results.
################################################################################
output_df <- duplicate(tfs_enrichment)
output_df <- output_df[!grepl('mouse', output_df$tf),]
output_df$tf <- unlist(lapply(strsplit(output_df$tf, " "), head, n = 1L))
rownames(output_df) <- NULL
colnames <- c('tf', 'module', 'database', 'p.value', 'adj_p.value', 'overlap', "odds_ratio")
output_df <- output_df[, colnames]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.table(output_df, paste(output_dir, 'tfs_enrichment.tsv', sep=''), sep='\t', 
            quote = FALSE, row.names = FALSE)
