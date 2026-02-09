suppressMessages(library(WGCNA)) #1.73 (cutreeDynamic, moduleEigengenes, mergeCloseModules)

################################################################################ 
# Description
################################################################################
# This script needs to run for different parameters (deep_split and min_size) 
# of WGCNA. A specific folder will be generated for each parameter set to save
# the created modules and eigen-genes. In the next steps of the analysis, the
# derived module sets will be validated to find the optimal one, regarding their
# ability to model the MASLD variables.
#
# Examined parameters in the study:
# deep split values: 2, 3, 4
# min_size: 30, 60, 90, 120, 150, 180, 210, 240, 270, 300
# Outputs:
# - modules.tsv: tsv file for the final mapping of genes to modules 
# - eigengenes.RData: This object will include the expression profiles of the
# eigen-genes of modules
################################################################################

################################################################################
# Inputs
################################################################################
args = commandArgs(trailingOnly=TRUE)
deep_split = as.integer(args[1]) # for example: 2
min_size = as.integer(args[2]) # for example: 30
soft_power = 4
key <- paste(as.character(deep_split), as.character(min_size), sep='_')
main_dir <- '../../results/ucam_sanyal/wgcna_and_linear_modelling/'
output_dir = paste(main_dir, 'grid_params/', key, '/', sep='')
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

################################################################################
# 1. Set up the parallel calculation back-end
################################################################################
enableWGCNAThreads()

################################################################################
# 2. Load the objects from the previous step
################################################################################
load(paste(main_dir, "soft_threshold_results.RData", sep=''))
load(paste(main_dir, "TOM_distance_matrix.RData", sep=''))
load(paste(main_dir, "gene_tree.RData", sep=''))

################################################################################
# 3. Cut the hierarchical tree and define the modules
################################################################################
# deepSplit: values in [1-4].
# 1: Few clusters with big sizes 
# 4: Many clusters with small sizes
# minClusterSize: minimum size of modules
dynamicMods = cutreeDynamic(dendro = gene_tree, distM = TOM_distance_matrix, 
                            deepSplit = deep_split, pamRespectsDendro = FALSE, 
                            minClusterSize = min_size)
dynamicColors = labels2colors(dynamicMods)

################################################################################
# 4. Calculate EigenGenes and merge those which are very close
################################################################################
MEList= moduleEigengenes(GeneXData, colors = dynamicColors, softPower = soft_power)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)
thr = 0.1 # means correlation > 0.9
merge = mergeCloseModules(GeneXData, colors = dynamicColors, cutHeight = thr)
mergedColors = merge$colors
mergedMEs = merge$newMEs

################################################################################
# 5. Write modules membership in a tsv file
################################################################################
moduleColors = mergedColors
colorOrder = c("grey", standardColors(200))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
genes <- colnames(GeneXData)
gene_mod <- cbind(ensembl_gene_id=genes, module_color=moduleColors, 
                  module_label=moduleLabels)
modules_df <- merge(gene_mod, ensembl_mapping, by.x='ensembl_gene_id', 
                    by.y='ensembl_gene_id', no.dups=FALSE)
names(modules_df)[names(modules_df) == 'external_gene_name'] <- 'gene_symbol'
modules_df <- modules_df[,c('ensembl_gene_id', 'gene_symbol', 'module_label', 'module_color')]
modules_df$module_color <- paste("ME", modules_df$module_color, sep="")
modules_df <- modules_df[order(modules_df$module_label),]
write.table(x = modules_df, file = paste(output_dir, 'modules.tsv', sep=''), 
            sep='\t', quote = FALSE, row.names = FALSE)

################################################################################
# 6. Save geneTree and MElist
################################################################################
save(MEList, file=paste(output_dir, "eigengenes.RData", sep=''))
