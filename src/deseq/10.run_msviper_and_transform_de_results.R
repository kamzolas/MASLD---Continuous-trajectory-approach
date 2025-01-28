suppressMessages(library(dplyr)) # 1.1.2
suppressMessages(library(tidyr)) # 1.3.0
suppressMessages(library(viper)) # 1.34.0
suppressMessages(library(purrr)) # 1.0.1
source("library.R")


################################################################################
#
# A function to calculate an average the p-value when there are 
# two or more results for a gene in the differential expression analysis
#
################################################################################
pvalue_transformation <- function(p) {
  10**mean(log10(p))
}


################################################################################
#
# A function to transform regulon data frame to regulon object (created by Giannis)
#
################################################################################
df2regulon <- function(df) {
  regulon = df %>%
    split(.$tf) %>%
    map(function(dat) {
      tf = dat %>% distinct(tf) %>% pull()
      targets = setNames(dat$mor, dat$target)
      likelihood = dat$likelihood
      list(tfmode =targets, likelihood = likelihood)
    })
  return(regulon)
}


################################################################################
#
# A worker function to calculate TF activities given a specific signature of gene
# scores
#
################################################################################
calculate_tf_activities <- function(signature) {
  mrs <- msviper(ges = signature, regulon = df2regulon(Regulon_file),
                ges.filter = F, minsize = 4)
  mrs_df <- data.frame('size' = mrs$es$size,
                       'p.value' = mrs$es$p.value,
                       'fdr' = p.adjust(mrs$es$p.value, method = 'fdr'),
                       'nes' = mrs$es$nes)
  #mrs_df$score <- -log10(mrs_df$fdr)*sign(mrs_df$nes)
  #mrs_df$score_plus_one <- (-log10(mrs_df$fdr)+1)*sign(mrs_df$nes)
  return(mrs_df)
}



################################################################################
#
# Output files:
#
# There are three output files, related to the DE and TF activity analyses. All 
# of them will be stores in one new folder:
# ../results/deep_split+'_'+min_size/networks/de_and_msviper
#
# 1. msviper_results.tsv: TF activation scores in each SW
# 2. de_results.tsv: DE scores in each SW
# 3. msviper_results.RData: All the results of TF activity analysis
#
################################################################################



################################################################################
#
# Input arguments
#
################################################################################
args = commandArgs(trailingOnly=TRUE)
deep_split = args[1] # 4
min_size = args[2] # 180
key = paste(deep_split, min_size, sep='_')
dir = paste('../results/', key, '/', sep='')



################################################################################
#
# 1. Loading of the modules gene sets. 
# 
# Only the significant modules (based on regression models) will be kept. 
# The respective filtered gene sets will be retrieved to create the 
# 'modules_gene_sets' list. Finally, a unified pool of these genes will be 
# created ('modules_genes') and only these genes will be used to calculate the 
# TF activities. 
#
################################################################################

# 1A. Selection of the significant modules (regression models)
filename <- paste(dir, 'module_variable_coefficients_FINAL.tsv', sep='')
modules_per_var_df <- read.table(filename, sep='\t', header=TRUE)
modules_per_var_df <- modules_per_var_df[modules_per_var_df$p.value_adjust < 0.05,]
modules_per_var_df <- modules_per_var_df[!modules_per_var_df$term %in% c('age', 'inflammation'),] 

# 1B. Retrieval of module gene sets (filtered with correlation score) and creation
# of the unified gene set
filename <- paste(dir, 'filtered_modules.tsv',sep='')
modules_df <- read.table(filename, sep='\t', header=TRUE)
modules_gene_sets <- split(x=modules_df$gene_symbol, f=modules_df$module_color)
for (module in names(modules_gene_sets)) {
  if (!(module %in% modules_per_var_df$term)) {
    modules_gene_sets[[module]] <- NULL
  }
}
modules_genes <- unique(unlist(modules_gene_sets))



################################################################################
#
# 2.Loading and transformation of the results from DE analysis. This data frame
# contains 3 columns for each SWs comparison (log2FC values, p-values and FDRs).
#
# 2A. A data frame which contains only the log2FC values will be created, 
# removing duplicated gene symbols by averaging the respective values.
#
# 2B. A data frame which contains only the FDR values will be created, 
# removing duplicated gene symbols by averaging the respective values.
#
################################################################################
de_analysis_df <- read.csv('../data/deseq_results_SW_vs_previous.csv')
columns <- colnames(de_analysis_df)
#columns <- c('ensembl_gene_id', columns)
#columns <- columns[1: length(columns)-1]
#colnames(de_analysis_df) <- columns
row.names(de_analysis_df) <- de_analysis_df$ensembl_gene_id
de_analysis_df$ensembl_gene_id <- NULL
de_analysis_df <- de_analysis_df[!de_analysis_df$external_gene_name == "",]
de_analysis_df <- de_analysis_df[!is.na(de_analysis_df$external_gene_name),]
dim(de_analysis_df)

# 2A. log2FC data frame
columns <- colnames(de_analysis_df)[grepl(x=colnames(de_analysis_df), pattern='log2')]
logfc_df <- de_analysis_df[,c(columns,'external_gene_name')]
logfc_df[is.na(logfc_df)] <- 0 
logfc_df <- logfc_df[!logfc_df$external_gene_name == "",]
dim(logfc_df)
logfc_df <- logfc_df %>% 
  group_by(external_gene_name) %>% 
  select(-external_gene_name) %>% 
  summarise_if(is.numeric, mean)
logfc_df <- as.data.frame(logfc_df)
rownames(logfc_df) <- logfc_df$external_gene_name
logfc_df$external_gene_name <- NULL
colnames(logfc_df) <- paste('SW', seq(2,dim(logfc_df)[2]+1), sep='')
dim(logfc_df)

# 2B. FDR data frame
columns <- colnames(de_analysis_df)[grepl(x=colnames(de_analysis_df), pattern='padj')]
adj_pvalues_df <- de_analysis_df[,c(columns,'external_gene_name')]
adj_pvalues_df[is.na(adj_pvalues_df)] <- 1
adj_pvalues_df <- adj_pvalues_df[!adj_pvalues_df$external_gene_name == "",]
dim(adj_pvalues_df)
adj_pvalues_df <- adj_pvalues_df %>% 
  group_by(external_gene_name) %>% 
  select(-external_gene_name) %>% 
  summarise_if(is.numeric, pvalue_transformation)
adj_pvalues_df <- as.data.frame(adj_pvalues_df)
rownames(adj_pvalues_df) <- adj_pvalues_df$external_gene_name
adj_pvalues_df$external_gene_name <- NULL
colnames(adj_pvalues_df) <- paste('SW', seq(2,dim(adj_pvalues_df)[2]+1), sep='')
dim(adj_pvalues_df)



################################################################################
#
# 3. Loading of CollectTRI data and creation of the 'background_genes' variable
# which will be used as reference to calculate the TF activites.
#
# The above data frame (log2FC and FDR) will be filtered to keep the results only
# for genes in the 'background_genes' variable.
#
################################################################################
Regulon_file<- read.csv("../data/collectTRI_network.tsv", sep='\t', header=T)
#Regulon_file<- read.csv("../data/obsolete/human_network_dorothea.csv", sep=',', header=T)
#Regulon_file<- Regulon_file[Regulon_file$confidence=='A'| Regulon_file$confidence=='B',]
background_genes <- intersect(modules_genes, rownames(adj_pvalues_df))
length(background_genes)
final_adj_pvalues_df <- adj_pvalues_df[background_genes, ]
final_logfc_df <- logfc_df[background_genes, ]



################################################################################
#
# 4. Perform the TF activity analysis for each sw and store all the results in 
# three different formats:
# 
# A. A list of dataframes, one for each sw
# B. A unified dataframe, which will be used later for the network analysis
################################################################################

# 4A
SWs = paste('SW', seq(2,dim(adj_pvalues_df)[2]+1,1), sep='')
msviper_results_list <- list()
for (sw in SWs) {
    # Step 1: construct the gene signature for msviper
    sw_df <- merge(final_adj_pvalues_df[,sw,drop=FALSE], 
                   final_logfc_df[,sw,drop=FALSE], 
                   by.x=0, by.y=0)
    colnames(sw_df) <- c('gene_symbol', 'padj', 'logfc')
    rownames(sw_df) <- sw_df$gene_symbol
    sw_df$gene_symbol <- NULL
    signature <- sw_df[,'padj',drop=TRUE]
    signature <- qnorm(signature/2, lower.tail = FALSE)*sign(sw_df[,'logfc',drop=TRUE])
    names(signature) <- rownames(sw_df)
    # Step 2: Run msviper and save the results
    msviper_results_df <- calculate_tf_activities(signature)
    msviper_results_list[[sw]] <- msviper_results_df
}

# 4B
msviper_results_df <- data.frame(tf=unique(Regulon_file$tf))
rownames(msviper_results_df) <- msviper_results_df$tf
for (sw in names(msviper_results_list)) {
    tmp <-  msviper_results_list[[sw]][,c('nes', 'fdr'),drop=FALSE]
    colnames(tmp) <- paste(colnames(tmp), c(sw), sep='_')
    msviper_results_df <- merge(msviper_results_df, tmp, 
                                by.x=0, by.y=0, all=TRUE)
    rownames(msviper_results_df) <- msviper_results_df$Row.names
    msviper_results_df$Row.names <- NULL
    
}
msviper_results_df[is.na(msviper_results_df)] <- 1
msviper_results_df$tf <- NULL




################################################################################
#
# 5. Transformation of the DE analysis results.
#
################################################################################
colnames(final_logfc_df) <- paste('logfc', colnames(final_logfc_df), sep='_')
colnames(final_adj_pvalues_df) <- paste('fdr', colnames(final_adj_pvalues_df), sep='_') 
de_df = merge(final_logfc_df, final_adj_pvalues_df, by.x=0, by.y=0, all = TRUE)

rownames(de_df) <- de_df$Row.names
de_df$Row.names <- NULL
sorted_columns <- c()
for (sw in SWs){
  sorted_columns <- c(sorted_columns, c(paste('logfc_',sw,sep=''), paste('fdr_',sw,sep='')))
}
de_df <- de_df[, sorted_columns]

################################################################################
#
# 6. Saving of the results (both TF activities & DEGs).
#
################################################################################
output_dir <- paste(dir, 'networks/de_and_msviper', sep='')
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

filename <- paste(output_dir, '/', 'msviper_results.tsv', sep='')
write.table(msviper_results_df, file = filename, quote = FALSE, sep='\t')

filename <- paste(output_dir, '/', 'de_results.tsv', sep='')
write.table(de_df, file = filename, quote = FALSE, sep='\t')

filename <- paste(output_dir, '/', 'msviper_results.RData', sep='')
save(msviper_results_list, file = filename)



