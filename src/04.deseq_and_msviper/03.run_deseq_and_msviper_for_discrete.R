suppressMessages(library(DESeq2)) # 1.46.0
suppressMessages(library(dplyr)) # 1.1.2
suppressMessages(library(tidyr)) # 1.3.0
suppressMessages(library(viper)) # 1.34.0
suppressMessages(library(purrr)) # 1.0.1

################################################################################
# Description
################################################################################
# Application of DESeq2 and msVIPER for discrete groups of patients, based on 
# SAF and NAS scores (check metadata file).
# SAF score division: Healthy, Mild, Moderate and Severe
# NAS score division: values from 0 to 7.
# Outputs:
# - Tsv files which contains the results for all the comparisons along the 
# trajectory (SAF- and NAS-based analysis).
################################################################################


################################################################################
# A simple function to adjust the batches if they contain only one value
################################################################################
adjust_batch <- function(template) {
  b1 <- length(unique(template$Dataset))
  b2 <- length(unique(template$Sex))
  if (b1 == 1 && b2 == 1) {
    for (i in 1:dim(template)[1]) {
      template[i, 'Sex'] = as.character(i)
    }
  } else {
    #pass
  }
  return(template)
}


################################################################################
# This small function is used to concatenate the duplicated results in 
# differential expression analysis, by averaging the p-value
################################################################################
pvalue_transformation <- function(p) {
  10**mean(log10(p))
}


################################################################################
# A function to transform regulon data frame to regulon object
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
# A function to calculate differential TF-activities given a specific signature 
# of gene scores
################################################################################
calculate_tf_activities <- function(signature) {
  mrs <- msviper(ges = signature, regulon = df2regulon(Regulon_file),
                 ges.filter = F, minsize = 30)
  mrs_df <- data.frame('size' = mrs$es$size,
                       'p.value' = mrs$es$p.value,
                       'fdr' = p.adjust(mrs$es$p.value, method = 'fdr'),
                       'nes' = mrs$es$nes)
  return(mrs_df)
}


################################################################################
# This function runs deseq2 for a pre-defined sequence of conditions
################################################################################
deseq_worker <- function(sorted_conditions, samples_in_conditions, output_dir) {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # DESeq2
  deseq_results <- list()
  idx1 = 1
  for (idx1 in seq(length(sorted_conditions)-1)) {
    idx2 <- idx1 + 1
    samples1 <- samples_in_conditions[[sorted_conditions[idx1]]]
    samples2 <- samples_in_conditions[[sorted_conditions[idx2]]]
    
    # Subset the counts matrix based on the input samples
    counts_matrix1 <- counts_matrix[,samples1,drop=FALSE]
    counts_matrix2 <- counts_matrix[,samples2,drop=FALSE]
    
    # Customize the template and cts for sw1_samples
    template_df1 <- template_df[samples1, c('Dataset', 'Sex'), drop=FALSE]
    key1 <- 'SW1'
    template_df1$sample_id <- paste(key1, seq(1,length(samples1)), sep='_')
    template_df1$condition <- rep(key1, length(samples1))
    template_df1 <- adjust_batch(template_df1)
    colnames(counts_matrix1) = template_df1$sample_id
    
    # Customize the template and cts for sw2_samples
    template_df2 <- template_df[samples2, c('Dataset', 'Sex'), drop=FALSE]
    key2 <- 'SW2'
    template_df2$sample_id <- paste(key2, seq(1,length(samples2)), sep='_')
    template_df2$condition <- rep(key2, length(samples2))
    template_df2 <- adjust_batch(template_df2)
    colnames(counts_matrix2) = template_df2$sample_id
    
    # Concatenate the count matrices and templates
    combined_counts_matrix = cbind(counts_matrix1, counts_matrix2)
    coldata = rbind(template_df1, template_df2)
    rownames(coldata) <- NULL
    coldata$condition <- as.factor(coldata$condition)
    
    # Factorization of variables
    coldata$batch_id <- paste0(coldata$Dataset, coldata$Sex)
    coldata$batch_id <- as.factor(coldata$batch_id)
    
    # Run DESeq and get the results
    dds <- DESeqDataSetFromMatrix(countData = combined_counts_matrix,
                                  colData = coldata,
                                  design = ~ batch_id + condition)
    dds$condition <- relevel(dds$condition, ref = key1)
    dds <- DESeq(dds, parallel=FALSE)
    res <- results(dds, cooksCutoff = FALSE)
    res <- as.data.frame(res)
    res <- res[!is.na(res$padj),]
    res <- res[,c('log2FoldChange', 'pvalue', 'padj')]
    key <- paste(sorted_conditions[idx2], '_vs_', sorted_conditions[idx1], sep='')
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
  
  ensembl_mapping <- read.table('../../data/ensembl_mapping.tsv', sep='\t')
  ensembl_mapping[,1] <- NULL
  colnames(ensembl_mapping) <- c('ensembl_gene_id', 'external_gene_name')
  
  final_df <- merge(final_df, ensembl_mapping, by.x=0, by.y='ensembl_gene_id',
                    all.x=TRUE, all.y=FALSE)
  colnames(final_df)[1] <- 'ensembl_gene_id'
  write.table(final_df, file = paste(output_dir, '/de_results.csv', sep=''), sep=',')
}


################################################################################
# This function applies msviper using the results of differential expression analysis
################################################################################
msviper_worker <- function(adj_pvalues_df, logfc_df, output_dir) {
  
  background_genes <- intersect(modules_genes, rownames(adj_pvalues_df))
  length(background_genes)
  final_adj_pvalues_df <- adj_pvalues_df[background_genes, ]
  final_logfc_df <- logfc_df[background_genes, ]
  
  discrete_stages = colnames(final_logfc_df)
  msviper_results_list <- list()
  for (stage in discrete_stages) {
    stage_df <- merge(final_adj_pvalues_df[,stage, drop=FALSE], 
                      final_logfc_df[,stage, drop=FALSE],
                      by.x=0, by.y=0)
    colnames(stage_df) <- c('gene_symbol', 'padj', 'logfc')
    rownames(stage_df) <- stage_df$gene_symbol
    stage_df$gene_symbol <- NULL
    signature <- stage_df[,'padj',drop=TRUE]
    signature <- qnorm(signature/2, lower.tail = FALSE)*sign(stage_df[,'logfc',drop=TRUE])
    names(signature) <- rownames(stage_df)
    msviper_results_df <- calculate_tf_activities(signature)
    msviper_results_list[[stage]] <- msviper_results_df
  }
  
  msviper_results_df <- data.frame(tf=unique(Regulon_file$tf))
  rownames(msviper_results_df) <- msviper_results_df$tf
  for (stage in names(msviper_results_list)) {
    
    tmp <- msviper_results_list[[stage]][,c('nes', 'fdr'),drop=FALSE]
    colnames(tmp) <- paste(colnames(tmp), c(stage), sep='_')
    msviper_results_df <- merge(msviper_results_df, tmp, 
                                by.x=0, by.y=0, all=TRUE)
    rownames(msviper_results_df) <- msviper_results_df$Row.names
    msviper_results_df$Row.names <- NULL
  }
  
  msviper_results_df[is.na(msviper_results_df)] <- 1
  msviper_results_df$tf <- NULL
  filename <- paste(output_dir, 'tf_results.tsv', sep='')
  write.table(msviper_results_df, file = filename, quote = FALSE, sep='\t')
  
}


################################################################################
# Inputs
################################################################################
modules_dir <- '../../results/ucam_sanyal/wgcna_and_linear_modelling/grid_params/4_180/'
lm_filename <- paste(modules_dir, 'module_variable_coefficients_FINAL.tsv', sep='')
modules_filename <- paste(modules_dir, 'filtered_modules.tsv',sep='')
Regulon_file<- read.csv("../../data/collectTRI_network.tsv", sep='\t', header=T)
main_dir <- "../../data/ucam_sanyal/"
sorted_samples_filename <- paste(main_dir, "PC1_sorted_samples.csv", sep='')
template_filename <- paste(main_dir, 'metadata.csv', sep='')
counts_filename <- paste(main_dir, "counts_matrix.csv", sep='')
sample_names_col <- 'Sample.name'
batches <- c("Dataset", "Sex")


################################################################################
# 1. Load the modules, keep only the significant ones based on regression 
# analysis results and then retrieve their filtered gene sets. Finally, create
# a pool of genes by unifying the module gene sets. Only this pool of genes will 
# be used in the downstream analysis 
################################################################################

# 1A. Get the significant modules
filename <- paste(modules_dir, 'module_variable_coefficients_FINAL.tsv', sep='')
modules_per_var_df <- read.table(filename, sep='\t', header=TRUE)
modules_per_var_df <- modules_per_var_df[modules_per_var_df$p.value_adjust < 0.05,]
modules_per_var_df <- modules_per_var_df[!modules_per_var_df$term %in% c('age', 'inflammation'),] 
# 1B. Get their gene sets
filename <- paste(modules_dir, 'filtered_modules.tsv',sep='')
modules_df <- read.table(filename, sep='\t', header=TRUE)
modules_gene_sets <- split(x=modules_df$gene_symbol, f=modules_df$module_color)
for (module in names(modules_gene_sets)) {
  if (!(module %in% modules_per_var_df$term)) {
    modules_gene_sets[[module]] <- NULL
  }
}
modules_genes <- unique(unlist(modules_gene_sets))
length(modules_genes)


################################################################################
# 2. Loading samples ranking, metadata and counts matrix
################################################################################
# Samples ranking
sorted_samples_df = read.table(sorted_samples_filename, sep=',', header=TRUE)
rownames(sorted_samples_df) = sorted_samples_df[,1] 
sorted_samples_df[,1] = NULL
sorted_samples_df[,2] = rownames(sorted_samples_df)
colnames(sorted_samples_df) <- c('Sorting_axis', 'Sample_name')

template_df = read.csv(file = template_filename, sep=',')
rownames(template_df) = template_df[,sample_names_col]
template_df <- template_df[rownames(sorted_samples_df),]
counts_matrix_df <- read.csv(file = counts_filename, sep = ",", check.names=FALSE,
                             row.names=1)
counts_matrix_df <- counts_matrix_df[, rownames(sorted_samples_df)]
counts_matrix_df <- counts_matrix_df[(rowSums(counts_matrix_df)>dim(counts_matrix_df)[2]),] #Exclude low expressed counts
counts_matrix = data.matrix(counts_matrix_df)
colnames(counts_matrix) == row.names(template_df)


################################################################################
# 3. Stratification using the the SAF score
################################################################################
template_df[,'SAF_score_division'] = "NASH"

NAFL = which(template_df$Inflammation < 2 & template_df$Ballooning == 0 & template_df$Fibrosis == 0)
template_df[NAFL, 'SAF_score_division'] = "NAFL-NASHF0"

NASH_F0 = which(template_df[,'SAF_score_division'] == "NASH" & template_df$Fibrosis == 0)
template_df[NASH_F0, 'SAF_score_division'] = "NAFL-NASHF0"

NASH_F1 = which(template_df[,'SAF_score_division'] == "NASH" & template_df$Fibrosis == 1)
template_df[NASH_F1, 'SAF_score_division'] = "NASH_F12"

NASH_F2 = which(template_df[,'SAF_score_division'] == "NASH" & template_df$Fibrosis == 2)
template_df[NASH_F2, 'SAF_score_division'] = "NASH_F12"

NASH_F3 = which(template_df[,'SAF_score_division'] == "NASH" & template_df$Fibrosis == 3)
template_df[NASH_F3, 'SAF_score_division'] = "NASH_F34"

NASH_F4 = which(template_df[,'SAF_score_division'] == "NASH" & template_df$Fibrosis == 4)
template_df[NASH_F4, 'SAF_score_division']= "NASH_F34"

CTRL = which(template_df$SAF.Score..check..Borderline..NAFL. == "CTRL")
template_df[CTRL, 'SAF_score_division'] = "CONTROL"

# CONTROL NAFL-NASHF0    NASH_F12    NASH_F34 
# 4          28          73          30 


################################################################################
# 4. Run the analysis for the SAF-based stratification
################################################################################
conditions <- list(
  'Control'= c('CONTROL'),
  'Mild'= c('NAFL-NASHF0'),
  'Moderate'= c('NASH_F12'),
  'Severe'= c('NASH_F34')
)

sorted_conditions <- c('Control', 'Mild', 'Moderate', 'Severe')
samples_in_conditions <- list()
for (c in names(conditions)) {
  tmp <- template_df[template_df[,'SAF_score_division'] %in% conditions[[c]],]
  samples_in_conditions[[c]] <- rownames(tmp)
}

output_dir = '../../results/ucam_sanyal/discrete_analysis/SAF_score/'
deseq_worker(sorted_conditions, samples_in_conditions, output_dir)

de_analysis_df <- read.csv(paste(output_dir,'de_results.csv', sep=''), row.names=1, sep=',')
de_analysis_df <- de_analysis_df[!de_analysis_df$external_gene_name == "",]
dim(de_analysis_df)

columns <- colnames(de_analysis_df)[grepl(x=colnames(de_analysis_df), pattern='log2')]
logfc_df <- de_analysis_df[,c(columns,'external_gene_name')]
logfc_df[is.na(logfc_df)] <- 0 
logfc_df <- logfc_df[!logfc_df$external_gene_name == "",]
logfc_df <- logfc_df %>% 
  group_by(external_gene_name) %>% 
  select(-external_gene_name) %>% 
  summarise_if(is.numeric, mean)
logfc_df <- as.data.frame(logfc_df)
rownames(logfc_df) <- logfc_df$external_gene_name
logfc_df$external_gene_name <- NULL
colnames(logfc_df) <- c('Mild_vs_Control', 'Moderate_vs_Mild', 'Severe_vs_Moderate')
dim(logfc_df)

columns <- colnames(de_analysis_df)[grepl(x=colnames(de_analysis_df), pattern='padj')]
adj_pvalues_df <- de_analysis_df[,c(columns,'external_gene_name')]
adj_pvalues_df[is.na(adj_pvalues_df)] <- 1
adj_pvalues_df <- adj_pvalues_df[!adj_pvalues_df$external_gene_name == "",]
adj_pvalues_df <- adj_pvalues_df %>% 
  group_by(external_gene_name) %>% 
  select(-external_gene_name) %>% 
  summarise_if(is.numeric, pvalue_transformation)
adj_pvalues_df <- as.data.frame(adj_pvalues_df)
rownames(adj_pvalues_df) <- adj_pvalues_df$external_gene_name
adj_pvalues_df$external_gene_name <- NULL
colnames(adj_pvalues_df) <- c('Mild_vs_Control', 'Moderate_vs_Mild', 'Severe_vs_Moderate')
dim(adj_pvalues_df)

msviper_worker(adj_pvalues_df, logfc_df, output_dir)


################################################################################
# 5. Stratification using the the NAS score
################################################################################
template_df[,'NAS_score_division'] = template_df[,'NAS']

conditions <- list(
  'NAS_0'= c(0),
  'NAS_1'= c(1),
  'NAS_2'= c(2),
  'NAS_3'= c(3),
  'NAS_4'= c(4),
  'NAS_5'= c(5),
  'NAS_6'= c(6),
  'NAS_7'= c(7)
)

sorted_conditions <- paste('NAS', seq(0,7), sep='_')
samples_in_conditions <- list()
for (c in names(conditions)) {
  tmp <- template_df[template_df[,'NAS_score_division'] %in% conditions[[c]],]
  samples_in_conditions[[c]] <- rownames(tmp)
}


################################################################################
# 6. Run the analysis for the SAF-based stratification
################################################################################
output_dir = '../../results/ucam_sanyal/discrete_analysis/NAS_score/'
deseq_worker(sorted_conditions, samples_in_conditions, output_dir)

de_analysis_df <- read.csv(paste(output_dir,'de_results.csv', sep=''), row.names=1, sep=',')
de_analysis_df <- de_analysis_df[!de_analysis_df$external_gene_name == "",]
dim(de_analysis_df)

columns <- colnames(de_analysis_df)[grepl(x=colnames(de_analysis_df), pattern='log2')]
logfc_df <- de_analysis_df[,c(columns,'external_gene_name')]
logfc_df[is.na(logfc_df)] <- 0 
logfc_df <- logfc_df[!logfc_df$external_gene_name == "",]
logfc_df <- logfc_df %>% 
  group_by(external_gene_name) %>% 
  select(-external_gene_name) %>% 
  summarise_if(is.numeric, mean)
logfc_df <- as.data.frame(logfc_df)
rownames(logfc_df) <- logfc_df$external_gene_name
logfc_df$external_gene_name <- NULL
colnames(logfc_df) <- sapply(colnames(logfc_df), function(c) {
  gsub('_log2FoldChange', '', c)
}, USE.NAMES=FALSE)

columns <- colnames(de_analysis_df)[grepl(x=colnames(de_analysis_df), pattern='padj')]
adj_pvalues_df <- de_analysis_df[,c(columns,'external_gene_name')]
adj_pvalues_df[is.na(adj_pvalues_df)] <- 1
adj_pvalues_df <- adj_pvalues_df[!adj_pvalues_df$external_gene_name == "",]
adj_pvalues_df <- adj_pvalues_df %>% 
  group_by(external_gene_name) %>% 
  select(-external_gene_name) %>% 
  summarise_if(is.numeric, pvalue_transformation)
adj_pvalues_df <- as.data.frame(adj_pvalues_df)
rownames(adj_pvalues_df) <- adj_pvalues_df$external_gene_name
adj_pvalues_df$external_gene_name <- NULL
colnames(adj_pvalues_df) <- sapply(colnames(adj_pvalues_df), function(c) {
  gsub('_padj', '', c)
}, USE.NAMES=FALSE)

msviper_worker(adj_pvalues_df, logfc_df, output_dir)






