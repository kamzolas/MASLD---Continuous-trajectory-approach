suppressMessages(library(dplyr)) # 1.1.2
suppressMessages(library(tidyr)) # 1.3.0
suppressMessages(library(viper)) # 1.34.0
suppressMessages(library(DESeq2))
suppressMessages(library(purrr)) # 1.0.1
source("library.R")


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
# Transform regulon data frame to regulon object (create by Giannis)
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
# Calculate transcription factor activities given a specific signature of gene
# scores
################################################################################
calculate_tf_activities <- function(signature) {
  mrs <- msviper(ges = signature, regulon = df2regulon(Regulon_file),
                ges.filter = F, minsize = 4)
  mrs_df <- data.frame('size' = mrs$es$size,
                       'p.value' = mrs$es$p.value,
                       'fdr' = p.adjust(mrs$es$p.value, method = 'fdr'),
                       'nes' = mrs$es$nes)
  return(mrs_df)
}





################################################################################
# Data and results directory
################################################################################
args = commandArgs(trailingOnly=TRUE)
deep_split = 4#args[1] # 2
min_size = 180#args[2] # 60
key = paste(deep_split, min_size, sep='_')
modules_dir = paste('../results/', key, '/', sep='')

output_dir = '../results/discrete_analysis/'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)



# Samples metadata
template_df = read.csv(file = "../data/metadata.csv")[, -c(1,3)]
rownames(template_df) = template_df$Sample.name

# Counts matrix
cts_ucamsanyal <- read.csv(file = "../data/merged_counts.csv", sep = ",")
rownames(cts_ucamsanyal) = as.character(cts_ucamsanyal$X)
cts_ucamsanyal = cts_ucamsanyal[,-c(1)]
cts_ucamsanyal[,'Sample.5'] <- NULL
colnames(cts_ucamsanyal) = template_df$Sample.name
cts_ucamsanyal <- cts_ucamsanyal[(rowSums(cts_ucamsanyal)>dim(cts_ucamsanyal)[2]),] #Exclude low expressed counts
unique(colnames(cts_ucamsanyal) == template_df$Sample.name)
cts = data.matrix(cts_ucamsanyal)
colnames(cts) == row.names(template_df)



template <- template_df
template$SAF.score2 = "NASH"
colnames(template)[13] = "NEW DIVISION GROUP"

NAFL = which(template$Inflammation < 2 & template$Ballooning == 0 & template$Fibrosis == 0)
template$'NEW DIVISION GROUP'[NAFL] = "NAFL-NASHF0"

NASH_F0 = which(template$'NEW DIVISION GROUP' == "NASH" & template$Fibrosis == 0)
template$'NEW DIVISION GROUP'[NASH_F0] = "NAFL-NASHF0"
NASH_F1 = which(template$'NEW DIVISION GROUP' == "NASH" & template$Fibrosis == 1)
template$'NEW DIVISION GROUP'[NASH_F1] = "NASH_F12"
NASH_F2 = which(template$'NEW DIVISION GROUP' == "NASH" & template$Fibrosis == 2)
template$'NEW DIVISION GROUP'[NASH_F2] = "NASH_F12"
NASH_F3 = which(template$`NEW DIVISION GROUP` == "NASH" & template$Fibrosis == 3)
template$'NEW DIVISION GROUP'[NASH_F3] = "NASH_F34"
NASH_F4 = which(template$'NEW DIVISION GROUP' == "NASH" & template$Fibrosis == 4)
template$'NEW DIVISION GROUP'[NASH_F4] = "NASH_F34"

CTRL = which(template$SAF.Score..check..Borderline..NAFL. == "CTRL")
template$`NEW DIVISION GROUP`[CTRL] = "CONTROL"



# Patients stratification

patient_classes <- list(
  'Control'= c('CONTROL'),
  'Mild'= c('NAFL-NASHF0'),
  'Moderate'= c('NASH_F12'),
  'Severe'= c('NASH_F34')
)

sorted_classes <- c('Control', 'Mild', 'Moderate', 'Severe')

patient_groups <- list()
for (c in names(patient_classes)) {
  tmp <- template_df[template$'NEW DIVISION GROUP' %in% patient_classes[[c]],]
  patient_groups[[c]] <- rownames(tmp)
}



# DESeq2
deseq_results <- list()

for (idx1 in seq(length(sorted_classes)-1)) {
  idx2 <- idx1 + 1
  samples1 <- patient_groups[[sorted_classes[idx1]]]
  samples2 <- patient_groups[[sorted_classes[idx2]]]

  # Subset the counts matrix based on the input samples
  cts1 <- cts[,samples1,drop=FALSE]
  cts2 <- cts[,samples2,drop=FALSE]
  
  # Customize the template and cts for sw1_samples
  template_df1 <- template[samples1, c('Dataset', 'Sex'), drop=FALSE]
  key1 <- 'SW1'
  template_df1$sample_id <- paste(key1, seq(1,length(samples1)), sep='_')
  template_df1$condition <- rep(key1, length(samples1))
  template_df1 <- adjust_batch(template_df1)
  colnames(cts1) = template_df1$sample_id
  
  # Customize the template and cts for sw2_samples
  template_df2 <- template[samples2, c('Dataset', 'Sex'), drop=FALSE]
  key2 <- 'SW2'
  template_df2$sample_id <- paste(key2, seq(1,length(samples2)), sep='_')
  template_df2$condition <- rep(key2, length(samples2))
  template_df2 <- adjust_batch(template_df2)
  colnames(cts2) = template_df2$sample_id
  
  # Concatenate the count matrices and templates
  counts_matrix = cbind(cts1, cts2)
  coldata = rbind(template_df1, template_df2)
  rownames(coldata) <- NULL
  coldata$condition <- as.factor(coldata$condition)
  
  # Factorization of variables
  coldata$batch_id <- paste0(coldata$Dataset, coldata$Sex)
  coldata$batch_id <- as.factor(coldata$batch_id)
  
  # Run DESeq and get the results
  dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                colData = coldata,
                                design = ~ batch_id + condition)
  dds$condition <- relevel(dds$condition, ref = key1)
  dds <- DESeq(dds, parallel=FALSE)
  res <- results(dds, cooksCutoff = FALSE)
  res <- as.data.frame(res)
  res <- res[!is.na(res$padj),]
  res <- res[,c('log2FoldChange', 'pvalue', 'padj')]
  key <- paste(sorted_classes[idx2], '_vs_', sorted_classes[idx1], sep='')
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

ensembl_mapping <- read.table('../data/ensembl_mapping.tsv', sep='\t')
ensembl_mapping[,1] <- NULL
colnames(ensembl_mapping) <- c('ensembl_gene_id', 'external_gene_name')

final_df <- merge(final_df, ensembl_mapping, by.x=0, by.y='ensembl_gene_id',
                  all.x=TRUE, all.y=FALSE)
colnames(final_df)[1] <- 'ensembl_gene_id'

write.table(final_df, file = '../results/discrete_analysis/deseq_results.csv', sep=',')

















################################################################################
#
# 1. Load the modules, keep only the significant ones based on regression 
# analysis results and then retrieve their filtered gene sets. Finally, create
# a pool of genes, unifying the module gene sets. Only this pool of genes will 
# be used in the downstream analysis 
#
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




################################################################################
#
# 2.Load the results from differential expression analysis and transform them
# appropriately
#
#
# 2A. Get only the log2fc values and remove duplicated gene symbols by averaging
# their logfcs
#
# 2B. Get only the adjusted pvalues and remove duplicated gene symbols. To 
# do that, calculate the mean log10(adj.pvalue) and create a new row
#
################################################################################
de_analysis_df <- read.csv('../results/discrete_analysis/deseq_results.csv', 
                           row.names=1, sep=',')
de_analysis_df <- de_analysis_df[!de_analysis_df$external_gene_name == "",]


# 2A
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

# 2B
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




################################################################################
#
# 3. Load and create some necessary objects for the TF activity analysis
#
# 3A. The Regulons' object
#
# 3B. Filter the above logfc and adj_pvalues dataframes, keeping the results only
# for the modules_genes
#
################################################################################

# 3A
Regulon_file<- read.csv("../data/collectTRI_network.tsv", sep='\t', header=T)
background_genes <- intersect(modules_genes, rownames(adj_pvalues_df))
length(background_genes)
final_adj_pvalues_df <- adj_pvalues_df[background_genes, ]
final_logfc_df <- logfc_df[background_genes, ]




################################################################################
#
# 4. Perform the TF analysis for each sw and store all the results in three formats:
# 
# A. A list of dataframes, one for each sw
# B. A unified dataframe, which will be used later for the network analysis
# C. A unified dataframe for the FDRs of TF activities
################################################################################

# 4A
discrete_stages = colnames(final_logfc_df)
msviper_results_list <- list()
for (stage in discrete_stages) {
    # Step 1: construct the gene signature for msviper
    stage_df <- merge(final_adj_pvalues_df[,stage,drop=FALSE], 
                      final_logfc_df[,stage,drop=FALSE],
                      by.x=0, by.y=0)
    colnames(stage_df) <- c('gene_symbol', 'padj', 'logfc')
    rownames(stage_df) <- stage_df$gene_symbol
    stage_df$gene_symbol <- NULL
    signature <- stage_df[,'padj',drop=TRUE]
    signature <- qnorm(signature/2, lower.tail = FALSE)*sign(stage_df[,'logfc',drop=TRUE])
    names(signature) <- rownames(stage_df)
    # Step 2: Run msviper and save the results
    msviper_results_df <- calculate_tf_activities(signature)
    msviper_results_list[[stage]] <- msviper_results_df
}

# 4B
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

filename <- paste('../results/discrete_analysis/', 'msviper_results.tsv', sep='')
write.table(msviper_results_df, file = filename, quote = FALSE, sep='\t')


