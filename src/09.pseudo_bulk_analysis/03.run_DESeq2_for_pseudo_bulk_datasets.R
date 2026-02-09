library(DESeq2) # 1.46.0

################################################################################
# Description
################################################################################
# Apply DESeq2 for the SW-based groups of pseudo bulk patients.
# Outputs:
# - A csv file which contains the results for all the comparisons along the SW-based
# trajectory
################################################################################


data_dir <- '../../results/ucam_sanyal/pseudo_bulk_analysis/pseudo_bulk_RNAseq_datasets/'
files <- list.files(data_dir)
files <- files[grepl('raw_counts', files)]
Ntrials <- length(files)-1

for (N in seq(0, Ntrials, 1)) {

  output_dir <- paste('../../results/ucam_sanyal/pseudo_bulk_analysis/all_de_results/', 
                      as.character(N), '/', sep='')
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  filename <- paste(data_dir, 'raw_counts_', as.character(N), '.tsv', sep='')
  counts_matrix <- read.csv(filename, sep='\t', row.names = 1)
  
  filename <- paste(data_dir, 'metadata_', as.character(N), '.tsv', sep='')
  template <-  read.csv(filename, sep='\t', row.names = 1)
  
  L <- length(unique(template$condition))
  comparisons <- list()
  for (i in 1:(L-1)) {
    sw1 <- paste('SW_', as.character(i), sep='')
    sw2 <- paste('SW_', as.character(i+1), sep='')
    comparisons[[i]] <- c(sw1, sw2)
  }
  
  comparison <- comparisons[[1]]
  for (comparison in comparisons) {
    
    sw1 <- comparison[1]
    sw1_template <- template[template$condition == sw1,]
    sw1_cts <- counts_matrix[,sw1_template$sample,drop=FALSE]
    
    sw2 <- comparison[2]
    sw2_template <- template[template$condition == sw2,]
    sw2_cts <- counts_matrix[,sw2_template$sample,drop=FALSE]
    
    # Concatenate the count matrices and templates
    cts = cbind(sw1_cts, sw2_cts)
    coldata = rbind(sw1_template, sw2_template)
    rownames(coldata) <- NULL
    coldata$condition <- as.factor(coldata$condition)
    
    print(paste(N, sw1, sw2))
    
    # Run DESeq and save the results
    dds <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = coldata,
                                  design = ~ condition)
    dds$condition <- relevel(dds$condition, ref = sw1)
    dds <- DESeq(dds, parallel=FALSE)
    res <- results(dds)
    res <- as.data.frame(res)
    #res <- res[!is.na(res$padj),]
    res$gene_symbol <- row.names(res)
    row.names(res) <- NULL
    res <- res[,c("gene_symbol","baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
    res <- res[order(res$padj, decreasing = FALSE),]
    comparison_str <- paste(comparison[1], comparison[2], sep='_vs_')
    filename <- paste(output_dir, comparison_str, '.tsv', sep='')
    write.table(res, file = filename, row.names=FALSE)
  }
}



