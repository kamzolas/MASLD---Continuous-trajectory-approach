suppressMessages(library(Rcpp)) # ‘1.0.12’
suppressMessages(library(igraph)) # ‘1.4.3’
suppressMessages(library(GOSemSim)) # ‘2.29.2’

################################################################################
#
# 1. Input arguments
#
################################################################################
main_dir = '../../results/networks/'
output_dir = '../../results/networks/semantic_similarities/'
dir.create(output_dir, showWarnings = FALSE)


################################################################################
#
# 2. Calculation of three different semantic similarity scores, for all the genes
# included in the file 'MASLD_nodes_without_semantics_filtered.csv'.
#
################################################################################

# A. Loading the the gene names (entrez ids).
df <- read.csv(paste(main_dir, 'MASLD_nodes_without_semantics_filtered.csv', sep=''))
ids <- as.character(df$entrez_id)

# B. Calculation of similarities. The results are stored in different squared 
# matrices.
ontologies <- c("BP")
measures <- c("Resnik", "Lin", 'Wang')
for (ont in ontologies) {
	hsGO <- godata('org.Hs.eg.db', ont=ont)
	for (measure in measures) {
		print(paste(ont, measure))
		matrix = mgeneSim(ids, hsGO, measure=measure)
		matrix <- as.data.frame(matrix)
		write.table(matrix, file = paste(output_dir, measure, '_', ont, '.tsv', sep=''), sep=',')
	}
}
