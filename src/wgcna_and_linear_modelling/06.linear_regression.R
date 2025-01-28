suppressMessages(library(broom)) # 1.0.4. tidy


################################################################################
# Get the input arguments
################################################################################
args = commandArgs(trailingOnly=TRUE)
deep_split = args[1] #2
min_size = args[2] #30

key <- paste(as.character(deep_split), as.character(min_size), sep='_')
output_dir <- paste('../../results/wgcna_and_linear_modelling/grid_params/', key, '/', sep='')

################################################################################
# Load geneTree and eigen-genes objects
################################################################################
load(paste(output_dir, 'eigengenes.RData', sep=''))


################################################################################
# Create the appropriate data frame from the regression
################################################################################
MEs <- MEList$eigengenes
patients_metadata <- read.csv("../../data/metadata.csv")
pc_df <- read.delim("../../data/PC1_sorted_samples.csv", header = T, sep = ",")
patients_metadata <- patients_metadata[patients_metadata$Sample.name %in% pc_df$X,]
regression_data<-data.frame(MEs, pc1 = pc_df$PC1[match(rownames(MEs), pc_df$X)], patients_metadata)



################################################################################
# Run linear regression for each variable
################################################################################
coefficient_sum <- data.frame()
variables <- c('Fibrosis', 'Steatosis', 'Ballooning', 'NAS')
for (variable in variables) {
  print(variable)
  tmp_df <- data.frame(regression_data[,grep("ME", colnames(regression_data))], 
                       dependent=regression_data[,variable],
                       age=regression_data$Age,
                       inflammation = regression_data$Inflammation)
  tmp_df <- data.frame(scale(tmp_df))
  
  lm_total <- lm(dependent~., data = tmp_df)
  step_total <- step(lm_total, direction = "backward", trace = F)
  lm_improved <- lm(dependent~., data = step_total$model)
  results_lm <- summary(lm_improved)
  tidy_lm<-data.frame(tidy(lm_improved, conf.int = T))
  coefficient_sum<-rbind(coefficient_sum, data.frame(tidy_lm, variable = variable))
}
module_coefs_df <- coefficient_sum[coefficient_sum$term != "(Intercept)",]
module_coefs_df$p.value_adjust <- p.adjust(module_coefs_df$p.value, method="BH")


################################################################################
# Save the models (module coefficients)
################################################################################
filename = paste(output_dir, 'module_variable_coefficients_FINAL.tsv', sep='')
write.table(x = module_coefs_df, file = filename, sep = '\t', 
	    quote = FALSE, row.names = FALSE)



