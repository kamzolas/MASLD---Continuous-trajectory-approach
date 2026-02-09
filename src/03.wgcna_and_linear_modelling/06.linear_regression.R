suppressMessages(library(broom)) # 1.0.10 (tidy)

################################################################################ 
# Description
################################################################################
# In this step the final linear models of MASLD variables are created. The 
# optimal set of WGCNA parameters has been defined in the previous steps. These
# optimal values of deep_split and min_size need to be defined in the input
# arguments by the user in order to create the final models.
# Outputs:
# - A tsv file with the coefficients of linear models.
################################################################################


################################################################################
# Inputs
################################################################################
args = commandArgs(trailingOnly=TRUE)
deep_split = args[1] # for example 2
min_size = args[2] # for example 30
key <- paste(as.character(deep_split), as.character(min_size), sep='_')
output_dir <- paste('../../results/ucam_sanyal/wgcna_and_linear_modelling/grid_params/', 
                    key, '/', sep='')
patients_metadata <- read.csv("../../data/ucam_sanyal/metadata.csv")
patients_df <- read.delim("../../data/ucam_sanyal/PC1_sorted_samples.csv", 
                          header=T, sep=",")
colnames(patients_df) <- c('X', 'PC1', 'PC2')
samples_name <- 'Sample.name'


################################################################################
# 1. Load geneTree and eigen-genes objects
################################################################################
load(paste(output_dir, 'eigengenes.RData', sep=''))


################################################################################
# 2. Create the appropriate data frame for the linear modeling by merging the
# expression profile of eigen-genes and patients MASLD phenotype measurements
################################################################################
MEs <- MEList$eigengenes
patients_metadata <- patients_metadata[patients_metadata[,samples_name] %in% patients_df$X,]
regression_data <- data.frame(MEs, pc1 = patients_df$PC1[match(rownames(MEs), patients_df$X)], patients_metadata)


################################################################################
# 3. Run linear regression for each variable
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
# 4.Save the models (module coefficients)
################################################################################
filename = paste(output_dir, 'module_variable_coefficients_FINAL.tsv', sep='')
write.table(x = module_coefs_df, file = filename, sep = '\t', 
	    quote = FALSE, row.names = FALSE)



