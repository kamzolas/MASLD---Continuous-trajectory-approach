suppressMessages(library(broom)) # 1.0.4 (tidy)
suppressMessages(library(modelr)) # 0.1.11 (crossv_kfold)


################################################################################
# Get the input arguments
################################################################################
args = commandArgs(trailingOnly=TRUE)
deep_split = args[1] # for example 2
min_size = args[2] # for example 30

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
regression_data <- data.frame(MEs, pc1 = pc_df$PC1[match(rownames(MEs), pc_df$X)], patients_metadata)


################################################################################
# Split the dataset into training and test sets
################################################################################
set.seed(1234)
samples <- pc_df$X
samples <- sample(samples)

kfolds = seq(1,4)
cv <- crossv_kfold(regression_data, k=4) # train = 0.75 and test = 0.25

################################################################################
# Run the cross validation (k=4)
################################################################################
for (fold in kfolds) {

	indexes <- as.integer(cv$train[[fold]])
	training_samples <- samples[indexes]

	indexes <- as.integer(cv$test[[fold]])
	testing_samples <- samples[indexes]

	##############################################################################
	# Run linear regression for each variable
	##############################################################################
	coefficient_sum <- data.frame()
	variables <- c('Fibrosis', 'Steatosis', 'Ballooning', 'NAS')
	training_residuals_list <- list()
	testing_residuals_list <- list()
	for (variable in variables) {

		tmp_df <- data.frame(regression_data[,grep("ME", colnames(regression_data))], 
		                     dependent = regression_data[,variable],
		                     age = regression_data$Age,
		                     inflammation = regression_data$Inflammation)
		tmp_df <- data.frame(scale(tmp_df))
  
		training_data <- tmp_df[training_samples,]
		testing_data <- tmp_df[testing_samples,]

		lm_total <- lm(dependent~., data = training_data)
		step_total <- step(lm_total, direction = "backward", trace = F)
		lm_improved <- lm(dependent~., data = step_total$model)
		results_lm <- summary(lm_improved)
		tidy_lm <- data.frame(tidy(lm_improved, conf.int = T))
		coefficient_sum <- rbind(coefficient_sum, data.frame(tidy_lm, variable = variable))
  
		# training error
		errors_df <- as.data.frame((lm_improved$residuals)^2)
		colnames(errors_df) <- c(variable)
		training_residuals_list[[variable]] <- errors_df

		# test error
		validation_df <- data.frame(predictions = predict(lm_improved, testing_data), 
		                            true_value=testing_data$dependent)
		errors_df <- as.data.frame((validation_df$true_value - validation_df$predictions)^2)
		row.names(errors_df) <- row.names(validation_df)
		colnames(errors_df) <- c(variable)
		testing_residuals_list[[variable]] <- errors_df
	}

	module_coefs_df <- coefficient_sum[coefficient_sum$term != "(Intercept)",]
	module_coefs_df$p.value_adjust <- p.adjust(module_coefs_df$p.value, method="BH")


	################################################################################
	# Save the models (module coefficients)
	################################################################################
	filename = paste(output_dir, 'module_variable_coefficients_', as.character(fold), '.tsv', sep='')
	write.table(x = module_coefs_df, file = filename, sep = '\t', quote = FALSE, row.names = FALSE)


	##############################################################################
	# Merge the errors from the 4 different models and save them into two data 
	# frame, one for the training and another one for the test
	##############################################################################
	merge_dfs <- function(x,y) {
		df <- merge(x, y, by=0, all=TRUE)
		row.names(df) <- df$Row.names
		df <- df[, !names(df) %in% c("Row.names")]
	}

	# training
	residuals_df <- Reduce(function(x,y) {merge_dfs(x,y)}, training_residuals_list)
	filename = paste(output_dir, 'training_squared_residuals_', as.character(fold), '.tsv', sep='')
	write.table(x = residuals_df, file = filename, sep = '\t', quote = FALSE, row.names = TRUE)

	# test
	residuals_df <- Reduce(function(x,y) {merge_dfs(x,y)}, testing_residuals_list)
	filename = paste(output_dir, 'testing_squared_residuals_', as.character(fold), '.tsv', sep='')
	write.table(x = residuals_df, file = filename, sep = '\t', quote = FALSE, row.names = TRUE)

}
