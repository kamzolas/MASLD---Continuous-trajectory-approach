library(readxl)
library(randomForest)
library(caret)
library(biomaRt)

# Load and preprocess data (similar to your original script)
template <- read.csv(file = "data/metadata.csv")[-5, -c(1,3)]
corrected_ucamsanyal <- read.csv("data/z_scores.csv")

# Preprocess the data as you did in the original script
rownames(corrected_ucamsanyal) <- as.character(corrected_ucamsanyal$X)
corrected_ucamsanyal <- corrected_ucamsanyal[,-c(1)]
colnames(corrected_ucamsanyal) <- template$Sample.name

# Load your 57 gene signature
gov1= read.csv('data/57_BMs.csv')[,-1]
gov = gov1$ensembl_gene_id[is.element(gov1$ensembl_gene_id, rownames(corrected_ucamsanyal))]
corrected_ucamsanyal <- corrected_ucamsanyal[gov,]

# Add PC1 values to your dataset
order_ = read.csv("data/PC1_sorted_samples.csv")
order_[,2] <- (order_[,2] - min(order_[,2])) / (max(order_[,2]) - min(order_[,2])) #0-1 normalisation the PC1 column
rownames(order_) = order_$X
corrected_ucamsanyal = corrected_ucamsanyal[ , rownames(order_)]
rownames(template) = template$Sample.name
template <- template[colnames(corrected_ucamsanyal),]
template$PC1 = order_$PC1

data <- corrected_ucamsanyal
data[dim(data)[1] + 1,] <- template$PC1
rownames(data)[dim(data)[1]] <- "PC1"
data <- t(data)
write.table(data, "RF_data.txt")
data <- read.table(file = "RF_data.txt")

# Convert PC1 to a numeric variable
data$PC1 <- as.numeric(data$PC1)

# Initialize vectors for storing results
mse_list <- c()
pval_list <- c()
r2_list <- c()
importance_list <- list()

predicted_positions = matrix(NA, nrow = dim(corrected_ucamsanyal)[2], ncol = 1000)
rownames(predicted_positions) = colnames(corrected_ucamsanyal)

for (variable in 1:1000) {
  print(variable)
  set.seed(Sys.time())
  
  # Data partition: 80% training and 20% testing
  sample <- sort(sample.int(n = nrow(data), size = floor(.8 * nrow(data)), replace = F))
  train <- data[sample, ]
  test <- data[-sample, ]
  
  # Random Forest regression model
  rf <- randomForest(PC1 ~ ., data = train, importance = TRUE)
  
  # Store variable importance
  importance_list[[variable]] <- importance(rf, type = 1)  # Using IncNodePurity for regression
  
  # Predict on the test set
  pred <- predict(rf, test)
  predicted_positions[ rownames(test), variable] = pred
  
  # Calculate and store MSE and R²
  mse <- mean((pred - test$PC1)^2)
  
  cor_test <- cor.test(pred, test$PC1, method = "pearson")
  
  r <- cor_test$estimate   # correlation coefficient
  r2 <- r^2                # R²
  pval <- cor_test$p.value # p-value
  
  
  
  mse_list <- c(mse_list, mse)
  r2_list <- c(r2_list, r2)
  pval_list <- c(pval_list, pval)
}

# Combine the importance values across all iterations
importance_df <- do.call(cbind, importance_list)
mean_importance <- rowMeans(importance_df, na.rm = TRUE)

# Sort the genes by their mean importance
sorted_importance <- sort(mean_importance, decreasing = TRUE)

# Output the top genes based on importance
top_genes <- head(sorted_importance, 10)
print("Top 10 Genes by Importance:")
print(top_genes)

# Optional: Write the top genes to a file
write.table(sorted_importance, "Ranked_Genes_by_Importance.txt", sep = "\t", col.names = NA)


temp = predicted_positions
temp = rowMeans(temp, na.rm = TRUE)



metadata = template

# Update metadata
metadata[ , 1] = 1-metadata$PC1
metadata[ , 2] = 1-temp
colnames(metadata)[c(1,2)] = c("real trajectory", "predicted trajectory")

# Sort by predicted values and write to CSV
sorted_names = (order(metadata$`predicted trajectory`))
metadata = metadata[ sorted_names, ]
write.csv(metadata, "ucam_vcu predicted order.csv")

#Group metadata in CTRL, MASL, MASHF01, MASHF2, MASHF34
metadata$SAF.score[metadata$SAF.score == "MASH F0" | metadata$SAF.score == "MASH F1"] <- "MASHF01"
metadata$SAF.score[metadata$SAF.score == "MASH F2"] <- "MASHF2"
metadata$SAF.score[metadata$SAF.score == "MASH F3" | metadata$SAF.score == "MASH F4"] <- "MASHF34"

metadata$`predicted trajectory` <- (metadata$`predicted trajectory` - min(metadata$`predicted trajectory`)) / (max(metadata$`predicted trajectory`) - min(metadata$`predicted trajectory`)) #0-1 normalisation the predicted trajectory

# Define colors for each disease group
group_colors <- c("CTRL" = "#009E73",    # Green
                  "MASL" = "#0072B2",    # Blue
                  "MASHF01" = "#E69F00", # Orange
                  "MASHF2" = "#D55E00",  # Red
                  "MASHF34" = "#CC79A7") # Magenta/Pink

# Create a vector of colors based on the disease group for each patient
point_colors <- group_colors[metadata$SAF.score]

pdf("PredictedvsReal_trajectory_.pdf")
# Create the plot with colored points
plot(metadata$`predicted trajectory`, metadata$`real trajectory`, 
     xlim = c(0, 1),    # Set x-axis from 0 to 1
     ylim = c(0, 1),    # Set y-axis from 0 to 1
     pch = 16,          # Solid circle points
     col = point_colors, # Color points based on disease group
     xlab = "Predicted", # X-axis label
     ylab = "Real",    # Y-axis label
     main = "Predicted vs Real")  # Title of the plot

# Add the red diagonal line
abline(a = 0, b = 1, col = "black", lwd = 2)

# # Add a legend to the plot
# legend("topright", legend = names(group_colors), 
#        col = group_colors, pch = 16, title = "Disease Group")
dev.off()





# Calculate MSE for each predicted value
mse_per_patient <- (metadata$`predicted trajectory` - metadata$`real trajectory`)^2
metadata$`MSE` = mse_per_patient

# Create a sequence of patient names
patient_names <- paste0("Patient_", seq_len(nrow(metadata)))

# Create a bar plot with smaller bars to fit all 57 patients without labels
par(mar = c(8, 5, 4, 2))  # Adjust margins: more space at the bottom for names

pdf("Predicted_trajectory_MSEs.pdf")
bar_positions <- barplot(metadata$MSE, 
                         col = point_colors[order(metadata$`predicted trajectory`)],  # Color bars by disease group
                         names.arg = NULL,  # No patient names for now
                         xlab = "", 
                         #xlab = "Patients (sorted by predicted trajectory)", 
                         ylab = "Mean Squared Error (MSE)", 
                         main = "MSE for Each Predicted Trajectory",
                         ylim = c(0, 0.5),
                         las = 2,            # Rotate x-axis labels vertically for better readability
                         border = NA,        # No borders around bars
                         cex.names = 0.5,    # Shrink x-axis labels to fit more names
                         width = 0.5)        # Reduce the bar width to fit all patients

# Add colored rownames (x-axis labels) manually, using the patient names
text(x = bar_positions, 
     y = par("usr")[3] - 0.05,  # Position below the bars
     labels = patient_names,    # Use the patient names (Patient_1, Patient_2, ...)
     srt = 90,                  # Rotate labels 90 degrees
     adj = 1,                   # Right align
     xpd = TRUE,                # Allow text to be drawn outside the plot region
     cex = 0.25,                 # Shrink label size
     col = point_colors[order(metadata$`predicted trajectory`)]  # Color the labels
)
dev.off()
# # Add a legend to the plot
# legend("topright", legend = names(group_colors), 
#        col = group_colors, pch = 16, title = "Disease Group")

# Save the sorted metadata to CSV
write.csv(metadata, "ucam_vcu_predicted_order_with_mse.csv")





# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggpubr)

# Sample metadata data frame (replace this with your actual metadata)
metadata_df <- as.data.frame(metadata[, c("predicted trajectory", "SAF.score")])
colnames(metadata_df) <- c("predicted_trajectory", "SAF_score")

# Set the order of the groups
metadata_df$SAF_score <- factor(metadata_df$SAF_score, 
                                levels = c("CTRL", "MASL", "MASHF01", "MASHF2", "MASHF34"))

# Define the colors for each group
group_colors <- c("CTRL" = "#009E73",    # Green
                  "MASL" = "#0072B2",    # Blue
                  "MASHF01" = "#E69F00", # Orange
                  "MASHF2" = "#D55E00",  # Red
                  "MASHF34" = "#CC79A7") # Magenta/Pink

# Create the horizontal violin plot with statistical comparisons
ggplot(metadata_df, aes(x = predicted_trajectory, y = SAF_score, fill = SAF_score)) +
  geom_violin(trim = TRUE, show.legend = FALSE) +  # Create violin plot without legend
  geom_jitter(width = 0.1, alpha = 0.5, size = 1, show.legend = FALSE) +  # Add jittered points without legend
  scale_fill_manual(values = group_colors) +  # Use the specified colors for groups
  labs(title = "",
       x = "",
       y = "") +
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1),  # Adjust y-axis labels
        plot.margin = margin(1, 1, 1, 1, "cm")) +  # Adjust margins
  stat_compare_means(comparisons = list(#c("CTRL", "MASL"), 
                                        c("CTRL", "MASHF01"), 
                                        c("CTRL", "MASHF2"), 
                                        c("CTRL", "MASHF34"),
                                        c("MASL", "MASHF01"),
                                        c("MASL", "MASHF2"),
                                        c("MASL", "MASHF34"),
                                        #c("MASHF01", "MASHF2"),
                                        c("MASHF01", "MASHF34")),
                                        #c("MASHF2", "MASHF34")), 
                     label = "p.signif", hide.ns = TRUE, remove = "non_sig")  # Add significance labels and remove non-significant lines

# Save the plot to a PDF
ggsave("violin_plot_predicted_trajectory_by_group_horizontal.pdf", width = 6.75, height = 2.2)



# Plot MSE distribution
pdf("MSE_Distribution_PC1.pdf")
hist(mse_list, main = "Distribution of MSE (PC1 Prediction)", xlab = paste0("Mean Squared Error: ", round(mean(mse_list), 3)), col = "lightblue", breaks = 10)
dev.off()

# Plot R² distribution
pdf("R2_Distribution_PC1.pdf")
hist(r2_list, main = "R² Distribution (Predicted vs Real patient position on the trajectory)", xlab = paste0("Mean R²: ", round(mean(cor(metadata$`real trajectory`, metadata$`predicted trajectory`))*100, 1), " %"), col = "lightgray", breaks = 10)
dev.off()

# Print summary statistics
cat("Mean MSE:", mean(mse_list), "\n")
cat("Mean R²:", mean(cor(metadata$`real trajectory`, metadata$`predicted trajectory`)), "\n")

