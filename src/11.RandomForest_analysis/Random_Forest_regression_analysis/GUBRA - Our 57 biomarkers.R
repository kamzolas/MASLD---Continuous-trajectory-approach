library(readxl)
library(randomForest)
library(caret)

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

#Gubra data
gubra=read.csv("data/Gubra_dataset/Normalised_Counts_GUBRA.txt", sep = "\t")
common_ = intersect(gov, rownames(gubra))

gubra = gubra[ common_, ]
gubra <- scale(gubra, center = TRUE, scale = TRUE)

corrected_ucamsanyal = corrected_ucamsanyal[ common_, ]


gubra = gubra[common_ ,]

#data variable for the gubra data
order_gubra = read.csv("data/Gubra_dataset/patient position on x-axis gubra's trajectory - 145 genes trajectory.csv")
order_gubra[,2] <- (order_gubra[,2] - min(order_gubra[,2])) / (max(order_gubra[,2]) - min(order_gubra[,2])) #0-1 normalisation the PC1 column

data2 = gubra
data2 = rbind(data2, order_gubra$x)
rownames(data2)[dim(data2)[1]] <- "PC1"
data2 <- t(data2)
write.table(data2, "RF_data.txt")
data2 <- read.table(file = "RF_data.txt")
  
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

for (variable in 1:1) { #We test in the independent gubra dataset, so here we need to run only once
  print(variable)
  #set.seed(Sys.time())
  set.seed(3)
  
  #sample <- sort(sample.int(n = nrow(data), size = floor(1 * nrow(data)), replace = F)) #100% of our training dataset is used to train the RF model
  #train <- data[sample, ]
  train = data
  #test <- data[-sample, ]
  test = data2[, ]
  
  # Random Forest regression model
  rf <- randomForest(PC1 ~ ., data = train, importance = TRUE)
  
  # Store variable importance
  importance_list[[variable]] <- importance(rf, type = 1)  # Using IncNodePurity for regression
  
  # Predict on the test set
  pred <- predict(rf, test)
  
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
write.table(top_genes, "Top_Genes_by_Importance.txt", sep = "\t", col.names = NA)





gubra_metadata = read.csv("data/Gubra_dataset/GSE126848_meta.txt")
rownames(gubra_metadata) = gubra_metadata$X
#gubra_metadata = gubra_metadata[ -which(gubra_metadata$X == "X0874"), ] #exclude trajectory outlier

# Update metadata
gubra_metadata[ , 1] = 1-test$PC1
gubra_metadata[ , 3] = pred
colnames(gubra_metadata)[c(1,3)] = c("real trajectory", "predicted trajectory")

# Sort by predicted values and write to CSV
pred = rev(sort(pred))
sorted_names = names(pred)
gubra_metadata = gubra_metadata[ sorted_names, ]
write.csv(gubra_metadata, "gubra predicted order.csv")




#pred = pred[ -which(names(pred) == "X0874")]
#test = test[-which(test$PC1 == 1), ]
gubra_metadata$`predicted trajectory` <- (gubra_metadata$`predicted trajectory` - min(gubra_metadata$`predicted trajectory`)) / (max(gubra_metadata$`predicted trajectory`) - min(gubra_metadata$`predicted trajectory`)) #0-1 normalisation the predicted trajectory

# Define colors for each disease group
group_colors <- c("healthy" = "#009E73", 
                  "obese" = "#0072B2", 
                  "NAFLD" = "#E69F00", 
                  "NASH" = "#D55E00")

# Create a vector of colors based on the disease group for each patient
point_colors <- group_colors[gubra_metadata$dis]

pdf("PredictedvsReal_trajectory_.pdf")
# Create the plot with colored points
plot(gubra_metadata$`predicted trajectory`, gubra_metadata$`real trajectory`, 
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
mse_per_patient <- (gubra_metadata$`predicted trajectory` - gubra_metadata$`real trajectory`)^2
gubra_metadata$`MSE` = mse_per_patient

# Create a sequence of patient names
patient_names <- paste0("Patient_", seq_len(nrow(gubra_metadata)))

# Create a bar plot with smaller bars to fit all 57 patients without labels
par(mar = c(8, 5, 4, 2))  # Adjust margins: more space at the bottom for names

pdf("Predicted_trajectory_MSEs.pdf")
bar_positions <- barplot(gubra_metadata$MSE, 
                         col = point_colors[order(gubra_metadata$`predicted trajectory`)],  # Color bars by disease group
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
     cex = 0.5,                 # Shrink label size
     col = point_colors[order(gubra_metadata$`predicted trajectory`)]  # Color the labels
)
dev.off()
# # Add a legend to the plot
# legend("topright", legend = names(group_colors), 
#        col = group_colors, pch = 16, title = "Disease Group")

# Save the sorted metadata to CSV
write.csv(gubra_metadata, "gubra_predicted_order_with_mse.csv")







# Plot MSE distribution
pdf("MSE_Distribution_PC1.pdf")
hist(mse_list, main = "Distribution of MSE (PC1 Prediction)", xlab = paste0("Mean Squared Error: ", round(mean(mse_list), 3)), col = "lightblue", breaks = 10)
dev.off()

# Plot R² distribution
pdf("R2_Distribution_PC1.pdf")
hist(cor(gubra_metadata$`predicted trajectory`, gubra_metadata$`real trajectory`), main = "Distribution of R² (PC1 Prediction)", xlab = paste0("R²: ", round(mean(cor(gubra_metadata$`predicted trajectory`, gubra_metadata$`real trajectory`))*100, 1), " %"), col = "lightgreen", breaks = 10)
dev.off()

# Print summary statistics
cat("Mean MSE:", mean(mse_list), "\n")
cat("Mean R²:", mean(cor(gubra_metadata$`predicted trajectory`, gubra_metadata$`real trajectory`)), "\n")




