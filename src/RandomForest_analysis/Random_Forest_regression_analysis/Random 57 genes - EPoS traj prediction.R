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

#Load our top145
#top145 = read.csv("/Users/kamzolas_macbookpro/Desktop/NAFLD -> T2DM/Slingshot - Pseudotemporal ordering/SANYAL:UCAM Slingshot after correction/Random Forest to Predict NAS components/Best gene signature - SAF/145 (from200) most important genes (variable importance with MeanDecreaseAccuracy>0).csv")[,-1]

# Load your gene signature
set.seed(123)  # Set seed for reproducibility, optional
gov_194 <- sample(1:nrow(corrected_ucamsanyal), 57)
#gov_194 <- read_excel("../../../../RF classifier on distinct groups/Our Top145 vs Govaeres Top_Proteomics-Transcriptomics/Govaeres Top_Proteomics-Transcriptomics.xlsx")
#gov_194 <- gov_194$ENSID_TOP
gov_194 <- rownames(corrected_ucamsanyal)[gov_194]
# corrected_ucamsanyal <- corrected_ucamsanyal[gov,]
#gov = top145$Sample.1
#gov <- gov[is.element(gov, rownames(corrected_ucamsanyal)) & !duplicated(gov)]

# Add PC1 values to your dataset
order_ = read.csv("data/PC1_sorted_samples.csv")
order_[,2] <- (order_[,2] - min(order_[,2])) / (max(order_[,2]) - min(order_[,2])) #0-1 normalisation the PC1 column
rownames(order_) = order_$X
corrected_ucamsanyal = corrected_ucamsanyal[ , rownames(order_)]
rownames(template) = template$Sample.name
template <- template[colnames(corrected_ucamsanyal),]
template$PC1 = order_$PC1

#EPoS metadata
epos_metadata = read.csv("data/EPoS_dataset/epos_metadata.csv")
rownames(epos_metadata) = epos_metadata$X

epos_metadata$NEW.DIVISION.GROUP[epos_metadata$HistoGroup == "control"] <- "CTRL"
epos_metadata$NEW.DIVISION.GROUP[epos_metadata$HistoGroup == "Steatosis"] <- "MASL"
epos_metadata$NEW.DIVISION.GROUP[epos_metadata$HistoGroup == "NASH_F0" | epos_metadata$HistoGroup == "NASH_F1"] <- "MASHF01"
epos_metadata$NEW.DIVISION.GROUP[epos_metadata$HistoGroup == "NASH_F2"] <- "MASHF2"
epos_metadata$NEW.DIVISION.GROUP[epos_metadata$HistoGroup == "NASH_F3" | epos_metadata$HistoGroup == "cirrhosis"] <- "MASHF34"


#EPoS data
epos=read.csv("data/EPoS_dataset/Corrected_NormalisedCounts_EPOS2022.csv")
rownames(epos) = epos$X
epos = epos[,-1]
colnames(epos) = rownames(epos_metadata)

epos_metadata = epos_metadata[ -which(epos_metadata$NEW.DIVISION.GROUP == "CTRL") , ]
epos = epos[ , rownames(epos_metadata)]

#Do the scale (z-score) based only on the 194 genes
gov_194 = gov_194[is.element(gov_194, rownames(epos))]
epos_all = epos
epos = epos[ gov_194, ]
table(is.element(rownames(epos), gov_194))
epos <- scale(epos, center = TRUE, scale = TRUE)
common_ = intersect(gov_194, rownames(epos))
epos = epos[ common_, ]
corrected_ucamsanyal = corrected_ucamsanyal[ common_, ]


#data variable for the epos data
pca = prcomp(t(epos_all))
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1) #How much variation in the original data each PC accounts for

order_epos = pca$x[,1]
order_epos <- (order_epos - min(order_epos)) / (max(order_epos) - min(order_epos)) #0-1 normalisation the PC1 column

data2 = epos
data2 = rbind(data2, order_epos)
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

for (variable in 1:1) { #We test in the independent epos dataset, so here we need to run only once
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

# Plot MSE distribution
pdf("MSE_Distribution_PC1.pdf")
hist(mse_list, main = "Distribution of MSE (PC1 Prediction)", xlab = paste0("Mean Squared Error: ", round(mean(mse_list), 3)), col = "lightblue", breaks = 10)
dev.off()

# Plot R² distribution
pdf("R2_Distribution_PC1.pdf")
hist(r2_list, main = "Distribution of R² (PC1 Prediction)", xlab = paste0("R²: ", round(mean(r2_list)*100, 1), " %"), col = "lightgreen", breaks = 10)
dev.off()

# Print summary statistics
cat("Mean MSE:", mean(mse_list), "\n")
cat("Mean R²:", mean(r2_list), "\n")







# Update metadata
epos_metadata[ , 1] = test$PC1
epos_metadata[ , 2] = 1-pred
colnames(epos_metadata)[c(1,2)] = c("real trajectory", "predicted trajectory")

# Sort by predicted values and write to CSV
pred = rev(sort(pred))
sorted_names = names(pred)
epos_metadata = epos_metadata[ sorted_names, ]
write.csv(epos_metadata, "epos predicted order.csv")




#pred = pred[ -which(names(pred) == "X0874")]
#test = test[-which(test$PC1 == 1), ]
epos_metadata$`predicted trajectory` <- (epos_metadata$`predicted trajectory` - min(epos_metadata$`predicted trajectory`)) / (max(epos_metadata$`predicted trajectory`) - min(epos_metadata$`predicted trajectory`)) #0-1 normalisation the predicted trajectory

# Define colors for each disease group
group_colors <- c(#"CTRL" = "#009E73",    # Green
  "MASL" = "#0072B2",    # Blue
  "MASHF01" = "#E69F00", # Orange
  "MASHF2" = "#D55E00",  # Red
  "MASHF34" = "#CC79A7") # Magenta/Pink

# Create a vector of colors based on the disease group for each patient
point_colors <- group_colors[epos_metadata$NEW.DIVISION.GROUP]

pdf("PredictedvsReal_trajectory_.pdf")
# Create the plot with colored points
plot(epos_metadata$`predicted trajectory`, epos_metadata$`real trajectory`, 
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
mse_per_patient <- (epos_metadata$`predicted trajectory` - epos_metadata$`real trajectory`)^2
epos_metadata$`MSE` = mse_per_patient

# Create a sequence of patient names
patient_names <- paste0("Patient_", seq_len(nrow(epos_metadata)))

# Create a bar plot with smaller bars to fit all 57 patients without labels
par(mar = c(8, 5, 4, 2))  # Adjust margins: more space at the bottom for names

pdf("Predicted_trajectory_MSEs.pdf")
bar_positions <- barplot(epos_metadata$MSE, 
                         col = point_colors[order(epos_metadata$`predicted trajectory`)],  # Color bars by disease group
                         names.arg = NULL,  # No patient names for now
                         xlab = "", 
                         #xlab = "Patients (sorted by predicted trajectory)", 
                         ylab = "Mean Squared Error (MSE)", 
                         main = "MSE for Each Predicted Trajectory",
                         ylim = c(0, 0.5),
                         las = 2,            # Rotate x-axis labels vertically for better readability
                         border = NA,        # No borders around bars
                         cex.names = 0.3,    # Shrink x-axis labels to fit more names
                         width = 0.5)        # Reduce the bar width to fit all patients

# Add colored rownames (x-axis labels) manually, using the patient names
text(x = bar_positions, 
     y = par("usr")[3] - 0.05,  # Position below the bars
     labels = patient_names,    # Use the patient names (Patient_1, Patient_2, ...)
     srt = 90,                  # Rotate labels 90 degrees
     adj = 1,                   # Right align
     xpd = TRUE,                # Allow text to be drawn outside the plot region
     cex = 0.2,                 # Shrink label size
     col = point_colors[order(epos_metadata$`predicted trajectory`)]  # Color the labels
)
dev.off()
# # Add a legend to the plot
# legend("topright", legend = names(group_colors), 
#        col = group_colors, pch = 16, title = "Disease Group")

# Save the sorted metadata to CSV
write.csv(epos_metadata, "epos_predicted_order_with_mse.csv")







# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggpubr)

# Use the new metadata table (epos_metadata)
metadata_df <- as.data.frame(epos_metadata[, c("predicted trajectory", "NEW.DIVISION.GROUP")])
colnames(metadata_df) <- c("predicted_trajectory", "disease.ch1")

# Set the order of the groups
metadata_df$disease.ch1 <- factor(metadata_df$disease.ch1, 
                                  levels = c("MASL", "MASHF01", "MASHF2", "MASHF34"))

# Define colors for each disease group
group_colors <- c(#"CTRL" = "#009E73",    # Green
  "MASL" = "#0072B2",    # Blue
  "MASHF01" = "#E69F00", # Orange
  "MASHF2" = "#D55E00",  # Red
  "MASHF34" = "#CC79A7") # Magenta/Pink

# Create the horizontal violin plot with statistical comparisons
ggplot(metadata_df, aes(x = predicted_trajectory, y = disease.ch1, fill = disease.ch1)) +
  geom_violin(trim = TRUE, show.legend = FALSE) +  # Create violin plot without legend
  geom_jitter(width = 0.1, alpha = 0.5, size = 1, show.legend = FALSE) +  # Add jittered points without legend
  scale_fill_manual(values = group_colors) +  # Use the specified colors for groups
  labs(title = "",
       x = "",
       y = "") +
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1),  # Adjust y-axis labels
        plot.margin = margin(1, 1, 1, 1, "cm")) +  # Adjust margins
  stat_compare_means(comparisons = list(#c("healthy", "obese"), 
    #c("MASL", "MASHF01"),
    c("MASL", "MASHF2"), 
    c("MASL", "MASHF34"), 
    c("MASHF01", "MASHF2"),
    c("MASHF01", "MASHF34"),
    c("MASHF2", "MASHF34")
  ), 
  label = "p.signif", hide.ns = TRUE)  # Only show significant comparisons, hide "ns" labels

# Save the plot to a PDF
ggsave("epos_violin_plot_predicted_trajectory_by_disease_group_horizontal.pdf", width = 6, height = 2.2)



#Create the correlations of the different MASLD variable along with the predicted trajectory
metadata_df = as.data.frame(epos_metadata[,c("predicted trajectory","Grade.of.Steatosis", "Hepatocyte.Ballooning", 
                                             "Lobular.inflammation","Fibrosis.stage", "NAFLD.Activity.Score..NAS..Kleiner")])
colnames(metadata_df) = c("Predicted_trajectory_position","Steatosis", "Ballooning", 
                          "Inflammation","Fibrosis", "NAS")
library("ggpubr")
library("cowplot")
xaxis = metadata_df$Predicted_trajectory_position
pdf("EPoS MASLD chars along predicted trajectory.pdf")

mydata = data.frame("xaxis" = xaxis, "Steatosis" = metadata_df$Steatosis)
STEATOSIS <- ggscatter(mydata, x = "xaxis", y = "Steatosis",
                       color = "#a6611a", shape = 20, size = 3, # Points color, shape and size
                       add = "reg.line",  # Add regression line
                       add.params = list(color = "black", fill = "lightgray"), conf.int = TRUE, 
                       cor.coef = TRUE, cor.method = "pearson", 
                       xlab = "Predicted Trajectory Position", ylab = "mean Steatosis", main = "Steatosis")

mydata = data.frame("xaxis" = xaxis, "Inflammation" = metadata_df$Inflammation)
INFLAMMATION <- ggscatter(mydata, x = "xaxis", y = "Inflammation",
                          color = "#dfc27d", shape = 20, size = 3, # Points color, shape and size
                          add = "reg.line",  # Add regression line
                          add.params = list(color = "black", fill = "lightgray"), conf.int = TRUE, 
                          cor.coef = TRUE, cor.method = "pearson", 
                          xlab = "Predicted Trajectory Position", ylab = "mean Inflammation", main = "Inflammation")

mydata = data.frame("xaxis" = xaxis, "Ballooning" = metadata_df$Ballooning)
BALLOONING <- ggscatter(mydata, x = "xaxis", y = "Ballooning",
                        color = "#80cdc1", shape = 20, size = 3, # Points color, shape and size
                        add = "reg.line",  # Add regression line
                        add.params = list(color = "black", fill = "lightgray"), conf.int = TRUE, 
                        cor.coef = TRUE, cor.method = "pearson", 
                        xlab = "Predicted Trajectory Position", ylab = "mean Ballooning", main = "Ballooning")

mydata = data.frame("xaxis" = xaxis, "Fibrosis" = metadata_df$Fibrosis)
FIBROSIS <- ggscatter(mydata, x = "xaxis", y = "Fibrosis",
                      color = "#018571", shape = 20, size = 3, # Points color, shape and size
                      add = "reg.line",  # Add regression line
                      add.params = list(color = "black", fill = "lightgray"), conf.int = TRUE, 
                      cor.coef = TRUE, cor.method = "pearson", 
                      xlab = "Predicted Trajectory Position", ylab = "mean Fibrosis", main = "Fibrosis")

plot_grid(STEATOSIS, INFLAMMATION, BALLOONING, FIBROSIS,
          nrow = 2, ncol = 2, align = "hv")
dev.off()




#Create the correlations of the different MASLD variable along with the real trajectory
metadata_df = as.data.frame(epos_metadata[,c("real trajectory","Grade.of.Steatosis", "Hepatocyte.Ballooning",
                                             "Lobular.inflammation","Fibrosis.stage", "NAFLD.Activity.Score..NAS..Kleiner")])
colnames(metadata_df) = c("Real_trajectory_position","Steatosis", "Ballooning", 
                          "Inflammation","Fibrosis", "NAS")

xaxis = metadata_df$Real_trajectory_position
pdf("EPoS MASLD chars along real trajectory.pdf")

mydata = data.frame("xaxis" = xaxis, "Steatosis" = metadata_df$Steatosis)
STEATOSIS <- ggscatter(mydata, x = "xaxis", y = "Steatosis",
                       color = "#a6611a", shape = 20, size = 3, # Points color, shape and size
                       add = "reg.line",  # Add regression line
                       add.params = list(color = "black", fill = "lightgray"), conf.int = TRUE, 
                       cor.coef = TRUE, cor.method = "pearson", 
                       xlab = "Predicted Trajectory Position", ylab = "mean Steatosis", main = "Steatosis")

mydata = data.frame("xaxis" = xaxis, "Inflammation" = metadata_df$Inflammation)
INFLAMMATION <- ggscatter(mydata, x = "xaxis", y = "Inflammation",
                          color = "#dfc27d", shape = 20, size = 3, # Points color, shape and size
                          add = "reg.line",  # Add regression line
                          add.params = list(color = "black", fill = "lightgray"), conf.int = TRUE, 
                          cor.coef = TRUE, cor.method = "pearson", 
                          xlab = "Predicted Trajectory Position", ylab = "mean Inflammation", main = "Inflammation")

mydata = data.frame("xaxis" = xaxis, "Ballooning" = metadata_df$Ballooning)
BALLOONING <- ggscatter(mydata, x = "xaxis", y = "Ballooning",
                        color = "#80cdc1", shape = 20, size = 3, # Points color, shape and size
                        add = "reg.line",  # Add regression line
                        add.params = list(color = "black", fill = "lightgray"), conf.int = TRUE, 
                        cor.coef = TRUE, cor.method = "pearson", 
                        xlab = "Predicted Trajectory Position", ylab = "mean Ballooning", main = "Ballooning")

mydata = data.frame("xaxis" = xaxis, "Fibrosis" = metadata_df$Fibrosis)
FIBROSIS <- ggscatter(mydata, x = "xaxis", y = "Fibrosis",
                      color = "#018571", shape = 20, size = 3, # Points color, shape and size
                      add = "reg.line",  # Add regression line
                      add.params = list(color = "black", fill = "lightgray"), conf.int = TRUE, 
                      cor.coef = TRUE, cor.method = "pearson", 
                      xlab = "Predicted Trajectory Position", ylab = "mean Fibrosis", main = "Fibrosis")

plot_grid(STEATOSIS, INFLAMMATION, BALLOONING, FIBROSIS,
          nrow = 2, ncol = 2, align = "hv")
dev.off()





# Create the Violin Plots for the different MASLD variables along the predicted trajectory
metadata_df = as.data.frame(epos_metadata[,c("predicted trajectory","Grade.of.Steatosis", 
                                             "Hepatocyte.Ballooning", "Lobular.inflammation",
                                             "Fibrosis.stage", "NAFLD.Activity.Score..NAS..Kleiner")])
colnames(metadata_df) = c("Predicted_trajectory_position","Steatosis", "Ballooning", 
                          "Inflammation","Fibrosis", "NAS")

library("ggpubr")
library("cowplot")

xaxis = metadata_df$Predicted_trajectory_position
metadata_df$Predicted_trajectory_position <- as.factor(metadata_df$Predicted_trajectory_position)

pdf("EPoS_MASLD_Violin_plots_along_predicted_trajectory.pdf")

metadata_df$Steatosis = as.factor(metadata_df$Steatosis)
mydata = data.frame("Trajectory" = xaxis, "Steatosis" = metadata_df$Steatosis)
STEATOSIS <- ggplot(mydata, aes(x=Trajectory, y=Steatosis)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="#a6611a") +
  geom_boxplot(width=0.1) + theme_minimal()

metadata_df$Inflammation = as.factor(metadata_df$Inflammation)
mydata = data.frame("Trajectory" = xaxis, "Inflammation" = metadata_df$Inflammation)
INFLAMMATION <- ggplot(mydata, aes(x=Trajectory, y=Inflammation)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="#dfc27d") +
  geom_boxplot(width=0.1) + theme_minimal()

metadata_df$Ballooning = as.factor(metadata_df$Ballooning)
mydata = data.frame("Trajectory" = xaxis, "Ballooning" = metadata_df$Ballooning)
BALLOONING <- ggplot(mydata, aes(x=Trajectory, y=Ballooning)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="#80cdc1") +
  geom_boxplot(width=0.1) + theme_minimal()

metadata_df$Fibrosis = as.factor(metadata_df$Fibrosis)
mydata = data.frame("Trajectory" = xaxis, "Fibrosis" = metadata_df$Fibrosis)
FIBROSIS <- ggplot(mydata, aes(x=Trajectory, y=Fibrosis)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="#018571") +
  geom_boxplot(width=0.1) + theme_minimal()

metadata_df$NAS = as.factor(metadata_df$NAS)
mydata = data.frame("Trajectory" = xaxis, "NAS" = metadata_df$NAS)
NAS <- ggplot(mydata, aes(x=Trajectory, y=NAS)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="#018571") +
  geom_boxplot(width=0.1) + theme_minimal()

plot_grid(STEATOSIS, INFLAMMATION, BALLOONING, FIBROSIS, NAS,
          nrow = 3, ncol = 2, align = "hv")
dev.off()



