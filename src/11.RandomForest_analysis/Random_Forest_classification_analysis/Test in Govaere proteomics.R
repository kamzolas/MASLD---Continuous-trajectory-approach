

library(readxl)
library(dplyr)

### === SETUP === ###
library(randomForest)
library(caret)
library(ggplot2)
library(pROC)
library(reshape2)
library(readxl)
library(tidyverse)

### === LOAD TRAINING DATA AND RF MODEL === ###
# Load patient metadata
template <- read.csv("data/metadata.csv")[-5, -c(1,3)]

# Load z-score normalised gene expression matrix
corrected_ucamsanyal <- read.csv("data/z_scores.csv", check.names = FALSE)

# Load 57 selected gene mapping table (external names + ENSG IDs)
selected_57 <- read.csv("data/57_BMs.csv")[,-1]
selected_57 = selected_57[ is.element(selected_57$ensembl_gene_id, corrected_ucamsanyal[,1]), ]

# Preprocess training matrix
rownames(corrected_ucamsanyal) <- corrected_ucamsanyal[,1]
corrected_ucamsanyal <- corrected_ucamsanyal[,-1]
corrected_ucamsanyal <- corrected_ucamsanyal[selected_57$ensembl_gene_id, ]
template <- template[match(colnames(corrected_ucamsanyal), template$Sample.name), ]
stopifnot(all(template$Sample.name == colnames(corrected_ucamsanyal)))
expr_57 <- as.data.frame(t(corrected_ucamsanyal))

# Define binary label
expr_57$fibrosis_bin <- factor(ifelse(as.numeric(as.character(template$Fibrosis)) >= 3, "Advanced", "Early"),
                               levels = c("Early", "Advanced"))

# Train RF model with 5-fold CV
set.seed(123)
ctrl <- trainControl(method = "cv",
                     number = 5,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary,
                     savePredictions = "final")
rf_cv <- train(fibrosis_bin ~ ., data = expr_57,
               method = "rf",
               trControl = ctrl,
               metric = "ROC")



### === LOAD GOVAERE DATASET === ###
govaere <- read_excel("data/Plasma_Proteomics_dataset/42255_2023_775_MOESM3_ESM.xlsx")[-c(1,2,3,4), -c(1:6,8)]
govaere <- govaere %>%
   distinct(.[[7]], .keep_all = TRUE) %>%
   rename(sample_id = 1) %>%
   filter(!duplicated(sample_id), !is.na(sample_id)) %>%
   mutate(across(-c(sample_id,2), as.double))
govaere <- as.data.frame(govaere[ ,-c(dim(govaere)[2])])
rownames(govaere) <- govaere$sample_id
govaere <- govaere[ , -c(1,2)]

# Metadata for fibrosis
metadata <- read_excel("data/Plasma_Proteomics_dataset/42255_2023_775_MOESM3_ESM.xlsx")[ -c(1,2) , -c(1:9)]
metadata <- metadata[, colnames(govaere)]
metadata <- as.double(metadata[1,])

### === PREPROCESS GOVAERE DATA === ###
# 1. Identify existing/missing genes
selected_genes <- selected_57$external_gene_name
existing_genes <- intersect(selected_genes, rownames(govaere))
missing_genes <- setdiff(selected_genes, rownames(govaere))

# 2. Subset existing genes
govaere_sub <- govaere[existing_genes, , drop = FALSE]

# 3. Add missing genes as rows of 0
if(length(missing_genes) > 0){
   missing_mat <- matrix(0, nrow = length(missing_genes), ncol = ncol(govaere_sub),
                         dimnames = list(missing_genes, colnames(govaere_sub)))
   govaere_sub <- rbind(govaere_sub, missing_mat)
}

# 4. Z-score normalize, avoiding NaNs for missing genes
govaere_sub_z <- t(apply(govaere_sub, 1, function(x) {
   if (sd(x) == 0) {
      rep(0, length(x))  # keep as 0 if all values identical
   } else {
      scale(x)
   }
}))


# 5. Transpose to samples-as-rows
expr_df <- as.data.frame(t(govaere_sub_z))

# 6. Add fibrosis_bin factor
expr_df$fibrosis_bin <- factor(ifelse(metadata == 1, "Advanced", "Early"),
                               levels = c("Early", "Advanced"))

# 7. Map column names to ENSG IDs
name2ensg <- setNames(selected_57$ensembl_gene_id, selected_57$external_gene_name)
expr_genes_no_fib <- setdiff(colnames(expr_df), "fibrosis_bin")
expr_df <- expr_df %>%
   rename_with(.cols = all_of(expr_genes_no_fib), .fn = ~ name2ensg[.])

# 8. Reorder columns to match training set
train_genes <- setdiff(colnames(expr_57), "fibrosis_bin")
expr_df <- expr_df[, c(train_genes, "fibrosis_bin")]

### === PREDICTIONS AND EVALUATION === ###
pred_rf <- predict(rf_cv$finalModel, newdata = expr_df)
prob_rf <- predict(rf_cv$finalModel, newdata = expr_df, type = "prob")[, "Advanced"]

# Confusion matrix
conf_rf <- confusionMatrix(pred_rf, expr_df$fibrosis_bin)

# ROC curve and AUC
roc_rf <- roc(expr_df$fibrosis_bin, prob_rf)
auc_rf <- auc(roc_rf)

cat("Govaere dataset RF AUC:", auc_rf, "\n")
print(conf_rf)

plot(roc_rf, print.auc = TRUE, col = "darkgreen", lwd = 2,
     main = "Govaere Dataset ROC: 57-Gene Classifier")











### === OPTIMAL THRESHOLD FOR HIGHER SENSITIVITY === ###
# Get ROC with correct positive/negative classes
roc_rf <- roc(expr_df$fibrosis_bin, prob_rf, levels = c("Early", "Advanced"))

# Option 1: Youden's index (balances sensitivity & specificity)
opt_coords <- coords(roc_rf, x = "best", best.method = "youden")
opt_threshold <- opt_coords["threshold"]
opt_threshold <- opt_threshold$threshold  # extract as numeric

# Option 2: Force high sensitivity (≥ 0.9) — uncomment to use this instead
# opt_coords <- coords(roc_rf, x = 0.9, input = "sensitivity")
# opt_threshold <- opt_coords["threshold"]

cat("Optimal threshold for higher sensitivity:", opt_threshold, "\n")

# Predict using the new threshold
pred_rf_opt <- factor(ifelse(prob_rf >= opt_threshold, "Advanced", "Early"),
                      levels = c("Early", "Advanced"))

# Confusion matrix with new threshold
conf_rf_opt <- confusionMatrix(pred_rf_opt, expr_df$fibrosis_bin)
print(conf_rf_opt)

# Plot ROC with optimal threshold point
plot(roc_rf, print.auc = TRUE, col = "darkgreen", lwd = 2,
     main = paste0("Govaere Dataset ROC: 57-Gene Classifier\nOptimal threshold = ", round(opt_threshold, 3)))
points(opt_coords["threshold"], opt_coords["sensitivity"], col = "red", pch = 19)







#####
#####
#####
library(ggplot2)
library(cowplot)
library(pROC)
library(caret)
library(PRROC)
library(pheatmap)

### === Panel A: ROC Curve === ###
roc_plot <- ggroc(roc_rf, colour = "darkgreen", size = 1) +
   geom_abline(linetype = "dashed", color = "grey50") +
   labs(title = paste0("ROC Curve (AUC = ", round(auc_rf, 3), ")"),
        x = "False Positive Rate", y = "True Positive Rate") +
   theme_minimal(base_size = 12)

### === Panel B: Publication-ready Sensitivity & Specificity Heatmap === ###
# Compute metrics for a range of thresholds
thresholds <- seq(0, 1, by = 0.01)
sens <- spec <- numeric(length(thresholds))

for(i in seq_along(thresholds)) {
   pred <- factor(ifelse(prob_rf >= thresholds[i], "Advanced", "Early"),
                  levels = c("Early", "Advanced"))
   cm <- confusionMatrix(pred, expr_df$fibrosis_bin)
   sens[i] <- cm$byClass["Sensitivity"]
   spec[i] <- cm$byClass["Specificity"]
}

# Build data frame for heatmap
heat_df <- data.frame(
   Threshold = thresholds,
   Sensitivity = sens,
   Specificity = spec
) %>%
   pivot_longer(cols = c("Sensitivity", "Specificity"),
                names_to = "Metric", values_to = "Value")

# Get values at optimal threshold
sens_opt <- heat_df$Value[heat_df$Threshold == round(opt_threshold, 2) & heat_df$Metric == "Sensitivity"]
spec_opt <- heat_df$Value[heat_df$Threshold == round(opt_threshold, 2) & heat_df$Metric == "Specificity"]





### === Panel C: Probability Distributions === ###
prob_df <- data.frame(Probability = prob_rf,
                      Class = expr_df$fibrosis_bin)
prob_plot <- ggplot(prob_df, aes(x = Probability, fill = Class)) +
   geom_density(alpha = 0.5) +
   geom_vline(xintercept = opt_threshold, linetype = "dashed", color = "red") +
   labs(title = "Predicted Probability Distributions",
        x = "RF Predicted Probability (Advanced)", y = "Density") +
   theme_minimal(base_size = 12)





### === Panel D: Heatmap of Misclassified Samples === ###
# ---- Misclassified heatmap (fixed) ----
library(pheatmap)

# 1) pick top N important genes that actually exist in expr_df
imp <- importance(rf_cv$finalModel)
imp_df <- data.frame(Gene = rownames(imp), Importance = imp[,1])
imp_df <- imp_df %>% arrange(desc(Importance)) %>% slice(1:10)

top_n <- 20
top_genes <- head(imp_df$Gene, top_n)
top_genes <- intersect(top_genes, colnames(expr_df))  # ensure present

# 2) index misclassified samples
mis_idx <- which(pred_rf != expr_df$fibrosis_bin)

if (length(mis_idx) > 0 && length(top_genes) > 1) {
   misclassified <- expr_df[mis_idx, , drop = FALSE]
   
   # 3) build annotation with matching rownames = sample IDs (i.e., columns of heatmap)
   ann <- data.frame(
      True = expr_df$fibrosis_bin[mis_idx],
      Pred = pred_rf[mis_idx],
      Prob_Advanced = round(prob_rf[mis_idx], 3)
   )
   rownames(ann) <- rownames(misclassified)
   
   # 4) matrix for heatmap: genes x samples (numeric only)
   mat <- t(as.matrix(misclassified[, top_genes, drop = FALSE]))
   
   # 5) (optional) reorder samples by predicted class, then by probability
   ord <- order(ann$Pred, -ann$Prob_Advanced)
   mat <- mat[, ord, drop = FALSE]
   ann <- ann[ord, , drop = FALSE]
   
   pheat <- pheatmap(
      mat,
      annotation_col = ann,
      main = "Misclassified Samples (Top Important Genes)",
      cluster_cols = TRUE,
      cluster_rows = TRUE,
      show_colnames = FALSE,
      silent = TRUE
   )
} else {
   message("No misclassified samples or too few genes to plot.")
}















### === Panel E: Sensitivity–Specificity Curves === ###
thresholds <- seq(0, 1, by = 0.01)

metrics <- data.frame(
   Threshold = thresholds,
   Sensitivity = sapply(thresholds, function(t) {
      tp <- sum(prob_rf >= t & expr_df$fibrosis_bin == "Advanced")
      fn <- sum(prob_rf <  t & expr_df$fibrosis_bin == "Advanced")
      tp / (tp + fn)
   }),
   Specificity = sapply(thresholds, function(t) {
      tn <- sum(prob_rf <  t & expr_df$fibrosis_bin == "Early")
      fp <- sum(prob_rf >= t & expr_df$fibrosis_bin == "Early")
      tn / (tn + fp)
   })
)

# Find crossing point (optimal threshold)
metrics$Diff <- abs(metrics$Sensitivity - metrics$Specificity)
optimal <- metrics[which.min(metrics$Diff), ]

sens_spec_plot <- ggplot(metrics, aes(x = Threshold)) +
   geom_line(aes(y = Sensitivity, color = "Sensitivity"), size = 0.6) +
   geom_line(aes(y = Specificity, color = "Specificity"), size = 0.6) +
   geom_vline(xintercept = optimal$Threshold, linetype = "dashed", color = "black") +
   geom_point(aes(x = optimal$Threshold, y = optimal$Sensitivity), 
              color = "red", size = 2.2) +
   labs(y = "Rate", x = "Threshold") +
   scale_color_manual(values = c("Sensitivity" = "#1f77b4",   # blue-ish
                                 "Specificity" = "#2ca02c")) + # green-ish
   theme_minimal(base_size = 9) +
   theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5)
   )

print(sens_spec_plot)

# Save smaller for 3x3 figure
ggsave("Sensitivity_Specificity_Curves.pdf", plot = sens_spec_plot,
       width = 2.5, height = 2.5, dpi = 300)







library(ggplot2)
library(pROC)

# Create data frame from roc object
roc_df <- data.frame(
   FPR = 1 - roc_rf$specificities,
   TPR = roc_rf$sensitivities
)

# AUC rounded
auc_label <- paste0("AUC = ", round(auc_rf, 3))

# === PUBLISHABLE ROC CURVE === #
p_roc_pub <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
   geom_line(color = "#6A1B9A", size = 1.2) +               # deep purple
   geom_abline(slope = 1, intercept = 0, 
               linetype = "dashed", color = "grey50", size = 0.8) +
   annotate("text", x = 0.65, y = 0.1, label = auc_label, 
            size = 4.5, hjust = 0, color = "black") +
   labs(
      title = "Proteomics dataset",
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
   ) +
   coord_fixed() +  # square aspect ratio
   theme_minimal(base_size = 14) +
   theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(size = 0.3),
      axis.title = element_text(face = "bold"),
      text = element_text(color = "black")
   )

print(p_roc_pub)

# === SAVE HIGH-RES IMAGES === #
ggsave("ROC_Govaere_Proteomics_RF.pdf", 
       plot = p_roc_pub, width = 6, height = 6)

df_57 <- data.frame(
   FPR = 1 - roc_rf$specificities,
   TPR = roc_rf$sensitivities,
   Model = "57 genes"
)















### === LOAD TRAINING DATA AND RF MODEL === ###
# Load patient metadata
template <- read.csv("data/metadata.csv")[-5, -c(1,3)]

# Load z-score normalised gene expression matrix
corrected_ucamsanyal <- read.csv("data/z_scores.csv", check.names = FALSE)

# Load 57 selected gene mapping table (external names + ENSG IDs)
selected_57 <- read.csv("data/top_15_BMs_from_ELBOW.csv")[,-1]
selected_57 = selected_57[ is.element(selected_57$ensembl_gene_id, corrected_ucamsanyal[,1]), ]

# Preprocess training matrix
rownames(corrected_ucamsanyal) <- corrected_ucamsanyal[,1]
corrected_ucamsanyal <- corrected_ucamsanyal[,-1]
corrected_ucamsanyal <- corrected_ucamsanyal[selected_57$ensembl_gene_id, ]
template <- template[match(colnames(corrected_ucamsanyal), template$Sample.name), ]
stopifnot(all(template$Sample.name == colnames(corrected_ucamsanyal)))
expr_57 <- as.data.frame(t(corrected_ucamsanyal))

# Define binary label
expr_57$fibrosis_bin <- factor(ifelse(as.numeric(as.character(template$Fibrosis)) >= 3, "Advanced", "Early"),
                               levels = c("Early", "Advanced"))

# Train RF model with 5-fold CV
set.seed(123)
ctrl <- trainControl(method = "cv",
                     number = 5,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary,
                     savePredictions = "final")
rf_cv <- train(fibrosis_bin ~ ., data = expr_57,
               method = "rf",
               trControl = ctrl,
               metric = "ROC")



### === LOAD GOVAERE DATASET === ###
govaere <- read_excel("data/Plasma_Proteomics_dataset/42255_2023_775_MOESM3_ESM.xlsx")[-c(1,2,3,4), -c(1:6,8)]
govaere <- govaere %>%
   distinct(.[[7]], .keep_all = TRUE) %>%
   rename(sample_id = 1) %>%
   filter(!duplicated(sample_id), !is.na(sample_id)) %>%
   mutate(across(-c(sample_id,2), as.double))
govaere <- as.data.frame(govaere[ ,-c(dim(govaere)[2])])
rownames(govaere) <- govaere$sample_id
govaere <- govaere[ , -c(1,2)]

# Metadata for fibrosis
metadata <- read_excel("data/Plasma_Proteomics_dataset/42255_2023_775_MOESM3_ESM.xlsx")[ -c(1,2) , -c(1:9)]
metadata <- metadata[, colnames(govaere)]
metadata <- as.double(metadata[1,])

### === PREPROCESS GOVAERE DATA === ###
# 1. Identify existing/missing genes
selected_genes <- selected_57$external_gene_name
existing_genes <- intersect(selected_genes, rownames(govaere))
missing_genes <- setdiff(selected_genes, rownames(govaere))

# 2. Subset existing genes
govaere_sub <- govaere[existing_genes, , drop = FALSE]

# 3. Add missing genes as rows of 0
if(length(missing_genes) > 0){
   missing_mat <- matrix(0, nrow = length(missing_genes), ncol = ncol(govaere_sub),
                         dimnames = list(missing_genes, colnames(govaere_sub)))
   govaere_sub <- rbind(govaere_sub, missing_mat)
}

# 4. Z-score normalize, avoiding NaNs for missing genes
govaere_sub_z <- t(apply(govaere_sub, 1, function(x) {
   if (sd(x) == 0) {
      rep(0, length(x))  # keep as 0 if all values identical
   } else {
      scale(x)
   }
}))


# 5. Transpose to samples-as-rows
expr_df <- as.data.frame(t(govaere_sub_z))

# 6. Add fibrosis_bin factor
expr_df$fibrosis_bin <- factor(ifelse(metadata == 1, "Advanced", "Early"),
                               levels = c("Early", "Advanced"))

# 7. Map column names to ENSG IDs
name2ensg <- setNames(selected_57$ensembl_gene_id, selected_57$external_gene_name)
expr_genes_no_fib <- setdiff(colnames(expr_df), "fibrosis_bin")
expr_df <- expr_df %>%
   rename_with(.cols = all_of(expr_genes_no_fib), .fn = ~ name2ensg[.])

# 8. Reorder columns to match training set
train_genes <- setdiff(colnames(expr_57), "fibrosis_bin")
expr_df <- expr_df[, c(train_genes, "fibrosis_bin")]

### === PREDICTIONS AND EVALUATION === ###
pred_rf <- predict(rf_cv$finalModel, newdata = expr_df)
prob_rf <- predict(rf_cv$finalModel, newdata = expr_df, type = "prob")[, "Advanced"]

# Confusion matrix
conf_rf <- confusionMatrix(pred_rf, expr_df$fibrosis_bin)

# ROC curve and AUC
roc_rf <- roc(expr_df$fibrosis_bin, prob_rf)
auc_rf <- auc(roc_rf)

cat("Govaere dataset RF AUC:", auc_rf, "\n")
print(conf_rf)

plot(roc_rf, print.auc = TRUE, col = "darkgreen", lwd = 2,
     main = "Govaere Dataset ROC: 57-Gene Classifier")











### === OPTIMAL THRESHOLD FOR HIGHER SENSITIVITY === ###
# Get ROC with correct positive/negative classes
roc_rf <- roc(expr_df$fibrosis_bin, prob_rf, levels = c("Early", "Advanced"))

# Option 1: Youden's index (balances sensitivity & specificity)
opt_coords <- coords(roc_rf, x = "best", best.method = "youden")
opt_threshold <- opt_coords["threshold"]
opt_threshold <- opt_threshold$threshold  # extract as numeric

# Option 2: Force high sensitivity (≥ 0.9) — uncomment to use this instead
# opt_coords <- coords(roc_rf, x = 0.9, input = "sensitivity")
# opt_threshold <- opt_coords["threshold"]

cat("Optimal threshold for higher sensitivity:", opt_threshold, "\n")

# Predict using the new threshold
pred_rf_opt <- factor(ifelse(prob_rf >= opt_threshold, "Advanced", "Early"),
                      levels = c("Early", "Advanced"))

# Confusion matrix with new threshold
conf_rf_opt <- confusionMatrix(pred_rf_opt, expr_df$fibrosis_bin)
print(conf_rf_opt)

# Plot ROC with optimal threshold point
plot(roc_rf, print.auc = TRUE, col = "darkgreen", lwd = 2,
     main = paste0("Govaere Dataset ROC: 57-Gene Classifier\nOptimal threshold = ", round(opt_threshold, 3)))
points(opt_coords["threshold"], opt_coords["sensitivity"], col = "red", pch = 19)







#####
#####
#####
library(ggplot2)
library(cowplot)
library(pROC)
library(caret)
library(PRROC)
library(pheatmap)

### === Panel A: ROC Curve === ###
roc_plot <- ggroc(roc_rf, colour = "darkgreen", size = 1) +
   geom_abline(linetype = "dashed", color = "grey50") +
   labs(title = paste0("ROC Curve (AUC = ", round(auc_rf, 3), ")"),
        x = "False Positive Rate", y = "True Positive Rate") +
   theme_minimal(base_size = 12)

### === Panel B: Publication-ready Sensitivity & Specificity Heatmap === ###
# Compute metrics for a range of thresholds
thresholds <- seq(0, 1, by = 0.01)
sens <- spec <- numeric(length(thresholds))

for(i in seq_along(thresholds)) {
   pred <- factor(ifelse(prob_rf >= thresholds[i], "Advanced", "Early"),
                  levels = c("Early", "Advanced"))
   cm <- confusionMatrix(pred, expr_df$fibrosis_bin)
   sens[i] <- cm$byClass["Sensitivity"]
   spec[i] <- cm$byClass["Specificity"]
}

# Build data frame for heatmap
heat_df <- data.frame(
   Threshold = thresholds,
   Sensitivity = sens,
   Specificity = spec
) %>%
   pivot_longer(cols = c("Sensitivity", "Specificity"),
                names_to = "Metric", values_to = "Value")

# Get values at optimal threshold
sens_opt <- heat_df$Value[heat_df$Threshold == round(opt_threshold, 2) & heat_df$Metric == "Sensitivity"]
spec_opt <- heat_df$Value[heat_df$Threshold == round(opt_threshold, 2) & heat_df$Metric == "Specificity"]





### === Panel C: Probability Distributions === ###
prob_df <- data.frame(Probability = prob_rf,
                      Class = expr_df$fibrosis_bin)
prob_plot <- ggplot(prob_df, aes(x = Probability, fill = Class)) +
   geom_density(alpha = 0.5) +
   geom_vline(xintercept = opt_threshold, linetype = "dashed", color = "red") +
   labs(title = "Predicted Probability Distributions",
        x = "RF Predicted Probability (Advanced)", y = "Density") +
   theme_minimal(base_size = 12)





### === Panel D: Heatmap of Misclassified Samples === ###
# ---- Misclassified heatmap (fixed) ----
library(pheatmap)

# 1) pick top N important genes that actually exist in expr_df
imp <- importance(rf_cv$finalModel)
imp_df <- data.frame(Gene = rownames(imp), Importance = imp[,1])
imp_df <- imp_df %>% arrange(desc(Importance)) %>% slice(1:10)

top_n <- 20
top_genes <- head(imp_df$Gene, top_n)
top_genes <- intersect(top_genes, colnames(expr_df))  # ensure present

# 2) index misclassified samples
mis_idx <- which(pred_rf != expr_df$fibrosis_bin)

if (length(mis_idx) > 0 && length(top_genes) > 1) {
   misclassified <- expr_df[mis_idx, , drop = FALSE]
   
   # 3) build annotation with matching rownames = sample IDs (i.e., columns of heatmap)
   ann <- data.frame(
      True = expr_df$fibrosis_bin[mis_idx],
      Pred = pred_rf[mis_idx],
      Prob_Advanced = round(prob_rf[mis_idx], 3)
   )
   rownames(ann) <- rownames(misclassified)
   
   # 4) matrix for heatmap: genes x samples (numeric only)
   mat <- t(as.matrix(misclassified[, top_genes, drop = FALSE]))
   
   # 5) (optional) reorder samples by predicted class, then by probability
   ord <- order(ann$Pred, -ann$Prob_Advanced)
   mat <- mat[, ord, drop = FALSE]
   ann <- ann[ord, , drop = FALSE]
   
   pheat <- pheatmap(
      mat,
      annotation_col = ann,
      main = "Misclassified Samples (Top Important Genes)",
      cluster_cols = TRUE,
      cluster_rows = TRUE,
      show_colnames = FALSE,
      silent = TRUE
   )
} else {
   message("No misclassified samples or too few genes to plot.")
}















### === Panel E: Sensitivity–Specificity Curves === ###
thresholds <- seq(0, 1, by = 0.01)

metrics <- data.frame(
   Threshold = thresholds,
   Sensitivity = sapply(thresholds, function(t) {
      tp <- sum(prob_rf >= t & expr_df$fibrosis_bin == "Advanced")
      fn <- sum(prob_rf <  t & expr_df$fibrosis_bin == "Advanced")
      tp / (tp + fn)
   }),
   Specificity = sapply(thresholds, function(t) {
      tn <- sum(prob_rf <  t & expr_df$fibrosis_bin == "Early")
      fp <- sum(prob_rf >= t & expr_df$fibrosis_bin == "Early")
      tn / (tn + fp)
   })
)

# Find crossing point (optimal threshold)
metrics$Diff <- abs(metrics$Sensitivity - metrics$Specificity)
optimal <- metrics[which.min(metrics$Diff), ]

sens_spec_plot <- ggplot(metrics, aes(x = Threshold)) +
   geom_line(aes(y = Sensitivity, color = "Sensitivity"), size = 0.6) +
   geom_line(aes(y = Specificity, color = "Specificity"), size = 0.6) +
   geom_vline(xintercept = optimal$Threshold, linetype = "dashed", color = "black") +
   geom_point(aes(x = optimal$Threshold, y = optimal$Sensitivity), 
              color = "red", size = 2.2) +
   labs(y = "Rate", x = "Threshold") +
   scale_color_manual(values = c("Sensitivity" = "#1f77b4",   # blue-ish
                                 "Specificity" = "#2ca02c")) + # green-ish
   theme_minimal(base_size = 9) +
   theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5)
   )

print(sens_spec_plot)

# Save smaller for 3x3 figure
ggsave("Sensitivity_Specificity_Curves_15.pdf", plot = sens_spec_plot,
       width = 2.5, height = 2.5, dpi = 300)







library(ggplot2)
library(pROC)

# Create data frame from roc object
roc_df <- data.frame(
   FPR = 1 - roc_rf$specificities,
   TPR = roc_rf$sensitivities
)

# AUC rounded
auc_label <- paste0("AUC = ", round(auc_rf, 3))

# === PUBLISHABLE ROC CURVE === #
p_roc_pub <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
   geom_line(color = "#6A1B9A", size = 1.2) +               # deep purple
   geom_abline(slope = 1, intercept = 0, 
               linetype = "dashed", color = "grey50", size = 0.8) +
   annotate("text", x = 0.65, y = 0.1, label = auc_label, 
            size = 4.5, hjust = 0, color = "black") +
   labs(
      title = "Proteomics dataset",
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
   ) +
   coord_fixed() +  # square aspect ratio
   theme_minimal(base_size = 14) +
   theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(size = 0.3),
      axis.title = element_text(face = "bold"),
      text = element_text(color = "black")
   )

print(p_roc_pub)

# === SAVE HIGH-RES IMAGES === #
ggsave("ROC_Govaere_Proteomics_RF_15.pdf", 
       plot = p_roc_pub, width = 6, height = 6)


df_15 <- data.frame(
   FPR = 1 - roc_rf$specificities,
   TPR = roc_rf$sensitivities,
   Model = "15 genes"
)




# Plot both for the 57 and 15 in one plot
df_all <- rbind(df_57, df_15)


### --- PLOTTING --- ###
combined_sens_spec_plot <- ggplot(df_all, aes(x = FPR, y = TPR, color = Model)) +
   geom_line(size = 1.4) +
   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
   scale_color_manual(values = c(
      "57 genes"  = "purple4",   # deep purple
      "15 genes"  = "mediumpurple"    # lighter purple
   )) +
   labs(title = "ROC curves: RF models (57 vs 15 genes)",
        x = "False Positive Rate",
        y = "True Positive Rate") +
   theme_minimal(base_size = 14) +
   theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_blank()
   )

ggsave("ROC_Govaere_Proteomics_RF_combined.pdf", 
       plot = combined_sens_spec_plot, width = 6, height = 6)



