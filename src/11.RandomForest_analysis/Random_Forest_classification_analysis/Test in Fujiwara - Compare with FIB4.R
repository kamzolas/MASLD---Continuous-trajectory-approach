### === SETUP === ###
library(randomForest)
library(caret)
library(ggplot2)
library(pROC)
library(reshape2)
library(readxl)
library(tidyverse)

### === LOAD DATA === ###
# Load patient metadata (remove row 5 and unnecessary columns)
template <- read.csv("data/metadata.csv")[-5, -c(1,3)]

# Load z-score normalised gene expression matrix
corrected_ucamsanyal <- read.csv(
   "data/z_scores.csv",
   check.names = FALSE
)

# Load 57 selected gene ENSEMBL IDs
selected_57 <- read.csv("data/57_BMs.csv")

### === PREPROCESS EXPRESSION DATA === ###
# Set rownames and remove gene ID column
rownames(corrected_ucamsanyal) <- corrected_ucamsanyal[,1]
corrected_ucamsanyal <- corrected_ucamsanyal[,-1]

# Subset expression matrix to 24 selected genes
selected_57 = selected_57[is.element(selected_57$ensembl_gene_id, rownames(corrected_ucamsanyal)),]
corrected_ucamsanyal <- corrected_ucamsanyal[selected_57$ensembl_gene_id, ]

# Ensure sample order matches metadata
template <- template[match(colnames(corrected_ucamsanyal), template$Sample.name), ]
stopifnot(all(template$Sample.name == colnames(corrected_ucamsanyal)))

# Transpose expression matrix: patients as rows
expr_24 <- as.data.frame(t(corrected_ucamsanyal))

### === DEFINE BINARY LABEL: EARLY (F0–F2) VS ADVANCED (F3–F4) === ###
expr_24$fibrosis_bin <- ifelse(as.numeric(as.character(template$Fibrosis)) >= 3, "Advanced", "Early")
expr_24$fibrosis_bin <- factor(expr_24$fibrosis_bin, levels = c("Early", "Advanced"))

### === RANDOM FOREST CLASSIFIER WITH CROSS-VALIDATION === ###
set.seed(123)
ctrl <- trainControl(method = "cv",
                     number = 5,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary,
                     savePredictions = "final")

# Train RF model using caret with CV
rf_cv <- train(fibrosis_bin ~ ., data = expr_24,
               method = "rf",
               trControl = ctrl,
               metric = "ROC")

# Print model performance summary
print(rf_cv)

### === EVALUATE CV PERFORMANCE === ###
# Confusion matrix from CV predictions
cv_conf <- confusionMatrix(rf_cv$pred$pred, rf_cv$pred$obs)
print(cv_conf)

### === ROC CURVE FROM CV === ###
roc_cv <- roc(rf_cv$pred$obs, rf_cv$pred$Advanced)
auc_cv <- auc(roc_cv)
cat("AUC from CV =", auc_cv, "\n")

plot(roc_cv, print.auc = TRUE, col = "darkblue", lwd = 2,
     main = "5-Fold CV ROC: 24-Gene Classifier (F3–F4 vs F0–F2)")

### === VARIABLE IMPORTANCE (from final RF model) === ###
varImpPlot(rf_cv$finalModel, main = "Gene Importance: Early vs Advanced Fibrosis")

### === EXPORT PREDICTIONS AND IMPORTANCE === ###
# Save CV predictions
#write.csv(rf_cv$pred, "RF_24gene_CV_Predictions.csv", row.names = FALSE)

# Export variable importance
importance_df <- as.data.frame(importance(rf_cv$finalModel))
importance_df$Gene <- rownames(importance_df)
importance_df <- importance_df[order(importance_df$MeanDecreaseGini, decreasing = TRUE), ]
#write.csv(importance_df, "GeneImportance_RF_24gene_CV.csv", row.names = FALSE)

# Optional: Heatmap-style confusion matrix
cm_table <- as.table(cv_conf$table)
cm_df <- as.data.frame(cm_table)
ggplot(cm_df, aes(Prediction, Reference, fill = Freq)) +
   geom_tile(color = "white") +
   geom_text(aes(label = Freq), vjust = 0.5, size = 5) +
   scale_fill_gradient(low = "white", high = "#d7191c") +
   labs(title = "Confusion Matrix (CV): 24-Gene Classifier", fill = "Count") +
   theme_minimal()











#External dataset for validation - FUJIWARA
fujiwara=read.csv("data/Fujiwara_dataset/raw_counts_fujiwara.csv", sep = ",")
rownames(fujiwara) = fujiwara$ENSID
fujiwara = fujiwara[ ,-c(1)]
fujiwara = as.matrix(fujiwara)
dim(fujiwara)


# Read sheets
df_bx1 <- read_xlsx("data/Fujiwara_dataset/Both biopsies - refined dataset.xlsx", sheet = 1)
df_bx2 <- read_xlsx("data/Fujiwara_dataset/Both biopsies - refined dataset.xlsx", sheet = 2)

# Select relevant columns from first biopsy
df_bx1_sub <- df_bx1 %>%
   select(gct.name, Patient, Biopsy, Age, Platelet, ALB, AST, ALT, BMI, Diabetes, Histology.fibrosis)

# Select relevant columns from second biopsy, rename *_2 columns to base names
df_bx2_sub <- df_bx2 %>%
   select(gct.name = sample.names, Patient, Biopsy, Age, PLT_2, ALB_2, AST_2, ALT_2, BMI_2, Diabetes, Histology.fibrosis_2) %>%
   rename(
      Platelet = PLT_2,
      ALB = ALB_2,
      AST = AST_2,
      ALT = ALT_2,
      BMI = BMI_2,
      Histology.fibrosis = Histology.fibrosis_2
   )

# Add biopsy number for clarity
df_bx1_sub$biopsy_num <- 1
df_bx2_sub$biopsy_num <- 2

# Convert AST and other columns to numeric in both dfs before combining
cols_to_convert <- c("AST", "ALT", "Platelet", "ALB", "BMI", "Age", "Diabetes", "Histology.fibrosis")

# Helper function to convert columns to numeric safely
convert_to_numeric <- function(df, cols) {
   for (col in cols) {
      if (col %in% colnames(df)) {
         df[[col]] <- as.numeric(as.character(df[[col]]))
      }
   }
   return(df)
}

# Apply to both data frames
df_bx1_sub <- convert_to_numeric(df_bx1_sub, cols_to_convert)
df_bx2_sub <- convert_to_numeric(df_bx2_sub, cols_to_convert)

# Now bind rows
df_combined <- bind_rows(df_bx1_sub, df_bx2_sub)

#write.csv(df_combined, "Fujiwara_Metadata_bothbiopsies_refined_and_combined.csv")





# Assuming your 57 genes list:
selected_genes <- selected_57$ensembl_gene_id

# Subset expression matrix to 57 genes
fujiwara_57 <- fujiwara[selected_genes, , drop=FALSE]

# Get samples common to expression and metadata
common_samples <- intersect(colnames(fujiwara_57), df_combined$gct.name)

# Subset expression and metadata
fujiwara_57 <- fujiwara_57[, common_samples, drop=FALSE]
df_combined_sub <- df_combined %>% filter(gct.name %in% common_samples)




# z-score normalize gene expression (genes in rows, samples in columns)
fujiwara_57_z <- t(scale(t(fujiwara_57)))

# Convert to data.frame with samples as rows (for prediction)
expr_df <- as.data.frame(t(fujiwara_57_z))

# Add fibrosis label factor consistent with training
expr_df$fibrosis_bin <- ifelse(df_combined_sub$Histology.fibrosis >= 3, "Advanced", "Early")
expr_df$fibrosis_bin <- factor(expr_df$fibrosis_bin, levels = c("Early", "Advanced"))




# Predict fibrosis class
pred_fujiwara <- predict(rf_cv$finalModel, newdata = expr_df)

# Predict probabilities for ROC etc.
prob_fujiwara <- predict(rf_cv$finalModel, newdata = expr_df, type = "prob")[, "Advanced"]

# Evaluate performance
library(caret)
conf_fujiwara <- confusionMatrix(pred_fujiwara, expr_df$fibrosis_bin)
print(conf_fujiwara)





library(pROC)

roc_fujiwara <- roc(expr_df$fibrosis_bin, prob_fujiwara)
auc_fujiwara <- auc(roc_fujiwara)
cat("FUJIWARA AUC =", auc_fujiwara, "\n")

plot(roc_fujiwara, print.auc = TRUE, col = "darkgreen", lwd = 2,
     main = "FUJIWARA Dataset ROC: 57-Gene Classifier (F3–F4 vs F0–F2)")





library(dplyr)

AST_ULN <- 40

df_combined_sub <- df_combined_sub %>%
   mutate(
      Age = as.numeric(Age),
      AST = as.numeric(AST),
      ALT = as.numeric(ALT),
      Platelet = as.numeric(Platelet),
      BMI = as.numeric(BMI),
      Diabetes_num = as.numeric(Diabetes),
      ALB = as.numeric(ALB),
      
      # Avoid zeros or negative in denominators (optional)
      AST = ifelse(is.na(AST) | AST <= 0, NA, AST),
      ALT = ifelse(is.na(ALT) | ALT <= 0, NA, ALT),
      Platelet = ifelse(is.na(Platelet) | Platelet <= 0, NA, Platelet),
      Age = ifelse(is.na(Age) | Age <= 0, NA, Age),
      
      FIB4 = (Age * AST) / (Platelet * sqrt(ALT)),
      APRI = ((AST / AST_ULN) / Platelet) * 100,
      NFS = -1.675 + 0.037 * Age + 0.094 * BMI + 1.13 * Diabetes_num +
         0.99 * (AST / ALT) - 0.013 * Platelet - 0.66 * ALB
   )


write.csv(df_combined_sub, "Fujiwara_Metadata_bothbiopsies_refined_and_combined_57.csv")










library(dplyr)
library(caret)
library(pROC)
library(tidyr)
library(ggplot2)

# --- Step 0: Define AST ULN ---
AST_ULN <- 40

# --- Step 1: Calculate NIT scores (if not done yet) ---
df_combined_sub <- df_combined_sub %>%
   mutate(
      Age = as.numeric(Age),
      AST = as.numeric(AST),
      ALT = as.numeric(ALT),
      Platelet = as.numeric(Platelet),
      BMI = as.numeric(BMI),
      Diabetes_num = as.numeric(Diabetes),
      ALB = ifelse("ALB" %in% colnames(df_combined_sub), as.numeric(ALB), NA_real_),
      
      AST = ifelse(is.na(AST) | AST <= 0, NA, AST),
      ALT = ifelse(is.na(ALT) | ALT <= 0, NA, ALT),
      Platelet = ifelse(is.na(Platelet) | Platelet <= 0, NA, Platelet),
      Age = ifelse(is.na(Age) | Age <= 0, NA, Age),
      
      FIB4 = (Age * AST) / (Platelet * sqrt(ALT)),
      APRI = ((AST / AST_ULN) / Platelet) * 100,
      NFS = ifelse(!is.na(ALB),
                   -1.675 + 0.037 * Age + 0.094 * BMI + 1.13 * Diabetes_num +
                      0.99 * (AST / ALT) - 0.013 * Platelet - 0.66 * ALB,
                   NA_real_)
   )

# --- Step 2: Create fibrosis_bin factor ---
df_combined_sub$fibrosis_bin <- ifelse(df_combined_sub$Histology.fibrosis >= 3, "Advanced", "Early")
df_combined_sub$fibrosis_bin <- factor(df_combined_sub$fibrosis_bin, levels = c("Early", "Advanced"))

# --- Step 3: RF predictions ---
# Ensure expr_df contains gene expression data for RF model predictions with fibrosis_bin
pred_rf <- predict(rf_cv$finalModel, newdata = expr_df)
prob_rf <- predict(rf_cv$finalModel, newdata = expr_df, type = "prob")[, "Advanced"]

# --- Step 4: NITs binary predictions with cutoffs ---
pred_fib4 <- factor(ifelse(df_combined_sub$FIB4 >= 1.45, "Advanced", "Early"), levels = c("Early", "Advanced"))
pred_apri <- factor(ifelse(df_combined_sub$APRI >= 0.7, "Advanced", "Early"), levels = c("Early", "Advanced"))
pred_nfs  <- factor(ifelse(df_combined_sub$NFS  >= -1.455, "Advanced", "Early"), levels = c("Early", "Advanced"))

# --- Step 5: Function to compute metrics (with AUC) ---
get_metrics <- function(pred, true_label, prob = NULL) {
   cm <- confusionMatrix(pred, true_label, positive = "Advanced")
   auc_val <- if (!is.null(prob) && !all(is.na(prob))) {
      tryCatch({
         roc_obj <- roc(true_label, prob, levels = c("Early", "Advanced"))
         auc(roc_obj)
      }, error = function(e) NA)
   } else {
      NA
   }
   data.frame(
      Accuracy = cm$overall["Accuracy"],
      Sensitivity = cm$byClass["Sensitivity"],
      Specificity = cm$byClass["Specificity"],
      AUC = auc_val
   )
}

# --- Step 6: Calculate metrics for all methods ---
metrics_rf <- get_metrics(pred_rf, expr_df$fibrosis_bin, prob = prob_rf)
metrics_fib4 <- get_metrics(pred_fib4, df_combined_sub$fibrosis_bin, prob = df_combined_sub$FIB4)
metrics_apri <- get_metrics(pred_apri, df_combined_sub$fibrosis_bin, prob = df_combined_sub$APRI)
metrics_nfs  <- get_metrics(pred_nfs,  df_combined_sub$fibrosis_bin, prob = df_combined_sub$NFS)

# --- Step 7: Combine results ---
metrics_df <- rbind(
   RF_Model = metrics_rf,
   FIB4 = metrics_fib4,
   APRI = metrics_apri,
   NFS = metrics_nfs
)

metrics_df$Method <- rownames(metrics_df)
metrics_df <- metrics_df[, c("Method", "Accuracy", "Sensitivity", "Specificity", "AUC")]

print(metrics_df)

# --- Step 8: Plot metrics comparison ---
# Convert AUC column to numeric (from 'auc' S3 class)
metrics_df$AUC <- as.numeric(metrics_df$AUC)

# Pivot to long format
metrics_long_57 <- metrics_df %>%
   pivot_longer(cols = Accuracy:AUC, names_to = "Metric", values_to = "Value")


ggplot(metrics_long_57, aes(x = Metric, y = Value, fill = Method)) +
   geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
   scale_fill_manual(values = c("purple", "gray40", "orange", "forestgreen")) +
   scale_y_continuous(breaks = seq(0.4, .9, 0.05)) +
   coord_cartesian(ylim = c(0.4, .9)) +   # zoom without dropping data
   labs(title = "Performance Comparison: RF Model vs NITs on Fujiwara Dataset",
        y = "Metric Value", x = "") +
   theme_minimal(base_size = 14) +
   theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, face = "bold")
   )

# --- Save the plot with white background ---
ggsave("Fujiwara_RF_vs_NITs_Performance_57.pdf", width = 8, height = 6, dpi = 300, bg = "white")

















#Same for the top 15 biomarkers


### === LOAD DATA === ###
# Load patient metadata (remove row 5 and unnecessary columns)
template <- read.csv("data/metadata.csv")[-5, -c(1,3)]

# Load z-score normalised gene expression matrix
corrected_ucamsanyal <- read.csv(
   "data/z_scores.csv",
   check.names = FALSE
)

# Load 15 selected gene ENSEMBL IDs
selected_15 <- read.csv("data/top_15_BMs_from_ELBOW.csv")

### === PREPROCESS EXPRESSION DATA === ###
# Set rownames and remove gene ID column
rownames(corrected_ucamsanyal) <- corrected_ucamsanyal[,1]
corrected_ucamsanyal <- corrected_ucamsanyal[,-1]

# Subset expression matrix to 15 selected genes
selected_15 = selected_15[is.element(selected_15$ensembl_gene_id, rownames(corrected_ucamsanyal)),]
corrected_ucamsanyal <- corrected_ucamsanyal[selected_15$ensembl_gene_id, ]

# Ensure sample order matches metadata
template <- template[match(colnames(corrected_ucamsanyal), template$Sample.name), ]
stopifnot(all(template$Sample.name == colnames(corrected_ucamsanyal)))

# Transpose expression matrix: patients as rows
expr_15 <- as.data.frame(t(corrected_ucamsanyal))

### === DEFINE BINARY LABEL: EARLY (F0–F2) VS ADVANCED (F3–F4) === ###
expr_15$fibrosis_bin <- ifelse(as.numeric(as.character(template$Fibrosis)) >= 3, "Advanced", "Early")
expr_15$fibrosis_bin <- factor(expr_15$fibrosis_bin, levels = c("Early", "Advanced"))

### === RANDOM FOREST CLASSIFIER WITH CROSS-VALIDATION === ###
set.seed(123)
ctrl <- trainControl(method = "cv",
                     number = 5,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary,
                     savePredictions = "final")

# Train RF model using caret with CV
rf_cv <- train(fibrosis_bin ~ ., data = expr_15,
               method = "rf",
               trControl = ctrl,
               metric = "ROC")

# Print model performance summary
print(rf_cv)

### === EVALUATE CV PERFORMANCE === ###
# Confusion matrix from CV predictions
cv_conf <- confusionMatrix(rf_cv$pred$pred, rf_cv$pred$obs)
print(cv_conf)

### === ROC CURVE FROM CV === ###
roc_cv <- roc(rf_cv$pred$obs, rf_cv$pred$Advanced)
auc_cv <- auc(roc_cv)
cat("AUC from CV =", auc_cv, "\n")

plot(roc_cv, print.auc = TRUE, col = "darkblue", lwd = 2,
     main = "5-Fold CV ROC: 15-Gene Classifier (F3–F4 vs F0–F2)")

### === VARIABLE IMPORTANCE (from final RF model) === ###
varImpPlot(rf_cv$finalModel, main = "Gene Importance: Early vs Advanced Fibrosis")

### === EXPORT PREDICTIONS AND IMPORTANCE === ###
# Save CV predictions
#write.csv(rf_cv$pred, "RF_15gene_CV_Predictions.csv", row.names = FALSE)

# Export variable importance
importance_df <- as.data.frame(importance(rf_cv$finalModel))
importance_df$Gene <- rownames(importance_df)
importance_df <- importance_df[order(importance_df$MeanDecreaseGini, decreasing = TRUE), ]
#write.csv(importance_df, "GeneImportance_RF_15gene_CV.csv", row.names = FALSE)

# Optional: Heatmap-style confusion matrix
cm_table <- as.table(cv_conf$table)
cm_df <- as.data.frame(cm_table)
ggplot(cm_df, aes(Prediction, Reference, fill = Freq)) +
   geom_tile(color = "white") +
   geom_text(aes(label = Freq), vjust = 0.5, size = 5) +
   scale_fill_gradient(low = "white", high = "#d7191c") +
   labs(title = "Confusion Matrix (CV): 15-Gene Classifier", fill = "Count") +
   theme_minimal()











#External dataset for validation - FUJIWARA
fujiwara=read.csv("data/Fujiwara_dataset/raw_counts_fujiwara.csv", sep = ",")
rownames(fujiwara) = fujiwara$ENSID
fujiwara = fujiwara[ ,-c(1)]
fujiwara = as.matrix(fujiwara)
dim(fujiwara)


# Read sheets
df_bx1 <- read_xlsx("data/Fujiwara_dataset/Both biopsies - refined dataset.xlsx", sheet = 1)
df_bx2 <- read_xlsx("data/Fujiwara_dataset/Both biopsies - refined dataset.xlsx", sheet = 2)

# Select relevant columns from first biopsy
df_bx1_sub <- df_bx1 %>%
   select(gct.name, Patient, Biopsy, Age, Platelet, ALB, AST, ALT, BMI, Diabetes, Histology.fibrosis)

# Select relevant columns from second biopsy, rename *_2 columns to base names
df_bx2_sub <- df_bx2 %>%
   select(gct.name = sample.names, Patient, Biopsy, Age, PLT_2, ALB_2, AST_2, ALT_2, BMI_2, Diabetes, Histology.fibrosis_2) %>%
   rename(
      Platelet = PLT_2,
      ALB = ALB_2,
      AST = AST_2,
      ALT = ALT_2,
      BMI = BMI_2,
      Histology.fibrosis = Histology.fibrosis_2
   )

# Add biopsy number for clarity
df_bx1_sub$biopsy_num <- 1
df_bx2_sub$biopsy_num <- 2

# Convert AST and other columns to numeric in both dfs before combining
cols_to_convert <- c("AST", "ALT", "Platelet", "ALB", "BMI", "Age", "Diabetes", "Histology.fibrosis")

# Helper function to convert columns to numeric safely
convert_to_numeric <- function(df, cols) {
   for (col in cols) {
      if (col %in% colnames(df)) {
         df[[col]] <- as.numeric(as.character(df[[col]]))
      }
   }
   return(df)
}

# Apply to both data frames
df_bx1_sub <- convert_to_numeric(df_bx1_sub, cols_to_convert)
df_bx2_sub <- convert_to_numeric(df_bx2_sub, cols_to_convert)

# Now bind rows
df_combined <- bind_rows(df_bx1_sub, df_bx2_sub)

#write.csv(df_combined, "Fujiwara_Metadata_bothbiopsies_refined_and_combined.csv")





# Assuming your 15 genes list:
selected_genes <- selected_15$ensembl_gene_id

# Subset expression matrix to 15 genes
fujiwara_15 <- fujiwara[selected_genes, , drop=FALSE]

# Get samples common to expression and metadata
common_samples <- intersect(colnames(fujiwara_15), df_combined$gct.name)

# Subset expression and metadata
fujiwara_15 <- fujiwara_15[, common_samples, drop=FALSE]
df_combined_sub <- df_combined %>% filter(gct.name %in% common_samples)




# z-score normalize gene expression (genes in rows, samples in columns)
fujiwara_15_z <- t(scale(t(fujiwara_15)))

# Convert to data.frame with samples as rows (for prediction)
expr_df <- as.data.frame(t(fujiwara_15_z))

# Add fibrosis label factor consistent with training
expr_df$fibrosis_bin <- ifelse(df_combined_sub$Histology.fibrosis >= 3, "Advanced", "Early")
expr_df$fibrosis_bin <- factor(expr_df$fibrosis_bin, levels = c("Early", "Advanced"))




# Predict fibrosis class
pred_fujiwara <- predict(rf_cv$finalModel, newdata = expr_df)

# Predict probabilities for ROC etc.
prob_fujiwara <- predict(rf_cv$finalModel, newdata = expr_df, type = "prob")[, "Advanced"]

# Evaluate performance
library(caret)
conf_fujiwara <- confusionMatrix(pred_fujiwara, expr_df$fibrosis_bin)
print(conf_fujiwara)





library(pROC)

roc_fujiwara <- roc(expr_df$fibrosis_bin, prob_fujiwara)
auc_fujiwara <- auc(roc_fujiwara)
cat("FUJIWARA AUC =", auc_fujiwara, "\n")

plot(roc_fujiwara, print.auc = TRUE, col = "darkgreen", lwd = 2,
     main = "FUJIWARA Dataset ROC: 15-Gene Classifier (F3–F4 vs F0–F2)")





library(dplyr)

AST_ULN <- 40

df_combined_sub <- df_combined_sub %>%
   mutate(
      Age = as.numeric(Age),
      AST = as.numeric(AST),
      ALT = as.numeric(ALT),
      Platelet = as.numeric(Platelet),
      BMI = as.numeric(BMI),
      Diabetes_num = as.numeric(Diabetes),
      ALB = as.numeric(ALB),
      
      # Avoid zeros or negative in denominators (optional)
      AST = ifelse(is.na(AST) | AST <= 0, NA, AST),
      ALT = ifelse(is.na(ALT) | ALT <= 0, NA, ALT),
      Platelet = ifelse(is.na(Platelet) | Platelet <= 0, NA, Platelet),
      Age = ifelse(is.na(Age) | Age <= 0, NA, Age),
      
      FIB4 = (Age * AST) / (Platelet * sqrt(ALT)),
      APRI = ((AST / AST_ULN) / Platelet) * 100,
      NFS = -1.675 + 0.037 * Age + 0.094 * BMI + 1.13 * Diabetes_num +
         0.99 * (AST / ALT) - 0.013 * Platelet - 0.66 * ALB
   )


write.csv(df_combined_sub, "Fujiwara_Metadata_bothbiopsies_refined_and_combined_top15.csv")










library(dplyr)
library(caret)
library(pROC)
library(tidyr)
library(ggplot2)

# --- Step 0: Define AST ULN ---
AST_ULN <- 40

# --- Step 1: Calculate NIT scores (if not done yet) ---
df_combined_sub <- df_combined_sub %>%
   mutate(
      Age = as.numeric(Age),
      AST = as.numeric(AST),
      ALT = as.numeric(ALT),
      Platelet = as.numeric(Platelet),
      BMI = as.numeric(BMI),
      Diabetes_num = as.numeric(Diabetes),
      ALB = ifelse("ALB" %in% colnames(df_combined_sub), as.numeric(ALB), NA_real_),
      
      AST = ifelse(is.na(AST) | AST <= 0, NA, AST),
      ALT = ifelse(is.na(ALT) | ALT <= 0, NA, ALT),
      Platelet = ifelse(is.na(Platelet) | Platelet <= 0, NA, Platelet),
      Age = ifelse(is.na(Age) | Age <= 0, NA, Age),
      
      FIB4 = (Age * AST) / (Platelet * sqrt(ALT)),
      APRI = ((AST / AST_ULN) / Platelet) * 100,
      NFS = ifelse(!is.na(ALB),
                   -1.675 + 0.037 * Age + 0.094 * BMI + 1.13 * Diabetes_num +
                      0.99 * (AST / ALT) - 0.013 * Platelet - 0.66 * ALB,
                   NA_real_)
   )

# --- Step 2: Create fibrosis_bin factor ---
df_combined_sub$fibrosis_bin <- ifelse(df_combined_sub$Histology.fibrosis >= 3, "Advanced", "Early")
df_combined_sub$fibrosis_bin <- factor(df_combined_sub$fibrosis_bin, levels = c("Early", "Advanced"))

# --- Step 3: RF predictions ---
# Ensure expr_df contains gene expression data for RF model predictions with fibrosis_bin
pred_rf <- predict(rf_cv$finalModel, newdata = expr_df)
prob_rf <- predict(rf_cv$finalModel, newdata = expr_df, type = "prob")[, "Advanced"]

# --- Step 4: NITs binary predictions with cutoffs ---
pred_fib4 <- factor(ifelse(df_combined_sub$FIB4 >= 1.45, "Advanced", "Early"), levels = c("Early", "Advanced"))
pred_apri <- factor(ifelse(df_combined_sub$APRI >= 0.7, "Advanced", "Early"), levels = c("Early", "Advanced"))
pred_nfs  <- factor(ifelse(df_combined_sub$NFS  >= -1.455, "Advanced", "Early"), levels = c("Early", "Advanced"))

# --- Step 5: Function to compute metrics (with AUC) ---
get_metrics <- function(pred, true_label, prob = NULL) {
   cm <- confusionMatrix(pred, true_label, positive = "Advanced")
   auc_val <- if (!is.null(prob) && !all(is.na(prob))) {
      tryCatch({
         roc_obj <- roc(true_label, prob, levels = c("Early", "Advanced"))
         auc(roc_obj)
      }, error = function(e) NA)
   } else {
      NA
   }
   data.frame(
      Accuracy = cm$overall["Accuracy"],
      Sensitivity = cm$byClass["Sensitivity"],
      Specificity = cm$byClass["Specificity"],
      AUC = auc_val
   )
}

# --- Step 6: Calculate metrics for all methods ---
metrics_rf <- get_metrics(pred_rf, expr_df$fibrosis_bin, prob = prob_rf)
metrics_fib4 <- get_metrics(pred_fib4, df_combined_sub$fibrosis_bin, prob = df_combined_sub$FIB4)
metrics_apri <- get_metrics(pred_apri, df_combined_sub$fibrosis_bin, prob = df_combined_sub$APRI)
metrics_nfs  <- get_metrics(pred_nfs,  df_combined_sub$fibrosis_bin, prob = df_combined_sub$NFS)

# --- Step 7: Combine results ---
metrics_df <- rbind(
   RF_Model = metrics_rf,
   FIB4 = metrics_fib4,
   APRI = metrics_apri,
   NFS = metrics_nfs
)

metrics_df$Method <- rownames(metrics_df)
metrics_df <- metrics_df[, c("Method", "Accuracy", "Sensitivity", "Specificity", "AUC")]

print(metrics_df)

# --- Step 8: Plot metrics comparison ---
# Convert AUC column to numeric (from 'auc' S3 class)
metrics_df$AUC <- as.numeric(metrics_df$AUC)

# Pivot to long format
metrics_long_15 <- metrics_df %>%
   pivot_longer(cols = Accuracy:AUC, names_to = "Metric", values_to = "Value")


ggplot(metrics_long_15, aes(x = Metric, y = Value, fill = Method)) +
   geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
   scale_fill_manual(values = c("purple", "gray40", "orange", "forestgreen")) +
   scale_y_continuous(breaks = seq(0.4, .9, 0.05)) +
   coord_cartesian(ylim = c(0.4, .9)) +   # zoom without dropping data
   labs(title = "Performance Comparison: RF Model vs NITs on Fujiwara Dataset",
        y = "Metric Value", x = "") +
   theme_minimal(base_size = 14) +
   theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, face = "bold")
   )

# --- Save the plot with white background ---
ggsave("Fujiwara_RF_vs_NITs_Performance_15.pdf", width = 8, height = 6, dpi = 300, bg = "white")







# --- Step 9: Combine 57-gene and 15-gene metrics ---
metrics_long_57$GeneSet <- "57 Genes"
metrics_long_15$GeneSet <- "Top 15 Genes"

combined_metrics <- bind_rows(metrics_long_57, metrics_long_15)[-c(5:16) , ]
combined_metrics$Method[1:4] = paste0(combined_metrics$Method[1:4], "_57")
combined_metrics$Method[5:8] = paste0(combined_metrics$Method[5:8], "_15")

# Optional: keep Method as factor for consistent ordering
combined_metrics$Method <- factor(combined_metrics$Method, levels = c("RF_Model_15", "RF_Model_57", "FIB4", "APRI", "NFS"))

# --- Step 10: Plot combined barplot ---
ggplot(combined_metrics, aes(x = Metric, y = Value, fill = Method)) +
   geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
   scale_fill_manual(values = c("mediumpurple", "purple4", "gray70", "gray50", "gray30")) +
   scale_y_continuous(breaks = seq(0.4, 0.95, 0.05)) +
   coord_cartesian(ylim = c(0.4, 0.95)) +
   labs(title = "Fujiwara: Our RF model vs NITs",
        y = "Metric Value", x = "") +
   theme_minimal(base_size = 14) +
   theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, face = "bold")
   )

# --- Save the combined plot ---
ggsave("Fujiwara_RF_vs_NITs_Performance_Combined.pdf", width = 10, height = 6, dpi = 300, bg = "white")


#Keep separately also the patients with paired biopsies
paired_patients = df_combined_sub[duplicated(df_combined_sub$Patient), ]$Patient
unique(paired_patients)
df_combined_sub1 = df_combined_sub[ is.element(df_combined_sub$Patient, paired_patients), ]
write.csv(df_combined_sub1, "Fujiwara_PairedBiopsies.csv")
