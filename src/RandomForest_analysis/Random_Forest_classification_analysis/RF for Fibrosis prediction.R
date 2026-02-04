### === SETUP === ###
library(randomForest)
library(caret)
library(ggplot2)
library(pROC)
library(reshape2)
library(readxl)
library(biomaRt)

### === LOAD DATA === ###
# Load patient metadata (remove row 5 and unnecessary columns)
template <- read.csv("data/metadata.csv")[-5, -c(1,3)]

# Load z-score normalised gene expression matrix
corrected_ucamsanyal <- read.csv("data/z_scores.csv", check.names = FALSE)
rownames(corrected_ucamsanyal) = corrected_ucamsanyal[,1]
corrected_ucamsanyal = corrected_ucamsanyal[,-1]

# Load 57  bms - after changing some parameters in the MASLD network
bms <- read.csv("new_markers.tsv", sep = "\t", stringsAsFactors = FALSE)

# Connect to Ensembl (human genes)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Map gene symbols to Ensembl gene IDs
mapping <- getBM(
   attributes = c("ensembl_gene_id", "external_gene_name"),
   filters = "external_gene_name",
   values = bms$gene_symbol,
   mart = ensembl
)

# Check results
head(mapping)

# Optionally merge back to your original list to maintain order
bms_merged <- merge(bms, mapping, by.x = "gene_symbol", by.y = "external_gene_name", all.x = TRUE)

# Reorder columns if you prefer
selected_57 <- bms_merged[, c("ensembl_gene_id", "gene_symbol")]
rm(bms, bms_merged, ensembl, mapping)
selected_57 = selected_57[!duplicated(selected_57$gene_symbol),]
colnames(selected_57) = c("ensembl_gene_id", "external_gene_name")
write.csv(selected_57, "57_BMs.csv")

### === PREPROCESS EXPRESSION DATA === ###
# Subset expression matrix to 57 selected genes
selected_57 = selected_57[ is.element(selected_57$ensembl_gene_id, rownames(corrected_ucamsanyal)), ]
corrected_ucamsanyal = corrected_ucamsanyal[ selected_57$ensembl_gene_id, ]


# Ensure sample order matches metadata
template <- template[match(colnames(corrected_ucamsanyal), template$Sample.name), ]
stopifnot(all(template$Sample.name == colnames(corrected_ucamsanyal)))

# Transpose expression matrix: patients as rows
expr_57 <- as.data.frame(t(corrected_ucamsanyal))

### === DEFINE BINARY LABEL: EARLY (F0–F2) VS ADVANCED (F3–F4) === ###
expr_57$fibrosis_bin <- ifelse(as.numeric(as.character(template$Fibrosis)) >= 3, "Advanced", "Early")
expr_57$fibrosis_bin <- factor(expr_57$fibrosis_bin, levels = c("Early", "Advanced"))

### === RANDOM FOREST CLASSIFIER WITH CROSS-VALIDATION === ###
set.seed(123)
ctrl <- trainControl(method = "cv",
                     number = 5,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary,
                     savePredictions = "final")

# Train RF model using caret with CV
rf_cv <- train(fibrosis_bin ~ ., data = expr_57,
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
     main = "5-Fold CV ROC: 57-Gene Classifier (F3–F4 vs F0–F2)")

### === VARIABLE IMPORTANCE (from final RF model) === ###
varImpPlot(rf_cv$finalModel, main = "Gene Importance: Early vs Advanced Fibrosis")

### === EXPORT PREDICTIONS AND IMPORTANCE === ###
# Save CV predictions
write.csv(rf_cv$pred, "RF_57gene_CV_Predictions.csv", row.names = FALSE)

# Export variable importance
importance_df <- as.data.frame(importance(rf_cv$finalModel))
importance_df$Gene <- rownames(importance_df)
importance_df <- importance_df[order(importance_df$MeanDecreaseGini, decreasing = TRUE), ]
# Ensure both data frames have a common column name for the merge - Merge by ENS_ID
importance_df$ensembl_gene_id <- importance_df$Gene
importance_with_symbols <- merge(importance_df, selected_57, by = "ensembl_gene_id", all.x = TRUE)
importance_with_symbols <- importance_with_symbols[, c("ensembl_gene_id", "external_gene_name", "MeanDecreaseGini")]
write.csv(importance_with_symbols, "GeneImportance_RF_57gene_CV.csv", row.names = FALSE)

# Optional: Heatmap-style confusion matrix
cm_table <- as.table(cv_conf$table)
cm_df <- as.data.frame(cm_table)
ggplot(cm_df, aes(Prediction, Reference, fill = Freq)) +
   geom_tile(color = "white") +
   geom_text(aes(label = Freq), vjust = 0.5, size = 5) +
   scale_fill_gradient(low = "white", high = "#d7191c") +
   labs(title = "Confusion Matrix (CV): 57-Gene Classifier", fill = "Count") +
   theme_minimal()









#External dataset for validation - FUJIWARA
fujiwara=read.csv("data/Fujiwara_dataset/raw_counts_fujiwara.csv", sep = ",")
rownames(fujiwara) = fujiwara$ENSID
fujiwara = fujiwara[ ,-c(1)]
fujiwara = as.matrix(fujiwara)



#First biopsy
template1 <- read_xlsx("data/Fujiwara_dataset/MASLD_data.xlsx", sheet = 1)
# Convert tibble to data frame
template1 <- as.data.frame(template1)
# Set first column as row names
rownames(template1) <- template1$gct.name
template1[,1] = 1
colnames(template1)[1] = "Biopsy"

#Second biopsy
template2 <- read_xlsx("data/Fujiwara_dataset/MASLD_data.xlsx", sheet = 2)
# Convert tibble to data frame
template2 <- as.data.frame(template2)
# Set first column as row names
rownames(template2) <- template2$sample.names
template2[,1] = 2
colnames(template2)[1] = "Biopsy"


#merge both templates
# Select relevant columns from each template
template1_selected <- template1[, c("Biopsy", "Patient", "Histology.fibrosis", "Histology.steatosis", 
                                    "Histology.ballooning", "Histology.inflammation", "Histology.NAS")]
template2_selected <- template2[, c("Biopsy", "Patient", "Histology.fibrosis_2", "Histology.steatosis_2", 
                                    "Histology.ballooning_2", "Histology.inflammation_2", "NAS_2")]
colnames(template2_selected) = colnames(template1_selected)

# Merge the two templates by Patient
template_fuji <- rbind(template1_selected, template2_selected)

# Print the first few rows to check the result
head(template_fuji)
rm(template1, template1_selected, template2, template2_selected)

dim(fujiwara)
unique(is.element(colnames(fujiwara), rownames(template_fuji)))




#Validate in Fujiwara
### === VALIDATION ON FUJIWARA ET AL. 2022 === ###

# 1. Load & preprocess Fujiwara data
# Assuming 'fujiwara' is a gene (rows) × sample (columns) matrix like the original
# And rownames are ENSEMBL IDs

# Subset to the same 57 genes
fuji_57 <- fujiwara[selected_57$ensembl_gene_id, ]

# Z-score normalisation (per gene across samples)
fuji_57_z <- t(scale(t(as.matrix(fuji_57))))  # genes in rows, samples in columns → standardize by row

# Transpose to samples as rows (as required by RF)
fuji_57_z <- as.data.frame(t(fuji_57_z))

# Match sample names in metadata
# Ensure rownames(fuji_57_z) match rownames(template_fuji) or vice versa
common_samples <- intersect(rownames(fuji_57_z), rownames(template_fuji))
fuji_57_z <- fuji_57_z[common_samples, ]
template_fuji <- template_fuji[common_samples, ]

# 2. Define early vs advanced fibrosis in Fujiwara
fuji_57_z$fibrosis_bin <- ifelse(template_fuji$Histology.fibrosis >= 3, "Advanced", "Early")
fuji_57_z$fibrosis_bin <- factor(fuji_57_z$fibrosis_bin, levels = c("Early", "Advanced"))

# 3. Predict using trained model (from UCAM)
# Use rf_cv$finalModel from the previous part (it’s trained on UCAM)

pred_fuji <- predict(rf_cv$finalModel, newdata = fuji_57_z)
prob_fuji <- predict(rf_cv$finalModel, newdata = fuji_57_z, type = "prob")[,"Advanced"]

# 4. Evaluate performance
conf_fuji <- confusionMatrix(pred_fuji, fuji_57_z$fibrosis_bin)
print(conf_fuji)

# ROC for Fujiwara
roc_fuji <- roc(fuji_57_z$fibrosis_bin, prob_fuji)
auc_fuji <- auc(roc_fuji)
cat("Fujiwara AUC =", auc_fuji, "\n")

plot(roc_fuji, print.auc = TRUE, col = "darkgreen", lwd = 2,
     main = "Fujiwara ROC: 57-Gene Classifier (F3–F4 vs F0–F2)")

# 5. Export predictions
fuji_output <- data.frame(
   Sample = rownames(fuji_57_z),
   True_Label = fuji_57_z$fibrosis_bin,
   Predicted_Label = pred_fuji,
   Prob_Advanced = prob_fuji,
   True_Fibrosis = template_fuji$Histology.fibrosis
)
write.csv(fuji_output, "Predictions_RF_57gene_Fujiwara.csv", row.names = FALSE)


# Heatmap-style confusion matrix
conf_fuji <- confusionMatrix(pred_fuji, fuji_57_z$fibrosis_bin)

# === HEATMAP-STYLE CONFUSION MATRIX FOR FUJIWARA ===
cm_fuji_table <- as.table(conf_fuji$table)
cm_fuji_df <- as.data.frame(cm_fuji_table)

ggplot(cm_fuji_df, aes(Prediction, Reference, fill = Freq)) +
   geom_tile(color = "white") +
   geom_text(aes(label = Freq), vjust = 0.5, size = 5) +
   scale_fill_gradient(low = "white", high = "#1a9641") +
   labs(title = "Fujiwara: Confusion Matrix (57-Gene Classifier)", fill = "Count") +
   theme_minimal()














# Load patchwork for combining plots
library(patchwork)

### === 1. ROC CURVES === ###
# Training (CV)
roc_cv <- roc(rf_cv$pred$obs, rf_cv$pred$Advanced)
auc_cv <- auc(roc_cv)

p_roc_cv <- ggplot(data = data.frame(
   FPR = 1 - roc_cv$specificities,
   TPR = roc_cv$sensitivities
), aes(x = FPR, y = TPR)) +
   geom_line(color = "darkblue", size = 1.2) +
   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +  # random classifier line
   annotate("text", x = 0.65, y = 0.1, label = paste0("AUC = ", round(auc_cv, 3)), size = 5) +
   labs(title = "CV ROC Curve (Training Set)", x = "1 - Specificity", y = "Sensitivity") +
   theme_minimal() +
   coord_cartesian(xlim = c(0,1), ylim = c(0,1))

# Fujiwara
roc_fuji <- roc(fuji_57_z$fibrosis_bin, prob_fuji)
auc_fuji <- auc(roc_fuji)

p_roc_fuji <- ggplot(data = data.frame(
   FPR = 1 - roc_fuji$specificities,
   TPR = roc_fuji$sensitivities
), aes(x = FPR, y = TPR)) +
   geom_line(color = "darkgreen", size = 1.2) +
   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +  # random classifier line
   annotate("text", x = 0.65, y = 0.1, label = paste0("AUC = ", round(auc_fuji, 3)), size = 5) +
   labs(title = "ROC Curve (Fujiwara Validation)", x = "1 - Specificity", y = "Sensitivity") +
   theme_minimal() +
   coord_cartesian(xlim = c(0,1), ylim = c(0,1))

### === 2. CONFUSION MATRICES === ###
# Training confusion matrix
cm_cv_df <- as.data.frame(as.table(cv_conf$table))
p_cm_cv <- ggplot(cm_cv_df, aes(Prediction, Reference, fill = Freq)) +
   geom_tile(color = "white") +
   geom_text(aes(label = Freq), size = 5) +
   scale_fill_gradient(low = "white", high = "#d7191c") +
   labs(title = "Confusion Matrix (Training CV)", fill = "Count") +
   theme_minimal()

# Fujiwara confusion matrix
cm_fuji_df <- as.data.frame(as.table(conf_fuji$table))
p_cm_fuji <- ggplot(cm_fuji_df, aes(Prediction, Reference, fill = Freq)) +
   geom_tile(color = "white") +
   geom_text(aes(label = Freq), size = 5) +
   scale_fill_gradient(low = "white", high = "#1a9641") +
   labs(title = "Confusion Matrix (Fujiwara)", fill = "Count") +
   theme_minimal()

### === 3. VARIABLE IMPORTANCE === ###
imp_df <- importance_with_symbols  # from earlier code
imp_df = imp_df[rev(order(imp_df$MeanDecreaseGini)), ]
#imp_df <- imp_df[1:20, ]  # top 20 genes for clean plotting

p_importance <- ggplot(imp_df, aes(x = reorder(external_gene_name, MeanDecreaseGini), y = MeanDecreaseGini)) +
   geom_col(fill = "steelblue") +
   coord_flip() +
   labs(title = "Gene Importance", x = "Gene", y = "Mean Decrease Gini") +
   theme(
      axis.text.y = element_text(size = 6)   # change y-axis text size
      
   )

### === COMBINE ALL PLOTS INTO ONE FIGURE === ###
final_plot <- (p_roc_cv | p_roc_fuji) /
   (p_cm_cv | p_cm_fuji) /
   p_importance +
   plot_layout(heights = c(1, 1, 1.3))  # more space for importance plot

### === SAVE FINAL FIGURE === ###
ggsave("Classifier_Summary_Figure.png", final_plot, width = 14, height = 12, dpi = 300)






#Selection of the top genes according to MeanDecreaseGini score - based on the ELBOW method
top15 = imp_df[1:15,]
top15 = top15[ ,c("ensembl_gene_id", "external_gene_name") ]
rownames(top15) = 1:15
write.csv(top15, "top_15_BMs_from_ELBOW.csv")





