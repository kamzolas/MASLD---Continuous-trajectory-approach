# === LIBRARIES === #
library(randomForest)
library(caret)
library(ggplot2)
library(pROC)
library(dplyr)
library(tidyr)
library(readxl)

# === PATHS === #
template_path <- "data/metadata.csv"
zscores_path <- "data/z_scores.csv"
epos_meta_path <- "data/epos_metadata.csv"
epos_expr_path <- "data/Corrected_NormalisedCounts_EPOS2022.csv"

# === HELPER FUNCTION === #
run_rf_analysis <- function(selected_genes_path, gene_set_name,
                            template_path, zscores_path,
                            epos_meta_path, epos_expr_path) {
   
   # ---- Load internal data ----
   template <- read.csv(template_path)[-5, -c(1,3)]
   corrected <- read.csv(zscores_path, check.names = FALSE)
   rownames(corrected) <- corrected[,1]
   corrected <- corrected[,-1]
   
   selected_genes <- read.csv(selected_genes_path)$ensembl_gene_id
   selected_genes <- selected_genes[selected_genes %in% rownames(corrected)]
   
   expr <- corrected[selected_genes, ]
   template <- template[match(colnames(expr), template$Sample.name), ]
   stopifnot(all(template$Sample.name == colnames(expr)))
   
   # ---- Transpose and label ----
   data_train <- as.data.frame(t(expr))
   data_train$fibrosis_bin <- factor(ifelse(as.numeric(as.character(template$Fibrosis)) >= 3,
                                            "Advanced", "Early"),
                                     levels = c("Early", "Advanced"))
   
   # ---- RF Training ----
   set.seed(123)
   ctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE,
                        summaryFunction = twoClassSummary, savePredictions = "final")
   rf_cv <- train(fibrosis_bin ~ ., data = data_train, method = "rf",
                  trControl = ctrl, metric = "ROC")
   
   # ---- CV Performance ----
   cv_conf <- confusionMatrix(rf_cv$pred$pred, rf_cv$pred$obs)
   roc_cv <- roc(rf_cv$pred$obs, rf_cv$pred$Advanced)
   auc_cv <- auc(roc_cv)
   cat(gene_set_name, "- CV AUC =", auc_cv, "\n")
   
   # ---- Load EPoS data ----
   epos_meta <- read.csv(epos_meta_path)
   epos <- read.csv(epos_expr_path)
   rownames(epos) <- epos$X
   epos <- epos[,-1]
   colnames(epos) <- paste0(colnames(epos), " ")
   epos <- epos[, epos_meta$X]
   
   expr_epos <- epos[selected_genes, ]
   z_epos <- t(scale(t(as.matrix(expr_epos))))
   z_epos <- as.data.frame(t(z_epos))
   
   common_samples <- intersect(rownames(z_epos), epos_meta$X)
   z_epos <- z_epos[common_samples, ]
   epos_meta <- epos_meta[match(common_samples, epos_meta$X), ]
   
   z_epos$fibrosis_bin <- factor(ifelse(epos_meta$Fibrosis.stage >= 3,
                                        "Advanced", "Early"),
                                 levels = c("Early", "Advanced"))
   
   # ---- Align & Impute ----
   expected_genes <- setdiff(colnames(data_train), "fibrosis_bin")
   missing <- setdiff(expected_genes, colnames(z_epos))
   for (g in missing) z_epos[[g]] <- NA
   z_epos <- z_epos[, c(expected_genes, "fibrosis_bin")]
   
   means_train <- colMeans(data_train[, expected_genes])
   for (g in expected_genes) {
      if (anyNA(z_epos[[g]])) z_epos[[g]][is.na(z_epos[[g]])] <- means_train[g]
   }
   
   # ---- Predict on EPoS ----
   pred_epos <- predict(rf_cv$finalModel, newdata = z_epos)
   prob_epos <- predict(rf_cv$finalModel, newdata = z_epos, type = "prob")[,"Advanced"]
   conf_epos <- confusionMatrix(pred_epos, z_epos$fibrosis_bin)
   roc_epos <- roc(z_epos$fibrosis_bin, prob_epos)
   auc_epos <- auc(roc_epos)
   cat(gene_set_name, "- EPoS AUC =", auc_epos, "\n")
   
   # ---- ROC Comparison with FIB-4 ----
   roc_fib4 <- roc(z_epos$fibrosis_bin, epos_meta$FIB4)
   auc_fib4 <- auc(roc_fib4)
   
   roc_df <- bind_rows(
      data.frame(FPR = 1 - roc_epos$specificities, TPR = roc_epos$sensitivities,
                 Method = paste0("RF ", gene_set_name, " (AUC=", round(auc_epos,2), ")")),
      data.frame(FPR = 1 - roc_fib4$specificities, TPR = roc_fib4$sensitivities,
                 Method = paste0("FIB-4 (AUC=", round(auc_fib4,2), ")"))
   )
   
   # ---- Metric comparison ----
   fib4_pred <- factor(ifelse(epos_meta$FIB4 >= 1.45, "Advanced", "Early"),
                       levels = c("Early", "Advanced"))
   
   comparison_df <- data.frame(
      Metric = c("Accuracy", "Sensitivity", "Specificity", "AUC"),
      RF_Model = c(
         round(conf_epos$overall["Accuracy"], 3),
         round(conf_epos$byClass["Sensitivity"], 3),
         round(conf_epos$byClass["Specificity"], 3),
         round(auc_epos, 3)
      ),
      FIB4 = c(
         round(confusionMatrix(fib4_pred, z_epos$fibrosis_bin)$overall["Accuracy"], 3),
         round(confusionMatrix(fib4_pred, z_epos$fibrosis_bin)$byClass["Sensitivity"], 3),
         round(confusionMatrix(fib4_pred, z_epos$fibrosis_bin)$byClass["Specificity"], 3),
         round(auc_fib4, 3)
      )
   )
   
   comparison_long <- pivot_longer(comparison_df, cols = c("RF_Model", "FIB4"),
                                   names_to = "Model", values_to = "Value")
   comparison_long$RF_Set <- gene_set_name
   
   return(list(rf_model = rf_cv,
               roc_df = roc_df,
               metrics_df = comparison_long))
}

# === RUN FOR MULTIPLE GENE SETS === #
gene_sets <- list(
   "57-gene" = "data/57_BMs.csv",
   "15-gene" = "data/top_15_BMs_from_ELBOW.csv"
)

results_list <- lapply(names(gene_sets), function(gs_name) {
   run_rf_analysis(selected_genes_path = gene_sets[[gs_name]],
                   gene_set_name = gs_name,
                   template_path = template_path,
                   zscores_path = zscores_path,
                   epos_meta_path = epos_meta_path,
                   epos_expr_path = epos_expr_path)
})
names(results_list) <- names(gene_sets)

# === COMBINE ROC CURVES & PLOT === #
roc_all <- do.call(rbind, lapply(results_list, function(x) x$roc_df))
roc_all <- roc_all %>% filter(Method != "FIB-4" | !duplicated(FPR))  # keep one FIB-4

p_roc_all <- ggplot(roc_all, aes(FPR, TPR, color = Method)) +
   geom_line(size = .8) +
   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
   scale_color_manual(values = c("gray70", "mediumpurple", "purple4")) +
   labs(title = "ROC curves: 57 BMs vs Top15 BMs vs FIB-4",
        x = "1 - Specificity", y = "Sensitivity") +
   coord_fixed() +
   theme_minimal(base_size = 12) +
   theme(legend.position = "bottom",
         panel.background = element_rect(fill = "white", color = NA),
         plot.background = element_rect(fill = "white", color = NA))
print(p_roc_all)
ggsave("ROC_Comparison_All_RF_vs_FIB4.png", plot = p_roc_all,
       width = 7, height = 6, dpi = 300)

# === COMBINE METRICS & PLOT === #
metrics_all <- do.call(rbind, lapply(results_list, function(x) x$metrics_df))
metrics_all = metrics_all[ -c(10,12,14,16), ]
metrics_all$Model <- ifelse(metrics_all$Model == "FIB4", "FIB-4", metrics_all$Model)

plot_bar_all <- ggplot(metrics_all, aes(x = Metric, y = Value, fill = interaction(RF_Set, Model))) +
   geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
   geom_text(aes(label = round(Value,2)), position = position_dodge(width=0.7), vjust=-0.5, size=3) +
   scale_fill_manual(
      values = c("gray70","mediumpurple", "purple4"),
      labels = c("FIB-4", "Top15 BMs", "57 BMs")  # rename legend
   ) +
   scale_y_continuous(breaks = seq(0.6, 1, 0.05)) +
   coord_cartesian(ylim = c(0.6,1)) +
   labs(title = "Performance: 57 BMs vs Top15 BMs vs FIB-4", y="Metric Value", x="") +
   theme_minimal(base_size = 14) +
   theme(legend.title = element_blank(),
         plot.title = element_text(hjust=0.5, face="bold"))

print(plot_bar_all)

ggsave("Metrics_Comparison_All_RF_vs_FIB4.png", plot = plot_bar_all,
       width=7, height=5, dpi=300, bg="white")

