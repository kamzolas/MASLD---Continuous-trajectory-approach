### ===============================
### === SETUP =====================
### ===============================
library(randomForest)
library(caret)
library(ggplot2)
library(pROC)
library(readxl)
library(tidyverse)

set.seed(123)

ctrl <- trainControl(
   method = "cv",
   number = 5,
   classProbs = TRUE,
   summaryFunction = twoClassSummary,
   savePredictions = "final"
)

AST_ULN <- 40

### ===============================
### === HELPER FUNCTIONS ==========
### ===============================

# ---- Load and preprocess training data ----
load_training_data <- function(template_path, zscore_path, gene_list_path) {
   
   template <- read.csv(template_path)[-5, -c(1,3)]
   
   expr <- read.csv(zscore_path, check.names = FALSE)
   rownames(expr) <- expr[,1]
   expr <- expr[,-1]
   
   genes <- read.csv(gene_list_path)$ensembl_gene_id
   genes <- genes[genes %in% rownames(expr)]
   
   expr <- expr[genes, ]
   
   template <- template[match(colnames(expr), template$Sample.name), ]
   stopifnot(all(template$Sample.name == colnames(expr)))
   
   expr <- as.data.frame(t(expr))
   
   expr$fibrosis_bin <- factor(
      ifelse(as.numeric(as.character(template$Fibrosis)) >= 3, "Advanced", "Early"),
      levels = c("Early","Advanced")
   )
   
   return(expr)
}

# ---- Train RF model ----
train_rf <- function(expr) {
   train(fibrosis_bin ~ ., data = expr,
         method = "rf",
         trControl = ctrl,
         metric = "ROC")
}

# ---- Load Fujiwara metadata ----
load_fujiwara_meta <- function(path) {
   
   df1 <- read_xlsx(path, sheet = 1) %>%
      select(gct.name, Patient, Biopsy, Age, Platelet, ALB, AST, ALT, BMI, Diabetes, Histology.fibrosis) %>%
      mutate(across(c(Age, Platelet, ALB, AST, ALT, BMI, Diabetes), as.numeric))
   
   df2 <- read_xlsx(path, sheet = 2) %>%
      select(gct.name = sample.names, Patient, Biopsy, Age, PLT_2, ALB_2, AST_2, ALT_2, BMI_2, Diabetes, Histology.fibrosis_2) %>%
      rename(
         Platelet = PLT_2,
         ALB = ALB_2,
         AST = AST_2,
         ALT = ALT_2,
         BMI = BMI_2,
         Histology.fibrosis = Histology.fibrosis_2
      ) %>%
      mutate(across(c(Age, Platelet, ALB, AST, ALT, BMI, Diabetes), as.numeric))
   
   df <- bind_rows(df1, df2)
   
   return(df)
}

# ---- Load Fujiwara expression ----
load_fujiwara_expr <- function(path) {
   expr <- read.csv(path)
   rownames(expr) <- expr$ENSID
   expr <- expr[,-1]
   as.matrix(expr)
}

# ---- Prepare Fujiwara dataset ----
prepare_fujiwara <- function(expr, meta, genes) {
   
   expr <- expr[genes, , drop=FALSE]
   common <- intersect(colnames(expr), meta$gct.name)
   
   expr <- expr[, common]
   meta <- meta %>% filter(gct.name %in% common)
   
   expr_z <- t(scale(t(expr)))
   expr_df <- as.data.frame(t(expr_z))
   
   expr_df$fibrosis_bin <- factor(
      ifelse(meta$Histology.fibrosis >= 3, "Advanced", "Early"),
      levels = c("Early","Advanced")
   )
   
   list(expr = expr_df, meta = meta)
}

# ---- Add NIT scores ----
add_nits <- function(meta) {
   meta %>%
      mutate(
         FIB4 = (Age * AST) / (Platelet * sqrt(ALT)),
         APRI = ((AST / AST_ULN) / Platelet) * 100,
         NFS = -1.675 + 0.037 * Age + 0.094 * BMI + 1.13 * Diabetes +
            0.99 * (AST / ALT) - 0.013 * Platelet - 0.66 * ALB,
         fibrosis_bin = factor(ifelse(Histology.fibrosis >= 3, "Advanced","Early"),
                               levels = c("Early","Advanced"))
      )
}

# ---- Metrics ----
get_metrics <- function(pred, truth, prob=NULL) {
   cm <- confusionMatrix(pred, truth, positive="Advanced")
   
   auc_val <- if (!is.null(prob)) {
      tryCatch(auc(roc(truth, prob)), error=function(e) NA)
   } else NA
   
   data.frame(
      Accuracy = cm$overall["Accuracy"],
      Sensitivity = cm$byClass["Sensitivity"],
      Specificity = cm$byClass["Specificity"],
      AUC = as.numeric(auc_val)
   )
}

# ---- Full evaluation ----
run_pipeline <- function(gene_path, label, training_paths, fujiwara_expr, fujiwara_meta) {
   
   expr_train <- load_training_data(training_paths$template,
                                    training_paths$zscore,
                                    gene_path)
   
   rf <- train_rf(expr_train)
   
   genes <- colnames(expr_train)[-ncol(expr_train)]
   
   fuji <- prepare_fujiwara(fujiwara_expr, fujiwara_meta, genes)
   fuji$meta <- add_nits(fuji$meta)
   
   pred_rf  <- predict(rf$finalModel, fuji$expr)
   prob_rf  <- predict(rf$finalModel, fuji$expr, type="prob")[,"Advanced"]
   
   pred_fib4 <- factor(ifelse(fuji$meta$FIB4 >= 1.45,"Advanced","Early"), levels=c("Early","Advanced"))
   pred_apri <- factor(ifelse(fuji$meta$APRI >= 0.7,"Advanced","Early"), levels=c("Early","Advanced"))
   pred_nfs  <- factor(ifelse(fuji$meta$NFS >= -1.455,"Advanced","Early"), levels=c("Early","Advanced"))
   
   metrics <- rbind(
      RF_Model = get_metrics(pred_rf, fuji$expr$fibrosis_bin, prob_rf),
      FIB4 = get_metrics(pred_fib4, fuji$meta$fibrosis_bin, fuji$meta$FIB4),
      APRI = get_metrics(pred_apri, fuji$meta$fibrosis_bin, fuji$meta$APRI),
      NFS  = get_metrics(pred_nfs, fuji$meta$fibrosis_bin, fuji$meta$NFS)
   )
   
   metrics$Method <- rownames(metrics)
   metrics$GeneSet <- label
   
   return(list(
      metrics = metrics,
      truth = fuji$expr$fibrosis_bin,
      prob_rf = prob_rf,
      meta = fuji$meta
   ))
}

### ===============================
### === RUN ANALYSIS ==============
### ===============================

paths <- list(
   template = "data/metadata.csv",
   zscore   = "data/z_scores.csv"
)

fujiwara_expr <- load_fujiwara_expr("data/raw_counts_fujiwara.csv")
fujiwara_meta <- load_fujiwara_meta("data/Both biopsies - refined dataset.xlsx")

res_57 <- run_pipeline("data/57_BMs.csv", "57 Genes", paths, fujiwara_expr, fujiwara_meta)
res_15 <- run_pipeline("data/top_15_BMs_from_ELBOW.csv", "Top15 Genes", paths, fujiwara_expr, fujiwara_meta)
res_3  <- run_pipeline("data/External_3_gene_panel_for_additional_validation.csv", "3-gene", paths, fujiwara_expr, fujiwara_meta)

results <- bind_rows(
   res_57$metrics,
   res_15$metrics,
   res_3$metrics
)


truth   <- res_57$truth   # same truth for all models

prob_57 <- res_57$prob_rf
prob_15 <- res_15$prob_rf
prob_3  <- res_3$prob_rf

meta    <- res_57$meta    # contains FIB4/APRI/NFS


# Separate RF models
rf_15 <- results %>% filter(GeneSet == "Top15 Genes", Method == "RF_Model")
rf_57 <- results %>% filter(GeneSet == "57 Genes", Method == "RF_Model")
rf_3 <- results %>% filter(GeneSet == "3-gene", Method == "RF_Model")

# Take NITs only once (they are identical across GeneSets)
nits <- results %>%
   filter(GeneSet == "57 Genes", Method %in% c("FIB4","APRI","NFS"))

# Rename RF models
rf_15$Method <- "RF_Model_15"
rf_57$Method <- "RF_Model_57"
rf_3$Method <- "3-Gene_Published"

# Combine everything
metrics_clean <- bind_rows(
   rf_15,
   rf_57,
   rf_3,    # <-- NEW
   nits
)




metrics_all <- metrics_clean %>%
   mutate(
      RF_Set = case_when(
         Method == "RF_Model_15" ~ "15-gene",
         Method == "RF_Model_57" ~ "57-gene",
         Method == "3-Gene_Published"  ~ "3-gene published panel",
         TRUE ~ "57-gene"
      ),
      Model = case_when(
         grepl("RF_Model", Method) ~ "RF_Model",
         Method == "3-Gene_Published" ~ "RF_Model",
         Method == "FIB4" ~ "FIB-4",
         Method == "APRI" ~ "APRI",
         Method == "NFS" ~ "NFS"
      )
   )



delong_test <- function(truth, prob1, prob2) {
   roc1 <- roc(truth, prob1)
   roc2 <- roc(truth, prob2)
   
   test <- roc.test(roc1, roc2, method = "delong")
   
   data.frame(
      AUC_1 = as.numeric(auc(roc1)),
      AUC_2 = as.numeric(auc(roc2)),
      p_value = test$p.value
   )
}



delong_results <- bind_rows(
   
   data.frame(
      Comparison = "RF_57 vs RF_15",
      delong_test(truth, prob_57, prob_15)
   ),
   
   data.frame(
      Comparison = "RF_57 vs 3-gene",
      delong_test(truth, prob_57, prob_3)
   ),
   
   data.frame(
      Comparison = "RF_57 vs FIB-4",
      delong_test(truth, prob_57, meta$FIB4)
   ),
   
   data.frame(
      Comparison = "RF_57 vs APRI",
      delong_test(truth, prob_57, meta$APRI)
   ),
   
   data.frame(
      Comparison = "RF_57 vs NFS",
      delong_test(truth, prob_57, meta$NFS)
   )
)



get_auc_ci <- function(truth, prob) {
   roc_obj <- roc(truth, prob)
   ci_vals <- ci.auc(roc_obj, method = "delong")
   
   data.frame(
      AUC = as.numeric(auc(roc_obj)),
      CI_low = ci_vals[1],
      CI_high = ci_vals[3]
   )
}



auc_ci_df <- bind_rows(
   
   data.frame(RF_Set="15-gene", Model="RF_Model",
              Metric="AUC",
              get_auc_ci(truth, prob_15)),
   
   data.frame(RF_Set="57-gene", Model="RF_Model",
              Metric="AUC",
              get_auc_ci(truth, prob_57)),
   
   data.frame(RF_Set="3-gene published panel", Model="RF_Model",
              Metric="AUC",
              get_auc_ci(truth, prob_3)),
   
   data.frame(RF_Set="57-gene", Model="FIB-4",
              Metric="AUC",
              get_auc_ci(truth, meta$FIB4)),
   
   data.frame(RF_Set="57-gene", Model="APRI",
              Metric="AUC",
              get_auc_ci(truth, meta$APRI)),
   
   data.frame(RF_Set="57-gene", Model="NFS",
              Metric="AUC",
              get_auc_ci(truth, meta$NFS))
)




metrics_all <- metrics_all %>%
   pivot_longer(cols = Accuracy:AUC,
                names_to = "Metric",
                values_to = "Value")

metrics_all <- left_join(metrics_all, auc_ci_df,
                         by = c("RF_Set","Model","Metric"))


metrics_all$fill_group <- interaction(metrics_all$RF_Set, metrics_all$Model)

metrics_all$fill_group <- factor(
   metrics_all$fill_group,
   levels = c(
      "15-gene.RF_Model",
      "57-gene.RF_Model",
      "57-gene.FIB-4",
      "57-gene.APRI",
      "57-gene.NFS",
      "3-gene published panel.RF_Model"
   )
)



metrics_all$Metric <- factor(
   metrics_all$Metric,
   levels = c("AUC", "Sensitivity", "Specificity", "Accuracy")
)


plot_bar_all <- ggplot(metrics_all,
                       aes(x = Metric,
                           y = Value,
                           fill = fill_group)) +
   
   geom_bar(stat = "identity",
            position = position_dodge(width = 0.7),
            width = 0.6) +
   
   geom_text(
      aes(
         label = round(Value, 2),
         y = ifelse(Metric == "AUC", Value + 0.12, Value + 0.02)
      ),
      position = position_dodge(width = 0.7),
      size = 2.2
   ) +
   
   geom_errorbar(
      data = subset(metrics_all, Metric == "AUC"),
      aes(ymin = CI_low, ymax = CI_high),
      position = position_dodge(width = 0.7),
      width = 0.2,
      linewidth = 0.7
   ) +
   
   scale_fill_manual(
      values = c(
         "15-gene.RF_Model" = "mediumpurple",
         "57-gene.RF_Model" = "purple4",
         "57-gene.FIB-4" = "gray70",
         "57-gene.APRI" = "gray50",
         "57-gene.NFS" = "gray30",
         "3-gene published panel.RF_Model" = "black"
      ),
      labels = c(
         "Top15 BMs",
         "57 BMs",
         "FIB-4",
         "APRI",
         "NFS",
         "Verschuren et al"
      )
   ) +
   
   labs(title = "Performance comparison of biomarker panels and clinical scores",
        y = "Metric Value",
        x = "") +
   
   theme_minimal(base_size = 14) +
   
   theme(
      legend.title = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 10.2),       # smaller text
      legend.key.height = unit(0.5, "lines"),     # reduce key height
      legend.key.width = unit(1.5, "lines"),      # reduce key width
      plot.title = element_text(hjust = 0.5, face = "bold")
   ) +
   
   guides(fill = guide_legend(nrow = 1))  # force single-line legend

print(plot_bar_all)

ggsave("Metrics_Comparison_All.pdf",
       plot = plot_bar_all,
       width = 7, height = 5.5, dpi = 300, bg = "white")


