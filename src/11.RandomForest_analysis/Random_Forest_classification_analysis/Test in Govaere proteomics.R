library(tidyverse)
library(caret)
library(randomForest)
library(pROC)
library(readxl)

### === 1. TRAIN RF MODEL FUNCTION === ###
train_rf_model <- function(gene_file, expr_file, template_file) {
   
   set.seed(123)
   
   template <- read.csv(template_file)[-5, -c(1,3)]
   expr <- read.csv(expr_file, check.names = FALSE)
   genes <- read.csv(gene_file)[,-1]
   
   rownames(expr) <- expr[,1]
   expr <- expr[,-1]
   genes <- genes[genes$ensembl_gene_id %in% rownames(expr), ]
   
   expr <- expr[genes$ensembl_gene_id, ]
   template <- template[match(colnames(expr), template$Sample.name), ]
   
   df <- as.data.frame(t(expr))
   df$fibrosis_bin <- factor(ifelse(as.numeric(template$Fibrosis) >= 3,
                                    "Advanced", "Early"),
                             levels = c("Early","Advanced"))
   
   trainControl_obj <- trainControl(method="cv", number=5,
                                    classProbs=TRUE,
                                    summaryFunction=twoClassSummary)
   
   train(fibrosis_bin ~ ., data=df,
         method="rf",
         trControl=trainControl_obj,
         metric="ROC")
}

### === 2. PREPROCESS GOVAERE FUNCTION === ###
prepare_govaere <- function(gene_map) {
   
   govaere <- read_excel("data/Plasma_Proteomics_dataset/42255_2023_775_MOESM3_ESM.xlsx")[-c(1:4), -c(1:6,8)] %>%
      rename(sample_id = 1) %>%
      distinct(sample_id, .keep_all = TRUE) %>%
      mutate(across(-c(sample_id,2), as.double))
   
   govaere <- as.data.frame(govaere)
   govaere <- govaere[!duplicated(govaere$sample_id) & !is.na(govaere$sample_id), ]
   rownames(govaere) <- govaere$sample_id
   govaere <- govaere[, -c(1,2)]
   
   metadata <- read_excel("data/Plasma_Proteomics_dataset/42255_2023_775_MOESM3_ESM.xlsx")[-c(1,2), -c(1:9)]
   metadata <- metadata[, colnames(govaere)]
   metadata <- as.double(metadata[1,])
   
   genes <- gene_map$external_gene_name
   
   # Subset + add missing
   govaere_sub <- govaere[intersect(genes, rownames(govaere)), , drop=FALSE]
   missing <- setdiff(genes, rownames(govaere_sub))
   
   if(length(missing) > 0){
      govaere_sub <- rbind(
         govaere_sub,
         matrix(0, length(missing), ncol(govaere_sub),
                dimnames = list(missing, colnames(govaere_sub)))
      )
   }
   
   # Z-score
   govaere_z <- t(scale(t(govaere_sub)))
   govaere_z[is.na(govaere_z)] <- 0
   
   df <- as.data.frame(t(govaere_z))
   df$fibrosis_bin <- factor(ifelse(metadata == 1, "Advanced", "Early"),
                             levels = c("Early","Advanced"))
   
   # --- SAFE RENAMING ---
   name2ensg <- setNames(gene_map$ensembl_gene_id,
                         gene_map$external_gene_name)
   
   gene_cols <- setdiff(colnames(df), "fibrosis_bin")
   mapped <- name2ensg[gene_cols]
   
   # Remove genes that failed mapping
   valid <- !is.na(mapped)
   df <- df[, c(gene_cols[valid], "fibrosis_bin")]
   colnames(df)[1:sum(valid)] <- mapped[valid]
   
   return(df)
}

### === 3. EVALUATE MODEL FUNCTION === ###
evaluate_model <- function(model, df, train_genes) {
   
   # Add missing genes
   missing <- setdiff(train_genes, colnames(df))
   if(length(missing) > 0){
      df[missing] <- 0
   }
   
   # Align order
   df <- df[, c(train_genes, "fibrosis_bin")]
   
   prob <- predict(model$finalModel, df, type="prob")[,"Advanced"]
   pred <- factor(ifelse(prob > 0.5, "Advanced", "Early"),
                  levels = c("Early","Advanced"))
   
   roc_obj <- roc(df$fibrosis_bin, prob)
   
   list(
      roc = roc_obj,
      auc = auc(roc_obj),
      prob = prob,
      pred = pred,
      cm = confusionMatrix(pred, df$fibrosis_bin)
   )
}

### === 4. RUN PIPELINE FOR BOTH MODELS AND 3-PUBLISHED GENES BMs === ###
expr_file <- "data/z_scores.csv"
template_file <- "data/metadata.csv"


# 57 genes
genes_57 <- read.csv("data/57_BMs.csv")[,-1]
rf_57 <- train_rf_model("data/57_BMs.csv", expr_file, template_file)
govaere_df_57 <- prepare_govaere(genes_57)
res_57 <- evaluate_model(rf_57, govaere_df_57,
                         setdiff(colnames(rf_57$trainingData),"fibrosis_bin"))

# 15 genes
genes_15 <- read.csv("data/top_15_BMs_from_ELBOW.csv")[,-1]
rf_15 <- train_rf_model("data/top_15_BMs_from_ELBOW.csv", expr_file, template_file)
govaere_df_15 <- prepare_govaere(genes_15)
res_15 <- evaluate_model(rf_15, govaere_df_15,
                         setdiff(colnames(rf_15$trainingData),"fibrosis_bin"))

# 3 published BMs
genes_3 <- read.csv("data/External_3_gene_panel_for_additional_validation.csv")[,-1]
rf_3 <- train_rf_model("data/External_3_gene_panel_for_additional_validation.csv", expr_file, template_file)
govaere_df_3 <- prepare_govaere(genes_3)
res_3 <- evaluate_model(rf_3, govaere_df_3,
                         setdiff(colnames(rf_3$trainingData),"fibrosis_bin"))

### === 5. COMBINED ROC PLOT === ###
# Smooth ROC objects
roc_57_smooth <- smooth(res_57$roc, n=1000)
roc_15_smooth <- smooth(res_15$roc, n=1000)
roc_3_smooth  <- smooth(res_3$roc,  n=1000)

df_all <- bind_rows(
   data.frame(FPR = 1 - roc_57_smooth$specificities,
              TPR = roc_57_smooth$sensitivities,
              Model = "57 genes"),
   data.frame(FPR = 1 - roc_15_smooth$specificities,
              TPR = roc_15_smooth$sensitivities,
              Model = "15 genes"),
   data.frame(FPR = 1 - roc_3_smooth$specificities,
              TPR = roc_3_smooth$sensitivities,
              Model = "3 genes")
)

df_all <- df_all %>%
   group_by(Model) %>%
   arrange(FPR, .by_group = TRUE)

combined_sens_spec_plot <- ggplot(df_all,
                                  aes(FPR, TPR, color=Model, group=Model)) +
   # 57 genes
   geom_line(data = df_all[df_all$Model=="57 genes", ], 
             aes(FPR, TPR, color=Model), size=1) +
   # 15 genes
   geom_line(data = df_all[df_all$Model=="15 genes", ], 
             aes(FPR, TPR, color=Model), size=1) +
   # 3 genes
   geom_line(data = df_all[df_all$Model=="3 genes", ], 
             aes(FPR, TPR, color=Model), size=0.5) +
   geom_abline(linetype="dashed") +
   scale_color_manual(values=c(
      "57 genes"="purple4",
      "15 genes"="mediumpurple",
      "3 genes"="black"
   )) +
   theme_minimal(base_size=14) +
   coord_fixed(ratio = 1)


ggsave("ROC_Govaere_Proteomics_RF_combined.pdf", 
       plot = combined_sens_spec_plot, width = 6, height = 6)

cat("AUC (57 genes):", round(as.numeric(res_57$auc), 2), "\n")
cat("AUC (15 genes):", round(as.numeric(res_15$auc), 2), "\n")
cat("AUC (3 genes):", round(as.numeric(res_3$auc), 2), "\n")


# DeLong statistical tests
delong_57_vs_15 <- roc.test(res_57$roc, res_15$roc, method="delong")
delong_57_vs_3  <- roc.test(res_57$roc, res_3$roc,  method="delong")
delong_15_vs_3  <- roc.test(res_15$roc, res_3$roc,  method="delong")

# Print p-values only
cat("DeLong p-value 57 vs 15:", delong_57_vs_15$p.value, "\n")
cat("DeLong p-value 57 vs 3: ", delong_57_vs_3$p.value, "\n")
cat("DeLong p-value 15 vs 3: ", delong_15_vs_3$p.value, "\n")


