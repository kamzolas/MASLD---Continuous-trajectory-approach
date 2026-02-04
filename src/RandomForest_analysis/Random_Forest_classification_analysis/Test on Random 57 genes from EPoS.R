# === LIBRARIES === #
library(randomForest)
library(caret)
library(ggplot2)
library(pROC)
library(readxl)
library(dplyr)
library(tidyr)


# === PARAMETERS === #
n_random_sets <- 300
set.seed(42)

# === PATHS === #
template_path <- "data/metadata.csv"
zscores_path <- "data/z_scores.csv"
selected_genes_path <- "data/57_BMs.csv"
epos_meta_path <- "data/EPoS_dataset/epos_metadata.csv"
epos_expr_path <- "data/EPoS_dataset/Corrected_NormalisedCounts_EPOS2022.csv"


# === LOAD INTERNAL DATA === #
template <- read.csv(template_path)[-5, -c(1,3)]
corrected <- read.csv(zscores_path, check.names = FALSE)
rownames(corrected) <- corrected[,1]
corrected <- corrected[,-1]

selected_genes <- read.csv(selected_genes_path)[,-1]$ensembl_gene_id

gov <- read_excel("data/Plasma_Proteomics_dataset/Govaeres Top_Proteomics-Transcriptomics.xlsx")
gov <- gov$ENSID_TOP
top145 = read.csv("data/145 (from200) most important genes (variable importance with MeanDecreaseAccuracy>0).csv")[,-1]
top145 = top145$Sample.1

selected_genes = union(selected_genes, top145)
selected_genes = selected_genes[ !duplicated(selected_genes)]
selected_genes = intersect(selected_genes, gov)



expr <- corrected[selected_genes, ]
template <- template[match(colnames(expr), template$Sample.name), ]
stopifnot(all(template$Sample.name == colnames(expr)))

# Transpose and label
data_train <- as.data.frame(t(expr))
data_train$fibrosis_bin <- factor(ifelse(as.numeric(as.character(template$Fibrosis)) >= 3, "Advanced", "Early"),
                                  levels = c("Early", "Advanced"))


# === RF TRAINING === #
ctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE,
                     summaryFunction = twoClassSummary, savePredictions = "final")
rf_cv <- train(fibrosis_bin ~ ., data = data_train, method = "rf",
               trControl = ctrl, metric = "ROC")

# === CV PERFORMANCE === #
cv_conf <- confusionMatrix(rf_cv$pred$pred, rf_cv$pred$obs)
roc_cv <- roc(rf_cv$pred$obs, rf_cv$pred$Advanced)
auc_cv <- auc(roc_cv)
cat("AUC (CV) =", auc_cv, "\n")


# === LOAD EPOS DATA === #
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

z_epos$fibrosis_bin <- factor(ifelse(epos_meta$Fibrosis.stage >= 3, "Advanced", "Early"),
                              levels = c("Early", "Advanced"))

expected_genes <- setdiff(colnames(data_train), "fibrosis_bin")
missing <- setdiff(expected_genes, colnames(z_epos))
for (g in missing) z_epos[[g]] <- NA
z_epos <- z_epos[, c(expected_genes, "fibrosis_bin")]
means_train <- colMeans(data_train[, expected_genes])
for (g in expected_genes) {
   if (anyNA(z_epos[[g]])) z_epos[[g]][is.na(z_epos[[g]])] <- means_train[g]
}

# === PREDICT ON EPOS === #
pred_epos <- predict(rf_cv$finalModel, newdata = z_epos)
prob_epos <- predict(rf_cv$finalModel, newdata = z_epos, type = "prob")[,"Advanced"]
conf_epos <- confusionMatrix(pred_epos, z_epos$fibrosis_bin)
roc_epos <- roc(z_epos$fibrosis_bin, prob_epos)
auc_epos <- auc(roc_epos)
cat("EPoS AUC =", auc_epos, "\n")

# === RANDOM GENE SIGNATURE TESTING === #
random_auc_results <- data.frame(
   Set = character(n_random_sets + 1),
   AUC = numeric(n_random_sets + 1),
   Type = character(n_random_sets + 1),
   stringsAsFactors = FALSE
)
random_auc_results[1, ] <- c("Selected_57", auc_epos, "Original")

all_genes <- rownames(corrected)
all_genes <- all_genes[all_genes %in% rownames(corrected)]

for (i in 1:n_random_sets) {
   cat("Running random gene set", i, "\n")
   rand_genes <- sample(setdiff(all_genes, selected_genes), 57)
   
   rand_expr_train <- corrected[rand_genes, ]
   data_rand_train <- as.data.frame(t(rand_expr_train))
   data_rand_train$fibrosis_bin <- data_train$fibrosis_bin
   
   rf_rand <- train(fibrosis_bin ~ ., data = data_rand_train,
                    method = "rf", trControl = ctrl, metric = "ROC")
   
   expr_epos_rand <- epos[rand_genes, ]
   z_epos_rand <- t(scale(t(as.matrix(expr_epos_rand))))
   z_epos_rand <- as.data.frame(t(z_epos_rand))
   z_epos_rand <- z_epos_rand[common_samples, ]
   z_epos_rand$fibrosis_bin <- z_epos$fibrosis_bin
   
   missing <- setdiff(rand_genes, colnames(z_epos_rand))
   for (g in missing) z_epos_rand[[g]] <- NA
   z_epos_rand <- z_epos_rand[, c(rand_genes, "fibrosis_bin")]
   train_means <- colMeans(data_rand_train[, rand_genes], na.rm = TRUE)
   for (g in rand_genes) {
      if (anyNA(z_epos_rand[[g]])) {
         z_epos_rand[[g]][is.na(z_epos_rand[[g]])] <- train_means[g]
      }
   }
   
   prob_rand <- predict(rf_rand$finalModel, newdata = z_epos_rand, type = "prob")[,"Advanced"]
   roc_rand <- roc(z_epos_rand$fibrosis_bin, prob_rand)
   auc_rand <- auc(roc_rand)
   
   random_auc_results[i + 1, ] <- c(paste0("Random_", i), auc_rand, "Random")
}



# === PLOT AUC Violin Plots === #
library(ggplot2)
library(dplyr)

# Ensure proper formatting
random_auc_results$AUC <- as.numeric(random_auc_results$AUC)
random_auc_results$Type <- factor(random_auc_results$Type, levels = c("Random", "Original"))

# === EMPIRICAL P-VALUE CALCULATION === #
original_auc <- as.numeric(random_auc_results$AUC[random_auc_results$Type == "Original"])
random_aucs <- as.numeric(random_auc_results$AUC[random_auc_results$Type == "Random"])
pval_auc <- mean(random_aucs >= original_auc)

# Labels
auc_label <- paste0("AUC (original)")
pval_label <- paste0("P-value ", format.pval(pval_auc, digits = 2, eps = 0.1^64))

# === PLOT === #
p <- ggplot(random_auc_results, aes(x = Type, y = AUC, fill = Type)) +
   geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
   geom_jitter(data = random_auc_results %>% filter(Type == "Random"),
               aes(x = Type, y = AUC),
               width = 0.15, alpha = 0.4, size = 1, color = "gray40") +
   geom_point(data = random_auc_results %>% filter(Type == "Original"),
              aes(x = Type, y = AUC),
              shape = 18, color = "red", size = 4) +
   # AUC label next to red dot
   geom_text(data = random_auc_results %>% filter(Type == "Original"),
             aes(x = 2, y = original_auc, label = auc_label),
             hjust = -0.2, vjust = 0.5, size = 4.5, fontface = "italic") +
   # p-value label at top center
   annotate("text", x = 1.5, y = max(random_auc_results$AUC) + 0.01,
            label = pval_label, size = 5, fontface = "bold") +
   scale_fill_manual(values = c("gray80", "darkblue")) +
   labs(title = paste0("AUC Distribution from ", n_random_sets, " Random Same-Size Gene Sets"),
        y = "AUC on EPoS", x = NULL) +
   theme_minimal(base_size = 14) +
   theme(legend.position = "none",
         axis.text.x = element_text(size = 12, face = "bold"),
         plot.title = element_text(size = 15, face = "bold", hjust = 0.5))

# Save to PDF
ggsave("AUC_violin_plot.pdf", plot = p, width = 6.5, height = 5, dpi = 300)


saveRDS(random_auc_results, "random_for_57_genes_EPoS.RDS")
