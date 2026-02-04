
corrected = read.csv("data/batch_corrected_counts_(dataset+gender).csv")
rownames(corrected) = corrected$X
corrected = corrected[,-1]
normalised_beforecorrection = read.csv("data/Human_Normalised_counts.csv")
rownames(normalised_beforecorrection) = normalised_beforecorrection$X
normalised_beforecorrection = normalised_beforecorrection[,-1]


colnames(corrected)
colnames(normalised_beforecorrection)
corrected = corrected[  intersect(rownames(normalised_beforecorrection),rownames(corrected)), intersect(colnames(normalised_beforecorrection),colnames(corrected))]
normalised_beforecorrection = normalised_beforecorrection[ intersect(rownames(normalised_beforecorrection),rownames(corrected)) , intersect(colnames(normalised_beforecorrection),colnames(corrected))]

codes = read.csv('data/metadata.csv')[,-1]
codes = codes[ -5 , ]
codes[,1] = colnames(corrected)

normalised_SANYAL_CTRL = rowMeans(subset(normalised_beforecorrection, select = which(codes$Dataset == "SANYAL" & codes$SAF.score == "CTRL"), na.rm = TRUE))
corrected_SANYAL_CTRL = rowMeans(subset(corrected, select = which(codes$Dataset == "SANYAL" & codes$SAF.score == "CTRL"), na.rm = TRUE))

normalised_SANYAL_MASL = rowMeans(subset(normalised_beforecorrection, select = which(codes$Dataset == "SANYAL" & codes$SAF.score == "MASL"), na.rm = TRUE))
corrected_SANYAL_MASL = rowMeans(subset(corrected, select = which(codes$Dataset == "SANYAL" & codes$SAF.score == "MASL"), na.rm = TRUE))

normalised_SANYAL_MASH_F0 = rowMeans(subset(normalised_beforecorrection, select = which(codes$Dataset == "SANYAL" & codes$SAF.score == "MASH F0"), na.rm = TRUE))
corrected_SANYAL_MASH_F0 = rowMeans(subset(corrected, select = which(codes$Dataset == "SANYAL" & codes$SAF.score == "MASH F0"), na.rm = TRUE))

normalised_SANYAL_MASH_F1 = rowMeans(subset(normalised_beforecorrection, select = which(codes$Dataset == "SANYAL" & codes$SAF.score == "MASH F1"), na.rm = TRUE))
corrected_SANYAL_MASH_F1 = rowMeans(subset(corrected, select = which(codes$Dataset == "SANYAL" & codes$SAF.score == "MASH F1"), na.rm = TRUE))

normalised_SANYAL_MASH_F2 = rowMeans(subset(normalised_beforecorrection, select = which(codes$Dataset == "SANYAL" & codes$SAF.score == "MASH F2"), na.rm = TRUE))
corrected_SANYAL_MASH_F2 = rowMeans(subset(corrected, select = which(codes$Dataset == "SANYAL" & codes$SAF.score == "MASH F2"), na.rm = TRUE))

normalised_SANYAL_MASH_F34 = rowMeans(subset(normalised_beforecorrection, select = which(codes$Dataset == "SANYAL" & (codes$SAF.score == "MASH F3" | codes$SAF.score == "MASH F4")), na.rm = TRUE))
corrected_SANYAL_MASH_F34 = rowMeans(subset(corrected, select = which(codes$Dataset == "SANYAL" & (codes$SAF.score == "MASH F3" | codes$SAF.score == "MASH F4")), na.rm = TRUE))




#plot meanNAS in each SW
library("ggpubr")
library(cowplot)
pdf("before-after correction comparison SANYAL.pdf")
mydata = data.frame("normalised_SANYAL_CTRL" = normalised_SANYAL_CTRL, "corrected_SANYAL_CTRL" = corrected_SANYAL_CTRL)
CTRL_SANYAL <- ggscatter(mydata, x = "normalised_SANYAL_CTRL", y = "corrected_SANYAL_CTRL",
                         color = "black", shape = 20, size = 1, alpha = 0.3, # Points color, shape, size and transparency
                         add = "reg.line",  # Add regressin line
                         add.params = list(color = "black", fill = "lightgray"), conf.int = TRUE, 
                         cor.coef = TRUE, cor.method = "pearson",
                         #ylim = c(1,6),
                         xlab = "Normalised CTRL\n(SANYAL)", ylab = "Corrected CTRL\n(SANYAL)", main = "CTRL SANYAL")

mydata = data.frame("normalised_SANYAL_MASL" = normalised_SANYAL_MASL, "corrected_SANYAL_MASL" = corrected_SANYAL_MASL)
MASL_SANYAL <- ggscatter(mydata, x = "normalised_SANYAL_MASL", y = "corrected_SANYAL_MASL",
                         color = "#f3acac", shape = 20, size = 1, alpha = 0.3, # Points color, shape, size and transparency
                         add = "reg.line",  # Add regressin line
                         add.params = list(color = "black", fill = "lightgray"), conf.int = TRUE, 
                         cor.coef = TRUE, cor.method = "pearson",
                         #ylim = c(1,6),
                         xlab = "Normalised MASL\n(SANYAL)", ylab = "Corrected MASL\n(SANYAL)", main = "MASL SANYAL")

mydata = data.frame("normalised_SANYAL_MASH_F0" = normalised_SANYAL_MASH_F0, "corrected_SANYAL_MASH_F0" = corrected_SANYAL_MASH_F0)
MASH_F0_SANYAL <- ggscatter(mydata, x = "normalised_SANYAL_MASH_F0", y = "corrected_SANYAL_MASH_F0",
                         color = "#ea6564", shape = 20, size = 1, alpha = 0.3, # Points color, shape, size and transparency
                         add = "reg.line",  # Add regressin line
                         add.params = list(color = "black", fill = "lightgray"), conf.int = TRUE, 
                         cor.coef = TRUE, cor.method = "pearson",
                         #ylim = c(1,6),
                         xlab = "Normalised MASHF0\n(SANYAL)", ylab = "Corrected MASHF0\n(SANYAL)", main = "MASHF0 SANYAL")

mydata = data.frame("normalised_SANYAL_MASH_F1" = normalised_SANYAL_MASH_F1, "corrected_SANYAL_MASH_F1" = corrected_SANYAL_MASH_F1)
MASH_F1_SANYAL <- ggscatter(mydata, x = "normalised_SANYAL_MASH_F1", y = "corrected_SANYAL_MASH_F1",
                            color = "#de1e1d", shape = 20, size = 1, alpha = 0.3, # Points color, shape, size and transparency
                            add = "reg.line",  # Add regressin line
                            add.params = list(color = "black", fill = "lightgray"), conf.int = TRUE, 
                            cor.coef = TRUE, cor.method = "pearson",
                            #ylim = c(1,6),
                            xlab = "Normalised MASHF1\n(SANYAL)", ylab = "Corrected MASHF1\n(SANYAL)", main = "MASHF1 SANYAL")

mydata = data.frame("normalised_SANYAL_MASH_F2" = normalised_SANYAL_MASH_F2, "corrected_SANYAL_MASH_F2" = corrected_SANYAL_MASH_F2)
MASH_F2_SANYAL <- ggscatter(mydata, x = "normalised_SANYAL_MASH_F2", y = "corrected_SANYAL_MASH_F2",
                            color = "#971414", shape = 20, size = 1, alpha = 0.3, # Points color, shape, size and transparency
                            add = "reg.line",  # Add regressin line
                            add.params = list(color = "black", fill = "lightgray"), conf.int = TRUE, 
                            cor.coef = TRUE, cor.method = "pearson",
                            #ylim = c(1,6),
                            xlab = "Normalised MASHF2\n(SANYAL)", ylab = "Corrected MASHF2\n(SANYAL)", main = "MASHF2 SANYAL")

mydata = data.frame("normalised_SANYAL_MASH_F34" = normalised_SANYAL_MASH_F34, "corrected_SANYAL_MASH_F34" = corrected_SANYAL_MASH_F34)
MASH_F34_SANYAL <- ggscatter(mydata, x = "normalised_SANYAL_MASH_F34", y = "corrected_SANYAL_MASH_F34",
                            color = "#4f0a0a", shape = 20, size = 1, alpha = 0.3, # Points color, shape, size and transparency
                            add = "reg.line",  # Add regressin line
                            add.params = list(color = "black", fill = "lightgray"), conf.int = TRUE, 
                            cor.coef = TRUE, cor.method = "pearson",
                            #ylim = c(1,6),
                            xlab = "Normalised MASHF34\n(SANYAL)", ylab = "Corrected MASHF34\n(SANYAL)", main = "MASHF34 SANYAL")

plot_grid(CTRL_SANYAL, MASL_SANYAL, MASH_F0_SANYAL, MASH_F1_SANYAL, MASH_F2_SANYAL, MASH_F34_SANYAL, nrow = 3, ncol = 2, align = "hv")

dev.off()













#same (similar) for UCAM
normalised_UCAM_MASL = rowMeans(subset(normalised_beforecorrection, select = which(codes$Dataset == "UCAM" & codes$SAF.score == "MASL"), na.rm = TRUE))
corrected_UCAM_MASL = rowMeans(subset(corrected, select = which(codes$Dataset == "UCAM" & codes$SAF.score == "MASL"), na.rm = TRUE))

normalised_UCAM_MASH_F0 = rowMeans(subset(normalised_beforecorrection, select = which(codes$Dataset == "UCAM" & codes$SAF.score == "MASH F0"), na.rm = TRUE))
corrected_UCAM_MASH_F0 = rowMeans(subset(corrected, select = which(codes$Dataset == "UCAM" & codes$SAF.score == "MASH F0"), na.rm = TRUE))

normalised_UCAM_MASH_F1 = rowMeans(subset(normalised_beforecorrection, select = which(codes$Dataset == "UCAM" & codes$SAF.score == "MASH F1"), na.rm = TRUE))
corrected_UCAM_MASH_F1 = rowMeans(subset(corrected, select = which(codes$Dataset == "UCAM" & codes$SAF.score == "MASH F1"), na.rm = TRUE))

normalised_UCAM_MASH_F2 = rowMeans(subset(normalised_beforecorrection, select = which(codes$Dataset == "UCAM" & codes$SAF.score == "MASH F2"), na.rm = TRUE))
corrected_UCAM_MASH_F2 = rowMeans(subset(corrected, select = which(codes$Dataset == "UCAM" & codes$SAF.score == "MASH F2"), na.rm = TRUE))

normalised_UCAM_MASH_F3 = rowMeans(subset(normalised_beforecorrection, select = which(codes$Dataset == "UCAM" & codes$SAF.score == "MASH F3"), na.rm = TRUE))
corrected_UCAM_MASH_F3 = rowMeans(subset(corrected, select = which(codes$Dataset == "UCAM" & codes$SAF.score == "MASH F3"), na.rm = TRUE))

normalised_UCAM_MASH_F4 = rowMeans(subset(normalised_beforecorrection, select = which(codes$Dataset == "UCAM" & codes$SAF.score == "MASH F4"), na.rm = TRUE))
corrected_UCAM_MASH_F4 = rowMeans(subset(corrected, select = which(codes$Dataset == "UCAM" & codes$SAF.score == "MASH F4"), na.rm = TRUE))



#plot meanNAS in each SW
library("ggpubr")
library(cowplot)
pdf("before-after correction comparison UCAM.pdf")

mydata = data.frame("normalised_UCAM_MASL" = normalised_UCAM_MASL, "corrected_UCAM_MASL" = corrected_UCAM_MASL)
MASL_UCAM <- ggscatter(mydata, x = "normalised_UCAM_MASL", y = "corrected_UCAM_MASL",
                         color = "#f3acac", shape = 20, size = 1, alpha = 0.3, # Points color, shape, size and transparency
                         add = "reg.line",  # Add regressin line
                         add.params = list(color = "black", fill = "lightgray"), conf.int = TRUE, 
                         cor.coef = TRUE, cor.method = "pearson",
                         #ylim = c(1,6),
                         xlab = "Normalised MASL\n(UCAM)", ylab = "Corrected MASL\n(UCAM)", main = "MASL UCAM")

mydata = data.frame("normalised_UCAM_MASH_F0" = normalised_UCAM_MASH_F0, "corrected_UCAM_MASH_F0" = corrected_UCAM_MASH_F0)

mydata = data.frame("normalised_UCAM_MASH_F1" = normalised_UCAM_MASH_F1, "corrected_UCAM_MASH_F1" = corrected_UCAM_MASH_F1)
MASH_F1_UCAM <- ggscatter(mydata, x = "normalised_UCAM_MASH_F1", y = "corrected_UCAM_MASH_F1",
                            color = "#de1e1d", shape = 20, size = 1, alpha = 0.3, # Points color, shape, size and transparency
                            add = "reg.line",  # Add regressin line
                            add.params = list(color = "black", fill = "lightgray"), conf.int = TRUE, 
                            cor.coef = TRUE, cor.method = "pearson",
                            #ylim = c(1,6),
                            xlab = "Normalised MASHF1\n(UCAM)", ylab = "Corrected MASHF1\n(UCAM)", main = "MASHF1 UCAM")

mydata = data.frame("normalised_UCAM_MASH_F2" = normalised_UCAM_MASH_F2, "corrected_UCAM_MASH_F2" = corrected_UCAM_MASH_F2)
MASH_F2_UCAM <- ggscatter(mydata, x = "normalised_UCAM_MASH_F2", y = "corrected_UCAM_MASH_F2",
                            color = "#971414", shape = 20, size = 1, alpha = 0.3, # Points color, shape, size and transparency
                            add = "reg.line",  # Add regressin line
                            add.params = list(color = "black", fill = "lightgray"), conf.int = TRUE, 
                            cor.coef = TRUE, cor.method = "pearson",
                            #ylim = c(1,6),
                            xlab = "Normalised MASHF2\n(UCAM)", ylab = "Corrected MASHF2\n(UCAM)", main = "MASHF2 UCAM")

mydata = data.frame("normalised_UCAM_MASH_F3" = normalised_UCAM_MASH_F3, "corrected_UCAM_MASH_F3" = corrected_UCAM_MASH_F3)
MASH_F3_UCAM <- ggscatter(mydata, x = "normalised_UCAM_MASH_F3", y = "corrected_UCAM_MASH_F3",
                             color = "#4f0a0a", shape = 20, size = 1, alpha = 0.3, # Points color, shape, size and transparency
                             add = "reg.line",  # Add regressin line
                             add.params = list(color = "black", fill = "lightgray"), conf.int = TRUE, 
                             cor.coef = TRUE, cor.method = "pearson",
                             #ylim = c(1,6),
                             xlab = "Normalised MASHF3\n(UCAM)", ylab = "Corrected MASHF3\n(UCAM)", main = "MASHF3 UCAM")

mydata = data.frame("normalised_UCAM_MASH_F4" = normalised_UCAM_MASH_F4, "corrected_UCAM_MASH_F4" = corrected_UCAM_MASH_F4)
MASH_F4_UCAM <- ggscatter(mydata, x = "normalised_UCAM_MASH_F4", y = "corrected_UCAM_MASH_F4",
                          color = "#400808", shape = 20, size = 1, alpha = 0.3, # Points color, shape, size and transparency
                          add = "reg.line",  # Add regressin line
                          add.params = list(color = "black", fill = "lightgray"), conf.int = TRUE, 
                          cor.coef = TRUE, cor.method = "pearson",
                          #ylim = c(1,6),
                          xlab = "Normalised MASHF4\n(UCAM)", ylab = "Corrected MASHF4\n(UCAM)", main = "MASHF4 UCAM")

plot_grid(MASL_UCAM, MASH_F1_UCAM, MASH_F2_UCAM, MASH_F3_UCAM, MASH_F4_UCAM, nrow = 3, ncol = 2, align = "hv")

dev.off()












# --- Stratified correlation analysis ---
library(ggplot2)
library(dplyr)

# Compute mean expression across all samples (before correction)
mean_expr <- rowMeans(normalised_beforecorrection, na.rm = TRUE)

# Divide genes into three expression bins (low, medium, high)
expr_bin <- cut(mean_expr,
                breaks = quantile(mean_expr, probs = c(0, 1/3, 2/3, 1)),
                labels = c("Low", "Medium", "High"),
                include.lowest = TRUE)

# Compute per-gene correlation between before and after correction
# Using fast vectorized row-wise correlation
rowCor <- function(x, y) {
  xm <- rowMeans(x)
  ym <- rowMeans(y)
  num <- rowSums((x - xm) * (y - ym))
  den <- sqrt(rowSums((x - xm)^2) * rowSums((y - ym)^2))
  num / den
}

cor_values <- rowCor(as.matrix(normalised_beforecorrection), as.matrix(corrected))

# Combine results
df_cor <- data.frame(
  gene = names(mean_expr),
  mean_expr = mean_expr,
  expr_bin = expr_bin,
  cor = cor_values
)

# Summarise mean correlation per expression bin
summary_df <- df_cor %>%
  group_by(expr_bin) %>%
  summarise(mean_correlation = mean(cor, na.rm = TRUE))

# Plot
pdf("Expression_Stratified_Correlation.pdf", width = 2.5, height = 2)
ggplot(summary_df, aes(x = expr_bin, y = mean_correlation, fill = expr_bin)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  scale_fill_manual(values = c("#B8E186", "#4DAC26", "#0571B0")) +
  labs(x = "Expression level", y = "Mean Pearson correlation") +
  coord_cartesian(ylim = c(0, 1)) +   # <-- y-axis from 0.4 to 1
  theme_minimal(base_size = 9) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9)
  )
dev.off()


# Print numeric summary
print(summary_df)

















# --- Stratified correlation per disease stage ---
library(dplyr)
library(ggplot2)

# Define a helper function to compute row-wise Pearson correlations quickly
rowCor <- function(x, y) {
  xm <- rowMeans(x)
  ym <- rowMeans(y)
  num <- rowSums((x - xm) * (y - ym))
  den <- sqrt(rowSums((x - xm)^2) * rowSums((y - ym)^2))
  num / den
}

# Prepare metadata
codes <- codes %>%
  mutate(SAF.score = factor(SAF.score, levels = c("CTRL", "MASL", "MASH F0", "MASH F1", "MASH F2", "MASH F3", "MASH F4")))

# Prepare expression bins (based on all samples)
mean_expr <- rowMeans(normalised_beforecorrection, na.rm = TRUE)
expr_bin <- cut(mean_expr,
                breaks = quantile(mean_expr, probs = c(0, 1/3, 2/3, 1)),
                labels = c("Low", "Medium", "High"),
                include.lowest = TRUE)

# Create list of disease stages
stages <- levels(codes$SAF.score)
stages <- stages[!is.na(stages)]

# Initialize results storage
results_list <- list()

# Loop through disease stages
for (stage in stages) {
  samples_stage <- codes$Sample.name[codes$SAF.score == stage]
  
  if (length(samples_stage) < 3) next  # skip if not enough samples
  
  expr_pre_stage <- normalised_beforecorrection[, colnames(normalised_beforecorrection) %in% samples_stage, drop = FALSE]
  expr_post_stage <- corrected[, colnames(corrected) %in% samples_stage, drop = FALSE]
  
  cor_stage <- rowCor(as.matrix(expr_pre_stage), as.matrix(expr_post_stage))
  
  results_list[[stage]] <- data.frame(
    gene = rownames(expr_pre_stage),
    expr_bin = expr_bin,
    cor = cor_stage,
    stage = stage
  )
}

# Combine all into one data frame
df_stage_cor <- do.call(rbind, results_list)

# Compute mean correlation per bin and stage
summary_stage_df <- df_stage_cor %>%
  group_by(stage, expr_bin) %>%
  summarise(mean_correlation = mean(cor, na.rm = TRUE), .groups = "drop")

# --- Plot grouped barplot ---
pdf("Expression_Stratified_Correlation_by_Stage.pdf", width = 4.5, height = 3)
ggplot(summary_stage_df, aes(x = expr_bin, y = mean_correlation, fill = expr_bin)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  facet_wrap(~stage, nrow = 2) +
  scale_fill_manual(values = c("#B8E186", "#4DAC26", "#0571B0")) +
  labs(x = "Expression level", y = "Mean Pearson correlation") +
  coord_cartesian(ylim = c(0.4, 1)) +
  theme_minimal(base_size = 9) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    strip.text = element_text(size = 8, face = "bold")
  )
dev.off()

# Print numeric summary
print(summary_stage_df)




# --- Overall correlation distribution across all stages ---
pdf("Expression_Stratified_Correlation_AllStages.pdf", width = 2.5, height = 2)

ggplot(df_stage_cor, aes(x = expr_bin, y = cor, fill = expr_bin)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6, color = "black") +
  scale_fill_manual(values = c("#B8E186", "#4DAC26", "#0571B0")) +
  labs(x = "Expression level", y = "Pearson correlation") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 9) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9)
  )

dev.off()






library(ggplot2)
library(dplyr)

# --- Compute summary statistics ---
# df_stage_cor must contain: expr_bin, stage, cor
summary_stage_df <- df_stage_cor %>%
  group_by(expr_bin, stage) %>%
  summarise(mean_correlation = mean(cor, na.rm = TRUE)) %>%
  ungroup()

# overall mean per expr_bin (for the bars)
overall_means <- summary_stage_df %>%
  group_by(expr_bin) %>%
  summarise(overall_mean = mean(mean_correlation))

# --- Barplot with stage means as dots ---
pdf("Expression_Stratified_Correlation_AllStages.pdf", width = 2.5, height = 2)

ggplot() +
  # Bars = overall mean correlation per expression bin
  geom_bar(data = overall_means, aes(x = expr_bin, y = overall_mean, fill = expr_bin),
           stat = "identity", color = "black", width = 0.6) +
  # Points = mean per disease stage within each expression bin
  geom_point(data = summary_stage_df, 
             aes(x = expr_bin, y = mean_correlation, color = stage),
             position = position_jitter(width = 0.05, height = 0), 
             size = .3) +
  scale_fill_manual(values = c("#B8E186", "#4DAC26", "#0571B0")) +
  scale_color_manual(values = c(
    "CTRL" = "black",
    "MASL" = "#f3acac",
    "MASH F0" = "#ea6564",
    "MASH F1" = "#de1e1d",
    "MASH F2" = "#971414",
    "MASH F3" = "#4f0a0a",
    "MASH F4" = "#4f0a0a"
  )) +
  labs(x = "Expression level", y = "Mean Pearson correlation") +
  coord_cartesian(ylim = c(0.4, 1)) +
  theme_minimal(base_size = 9) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9)
  )

dev.off()











# Add standard deviation for the error bars
overall_means <- summary_stage_df %>%
  group_by(expr_bin) %>%
  summarise(
    overall_mean = mean(mean_correlation),
    sd = sd(mean_correlation),
    se = sd / sqrt(n())
  )

pdf("Expression_Stratified_Correlation_AllStages.pdf", width = 2.5, height = 2)

ggplot() +
  # Bars = overall mean correlation per expression bin
  geom_bar(data = overall_means, aes(x = expr_bin, y = overall_mean, fill = expr_bin),
           stat = "identity", color = "black", width = 0.6) +
  # Error bars showing variability
  geom_errorbar(data = overall_means, 
                aes(x = expr_bin, ymin = overall_mean - sd, ymax = overall_mean + sd),
                width = 0.15, color = "black", linewidth = 0.3) +
  # Dots per disease stage
  geom_point(data = summary_stage_df, 
             aes(x = expr_bin, y = mean_correlation, color = stage),
             position = position_jitter(width = 0.05, height = 0), 
             size = 0.3) +
  scale_fill_manual(values = c("#B8E186", "#4DAC26", "#0571B0")) +
  scale_color_manual(values = c(
    "CTRL" = "black",
    "MASL" = "#f3acac",
    "MASH F0" = "#ea6564",
    "MASH F1" = "#de1e1d",
    "MASH F2" = "#971414",
    "MASH F3" = "#4f0a0a",
    "MASH F4" = "#4f0a0a"
  )) +
  labs(x = "Expression level", y = "Mean Pearson correlation") +
  coord_cartesian(ylim = c(0.4, 1)) +
  theme_minimal(base_size = 8) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10)
  )

dev.off()

