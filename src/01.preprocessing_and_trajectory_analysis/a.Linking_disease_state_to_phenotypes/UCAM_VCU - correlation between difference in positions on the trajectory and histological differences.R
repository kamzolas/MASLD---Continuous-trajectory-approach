# Load data
template <- read.csv(file = "data/metadata.csv",
                     header = TRUE, stringsAsFactors = FALSE, row.names = 1)
template <- template[template$Sample.name != "Sample 5", ]

sorted_samples <- read.csv("data/PC1_sorted_samples.csv")

# Merge on Sample name
merged_data <- merge(template, sorted_samples, by.x = "Sample.name", by.y = "X")
merged_data$PC1 <- -merged_data$PC1

# Get only the relevant columns
df <- merged_data[, c("Sample.name", "PC1", "NAS")]

# Generate all pairwise combinations
combinations <- expand.grid(1:nrow(df), 1:nrow(df))
combinations <- combinations[combinations$Var1 < combinations$Var2, ]  # Remove self-pairs and duplicates

# Compute pairwise differences
pairwise_df <- data.frame(
  Sample1 = df$Sample.name[combinations$Var1],
  Sample2 = df$Sample.name[combinations$Var2],
  PC1_diff = abs(df$PC1[combinations$Var1] - df$PC1[combinations$Var2]),
  NAS_diff = abs(df$NAS[combinations$Var1] - df$NAS[combinations$Var2])
)

library(ggplot2)
library(gridExtra)

# Function to make individual plots
make_plot <- function(y, y_label) {
  cor_res <- cor.test(pairwise_df$PC1_diff, pairwise_df[[y]], method = "pearson")
  
  # Format p-value for very small values
  p_value_formatted <- ifelse(cor_res$p.value < 1e-300, "< 1e-300", 
                              formatC(cor_res$p.value, format = "e", digits = 2))
  
  ggplot(pairwise_df, aes(x = PC1_diff, y = !!as.name(y))) +
    geom_point(alpha = 0.3, color = "darkblue") +
    geom_smooth(method = "lm", se = TRUE, color = "orange") +
    labs(
      title = paste0(y_label, "\nR = ", round(cor_res$estimate, 2),
                     ", p = ", p_value_formatted),
      x = "Pairwise difference on trajectory",
      y = paste0(y_label, " Difference")
    ) +
    theme_minimal(base_size = 12)
}

# Compute pairwise differences for all histology scores
# Signed differences (Var2 - Var1)
pairwise_df$PC1_diff <- merged_data$PC1[combinations$Var2] - merged_data$PC1[combinations$Var1]
pairwise_df$NAS_diff <- merged_data$NAS[combinations$Var2] - merged_data$NAS[combinations$Var1]
pairwise_df$Steatosis_diff <- merged_data$Fat[combinations$Var2] - merged_data$Fat[combinations$Var1]
pairwise_df$Inflammation_diff <- merged_data$Inflammation[combinations$Var2] - merged_data$Inflammation[combinations$Var1]
pairwise_df$Ballooning_diff <- merged_data$Ballooning[combinations$Var2] - merged_data$Ballooning[combinations$Var1]
pairwise_df$Fibrosis_diff <- merged_data$Fibrosis[combinations$Var2] - merged_data$Fibrosis[combinations$Var1]

# Create all plots
p1 <- make_plot("NAS_diff", "NAS")
p2 <- make_plot("Steatosis_diff", "Steatosis")
p3 <- make_plot("Inflammation_diff", "Inflammation")
p4 <- make_plot("Ballooning_diff", "Ballooning")
p5 <- make_plot("Fibrosis_diff", "Fibrosis")

# Combine into a multi-panel figure (2 rows)
grid.arrange(p1, p2, p3, p4, p5, ncol = 2)
