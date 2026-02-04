

### GWAS analysis

library(dplyr)
library(readr)
library(biomaRt)
library(tidyr)



# Define trait category mapping
trait_categories <- list(
  Metabolic = c("BMI", "body mass", "cholesterol", "lipid", "obesity", "insulin", "glucose",
                "triglyceride", "metabolic", "diabetes", "adipose", "HDL", "LDL"),
  Cardiovascular = c("blood pressure", "hypertension", "heart", "coronary", "stroke", "cardio",
                     "atherosclerosis", "vascular"),
  Inflammatory = c("inflammation", "cytokine", "C-reactive", "immune", "autoimmune", "lupus",
                   "psoriasis", "Crohn", "colitis", "asthma"),
  Liver = c("liver", "steatosis", "fibrosis", "ALT", "AST", "bilirubin", "gamma-glutamyl",
            "transaminase", "fatty liver", "hepatitis"),
  Neurological = c("schizophrenia", "cognition", "autism", "Alzheimer", "Parkinson", "depression",
                   "neuro", "brain"),
  Cancer = c("cancer", "carcinoma", "tumor", "neoplasm", "melanoma", "leukemia")
)



assign_trait_category <- function(trait, categories) {
  assigned <- c()
  for (cat in names(categories)) {
    if (any(sapply(categories[[cat]], function(keyword)
      grepl(keyword, trait, ignore.case = TRUE)))) {
      assigned <- c(assigned, cat)
    }
  }
  if (length(assigned) == 0) assigned <- "Other"
  paste(unique(assigned), collapse = ", ")
}



gwas_summary <- read.csv("data/GWAS_catalog/gwas_summary.csv")

# Assign categories
gwas_summary$Category <- sapply(gwas_summary$Trait, assign_trait_category, categories = trait_categories)

# Count number of genes per category
category_summary <- gwas_summary %>%
  separate_rows(Category, sep = ", ") %>%
  group_by(Category) %>%
  summarise(n_genes = n_distinct(MAPPED_GENE)) %>%
  arrange(desc(n_genes))

write.csv(category_summary, "gwas_trait_categories.csv", row.names = FALSE)




library(ggplot2)

p <- ggplot(category_summary, aes(x = reorder(Category, -n_genes), y = n_genes, fill = "blue")) +
  geom_bar(stat = "identity") +
  labs(x = "Trait Category", y = "Number of Biomarker Genes",
       title = "GWAS Trait Categories for MASLD Biomarkers") +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_flip()

ggsave("GWAS_Trait_Categories.pdf", plot = p, width = 8, height = 6, units = "in")


# Summary percentages and coverage
total_genes <- length(unique(gwas_summary$MAPPED_GENE))
category_summary <- category_summary %>%
  mutate(Percentage = round(100 * n_genes / total_genes, 1))

print(category_summary)



#Multi-system biomarkers
multi_category_genes <- gwas_summary %>%
  mutate(n_categories = sapply(strsplit(Category, ", "), length)) %>%
  filter(n_categories > 1)

write.csv(multi_category_genes, "gwas_multi_category_genes.csv", row.names = FALSE)





# Visualise gene category relationships
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggtext)  # needed for element_markdown

# Load Top 15 gene names
top15 <- read.csv("data/top_15_BMs_from_ELBOW.csv")[,-1]
top15 <- top15$external_gene_name

# Mark top15 genes
gene_category_matrix <- gwas_summary %>%
  separate_rows(Category, sep = ", ") %>%
  distinct(MAPPED_GENE, Category) %>%
  mutate(
    Highlight = ifelse(MAPPED_GENE %in% top15, "Top 15 Biomarkers", "57 Biomarkers")
  )

# Order genes to control y-axis order
gene_order <- gene_category_matrix %>%
  distinct(MAPPED_GENE, Highlight) %>%
  arrange(desc(Highlight), MAPPED_GENE) %>%
  pull(MAPPED_GENE)

gene_category_matrix$MAPPED_GENE <- factor(gene_category_matrix$MAPPED_GENE, levels = rev(gene_order))

# Create labels: bold only for Top 15 genes
gene_labels <- sapply(levels(gene_category_matrix$MAPPED_GENE), function(g) {
  if (g %in% top15) {
    paste0("**", g, "**")  # Markdown bold
  } else {
    g
  }
})

# Color palette
cols <- c("57 Biomarkers" = "grey70",
          "Top 15 Biomarkers" = "purple3")

# Plot
p <- ggplot(gene_category_matrix, aes(x = Category, y = MAPPED_GENE)) +
  geom_point(aes(color = Highlight, size = Highlight)) +
  scale_color_manual(values = cols) +
  scale_size_manual(values = c("Top 15 Biomarkers" = 2, "57 Biomarkers" = 2)) +
  scale_y_discrete(labels = gene_labels) +  # apply bold labels
  theme_minimal() +
  labs(x = "Trait Category", y = "Biomarker Gene",
       title = "Gene–Trait Category Associations (GWAS Catalog)") +
  theme(
    axis.text.y = element_markdown(size = 7),  # render Markdown for bold
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

ggsave("Gene_category_matrix_boldTop15.pdf", plot = p, width = 8, height = 6, units = "in")








# Fisher enrichment per category
# required libraries
library(dplyr)
library(tidyr)
library(readr)

# Load your GWAS summary with trait categories
gwas_summary <- read.csv("data/GWAS_catalog/gwas_summary.csv")


gwas_summary$Category <- sapply(gwas_summary$Trait, assign_trait_category, categories = trait_categories)

# Ensure Category column exists
if(!"Category" %in% colnames(gwas_summary)) {
  stop("Please assign trait categories first (as in your previous script).")
}

# Split multi-category rows into multiple rows
gene_cat <- gwas_summary %>%
  dplyr::select(MAPPED_GENE, Category) %>%
  mutate(Category = ifelse(is.na(Category) | Category == "", "Other", Category)) %>%
  tidyr::separate_rows(Category, sep = ",\\s*") %>%
  distinct()




# Your 35 MASLD biomarker genes
masld_genes <- gwas_summary$MAPPED_GENE

# Background genes (all genes in GWAS summary)
background_genes <- unique(gwas_summary$MAPPED_GENE)




enrichment_results <- gene_cat %>%
  group_by(Category) %>%
  summarise(
    n_in_masld = sum(MAPPED_GENE %in% masld_genes),
    n_not_in_masld = sum(!MAPPED_GENE %in% masld_genes)
  ) %>%
  rowwise() %>%
  mutate(
    # Build 2x2 table: [MASLD, Not MASLD] x [In category, Not in category]
    fisher_test = list(fisher.test(matrix(
      c(n_in_masld,
        length(masld_genes) - n_in_masld,
        n_not_in_masld,
        length(background_genes) - length(masld_genes) - n_not_in_masld + n_in_masld),
      nrow = 2
    ))),
    OR = fisher_test$estimate,
    p_value = fisher_test$p.value
  ) %>%
  ungroup() %>%
  mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  dplyr::select(Category, OR, p_value, FDR) %>%
  arrange(FDR)


write.csv(enrichment_results, "gwas_category_enrichment.csv", row.names = FALSE)





library(ggplot2)

# Define significance categories
enrichment_results <- enrichment_results %>%
  mutate(
    sig_status = case_when(
      FDR < 0.05 ~ "FDR_significant",
      p_value < 0.05 & FDR >= 0.05 ~ "p_significant",
      TRUE ~ "Not_significant"
    ),
    logFDR = -log10(FDR)
  )

# Define custom colors
sig_colors <- c(
  "FDR_significant" = "#D55E00",
  "p_significant" = "#E69F00",
  "Not_significant" = "grey70"
)

# Plot
p <- ggplot(enrichment_results, aes(x = reorder(Category, logFDR), y = logFDR, fill = sig_status)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = sig_colors) +
  labs(
    x = "Trait Category",
    y = "-log10(FDR)",
    title = "Enrichment of MASLD Biomarker Genes Across GWAS Trait Categories",
    fill = "Significance"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

ggsave("Enrichment_Results.pdf", plot = p, width = 6, height = 5, units = "in")



