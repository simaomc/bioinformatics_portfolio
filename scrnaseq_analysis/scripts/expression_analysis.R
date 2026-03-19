############################################################
# Single-cell RNA-seq Expression Analysis
# Author: Maryanna Simao
# Location: São Paulo, Brazil
# Description: Analysis of gene expression 
# in Drosophila spermatogenesis
############################################################
############################################################
# 📦 1. PACKAGES
############################################################

library(dplyr)
library(ggplot2)
library(patchwork)
library(pheatmap)

############################################################
# 📂 2. LOAD DATA
############################################################

df <- read.csv("../data/avgexp.txt", skip = 5)
df <- df[, 1:10]

cat("Total genes:", nrow(df), "\n")

############################################################
# 📊 3. EXPRESSION THRESHOLD ANALYSIS
############################################################

thresholds <- c(0.5, 0.25, 0.1, 0.05, 0.01)

for(t in thresholds){
  
  early <- sum(df$RNA.Early.spermatids >= t, na.rm = TRUE)
  late  <- sum(df$RNA.Late.spermatids >= t, na.rm = TRUE)
  
  cat("\nThreshold:", t, "\n")
  cat("Early:", early, "\n")
  cat("Late :", late, "\n")
}

############################################################
# 🧬 4. SAVE GENE LISTS
############################################################

genes_early <- df$Gene[df$RNA.Early.spermatids >= 0.5]
genes_late  <- df$Gene[df$RNA.Late.spermatids >= 0.5]

write.table(genes_early, "../results/genes_early.txt",
            row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(genes_late, "../results/genes_late.txt",
            row.names=FALSE, col.names=FALSE, quote=FALSE)

############################################################
# 📈 5. DISTRIBUTION PLOT
############################################################

df_long <- data.frame(
  Expression = c(df$RNA.Early.spermatids, df$RNA.Late.spermatids),
  Stage = rep(c("Early","Late"), each = nrow(df))
)

ggplot(df_long, aes(x = log10(Expression + 1e-4), fill = Stage)) +
  geom_histogram(bins = 80, alpha = 0.6) +
  theme_classic() +
  scale_fill_manual(values = c("#8B0000", "#1F4E79"))

############################################################
# 📉 6. DETECTION CURVE
############################################################

thresholds <- seq(0, 1, by = 0.01)

early_counts <- sapply(thresholds, function(t){
  sum(df$RNA.Early.spermatids >= t)
})

late_counts <- sapply(thresholds, function(t){
  sum(df$RNA.Late.spermatids >= t)
})

df_curve <- data.frame(
  threshold = rep(thresholds, 2),
  genes = c(early_counts, late_counts),
  stage = rep(c("Early","Late"), each = length(thresholds))
)

ggplot(df_curve, aes(threshold, genes, color = stage)) +
  geom_line(size = 1.2) +
  theme_classic()

############################################################
# 🔬 7. ENRICHMENT ANALYSIS
############################################################

other_cells <- df[, c(
  "RNA.Late.spermatocytes",
  "RNA.Early.spermatocytes",
  "RNA.Hub.cells",
  "RNA.Cyst",
  "RNA.Late.spermatogonia",
  "RNA.GSC..Early.spermatogonia",
  "RNA.Epithelial.cells"
)]

max_other <- apply(other_cells, 1, max)

classify_fc <- function(fc, threshold = 0.25){
  
  early <- df$RNA.Early.spermatids >= fc * max_other &
           df$RNA.Early.spermatids >= threshold
  
  late  <- df$RNA.Late.spermatids >= fc * max_other &
           df$RNA.Late.spermatids >= threshold
  
  early_genes <- df$Gene[early]
  late_genes  <- df$Gene[late]
  
  shared <- intersect(early_genes, late_genes)
  
  cat("\nFC:", fc, "\n")
  cat("Early:", length(early_genes), "\n")
  cat("Late :", length(late_genes), "\n")
  cat("Shared:", length(shared), "\n")
}

classify_fc(1.5)
classify_fc(2)
classify_fc(3)

############################################################
# 🧪 8. SCATTER PLOTS (PUBLICATION STYLE)
############################################################

t <- 0.25

df_plot <- df %>%
  mutate(
    FC = RNA.Late.spermatids / (max_other + 1e-6),
    category = case_when(
      RNA.Late.spermatids < t ~ "low",
      FC >= 3   ~ "3x",
      FC >= 2   ~ "2x",
      FC >= 1.5 ~ "1.5x",
      TRUE ~ "no"
    )
  )

cell_types <- colnames(other_cells)

plot_scatter <- function(cell){
  
  ggplot(df_plot, aes(x = .data[[cell]] + 1e-6,
                      y = RNA.Late.spermatids + 1e-6)) +
    
    geom_point(alpha = 0.3, color = "grey") +
    
    geom_point(data = df_plot %>% filter(category != "no"),
               aes(color = category)) +
    
    geom_abline(slope = 1, linetype = "dashed") +
    
    scale_x_log10() +
    scale_y_log10() +
    
    theme_classic() +
    labs(title = cell)
}

plots <- lapply(cell_types, plot_scatter)

final_plot <- wrap_plots(plots)

ggsave("../results/scatter_panel.png",
       final_plot, width = 12, height = 10, dpi = 300)

############################################################
# 🔥 9. HEATMAP (OPTIONAL)
############################################################

mat <- df[, grep("RNA", colnames(df))]
rownames(mat) <- df$Gene

pheatmap(log2(mat + 1))
