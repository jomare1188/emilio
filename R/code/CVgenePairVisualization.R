library(ggplot2)

# *** Reset R variables ***
rm(list = ls())

DIR="/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/correr/AnalysisCV/top20CV"

setwd(DIR)

# Read coefficient of variation file
genes <- read.table(file.path(DIR,"Correr2020_top20CV.txt"), header = TRUE)

# Read quantification matrix
quantification_matrix <- read.table(file.path(DIR,'Correr2020_counts_filters_VST_top20CV.txt'), header = TRUE, row.names = 1, check.names = FALSE)

# Split genes in quartiles
genes$Quartil <- cut(genes$CV, breaks = quantile(genes$CV, probs = seq(0, 1, 0.25), na.rm = TRUE), include.lowest = TRUE, labels = FALSE)

# Quantile analysis
quantile(genes$CV, probs = seq(0, 1, 0.25), na.rm = TRUE)

quartis <- quantile(genes$CV, probs = seq(0, 1, 0.25), na.rm = TRUE)
percent_quartis <- format(quartis, digits = 2) 
names(percent_quartis) <- c("Min CV", "25%", "50%", "75%", "Max CV")
print(percent_quartis)

# Def gene pairs per plot
n_pairs_per_quartil <- 6

# Get random pairs for each quantile
random_gene_pairs <- lapply(split(genes, genes$Quartil), function(subset) {
  sample(subset$Gene, size = 2 * n_pairs_per_quartil, replace = TRUE)
})

percent_range <- c("0-25%", "25-50%", "50-75%", "75-100%")

for (i in seq_along(random_gene_pairs)) {
  matrix_gene_pairs <- matrix(random_gene_pairs[[i]], ncol = 2, byrow = TRUE)
  
  plots <- list()
  
  for (j in seq_len(n_pairs_per_quartil)) {
    gene_pair <- matrix_gene_pairs[j, ]
    gene_pair_quantification <- quantification_matrix[gene_pair, ]
    
    gene_labels <- rownames(gene_pair_quantification)
   
    plot_data <- data.frame(
      Condition = rep(colnames(gene_pair_quantification), each = 1),
      Counts_Gene1 = as.vector(t(gene_pair_quantification[1, ])),
      Counts_Gene2 = as.vector(t(gene_pair_quantification[2, ])),
      Gene_Color1 = gene_labels[1],
      Gene_Color2 = gene_labels[2]
    )
    
    plot <- ggplot(plot_data, aes(x = Condition)) +
      geom_point(aes(y = Counts_Gene1, color = Gene_Color1), shape = 1) +
      geom_point(aes(y = Counts_Gene2, color = Gene_Color2), shape = 2) +
      labs(title = paste("Q", i, "(", percent_range[i], ")", '- Gene Pair', j, "(top 50% CV)"),
           x = 'Samples',
           y = 'Counts (VST)',
           color = "Genes") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    plots[[length(plots) + 1]] <- plot
  }
  
  # Save one graph per quantile (one quantile and 6 pairs of genes)
  ggsave(paste0("gene_pairs_quartil_", i, ".png"), do.call(gridExtra::arrangeGrob, plots), width = 15, height = 10)
}
