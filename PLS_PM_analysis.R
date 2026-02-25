############################################
# PLS-PM Analysis of EMF
# Description: Complete PLS-PM workflow including:
# - Data import
# - SEM measurement and structural model
# - Path coefficients, effects, and GOF
# - Visualization
# - Network table
############################################

# Load required packages
library(plspm)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(stringr)

# ======================================================
# 1. Data Preparation
# ======================================================

# Read data
data <- read.csv("E:/aaa2019zydna/plspm2/all_emf8.csv", 
                 header = TRUE, check.names = FALSE, row.names = 1)
data$sample <- rownames(data)

# Read structural matrix (path relationships)
struct <- read.csv("E:/aaa2019zydna/plspm2/structure_crop.csv", 
                   header = TRUE, check.names = FALSE)

# Read measurement model (factor definitions)
measure <- read.csv("E:/aaa2019zydna/plspm2/factor_crop.csv", 
                    header = TRUE, check.names = FALSE)

# ======================================================
# 2. Build PLS-PM Model
# ======================================================

# Prepare model parameters
group <- data.frame(factors = measure$target, group = measure$source)
sm <- as.matrix(struct)
mm <- as.matrix(measure)

# Check if latent variables exist in data
latent <- unique(as.vector(sm))
latent[latent %in% colnames(data)]

# Build model
model <- plsm(data, sm, mm)
path <- t(model$D)

# Save path matrix
write.table(path, "E:/aaa2019zydna/plspm2/pls_sj/step57.SEM.PLSPM_PNA.structural_matrix.xls", 
            sep = "\t", col.names = NA, quote = FALSE)

# ======================================================
# 3. Define Blocks and Modes
# ======================================================

# Order by group
group$group <- factor(group$group, levels = rownames(path))
group <- group[order(group$group), ]

# Generate block list
temp <- lapply(1:length(unique(group[, 2])), function(x) {
  group[which(group[, 2] %in% unique(group[, 2])[x]), 1]
})
block <- lapply(1:length(temp), function(x) which(colnames(data) %in% temp[[x]]))
mod <- rep("A", length(block))

# ======================================================
# 4. Run PLS-PM (with bootstrap validation)
# ======================================================

set.seed(123)  # Set random seed for reproducibility
pls <- plspm(data, path, block, modes = mod, boot.val = TRUE, br = 1000)

# Print model fit
print(pls$gof)

# ======================================================
# 5. Save Model Results
# ======================================================

# 5.1 Latent variable scores
latent_scores <- pls$scores
write.csv(latent_scores, "E:/aaa2019zydna/plspm2/pls_sj/latent_scores.csv", quote = FALSE)

# 5.2 Outer model (loadings)
write.table(pls$outer_model, "E:/aaa2019zydna/plspm2/pls_sj/step57.SEM.outer_model.txt", 
            row.names = FALSE, sep = '\t', quote = FALSE)

# 5.3 Effects decomposition
write.table(pls$effects, "E:/aaa2019zydna/plspm2/pls_sj/step57.SEM.effects.txt", 
            row.names = FALSE, sep = '\t', quote = FALSE)

# 5.4 Inner model summary
write.table(pls$inner_summary, "E:/aaa2019zydna/plspm2/pls_sj/step57.SEM.PLSPM_PNA.inner_summary.xls", 
            sep = "\t", col.names = NA, quote = FALSE)

# 5.5 Path coefficients matrix
write.table(pls$path_coefs, "E:/aaa2019zydna/plspm2/pls_sj/step57.SEM.PLSPM_PNA.path_coefs.xls", 
            sep = "\t", col.names = NA, quote = FALSE)

# ======================================================
# 6. Process Inner Model Results (for network visualization)
# ======================================================

# Extract inner model
inner_model_trans <- do.call(rbind, lapply(pls$inner_model, data.frame))

# Remove intercept terms
intercept <- inner_model_trans[grep(pattern = "Intercept", rownames(inner_model_trans)), ]
inner_model_trans <- inner_model_trans[-which(rownames(inner_model_trans) %in% rownames(intercept)), ]

# Rename columns
colnames(inner_model_trans) <- c("Estimate", "Std.Error", "t value", "P")

# Extract path names
p_rowname <- strsplit(rownames(inner_model_trans), "[.]")
p_rowname <- data.frame(matrix(unlist(p_rowname), nrow = length(p_rowname), byrow = TRUE))
inner_model_trans <- cbind(p_rowname, inner_model_trans)
colnames(inner_model_trans) <- c("target", "source", "coef", "std.error", "t value", "P")

# Save inner model
write.table(inner_model_trans, "E:/aaa2019zydna/plspm2/pls_sj/step57.SEM.PLSPM_PNA.inner_model.xls", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# ======================================================
# 7. Prepare Network Data (Cytoscape format)
# ======================================================

inner_model_trans$coef <- round(inner_model_trans$coef, digits = 4)
inner_model_trans$linewidth <- abs(inner_model_trans$coef)
inner_model_trans$po_ne <- ifelse(inner_model_trans$coef > 0, 'positive', 'negative')

# Add significance markers
sig <- rep("", nrow(inner_model_trans))
sig[which(inner_model_trans$P < 0.05)] <- "*"
sig[which(inner_model_trans$P < 0.01)] <- "**"
sig[which(inner_model_trans$P < 0.001)] <- "***"
inner_model_trans$sig <- sig

# Save complete network table
write.table(inner_model_trans, "E:/aaa2019zydna/plspm2/pls_sj/step57.SEM.PLSPM_PNA.inner_model_for_cytoscape.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Filter and prepare network data (keep positive coefficients)
inner_filter <- inner_model_trans[inner_model_trans$coef > 0.0, ]
network <- inner_filter
network$edge_label <- paste0(round(network$coef, 2), network$sig)

# Set edge colors and widths
edge_color <- rep("black", length(rownames(network)))
edge_color[which(network$coef > 0)] <- "red"
edge_color[which(network$coef < 0)] <- "blue"
network$edge_color <- edge_color
network$edge_width <- abs(network$coef)

write.table(network, "E:/aaa2019zydna/plspm2/pls_sj/step57.SEM.PLSPM_PNA.network_table.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# ======================================================
# 8. Effects Decomposition Visualization (Effects on EMFunctions)
# ======================================================

# Read effects file
df <- read.delim('E:/aaa2019zydna/plspm2/pls_sj/step57.SEM.effects.txt', header = TRUE)

# Split paths
temp <- str_split_fixed(df$relationships, ' -> ', n = 2)
colnames(temp) <- c('Source', 'Target')
df2 <- as.data.frame(cbind(df, temp))

# Filter for target variable EMFunctions
df3 <- df2 %>% dplyr::filter(Target == 'EMFunctions') %>% 
  dplyr::select(Source, direct, indirect, total)

# Reshape to long format
df4 <- melt(df3, id.vars = 'Source')

# Plot effects decomposition
ggplot(data = df4, aes(x = Source, y = value, fill = variable)) +
  geom_bar(stat = 'identity', position = position_dodge(), color = 'black') +
  scale_fill_manual(values = c('#f09461', '#fff1c3', '#90c18e')) +
  geom_hline(yintercept = 0) +
  labs(x = NULL, y = 'Standardized effects on EMFunctions (PLS-PM)') +
  theme_classic() +
  theme(axis.line.x = element_blank(),
        legend.position = 'top',
        legend.title = element_blank(),
        axis.ticks.x = element_blank())

# ggsave('step57.SEM.direct_indirect_total.pdf', width = 6, height = 4)

# ======================================================
# 9. Plot PLS-PM Path Diagram
# ======================================================

# pdf("step56.SEM.plot_SEM_PLSPM_PNA.pdf", width = 20, height = 20)
innerplot(pls, 
          colpos = "red", colneg = "blue", lcol = "black", 
          box.col = "white", box.size = 0.06, box.cex = 2, 
          arr.pos = 0.2, arr.width = 0.1, arr.lwd = 1,
          txt.col = "black", cex.txt = 1.5, dtext = -1, show.values = TRUE)
# dev.off()

# ======================================================
# 10. Recalculate P-values Based on Bootstrap (More Accurate)
# ======================================================

# Extract bootstrap results from pls object
paths_boot <- pls$boot$paths

# Calculate p-values based on normal approximation
paths_df <- data.frame(
  Path = rownames(paths_boot),
  Coefficient = paths_boot[, "Original"],
  Boot_SE = paths_boot[, "Std.Error"],
  t_value = paths_boot[, "Original"] / paths_boot[, "Std.Error"],
  P_value = 2 * (1 - pnorm(abs(paths_boot[, "Original"] / paths_boot[, "Std.Error"])))
)

# Add significance markers
paths_df$Significance <- ifelse(
  paths_df$P_value < 0.001, "***",
  ifelse(paths_df$P_value < 0.01, "**",
         ifelse(paths_df$P_value < 0.05, "*", "ns"))
)

# View results
print(paths_df)

# ======================================================
# 11. Save Bootstrap Validation Results
# ======================================================

write.csv(paths_df, "E:/aaa2019zydna/plspm2/pls_sj/PLS_path_coefficients_bootstrap.csv", 
          row.names = FALSE)

