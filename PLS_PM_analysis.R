PLS_PM_analysis
############################################
# PLS-PM Analysis of EMF
# Description: Complete PLS-PM workflow including:
# - Data import
# - SEM measurement and structural model
# - Path coefficients, effects, and GOF
# - Visualization
# - Network table
############################################

# ----------------------------
# 1. Load required packages
# ----------------------------
library(devtools)
library(lattice)
library(semPLS)
library(plspm)
library(reshape2)
library(dplyr)
library(ggplot2)
library(stringr)

# Clear environment
rm(list = ls())

# ----------------------------
# 2. Data Import
# ----------------------------
# Main data (EMFunctions)
data1 <- read.csv("E:/aaa2019zydna/plspm2/all_emf8.csv", 
                  header = TRUE, check.names = FALSE, row.names = 1)
data1$sample <- rownames(data1)
data <- data1

# Structural model
struct <- read.csv("E:/aaa2019zydna/plspm2/structure_nsc1.csv", 
                   header = TRUE, check.names = FALSE)

# Measurement model
measure <- read.csv("E:/aaa2019zydna/plspm2/factor_nsc1.csv", 
                    header = TRUE, check.names = FALSE)

# ----------------------------
# 3. Prepare PLS-PM input
# ----------------------------
group <- data.frame(factors = measure$target, group = measure$source)
sm <- as.matrix(struct)
mm <- as.matrix(measure)

# Identify latent variables in data
latent <- unique(as.vector(sm))
latent[latent %in% colnames(data)]

# ----------------------------
# 4. SEM Measurement Model
# ----------------------------
model <- plsm(data, sm, mm)
path <- t(model$D)

# Export structural matrix
write.table(path, 
            "E:/aaa2019zydna/plspm2/pls_sj/step57.SEM.PLSPM_PNA.structural_matrix.xls", 
            sep = "\t", col.names = NA, quote = FALSE)

# ----------------------------
# 5. Prepare blocks for plspm()
# ----------------------------
group$group <- factor(group$group, levels = rownames(path))
group <- group[order(group$group), ]

temp <- lapply(1:length(unique(group[, 2])), function(x) {
  group[which(group[, 2] %in% unique(group[, 2])[x]), 1]
})

block <- lapply(1:length(temp), function(x) {
  which(colnames(data) %in% temp[[x]])
})

# All latent variables are reflective
mod <- rep("A", length(block))

# ----------------------------
# 6. Run PLS-PM
# ----------------------------
pls <- plspm(data, path, block, modes = mod, boot.val = FALSE)

# Goodness-of-fit
pls$gof

# Export outer model
write.table(pls$outer_model,
            "E:/aaa2019zydna/plspm2/pls_sj/step57.SEM.outer_model.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)

# ----------------------------
# 7. Effects
# ----------------------------
write.table(pls$effects,
            "E:/aaa2019zydna/plspm2/pls_sj/step57.SEM.effects.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)

df <- read.delim("E:/aaa2019zydna/plspm2/pls_sj/step57.SEM.effects.txt", header = TRUE)
temp <- str_split_fixed(df$relationships, " -> ", n = 2)
colnames(temp) <- c("Source", "Target")
df2 <- cbind(df, temp)

# Select EMFunctions target
df3 <- df2 %>% filter(Target == "EMFunctions") %>% select(Source, direct, indirect, total)
df4 <- melt(df3, id.vars = "Source")

# Bar plot: direct, indirect, total effects
ggplot(df4, aes(x = Source, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  scale_fill_manual(values = c('#f09461', '#fff1c3', '#90c18e')) +
  geom_hline(yintercept = 0) +
  labs(x = NULL, y = "Standardized on EMFunctions from PLS_PM") +
  theme_classic() +
  theme(axis.line.x = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        axis.ticks.x = element_blank())

# ----------------------------
# 8. Plot internal model
# ----------------------------
innerplot(pls,
          colpos = "red", colneg = "blue", lcol = "black", 
          box.col = "white", box.size = 0.06, box.cex = 2, 
          arr.pos = 0.2, arr.width = 0.1, arr.lwd = 1,
          txt.col = "black", cex.txt = 1.5, dtext = -1,
          show.values = TRUE)

# Export summary and path coefficients
write.table(pls$inner_summary,
            "E:/aaa2019zydna/plspm2/pls_sj/step57.SEM.PLSPM_PNA.inner_summary.xls",
            sep = "\t", col.names = NA, quote = FALSE)

write.table(pls$path_coefs,
            "E:/aaa2019zydna/plspm2/pls_sj/step57.SEM.PLSPM_PNA.path_coefs.xls",
            sep = "\t", col.names = NA, quote = FALSE)

# ----------------------------
# 9. Transform inner model for Cytoscape/network
# ----------------------------
inner_model_trans <- do.call(rbind, lapply(pls$inner_model, data.frame))
intercept <- inner_model_trans[grep("Intercept", rownames(inner_model_trans)), ]
inner_model_trans <- inner_model_trans[!rownames(inner_model_trans) %in% rownames(intercept), ]

colnames(inner_model_trans) <- c("Estimate", "Std.Error", "t.value", "P")
p_rowname <- strsplit(rownames(inner_model_trans), "[.]")
p_rowname <- data.frame(matrix(unlist(p_rowname), nrow = nrow(inner_model_trans), byrow = TRUE))
inner_model_trans <- cbind(p_rowname, inner_model_trans)
colnames(inner_model_trans) <- c("target", "source", "coef", "std.error", "t.value", "P")

# Round and add positive/negative labels
inner_model_trans$coef <- round(inner_model_trans$coef, 4)
inner_model_trans$po_ne <- ifelse(inner_model_trans$coef > 0, "positive", "negative")

# Significance
sig <- rep("", nrow(inner_model_trans))
sig[inner_model_trans$P < 0.05] <- "*"
sig[inner_model_trans$P < 0.01] <- "**"
sig[inner_model_trans$P < 0.001] <- "***"
inner_model_trans$sig <- sig

# Filter for positive edges
inner_filter <- inner_model_trans[inner_model_trans$coef > 0, ]

# Network table
network <- inner_filter
network$edge_label <- paste0(round(network$coef, 2), network$sig)
network$edge_color <- ifelse(network$coef > 0, "red", "blue")
network$edge_width <- abs(network$coef)

write.table(network,
            "E:/aaa2019zydna/plspm2/pls_sj/step57.SEM.PLSPM_PNA.network_table.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

