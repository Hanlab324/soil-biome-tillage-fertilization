community_analysis_PERMANOVA_dbRDA
# ------------------------------
# Load ASV data
# ------------------------------
asv1 <- read.csv("bacteria.csv", row.names = 1, header = TRUE)
asv <- asv1

# ------------------------------
# Load environmental data
# ------------------------------
env_data1 <- read.csv("envirpcoa.csv", row.names = 1)
env <- env_data1[, -c(1:3)] # remove unnecessary columns

# ------------------------------
# PERMANOVA
# ------------------------------
bray_dist <- vegdist(otu, method = "bray")


permanova <- adonis2(bray_dist ~ Tillage * Fertilization, 
                     data = env_data1, blocks = env_data1$block)
# PERMDISP
# Fertilization
disp_fert <- betadisper(bray_dist, env_data1$Fertilization)
permdisp_fert <- permutest(disp_fert, permutations = 999)

# Tillage
disp_till <- betadisper(bray_dist, env_data1$Tillage)
permdisp_till <- permutest(disp_till, permutations = 999)

# Interaction
env_data1$int <- interaction(env_data1$Tillage, env_data1$Fertilization)
disp_int <- betadisper(bray_dist, env_data1$int)
permdisp_int <- permutest(disp_int, permutations = 999)


# ------------------------------
# dbRDA

############################################################
library(vegan)
library(ggplot2)

# ----------------------------------------------------------
# 1. Import data
# ----------------------------------------------------------

asv <- read.csv("bacteria.csv", row.names = 1, header = TRUE)

env_data <- read.csv("envirpcoa.csv", row.names = 1)

# Remove categorical columns (first 3 columns assumed metadata)
env <- env_data[, -c(1:3)]

# ----------------------------------------------------------
# 2. dbRDA (Bray–Curtis distance)
# ----------------------------------------------------------

rda_full <- dbrda(asv ~ ., data = env, distance = "bray", scale = TRUE)

# Model significance
anova_full <- anova(rda_full, permutations = 999)
anova_axis <- anova(rda_full, by = "axis", permutations = 999)

# Adjust P values for axes
anova_axis$`Pr(>F)` <- p.adjust(anova_axis$`Pr(>F)`, method = "bonferroni")

# Variance inflation factors
vif_values <- vif.cca(rda_full)

# Adjusted R2
r2 <- RsquareAdj(rda_full)
R2_adj <- r2$adj.r.squared
R2_raw <- r2$r.squared

# ----------------------------------------------------------
# 3. Forward selection (based on adjusted R2)
# ----------------------------------------------------------

rda_null <- dbrda(asv ~ 1, data = env, distance = "bray")

rda_forward <- ordiR2step(
  rda_null,
  scope = formula(rda_full),
  R2scope = R2_adj,
  direction = "forward",
  permutations = 999
)

summary(rda_forward)

# ----------------------------------------------------------
# 4. Extract site and environmental scores
# ----------------------------------------------------------

scores_forward <- summary(rda_forward, scaling = 1)

site_scores <- as.data.frame(scores_forward$sites[, 1:2])
env_scores  <- as.data.frame(scores_forward$biplot[, 1:2])

site_scores$Sample <- rownames(site_scores)

# Add treatment information
site_scores$Fertilization <- factor(env_data$Fertilization)
site_scores$Tillage <- factor(env_data$Tillage)

# ----------------------------------------------------------
# 5. Calculate explained variance for axes
# ----------------------------------------------------------

constrained_eig <- eigenvals(rda_forward, model = "constrained")

axis1_prop <- constrained_eig[1] / sum(constrained_eig)
axis2_prop <- constrained_eig[2] / sum(constrained_eig)

axis1_total <- axis1_prop * RsquareAdj(rda_forward)$r.squared * 100
axis2_total <- axis2_prop * RsquareAdj(rda_forward)$r.squared * 100

# ----------------------------------------------------------
# 6. Plot dbRDA
# ----------------------------------------------------------
# 给 site_scores 添加列名
colnames(site_scores)[1:2] <- c("RDA1", "RDA2")
colnames(env_scores)[1:2] <- c("RDA1", "RDA2")

p <- ggplot(site_scores, aes(RDA1, RDA2)) +
  geom_point(aes(color = Tillage, shape = Fertilization), size = 3) +
  scale_color_manual(values = c("86", "19")) +
  scale_fill_manual(values = c("86", "19")) +
  scale_shape_manual(values = c(1, 0, 2, 7, 8, 12)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_segment(
    data = env_scores,
    aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
    arrow = arrow(length = unit(0.2, "cm")),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = env_scores,
    aes(RDA1, RDA2, label = rownames(env_scores)),
    inherit.aes = FALSE,
    vjust = -0.5
  ) +
  labs(
    x = paste0("dbRDA1 (", round(axis1_total, 2), "%)"),
    y = paste0("dbRDA2 (", round(axis2_total, 2), "%)"),
    title = "Bacterial community (dbRDA)"
  ) +
  theme_bw()

p
