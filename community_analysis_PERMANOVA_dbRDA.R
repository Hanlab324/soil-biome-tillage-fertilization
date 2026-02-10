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
# Transform ASV table
# ------------------------------
asv_hel <- decostand(asv, method = 'hellinger')

# ------------------------------
# PERMANOVA
# ------------------------------
adonis2(vegdist(sqrt(asv_hel), method = "bray") ~ 
          env_data1$Fertilization * env_data1$Tillage, 
        strata = env_data1$block)

# ------------------------------
# dbRDA
# ------------------------------
rda_db <- dbrda(asv_hel ~ ., env, scale = TRUE)
rda_db.scaling1 <- summary(rda_db, scaling = 1)

rda_site.scaling1 <- rda_db.scaling1$site[, 1:4]
rda_sp.scaling1 <- rda_db.scaling1$species[, 1:4]
rda_env.scaling1 <- rda_db.scaling1$biplot[, 1:4]

# Adjusted R2
r2 <- RsquareAdj(rda_db)
rda_noadj <- r2$r.squared
rda_adj <- r2$adj.r.squared

# Significance test
rda_db_test <- anova(rda_db, permutations = 999)
rda_db_test_axis <- anova(rda_db, by = "axis", permutations = 999)
rda_db_test_axis$`Pr(>F)` <- p.adjust(rda_db_test_axis$`Pr(>F)`, method = "bonferroni")
vif.cca(rda_db)

# ------------------------------
# Forward selection
# ------------------------------
rda_db_forward_p <- ordistep(rda(asv_hel ~ 1, env, scale = FALSE), 
                             scope = formula(rda_db), 
                             direction = 'forward', 
                             permutations = 999)

rda_db_forward_r <- ordiR2step(rda(asv_hel ~ 1, env, scale = FALSE), 
                               scope = formula(rda_db), 
                               R2scope = rda_adj, 
                               direction = 'forward', 
                               permutations = 999)

# Extract forward selection results
anova_forward <- rda_db_forward_r$anova
summary(rda_db_forward_r)
summary(rda_db_forward_p)

# ------------------------------
# Plot dbRDA
# ------------------------------
rda_db_forward_r.scaling1 <- summary(rda_db_forward_r, scaling = 1)
rda_db_forward_r.site <- data.frame(rda_db_forward_r.scaling1$sites)[1:2]
rda_db_forward_r.env <- data.frame(rda_db_forward_r.scaling1$biplot)[1:2]

rda_db_forward_r.site$sample <- rownames(rda_db_forward_r.site)
a_new2 <- cbind(rda_db_forward_r.site, asv1[, 2:3])
rda_db_forward_r.env$sample <- rownames(rda_db_forward_r.env)

a_new2$Fertilization <- factor(env_data1$Fertilization, levels = c("Ctrl", "N", "P", "NP", "M", "MNP"))
a_new2$Tillage <- factor(env_data1$Tillage, levels = c("CT", "NT"))

rda_plot <- ggplot(a_new2, aes(RDA1, RDA2)) +
  geom_point(aes(color = Tillage, shape = Fertilization), size = 3, stroke = 1.2) +
  scale_color_manual(values = c("86", "19")) +
  scale_fill_manual(values = c("86", "19")) +
  scale_shape_manual(values = c(1,0,2,7,8,12)) +     
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5, size = 18), 
        legend.key = element_rect(fill = 'transparent'),
        axis.text.x = element_text(size = rel(1.5), colour = "black"),
        axis.text.y = element_text(size = rel(1.5), colour = "black"),
        axis.title.x = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(1.5)),
        legend.position = "none",
        text = element_text(face = "bold")) +
  labs(title = "Bacteria",
       x = paste("dbRDA1", round(rda_db_forward_r$CCA$eig[1]/sum(rda_db_forward_r$CCA$eig)*100*RsquareAdj(rda_db_forward_p)$adj.r.squared, 2), "%"),
       y = paste("dbRDA2", round(rda_db_forward_r$CCA$eig[2]/sum(rda_db_forward_r$CCA$eig)*100*RsquareAdj(rda_db_forward_p)$adj.r.squared, 2), "%")) +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  geom_segment(data = rda_db_forward_r.env, aes(x = 0, y = 0, xend = RDA1*0.7, yend = RDA2*0.7), 
               arrow = arrow(length = unit(0.2, 'cm')), lwd = 0.8, size = 0.4, color = 'black') +
  geom_text(data = rda_db_forward_r.env, aes(RDA1*0.79, RDA2*0.79, label = sample), color = 'black', size = 6)

rda_plot

