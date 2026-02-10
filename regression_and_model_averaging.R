regression_and_model_averaging 
# ------------------------------
# Load libraries and data
# ------------------------------
library(MuMIn)
library(ggplot2)
library(dplyr)
library(cowplot)

alldata <- read.csv("E:/aaa2019zydna/functions_heat.csv", header = TRUE, row.names = 1)
soil <- read.csv("E:/aaa2019zydna/newsummary.csv", header = TRUE, row.names = 1)

# ------------------------------
# Prepare model dataset
# ------------------------------
temp <- data.frame(var = alldata$nc_emf, alldata[, 38:46], soil[, 11:15])
model_data <- data.frame(scale(temp))  # standardize
class(model_data)

# ------------------------------
# Fit full linear model
# ------------------------------
fit <- lm(var ~ ., data = model_data)
summary(fit)

options(na.action = 'na.fail')  # required for dredge

# Model selection via dredge
fit_select <- dredge(fit)
head(fit_select)

# Select models with delta AIC < 2
subset(fit_select, delta < 2)
de12 <- model.avg(fit_select, subset = delta < 2)
fit_avg_model <- summary(de12)

# ------------------------------
# Extract significant variables and refit
# ------------------------------
significant_vars <- names(which(!is.na(coef(de12)[-1])))  # exclude intercept
formula_avg <- as.formula(paste("var ~", paste(significant_vars, collapse = " + ")))
fit <- lm(formula_avg, data = model_data)
summary(fit)

# ------------------------------
# Prepare stats for plotting
# ------------------------------
stat.lm <- fit_avg_model$coefmat.subset[-1, ]  # remove intercept
stat.lm <- data.frame(stat.lm, check.names = FALSE)
stat.lm$env <- rownames(stat.lm)

# Assign significance stars
stat.lm$sig <- ifelse(stat.lm$'Pr(>|z|)' > 0.1, '',
                      ifelse(stat.lm$'Pr(>|z|)' > 0.05, '',
                             ifelse(stat.lm$'Pr(>|z|)' > 0.01, '*',
                                    ifelse(stat.lm$'Pr(>|z|)' > 0.001, '**', '***'))))
stat.lm$label <- paste(stat.lm$env, stat.lm$sig)

# Categorize variables
biotic_vars <- colnames(alldata[, 38:46])
abiotic_vars <- colnames(soil[, 11:15])

stat.lm$type <- NA
stat.lm[stat.lm$env %in% biotic_vars, 'type'] <- 'Biotic variables'
stat.lm[stat.lm$env %in% abiotic_vars, 'type'] <- 'Abiotic variables'

# Order and factor
stat.lm <- stat.lm[order(stat.lm$type, stat.lm$Estimate), ]
stat.lm$env <- factor(stat.lm$env, levels = stat.lm$env)
stat.lm$type <- factor(stat.lm$type, levels = c("Biotic variables", "Abiotic variables"))

# Compute relative effect percentages
stat.lm$env.ra <- abs(stat.lm$Estimate) / sum(abs(stat.lm$Estimate)) * 100
type_sum <- aggregate(env.ra ~ type, data = stat.lm, sum)

# ------------------------------
# Plot: parameter estimates
# ------------------------------
adjr2_val <- round(summary(fit)$adj.r.squared, 2)
p_model_estimate <- ggplot(stat.lm, aes(x = env, y = Estimate, color = type)) +
  geom_point(size = 6, shape = 15) +
  geom_errorbar(aes(ymin = Estimate - `Std. Error`, ymax = Estimate + `Std. Error`),
                width = 0.1, linewidth = 1, show.legend = FALSE) +
  scale_color_manual(values = c('#BD3786FF','#386CB0')) +
  scale_x_discrete(breaks = stat.lm$env, labels = stat.lm$label, position = 'top') +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = '', y = 'Parameter estimate', color = '') +
  theme_bw(base_size = 15) +
  coord_flip() +
  ggtitle(bquote("EMF (Crop) Adj.R"^2*" = "*.(adjr2_val))) +
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.position = "none")

# ------------------------------
# Plot: relative contributions by type
# ------------------------------
p_model_env <- ggplot(type_sum, aes(x = '', y = env.ra, fill = type)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(env.ra, 1), "%")),
            position = position_stack(vjust = 0.5),
            color = "white", size = 6, fontface = "bold") +
  labs(x = "", y = "Relative effects (%)", fill = "Subcategory") +
  scale_fill_manual(values = c('#BD3786FF','#386CB0')) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

# ------------------------------
# Combine plots
# ------------------------------
combined_plot <- plot_grid(
  p_model_env,
  p_model_estimate,
  align = "h",
  axis = "tb",
  rel_widths = c(2, 6)
)

combined_plot
