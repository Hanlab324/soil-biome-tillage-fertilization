# ============================================================
# Random Forest Analysis with rfPermute
# Purpose: Identify important predictors for NC_EMF and assess significance
# Data: emf15_ran.csv
# ============================================================

# Clear workspace
rm(list = ls())

# ============================================================
# 1. Load required packages
# ============================================================
library(randomForest)
library(rfPermute)
library(tidyverse)

# ============================================================
# 2. Data preparation
# ============================================================

# Read data
mydata <- read.csv("E:/emf_ran.csv", row.names = 1)

# Extract response and predictor variables
response <- mydata$nc_emf
predictors <- mydata[, !(names(mydata) %in% c("nc_emf", names(mydata)[1:5]))]

# Combine into RF-ready dataframe
rf_data <- data.frame(nc_emf = response, predictors)

# ============================================================
# 3. Manual parameter optimization (mtry selection)
# ============================================================

set.seed(123)

# Function to evaluate model error for a given mtry value
evaluate_mtry <- function(mtry_value, data, ntree = 1000, folds = 10, repeats = 3) {
  all_errors <- c()
  n <- nrow(data)
  
  for (r in 1:repeats) {
    fold_indices <- sample(rep(1:folds, length.out = n))
    for (f in 1:folds) {
      test_rows <- which(fold_indices == f)
      train_data <- data[-test_rows, ]
      test_data <- data[test_rows, ]
      
      rf_fold <- randomForest(nc_emf ~ ., data = train_data,
                              ntree = ntree, mtry = mtry_value)
      
      pred <- predict(rf_fold, newdata = test_data)
      mse <- mean((test_data$nc_emf - pred)^2)
      all_errors <- c(all_errors, mse)
    }
  }
  return(mean(all_errors))
}

# Define mtry range to test
mtry_range <- seq(1, ncol(predictors), by = 1)

# Calculate CV error for each mtry value
cv_errors <- sapply(mtry_range, function(m) {
  cat("Testing mtry =", m, "\n")
  evaluate_mtry(m, rf_data)
})

# Find optimal mtry
best_mtry <- mtry_range[which.min(cv_errors)]
cat("\nOptimal mtry =", best_mtry, "\n")

# Plot mtry vs error relationship
plot(mtry_range, cv_errors, type = "b", 
     xlab = "mtry", ylab = "CV Error (MSE)",
     main = "Cross-validation for mtry selection")
abline(v = best_mtry, col = "red", lty = 2)

# ============================================================
# 4. Fit final rfPermute model with optimal parameters
# ============================================================

rf_perm <- rfPermute(
  nc_emf ~ .,
  data = rf_data,
  ntree = 1000,
  mtry = best_mtry,
  importance = TRUE,
  nrep = 500,
  num.cores = 4
)

# ============================================================
# 5. Extract importance and p-values
# ============================================================

# Extract importance values
imp <- importance(rf_perm, scale = TRUE)
p_values <- rf_perm$pval

# Extract scaled %IncMSE p-values and reorder
p_values_scaled <- p_values[, "%IncMSE", "scaled"]
p_values_scaled <- p_values_scaled[rownames(imp)]

# Create importance dataframe
imp_df <- data.frame(
  Variable = rownames(imp),
  IncMSE = imp[, "%IncMSE"],
  p_value = p_values_scaled
)

# Order variables as in original predictors
var_order <- colnames(predictors)
imp_df$Variable <- factor(imp_df$Variable, levels = var_order)

# Add significance markers
imp_df <- imp_df %>%
  mutate(
    sig_level = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    color = ifelse(sig_level == "ns", "#2687bc", "#A2142F")
  )

# Check significance distribution
cat("\nSignificance distribution:\n")
print(table(imp_df$sig_level))

# ============================================================
# 6. Calculate cross-validated R2
# ============================================================

y_var <- var(rf_data$nc_emf)
best_mse <- min(cv_errors)
cv_r2 <- 1 - (best_mse / y_var)
cat("\nCross-validated R2 =", round(cv_r2, 3), "\n")

# ============================================================
# 7. Visualization
# ============================================================

ggplot(imp_df, aes(x = Variable, y = IncMSE, fill = color)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_identity() +
  geom_text(
    aes(label = ifelse(sig_level != "ns", sig_level, "")),
    vjust = -0.5,
    size = 5
  ) +
  labs(
    x = "",
    y = "Increase in MSE (%) - NC_EMF",
    title = "Variable Importance for NC_EMF"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_blank(),
    axis.ticks.x = element_line(color = "black"),
    axis.line.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5)
  ) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = paste0("R2 = ", round(cv_r2, 2)),
    hjust = -0.1,
    vjust = 1.5,
    size = 5
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

