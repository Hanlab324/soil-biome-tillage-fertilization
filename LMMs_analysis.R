LMM_analysis 
#---------------------------------------
# Soil multifunctionality LMMs analysis
#---------------------------------------

# 1. Load libraries
library(lme4)      # For linear mixed models
library(lmerTest)  # To get p-values in lmer
library(tidyverse) # Optional, for data handling

# 2. Load dataset
# Replace the path with your own CSV file location
a <- read.csv("E:/aaa2019zydna/newsummary.csv", header = TRUE, row.names = 1)

# 3. Set factor levels
# Fertilization: Control, N, P, NP, M, MNP
a$Fertilization <- factor(a$Fertilization, 
                          levels = c("Ctrl","N","P","NP","M","MNP"))

# Tillage: Conventional (CT), No-tillage (NT)
a$Tillage <- factor(a$Tillage, levels = c("CT","NT"))

# Block is assumed to be a column in your data
a$block <- factor(a$block)

# 4. Fit linear mixed model
# Response: nccrop (nutrient-cycling/crop productivity index, adjust as needed)
# Fixed effects: Tillage, Fertilization, and their interaction
# Random effects: block nested within Tillage
model <- lmer(nccrop ~ Tillage * Fertilization + (1 | block/Tillage), data = a)

# 5. ANOVA for fixed effects
anova_results <- anova(model)
print(anova_results)
