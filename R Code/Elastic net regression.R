
# July 2026

# Sarah Cordell et al.


### Getting setup -----------------------------------------------------------

library(tidyverse)
library(glmnet)
library(coefplot)
library(caret)
library(lmerTest)
library(ggtext)
library(patchwork)
library(grid)
library(readxl)

set.seed(1234567)


### Loading data ------------------------------------------------------------

# Load beta values 
beta <- read.csv('Data/data_sesame_40k.csv')
rownames(beta) <- beta$CGid
beta$CGid <- NULL

# Load metadata (sample information)
meta <- read.csv('Data/ET0835.DayTarongaPlatypus.m40k_metadata.csv')

# Load Universal Clock Coefficients
uniclk <- read.csv("Data/Universal_Coef.CombAnn.csv", stringsAsFactors = FALSE)


### Add useful columns -----------------------------------------------------

# Create AnimalID variable

# To be used for LOIOCV, to account for non-independent samples
meta <- meta %>%
  dplyr::mutate(
    AnimalID = stringr::str_extract(SampleID, "^(Plat|SBE)_[^_]+")
  )

# Sanity check, should = 0
sum(is.na(meta$AnimalID))

# Add common name to metadata
meta <- meta %>%
  mutate(
    Species = case_when(
      SpeciesLatinName == "Ornithorhynchus anatinus" ~ "Platypus",
      SpeciesLatinName == "Tachyglossus aculeatus" ~ "Short-beaked echidna",
      TRUE ~ NA_character_
    )
  )

# Add species origin to metadata

meta <- meta %>% 
  left_join(master %>% select(SampleID, individual.origin), by = "SampleID")


### Add life history information -------------------------------------------

# Information required by different clocks
species_params <- data.frame(
  SpeciesLatinName = c("Ornithorhynchus anatinus","Tachyglossus aculeatus"),
  AgeASM = c(3.5, 5),
  Gestation = c(0.044, 0.052),
  MaxAge = c(30, 50)
)

meta <- meta %>%
  left_join(species_params, by = "SpeciesLatinName")


### Quality control, remove outliers and other -----------------------------

# Remove outliers; hierarchical and probe fails
# Samples flagged as outliers using PCA only were not excluded due to aggressive nature of this analysis 
outliers <- c("ET083571", "ET083573", "ET0835118")

# Removing other questionable samples
questionable <- c("ET0835120", # cell culture sample
                  "ET08353", "ET083568", "ET083577") # age estimate high

sampletidy <- c(outliers, questionable)

beta <- beta[, !(colnames(beta) %in% sampletidy)]
meta <- meta[!(meta$SID %in% sampletidy), ]


### Create age transformation functions ------------------------------------

#### Relative age transformation for Uni Clk 2 -----------------------------

# Relative age transformation
fun_relage_trans <- Vectorize(function(age, Gestation, MaxAge) {
  if (is.na(age) || is.na(Gestation) || is.na(MaxAge)) {
    return(NA)
  }
  relativeAge <-
  (age + Gestation) / (MaxAge + Gestation)
  -log(-log(relativeAge))
  })

# Inverse age transformation
fun_relage_inv <- Vectorize(function(y, Gestation, MaxAge) {
  if (is.na(y) || is.na(Gestation) || is.na(MaxAge)) {
    return(NA)
  }
  relativeAge <- exp(-exp(-y))
  age <-
    relativeAge * (MaxAge + Gestation) - Gestation
  age
  })


#### Log/Linear age transformation for Uni Clk 3 and ENR -------------------

# Log age transformation
fun_llin_trans <- Vectorize(function(age, maturity) {
  if (is.na(age) || is.na(maturity)) return(NA)
  k <- 1.5
  if (age < maturity) {
    log((age + k) / (maturity + k))
  } else {(age - maturity) / (maturity + k)
  }
})


# Inverse age transformation
fun_llin_inv <- Vectorize(function(y, maturity) {
  if (is.na(y) || is.na(maturity)) return(NA)
  k <- 1.5
  if (y < 0) {
    (maturity + k) * exp(y) - k
  } else {(maturity + k) * y + maturity
  }
})


### Creating training and testing vectors ----------------------------------

knownAge <- !is.na(meta$Age)

beta_train <- beta[, knownAge]
beta_test  <- beta[, !knownAge]

meta_train <- meta[knownAge, ]
meta_test  <- meta[!knownAge, ]

age_train <- meta_train$Age

# Convert meta Age to log/linear age
meta_train$tAge <- fun_llin_trans(
  age = meta_train$Age,
  maturity = meta_train$AgeASM
)


### Create helper function for summary statistics and plot labels -----------

# Summary statistics
clock_summary_stats <- function(actual,
                             predicted,
                             individuals = NULL) {
  keep <- complete.cases(actual, predicted)
  actual <- actual[keep]
  predicted <- predicted[keep]
  if (!is.null(individuals))
    individuals <- individuals[keep]
  data.frame(
    n_samples = length(actual),
    n_individuals = if (is.null(individuals))
      NA_integer_
    else
      length(unique(individuals)),
    r = cor(actual, predicted),
    MAE = median(abs(predicted - actual)),
    MeanAE = mean(abs(predicted - actual)),
    RMSE = sqrt(mean((predicted - actual)^2))
  )
}

# Summary labels for plots
clock_summary_label <- function(actual, predicted, individuals = NULL) {
  stats <- clock_summary_stats(actual, predicted, individuals)
  paste0(
    "n = ", stats$n_samples, "(", stats$n_individuals, ")\n",
    "r = ", sprintf("%.2f", stats$r),
    ", MAE = ", sprintf("%.2f", stats$MAE)
  )
}


### Applying Universal Clocks 2 and 3 ---------------------------------------

beta_train_uniclks <- beta_train

#rownames(beta_train_uniclks) <- beta_train_uniclks$CGid
beta_train_uniclks$CGid <- NULL

intercept_uniclk2 <- uniclk$Coef.Universal2_RelativeAge[uniclk$var == "(Intercept)"]
intercept_uniclk3 <- uniclk$Coef.Universal3_Age.LogLinearRelAdult[uniclk$var == "(Intercept)"]

# Remove intercept row
coef_table <- subset(uniclk, var != "Intercept")

# Keep only common CpGs
common_cpgs <- intersect(rownames(beta_train_uniclks), uniclk$var)

length(common_cpgs)

beta_subset <- beta_train_uniclks[common_cpgs, , drop = FALSE]

coef_subset <- uniclk[match(common_cpgs, uniclk$var), ]

beta_matrix <- t(as.matrix(beta_subset))

### Universal Clock 2 -------------------------------------------------------

coef_uniclk2 <- coef_subset$Coef.Universal2_RelativeAge
coef_uniclk2[is.na(coef_uniclk2)] <- 0

pred_uniclk2_trans <- intercept_uniclk2 + beta_matrix %*% coef_uniclk2

pred_uniclk2_trans <- as.numeric(pred_uniclk2_trans)

pred_uniclk2 <- fun_relage_inv(
  y = pred_uniclk2_trans,
  Gestation = meta_train$Gestation,
  MaxAge = meta_train$MaxAge
)

# Store in metadata
meta_train$Universal2Age <- pred_uniclk2

#### Summary statistics Universal Clock 2 ----------------------------------- 

# Summary statistics Universal Clock 2 - All
summary_uniclk2 <- clock_summary_stats(
  actual = meta_train$Age,
  predicted = meta_train$Universal2Age,
  individuals = meta_train$AnimalID
)
summary_uniclk2

# Summary statistics Universal Clock 2 - Species
summary_uniclk2_species <- meta_train %>%
  group_by(SpeciesLatinName) %>%
  reframe(clock_summary_stats(
      actual = Age,
      predicted = Universal2Age,
      individuals = AnimalID
  ))
summary_uniclk2_species

# Summary statistics Universal Clock 2 - Tissue
summary_uniclk2_tissue <- meta_train %>%
  group_by(Tissue) %>%
  reframe(clock_summary_stats(
    actual = Age,
    predicted = Universal2Age,
    individuals = AnimalID
  ))
summary_uniclk2_tissue

# Summary statistics Universal Clock 2 - Species and Tissue
summary_uniclk2_species_tissue <- meta_train %>%
  group_by(SpeciesLatinName, Tissue) %>%
  reframe(clock_summary_stats(
    actual = Age,
    predicted = Universal2Age,
    individuals = AnimalID
  ))
summary_uniclk2_species_tissue


### Universal Clock 3 ------------------------------------------------------

coef_uniclk3 <- coef_subset$Coef.Universal3_Age.LogLinearRelAdult
coef_uniclk3[is.na(coef_uniclk3)] <- 0

pred_uniclk3_trans <- intercept_uniclk3 + beta_matrix %*% coef_uniclk3

pred_uniclk3_trans <- as.numeric(pred_uniclk3_trans)

pred_uniclk3 <- fun_llin_inv(
  y = pred_uniclk3_trans,
  maturity = meta_train$AgeASM
)

# Store in metadata
meta_train$Universal3Age <- pred_uniclk3

#### Summary statistics Universal Clock 3 ----------------------------------- 

# Summary statistics Universal Clock 3 - All
summary_uniclk3 <- clock_summary_stats(
  actual = meta_train$Age,
  predicted = meta_train$Universal3Age,
  individuals = meta_train$AnimalID
)
summary_uniclk3

# Summary statistics Universal Clock 2 - Species
summary_uniclk3_species <- meta_train %>%
  group_by(SpeciesLatinName) %>%
  reframe(clock_summary_stats(
    actual = Age,
    predicted = Universal3Age,
    individuals = AnimalID
  ))
summary_uniclk3_species

# Summary statistics Universal Clock 3 - Tissue
summary_uniclk3_tissue <- meta_train %>%
  group_by(Tissue) %>%
  reframe(clock_summary_stats(
    actual = Age,
    predicted = Universal3Age,
    individuals = AnimalID
  ))
summary_uniclk3_tissue

# Summary statistics Universal Clock 3 - Species and Tissue
summary_uniclk3_species_tissue <- meta_train %>%
  group_by(SpeciesLatinName, Tissue) %>%
  reframe(clock_summary_stats(
    actual = Age,
    predicted = Universal3Age,
    individuals = AnimalID
  ))
summary_uniclk3_species_tissue


### Create vectors for elastic net model construction -----------------------

#beta_train <- beta_train[, !(colnames(beta) %in% "CGid")]

common_samples <- intersect(colnames(beta_train), meta_train$SID)

beta_train <- beta_train[, common_samples]
meta_train <- meta_train[match(common_samples, meta_train$SID), ]

beta_data <- t(as.matrix(beta_train))
age_data <- meta_train$tAge

# Sanity check
stopifnot(
  nrow(beta_data) == length(age_data),
  all(rownames(beta_data) == meta_train$SID),
  !anyNA(age_data)
)


### Elastic net model -----------------------------------------------------

# Using alpha = 0.5 for elastic net regression
ALPH = 0.5

glmnet.Training.CV <- cv.glmnet(
  x = as(beta_data, "dgCMatrix"),
  y = as.numeric(age_data),
  nfolds = 10,
  alpha = ALPH,
  family = "gaussian"
)

lambda.training <- glmnet.Training.CV$lambda.1se


### LOIOCV -----------------------------------------------------------------

#Create column for DNAm
meta_train$DNAmoClockLOIO <- NA_real_

# Prepare for LOIOCV
individuals <- unique(meta_train$AnimalID)

for (id in individuals) {
  test_idx <- meta_train$AnimalID == id
  train_idx <- !test_idx
  beta_data_train <- beta_data[train_idx, ]
  beta_data_test  <- beta_data[test_idx, ]
  age_data_train <- age_data[train_idx]
  fit1 <- glmnet(
    x = beta_data_train,
    y = age_data_train,
    lambda = lambda.training,
    alpha = ALPH,
    family = "gaussian"
  )
  y.fitted <- predict(
    fit1,
    beta_data_test,
    type = "response"
  )
  meta_train$DNAmoClockLOIO[test_idx] <- as.numeric(y.fitted)
  print(id)
}

meta_train$PredictedAgeLOIOCV <- fun_llin_inv(
  y = meta_train$DNAmoClockLOIO,
  maturity = meta_train$AgeASM
)


#### Summary statistics Monotreme Epigenetic Clock ----------------------------------- 

# Summary statistics Monotreme Epigenetic Clock - All
summary_LOIOCV <- clock_summary_stats(
  actual = meta_train$Age,
  predicted = meta_train$PredictedAgeLOIOCV,
  individuals = meta_train$AnimalID
)
summary_LOIOCV

# Summary statistics Monotreme Epigenetic Clock- Species
summary_LOIOCV_species <- meta_train %>%
  group_by(SpeciesLatinName) %>%
  reframe(clock_summary_stats(
    actual = Age,
    predicted = PredictedAgeLOIOCV,
    individuals = AnimalID
  ))
summary_LOIOCV_species

# Summary statistics Monotreme Epigenetic Clock - Tissue
summary_LOIOCV_tissue <- meta_train %>%
  group_by(Tissue) %>%
  reframe(clock_summary_stats(
    actual = Age,
    predicted = PredictedAgeLOIOCV,
    individuals = AnimalID
  ))
summary_LOIOCV_tissue

# Summary statistics Monotreme Epigenetic Clock - Species and Tissue
summary_LOIOCV_species_tissue <- meta_train %>%
  group_by(SpeciesLatinName, Tissue) %>%
  reframe(clock_summary_stats(
    actual = Age,
    predicted = PredictedAgeLOIOCV,
    individuals = AnimalID
  ))
summary_LOIOCV_species_tissue


### Final LOIOCV Clock model ------------------------------------------------

beta_data_train <- t(beta_train)

age_data_train <- meta_train$tAge

final_model_LOIO <- glmnet(
  x = beta_data_train,
  y = age_data_train,
  lambda = lambda.training,
  alpha = ALPH,
  family = "gaussian"
)

beta_data_test <- t(beta_test)

meta_test$DNAmAge_t <- as.numeric(
  predict(
    final_model_LOIO,
    beta_data_test,
    type = "response"
  )
)

meta_test$PredictedAgeLOIO <- fun_llin_inv(
  y = meta_test$DNAmAge_t,
  maturity = meta_test$AgeASM
)


# CpG Count ----------------------------------------------------------------

coef_mat <- coef(glmnet.Training.CV, s = "lambda.1se")

sum(coef_mat != 0) # including intercept

sum(coef_mat[-1, ] != 0) # excluding intercept

# Extract CpG names

coef_df <- as.matrix(coef_mat)

selected_cpgs <- rownames(coef_df)[coef_df[,1] != 0]

selected_cpgs

selected_cpgs <- selected_cpgs[selected_cpgs != "(Intercept)"]

coef_table <- data.frame(
  CpG = rownames(coef_df),
  Coefficient = coef_df[,1]
)

coef_table <- subset(coef_table, Coefficient != 0)


### Assessing batch effects ------------------------------------------------

#Plat_C42012_sk_20_29052024_REPLICATE
#SIDs = ET083595 and ET0835103
cor(beta[, "ET083595"],
    beta[, "ET0835103"],
    method = "pearson")


#SBE_B50077_sk_80_28082025_REPLICATE
# SIDs =  ET083596 and ET0835104

cor(beta[, "ET083596"],
    beta[, "ET0835104"],
    method = "pearson")


### Analysis of treatment groups ----------------------------
# ANOVA test, not appropriate as data is not normally distributed
# Kruskal_Wallis test, not appropriate as samples not independent

treatment_df <- read.csv('Outputs_ENR/EpigeneticClock_KnownAgePredictions.csv')

group <- as.factor(treatment_df$Group)

#abs_error1 <- abs(errorLOOCV) # likely more appropriate to use to assess impact of accuracy on predicted ages
#abs_error2 <- abs(errorLOIOCV)

treatment_df$error <- treatment_df$Age - treatment_df$PredictedAgeLOIOCV
treatment_df$abs_error <- abs(treatment_df$Age - treatment_df$PredictedAgeLOIOCV)

# Treatment type

lmm_treatment1 <- lmer(
  error ~ group + Tissue + (1 | AnimalID),
  data = treatment_df
)

summary(lmm_treatment1)

coef_table_treatment1 <- as.data.frame(summary(lmm_treatment1)$coefficients)

lmm_treatment2 <- lmer(
  abs_error ~ group + Tissue + (1 | AnimalID),
  data = treatment_df
)

summary(lmm_treatment2)

coef_table_treatment2 <- as.data.frame(summary(lmm_treatment2)$coefficients)

# Sex

lmm_sex1 <- lmer(
  error ~ Sex + Tissue + (1 | AnimalID),
  data = treatment_df
)

summary(lmm_sex1)

coef_table_sex1 <- as.data.frame(summary(lmm_sex1)$coefficients)

lmm_sex2 <- lmer(
  abs_error ~ Sex + Tissue + (1 | AnimalID),
  data = treatment_df
)

summary(lmm_sex2)

coef_table_sex2 <- as.data.frame(summary(lmm_sex2)$coefficients)

# Species

lmm_species1 <- lmer(
  error ~ SpeciesLatinName + Tissue + (1 | AnimalID),
  data = treatment_df
)

summary(lmm_species1)

lmm_species2 <- lmer(
  abs_error ~ SpeciesLatinName + Tissue + (1 | AnimalID),
  data = treatment_df
)

summary(lmm_species2)
