
# July 2026

# Sarah Cordell et al.

### Getting setup -----------------------------------------------------------

library(tidyverse)
library(glmnet)
library(coefplot)
library(caret)
#library(lme4)
library(lmerTest)
library(ggtext)
library(patchwork)
library(grid)
library(readxl)


setwd("C:/Users/z5628859/UNSW/Platypus MSc - Sarah - General/003_Chapter 2/Chapter 2 Analysis/Chapter 2 Monotreme Epigenetics")

set.seed(1234567)


### Loading data ------------------------------------------------------------

# Load beta values 
beta <- read.csv('Data/data_sesame_40k.csv')
rownames(beta) <- beta$CGid
beta$CGid <- NULL

dim(beta)

# Load metadata (sample information)
meta <- read.csv('Data/ET0835.DayTarongaPlatypus.m40k_metadata.csv')

head(meta)

# Load master data sheet
# Add individual origin (wild/captive)

master <- read_excel('Chapter 2_DATA MASTER.xlsx')

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


#### Appendix S3. Summary of life history parameters -----------------------

write.csv(
  species_params,
  file = "Outputs_ENR/Species_Parameters.csv",
  row.names = FALSE
)


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


#### Table X Summary of samples analysed to develop clock model -------------

sample_summary1 <- meta %>%
  dplyr::mutate(
    KnownAge = !is.na(Age)
  ) %>%
  dplyr::group_by(
    SpeciesLatinName,
    KnownAge,
    Tissue,
    Sex
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    MinAge = ifelse(all(is.na(Age)), NA, min(Age, na.rm = TRUE)),
    MaxAge = ifelse(all(is.na(Age)), NA, max(Age, na.rm = TRUE)),
    MedianAge = ifelse(all(is.na(Age)), NA, median(Age, na.rm = TRUE)),
    .groups = "drop"
  )

write.csv(
  sample_summary1,
  file = "Outputs_ENR/Sample_Summary_Age_Sex.csv",
  row.names = FALSE
)


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


#### Appendix, potentially, duplicated info from S3 ------------------------

write.csv(meta |>
  dplyr::count(SpeciesLatinName, !is.na(meta$Age), Tissue),
  file = "Outputs_ENR/Sample_Summary.csv",
  row.names = FALSE
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

# Summary table Universal Clock 2 - All
write.csv(
  summary_uniclk2,
  file = "Outputs_ENR/Summary_Stats_UniversalClock2_All.csv",
  row.names = FALSE
)

# Summary statistics Universal Clock 2 - Species
summary_uniclk2_species <- meta_train %>%
  group_by(SpeciesLatinName) %>%
  reframe(clock_summary_stats(
      actual = Age,
      predicted = Universal2Age,
      individuals = AnimalID
  ))
summary_uniclk2_species

# Summary table Universal Clock 2 - Species
write.csv(
  summary_uniclk2_species,
  file = "Outputs_ENR/Summary_Stats_UniversalClock2_Species.csv",
  row.names = FALSE
)

# Summary statistics Universal Clock 2 - Tissue
summary_uniclk2_tissue <- meta_train %>%
  group_by(Tissue) %>%
  reframe(clock_summary_stats(
    actual = Age,
    predicted = Universal2Age,
    individuals = AnimalID
  ))
summary_uniclk2_tissue

# Summary table Universal Clock 2 - Tissue
write.csv(
  summary_uniclk2_tissue,
  file = "Outputs_ENR/Summary_Stats_UniversalClock2_Tissue.csv",
  row.names = FALSE
)

# Summary statistics Universal Clock 2 - Species and Tissue
summary_uniclk2_species_tissue <- meta_train %>%
  group_by(SpeciesLatinName, Tissue) %>%
  reframe(clock_summary_stats(
    actual = Age,
    predicted = Universal2Age,
    individuals = AnimalID
  ))
summary_uniclk2_species_tissue

# Summary table Universal Clock 2 - Species and Tissue
write.csv(
  summary_uniclk2_species_tissue,
  file = "Outputs_ENR/Summary_Stats_UniversalClock2_Species_Tissue.csv",
  row.names = FALSE
)

#### Plot Universal Clock 2 -------------------------------------------------

UniClock2_plot <- ggplot(meta_train, aes(x = Age, y = Universal2Age, 
                          colour = Species, shape = Tissue)) +
  geom_point(size = 2, stroke = 0.7) +
  geom_smooth(
    inherit.aes = FALSE,
    aes(x = Age, y = Universal2Age),
    method = "lm",
    se = FALSE,
    colour = "black",
    linewidth = 1
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_colour_manual(
    name = "Species",
    values = c(
      "Platypus" = "#293583",
      "Short-beaked echidna" = "#ff7254"
    )) +
  scale_shape_manual(
    name = "Tissue",
    values = c(
      "blood" = 21,
      "skin" = 24
    ),
    labels = c(
      "Blood",
      "Skin"
    )) +
  coord_equal(
    xlim = c(0, 40),
    ylim = c(0, 50)
  ) +
  theme_classic() +
  labs(
    x = "Chronological age (years)",
    y = "Predicted age (years)"
  )
UniClock2_plot

label_uniclk2 <- clock_summary_label(
  actual = meta_train$Age,
  predicted = meta_train$Universal2Age,
  individuals = meta_train$AnimalID
)

UniClock2_plot_labeled <- UniClock2_plot +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = label_uniclk2,
    hjust = 1.1,
    vjust = 1.1,
    size = 6,
    colour = "black"
  )
UniClock2_plot_labeled

ggsave(
  filename = "Outputs_ENR/UniClk2_plot.png",
  plot = UniClock2_plot_labeled,
  width = 6,
  height = 6,
  units = "in",
  dpi = 300
)

UniClock2_plot_clean <- UniClock2_plot_labeled +
  labs(
    x = NULL,
    y = "Predicted age (years)"
  )
UniClock2_plot_clean

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

# Summary table Universal Clock 2 - All
write.csv(
  summary_uniclk3,
  file = "Outputs_ENR/Summary_Stats_UniversalClock3_All.csv",
  row.names = FALSE
)

# Summary statistics Universal Clock 2 - Species
summary_uniclk3_species <- meta_train %>%
  group_by(SpeciesLatinName) %>%
  reframe(clock_summary_stats(
    actual = Age,
    predicted = Universal3Age,
    individuals = AnimalID
  ))
summary_uniclk3_species

# Summary table Universal Clock 3 - Species
write.csv(
  summary_uniclk3_species,
  file = "Outputs_ENR/Summary_Stats_UniversalClock3_Species.csv",
  row.names = FALSE
)

# Summary statistics Universal Clock 3 - Tissue
summary_uniclk3_tissue <- meta_train %>%
  group_by(Tissue) %>%
  reframe(clock_summary_stats(
    actual = Age,
    predicted = Universal3Age,
    individuals = AnimalID
  ))
summary_uniclk3_tissue

# Summary table Universal Clock 3 - Tissue
write.csv(
  summary_uniclk3_tissue,
  file = "Outputs_ENR/Summary_Stats_UniversalClock3_Tissue.csv",
  row.names = FALSE
)

# Summary statistics Universal Clock 3 - Species and Tissue
summary_uniclk3_species_tissue <- meta_train %>%
  group_by(SpeciesLatinName, Tissue) %>%
  reframe(clock_summary_stats(
    actual = Age,
    predicted = Universal3Age,
    individuals = AnimalID
  ))
summary_uniclk3_species_tissue

# Summary table Universal Clock 2 - Species and Tissue
write.csv(
  summary_uniclk3_species_tissue,
  file = "Outputs_ENR/Summary_Stats_UniversalClock3_Species_Tissue.csv",
  row.names = FALSE
)

#### Plot Universal Clock 3 -------------------------------------------------

UniClock3_plot <- ggplot(meta_train, aes(x = Age, y = Universal3Age, 
                                         colour = Species, shape = Tissue)) +
  geom_point(size = 2, stroke = 0.7) +
  geom_smooth(
    inherit.aes = FALSE,
    aes(x = Age, y = Universal3Age),
    method = "lm",
    se = FALSE,
    colour = "black",
    linewidth = 1
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_colour_manual(
    name = "Species",
    values = c(
      "Platypus" = "#293583",
      "Short-beaked echidna" = "#ff7254"
    )) +
  scale_shape_manual(
    name = "Tissue",
    values = c(
      "blood" = 21,
      "skin" = 24
    ),
    labels = c(
      "Blood",
      "Skin"
    )) +
  coord_equal(
    xlim = c(0, 40),
    ylim = c(0, 50)
  ) +
  theme_classic() +
  labs(
    x = "Chronological age (years)",
    y = "Predicted age (years)"
  )
UniClock3_plot

label_uniclk3 <- clock_summary_label(
  actual = meta_train$Age,
  predicted = meta_train$Universal3Age,
  individuals = meta_train$AnimalID
)

UniClock3_plot_labeled <- UniClock3_plot +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = label_uniclk2,
    hjust = 1.1,
    vjust = 1.1,
    size = 6,
    colour = "black"
  )
UniClock3_plot_labeled

ggsave(
  filename = "Outputs_ENR/UniClk3_plot.png",
  plot = UniClock3_plot_labeled,
  width = 6,
  height = 6,
  units = "in",
  dpi = 300
)

UniClock3_plot_clean <- UniClock3_plot_labeled +
  labs(
    x = NULL,
    y = NULL
  )
UniClock3_plot_clean


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
# Try lambda.min instead to see how accuracy is impacted
#lambda.training <- glmnet.Training.CV$lambda.min

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

#### Appendix S4. Summary of unique individuals ----------------------------

unique_individuals <- meta %>%
  count(AnimalID, Species, name = "n_samples") %>%
  filter(n_samples > 1) %>%
  arrange(desc(n_samples))

unique_individuals

write.csv(
  unique_individuals,
  file = "Outputs_ENR/Unique_Individual_IDs.csv",
  row.names = FALSE
)

unique_individuals1 <- meta_train %>%
  dplyr::group_by(SpeciesLatinName) %>%
  dplyr::summarise(
    UniqueIndividuals = dplyr::n_distinct(AnimalID),
    Samples = dplyr::n(),
    .groups = "drop"
  )

write.csv(
  unique_individuals1,
  file = "Outputs_ENR/Unique_Individuals.csv",
  row.names = FALSE
)

#### Appendix S5. Summary of unique individuals per tissue -----------------

unique_individuals2 <- meta_train %>%
  dplyr::group_by(SpeciesLatinName, Tissue) %>%
  dplyr::summarise(
    UniqueIndividuals = dplyr::n_distinct(AnimalID),
    Samples = dplyr::n(),
    .groups = "drop"
)

write.csv(
  unique_individuals2,
  file = "Outputs_ENR/Unique_Individuals_Species_Tissue.csv",
  row.names = FALSE
)

meta_test %>%
  dplyr::group_by(SpeciesLatinName) %>%
  dplyr::summarise(
    UniqueIndividuals = dplyr::n_distinct(AnimalID),
    Samples = dplyr::n(),
    .groups = "drop"
  )


#### Summary statistics Monotreme Epigenetic Clock ----------------------------------- 

# Summary statistics Monotreme Epigenetic Clock - All
summary_LOIOCV <- clock_summary_stats(
  actual = meta_train$Age,
  predicted = meta_train$PredictedAgeLOIOCV,
  individuals = meta_train$AnimalID
)
summary_LOIOCV

# Summary table Monotreme Epigenetic Clock - All
write.csv(
  summary_LOIOCV,
  file = "Outputs_ENR/Summary_Stats_MonotremeClock_All.csv",
  row.names = FALSE
)

# Summary statistics Monotreme Epigenetic Clock- Species
summary_LOIOCV_species <- meta_train %>%
  group_by(SpeciesLatinName) %>%
  reframe(clock_summary_stats(
    actual = Age,
    predicted = PredictedAgeLOIOCV,
    individuals = AnimalID
  ))
summary_LOIOCV_species

# Summary table Monotreme Epigenetic Clock - Species
write.csv(
  summary_LOIOCV_species,
  file = "Outputs_ENR/Summary_Stats_MonotremeClock_Species.csv",
  row.names = FALSE
)

# Summary statistics Monotreme Epigenetic Clock - Tissue
summary_LOIOCV_tissue <- meta_train %>%
  group_by(Tissue) %>%
  reframe(clock_summary_stats(
    actual = Age,
    predicted = PredictedAgeLOIOCV,
    individuals = AnimalID
  ))
summary_LOIOCV_tissue

# Summary table Monotreme Epigenetic Clock - Tissue
write.csv(
  summary_LOIOCV_tissue,
  file = "Outputs_ENR/Summary_Stats_MonotremeClock_Tissue.csv",
  row.names = FALSE
)

# Summary statistics Monotreme Epigenetic Clock - Species and Tissue
summary_LOIOCV_species_tissue <- meta_train %>%
  group_by(SpeciesLatinName, Tissue) %>%
  reframe(clock_summary_stats(
    actual = Age,
    predicted = PredictedAgeLOIOCV,
    individuals = AnimalID
  ))
summary_LOIOCV_species_tissue

# Summary table Monotreme Epigenetic Clock - Species and Tissue
write.csv(
  summary_LOIOCV_species_tissue,
  file = "Outputs_ENR/Summary_Stats_MonotremeClock_Species_Tissue.csv",
  row.names = FALSE
)

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

# Export results 
write.csv(
  meta_train,
  file = "Outputs_ENR/EpigeneticClock_KnownAgePredictions.csv",
  row.names = FALSE
)

write.csv(
  meta_test,
  file = "Outputs_ENR/EpigeneticClock_UnknownAgePredictions.csv",
  row.names = FALSE
)


#### Plot Monotreme Epigenetic Clock ----------------------------------------

MonotremeClock_plot <- ggplot(meta_train, aes(x = Age, y = PredictedAgeLOIOCV, 
                                         colour = Species, shape = Tissue)) +
  geom_point(size = 2, stroke = 0.7) +
  geom_smooth(
    inherit.aes = FALSE,
    aes(x = Age, y = PredictedAgeLOIOCV),
    method = "lm",
    se = FALSE,
    colour = "black",
    linewidth = 1
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_colour_manual(
    name = "Species",
    values = c(
      "Platypus" = "#293583",
      "Short-beaked echidna" = "#ff7254"
    )) +
  scale_shape_manual(
    name = "Tissue",
    values = c(
      "blood" = 21,
      "skin" = 24
    ),
    labels = c(
      "Blood",
      "Skin"
    )) +
  coord_equal(
    xlim = c(0, 40),
    ylim = c(0, 50)
  ) +
  theme_classic() +
  labs(
    x = "Chronological age (years)",
    y = "Predicted age (years)"
  )
MonotremeClock_plot

label_MonotremeClock <- clock_summary_label(
  actual = meta_train$Age,
  predicted = meta_train$PredictedAgeLOIOCV,
  individuals = meta_train$AnimalID
)

MonotremeClock_plot_labeled <- MonotremeClock_plot +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = label_MonotremeClock,
    hjust = 1.1,
    vjust = 1.1,
    size = 6,
    colour = "black"
  )
MonotremeClock_plot_labeled

ggsave(
  filename = "Outputs_ENR/MonoClk_plot.png",
  plot = MonotremeClock_plot_labeled,
  width = 6,
  height = 6,
  units = "in",
  dpi = 300
)

MonotremeClock_plot_clean <- MonotremeClock_plot_labeled +
labs(
  x = NULL,
  y = NULL
)
MonotremeClock_plot_clean


#### Plot Species - Monotreme Epigenetic Clock ------------------------------

MonoClk_species_plot <- ggplot(meta_train, aes(x = Age, y = PredictedAgeLOIOCV, 
                                              colour = Species, shape = Tissue)) +
  geom_point(size = 2, stroke = 0.8) +
  geom_smooth(
    inherit.aes = FALSE,
    aes(x = Age, y = PredictedAgeLOIOCV),
    method = "lm",
    se = FALSE,
    colour = "black",
    linewidth = 1
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  facet_wrap( ~ Species) +
  scale_colour_manual(
    name = "Species",
    values = c(
      "Platypus" = "#293583",
      "Short-beaked echidna" = "#ff7254"
    )) +
  scale_shape_manual(
    name = "Tissue",
    values = c(
      "blood" = 21,
      "skin" = 24
    ),
    labels = c(
      "Blood",
      "Skin"
    )) +
  coord_equal(
    xlim = c(0, 40),
    ylim = c(0, 50)
  ) +
  theme_classic() +
  labs(
    x = "Chronological age (years)",
    y = "Predicted age (years)"
  ) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14)
  )
MonoClk_species_plot

label_MonoClk_species <- meta_train %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    label = clock_summary_label(
      actual = Age,
      predicted = PredictedAgeLOIOCV,
      individuals = AnimalID
    ),
    .groups = "drop"
  )

MonoClk_species_plot_labeled <- MonoClk_species_plot +
  geom_text(
    data = label_MonoClk_species,
    aes(label = label),
    x = 39,
    y = 49,
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 1,
    size = 5
  )
MonoClk_species_plot_labeled

# Save plot
ggsave(
  filename = "Outputs_ENR/MonoClk_plot_species.png",
  plot = MonoClk_species_plot_labeled,
  width = 10,
  height = 4,
  units = "in",
  dpi = 300
)

MonoClk_species_plot_clean <- MonoClk_species_plot_labeled +
  labs(
    x = NULL,
    y = "Predicted age (years)"
  )
MonoClk_species_plot_clean


#### Plot Tissue - Monotreme Epigenetic Clock ------------------------------

MonoClk_tissue_plot <- ggplot(meta_train, aes(x = Age, y = PredictedAgeLOIOCV, 
                                               colour = Species, shape = Tissue)) +
  geom_point(size = 2, stroke = 0.8) +
  geom_smooth(
    inherit.aes = FALSE,
    aes(x = Age, y = PredictedAgeLOIOCV),
    method = "lm",
    se = FALSE,
    colour = "black",
    linewidth = 1
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  facet_wrap( ~ Tissue,
              labeller = labeller(
                Tissue = c(
                  "blood" = "Blood",
                  "skin" = "Skin"
                ))
  ) +
  scale_colour_manual(
    name = "Species",
    values = c(
      "Platypus" = "#293583",
      "Short-beaked echidna" = "#ff7254"
    )) +
  scale_shape_manual(
    name = "Tissue",
    values = c(
      "blood" = 21,
      "skin" = 24
    ),
    labels = c(
      "Blood",
      "Skin"
    )) +
  coord_equal(
    xlim = c(0, 40),
    ylim = c(0, 50)
  ) +
  theme_classic() +
  labs(
    x = "Chronological age (years)",
    y = "Predicted age (years)"
  ) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14)
  )
MonoClk_tissue_plot

label_MonoClk_tissue <- meta_train %>%
  dplyr::group_by(Tissue) %>%
  dplyr::summarise(
    label = clock_summary_label(
      actual = Age,
      predicted = PredictedAgeLOIOCV,
      individuals = AnimalID
    ),
    .groups = "drop"
  )

MonoClk_tissue_plot_labeled <- MonoClk_tissue_plot +
  geom_text(
    data = label_MonoClk_tissue,
    aes(label = label),
    x = 39,
    y = 49,
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 1,
    size = 5
  )
MonoClk_tissue_plot_labeled

# Save plot
ggsave(
  filename = "Outputs_ENR/MonoClk_plot_tissue.png",
  plot = MonoClk_tissue_plot_labeled,
  width = 10,
  height = 4,
  units = "in",
  dpi = 300
)

MonoClk_tissue_plot_clean <- MonoClk_tissue_plot_labeled +
  labs(
    x = NULL,
    y = NULL
  )
MonoClk_tissue_plot_clean


#### Plot Species and Tissue - Monotreme Epigenetic Clock ------------------

MonoClk_species_tissue_plot <- ggplot(meta_train, aes(x = Age, y = PredictedAgeLOIOCV, 
                                              colour = Species, shape = Tissue)) +
  geom_point(size = 2, stroke = 0.8) +
  geom_smooth(
    inherit.aes = FALSE,
    aes(x = Age, y = PredictedAgeLOIOCV),
    method = "lm",
    se = FALSE,
    colour = "black",
    linewidth = 1
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  facet_wrap( ~ Species + Tissue, nrow = 1,
              labeller = labeller(
                Tissue = c(
                  "blood" = "Blood",
                  "skin" = "Skin"
                ))
  ) +
  scale_colour_manual(
    name = "Species",
    values = c(
      "Platypus" = "#293583",
      "Short-beaked echidna" = "#ff7254"
    )) +
  scale_shape_manual(
    name = "Tissue",
    values = c(
      "blood" = 21,
      "skin" = 24
    ),
    labels = c(
      "Blood",
      "Skin"
    )) +
  coord_equal(
    xlim = c(0, 40),
    ylim = c(0, 50)
  ) +
  theme_classic() +
  labs(
    x = "Chronological age (years)",
    y = "Predicted age (years)"
  ) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14)
  )
MonoClk_species_tissue_plot

label_MonoClk_species_tissue <- meta_train %>%
  dplyr::group_by(Species, Tissue) %>%
  dplyr::summarise(
    label = clock_summary_label(
      actual = Age,
      predicted = PredictedAgeLOIOCV,
      individuals = AnimalID
    ),
    .groups = "drop"
  )

MonoClk_species_tissue_plot_labeled <- MonoClk_species_tissue_plot +
  geom_text(
    data = label_MonoClk_species_tissue,
    aes(label = label),
    x = 39,
    y = 49,
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 1,
    size = 5
  )
MonoClk_species_tissue_plot_labeled

# Save plot
ggsave(
  filename = "Outputs_ENR/MonoClk_plot_species_tissue.png",
  plot = MonoClk_species_tissue_plot_labeled,
  width = 10,
  height = 4,
  units = "in",
  dpi = 300
)

MonoClk_species_tissue_plot_clean <- MonoClk_species_tissue_plot_labeled +
  labs(
    x = NULL,
    y = "Predicted Age (years)"
  )
MonoClk_species_tissue_plot_clean

#### Figure 1. Comparison of the three clocks -----------------------------

Figure_1 <- (UniClock2_plot_clean | UniClock3_plot_clean | MonotremeClock_plot_clean) +
  plot_layout(guides = "collect") +
  plot_annotation(
    caption = "Chronological age (years)",
    tag_levels = list(c("(a)", "(b)", "(c)"))
  ) &
  theme(
    legend.position = "right",
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.caption = element_text(
      hjust = 0.5,
      size = 18
    ),
    plot.tag = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )
Figure_1

ggsave(
  filename = "Outputs_ENR/Chapter 2_Figure 1.png",
  plot = Figure_1,
  width = 14,
  height = 7,
  units = "in",
  dpi = 300
)


#### Figure 2. Monotreme Clock filtered ------------------------------------

Figure_2 <- (MonoClk_species_plot_clean | MonoClk_tissue_plot_clean) / 
  (MonoClk_species_tissue_plot_clean) +
  plot_layout(guides = "collect") +
  plot_annotation(
    caption = "Chronological age (years)", 
    tag_levels = list(c("(a)", "(b)", "(c)"))
  ) &
  theme(
    legend.position = "right",
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.caption = element_text(
      hjust = 0.5,
      size = 18
    ),
    plot.tag = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )
Figure_2

ggsave(
  filename = "Outputs_ENR/Chapter 2_Figure 2.png",
  plot = Figure_2,
  width = 14,
  height = 8,
  units = "in",
  dpi = 300
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

# Table for Supplementary data

write.csv(
  coef_table,
  file = "Outputs_ENR/CpGs_Retained.csv",
  row.names = FALSE
)    
    
# ????
results$nCpGs <- sapply(results$alpha, function(a) {
  cvfit <- cv.glmnet(
    x = as(beta_data, "dgCMatrix"),
    y = as.numeric(age_data),
    alpha = a,
    family = "gaussian"
  )
  coef_mat <- coef(cvfit, s = "lambda.1se")
  sum(coef_mat[-1, ] != 0)
})

alpha | lambda | nCpGs | MAE | r


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

write.csv(
  coef_table_treatment1,
  file = "Outputs_ENR/Storage treatment_lmm.csv",
  row.names = TRUE
)

lmm_treatment2 <- lmer(
  abs_error ~ group + Tissue + (1 | AnimalID),
  data = treatment_df
)

summary(lmm_treatment2)

coef_table_treatment2 <- as.data.frame(summary(lmm_treatment2)$coefficients)

write.csv(
  coef_table_treatment2,
  file = "Outputs_ENR/Storage treatment_abserror_lmm.csv",
  row.names = TRUE
)

# Individuals' origin

lmm_origin1 <- lmer(
  error ~ individual.origin + Tissue + (1 | AnimalID),
  data = treatment_df
)

summary(lmm_origin1)

coef_table_origin1 <- as.data.frame(summary(lmm_origin1)$coefficients)

write.csv(
  coef_table_origin1,
  file = "Outputs_ENR/Individuals origin_lmm.csv",
  row.names = TRUE
)

lmm_origin2 <- lmer(
  abs_error ~ individual.origin + Tissue + (1 | AnimalID),
  data = treatment_df
)

summary(lmm_origin2)

coef_table_origin2 <- as.data.frame(summary(lmm_origin2)$coefficients)

write.csv(
  coef_table_origin2,
  file = "Outputs_ENR/Individuals origin_abserror_lmm.csv",
  row.names = TRUE
)

# Sex

lmm_sex1 <- lmer(
  error ~ Sex + Tissue + (1 | AnimalID),
  data = treatment_df
)

summary(lmm_sex1)

coef_table_sex1 <- as.data.frame(summary(lmm_sex1)$coefficients)

write.csv(
  coef_table_sex1,
  file = "Outputs_ENR/Sex_lmm.csv",
  row.names = TRUE
)

lmm_sex2 <- lmer(
  abs_error ~ Sex + Tissue + (1 | AnimalID),
  data = treatment_df
)

summary(lmm_sex2)

coef_table_sex2 <- as.data.frame(summary(lmm_sex2)$coefficients)

write.csv(
  coef_table_sex2,
  file = "Outputs_ENR/Sex_abserror_lmm.csv",
  row.names = TRUE
)

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


# Table for supplementary data

treatment_summary <- meta_train %>%
  dplyr::group_by(SpeciesLatinName, Group, Tissue) %>%
  dplyr::summarise(
    Samples = dplyr::n(),
    .groups = "drop"
  )

write.csv(
  treatment_summary,
    file = "Outputs_ENR/Treatment_Summary.csv",
  row.names = FALSE
)


### Splitting ASM by sex ----------------------------------------------------

matAGEfem <- 7
matAGEmale <- 7
############# sexual maturity

SplitBySex <- split(x = info, f = info$Sex)
### Applies the LLin3 transformation to the input vector x
fun_llin_trans_Female <- Vectorize(function(x, maturity = matAGEfem) {
  if (is.na(x) | is.na(maturity)) {return(NA)}
  k <- 1.5
  y <- 0
  if (x < maturity) {y = log((x+k)/(maturity+k))}
  else {y = (x-maturity)/(maturity+k)}
  return(y)
})

fun_llin_trans_Male <- Vectorize(function(x, maturity = matAGEmale) {
  if (is.na(x) | is.na(maturity)) {return(NA)}
  k <- 1.5
  y <- 0
  if (x < maturity) {y = log((x+k)/(maturity+k))}
  else {y = (x-maturity)/(maturity+k)}
  return(y)
})

SplitBySex$Male$tAge = fun_llin_trans_Male(SplitBySex$Male$Age)
SplitBySex$Female$tAge = fun_llin_trans_Female(SplitBySex$Female$Age)


#inverse sex-specific

fun_llin_inv_Female <- Vectorize(function(y, maturity = matAGEfem) {
  if (is.na(y) | is.na(maturity)) {return(NA)}
  k <- 1.5
  x <- 0
  if (y < 0) {x = (maturity+k)*exp(y)-k}
  else {x = (maturity+k)*y+maturity}
  return(x)
})

fun_llin_inv_Male <- Vectorize(function(y, maturity = matAGEmale) {
  if (is.na(y) | is.na(maturity)) {return(NA)}
  k <- 1.5
  x <- 0
  if (y < 0) {x = (maturity+k)*exp(y)-k}
  else {x = (maturity+k)*y+maturity}
  return(x)
})

SplitBySex$Female$iAge = fun_llin_inv_Female(SplitBySex$Female$tAge)
SplitBySex$Male$iAge = fun_llin_inv_Male(SplitBySex$Male$tAge)

info <- unsplit(SplitBySex, info$Sex)


# Choosing alpha -----------------------------------------------------------

alphas <- seq(0.1, 1, by = 0.1)

alpha_results <- data.frame()

for (a in alphas) {
  cvfit <- cv.glmnet(
    x = as(beta_data, "dgCMatrix"),
    y = as.numeric(age_data),
    alpha = a,
    family = "gaussian",
    nfolds = 10
  )
  idx_1se <- which(cvfit$lambda == cvfit$lambda.1se)
  alpha_results <- rbind(
    alpha_results,
    data.frame(
      alpha = a,
      lambda = cvfit$lambda.1se,
      cvm = cvfit$cvm[idx_1se],
      nCpGs = sum(coef(cvfit, s = "lambda.1se")[-1, ] != 0)
    )
  )
}

alpha_results# choose lowest cvm (cross-validation error), = 0.6, however models 
# perform best overall with alpha = 0.5

# Summary alpha review ------------------------------------------------------
write.csv(
  alpha_results,
  file = "Outputs_ENR/Summary_alpha_review.csv",
  row.names = FALSE
)

