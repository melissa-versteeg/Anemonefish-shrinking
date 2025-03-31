

#---
title: "Individual clown anemonefish shrink to survive heat stress and social conflict"
authors: "Melissa Versteeg, Chancey MacDonald, Morgan Bennett-Smith, Peter Buston, Theresa Rueger"
year: '2025'

#---
  
# The code has been set up to align with the structure of the hypotheses. It contains the following sections:
#   
# Section 1: Libraries - With session and package information.
# Section 2: Data sets, transformations, subsetting, and data projections
# Section 3: Model results, plots, model fit, and effect extractions
#   Result 1a: Clown anemonefish shrinking is common during a heat stress event
#   Result 1b: Clown anemonefish growth and shrinking is related to size and social rank      
#   Result 2: Clown anemonefish growth and shrinking are linked
#   Result 3: Clown anemonefish growth is predicted by heat stress
#   Result 4: Clown anemonefish shrinking is also predicted by social conflict
#   Result 5: Clown anemonefish that shrink more often have higher chances of survival
# Section 4: Figures
#   Figure panel 1: Change in total length (TL) of clown anemonefish during a heat stress event.
#   Fig. 1A: The distribution of percent body length changes of individual rank 1 and rank 2 clown anemonefish.
#   Fig. 1B: The percent change in total length modelled in relation to initial total length of fish for rank 1.
#   Fig. 1C: The percent change in total length modelled in relation to initial total length of fish for rank 2.     
#   Figure 2: Percent change in total length (TL) for clown anemonefish with different shrinking response patterns.
#   Figure panel 3: Change in body length of clown anemonefish in relation to heat stress.
#   Fig. 3A: The percentage change in total length fitted against temperature exposure in the lunar month of length measurements for rank 1 and rank 2.
#   Fig. 3B: The distribution of the corresponding slope estimates for each rank, for current heat stress model.
#   Fig. 3C: The distribution of the corresponding slope contrasts between ranks for current heat stress model.
#   Fig. 3D: The percentage change in total length was fitted against temperature exposure in the preceding lunar month.
#   Fig. 3E: The distribution of the corresponding slope estimates for each rank, for preceding heat stress model.
#   Fig. 3F: The distribution of the corresponding slope contrasts between ranks for preceding heat stress model.  
#   Figure panel 4: Shrinking likelihood in relation to size ratio of both clown anemonefish within a breeding pair. 
#   Fig. 4A: The likelihood of shrinking was fitted against size ratios of rank 2 to rank 1 at the start of the lunar month.
#   Fig. 4B: The distribution of model estimates for shrinking probability among rank 1 and rank 2 individuals, at three discrete size ratios 0.6, 0.8, and 1.0
#   Fig. 4C: The corresponding contrasts in shrinking probabilities between the two ranks
#   Figure panel 5: Survival of clown anemonefish. 
#   Fig. 5A: Kaplan-Meier survival curves for heat stress at the anemone scale. 
#   Fig. 5B: Kaplan-Meier survival curves for growth patterns based on shrinking frequency. 
#   Fig. 5C: Kaplan-Meier survival curves for social factors of shrinking within a breeding pair.
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Section 1: Libraries ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
# R version 4.4.1 (2024-06-14 ucrt) -- "Race for Your Life"

library(tidyverse)    # ggplot, dplyr, %>%, and friends
library(readr)        # Importing csv files into R environment
library(brms)         # Bayesian regression modeling through Stan
library(tidybayes)    # Manipulate Stan objects in a tidy way
library(modelr)       # Modelling Functions that work with the Pipe
library(coda)         # Bayesian model convergence statistics and plotting functions
library(loo)          # Comparisons of Bayesian model fit parameters
library(emmeans)      # Calculate marginal effects in even fancier ways
library(posterior)    # Bayesian model functions for summarizing, visualising, and manipulating posterior distributions
library(bayesplot)    # Bayesian model assessment and diagnostic visualisation, using posterior distributions and MCMC diagnostics
library(bayestestR)   # Bayesian model analysis and interpretation, including calculation of credible intervals
library(ggeffects)    # Functions to compute and visualise predicted effects
library(parameters)   # Tools to extract and interpret effect sizes, confidence intervals and p-values
library(ggplot2)      # Plots and visualisation customisation
library(patchwork)    # Combine ggplot objects
library(forecast)     # ARIMA fore- and hindcasting of environmental data 
library(tree)         # Classification and regression trees
library(survival)     # Survival analysis
library(survminer)    # Drawing survival curves in ggplot2

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Section 2: Data import, transformations, subsetting, glimpses, and data projections ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Multiple observations of change in body length per anemonefish over entire study period ####
### Read in datasets ####
# Data: PercentChange_Versteeg_etal_2025
R1R2.bmrs <- read.csv("PercentChange_Versteeg_etal_2025.csv")

### Scaling of predictors in models ####
R1R2.bmrs$Init_TL_scl<-scale(R1R2.bmrs$Init_TL, center=TRUE, scale=TRUE)
R1R2.bmrs$M_Temp_LC_scl<-scale(R1R2.bmrs$M_Temp_LC, center=TRUE, scale=TRUE)
R1R2.bmrs$M_Temp_LC_INIT_scl<-scale(R1R2.bmrs$M_Temp_LC_INIT, center=TRUE, scale=TRUE)
R1R2.bmrs$Max_Temp_LC_scl<-scale(R1R2.bmrs$Max_Temp_LC, center=TRUE, scale=TRUE)
R1R2.bmrs$Max_Temp_LC_INIT_scl<-scale(R1R2.bmrs$Max_Temp_LC_INIT, center=TRUE, scale=TRUE)
R1R2.bmrs$Dev_Temp_LC_scl<-scale(R1R2.bmrs$Dev_Temp_LC, center=TRUE, scale=TRUE)
R1R2.bmrs$Dev_Temp_LC_INIT_scl<-scale(R1R2.bmrs$Dev_Temp_LC_INIT, center=TRUE, scale=TRUE)

### Adjusting variables for models ####
### Set as factors:
R1R2.bmrs$FishID <- as.factor(R1R2.bmrs$FishID)
R1R2.bmrs$GroupID <- as.factor(R1R2.bmrs$GroupID)
R1R2.bmrs$Rank <- as.factor(R1R2.bmrs$Rank)
R1R2.bmrs$LC_Adj <- as.factor(R1R2.bmrs$LC_Adj)
R1R2.bmrs$Shrink <- as.factor(R1R2.bmrs$Shrink)
R1R2.bmrs$Growth_ThreeGroups <- as.factor(R1R2.bmrs$Growth_ThreeGroups)

### Selecting of complete cases and subsets for rank ####
consistR1R2 <- R1R2.bmrs[complete.cases(R1R2.bmrs[, c("Perc_change_TL", "Dev_Temp_LC_scl", "Dev_Temp_LC_INIT_scl", 
                                                      "M_Temp_LC_scl", "M_Temp_LC_INIT_scl",
                                                      "Max_Temp_LC_scl", "Max_Temp_LC_INIT_scl", "Rank", "Shrink",
                                                      "Init_TL_scl")]), ]

### Rank 1 anemonefish:
R1 <- subset(consistR1R2, Rank == 1)
### Rank 2 anemonefish:
R2 <- subset(consistR1R2, Rank == 2) 

## Single observations of change in body length per anemonefish over entire study period ####
### Read in datasets ####
# Data: Survival_Versteeg_etal_2025.xlsx
R1R2.complete.g <- read.csv("Survival_Versteeg_etal_2025.csv")
R1R2.complete.h <- subset(R1R2.complete.g, GroupID != "B003")

### Adjusting variables for models ####
### Set as factors and scales:
R1R2.complete.h$GroupID <- as.factor(R1R2.complete.h$GroupID)
R1R2.complete.h$Rank <- as.factor(R1R2.complete.h$Rank)
R1R2.complete.h$LC_Adj <- as.factor(R1R2.complete.h$LC_Adj)
R1R2.complete.h$Br_Pair_Shr <- as.factor(R1R2.complete.h$Br_Pair_Shr)
R1R2.complete.h$Event_Died <- as.factor(R1R2.complete.h$Event_Died)
R1R2.complete.h$Growth_ThreeGroups <- as.factor(R1R2.complete.h$Growth_ThreeGroups)

R1R2.complete.h$MinOverall_FishID_scl <- scale(R1R2.complete.h$MinOverall_FishID, center=TRUE, scale=TRUE)
R1R2.complete.h$Init_TL_scl <- scale(R1R2.complete.h$Init_TL, center=TRUE, scale=TRUE)

## Multiple observations shrinking in body length over the entire study period ####
### Read in datasets ####
# Excel file: R1andR2Shrinking_Versteeg_etal_2025
binomR1R2.b <- read.csv("R1andR2Shrinking_Versteeg_etal_2025.csv")

### Adjusting variables for models ####
### Set as factors and scales:
binomR1R2.b$FishID <- as.factor(binomR1R2.b$FishID)
binomR1R2.b$GroupID <- as.factor(binomR1R2.b$GroupID)
binomR1R2.b$Rank <- as.factor(binomR1R2.b$Rank)
binomR1R2.b$LC_Adj <- as.factor(binomR1R2.b$LC_Adj)
binomR1R2.b$Shrinking <- as.factor(binomR1R2.b$Shrinking)

binomR1R2.b$Init_TL_scl <- scale(binomR1R2.b$Init_TL, center=TRUE, scale=TRUE)
binomR1R2.b$Size_R2R1_Initial_scl <- scale(binomR1R2.b$Size_R2R1_Initial, center=TRUE, scale=TRUE)
binomR1R2.b$Max_Temp_LC_INIT_scl <- scale(binomR1R2.b$Max_Temp_LC_INIT, center=TRUE, scale=TRUE)

### Select cases complete:
binomR1R2.c <- binomR1R2.b[complete.cases(binomR1R2.b[, c("Shrinking", "Size_R2R1_Initial_scl", "Max_Temp_LC_INIT_scl", "Init_TL_scl", 
                                                    "Rank","FishID", "GroupID", "LC_Adj")]), ]

## Temperature hindcasting and forecasting ####
### Read in datasets ####
# Excel file: ARIMA_Versteeg_etal_2025
R1R2_b <- read_csv("ARIMA_Versteeg_etal_2025.csv")

### Selecting of complete cases ####
Forecasting <- select(R1R2_b, "GroupID", "LC_Adj",  "Dev_Temp_LC", "Dev_Temp_LC_INIT", 
                      "M_Temp_LC", "M_Temp_LC_INIT",
                      "Max_Temp_LC", "Max_Temp_LC_INIT")

### Factors:
Forecasting$GroupID <- as.factor(Forecasting$GroupID)

### Maximum temperature exposure ####
### Visual inspection:
ggplot(data = Forecasting, aes(x = LC_Adj, y = Max_Temp_LC)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ GroupID) +
  ggtitle("Maximum temperature trends across lunar months 1-4 for each anemone") +
  theme_minimal()

### Temperature projections: Assessment of linear or polynomial approach for ARIMA
### Fit linear regression model: linear model to check best fit for temperature projection.
linear_model <- lm(Max_Temp_LC ~ LC_Adj, data = Forecasting)
summary(linear_model)

### Fit polynomial regression model: linear model to check best fit for temperature projection.
poly_model <- lm(Max_Temp_LC ~ poly(LC_Adj, degree = 2), data = Forecasting)
summary(poly_model)

### Comparison of linear_model and poly_model summary: Polynomial best approach based on Adjusted R-squared, 
# Residual Standard Error, and F-statistic.

### Check Residuals #### 
par(mfrow=c(2,2))
plot(linear_model, which=1:4, main="Linear Model Residuals")
plot(poly_model, which=1:4, main="Polynomial Model Residuals")
linear_resid <- residuals(linear_model)
poly_resid <- residuals(poly_model)

par(mfrow=c(1,2))
hist(linear_resid, main="Linear Model Residuals", xlab="Residuals", col="blue", breaks=30)
hist(poly_resid, main="Polynomial Model Residuals", xlab="Residuals", col="red", breaks=30)

### LC_Adj as numeric:
Forecasting$LC_Adj <- as.numeric(as.character(Forecasting$LC_Adj))
grouped_data <- split(Forecasting, Forecasting$GroupID)

### Polynomial model per GroupID:
fit_poly <- function(group_data) {
# Fit a polynomial regression model:
  poly_model <- lm(Max_Temp_LC ~ poly(LC_Adj, degree = 2), data = group_data)
  return(poly_model)
}

poly_models <- lapply(grouped_data, fit_poly)

### LC_Adj == 0 ####
## Predict Max_Temp_LC for LC_Adj == 0 for each GroupID:
LC_Adj_0_forecasts.Max <- lapply(poly_models, function(model) {
  # Predict Max_Temp_LC for LC_Adj == 0 using the polynomial regression model
  LC_Adj_0_prediction <- predict(model, newdata = data.frame(LC_Adj = 0))
  return(LC_Adj_0_prediction)
})

## Long format:
forecast_df_adj_0.Max <- map_dfr(names(LC_Adj_0_forecasts.Max), ~{
  data.frame(GroupID = .x, Max_Temp_LC = LC_Adj_0_forecasts.Max[[.x]])
})

print(forecast_df_adj_0.Max)

## Accuracy measure per GroupID:
calculate_accuracy <- function(group_data, forecast_data) {

merged_data <- merge(group_data, forecast_data, by = "GroupID", suffixes = c(".actual", ".forecast"))

accuracy_measures <- accuracy(merged_data$Max_Temp_LC.forecast, merged_data$Max_Temp_LC.actual)
  return(accuracy_measures)
}

accuracy_results <- lapply(names(grouped_data), function(group_id) {
  group_actual_data <- grouped_data[[group_id]]
  group_forecast_data <- forecast_df_adj_0.Max %>% filter(GroupID == group_id)
  accuracy_measures <- calculate_accuracy(group_actual_data, group_forecast_data)
  return(data.frame(GroupID = group_id, accuracy_measures))
})

## Combine:
accuracy_df_LC0.Max <- bind_rows(accuracy_results)
print(accuracy_df_LC0.Max)

### LC_Adj == 5 ####
## Predict Max_Temp_LC for LC_Adj == 5 for each GroupID:
LC_Adj_5_forecasts.Max <- lapply(poly_models, function(model) {
  # Predict Max_Temp_LC for LC_Adj == 5 using the polynomial regression model
  LC_Adj_5_prediction <- predict(model, newdata = data.frame(LC_Adj = 5))
  return(LC_Adj_5_prediction)
})

## Long format:
forecast_df_adj_5.Max <- map_dfr(names(LC_Adj_5_forecasts.Max), ~{
  data.frame(GroupID = .x, Max_Temp_LC = LC_Adj_5_forecasts.Max[[.x]])
})

print(forecast_df_adj_5.Max)

## Accuracy measures for each GroupID:
calculate_accuracy <- function(group_data, forecast_data) {
 
merged_data <- merge(group_data, forecast_data, by = "GroupID", suffixes = c(".actual", ".forecast"))

accuracy_measures <- accuracy(merged_data$Max_Temp_LC.forecast, merged_data$Max_Temp_LC.actual)
  return(accuracy_measures)
}

accuracy_results <- lapply(names(grouped_data), function(group_id) {
  group_actual_data <- grouped_data[[group_id]]
  group_forecast_data <- forecast_df_adj_5.Max %>% filter(GroupID == group_id)
  accuracy_measures <- calculate_accuracy(group_actual_data, group_forecast_data)
  return(data.frame(GroupID = group_id, accuracy_measures))
})

## Combine: 
accuracy_df_LC5.Max <- bind_rows(accuracy_results)
print(accuracy_df_LC5.Max)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Section 3: Model outputs, plots, model fit, and effect estimates ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Result 1a: Clown anemonefish shrinking is common during a heat stress event ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Descriptives of shrinking presence and patterns of shrinking events- considering rank and intitial total length of the anemonefish ####

### Observed shrinking among individuals:
table(R1R2.complete.h$Growth_ThreeGroups)
table(R1R2.complete.h$Rank, R1R2.complete.h$Growth_ThreeGroups)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Result 1b: Clown anemonefish growth and shrinking is related to size and social rank ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Regression trees for initial size threshold for percent change in size - considering rank of the anemonefish ####
### R1 only:
attach(R1)
threshR1 <- tree(Perc_change_TL ~ Init_TL)
print(threshR1)
tree_R1 <- tree(Perc_change_TL ~ Init_TL, data = R1)
plot(tree_R1, type = "uniform")
text(tree_R1, pretty = 0, cex = 0.8)
title("Regression tree for size thresholds and mean change in TL for rank 1 percula")
detach(R1)

### R2 only:
attach(R2)
threshR2 <- tree(Perc_change_TL ~ Init_TL)
print(threshR2)
tree_R2 <- tree(Perc_change_TL ~ Init_TL, data = R2)
plot(tree_R2, type = "uniform")
text(tree_R2, pretty = 0, cex = 0.8)
title("Regression tree for size thresholds and mean change in TL for rank 2 percula")
detach(R2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Result 2: Clown anemonefish growth and shrinking are linked ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Bayesian models of growth patterns- considering initial size and rank of the anemonefish and heat stress ####
### Growth pattern model 1a: 
Model.1a. <- brm(Perc_change_TL ~ Growth_ThreeGroups + 
              (1|GroupID) + (1|FishID) + (1|LC_Adj),
             family = gaussian(link = "identity"), 
             prior = c(set_prior('student_t(3, 0, 3)', class = 'Intercept'),
                       set_prior("normal(0, 4)", class = "b")),
             data = consistR1R2,
             cores = 4, 
             chains = 4,
             iter = 5000,
             warmup = 1000,
             control = list(adapt_delta = 0.999, max_treedepth = 20), 
             save_pars = save_pars(all = T), backend = 'cmdstanr'
)
### Model fit:
Model.1a. <- add_criterion(Model.1a., "loo",  save_psis = TRUE, reloo = T)
bayes_R2(Model.1a.) 
loo(Model.1a.)

### Model checks:
pp_check(Model.1a., ndraws = 50) 
pp_check(Model.1a.,  type = "scatter_avg", ndraws = 100) 
pp_check(Model.1a., type = "ecdf_overlay")
rstan::check_hmc_diagnostics(Model.1a.$fit)

mcmc_plot(Model.1a., type = "trace") 
mcmc_plot(Model.1a., type = "acf")   
mcmc_plot(Model.1a., type = 'rhat')  

### Model output:
plot(Model.1a.)
Model.1a. 
mcmc_plot(Model.1a., type = "dens") 

### EMmeans, contrasts and trends:
emmeans_results.Model.1a.2 <- emmeans(Model.1a., ~ Growth_ThreeGroups)
print(emmeans_results.Model.1a.2)
pairwise_comparisons.Model.1a.2 <- pairs(emmeans_results.Model.1a.2)
summary(pairwise_comparisons.Model.1a.2)

### Bayesian models of growth patterns- considering initial size and rank of the anemonefish and heat stress ####
### Growth pattern model 1b: Positive growth only:
### Subset:
consistR1R2.growth <- subset(consistR1R2, Shrink == 0)

Model.1b. <- brm(Perc_change_TL ~ Growth_ThreeGroups + 
                    (1|GroupID) + (1|FishID) + (1|LC_Adj),
                    family = gaussian(link = "identity"), 
                    prior = c(set_prior('student_t(3, 0, 3)', class = 'Intercept'),
                              set_prior("normal(0, 4)", class = "b")),
                    data = consistR1R2.growth,
                    cores = 4, 
                    chains = 4,
                    iter = 5000,
                    warmup = 1000,
                    control = list(adapt_delta = 0.999, max_treedepth = 20), 
                    save_pars = save_pars(all = T), backend = 'cmdstanr'
)

### Model fit:
Model.1b. <- add_criterion(Model.1b., "loo",  save_psis = TRUE, reloo = T)
bayes_R2(Model.1b.) 
loo(Model.1b.)

### Model checks:
pp_check(Model.1b., ndraws = 50) 
pp_check(Model.1b.,  type = "scatter_avg", ndraws = 100) 
pp_check(Model.1b., type = "ecdf_overlay")
rstan::check_hmc_diagnostics(Model.1b.$fit)

mcmc_plot(Model.1b., type = "trace") 
mcmc_plot(Model.1b., type = "acf")   
mcmc_plot(Model.1b., type = 'rhat')  

### Model output:
plot(Model.1b.)
Model.1b. 
mcmc_plot(Model.1b., type = "dens") 

### EMmeans, contrasts and trends:
emmeans_results.Model.1b.2 <- emmeans(Model.1b., ~ Growth_ThreeGroups)
print(emmeans_results.Model.1b.2)
pairwise_comparisons.Model.1b.2 <- pairs(emmeans_results.Model.1b.2)
summary(pairwise_comparisons.Model.1b.2)

# Model 1c: Shrinking only
### Subset:
consistR1R2.shrink <- subset(consistR1R2, Shrink == 1) %>% droplevels

Model.1c. <- brm(Perc_change_TL ~ Growth_ThreeGroups + 
                      (1|GroupID) + (1|FishID) + (1|LC_Adj),
                    family = gaussian(link = "identity"), 
                    prior = c(set_prior('student_t(3, 0, 3)', class = 'Intercept'),
                              set_prior("normal(0, 4)", class = "b")),
                    data = consistR1R2.shrink,
                    cores = 4, 
                    chains = 4,
                    iter = 5000,
                    warmup = 1000,
                    control = list(adapt_delta = 0.999, max_treedepth = 20), 
                    save_pars = save_pars(all = T), backend = 'cmdstanr'
)

### Model fit:
Model.1c. <- add_criterion(Model.1c., "loo",  save_psis = TRUE, reloo = T)
bayes_R2(Model.1c.) 
loo(Model.1c.)

### Model checks:
pp_check(Model.1c., ndraws = 50) 
pp_check(Model.1c.,  type = "scatter_avg", ndraws = 100) 
pp_check(Model.1c., type = "ecdf_overlay")
rstan::check_hmc_diagnostics(Model.1c.$fit)

mcmc_plot(Model.1c., type = "trace") 
mcmc_plot(Model.1c., type = "acf")   
mcmc_plot(Model.1c., type = 'rhat')  

### Model output:
plot(Model.1c.)
Model.1c. 
mcmc_plot(Model.1c., type = "dens") 

### EMmeans, contrasts and trends:
emmeans_results.Model.1c.2 <- emmeans(Model.1c., ~ Growth_ThreeGroups)
print(emmeans_results.Model.1c.2)
pairwise_comparisons.Model.1c.2 <- pairs(emmeans_results.Model.1c.2)
summary(pairwise_comparisons.Model.1c.2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Result 3: Clown anemonefish growth is predicted by heat stress ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Selecting heat stress predictors ####
### Model 2: Mean current temperature exposure:
Model.2. <- brm(Perc_change_TL ~ M_Temp_LC_scl + Rank 
                     + Init_TL_scl:Rank
                     + Init_TL_scl + (1|GroupID) + (1|LC_Adj),
                     family = gaussian(link = "identity"), 
                     prior = c(set_prior('student_t(3, 0, 3)', class = 'Intercept'),
                               set_prior("normal(0, 4)", class = "b")),
                     data = consistR1R2,
                     cores = 4, 
                     chains = 4,
                     iter = 5000,
                     warmup = 1000,
                     control = list(adapt_delta = 0.995, max_treedepth = 20), 
                     save_pars = save_pars(all = T), backend = 'cmdstanr'
)


### Model fit:
Model.2. <- add_criterion(Model.2., "loo",  save_psis = TRUE, reloo = T)
bayes_R2(Model.2.) 
loo(Model.2.)

### Model checks:
pp_check(Model.2., ndraws = 50) 
pp_check(Model.2.,  type = "scatter_avg", ndraws = 100) 
pp_check(Model.2., type = "ecdf_overlay")
rstan::check_hmc_diagnostics(Model.2.$fit)

mcmc_plot(Model.2., type = "trace") 
mcmc_plot(Model.2., type = "acf")   
mcmc_plot(Model.2., type = 'rhat')  

### Model output:
plot(Model.2.)
Model.2. 
mcmc_plot(Model.2., type = "dens") 

### Model 3: Mean preceding temperature exposure:
Model.3. <- brm(Perc_change_TL ~ M_Temp_LC_INIT_scl + Rank 
                      + Init_TL_scl:Rank
                      + Init_TL_scl + (1|GroupID) + (1|LC_Adj),
                      family = gaussian(link = "identity"), 
                      prior = c(set_prior('student_t(3, 0, 3)', class = 'Intercept'),
                                set_prior("normal(0, 4)", class = "b")),
                      data = consistR1R2,
                      cores = 4, 
                      chains = 4,
                      iter = 5000,
                      warmup = 1000,
                      control = list(adapt_delta = 0.995, max_treedepth = 20), 
                      save_pars = save_pars(all = T), backend = 'cmdstanr'
)


### Model fit:
Model.3. <- add_criterion(Model.3., "loo",  save_psis = TRUE, reloo = T)
bayes_R2(Model.3.) 
loo(Model.3.)

### Model checks:
pp_check(Model.3., ndraws = 50) 
pp_check(Model.3.,  type = "scatter_avg", ndraws = 100) 
pp_check(Model.3., type = "ecdf_overlay")
rstan::check_hmc_diagnostics(Model.3.$fit)

mcmc_plot(Model.3., type = "trace") 
mcmc_plot(Model.3., type = "acf")   
mcmc_plot(Model.3., type = 'rhat')  

### Model output:
plot(Model.3.)
Model.3. 
mcmc_plot(Model.3., type = "dens") 

### Model 4: Mean current & preceding temperature exposure:
Model.4. <- brm(Perc_change_TL ~ M_Temp_LC_INIT_scl + M_Temp_LC_scl
                      +  Rank 
                      + Init_TL_scl:Rank
                      + Init_TL_scl + (1|GroupID) + (1|LC_Adj),
                      family = gaussian(link = "identity"), 
                      prior = c(set_prior('student_t(3, 0, 3)', class = 'Intercept'),
                                set_prior("normal(0, 4)", class = "b")),
                      data = consistR1R2,
                      cores = 4, 
                      chains = 4,
                      iter = 5000,
                      warmup = 1000,
                      control = list(adapt_delta = 0.999, max_treedepth = 20), 
                      save_pars = save_pars(all = T), backend = 'cmdstanr'
)


### Model fit:
Model.4. <- add_criterion(Model.4., "loo",  save_psis = TRUE, reloo = T)
bayes_R2(Model.4.) 
loo(Model.4.)

### Model checks:
pp_check(Model.4., ndraws = 50) 
pp_check(Model.4.,  type = "scatter_avg", ndraws = 100) 
pp_check(Model.4., type = "ecdf_overlay")
rstan::check_hmc_diagnostics(Model.4.$fit)

mcmc_plot(Model.4., type = "trace") 
mcmc_plot(Model.4., type = "acf")   
mcmc_plot(Model.4., type = 'rhat')  

### Model output:
plot(Model.4.)
Model.4. 
mcmc_plot(Model.4., type = "dens") 

### Model 5: maximum current temperature exposure:
Model.5. <- brm(Perc_change_TL ~ Max_Temp_LC_scl + Rank 
                     + Init_TL_scl:Rank
                     + Init_TL_scl + (1|GroupID) + (1|LC_Adj),
                     family = gaussian(link = "identity"), 
                     prior = c(set_prior('student_t(3, 0, 3)', class = 'Intercept'),
                               set_prior("normal(0, 4)", class = "b")),
                     data = consistR1R2,
                     cores = 4, 
                     chains = 4,
                     iter = 5000,
                     warmup = 1000,
                     control = list(adapt_delta = 0.995, max_treedepth = 20), 
                     save_pars = save_pars(all = T), backend = 'cmdstanr'
)


### Model fit:
Model.5. <- add_criterion(Model.5., "loo",  save_psis = TRUE, reloo = T)
bayes_R2(Model.5.) 
loo(Model.5.)

### Model checks:
pp_check(Model.5., ndraws = 50) 
pp_check(Model.5.,  type = "scatter_avg", ndraws = 100) 
pp_check(Model.5., type = "ecdf_overlay")
rstan::check_hmc_diagnostics(Model.5.$fit)

mcmc_plot(Model.5., type = "trace") 
mcmc_plot(Model.5., type = "acf")   
mcmc_plot(Model.5., type = 'rhat')  

### Model output:
plot(Model.5.)
Model.5. 
mcmc_plot(Model.5., type = "dens") 

### Model 6: Maximum preceding temperature exposure:
Model.6. <- brm(Perc_change_TL ~ Max_Temp_LC_INIT_scl + Rank 
                      + Init_TL_scl:Rank
                      + Init_TL_scl + (1|GroupID) + (1|LC_Adj),
                      family = gaussian(link = "identity"), 
                      prior = c(set_prior('student_t(3, 0, 3)', class = 'Intercept'),
                                set_prior("normal(0, 4)", class = "b")),
                      data = consistR1R2,
                      cores = 4, 
                      chains = 4,
                      iter = 5000,
                      warmup = 1000,
                      control = list(adapt_delta = 0.99, max_treedepth = 20), 
                      save_pars = save_pars(all = T), backend = 'cmdstanr'
)


### Model fit:
Model.6. <- add_criterion(Model.6., "loo",  save_psis = TRUE, reloo = T)
bayes_R2(Model.6.) 
loo(Model.6.)

### Model checks:
pp_check(Model.6., ndraws = 50) 
pp_check(Model.6.,  type = "scatter_avg", ndraws = 100) 
pp_check(Model.6., type = "ecdf_overlay")
rstan::check_hmc_diagnostics(Model.6.$fit)

mcmc_plot(Model.6., type = "trace") 
mcmc_plot(Model.6., type = "acf")   
mcmc_plot(Model.6., type = 'rhat')  

### Model output:
plot(Model.6.)
Model.6. 
mcmc_plot(Model.6., type = "dens") 

### Model 7: Maximum current & preceding temperature exposure:
Model.7. <- brm(Perc_change_TL ~ Max_Temp_LC_INIT_scl + Max_Temp_LC_scl
                      +  Rank 
                      + Init_TL_scl:Rank
                      + Init_TL_scl + (1|GroupID) + (1|LC_Adj),
                      family = gaussian(link = "identity"), 
                      prior = c(set_prior('student_t(3, 0, 3)', class = 'Intercept'),
                                set_prior("normal(0, 4)", class = "b")),
                      data = consistR1R2,
                      cores = 4, 
                      chains = 4,
                      iter = 5000,
                      warmup = 1000,
                      control = list(adapt_delta = 0.99, max_treedepth = 20), 
                      save_pars = save_pars(all = T), backend = 'cmdstanr'
)


### Model fit:
Model.7. <- add_criterion(Model.7., "loo",  save_psis = TRUE, reloo = T)
bayes_R2(Model.7.) 
loo(Model.7.)

### Model checks:
pp_check(Model.7., ndraws = 50) 
pp_check(Model.7.,  type = "scatter_avg", ndraws = 100) 
pp_check(Model.7., type = "ecdf_overlay")
rstan::check_hmc_diagnostics(Model.7.$fit)

mcmc_plot(Model.7., type = "trace") 
mcmc_plot(Model.7., type = "acf")   
mcmc_plot(Model.7., type = 'rhat')  

### Model output:
plot(Model.7.)
Model.7. 
mcmc_plot(Model.7., type = "dens") 

### Effects estimates for this model extracted as this is the best model for heat stress. Decision based on loo_comparison
# results as displayed on line 874

### EMmeans, contrasts, and trends:
emmeans_results.Model.7. <- emmeans(Model.7., ~ Rank)
print(emmeans_results.Model.7.)
pairwise_comparisons.Model.7. <- pairs(emmeans_results.Model.7.)
summary(pairwise_comparisons.Model.7.)
trends_results.Model.7..1 <- emtrends(Model.7., pairwise ~ Rank, var = "Max_Temp_LC_INIT_scl")
summary(trends_results.Model.7..1)
trends_results.Model.7..2 <- emtrends(Model.7., pairwise ~ Rank, var = "Max_Temp_LC_scl")
summary(trends_results.Model.7..2)
trends_results_init.Model.7. <- emtrends(Model.7., pairwise ~ Rank, var = "Init_TL_scl")
summary(trends_results_init.Model.7.)

### Model 8: Range in current temperature exposure:
Model.8. <- brm(Perc_change_TL ~ Dev_Temp_LC_scl + Rank 
                     + Init_TL_scl:Rank
                     + Init_TL_scl + (1|GroupID) + (1|LC_Adj),
                     family = gaussian(link = "identity"), 
                     prior = c(set_prior('student_t(3, 0, 3)', class = 'Intercept'),
                               set_prior("normal(0, 4)", class = "b")),
                     data = consistR1R2,
                     cores = 4, 
                     chains = 4,
                     iter = 5000,
                     warmup = 1000,
                     control = list(adapt_delta = 0.995, max_treedepth = 20), 
                     save_pars = save_pars(all = T), backend = 'cmdstanr'
)


### Model fit:
Model.8. <- add_criterion(Model.8., "loo",  save_psis = TRUE, reloo = T)
bayes_R2(Model.8.) 
loo(Model.8.)

### Model checks:
pp_check(Model.8., ndraws = 50) 
pp_check(Model.8.,  type = "scatter_avg", ndraws = 100) 
pp_check(Model.8., type = "ecdf_overlay")
rstan::check_hmc_diagnostics(Model.8.$fit)

mcmc_plot(Model.8., type = "trace") 
mcmc_plot(Model.8., type = "acf")   
mcmc_plot(Model.8., type = 'rhat')  

### Model output:
plot(Model.8.)
Model.8. 
mcmc_plot(Model.8., type = "dens") 

### Model 9: Range in preceding temperature exposure:
Model.9. <- brm(Perc_change_TL ~ Dev_Temp_LC_INIT_scl + Rank 
                      + Init_TL_scl:Rank
                      + Init_TL_scl + (1|GroupID) + (1|LC_Adj),
                      family = gaussian(link = "identity"), 
                      prior = c(set_prior('student_t(3, 0, 3)', class = 'Intercept'),
                                set_prior("normal(0, 4)", class = "b")),
                      data = consistR1R2,
                      cores = 4, 
                      chains = 4,
                      iter = 5000,
                      warmup = 1000,
                      control = list(adapt_delta = 0.99, max_treedepth = 20), 
                      save_pars = save_pars(all = T), backend = 'cmdstanr'
)


### Model fit:
Model.9. <- add_criterion(Model.9., "loo",  save_psis = TRUE, reloo = T)
bayes_R2(Model.9.) 
loo(Model.9.)

### Model checks:
pp_check(Model.9., ndraws = 50) 
pp_check(Model.9.,  type = "scatter_avg", ndraws = 100) 
pp_check(Model.9., type = "ecdf_overlay")
rstan::check_hmc_diagnostics(Model.9.$fit)

mcmc_plot(Model.9., type = "trace") 
mcmc_plot(Model.9., type = "acf")   
mcmc_plot(Model.9., type = 'rhat')  

### Model output:
plot(Model.9.)
Model.9. 
mcmc_plot(Model.9., type = "dens") 

### Model 10: Range in current & preceding temperature exposure:
Model.10. <- brm(Perc_change_TL ~ Dev_Temp_LC_INIT_scl + Dev_Temp_LC_scl
                      +  Rank 
                      + Init_TL_scl:Rank
                      + Init_TL_scl + (1|GroupID) + (1|LC_Adj),
                      family = gaussian(link = "identity"), 
                      prior = c(set_prior('student_t(3, 0, 3)', class = 'Intercept'),
                                set_prior("normal(0, 4)", class = "b")),
                      data = consistR1R2,
                      cores = 4, 
                      chains = 4,
                      iter = 5000,
                      warmup = 1000,
                      control = list(adapt_delta = 0.999, max_treedepth = 20), 
                      save_pars = save_pars(all = T), backend = 'cmdstanr'
)


### Model fit:
Model.10. <- add_criterion(Model.10., "loo",  save_psis = TRUE, reloo = T)
bayes_R2(Model.10.) 
loo(Model.10.)

### Model checks:
pp_check(Model.10., ndraws = 50) 
pp_check(Model.10.,  type = "scatter_avg", ndraws = 100) 
pp_check(Model.10., type = "ecdf_overlay")
rstan::check_hmc_diagnostics(Model.10.$fit)

mcmc_plot(Model.10., type = "trace") 
mcmc_plot(Model.10., type = "acf")   
mcmc_plot(Model.10., type = 'rhat')  

### Model output:
plot(Model.10.)
Model.10. 
mcmc_plot(Model.10., type = "dens") 

### Model 11: Reduced model:
Model.11. <- brm(Perc_change_TL ~ Rank + Init_TL_scl +
                      (1|GroupID) + (1|LC_Adj),
                    family = gaussian(link = "identity"), 
                    prior = c(set_prior('student_t(3, 0, 3)', class = 'Intercept'),
                              set_prior("normal(0, 4)", class = "b")),
                    data = consistR1R2,
                    cores = 4, 
                    chains = 4,
                    iter = 5000,
                    warmup = 1000,
                    control = list(adapt_delta = 0.99, max_treedepth = 20), 
                    save_pars = save_pars(all = T), backend = 'cmdstanr'
)


### Model fit:
Model.11. <- add_criterion(Model.11., "loo",  save_psis = TRUE, reloo = T)
bayes_R2(Model.11.) 
loo(Model.11.)

### Model checks:
pp_check(Model.11., ndraws = 50) 
pp_check(Model.11.,  type = "scatter_avg", ndraws = 100) 
pp_check(Model.11., type = "ecdf_overlay")
rstan::check_hmc_diagnostics(Model.11.$fit)

mcmc_plot(Model.11., type = "trace") 
mcmc_plot(Model.11., type = "acf")   
mcmc_plot(Model.11., type = 'rhat')  

### Model output:
plot(Model.11.)
Model.11. 
mcmc_plot(Model.11., type = "dens") 

### Model 12: Null model:
Model.12. <- brm(Perc_change_TL ~ (1|GroupID) + (1|LC_Adj),
                  family = gaussian(link = "identity"), 
                  # prior = c(set_prior('student_t(3, 0, 3)', class = 'Intercept'),
                  #           set_prior("normal(0, 4)", class = "b")),
                  data = consistR1R2,
                  cores = 4, 
                  chains = 4,
                  iter = 5000,
                  warmup = 1000,
                  control = list(adapt_delta = 0.99, max_treedepth = 20), 
                  save_pars = save_pars(all = T), backend = 'cmdstanr'
)


### Model fit:
Model.12. <- add_criterion(Model.12., "loo",  save_psis = TRUE, reloo = T)
bayes_R2(Model.12.) 
loo(Model.12.)

### Model checks:
pp_check(Model.12., ndraws = 50) 
pp_check(Model.12.,  type = "scatter_avg", ndraws = 100) 
pp_check(Model.12., type = "ecdf_overlay")
rstan::check_hmc_diagnostics(Model.12.$fit)

mcmc_plot(Model.12., type = "trace") 
mcmc_plot(Model.12., type = "acf")   
mcmc_plot(Model.12., type = 'rhat')  

### Model output:
plot(Model.12.)
Model.12. 
mcmc_plot(Model.12., type = "dens") 

### Compare the models: LOO ####
loo_comparison.Model.Set.2 <- loo_compare(Model.2., Model.3., Model.4., Model.5., Model.6.,
                                          Model.7., Model.8., Model.9., Model.10., Model.12., Model.11.)
loo_comparison.Model.Set.2

### Combining heat stress models with social factors in anemonefish ####
### Model 13: Maximum current & preceding temperature exposure and rank interaction:
Model.13. <- brm(Perc_change_TL ~ Max_Temp_LC_INIT_scl + Max_Temp_LC_scl
                     + Rank 
                     + Max_Temp_LC_INIT_scl:Rank
                     + Max_Temp_LC_scl:Rank
                     + Init_TL_scl:Rank
                     + Init_TL_scl + (1|GroupID) + (1|LC_Adj),
                     family = gaussian(link = "identity"), 
                     prior = c(set_prior('student_t(3, 0, 3)', class = 'Intercept'),
                               set_prior("normal(0, 4)", class = "b")),
                     data = consistR1R2,
                     cores = 4, 
                     chains = 4,
                     iter = 5000,
                     warmup = 1000,
                     control = list(adapt_delta = 0.995, max_treedepth = 20), 
                     save_pars = save_pars(all = T), backend = 'cmdstanr'
)


### Model fit:
Model.13. <- add_criterion(Model.13., "loo",  save_psis = TRUE, reloo = T)
bayes_R2(Model.13.) 
loo(Model.13.)

### Model checks:
pp_check(Model.13., ndraws = 50) 
pp_check(Model.13.,  type = "scatter_avg", ndraws = 100) 
pp_check(Model.13., type = "ecdf_overlay")
rstan::check_hmc_diagnostics(Model.13.$fit)

mcmc_plot(Model.13., type = "trace") 
mcmc_plot(Model.13., type = "acf")   
mcmc_plot(Model.13., type = 'rhat')  

### Model output:
plot(Model.13.)
Model.13. 
mcmc_plot(Model.13., type = "dens") 

### EMmeans, contrasts, and trends:
emmeans_results.Model.13. <- emmeans(Model.13., ~ Rank)
print(emmeans_results.Model.13.)
pairwise_comparisons.Model.13. <- pairs(emmeans_results.Model.13.)
summary(pairwise_comparisons.Model.13.)
trends_results.Model.13. <- emtrends(Model.13., pairwise ~ Rank, var = "Max_Temp_LC_scl")
summary(trends_results.Model.13.)
trends_results_init.Model.13. <- emtrends(Model.13., pairwise ~ Rank, var = "Max_Temp_LC_INIT_scl")
summary(trends_results_init.Model.13.)
trends_results_size.Model.13. <- emtrends(Model.13., pairwise ~ Rank, var = "Init_TL_scl")
summary(trends_results_size.Model.13.)

### Compare the models: LOO ####
loo_comparison.Model.Set.3 <- loo_compare( Model.13., Model.12., Model.11.)
loo_comparison.Model.Set.3

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Result 4: Clown anemonefish shrinking is also predicted by social conflict ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Shrinking by the end of a lunar month - compared to size ratio at the start ####
### Model 14: Size ratio of rank 2: rank 1 at the start of the lunar month:
Model.14. <- brm(
  Shrinking ~ Size_R2R1_Initial_scl + Max_Temp_LC_INIT_scl + Rank +
    Size_R2R1_Initial_scl:Rank +
    Init_TL_scl + (1|GroupID) + (1|LC_Adj),
  family = bernoulli(link = "logit"),  
  prior = c(
    set_prior('student_t(3, 0, 3)', class = 'Intercept'),
    set_prior("normal(0, 4)", class = "b")
  ),
  data = binomR1R2.c,
  cores = 4,
  chains = 4,
  iter = 5000,
  warmup = 1000,
  control = list(adapt_delta = 0.995, max_treedepth = 20),
  # save_pars = save_pars(all = T), backend = 'cmdstanr'
)

### Model fit:
Model.14. <- add_criterion(Model.14., "loo",  save_psis = TRUE, reloo = T)
bayes_R2(Model.14.) 
loo(Model.14.)

### Model checks:
pp_check(Model.14., ndraws = 50) 
pp_check(Model.14.,  type = "scatter_avg", ndraws = 100) 
pp_check(Model.14., type = "ecdf_overlay")
rstan::check_hmc_diagnostics(Model.14.$fit)

mcmc_plot(Model.14., type = "trace") 
mcmc_plot(Model.14., type = "acf")   
mcmc_plot(Model.14., type = 'rhat')  

### Model output:
plot(Model.14.)
summary(Model.14.) 
mcmc_plot(Model.14., type = "dens") 

### Effects
ggpredict(Model.14., terms = c("Size_R2R1_Initial_scl", "Rank"))
model_parameters(Model.14., exponentiate = T)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Result 5: Clown anemonefish that shrink more often have higher chances of survival ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### For Survival models, factors need to be transformed to numeric:
### Set factors/ Numeric:
R1R2.complete.h$GroupID <- as.factor(R1R2.complete.h$GroupID)
R1R2.complete.h$Rank <- as.factor(R1R2.complete.h$Rank)
R1R2.complete.h$LC_Adj <- as.factor(R1R2.complete.h$LC_Adj)
R1R2.complete.h$Br_Pair_Shr <- as.factor(R1R2.complete.h$Br_Pair_Shr)
R1R2.complete.h$Growth_ThreeGroups <- as.factor(R1R2.complete.h$Growth_ThreeGroups)

R1R2.complete.h$TimeTo_Event <- as.numeric(R1R2.complete.h$TimeTo_Event)
R1R2.complete.h$Event_Died <- as.numeric(R1R2.complete.h$Event_Died)
R1R2.complete.h$Growth_ThreeGroups <- as.numeric(R1R2.complete.h$Growth_ThreeGroups)
R1R2.complete.h$Rank <- as.numeric(R1R2.complete.h$Rank)
R1R2.complete.h$Br_Pair_Shr <- as.numeric(R1R2.complete.h$Br_Pair_Shr)

### Calculation of temperature baseline at the level of an individual anemonefish:
R1R2.complete.h$MinOverall_FishID_scl <- scale(R1R2.complete.h$MinOverall_FishID, center=TRUE, scale=TRUE)

### Cox Hazard Model I- heat stress:
Model.15 <- coxph(Surv(TimeTo_Event, Event_Died) ~ MinOverall_FishID_scl , data = R1R2.complete.h)
summary(Model.15)
ph_test.Model.15 <- cox.zph(Model.15)
print(ph_test.Model.15)

### Cox Hazard Model II- initial size: 
Model.16 <- coxph(Surv(TimeTo_Event, Event_Died) ~ Init_TL_scl , data = R1R2.complete.h)
summary(Model.16)
ph_test.Model.16 <- cox.zph(Model.16)
print(ph_test.Model.16)

### Cox Hazard Model III- rank:
Model.17 <- coxph(Surv(TimeTo_Event, Event_Died) ~ Rank , data = R1R2.complete.h)
summary(Model.17)
ph_test.Model.17 <- cox.zph(Model.17)
print(ph_test.Model.17)

### Cox Hazard Model IV- shrinking within a breeding pair:
Model.18 <- coxph(Surv(TimeTo_Event, Event_Died) ~ Br_Pair_Shr, data = R1R2.complete.h)
summary(Model.18)
ph_test.Model.18 <- cox.zph(Model.18)
print(ph_test.Model.18)

### Cox Hazard Model V- growth patterns:
Model.19 <- coxph(Surv(TimeTo_Event, Event_Died) ~ Growth_ThreeGroups , data = R1R2.complete.h)
summary(Model.19)
ph_test.Model.19 <- cox.zph(Model.19)
print(ph_test.Model.19)

#~~~~~~~~~~~~~~~~~~~~~~~~
# Section 4: Figures ####
#~~~~~~~~~~~~~~~~~~~~~~~~

## Figure panel 1: Change in total length (TL) of clown anemonefish during a heat stress event #### 
## Fig. 1A: The distribution of percent body length changes of individual rank 1 and rank 2 clown anemonefish.
Fig1a <-  ggplot(R1R2.bmrs, aes(x=Perc_change_TL, fill=Rank, colour=Rank)) +
geom_vline(xintercept = 0, color = "black", linetype = "solid", size = .5) +
  geom_density( alpha = 0.7 ) +
  scale_fill_manual(values = c("1" = "#DF8F44", "2" = "#00A1D5" )) +
  scale_colour_manual(values = c( "1" = "#DF8F44", "2" = "#00A1D5")) +
  scale_x_continuous(breaks = seq(floor(-10), ceiling(15), by = 2.5)) +
  annotate("text", x = -6.5, y = .3, label = "Shrinking", color = "black", vjust = -0.5, hjust = 0, family = "sans", size = 5) +
  annotate("text", x = 3, y = .3, label = "Growing", color = "black", vjust = -0.5, hjust = 0, family = "sans", size = 5) +
  labs(
    x = "Percent change in Total Length",
    y = "Density",
    fill = "Rank"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "sans", size = 16),
    legend.position = "none",
    axis.title = element_text(family = "sans", size = 16),
    axis.text = element_text(family = "sans", size = 16),
    plot.margin = margin(t = 0, r = 0, b = 10, l = 10, unit = "pt"),
    legend.title = element_text(family = "sans", size = 16),
    legend.text = element_text(family = "sans", size = 16)
  )

## Fig. 1B: The percent change in total length modelled in relation to initial total length of fish for rank 1.
### Rank 1 anemonefish:
R1 <- subset(consistR1R2, Rank == 1)
### Zoom
x_limits.Fig1R1 <- c(40, 85)  
y_limits.Fig1R1 <- c(-3, 4)  

Fig1b <- ggplot(R1, aes(x = Init_TL, y = Perc_change_TL)) + 
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 1, alpha = 0.3) +
  geom_vline(xintercept = 79.92, color = "black", linetype = "dashed", linewidth = 1) +
  geom_point(aes(shape = ifelse(Perc_change_TL < 0, "Shrinking", "Growing")), 
             color = "#C7B3E5", size = 5) +
  geom_smooth(method = "loess", color = NA, fill = "#C7B3E5", se = TRUE, alpha = 0.1) +
  annotate("text", x = 79.92, y = -2.5, label = "Threshold", color = "black", angle = 90, vjust = -0.5, hjust = 0.5, family = "sans", size = 5) +
  annotate("text", x = 41, y = 0.5, label = "Growing", color = "black", hjust = 0.0, family = "sans", size = 5) +
  annotate("text", x = 41, y = -0.5, label = "Shrinking", color = "black", hjust = 0.0, family = "sans", size = 5) +
  ggtitle("Rank 1") + 
  labs(y = "Change in total length (%)",
       x = "Initial total length (mm)") +
  scale_shape_manual(values = c("Shrinking" = 16, "Growing" = 1)) +
  theme_classic() +
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(size = 16),
    legend.position = "none",   
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14)
  ) +
  coord_cartesian(xlim = x_limits.Fig1R1, ylim = y_limits.Fig1R1)

## Fig. 1C: The percent change in total length modelled in relation to initial total length of fish for rank 2.     
### Rank 2 anemonefish:
R2 <- subset(consistR1R2, Rank == 2) 
### Zoom
x_limits.Fig1R2 <- c(40, 85)  
y_limits.Fig1R2 <- c(-3, 4)  

Fig1c <- ggplot(R2, aes(x = Init_TL, y = Perc_change_TL)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 1, alpha = 0.3) +
  geom_vline(xintercept = 60.75, color = "black", linetype = "dashed", linewidth = 1) +
  geom_point(aes(shape = ifelse(Perc_change_TL < 0, "Shrinking", "Growing")), 
             color = "#B4E197", size = 5) +
  geom_smooth(method = "loess", color = NA, fill = "#B4E197", se = TRUE, alpha = 0.1) +
  annotate("text", x = 60.75, y = -2.5, label = "Threshold", color = "black", angle = 90, vjust = -0.5, hjust = 0.5, family = "sans", size = 5) +
  annotate("text", x = 41, y = 0.5, label = "Growing", color = "black", hjust = 0.0 , family = "sans", size = 5) +
  annotate("text", x = 41, y = -0.5, label = "Shrinking", color = "black", hjust = 0.0 , family = "sans", size = 5) +
  ggtitle("Rank 2") + 
  labs( y = "Change in total length (%)",
        x = "Initial total length (mm)"
  ) +
  scale_shape_manual(values = c("Shrinking" = 16, "Growing" = 1)) +
  theme_classic() +
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(size = 16),
    legend.position = "none",   
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14)
  ) +
  coord_cartesian(xlim = x_limits.Fig1R2, ylim = y_limits.Fig1R2)

## Figure panel 2: Percent change in total length (TL) for clown anemonefish with different shrinking response patterns ####
### Top row: The slope estimates of rank 1 and rank 2 clown anemonefish using each of different types of growth patterns, for all observations.
### All observations of Percent change in TL:
emmeans_results.Model.1a.2 <- emmeans(Model.1a., ~ Growth_ThreeGroups)
print(emmeans_results.Model.1a.2)
pairwise_comparisons.Model.1a.2 <- pairs(emmeans_results.Model.1a.2)
summary(pairwise_comparisons.Model.1a.2)

### Estimates
Fig2.1 <- consistR1R2 %>%
  modelr::data_grid(Growth_ThreeGroups) %>%
  add_epred_draws(Model.1a., re_formula = NA) %>%
  mutate(Mpred=mean(.epred)) %>%
  sample_draws(10000) %>%
  ggplot(aes(x = .epred, y = Growth_ThreeGroups, fill = Growth_ThreeGroups)) +
  # stat_pointinterval(.width = c(.89, .95))+
  stat_halfeye(.width = c(.89, .95)) +
  scale_fill_viridis_d() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic()+
  coord_cartesian(xlim=c(-3,4))+
  guides(fill="none")

### Middle row: The slope estimates of rank 1 and rank 2 clown anemonefish using each of different types of growth patterns, for postive growth only.
### Growth only:
emmeans_results.Model.1b.2 <- emmeans(Model.1b., ~ Growth_ThreeGroups)
print(emmeans_results.Model.1b.2)
pairwise_comparisons.Model.1b.2 <- pairs(emmeans_results.Model.1b.2)
summary(pairwise_comparisons.Model.1b.2)

### Estimates
Fig2.2 <- consistR1R2.growth %>%
  modelr::data_grid(Growth_ThreeGroups) %>%
  add_epred_draws(Model.1b., re_formula = NA) %>%
  mutate(Mpred=mean(.epred)) %>%
  sample_draws(10000) %>%
  ggplot(aes(x = .epred, y = Growth_ThreeGroups, fill = Growth_ThreeGroups)) +
  # stat_pointinterval(.width = c(.89, .95))+
  stat_halfeye(.width = c(.89, .95)) +
  scale_fill_viridis_d() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic()+
  coord_cartesian(xlim=c(-3,4))+
  guides(fill="none")

### Bottom row: The slope estimates of rank 1 and rank 2 clown anemonefish using each of different types of growth patterns, for shrinking only.
### Shrinking only:
emmeans_results.Model.1c.2 <- emmeans(Model.1c., ~ Growth_ThreeGroups)
print(emmeans_results.Model.1c.2)
pairwise_comparisons.Model.1c.2 <- pairs(emmeans_results.Model.1c.2)
summary(pairwise_comparisons.Model.1c.2)

### Estimates
Fig2.3 <- consistR1R2.shrink %>%
  modelr::data_grid(Growth_ThreeGroups) %>%
  add_epred_draws(Model.1c., re_formula = NA) %>%
  mutate(Mpred=mean(.epred)) %>%
  sample_draws(10000) %>%
  ggplot(aes(x = .epred, y = Growth_ThreeGroups, fill = Growth_ThreeGroups)) +
  # stat_pointinterval(.width = c(.89, .95))+
  stat_halfeye(.width = c(.89, .95)) +
  scale_fill_viridis_d() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic()+
  coord_cartesian(xlim=c(-3,4))+
  guides(fill="none")

### Combine 3 plots:
Fig2 <- Fig2.1 / Fig2.2 / Fig2.3
print(Fig2)

## Figure panel 3: Change in body length of clown anemonefish in relation to heat stress ####
## Fig. 3A: The percentage change in total length fitted against temperature exposure in the lunar month of length measurements for rank 1 and rank 2. 
Fig3a <- consistR1R2 %>%
  group_by(Rank) %>%
  modelr::data_grid(Max_Temp_LC_scl = seq_range(Max_Temp_LC_scl, n = 501),
                    Max_Temp_LC_INIT_scl = mean(Max_Temp_LC_INIT_scl),
                    Init_TL_scl = mean(Init_TL_scl)) %>%
  add_epred_draws(Model.13., re_formula = NA) %>%
  mutate(Mpred=mean(.epred)) %>%
  sample_draws(100) %>%
  ggplot(aes(x = Max_Temp_LC_scl, y = Perc_change_TL, color = ordered(Rank))) +
  geom_line(aes(y = .epred, group = paste(Rank, .draw)), alpha = 0.1) +
  geom_point(data = consistR1R2, alpha = 0.2)+
  geom_line(aes(y = Mpred, group = paste(Rank, .draw)), alpha = 1, size = 1) +
  scale_colour_manual(values = c( "1" = "#DF8F44", "2" = "#00A1D5")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()+
  theme(
    text = element_text(family = "sans", size = 12),
    legend.position = "none",
    axis.title = element_text(family = "sans", size = 12),
    axis.text = element_text(family = "sans", size = 12),
    plot.margin = margin(t = 0, r = 0, b = 10, l = 10, unit = "pt"),
    legend.title = element_text(family = "sans", size = 12),
    legend.text = element_text(family = "sans", size = 12)
  )

## Fig. 3B: The distribution of the corresponding slope estimates for each rank, for current heat stress model. 
### Trends ####
trends_results.Model.13..Plot1 <- emtrends(Model.13., ~ Rank, var = "Max_Temp_LC_scl")
summary(trends_results.Model.13..Plot1)
pd.currentTemp <- p_direction(trends_results.Model.13..Plot1)
pd.currentTemp.Plot <- plot(pd.currentTemp)
pd.currentTemp.Plot

Fig3b <- pd.currentTemp.Plot +
  annotate("text", x = 0.4, y = 2.1, label = "Rank 1 change in size", color = "black", family = "sans", size = 5) +
  annotate("text", x = 0.7, y = 1.1, label = "Rank 2 change in size", color = "black", family = "sans", size = 5) +
  labs(
    title = "Current heat stress exposure parameter estimates and contrasts between ranks",
    x = "Parameter values",
    y = "Rank",
    fill = "Value"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14), 
    legend.position = "none", 
    axis.title.x = element_blank(),
    axis.title.y = element_text(family = "sans", size = 12),
    axis.text = element_text(family = "sans"),
    legend.title = element_text(family = "sans", size = 12),
    legend.text = element_text(family = "sans", size = 10)
  )

## Fig. 3C: The distribution of the corresponding slope contrasts between ranks for current heat stress model.
### Contrast ####
summary(trends_results.Model.13..Plot1)
contrast_results_Plot1 <- contrast(trends_results.Model.13..Plot1, method = "pairwise")
summary(contrast_results_Plot1)
# Plot the contrasts:
pd.currentTempContr <- p_direction(contrast_results_Plot1)
pd.currentTempContr.Plot <- plot(pd.currentTempContr)
pd.currentTempContr.Plot

Fig3c <- pd.currentTempContr.Plot +
  annotate("text", x = -0.4, y = 0.05, label = "Contrast rank 1 - rank 2", color = "black", family = "sans", size = 5) +
  labs(
    title = "Current heat stress exposure contrast between rank 1 and rank 2",
    x = "Parameter values",
    y = "Posterior distribution",
    fill = "Value"
  ) +
  theme_classic() +
  theme(
    plot.title = element_blank(), 
    legend.position = "none", 
    axis.title.x = element_text(family = "sans", size = 12),
    axis.title.y = element_text(family = "sans", size = 12),
    axis.text = element_text(family = "sans"),
    legend.title = element_text(family = "sans", size = 12),
    legend.text = element_text(family = "sans", size = 10)
  )

## Fig. 3D: The percentage change in total length was fitted against temperature exposure in the preceding lunar month. 
Fig3d <- consistR1R2 %>%
  group_by(Rank) %>%
  modelr::data_grid(Max_Temp_LC_scl = mean(Max_Temp_LC_scl),
                    Max_Temp_LC_INIT_scl = seq_range(Max_Temp_LC_INIT_scl, n = 501),
                    Init_TL_scl = mean(Init_TL_scl)) %>%
  add_epred_draws(Model.13., re_formula = NA) %>%
  mutate(Mpred=mean(.epred)) %>%
  sample_draws(100) %>%
  ggplot(aes(x = Max_Temp_LC_INIT_scl, y = Perc_change_TL, color = ordered(Rank))) +
  geom_line(aes(y = .epred, group = paste(Rank, .draw)), alpha = 0.1) +
  geom_point(data = consistR1R2, alpha = 0.2)+
  geom_line(aes(y = Mpred, group = paste(Rank, .draw)), alpha = 1, size = 1) +
  scale_colour_manual(values = c( "1" = "#DF8F44", "2" = "#00A1D5")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()+
  theme(
    text = element_text(family = "sans", size = 12),
    legend.position = "none",
    axis.title = element_text(family = "sans", size = 12),
    axis.text = element_text(family = "sans", size = 12),
    plot.margin = margin(t = 0, r = 0, b = 10, l = 10, unit = "pt"),
    legend.title = element_text(family = "sans", size = 12),
    legend.text = element_text(family = "sans", size = 12)
  )

## Fig. 3E: The distribution of the corresponding slope estimates for each rank, for preceding heat stress model. 
### Trends ####
trends_results.Model.13..Plot2 <- emtrends(Model.13., ~ Rank, var = "Max_Temp_LC_INIT_scl")
summary(trends_results.Model.13..Plot2)
pd.precedingTemp <- p_direction(trends_results.Model.13..Plot2)
plot(pd.precedingTemp)
pd.precedingTemp
pd.precedingTemp.Plot <- plot(pd.precedingTemp)

Fig3e <- pd.precedingTemp.Plot +
  annotate("text", x = -0.3, y = 2.1, label = "Rank 1 change in size", color = "black", family = "sans", size = 5) +
  annotate("text", x = -0.5, y = 1.1, label = "Rank 2 change in size", color = "black", family = "sans", size = 5) +
  labs(
    title = "Preceding heat stress exposure parameter estimates and contrasts between ranks",
    x = "Parameter values",
    y = "Rank",
    fill = "Value"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5), 
    legend.position = "none", 
    axis.title.x = element_blank(),
    axis.title.y = element_text(family = "sans", size = 12),
    axis.text = element_text(family = "sans"),
    legend.title = element_text(family = "sans", size = 12),
    legend.text = element_text(family = "sans", size = 10)
  )

## Fig. 3F: The distribution of the corresponding slope contrasts between ranks for preceding heat stress model.
### Contrasts ####
summary(trends_results.Model.13..Plot2)
contrast_results_Plot2 <- contrast(trends_results.Model.13..Plot2, method = "pairwise")
summary(contrast_results_Plot2)
pd.precedingTempContr <- p_direction(contrast_results_Plot2)
plot(pd.precedingTempContr)
pd.precedingTempContr
pd.precedingTempContr.Plot <- plot(pd.precedingTempContr)

Fig3f <- pd.precedingTempContr.Plot +
  annotate("text", x = 0.2, y = 0.05, label = "Contrast rank 1 - rank 2", color = "black", family = "sans", size = 5) +
  labs(
    title = "Preceding heat stress exposure contrast between rank 1 and rank 2",
    x = "Parameter values",
    y = "Posterior distribution",
    fill = "Value"
  ) +
  theme_classic() +
  theme(
    plot.title = element_blank(), 
    legend.position = "none", 
    axis.title.x = element_text(family = "sans", size = 12),
    axis.title.y = element_text(family = "sans", size = 12),
    axis.text = element_text(family = "sans"),
    legend.title = element_text(family = "sans", size = 12),
    legend.text = element_text(family = "sans", size = 10)
  )

## Figure Panel 4: Shrinking likelihood in relation to size ratio of both clown anemonefish within a breeding pair ####
## Fig. 4A: The likelihood of shrinking was fitted against size ratios of rank 2 to rank 1 at the start of the lunar month.
(scale(binomR1R2.c$Size_R2R1_Initial, center = T, scale = T)) %>% str 
back_scale <- function(x){x * 0.0525 + 0.778} 
Fig4a <- binomR1R2.c %>%
group_by(Rank) %>%
modelr::data_grid(Size_R2R1_Initial_scl = seq_range(Size_R2R1_Initial_scl, n = 501),
        Max_Temp_LC_INIT_scl = 0,
        Init_TL_scl = 0) %>%
add_epred_draws(Model.14., re_formula = NA) %>%
sample_draws(1000) %>%
mutate(Mpred=mean(.epred),
Size_R2R1_Initial_bs = back_scale(Size_R2R1_Initial_scl) ) %>%
ggplot(aes(x = Size_R2R1_Initial_bs, y = .epred, color = ordered(Rank))) +
geom_line(aes(y = .epred, group = paste(Rank, .draw)), alpha = 0.05) +
geom_line(aes(y = Mpred, group = paste(Rank, .draw)), alpha = 1, size = 1) +
scale_colour_manual(values = c( "1" = "#DF8F44", "2" = "#00A1D5")) +
theme_classic()+
  labs(x = "Size ratio", y = "Probability of shrinking") +
theme(
text = element_text(family = "sans", size = 12),
  legend.position = "none",
  axis.title = element_text(family = "sans", size = 12),
  axis.text = element_text(family = "sans", size = 12),
  plot.margin = margin(t = 0, r = 0, b = 10, l = 10, unit = "pt"),
  legend.title = element_text(family = "sans", size = 12),
  legend.text = element_text(family = "sans", size = 12)
  )

## Fig. 4B: The distribution of model estimates for shrinking probability among rank 1 and rank 2 individuals, at three discrete size ratios 0.6, 0.8, and 1.0.
Fig4b <- binomR1R2.c %>%
group_by(Rank) %>%
modelr::data_grid(Size_R2R1_Initial_scl = (c(.6, .8, 1)- 0.778)/0.0525 ,
Max_Temp_LC_INIT_scl = 0,
Init_TL_scl = 0) %>%
add_epred_draws(Model.14., re_formula = NA) %>%
sample_draws(1000) %>%
mutate(Mpred=mean(.epred),
Size_R2R1_Initial_bs = back_scale(Size_R2R1_Initial_scl) ) %>%
ggplot(aes(x = .epred, y = Rank, fill = Rank)) +
# stat_pointinterval(.width = c(.89, .95))+
stat_halfeye(.width = c(.89, .95)) +
  scale_fill_manual(values = c( "1" = "#DF8F44", "2" = "#00A1D5")) +
geom_vline(xintercept = 0, linetype = "dashed") +
theme_classic()+
facet_wrap(~Size_R2R1_Initial_bs)+
# coord_cartesian(xlim=c(-3,4))+
# guides(fill="none")
labs(x = "Probability of shrinking", y = "Density") +
  theme(
  text = element_text(family = "sans", size = 12),
  legend.position = "none",
  axis.title = element_text(family = "sans", size = 12),
  axis.text = element_text(family = "sans", size = 12),
  plot.margin = margin(t = 0, r = 0, b = 10, l = 10, unit = "pt"),
  legend.title = element_text(family = "sans", size = 12),
  legend.text = element_text(family = "sans", size = 12)
  )

## Fig. 4C: The corresponding contrasts in shrinking probabilities between the two ranks.
Fig4c <- binomR1R2.c %>%
group_by(Rank) %>%
modelr::data_grid(Size_R2R1_Initial_scl = (c(.6, .8, 1)- 0.778)/0.0525 ,
Max_Temp_LC_INIT_scl = 0,
Init_TL_scl = 0) %>%
add_epred_draws(Model.14., re_formula = NA) %>%
sample_draws(1000) %>%
mutate(Size_R2R1_Initial_bs = back_scale(Size_R2R1_Initial_scl)) %>%
ungroup() %>%
select(.draw, .epred, Rank, Size_R2R1_Initial_bs) %>%
group_by(.draw, Rank) %>%
pivot_wider(names_from = Rank, names_prefix = "Rank_", values_from = .epred) %>%
mutate(Contrast = Rank_2 - Rank_1) %>%
ggplot(aes(x = Contrast)) +
# stat_pointinterval(.width = c(.89, .95))+
stat_halfeye(.width = c(.89, .95)) +
scale_fill_viridis_d() +
geom_vline(xintercept = 0, linetype = "dashed") +
theme_classic()+
facet_wrap(~Size_R2R1_Initial_bs)+
# coord_cartesian(xlim=c(-3,4))+
# guides(fill="none")
labs(x = "Between rank contrast of shrinking probability", y = "Density") +
  theme(
  text = element_text(family = "sans", size = 12),
  legend.position = "none",
  axis.title = element_text(family = "sans", size = 12),
  axis.text = element_text(family = "sans", size = 12),
  plot.margin = margin(t = 0, r = 0, b = 10, l = 10, unit = "pt"),
  legend.title = element_text(family = "sans", size = 12),
  legend.text = element_text(family = "sans", size = 12)
  )

### Combine 2 plots:
Fig4 <- Fig4a|(Fig4b / Fig4c)
print(Fig4)

## Figure panel 5: Survival of clown anemonefish ####
### Scale variables:
R1R2.complete.h$Init_TL_scl <- scale(R1R2.complete.h$Init_TL, center=TRUE, scale=TRUE)
R1R2.complete.h$MinOverall_FishID_scl <- scale(R1R2.complete.h$MinOverall_FishID, center=TRUE, scale=TRUE)
### Numeric (required for model package)
R1R2.complete.h$TimeTo_Event <- as.numeric(R1R2.complete.h$TimeTo_Event)
R1R2.complete.h$Event_Died <- as.numeric(R1R2.complete.h$Event_Died)
R1R2.complete.h$Growth_ThreeGroups <- as.numeric(R1R2.complete.h$Growth_ThreeGroups)
R1R2.complete.h$Rank <- as.numeric(R1R2.complete.h$Rank)
R1R2.complete.h$Br_Pair_Shr <- as.numeric(R1R2.complete.h$Br_Pair_Shr)

## Fig. 5A: Kaplan-Meier survival curves for heat stress at the anemone scale. 
Model.15 <- coxph(Surv(TimeTo_Event, Event_Died) ~ MinOverall_FishID_scl , data = R1R2.complete.h)
summary(Model.15)

### Temperature into quartiles
R1R2.complete.h$MinOverall_FishID_cat <- cut(R1R2.complete.h$MinOverall_FishID, 
                                             breaks = quantile(R1R2.complete.h$MinOverall_FishID, probs = seq(0, 1, 0.25)),
                                             include.lowest = TRUE, 
                                             labels = c("Low heat stress", "Below Average", "Above Average", "High heat stress"))

### Change colours
custom_temps <- c(
  "Low minimum" = "#00008B",          
  "Below average minimum" = "yellow", 
  "Above average minimum" = "orange", 
  "High minimum" = "#FF0000"           
)
### Change names
names(custom_temps) <- c("Low heat stress", "Below average", "Above average", "High heat stress")

### Plot KM curves ####
surv_fit.temp <- survfit(Surv(TimeTo_Event, Event_Died) ~ MinOverall_FishID_cat, data = R1R2.complete.h)

Fig5.a <- ggsurvplot(
  surv_fit.temp,                       
  data = R1R2.complete.h,         
  conf.int = TRUE,                
  legend.labs = names(custom_temps), 
  xlab = "Time (days)",         
  ylab = "Survival probability",  
  palette = custom_temps        
)
### Zoom
x_limits.SurvTemps <- c(75, 180)  
y_limits.SurvTemps <- c(0.4, 1)  
### Customise 
Fig5.a$plot <- Fig5.a$plot + 
  theme_classic(base_family = "sans") +  
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "none",
  ) +
  annotate("text", x = 200, y = 0.8, label = "Low temperatures", color = custom_temps["Low heat stress"], hjust = 0) +
  annotate("text", x = 200, y = 0.6, label = "Below average", color = custom_temps["Below average"], hjust = 0) +
  annotate("text", x = 200, y = 0.4, label = "Above average", color = custom_temps["Above average"], hjust = 0) +
  annotate("text", x = 200, y = 0.2, label = "High temperatures", color = custom_temps["High heat stress"], hjust = 0)  +
  coord_cartesian(xlim = x_limits.SurvTemps, ylim = y_limits.SurvTemps)

## Fig. 5B: Kaplan-Meier survival curves for growth patterns based on shrinking frequency. 
Model.19 <- coxph(Surv(TimeTo_Event, Event_Died) ~ Growth_ThreeGroups , data = R1R2.complete.h)
summary(Model.19)
### Change colours
custom_growth <- c(
  "1" = "grey",          
  "2" = "pink",
  "3" = "red"
)
### Change names
names(custom_growth) <- c("No shrinking", "Single shrinking", "Multiple shrinking")
### Plot KM curves ####
surv_fit_growth <- survfit(Surv(TimeTo_Event, Event_Died) ~ Growth_ThreeGroups, data = R1R2.complete.h)

Fig5.b <-ggsurvplot(
  surv_fit_growth,                       
  data = R1R2.complete.h,         
  conf.int = TRUE,                
  legend.labs = c("No shrinking", "Single shrinking", "Multiple shrinking"), 
  xlab = "Time (days)",         
  ylab = "Survival probability",
  palette = custom_growth       # set custom colors
)
### Zoom
x_limits.SurvGrowPatt <- c(75, 180)  # Adjust these values as needed
y_limits.SurvGrowPatt <- c(0.6, 1)  # Adjust these values as needed
### Customise 
Fig5.b$plot <- Fig5.b$plot + 
  theme_classic() +  
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "none",
  ) +
  annotate("text", x = 200, y = 0.8, label = "No shrinking", color = custom_growth["1"], hjust = 0) +
  annotate("text", x = 200, y = 0.6, label = "Single shrinking", color = custom_growth["2"], hjust = 0) +
  annotate("text", x = 200, y = 0.6, label = "Multiple shrinking", color = custom_growth["3"], hjust = 0)  +
  coord_cartesian(xlim = x_limits.SurvGrowPatt, ylim = y_limits.SurvGrowPatt)

## Fig. 5C: Kaplan-Meier survival curves for social factors of shrinking within a breeding pair.
Model.18 <- coxph(Surv(TimeTo_Event, Event_Died) ~ Br_Pair_Shr, data = R1R2.complete.h)
summary(Model.18)
### Change colours
custom_Br_Pair <- c(
  "1" = "purple",          
  "2" = "green"
)
### Change names
names(custom_Br_Pair) <- c("Not shrinking within breeding pairs", "Shrinking within breeding pairs")
### Plot KM curves ####
surv_fit_Br_Pair <- survfit(Surv(TimeTo_Event, Event_Died) ~ Br_Pair_Shr, data = R1R2.complete.h)

Fig5.c <-ggsurvplot(
  surv_fit_Br_Pair,                       
  data = R1R2.complete.h,         
  conf.int = TRUE,                
  legend.labs = c("Not shrinking within breeding pairs", "Shrinking within breeding pairs"), 
  xlab = "Time (days)",         
  ylab = "Survival probability",
  palette = custom_Br_Pair       # set custom colors
)
### Zoom
x_limits.SurvBrPair <- c(75, 180)  # Adjust these values as needed
y_limits.SurvBrPair <- c(0.7, 1)  # Adjust these values as needed
### Customise 
Fig5.c$plot <- Fig5.c$plot + 
  theme_classic() +  
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "none",
  ) +
  annotate("text", x = 200, y = 0.8, label = "No shrinking within breeding pair", color = custom_Br_Pair["1"], hjust = 0) +
  annotate("text", x = 200, y = 0.6, label = "Shrinking within breeding pair", color = custom_Br_Pair["2"], hjust = 0)  +
  coord_cartesian(xlim = x_limits.SurvBrPair, ylim = y_limits.SurvBrPair)

### Combine 3 plots:
Fig5 <- (Fig5.a$plot/Fig5.b$plot/Fig5.c$plot)
print(Fig5)