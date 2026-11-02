
# Lung Cancer Survival Analysis

install.packages("survival")

library(survival)

# Load the NCCTG lung cancer dataset
data("lung", package = "survival")

# Dataset dimensions
dim(lung)

# First six observations
head(lung)

# Variable names and data types
str(lung)

# Summary statistics
summary(lung)

# Missing values by variable
colSums(is.na(lung))

# Survival event status
# 1 = censored, 2 = death
table(lung$status)

# Open the dataset viewer
View(lung)


# --------------------------------------------------
# Step 2: Data preparation
# --------------------------------------------------

# Create a working copy of the original dataset
lung_data <- lung

# Convert survival status into a binary event indicator
# Original: 1 = censored, 2 = death
# New:      0 = censored, 1 = death
lung_data$event <- as.integer(lung_data$status == 2)

# Convert sex into a labelled factor
# Original: 1 = male, 2 = female
lung_data$sex <- factor(
  lung_data$sex,
  levels = c(1, 2),
  labels = c("Male", "Female")
)

# Check event and sex distributions
table(lung_data$event)
table(lung_data$sex)

# Check missing values
colSums(is.na(lung_data))

# Summary of key variables
summary(lung_data[, c("time", "event", "age", "sex", "ph.ecog")])

# Verify that event coding is correct
table(lung_data$status, lung_data$event)




# Step 3: Kaplan-Meier survival analysis



# Create the survival object
surv_object <- Surv(
  time = lung_data$time,
  event = lung_data$event
)

# Fit overall Kaplan-Meier survival curve
km_overall <- survfit(surv_object ~ 1)


# Overall Kaplan-Meier summary
print(km_overall)

# Median survival and 95% confidence interval
summary(km_overall)$table

# Estimated survival at selected time points
summary(km_overall, times = c(180, 365))


# Plot overall survival curve
plot(
  km_overall,
  conf.int = TRUE,
  col = "blue",
  lwd = 2,
  xlab = "Time (days)",
  ylab = "Estimated survival probability",
  main = "Kaplan-Meier Overall Survival"
)

# The dataset contains 228 patients, with 165 observed
# deaths (72.4%) and 63 censored observations (27.6%).
#
# The estimated median survival time was 310 days,
# with a 95% confidence interval of 285 to 363 days.
#
# At 180 days, the estimated survival probability was
# 72.2% (95% CI: 66.6% to 78.3%).
#
# At 365 days, the estimated survival probability was
# 40.9% (95% CI: 34.5% to 48.6%).
#
# The Kaplan-Meier curve shows declining survival
# probability over time. These findings describe overall
# survival and do not establish which patient characteristics
# are associated with survival.



# Step 4: Kaplan-Meier Survival Analysis by Sex


# Objective:
# Investigate whether survival differs between
# male and female lung cancer patients.

# Fit Kaplan-Meier survival curves by sex
km_sex <- survfit(
  Surv(time, event) ~ sex,
  data = lung_data
)

# Display median survival and confidence intervals
print(km_sex)

# Survival probabilities at 180 and 365 days
summary(km_sex, times = c(180, 365))

# Plot survival curves by sex
plot(
  km_sex,
  col = c("blue", "red"),
  lwd = 2,
  conf.int = FALSE,
  mark.time = TRUE,
  xlab = "Time (days)",
  ylab = "Estimated survival probability",
  main = "Kaplan-Meier Survival by Sex"
)

legend(
  "topright",
  legend = c("Male", "Female"),
  col = c("blue", "red"),
  lwd = 2,
  bty = "n"
)

# Perform log-rank test
logrank_sex <- survdiff(
  Surv(time, event) ~ sex,
  data = lung_data
)

print(logrank_sex)

# Calculate p-value
p_value_sex <- pchisq(
  logrank_sex$chisq,
  df = length(logrank_sex$n) - 1,
  lower.tail = FALSE
)

print(p_value_sex)



# Kaplan-Meier survival curves were compared between
# male and female patients.
#
# The female survival curve generally remained above
# the male survival curve, suggesting better survival
# outcomes among female patients in this cohort.
#
# A log-rank test was performed to compare the
# survival functions between the two groups.
#
# The test produced a p-value of approximately 0.00131.
# Since p < 0.05, we reject the null hypothesis of
# equal survival functions.
#
# These results provide statistically significant
# evidence of a survival difference between male
# and female patients.
#
# However, this is an unadjusted comparison and does
# not account for other potential explanatory variables
# such as age or clinical performance status.



# Step 5: Cox Proportional Hazards Regression


# Objective:
# Investigate associations between age, sex, ECOG
# performance status, and overall survival.

# Create a dataset containing model variables
cox_data <- lung_data[, c(
  "time", "event", "age", "sex", "ph.ecog"
)]

# Check missing values
colSums(is.na(cox_data))

# Remove observations missing model variables
cox_data <- na.omit(cox_data)

# Check the number of observations retained
nrow(cox_data)

# Fit multivariable Cox proportional hazards model
cox_model <- coxph(
  Surv(time, event) ~ age + sex + ph.ecog,
  data = cox_data
)

# Display model results
summary(cox_model)

# Hazard ratios
exp(coef(cox_model))

# 95% confidence intervals for hazard ratios
exp(confint(cox_model))

# Test the proportional hazards assumption
ph_test <- cox.zph(cox_model)

print(ph_test)

# Plot Schoenfeld residual diagnostics
plot(ph_test)


# A multivariable Cox proportional hazards model was
# fitted using age, sex, and ECOG performance status.
#
# The model included 227 patients and 164 observed deaths.
# One patient was excluded due to missing ECOG data.
#
# Age:
# HR = 1.011 (95% CI: 0.993 to 1.030), p = 0.232.
# Age was not significantly associated with survival
# after adjusting for sex and ECOG performance status.
#
# Sex:
# HR = 0.575 (95% CI: 0.414 to 0.799), p < 0.001.
# Female patients had approximately 42.5% lower
# estimated hazard of death compared with male patients,
# after adjusting for age and ECOG performance status.
#
# ECOG Performance Status:
# HR = 1.590 (95% CI: 1.273 to 1.986), p < 0.001.
# Each one-point increase in ECOG score was associated
# with approximately 59% higher estimated hazard of death,
# holding the other predictors constant.
#
# The model concordance index was 0.637, indicating
# modest apparent discrimination.
#
# Proportional Hazards Assumption:
# The global Schoenfeld residual test produced p = 0.22.
# Individual tests were also non-significant.
# No statistically significant evidence of a violation
# of the proportional hazards assumption was detected.
#
# These results represent associations rather than
# causal relationships.



# Step 6: Examine ECOG Performance Status


# Number of patients in each ECOG category
table(cox_data$ph.ecog)

# Create a categorical ECOG variable
cox_data$ecog_factor <- factor(cox_data$ph.ecog)

# Fit Cox model with categorical ECOG
cox_ecog_factor <- coxph(
  Surv(time, event) ~ age + sex + ecog_factor,
  data = cox_data
)

# Compare the two model specifications
anova(cox_model, cox_ecog_factor, test = "Chisq")

# Display categorical ECOG model
summary(cox_ecog_factor)



# Step 6: Model Comparison Results


# Two Cox regression models were compared:
# Model 1 included ECOG as a numeric predictor.
# Model 2 included ECOG as a categorical predictor.
#
# A likelihood ratio test produced:
# Chi-square = 0.3661, df = 2, p = 0.8327.
#
# The categorical ECOG model did not provide a
# statistically significant improvement in fit.
#
# The numeric ECOG model was retained because it
# is more parsimonious and its linear effect
# specification was not contradicted by this test.
#
# In the categorical model, higher ECOG scores were
# associated with higher estimated mortality hazards.
# However, the ECOG 3 estimate had a wide confidence
# interval, indicating substantial uncertainty.

# Number of patients by ECOG score
table(cox_data$ph.ecog)

# Number of observed deaths by ECOG score
table(cox_data$ph.ecog, cox_data$event)



# Step 7: ECOG Distribution and Interpretation


# ECOG performance status distribution:
# ECOG 0: 63 patients, 37 deaths
# ECOG 1: 113 patients, 82 deaths
# ECOG 2: 50 patients, 44 deaths
# ECOG 3: 1 patient, 1 death
#
# Most patients had ECOG scores of 0 or 1.
#
# Only one patient had an ECOG score of 3,
# explaining the substantial uncertainty in the
# categorical model's ECOG 3 hazard ratio.
#
# The categorical ECOG model did not significantly
# improve fit compared with the numeric model.
#
# Therefore, the numeric ECOG model was retained
# as the primary model.
#
# The sparse ECOG 3 category is a limitation,
# and the results should be interpreted cautiously.


## Step 8 — Model validation

#We've completed the main survival analyses. Next, let's evaluate whether our Cox model's discrimination holds up when tested on patients not used to fit it.

#Rather than using a single random train/test split on only 227 patients, I recommend **bootstrap internal validation**. It will help us estimate how optimistic our in-sample concordance of 0.637 might be.

#For now, we'll use 200 bootstrap samples as an initial assessment.


# Step 8: Bootstrap Internal Validation




set.seed(123)

B <- 200
n <- nrow(cox_data)

# Apparent concordance from the original model
apparent_c <- unname(
  summary(cox_model)$concordance[1]
)

optimism <- numeric(B)

for (b in seq_len(B)) {
  
  # Bootstrap sample with replacement
  sample_id <- sample(seq_len(n), n, replace = TRUE)
  boot_data <- cox_data[sample_id, ]
  
  # Fit model to bootstrap sample
  boot_model <- coxph(
    Surv(time, event) ~ age + sex + ph.ecog,
    data = boot_data
  )
  
  # Performance on bootstrap sample
  c_boot <- unname(
    concordance(boot_model)$concordance
  )
  
  # Apply bootstrap model to original dataset
  risk_original <- predict(
    boot_model,
    newdata = cox_data,
    type = "lp",
    reference = "zero"
  )
  
  c_original <- unname(
    concordance(
      Surv(time, event) ~ risk_original,
      data = cox_data,
      reverse = TRUE
    )$concordance
  )
  
  optimism[b] <- c_boot - c_original
}

# Average optimism
mean_optimism <- mean(optimism)

# Optimism-corrected concordance
corrected_c <- apparent_c - mean_optimism

cat("Apparent C-index:", round(apparent_c, 3), "\n")
cat("Mean optimism:", round(mean_optimism, 3), "\n")
cat("Optimism-corrected C-index:", round(corrected_c, 3), "\n")



# Step 9: Sensitivity Analysis


# Objective:
# Investigate whether the single patient with an
# ECOG score of 3 substantially influences the
# results of the original Cox regression model.

# Exclude the patient with ECOG score 3.
sensitivity_data <- subset(cox_data, ph.ecog <= 2)

# Check the number of remaining observations.
nrow(sensitivity_data)

# Fit the Cox model using the reduced dataset.
cox_sensitivity <- coxph(
  Surv(time, event) ~ age + sex + ph.ecog,
  data = sensitivity_data
)

# Display the sensitivity model results.
summary(cox_sensitivity)

# Compare hazard ratios from the original model
# and the sensitivity model.
comparison <- data.frame(
  Predictor = names(coef(cox_model)),
  Original_HR = round(exp(coef(cox_model)), 3),
  Sensitivity_HR = round(exp(coef(cox_sensitivity)), 3)
)

# Display the comparison table.
print(comparison)



# Step 9: Results and Interpretation


# A sensitivity analysis was conducted by excluding
# the single patient with an ECOG score of 3.
#
# The original model included 227 patients, while
# the sensitivity model included 226 patients.
#
# The estimated hazard ratios were very similar:
#
# Age:
# Original HR = 1.011
# Sensitivity HR = 1.011
#
# Sex (Female vs. Male):
# Original HR = 0.575
# Sensitivity HR = 0.578
#
# ECOG Performance Status:
# Original HR = 1.590
# Sensitivity HR = 1.573
#
# The sensitivity model had a concordance index
# of 0.635, compared with 0.637 in the original model.
#
# Excluding the single ECOG 3 patient did not
# materially change the estimated hazard ratios
# or the model's apparent discrimination.
#
# These findings suggest that the main conclusions
# are not driven by this single observation.
#
# However, this analysis does not eliminate other
# potential sources of bias or model uncertainty.



