
# Lung Cancer Survival Analysis in R

**Biostatistics | Survival Analysis | Cox Regression | R**

## Overview

This project investigates survival outcomes among patients with advanced lung cancer using statistical survival analysis methods in R.

The objective is to examine associations between patient characteristics and overall survival while accounting for right-censored observations.

The analysis includes Kaplan-Meier estimation, log-rank testing, multivariable Cox proportional hazards regression, model diagnostics, bootstrap internal validation, and sensitivity analysis.

## Dataset

The project uses the North Central Cancer Treatment Group (NCCTG) lung cancer dataset, available through the `survival` package in R.

| Characteristic | Value |
|---|---|
| Total patients | 228 |
| Observed deaths | 165 |
| Censored observations | 63 |
| Mean age | 62.45 years |
| Male patients | 138 |
| Female patients | 90 |

The primary variables include:

- `time`: Survival time in days
- `status`: Original censoring/event status
- `event`: Binary event indicator (1 = death, 0 = censored)
- `age`: Patient age
- `sex`: Patient sex
- `ph.ecog`: ECOG performance status

The dataset contains missing observations in several clinical variables. The primary Cox regression analysis uses 227 complete cases.

## Statistical Methods

### 1. Exploratory Data Analysis

Examined patient characteristics, survival times, event frequencies, and missing observations.

### 2. Kaplan-Meier Survival Estimation

Estimated the overall survival function and corresponding 95% confidence intervals.

Survival probabilities were evaluated at 180 and 365 days.

### 3. Log-Rank Hypothesis Testing

Compared survival distributions between male and female patients using Kaplan-Meier curves and the log-rank test.

### 4. Cox Proportional Hazards Regression

Fitted a multivariable Cox regression model using:

- Age
- Sex
- ECOG performance status

Estimated adjusted hazard ratios, confidence intervals, and statistical significance.

### 5. Model Diagnostics

Evaluated the proportional hazards assumption using Schoenfeld residuals and the `cox.zph()` function.

Compared numeric and categorical ECOG specifications using a likelihood ratio test.

### 6. Bootstrap Internal Validation

Performed 200 bootstrap resamples to estimate optimism in the model's concordance index.

### 7. Sensitivity Analysis

Refitted the Cox model after excluding the single patient with ECOG performance status 3.

Compared hazard ratios to assess the influence of this sparse category.

## Results

### Overall Survival

| Measure | Estimate |
|---|---|
| Median survival | 310 days |
| Median survival 95% CI | 285–363 days |
| 180-day survival probability | 72.2% |
| 365-day survival probability | 40.9% |

The Kaplan-Meier estimator indicated that approximately 40.9% of patients survived beyond one year.

### Survival Differences by Sex

Kaplan-Meier estimates indicated higher survival probabilities among female patients.

The log-rank test identified a statistically significant difference between the survival functions.

**Log-rank p-value: 0.00131**

This represents an unadjusted association and does not establish a causal effect.

### Multivariable Cox Regression

The primary model included 227 patients and 164 observed deaths.

| Predictor | Adjusted HR | 95% CI | p-value |
|---|---:|---|---:|
| Age (per year) | 1.011 | 0.993–1.030 | 0.232 |
| Female vs. male | 0.575 | 0.414–0.799 | 0.000986 |
| ECOG (per point) | 1.590 | 1.273–1.986 | 0.0000445 |

**Interpretation:**

- Female sex was associated with approximately 42.5% lower estimated mortality hazard after adjusting for age and ECOG performance status.
- Each one-point increase in ECOG performance status was associated with approximately 59% higher estimated mortality hazard, assuming a log-linear relationship.
- Age was not statistically significant in the adjusted model.

These findings describe statistical associations rather than causal relationships.

### Proportional Hazards Diagnostics

The Schoenfeld residual tests produced the following results:

| Variable | p-value |
|---|---:|
| Age | 0.66 |
| Sex | 0.13 |
| ECOG | 0.15 |
| Global test | 0.22 |

The global test did not identify statistically significant evidence against the proportional hazards assumption.

However, nonsignificant tests do not prove that the assumption holds.

### ECOG Model Comparison

Numeric and categorical ECOG specifications were compared using a likelihood ratio test.

**Likelihood ratio test p-value: 0.833**

The categorical specification did not provide statistically significant evidence of improved model fit.

The numeric ECOG specification was retained for parsimony.

### Bootstrap Internal Validation

The Cox model was internally validated using 200 bootstrap resamples.

| Performance Metric | Estimate |
|---|---:|
| Apparent C-index | 0.637 |
| Mean optimism | 0.006 |
| Optimism-corrected C-index | 0.631 |

The small estimated optimism suggests limited apparent overfitting.

However, the optimism-corrected C-index of 0.631 indicates modest discriminatory ability.

This procedure provides internal validation and does not establish performance on independent patient data.

### Sensitivity Analysis

The single patient with ECOG score 3 was excluded, and the Cox model was refitted using 226 observations.

| Predictor | Original HR | Sensitivity HR |
|---|---:|---:|
| Age | 1.011 | 1.011 |
| Female vs. male | 0.575 | 0.578 |
| ECOG | 1.590 | 1.573 |

The sensitivity model produced a concordance index of 0.635.

The hazard ratios remained similar, suggesting that the primary findings were not materially driven by the single ECOG 3 observation.

## Limitations

Several limitations should be considered:

- The dataset contains a relatively small historical cohort of patients with advanced lung cancer.
- Findings may not generalize to contemporary or broader lung cancer populations.
- The observational study design does not support causal conclusions.
- The primary model includes a limited number of predictors.
- The ECOG 3 category contains only one patient.
- Complete-case analysis excludes one observation from the primary regression model.
- Internal bootstrap validation does not replace external validation.
- Model calibration and potential nonlinear predictor effects were not formally evaluated.

The model is intended for statistical analysis and educational research rather than clinical prediction.

## Technologies

- **Programming Language:** R
- **Development Environment:** RStudio
- **Package:** `survival`
- **Methods:** Kaplan-Meier estimation, log-rank testing, Cox regression, Schoenfeld residual diagnostics, likelihood ratio testing, bootstrap resampling

## Repository Structure

```text
lung_cancer_survival_analysis/
├── cancer_research_rdataset.R
└── README.md
```

## Reproducibility

### Requirements

- R
- RStudio (recommended)
- `survival` package

### Installation

Install the required R package:

```r
install.packages("survival")
```

### Running the Analysis

1. Clone the repository:

```bash
git clone https://github.com/alaa-masaoud/lung_cancer_survival_analysis.git
```

2. Open `cancer_research_rdataset.R` in RStudio.

3. Run the script from beginning to end.

The dataset is included in the `survival` package, so no external data download is required.

The bootstrap analysis uses a fixed random seed for reproducibility.

## Future Work

Potential extensions include:

- Investigating nonlinear age effects
- Assessing survival probability calibration
- Evaluating additional clinical predictors
- Validating the model using an independent patient cohort
- Extending the analysis to larger cancer registry datasets

## Disclaimer

This project is intended for educational and statistical research purposes. It is not a clinically validated prognostic tool.
