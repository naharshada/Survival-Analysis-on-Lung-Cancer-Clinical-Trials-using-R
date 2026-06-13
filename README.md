# Survival-Analysis-on-Lung-Cancer-Clinical-Trials-using-R
#package installation
install.packages(c("survival","survminer","dplyr","ggplot2","gtsummary","car","rmarkdown","knitr"))
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(gtsummary)
#loading dataset
df=read.csv(file.choose(),header = T)
df
colnames(df)
dim(df)         
head(df)         
str(df)          
summary(df)      
#Cleaning data
colSums(is.na(df))

df <- df[, -which(names(df) == "X")]
df

df$sex=factor(df$sex,levels = c(1,2),labels = c("male","female"))
df$sex
df$ph.ecog <- factor(df$ph.ecog)
df[complete.cases(df), ]
df_clean <- df %>%
  select(time, status, age, sex, ph.ecog, ph.karno, pat.karno, wt.loss) %>%
  filter(complete.cases(.))

nrow(df_clean)  
df_clean

#Descriptive Statistics
table1 <- df %>%
  select(age, sex, ph.ecog, ph.karno, pat.karno,
         meal.cal, wt.loss, time, status) %>%
  tbl_summary(
    by = sex,
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = list(
      age      ~ "Age (years)",
      ph.ecog  ~ "ECOG Performance Score",
      ph.karno ~ "Karnofsky Score (Physician)",
      wt.loss  ~ "Weight Loss (lbs)"
    )
  ) %>%
  add_p() %>%         
  add_overall()      
table1   

#Survival function
surv_obj= Surv(time = df$time,
               event = df$status == 2)
head(surv_obj)
#fitting KM Curve
km_overall <- survfit(surv_obj ~ 1, data = df)
km_overall
ggsurvplot(km_overall,
           data          = df,
           conf.int      = TRUE,       # shaded 95% CI band
           risk.table    = TRUE,       # patients at risk table below
           xlab          = "Time (days)",
           ylab          = "Survival Probability",
           title         = "Overall Survival — NCCTG Lung Cancer",
           surv.median.line = "hv"    # draw median survival line
)
#KM stratified by sex
km_sex <- survfit(surv_obj ~ sex, data = df)
km_sex
# Log-rank test
log_rank <- survdiff(surv_obj ~ sex, data = df)
print(log_rank)
ggsurvplot(km_sex,
           data         = df,
           pval         = TRUE,         # show p-value on plot
           conf.int     = TRUE,
           risk.table   = TRUE,
           palette      = c("#378ADD", "#D4537E"),
           legend.labs  = c("Male", "Female"),
           title        = "Survival by Sex"
)

#univariate cox regression
surv_obj <- Surv(time   = df_clean$time,
                 event  = df_clean$status == 2)
surv_obj
library(broom)

vars <- c("age", "sex", "ph.ecog", "ph.karno", "pat.karno", "wt.loss")

univ_results <- lapply(vars, function(v) {
  formula <- as.formula(paste("surv_obj ~", v))
  model   <- coxph(formula, data = df_clean)
  tidy(model, exponentiate = TRUE, conf.int = TRUE)
})

univ_table <- bind_rows(univ_results)
print(univ_table)
# Multivariate Cox model
cox_model <- coxph(
  surv_obj ~ age + sex + ph.ecog + ph.karno + wt.loss,
  data = df_clean
)

# View detailed results
summary(cox_model)

# Get tidy HR table with confidence intervals
library(broom)
cox_results <- tidy(cox_model,
                    exponentiate = TRUE,  # gives HRs not log-HRs
                    conf.int     = TRUE)
print(cox_results)

# Create publication-ready table
tbl_regression(cox_model, exponentiate = TRUE)
# Forest plot using ggforest from survminer
ggforest(cox_model,
         data      = df_clean,
         main      = "Hazard Ratios — Multivariate Cox Model",
         cpositions = c(0.02, 0.22, 0.4),
         fontsize  = 0.8,
         refLabel  = "Reference",
         noDigits  = 2
)

# Test PH assumption using Schoenfeld residuals
ph_test <- cox.zph(cox_model)
print(ph_test)
# For each variable: p > 0.05 means PH assumption holds

# Visual check — residual plots
ggcoxzph(ph_test)
# Want: flat horizontal line (no trend)
# Bad: clear upward or downward slope

# Plot deviance residuals
ggcoxdiagnostics(cox_model,
                 type        = "deviance",
                 linear.predictions = FALSE
)
