# R Script to perform baseline model (linear regression of post on pre + age + sex) 
# in a cross-validation framework

# Libraries
library(fs)
library(dplyr)
library(stringr)
library(janitor)
library(readxl)
library(purrr)

processed_data_folder = "data"

# Use fs::path() to specify the data path robustly
p_load_hipc_immResp <- fs::path(processed_data_folder, "hipc_immResp.rds")
p_load_hipc_clinical <- fs::path(processed_data_folder, "hipc_clinical.rds")

# Load dataframes
hipc_immResp = readRDS(p_load_hipc_immResp) %>% 
  filter(study_accession == "SDY1276")

hipc_clinical = readRDS(p_load_hipc_clinical) %>% 
  filter(study_accession == "SDY1276")

# Join together the age and sex columns 
hipc_demographics = hipc_clinical %>% 
  mutate(genderMale = ifelse(gender == "Male", 1, 0)) %>% 
  select(participant_id, 
         age_imputed,
         genderMale) %>% 
  distinct()

hipc_immResp_merged = merge(x = hipc_demographics,
                            y = hipc_immResp,
                            by = "participant_id") 

# Cross-validation framework

# Parameters
set.seed(1234567)   # reproducible seed; change if you want different fold assignment
K <- 10           # number of folds

# Columns (as requested)
resp_col <- "immResp_MFC_nAb_log2_post_value"
pred_cols <- c("immResp_MFC_nAb_log2_pre_value", "age_imputed", "genderMale")

# 1) Prepare data: keep only complete cases for the model variables
df <- hipc_immResp_merged %>%
  dplyr::select(all_of(c("participant_id", resp_col, pred_cols))) %>%
  filter(!is.na(.data[[resp_col]])) %>%
  filter(complete.cases(.[pred_cols])) %>%
  # keep row index to map predictions back if needed
  mutate(.row = row_number())

n <- nrow(df)
if (n < K) stop("Number of complete cases is less than number of folds K.")

# 2) Create random folds
folds <- sample(rep(1:K, length.out = n))  # approximate equal sizes
df$cv_fold <- folds

# 3) Cross-validated predictions
pred_vec <- rep(NA_real_, n)  # placeholder for CV predictions

for (k in seq_len(K)) {
  train_idx <- which(df$cv_fold != k)
  test_idx  <- which(df$cv_fold == k)
  
  train_df <- df[train_idx, ]
  test_df  <- df[test_idx, ]
  
  # fit linear model on training data
  formula <- as.formula(paste(resp_col, "~", paste(pred_cols, collapse = " + ")))
  fit <- lm(formula, data = train_df)
  
  # predict on held-out fold
  pred_vec[test_idx] <- predict(fit, newdata = test_df)
}

# attach predictions to dataframe
df$pred_cv <- pred_vec
df_obs <- df[[resp_col]]
df_pred <- df$pred_cv

# 4) Compute metrics: RMSE, standardized RMSE (RMSE / sd(obs)), and CV R-squared
# RMSE
rmse <- sqrt(mean((df_obs - df_pred)^2, na.rm = TRUE))
# standardized RMSE: defined here as RMSE divided by sd(observed)
sd_obs <- sd(df_obs, na.rm = TRUE)
sRMSE <- rmse / sd_obs

# Cross-validated R^2: 1 - SSE / SST (using overall observed mean)
sse <- sum((df_obs - df_pred)^2, na.rm = TRUE)
sst <- sum((df_obs - mean(df_obs, na.rm = TRUE))^2, na.rm = TRUE)
R2_cv <- 1 - sse / sst

# Round stats for annotation
rmse_r <- signif(rmse, 3)
sRMSE_r <- signif(sRMSE, 3)
R2_r <- signif(R2_cv, 3)

# 5) Build observed vs predicted plot with identical scales and diagonal
# compute axis limits with a small margin
all_vals <- range(c(df_obs, df_pred), na.rm = TRUE)
margin <- 0.02 * diff(all_vals)
xlim_plot <- c(all_vals[1] - margin, all_vals[2] + margin)

p_cv <- ggplot(df, aes(x = .data[[resp_col]], y = pred_cv)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.6) +
  coord_equal(xlim = xlim_plot, ylim = xlim_plot) +   # identical scales and limits
  labs(
    x = "Observed (post) immResp_MFC_nAb_log2",
    y = "CV predicted immResp_MFC_nAb_log2",
    title = "Observed vs cross-validated predicted values"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  # annotate sRMSE and R2 in top-left corner of the plotting region
  annotate(
    "text",
    x = xlim_plot[1] + 0.02 * diff(xlim_plot),
    y = xlim_plot[2] - 0.02 * diff(xlim_plot),
    label = paste0("sRMSE = ", sRMSE_r, "\nR²_CV = ", R2_r),
    hjust = 0, vjust = 1, size = 4
  )

print(p_cv)

# 6) Print numeric summary to console
message("Cross-validated performance:")
message("  RMSE      = ", rmse_r)
message("  sRMSE     = ", sRMSE_r, "  (RMSE / sd(observed))")
message("  R2 (CV)   = ", R2_r)
