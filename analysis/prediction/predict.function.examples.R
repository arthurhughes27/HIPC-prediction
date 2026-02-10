sapply(list.files("R/", pattern = "\\.R$", full.names = TRUE), source)

spec = list(
  study_of_interest = "SDY1276",
  dataset = "all_noNorm",
  data.sel = "d1",
  include.cov = T,
  covariate.cols = c(
    # "genderMale",
    "age_imputed"
    # ,"immResp_mean_hai_log2_pre_value"
  ),
  response.col = "immResp_mean_hai_log2_post_value",
  seed = 21012026,
  K = 5,
  mod = "elasticnet",
  feat.eng.col = "none",
  feat.eng.row = "none",
  feat.select = "univariate",
  feat.select.metric = "sRMSE",
  feat.select.metric.threshold = 2,
  feat.select.model = "lm",
  feat.select.criterion = "relative.gain",
  n.cores = 1,
  gender.select = "Male"
)

list2env(spec, envir = .GlobalEnv)

res_specification <- do.call(cv.predict.specification, spec)

# spec$mod = "lm"
# 
# res_baseline <- do.call(cv.predict.baseline.specification, spec)

# p1 = res_specification$prediction.plot
# 
# p2 = res_baseline$prediction.plot
# 
# diagnosis_figures_folder = fs::path("output", "figures", "diagnosis")
# 
# ggsave(
#   filename = "elasticnet_nocov.pdf",
#   path = diagnosis_figures_folder,
#   plot = p1,
#   width = 25,
#   height = 15,
#   units = "cm"
# )
# 
# ggsave(
#   filename = "lm_gender.pdf",
#   path = diagnosis_figures_folder,
#   plot = p2,
#   width = 25,
#   height = 15,
#   units = "cm"
# )
