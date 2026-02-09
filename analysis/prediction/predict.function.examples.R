sapply(list.files("R/", pattern = "\\.R$", full.names = TRUE), source)

spec = list(
  study_of_interest = "SDY1276",
  dataset = "all_noNorm",
  data.sel = "d1",
  include.cov = TRUE,
  covariate.cols = c(
    "genderMale",
    "age_imputed",
    "immResp_MFC_nAb_log2_pre_value"
  ),
  response.col = "immResp_MFC_nAb_log2_post_value",
  seed = 21012026,
  K = 10,
  mod = "ranger",
  feat.eng.col = "none",
  feat.eng.row = "mean",
  feat.select = "none",
  feature.select.metric = "sRMSE",
  feature.select.metric.threshold = 0,
  feature.select.model = "lm",
  feature.select.criterion = "relative.gain",
  n.cores = 10
)

res_specification <- do.call(cv.predict.specification, spec)

res_baseline <- do.call(cv.predict.baseline.specification, spec)
