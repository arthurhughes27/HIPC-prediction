balance_folds <- function(df,
                          ind.col,
                          covariate.cols,
                          n.folds = 10,
                          n.continuous.split = 5) {
  # df: one row per individual
  # ind.col: column name of participant identifier
  # covariate.cols: character vector of covariates to balance on
  # returns: data.frame with ind.col and fold assignment
  
  participants <- df[[ind.col]]
  
  # if no covariates, simple random balanced assignment
  if (length(covariate.cols) == 0) {
    shuffled <- sample(participants)
    folds <- ((seq_along(shuffled) - 1) %% n.folds) + 1
    return(data.frame(participant = participants, fold = folds))
  }
  
  # create strata for each covariate
  strata_list <- lapply(covariate.cols, function(col) {
    vec <- df[[col]]
    if (is.numeric(vec) && length(unique(vec)) > 2) {
      # continuous -> quantile bins
      brks <- unique(quantile(
        vec,
        probs = seq(0, 1, length.out = n.continuous.split + 1),
        na.rm = TRUE
      ))
      if (length(brks) <= 1) {
        strata <- rep("1", length(vec))
      } else {
        strata <- as.character(cut(
          vec,
          breaks = brks,
          include.lowest = TRUE,
          labels = FALSE
        ))
      }
    } else {
      # categorical / binary
      strata <- as.character(vec)
    }
    ifelse(is.na(strata), "NA", strata)
  })
  
  # combine strata across covariates
  combined_strata <- apply(do.call(data.frame, strata_list), 1, paste, collapse = "___")
  
  # assign folds within each stratum
  fold_map <- integer(length(participants))
  for (stratum in unique(combined_strata)) {
    ids <- which(combined_strata == stratum)
    shuffled <- sample(ids)
    fold_map[shuffled] <- ((seq_along(shuffled) - 1) %% n.folds) + 1
  }
  
  res = data.frame(participant = participants, fold = fold_map)
  
  return(res)
}