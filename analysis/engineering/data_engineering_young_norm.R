# R script to engineer expression data

# Libraries
library(dplyr)
library(GSVA)
library(fs)   # used for path construction (fs::path used below)
library(tidyr)
library(devtools)
load_all("/home/ah3/Desktop/Work/PhD/dearseq/dearseq-devel/")


# Directory containing data files
processed_data_folder <- "data"

# Paths to processed gene-level data and gene-set objects
p_load_expr_young_norm <- fs::path(processed_data_folder, "hipc_merged_young_norm.rds")
p_load_btm           <- fs::path(processed_data_folder, "BTM_processed.rds")

# Load data objects
hipc_merged_young_norm <- readRDS(p_load_expr_young_norm)
BTM                  <- readRDS(p_load_btm)
idx <- which(BTM[["geneset.descriptions"]] == "TBA")

BTM[["genesets"]] <- BTM[["genesets"]][-idx]
BTM[["geneset.descriptions"]] <- BTM[["geneset.descriptions"]][-idx]
BTM[["geneset.names"]] <- BTM[["geneset.names"]][-idx]
BTM[["geneset.aggregates"]] = BTM[["geneset.aggregates"]][-idx]
BTM[["geneset.names.descriptions"]] = BTM[["geneset.names.descriptions"]][-idx]

# Timepoints of interest (numeric)
timepoints_of_interest <- c(0, 1, 3)

studies_of_interest = c("SDY1276")

# Filter to samples with non-missing immune response, Influenza vaccine,
# and collected at one of the specified timepoints.
hipc_merged_young_norm_filtered <- hipc_merged_young_norm %>%
  filter(
    !is.na(immResp_MFC_anyAssay_log2_MFC),
    vaccine_name == "Influenza (IN)",
    study_time_collected %in% timepoints_of_interest,
    study_accession %in% studies_of_interest
  )

gene_names = hipc_merged_young_norm_filtered %>%
  dplyr::select(a1cf:zzz3) %>%
  colnames()

# df = hipc_merged_young_norm_filtered
genesets = BTM[["genesets"]]
geneset_names = BTM[["geneset.names.descriptions"]]
id_col = "participant_id"
time_col = "study_time_collected"
gene_cols = gene_names
# transformation = "rank"
# timepoint = 1


compute_gene_fc = function(df, id_col, time_col, gene_cols, timepoint) {
  df = df %>%
    arrange(.data[[id_col]])
  
  df_filtered_baseline = df %>%
    filter(.data[[time_col]] == 0)
  
  df_filtered_post = df %>%
    filter(.data[[time_col]] == timepoint)
  
  participants = intersect(df_filtered_baseline[[id_col]], df_filtered_post[[id_col]])
  
  df_filtered_baseline = df_filtered_baseline %>%
    filter(.data[[id_col]] %in% participants) %>%
    arrange(.data[[id_col]])
  
  df_filtered_post = df_filtered_post %>%
    filter(.data[[id_col]] %in% participants) %>%
    arrange(.data[[id_col]])
  
  expr_baseline = df_filtered_baseline %>%
    select(any_of(gene_cols)) %>%
    as.matrix()
  
  expr_post = df_filtered_post %>%
    select(any_of(gene_cols)) %>%
    as.matrix()
  
  expr_fc = expr_post - expr_baseline
  
  result <- cbind(participant_id = df_filtered_post[[id_col]],
                  as.data.frame(expr_fc, stringsAsFactors = FALSE))
  
  return(result)
}



compute_col_transformation = function(df,
                                      id_col,
                                      time_col,
                                      gene_cols,
                                      timepoint,
                                      transformation = "none") {
  df_filtered = df %>%
    filter(.data[[time_col]] == timepoint)
  
  expr_mat = df_filtered %>%
    select(any_of(gene_cols)) %>%
    as.matrix()
  
  transformation_function = function(col, transformation) {
    if (transformation == "none") {
      return(col)
    } else if (transformation == "rank") {
      return(rank(-col, ties.method = "average"))
    } else if (transformation == "z") {
      return(as.numeric(scale(
        col, center = TRUE, scale = TRUE
      )))
    }
  }
  
  expr_mat_transformed <- apply(expr_mat, 2, FUN = transformation_function, transformation = transformation)
  
  result <- cbind(participant_id = df_filtered[[id_col]],
                  as.data.frame(expr_mat_transformed, stringsAsFactors = FALSE))
  
  return(result)
}


compute_row_transformation = function(df,
                                      id_col,
                                      time_col,
                                      gene_cols,
                                      genesets,
                                      geneset_names,
                                      timepoint,
                                      transformation = "mean") {
  df_filtered = df %>%
    filter(.data[[time_col]] == timepoint)
  
  expr_mat = df_filtered %>%
    select(any_of(gene_cols)) %>%
    as.matrix()
  
  if (transformation == "none") {
    transformed_mat = expr_mat
    
    return(transformed_mat)
  }
  
  if (transformation == "ssgsea") {
    names(genesets) <- geneset_names
      params <- ssgseaParam(
        exprData = t(expr_mat),
        geneSets = genesets,
        minSize = 1,
        maxSize = Inf,
        normalize = TRUE
      )
      transformed_mat <- t(gsva(params))
    
    
    transformed_mat <- as.data.frame(transformed_mat, stringsAsFactors = FALSE)
    transformed_mat <- cbind(participant_id = df_filtered[[id_col]], transformed_mat)
    
    return(transformed_mat)
  }
  
  
  transformation_function = function(row, transformation) {
    if (transformation == "mean") {
      return(mean(row))
    } else if (transformation == "median") {
      return(median(row))
    } else if (transformation == "max") {
      return(max(row))
    } else if (transformation == "cv") {
      return(ifelse(mean(row) == 0, NA_real_, sd(row) / mean(row)))
    } else if (transformation == "iqr") {
      return(IQR(row))
    }
  }
  
  transformed_mat = data.frame(matrix(
    NA_real_,
    nrow = nrow(expr_mat),
    ncol = length(genesets)
  ))
  
  colnames(transformed_mat) = geneset_names
  
  for (i in seq_along(genesets)) {
    gs.genes = intersect(colnames(expr_mat), genesets[[i]])
    
    expr_mat_gs = expr_mat[, gs.genes, drop = FALSE]
    
    if (ncol(expr_mat_gs) == 0) {
      transformed <- rep(NA_real_, nrow(expr_mat_gs))
    } else if (transformation == "pc1") {
      pca_res <- stats::prcomp(expr_mat_gs, center = TRUE, scale. = TRUE)
      transformed <- as.numeric(pca_res$x[, 1])
    } else {
      transformed <- apply(expr_mat_gs,
                           1,
                           FUN = transformation_function,
                           transformation = transformation)
      transformed <- as.numeric(transformed)
    }
    
    transformed_mat[[geneset_names[i]]] <- transformed
  }
  
  transformed_mat <- cbind(participant_id = df_filtered[[id_col]], transformed_mat)
  
  return(transformed_mat)
}


apply_transformations = function(df,
                                 id_col,
                                 time_col,
                                 gene_cols,
                                 genesets,
                                 geneset_names,
                                 timepoint,
                                 steps) {
  # start from the filtered data (rows to keep)
  df_filtered <- df %>% filter(.data[[time_col]] == timepoint)
  df_work <- df_filtered
  current_gene_cols <- gene_cols
  
  for (step in steps) {
    if (!is.list(step) ||
        is.null(step$type) || is.null(step$transformation)) {
      stop("Each step must be a list with elements $type and $transformation")
    }
    
    if (step$type == "col") {
      mat <- compute_col_transformation(
        df_work,
        id_col,
        time_col,
        current_gene_cols,
        timepoint,
        transformation = step$transformation
      )
      # mat already contains participant_id; ensure time_col present for subsequent filtering
      df_work <- as.data.frame(mat, stringsAsFactors = FALSE)
      df_work[[time_col]] <- timepoint
      current_gene_cols <- setdiff(colnames(df_work), c(id_col, time_col))
    } else if (step$type == "row") {
      res <- compute_row_transformation(
        df_work,
        id_col,
        time_col,
        current_gene_cols,
        genesets,
        geneset_names,
        timepoint,
        transformation = step$transformation
      )
      df_work <- as.data.frame(res, stringsAsFactors = FALSE)
      df_work[[time_col]] <- timepoint
      current_gene_cols <- setdiff(colnames(df_work), c(id_col, time_col))
    } else {
      stop("Step type must be 'col' or 'row'")
    }
  }
  
  # final data.frame with participant_id and time_col preserved
  final_df <- df_work[, c(id_col, time_col, current_gene_cols), drop = FALSE]
  return(final_df)
}


d0 = hipc_merged_young_norm_filtered %>%
  filter(study_time_collected == 0) %>%
  select(participant_id, any_of(gene_names))

d1 = hipc_merged_young_norm_filtered %>%
  filter(study_time_collected == 1) %>%
  select(participant_id, any_of(gene_names))

d3 = hipc_merged_young_norm_filtered %>%
  filter(study_time_collected == 3) %>%
  select(participant_id, any_of(gene_names))


col_transformations <- c("none", "z", "rank")   # matches compute_col_transformation's accepted strings
row_transformations <- c("none", "mean", "median", "max", "iqr", "cv", "ssgsea", "pc1")

# precompute fold-change tables (these lack the time column; we'll add it before passing to functions)
d1fc <- compute_gene_fc(
  df = hipc_merged_young_norm_filtered,
  id_col = id_col,
  time_col = time_col,
  gene_cols = gene_cols,
  timepoint = 1
)

d3fc <- compute_gene_fc(
  df = hipc_merged_young_norm_filtered,
  id_col = id_col,
  time_col = time_col,
  gene_cols = gene_cols,
  timepoint = 3
)

# add a time column so the existing functions (which filter by timepoint) work with these dfs
d1fc_time <- d1fc
d1fc_time[[time_col]] <- 1L

d3fc_time <- d3fc
d3fc_time[[time_col]] <- 3L

# mapping of top-level datasets to source data and timepoint
sources <- list(
  d0   = list(df = hipc_merged_young_norm_filtered, tp = 0L),
  d1   = list(df = hipc_merged_young_norm_filtered, tp = 1L),
  d3   = list(df = hipc_merged_young_norm_filtered, tp = 3L),
  d1fc = list(df = d1fc_time, tp = 1L),
  d3fc = list(df = d3fc_time, tp = 3L)
)

engineered <- list()

for (top_name in names(sources)) {
  src <- sources[[top_name]]$df
  tp  <- sources[[top_name]]$tp
  message("Processing dataset: ", top_name, " (timepoint ", tp, ")")
  engineered[[top_name]] <- list()
  
  for (col_tf in col_transformations) {
    message("  Column transformation: ", col_tf)
    engineered[[top_name]][[col_tf]] <- list()
    
    for (row_tf in row_transformations) {
      message("    Row transformation: ", row_tf)
      
      # build steps according to the rule:
      steps <- list()
      if (row_tf %in% c("ssgsea", "pc1")) {
        steps <- append(steps, list(list(
          type = "row", transformation = row_tf
        )))
        if (col_tf != "none")
          steps <- append(steps, list(list(
            type = "col", transformation = col_tf
          )))
      } else {
        if (col_tf != "none")
          steps <- append(steps, list(list(
            type = "col", transformation = col_tf
          )))
        if (row_tf != "none")
          steps <- append(steps, list(list(
            type = "row", transformation = row_tf
          )))
      }
      
      res_df <- apply_transformations(
        df = src,
        id_col = id_col,
        time_col = time_col,
        gene_cols = gene_cols,
        genesets = genesets,
        geneset_names = geneset_names,
        timepoint = tp,
        steps = steps
      )
      
      engineered[[top_name]][[col_tf]][[row_tf]] <- res_df
    }
  }
}



p_save <- fs::path(processed_data_folder,
                   "engineered_dataframes_influenzain_young_norm.rds")

# Save the data
saveRDS(engineered, p_save)
