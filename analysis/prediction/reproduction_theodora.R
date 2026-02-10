# Libraries
library(dplyr)
library(tidyverse)
library(readxl)
library(openxlsx)
library(BiocManager)
library(ggplot2)
library(ggrepel)
library(edgeR)
library(mixOmics)
library(FactoMineR)
library(patchwork)
library(xtable)
library(purrr)
library(gtsummary)
library(gt)
library(RColorBrewer)
library(stringr)
library(pbmcapply)
library(SurrogateRank)
library(caTools)
library(Rsurrogate)
library(ggvenn)
library(VennDiagram)
library(ggpubr)
library(caret)

# Loading the data
## Define the paths
raw_data_path = fs::path("data-raw")
path_ge = fs::path(raw_data_path, "all_norm_withResponse_eset.rds")
path_nab = fs::path(raw_data_path, "neut_ab_titer_2025-01-10_01-13-22.xlsx")

# Load the data
data = readRDS(path_ge)
neut.ab.assay <- read_excel(path_nab)

# Extract the clinical data
rawdata <- as.data.frame(data@phenoData@data)
# Extract the gene expression data
ge.data <- as.data.frame(data@assayData[["exprs"]])
rm(data)

# Preprocessing
## Select only females with 2008 TIV vaccine with gene expression info at baseline and 1 day, from SDY1276
rawdata <- rawdata %>%
  filter(
    pathogen == "Influenza",
    vaccine == "TIV (2008)",
    gender == "Female",
    study_time_collected %in% c(0, 1),
    study_accession == "SDY1276"
  ) %>%
  rename(expression_time = study_time_collected) %>%
  group_by(participant_id) %>% # keep only participants with Both 0 and 1
  filter(n_distinct(expression_time) == 2) %>%
  ungroup()

# Extract participant ids from filtered data
filtered_ids = rawdata %>%
  pull(participant_id) %>%
  unique()

# Filter particpants from neutralizing antibody dataframe
## Only choose participants with antibody levels at day 0 and 28 without missing values
nab_filtered <- neut.ab.assay %>%
  filter(
    `Participant ID` %in% filtered_ids,
    `Study Time Collected` %in% c(0, 28),
    !is.na(`Value Preferred`)
  )

# Grab the average of the antibody response
nab_composite <- nab_filtered %>%
  group_by(`Participant ID`, `Study Time Collected`) %>%
  summarize(response = mean(`Value Preferred`, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(`Participant ID`) %>%
  filter(n_distinct(`Study Time Collected`) == 2) %>%
  ungroup()

# Match expression days with response days
## We match d0 gene expression levels with d0 Ab levels, and d1 gene expression levels with d28 Ab levels
nab_final <- nab_composite %>%
  mutate(
    StudyTimeMatch = case_when(
      `Study Time Collected` == 0 ~ 0,
      `Study Time Collected` %in% c(28) ~ 1,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(StudyTimeMatch))

# Merge dataframes
eligible.participants <- nab_final %>%
  inner_join(
    rawdata %>% dplyr::select("participant_id", "study_accession", "uid", "expression_time"),
    by = c("Participant ID" = "participant_id", "StudyTimeMatch" = "expression_time")
  )

# Create a final data frame with the "treatment status"
filtered.assay.data <- eligible.participants %>%
  dplyr::select(`Participant ID`,
                uid,
                `Study Time Collected`,
                response,
                study_accession) %>%
  rename(`Antibody Value` = response) %>%
  mutate(Treatment = ifelse(`Study Time Collected` == 0, 0, 1))

# Change UID from columns to rows
ge.data.transposed <- ge.data %>%
  rownames_to_column(var = "uid") %>%  # Make row names a column
  pivot_longer(-uid, names_to = "gene", values_to = "expression") %>% # Convert wide to long
  pivot_wider(names_from = uid, values_from = expression) %>% # Convert back to wide with uid as row names
  rename(uid = gene)
rm(ge.data)

# Merge gene expression data and assay data
processed.data <- filtered.assay.data %>%
  left_join(ge.data.transposed, by = "uid")


## Pre-filtration with Wilcoxon-rank test

# Prepare a new data frame
df_expr <- processed.data
# Change column name
colnames(df_expr)[colnames(df_expr) == "Participant ID"] <- "ParticipantID"

# Identify metadata columns
meta_cols <- c(
  "ParticipantID",
  "uid",
  "Study Time Collected",
  "Antibody Value",
  "Treatment",
  "study_accession"
)

# Identify gene expression columns
gene_cols <- setdiff(colnames(df_expr), meta_cols)

# Initialize result table
wilcoxon.results <- data.frame(Gene = gene_cols, P_value = NA)

# Paired Wilcoxon test loop for each gene
for (g in gene_cols) {
  gene_data <- df_expr[, c("ParticipantID", "Treatment", g)]
  colnames(gene_data)[3] <- "GeneExpr"
  wide <- gene_data %>% pivot_wider(
    names_from = Treatment,
    values_from = GeneExpr,
    names_prefix = "T"
  )
  
  res <- wilcox.test(wide$T0, wide$T1, paired = TRUE, exact = FALSE)
  wilcoxon.results$P_value[wilcoxon.results$Gene == g] <- res$p.value
}


# Apply low-stringency threshold (e.g., p < 0.1)
filtered_genes <- wilcoxon.results[wilcoxon.results$P_value < 0.05, ]

# Remove unecessary vairables
rm(
  list = c(
    "df_expr",
    "meta_cols",
    "gene_cols",
    "wilcoxon.results",
    "gene_data",
    "g",
    "wide",
    "complete_cases",
    "res"
  )
)


# Prepare the data frame
Wilcoxon.filtered <- processed.data %>% filter(Treatment == 1) # filter for only vaccinated individuals

# Create a matrix with gene expression values
ge.values <- Wilcoxon.filtered %>%
  select(all_of(filtered_genes$Gene))  %>% #filter only for Wilcoxon genes
  as.matrix()

# Create a vector with antibody values
# antibody.values <- log(Wilcoxon.filtered$`Antibody Value` + 0.001) # create a vector for antibody values
antibody.values <- Wilcoxon.filtered$`Antibody Value`

# removing unnecessary  variables
rm(list = c("Wilcoxon.filtered"))


# Bootstrapping for variable selection
# Designate parameters
n_bootstraps <- 10
ncomp.grid <- 1:5
keepX.grid <- c(1, 2, 5, 10, 20, 50, 100)
n.inner.folds = 5

X <- ge.values
Y <- antibody.values

# Initiate dataframes to store results
selected.genes.list <- list()

# Set seed
set.seed(666)

for (i in 1:n_bootstraps) {
  cat("Bootstrap:", i, "\n") # progress message
  
  train.idx <- sample(1:nrow(X), replace = TRUE) # training ids
  test.idx <- setdiff(1:nrow(X), unique(train.idx)) # testing ids
  
  X.train <- X[train.idx, , drop = FALSE] # training data
  Y.train <- Y[train.idx]
  X.test <- X[test.idx, , drop = FALSE] # testing data
  Y.test <- Y[test.idx]
  
  best.cor <- -Inf # store best correlation value
  best.ncomp <- NULL # store best ncomp value
  best.keepX <- NULL # store best keepX
  
  for (ncomp in ncomp.grid) {
    # for each ncomp
    for (keepX in keepX.grid) {
      # for each keepX
      
      preds_all <- rep(NA_real_, length(train.idx)) # store predicttions
      inner.folds <- sample(rep(1:n.inner.folds, length.out = length(train.idx))) # create inner folds
      
      for (inner.fold in 1:n.inner.folds) {
        # inner loop for hyperparameter tuning
        inner.test.idx <- which(inner.folds == inner.fold) # inner testing ids
        inner.train.idx <- setdiff(1:length(train.idx), inner.test.idx) # inner training ids
        
        X.inner.train <- X.train[inner.train.idx, , drop = FALSE] # training data
        Y.inner.train <- Y.train[inner.train.idx]
        X.inner.test <- X.train[inner.test.idx, , drop = FALSE] # testing data
        
        model <- spls(
          X.inner.train,
          Y.inner.train,
          ncomp = ncomp,
          keepX = rep(keepX, ncomp)
        ) # build model with hyperparameter values
        pred <- predict(model, X.inner.test)$predict[, , ncomp] # predict on test data
        
        preds_all[inner.test.idx] <- as.vector(pred) # store predictions
      }
      
      pooled.cor <- cor(preds_all, Y.train, use = "complete.obs") # compute correlation as metric to select hyperparameters
      if (!is.na(pooled.cor) && pooled.cor > best.cor) {
        best.cor <- pooled.cor
        best.ncomp <- ncomp
        best.keepX <- keepX
      }
    }
  }
  
  final.model <- spls(X.train,
                      Y.train,
                      ncomp = best.ncomp,
                      keepX = rep(best.keepX, best.ncomp)) # build final model
  
  selected.variables <- unlist(lapply(1:best.ncomp, function(comp) {
    selectVar(final.model, comp = comp)$X$name
  })) # extract selected variables
  selected.genes.list[[i]] <- selected.variables # store selected variables
}

# Remove unnecessary variables
# rm(list = c("n_bootstraps", "ncomp.grid", "keepX.grid", "ge.values", "antibody.values", "i", "train.idx", "test.idx", "X.train", "Y.train", "X.test", "Y.test", "best.cor", "best.ncomp", "best.keepX", "ncomp", "keepX", "model", "pred", "corr", "final.model", "selected.variables"))

# Flatten all selected genes from all bootstrap iterations
all_selected_genes <- unlist(selected.genes.list)

# Count how often each gene was selected
gene_freq_table <- table(all_selected_genes)

# Convert to data frame and calculate proportions
gene_summary <- data.frame(
  gene = names(gene_freq_table),
  count = as.vector(gene_freq_table),
  proportion = as.vector(gene_freq_table) / length(selected.genes.list)
)

# Sort from most to least frequent
gene_summary <- gene_summary[order(-gene_summary$proportion), ]
gene_summary$rank <- seq_len(nrow(gene_summary))

# Select top 10 genes to show on graph
top10 <- head(gene_summary, 10)

# Plot
elbow.plot.predictors <- ggplot(gene_summary, aes(x = rank, y = proportion)) +
  geom_point(color = "darkblue") +
  geom_label_repel(
    data = top10,
    aes(label = gene),
    nudge_x = 20,
    direction = "y",
    hjust = 0,
    size = 3.2,
    box.padding = 0.35,
    label.padding = 0.15,
    label.size = 0.3,
    segment.color = "gray40",
    max.overlaps = Inf,
    fill = "white",
    color = "black"
  ) +
  geom_hline(yintercept = 0.22,
             linetype = "dashed",
             color = "red") +
  annotate(
    "text",
    x = Inf,
    y = 0.22,
    label = "Threshold: > 0.22",
    hjust = 1.05,
    vjust = -0.5,
    color = "black",
    size = 3.5
  ) +
  scale_x_continuous(breaks = seq(0, max(gene_summary$rank), by = 500)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  labs(x = "Gene Rank", y = "Proportion of Bootstraps Selected") +
  coord_cartesian(clip = "off")

print(elbow.plot.predictors)

predictor.genes <- gene_summary$gene[gene_summary$proportion >= 0.22]
# Parameters
n.outer.folds <- 10
n.inner.folds <- 10
ncomp.grid <- 1:5
keepX.grid <- c(1, 2, 5, 10, 20, 50)

# Delimiting x and y
X.filtered <- X[, predictor.genes]

# Store results
all.actual <- c()
all.predicted <- c()

# Set seed
set.seed(666)

# Designate outer folds
outer.folds <- sample(rep(1:n.outer.folds, length.out = nrow(X.filtered)))

# Cross-validation
for (outer.fold in 1:n.outer.folds) {
  cat("Outer Fold:", outer.fold, "\n")
  
  test.idx <- which(outer.folds == outer.fold)
  train.idx <- setdiff(1:nrow(X.filtered), test.idx)
  
  X.train <- X.filtered[train.idx, ]
  Y.train <- Y[train.idx]
  X.test <- X.filtered[test.idx, ]
  Y.test <- Y[test.idx]
  
  # Inner CV
  best.cor <- -Inf
  best.ncomp <- NULL
  best.keepX <- NULL
  
  for (ncomp in ncomp.grid) {
    for (keepX in keepX.grid) {
      preds_all <- rep(NA_real_, length(train.idx))
      inner.folds <- sample(rep(1:n.inner.folds, length.out = length(train.idx)))
      
      for (inner.fold in 1:n.inner.folds) {
        inner.test.idx <- which(inner.folds == inner.fold)
        inner.train.idx <- setdiff(1:length(train.idx), inner.test.idx)
        
        X.inner.train <- X.train[inner.train.idx, ]
        Y.inner.train <- Y.train[inner.train.idx]
        X.inner.test <- X.train[inner.test.idx, ]
        Y.inner.test <- Y.train[inner.test.idx]
        
        model <- spls(
          X.inner.train,
          Y.inner.train,
          ncomp = ncomp,
          keepX = rep(keepX, ncomp)
        )
        pred <- predict(model, X.inner.test)$predict[, , ncomp]
        
        preds_all[inner.test.idx] <- as.vector(pred)
      }
      
      pooled.cor <- cor(preds_all, Y.train, use = "complete.obs")
      
      if (!is.na(pooled.cor) && pooled.cor > best.cor) {
        best.cor <- pooled.cor
        best.ncomp <- ncomp
        best.keepX <- keepX
      }
    }
  }
  
  # Train final model
  final.model <- spls(X.train,
                      Y.train,
                      ncomp = best.ncomp,
                      keepX = rep(best.keepX, best.ncomp))
  final.pred <- predict(final.model, X.test)$predict[, , best.ncomp]
  
  all.actual <- c(all.actual, Y.test)
  all.predicted <- c(all.predicted, final.pred)
}

# Data frame of results
predictors.cv.plot.data <- data.frame(Actual = all.actual, Predicted = all.predicted)

# Determine limits from the combined actual + predicted range
axis.limits <- range(c(
  predictors.cv.plot.data$Actual,
  predictors.cv.plot.data$Predicted
),
na.rm = TRUE)

# Plot
predictors.cv.plot <- ggplot(predictors.cv.plot.data, aes(x = Actual, y = Predicted)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "red",
    size = 1
  ) +
  coord_fixed(
    ratio = 1,
    xlim = axis.limits,
    ylim = axis.limits,
    clip = "off"
  ) +
  theme_minimal(base_size = 18) +
  labs(
    #title = "Predicted vs Actual Values (sPLS)",
    subtitle = paste("Pearson Correlation:", round(
      cor(all.actual, all.predicted, method = "pearson"), 3
    )),
    x = expression("Actual " * log(Y^d28)),
    y = expression("Predicted " * log(Y^d28))
  )  +
  theme(axis.title = element_text(size = 20))

print(predictors.cv.plot)
