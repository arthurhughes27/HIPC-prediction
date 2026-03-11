for (assay_of_interest in c("nAb", "hai")) {
  for (gender_of_interest in c("none", "Male", "Female")) {
    for (model_of_interest in c("elasticnet")) {
      for (feat.eng.col in c("none", "z", "rank")) {
        for (metric_of_interest in c("R.spearman", "sRMSE")) {
          for (feature_selection_of_interest in c("none", "univariate", "variance")) {
            # R script to interpret results from the prediction task for HAI in SDY1276 in both males and females
            library(dplyr)
            library(ggplot2)
            library(tidyr)
            
            # Path to save the results
            results_folder = fs::path("output", "results")
            figures_folder = fs::path("output",
                                      "figures",
                                      "comparison",
                                      gender_of_interest)
            
            # Study of interest
            study_of_interest = "SDY1276"
            
            # Vaccine of interest
            vaccine_of_interest = "Influenza (IN)"
            
            # Dataset of interest
            dataset_of_interest = "all_noNorm"
            
            # Assay of interest
            # assay_of_interest = "hai"
            
            # Model of interest
            # model_of_interest = "elasticnet"
            
            # Gender of interest
            # gender_of_interest = "none"
            
            # Feature engineering column of interest
            # feat.eng.col = "none"
            
            # Metric of interest
            # metric_of_interest = "R.spearman"
            
            # Feature selection algorithm
            # feature_selection_of_interest = "univariate"
            
            file_name_male = paste0(
              "metrics_transformations_",
              gsub(
                "[[:space:]()]",
                "",
                tolower(vaccine_of_interest)
              ),
              "_",
              dataset_of_interest,
              "_",
              study_of_interest,
              "_",
              assay_of_interest,
              "_Male_",
              model_of_interest,
              ifelse(
                feature_selection_of_interest == "none",
                "",
                paste0("_", feature_selection_of_interest, "Filtration")
              ),
              ".rds"
            )
            
            file_name_female = paste0(
              "metrics_transformations_",
              gsub(
                "[[:space:]()]",
                "",
                tolower(vaccine_of_interest)
              ),
              "_",
              dataset_of_interest,
              "_",
              study_of_interest,
              "_",
              assay_of_interest,
              "_Female_",
              model_of_interest,
              ifelse(
                feature_selection_of_interest == "none",
                "",
                paste0("_", feature_selection_of_interest, "Filtration")
              ),
              ".rds"
            )
            
            file_name_none = paste0(
              "metrics_transformations_",
              gsub(
                "[[:space:]()]",
                "",
                tolower(vaccine_of_interest)
              ),
              "_",
              dataset_of_interest,
              "_",
              study_of_interest,
              "_",
              assay_of_interest,
              "_none_",
              model_of_interest,
              ifelse(
                feature_selection_of_interest == "none",
                "",
                paste0("_", feature_selection_of_interest, "Filtration")
              ),
              ".rds"
            )
            
            p_load_male <- fs::path(results_folder, file_name_male)
            p_load_female <- fs::path(results_folder, file_name_female)
            p_load_none <- fs::path(results_folder, file_name_none)
            
            # Read the results
            metrics_df_male = readRDS(p_load_male) %>%
              mutate(gender.select = "Male")
            metrics_df_female = readRDS(p_load_female) %>%
              mutate(gender.select = "Female")
            metrics_df_none = readRDS(p_load_none) %>%
              mutate(gender.select = "none")
            
            
            metrics_df = bind_rows(bind_rows(metrics_df_male, metrics_df_female),
                                   metrics_df_none) %>%
              filter(!(feature.engineering.row %in% c("cv", "none")),!(data.selection %in% c("d1fc", "d3fc")))
            
            # Generate a figure comparing a given metric across data selections for a given model, gender, and column pre-transformation
            
            # Tidy names for the plot
            model_tidy = if (model_of_interest == "elasticnet") {
              "Elastic Net"
            } else if (model_of_interest == "ranger") {
              "Random Forest"
            } else if (model_of_interest == "xgboost") {
              "XGboost"
            } else if (model_of_interest == "lm") {
              "Linear regression"
            } else if (model_of_interest == "spls") {
              "sPLS"
            }
            
            gender_tidy = if (gender_of_interest == "none") {
              "both genders"
            } else {
              (paste0(tolower(gender_of_interest), "s only"))
            }
            
            feature_engineering_col_tidy = if (feat.eng.col == "none") {
              "No gene pre-transformation"
            } else if (feat.eng.col == "z") {
              "Genes pre-transformed with z-score"
            } else if (feat.eng.col == "rank") {
              "Genes pre-transformed with rank"
            }
            
            metric_tidy = if (metric_of_interest == "R2") {
              "R²"
            } else if (metric_of_interest == "R.spearman") {
              "Spearman R"
            } else if (metric_of_interest == "sRMSE") {
              "sRMSE"
            } else if (metric_of_interest == "RMSE") {
              "RMSE"
            }
            
            featureSelection_tidy = if (feature_selection_of_interest == "none") {
              "No feature selection"
            } else if (feature_selection_of_interest == "univariate") {
              "Univariate feature selection"
            } else if (feature_selection_of_interest == "variance") {
              "Top 100 by variance feature selection"
            }
            
            vaccine_tidy = tolower(gsub("[^[:alnum:]]", "", vaccine_of_interest))
            
            metric_tidy <- gsub("\\.", "", metric_of_interest)
            
            levels_order <- c("d0",
                              "d1",
                              "d3",
                              "d0+d1",
                              "d0+d3",
                              "d1+d3",
                              "d0+d1+d3")
            row_order <- c(
              "mean",
              "median",
              "max",
              "pc1",
              "ssgsea",
              "iqr"
              # ,"none"
            )
            
            baseline_R2 <- metrics_df %>%
              filter(data.selection == "baseline",
                     gender.select == gender_of_interest) %>%
              pull(.data [[metric_of_interest]])
            
            df <- metrics_df %>%
              filter(
                data.selection %in% levels_order,
                model == model_of_interest,
                gender.select == gender_of_interest,
                feature.engineering.col == feat.eng.col
              ) %>%
              mutate(
                data.selection = factor(data.selection, levels = levels_order),
                feature.engineering.row = factor(feature.engineering.row, levels = row_order)
              )
            
            min_val = metrics_df %>%
              filter(gender.select == gender_of_interest) %>%
              pull(.data [[metric_of_interest]]) %>%
              min()
            
            y_lim <- c(ifelse(min_val < 0, min_val, 0), 1)
            
            
            title = paste0("Comparison of ",
                           metric_tidy,
                           " values across data selections")
            subtitle = paste0(
              "Data from ",
              gender_tidy,
              ", model = ",
              model_tidy,
              ", ",
              feature_engineering_col_tidy,
              ", ",
              featureSelection_tidy
            )
            
            p1 <- ggplot(
              df %>% filter(include.covariates == TRUE),
              aes(
                x = data.selection,
                y = .data[[metric_of_interest]],
                fill = feature.engineering.row
              )
            ) +
              geom_col(width = 0.6,
                       position = position_dodge2(width = 0.9, padding = 0.2)) +
              geom_hline(aes(yintercept = baseline_R2, linetype = "Baseline"),
                         color = "black") +
              scale_linetype_manual(name = "",
                                    values = c("Baseline" = "dashed")) +
              scale_fill_brewer(palette = "Set2", drop = FALSE) +
              coord_cartesian(ylim = y_lim) +
              labs(
                title = title,
                subtitle = subtitle,
                x = "Predictor Set",
                y = metric_tidy,
                fill = "Row transformation"
              ) +
              theme_minimal(base_size = 20) +
              theme(
                plot.title = element_text(size = 25, hjust = 0.5),
                plot.subtitle = element_text(size = 17, hjust = 0.5)
              )
            
            p1
            
            # Path to save the results
            file_name = paste0(
              "comparison_",
              dataset_of_interest,
              "_",
              study_of_interest,
              "_",
              vaccine_tidy,
              "_",
              assay_of_interest,
              "_",
              gender_of_interest,
              "_",
              feat.eng.col,
              "_",
              model_of_interest,
              "_",
              metric_of_interest,
              ifelse(
                feature_selection_of_interest == "none",
                "",
                paste0("_", feature_selection_of_interest, "Filtration")
              ),
              ".pdf"
            )
            
            ggsave(
              filename = file_name,
              path = figures_folder,
              plot = p1,
              width = 40,
              height = 15,
              units = "cm"
            )
            
          }
        }
      }
    }
  }
}
