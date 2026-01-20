# Function to produce cross-validation predicted vs observed plot

cv.plot = function(pred, 
                   obs){
  
  predictions = data.frame(obs = obs, pred = pred)
  
  R2 = cor(obs, pred)^2
  sRMSE = rmse(obs, pred)/sd(obs)
  
  # plotting limits (same for x and y), add tiny padding
  min_all <- min(c(obs, pred), na.rm = TRUE)
  max_all <- max(c(obs, pred), na.rm = TRUE)
  span <- max_all - min_all
  pad <- if (span == 0) 1e-6 else 0.02 * span
  lims <- c(min_all - pad, max_all + pad)
  
  # annotation text
  anno <- paste0("R² = ", formatC(R2, format = "f", digits = 3),
                 "\nsRMSE = ", formatC(sRMSE, format = "f", digits = 3))
  
  # location for annotation: top-left inside plot
  x_anno <- lims[1] + 0.02 * (lims[2] - lims[1])
  y_anno <- lims[2] - 0.02 * (lims[2] - lims[1])
  
  prediction.plot = ggplot(predictions, aes(x = obs, y = pred)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
    coord_fixed(ratio = 1, xlim = lims, ylim = lims, expand = FALSE) +
    annotate("text", x = x_anno, y = y_anno, label = anno, hjust = 0, vjust = 1,
             size = 5) +
    labs(title = "Cross-validation: predicted vs observed",
         x = "Observed",
         y = "Predicted") +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 18),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14)
    )
  
  return(prediction.plot)
}