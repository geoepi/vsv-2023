phylodynamic_process <- function(tree,
                                 mrsd = as.Date("2023-11-08"),
                                 x_min = as.Date("2023-01-01"),
                                 x_max = as.Date("2023-11-20"),
                                 ribbon_color = "steelblue",
                                 alpha_outer = 0.25,
                                 alpha_inner = 0.45) {
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(scales)
    library(lubridate)
  })
  
  coal_pref <- BNPR2(tree, 
                     lengthout = 500, 
                     prec_alpha = 0.001, 
                     prec_beta = 0.001,
                     beta1_prec = 0.0001, 
                     fns = NULL, 
                     log_fns = FALSE, 
                     simplify = TRUE,
                     derivative = FALSE, 
                     forward = TRUE)
  
  coal_pref_df <- as.data.frame(
    cbind(
      date = coal_pref$x,
      Ne = coal_pref$effpop,
      Ne.0.025 = coal_pref$effpop025,
      Ne.0.25 = coal_pref$effpop25,
      Ne.0.75 = coal_pref$effpop75,
      Ne.0.975 = coal_pref$effpop975))
  
  coal_pref_df <- coal_pref_df %>%
    mutate(
      date = mrsd - days(round(date * 365.25, 0))
    ) %>%
    arrange(desc(date))
  
  ymin <- floor(log10(min(coal_pref_df$Ne.0.025, na.rm = TRUE)))
  ymax <- ceiling(log10(max(coal_pref_df$Ne.0.975, na.rm = TRUE)))
  log_breaks <- 10^(ymin:ymax)
  
  ggplot(coal_pref_df, aes(date, Ne)) +
    geom_ribbon(
      aes(ymin = Ne.0.025, ymax = Ne.0.975),
      fill = ribbon_color,
      alpha = alpha_outer
    ) +
    geom_ribbon(
      aes(ymin = Ne.0.25, ymax = Ne.0.75),
      fill = ribbon_color,
      alpha = alpha_inner
    ) +
    geom_line(color = "black", linewidth = 1) +
    scale_x_date(
      date_breaks = "60 days",
      date_labels = "%b %Y",
      limits = c(x_min, x_max)
    ) +
    scale_y_continuous(
      trans = "log10",
      breaks = log_breaks,
      labels = trans_format("log10", math_format(10^.x))
    ) +
    ylab("Effective Population Size (Ne)") +
    xlab(" ") +
    theme_minimal() +
    theme(
      plot.margin = unit(c(2, 0.5, 2, 0.5), "cm"),
      legend.direction = "vertical",
      legend.position = c(0.9, 0.8),
      strip.text = element_text(size = 26, face = "bold"),
      strip.background = element_blank(),
      legend.key.size = unit(2, "line"),
      legend.key.width = unit(1, "line"),
      legend.text = element_text(size = 16, face = "bold"),
      legend.title = element_text(size = 18, face = "bold"),
      axis.title.x = element_text(size = 20, face = "bold", vjust = -2),
      axis.title.y = element_text(size = 20, face = "bold", vjust = 3),
      axis.text.x = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 28, face = "bold")
    )
}
