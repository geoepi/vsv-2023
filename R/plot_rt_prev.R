plot_rt_prev <- function(vsv_simulation,
                             death_rate,
                             epi_start,
                             epi_end,
                             burn_in = 2000,
                             alpha_outer = 0.25,
                             alpha_inner = 0.45) {
  
  stopifnot(
    all(c("birth_rate", "prevalence") %in% names(vsv_simulation))
  )
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(matrixStats)
  
  epi_start <- as.Date(epi_start)
  epi_end   <- as.Date(epi_end)
  
  summarize_mcmc <- function(mat, scale = 1) {
    mat_post <- mat[-(1:burn_in), , drop = FALSE]
    
    n_time <- ncol(mat_post)
    dates  <- seq(epi_start, epi_end, length.out = n_time)
    
    tibble(
      date  = dates,
      mean  = colMeans(mat_post, na.rm = TRUE) / scale,
      q025  = colQuantiles(mat_post, probs = 0.025, na.rm = TRUE) / scale,
      q25   = colQuantiles(mat_post, probs = 0.25,  na.rm = TRUE) / scale,
      q75   = colQuantiles(mat_post, probs = 0.75,  na.rm = TRUE) / scale,
      q975  = colQuantiles(mat_post, probs = 0.975, na.rm = TRUE) / scale
    )
  }
  
  rt_df <- summarize_mcmc(vsv_simulation$birth_rate, scale = death_rate) %>%
    mutate(panel = "Reproduction number")
  
  prev_df <- summarize_mcmc(vsv_simulation$prevalence, scale = 1) %>%
    mutate(panel = "Prevalence")
  
  plot_df <- bind_rows(rt_df, prev_df)
  
  ggplot(plot_df, aes(x = date)) +
    geom_ribbon(
      aes(ymin = q025, ymax = q975),
      fill = "steelblue",
      alpha = alpha_outer
    ) +
    geom_ribbon(
      aes(ymin = q25, ymax = q75),
      fill = "steelblue",
      alpha = alpha_inner
    ) +
    geom_line(aes(y = mean), linewidth = 0.9) +
    geom_hline(
      data = tibble(panel = "Reproduction number"),
      aes(yintercept = 1),
      color = "gray40",
      linetype = "dashed",
      linewidth = 0.6
    ) +
    facet_wrap(~ panel, ncol = 1, scales = "free_y") +
    scale_x_date(
      limits = c(epi_start, epi_end),
      expand = c(0, 0)
    ) +
    labs(
      x = " ",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      strip.text     = element_text(size = 18, face = "bold", color = "gray40"),
      axis.title.x   = element_text(size = 20, face = "bold"),
      axis.title.y   = element_text(size = 20, face = "bold"),
      axis.text.x    = element_text(size = 10, face = "bold"),
      axis.text.y    = element_text(size = 10, face = "bold"),
    )
}

