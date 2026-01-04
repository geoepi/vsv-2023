plot_migration_network <- function(tree,
                                   vsv_meta,
                                   ca_counties,
                                   epi_start,
                                   epi_end) {
  
  # tree data
  df <- as_tibble(tree) %>% as.data.frame()
  edges <- tree@phylo$edge
  df$node_num <- seq_len(nrow(df))
  
  migrations <- data.frame(
    from   = df$location[edges[, 1]],
    to     = df$location[edges[, 2]],
    parent = edges[, 1],
    child  = edges[, 2]
  ) %>%
    filter(from != to, !is.na(from), !is.na(to))
  
  # match coords
  locs_match <- vsv_meta %>%
    mutate(location = gsub(" ", "", county)) %>%
    select(location, x, y) %>%
    distinct()
  
  tmp_from <- migrations %>%
    select(from) %>%
    left_join(locs_match, by = c("from" = "location")) %>%
    dplyr::rename(x_from = x, y_from = y)
  
  tmp_to <- migrations %>%
    select(to) %>%
    left_join(locs_match, by = c("to" = "location")) %>%
    dplyr::rename(x_to = x, y_to = y)
  
  migrations_plot <- migrations %>%
    bind_cols(tmp_from[, -1], tmp_to[, -1])
  
  # heights to days
  latest_sample_date <- epi_end
  
  latest_decimal <- as.numeric(format(latest_sample_date, "%Y")) +
    (as.numeric(latest_sample_date - as.Date(paste0(format(latest_sample_date, "%Y"), "-01-01"))) / 365.25)
 
  migrations_plot <- migrations_plot %>%
    left_join(df %>% select(node, height), by = c("child" = "node")) %>%
    mutate(
      migration_time_days = as.numeric(height) * 365.25
    )
  
  # labels
  locs_labels <- vsv_meta %>%
    mutate(location = gsub(" ", "", county)) %>%
    select(location, county, x, y) %>%
    distinct()
  
  # map units
  us_states <- tigris::states(class = "sf", year = 2022)
  western_states <- c(
    "California", "Oregon", "Washington", "Nevada", "Idaho", "Texas",
    "Utah", "Arizona", "Colorado", "New Mexico", "Montana", "Wyoming",
    "Oklahoma", "Kansas", "Nebraska"
  )
  us_west_sf <- us_states %>% filter(NAME %in% western_states)
  
  mex_states <- ne_states(country = "Mexico", returnclass = "sf")
  
  
  p <- ggplot() +
    geom_sf(data = mex_states, fill = "gray95", color = "gray70") +
    geom_sf(data = us_west_sf, fill = "gray97", color = "gray40", linewidth = 0.25) +
    geom_sf(data = ca_counties, fill = "gray97", color = "gray40", linewidth = 0.25) +
    geom_point(
      data = locs_match,
      aes(x = x, y = y),
      size = 5,
      color = "tan"
    ) +
    geom_segment(
      data = migrations_plot,
      aes(
        x = x_from, y = y_from,
        xend = x_to, yend = y_to,
        color = migration_time_days
      ),
      arrow = arrow(length = unit(0.30, "cm"), type = "closed"),
      alpha = 0.5,
      linewidth = 1
    ) +
    scale_color_viridis_c(
      option = "D",
      direction = -1,
      name = "Days"
    ) +
    ggrepel::geom_text_repel(
      data = locs_labels,
      aes(x = x, y = y, label = county),
      color = "gray30",
      size = 3.5,
      box.padding = 0.3,
      point.padding = 0.2,
      segment.color = "gray50",
      segment.size = 0.3,
      max.overlaps = Inf
    )+
    coord_sf(
      xlim = c(-125, -89),
      ylim = c(26, 43),
      expand = FALSE
    ) +
    ggspatial::annotation_scale(
      location = "bl",
      width_hint = 0.25,
      text_cex = 0.9,
      pad_y = unit(0.6, "cm"),
      line_width = 0.2
    ) +
    ggspatial::annotation_north_arrow(
      location = "bl",
      pad_x = unit(3, "cm"),
      pad_y = unit(0.8, "cm"),
      which_north = "true",
      style = ggspatial::north_arrow_fancy_orienteering(
        text_size = 15,
        line_width = 1
      )
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = c(0.05, 0.3),
      legend.position.inside = TRUE,
      legend.title = element_text(size = 16, color="gray30"),
      legend.text = element_text(size = 10, color="gray30"),
      legend.key.width = unit(1, "line"),
      legend.key.height = unit(2, "line"),
      axis.title.x = element_text(size = 24, color="gray30"),
      axis.title.y = element_text(size = 24, color = "gray30"),
      axis.text.x = element_text(size = 10, face = "bold"),
      axis.text.y = element_text(size = 10, face = "bold"),
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5)
    ) +
    labs(
      x = "Longitude",
      y = "Latitude",
      title = " "
    )
  
  p
}
