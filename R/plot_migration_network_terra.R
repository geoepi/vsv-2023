plot_migration_network_terra <- function(tree,
                                   vsv_meta,
                                   ca_counties,
                                   epi_start,
                                   epi_end,
                                   loc_cols) {
  
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
  
  bbox_pts <- with(locs_match, c(
    xmin = min(x, na.rm = TRUE),
    ymin = min(y, na.rm = TRUE),
    xmax = max(x, na.rm = TRUE),
    ymax = max(y, na.rm = TRUE)
  ))
  
  bbox_pts <- as.numeric(bbox_pts)
  
  pad <- 2.0  # degrees expand
  
  bbox_coords <- c(
    left   = bbox_pts[1] - pad,
    bottom = bbox_pts[2] - pad,
    right  = bbox_pts[3] + pad,
    top    = bbox_pts[4] + pad
  )
  
  
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
  
  background_map <- get_map(location = bbox_coords, 
                            source = "stadia", maptype = "stamen_terrain")
  
  terra_p <- ggmap(background_map, darken = c(0.3, "white")) +
    geom_point(
      data = locs_match,
      aes(x = x, y = y, color = location),
      size = 7
    ) +
    scale_color_manual(
      values = loc_cols,
      name   = "Location",
      guide = guide_legend(
        title.position = "top",
        title.hjust = 0.5,
        label.position = "bottom",
        direction = "horizontal"
      )
    ) +
    ggnewscale::new_scale_color() +
    geom_segment(
      data = migrations_plot,
      aes(
        x = x_from, y = y_from,
        xend = x_to, yend = y_to,
        color = migration_time_days
      ),
      arrow = arrow(length = unit(0.30, "cm"), type = "closed"),
      alpha = 0.8,
      linewidth = 1
    ) +
    scale_color_gradient(
      low  = "gray65",
      high = "gray15",
      name = "Days",
      guide = guide_colourbar(
        title.position = "top",
        title.hjust = 0.5,
        direction = "horizontal"
      )
    ) + 
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      #legend.position = c(0.05, 0.3),
      #legend.position.inside = TRUE,
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 16, color="gray30"),
      legend.text = element_text(size = 10, color="gray30"),
      legend.key.width = unit(2, "line"),
      legend.key.height = unit(1, "line"),
      axis.title.x = element_text(size = 24, color="gray30"),
      axis.title.y = element_text(size = 24, color = "gray30"),
      axis.text.x = element_text(size = 20, face = "bold"),
      axis.text.y = element_text(size = 20, face = "bold"),
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5)
    ) +
    labs(
      x = "Longitude",
      y = "Latitude",
      title = " "
    )
  
  terra_p <- terra_p + ggspatial::annotation_scale(location = "bl", 
                                       width_hint = 0.25, 
                                       text_cex = 0.9, 
                                       pad_y = unit(0.6, "cm"), 
                                       line_width = 0.2
  ) + 
    ggspatial::annotation_north_arrow(location = "bl", 
                                      pad_x = unit(1.9, "cm"),
                                      pad_y = unit(0.8, "cm"), 
                                      which_north = "true", 
                                      style = ggspatial::north_arrow_fancy_orienteering(
                                        text_size = 15, line_width = 1 )
                                      ) +
    coord_sf(crs = 4326)

  terra_p
}
