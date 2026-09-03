library(readxl)
library(tidyverse)
library(ggbreak)

file_GP <- "GP relative abundance.xlsx"
file_RP <- "RP relative abundance.xlsx"
file_MD <- "MD relative abundance.xlsx"
file_RB <- "RB relative abundance.xlsx"

output_dir <- "."

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

ASV_SPACING <- 0.88

GROUP_OFFSET <- 0.17

BOX_WIDTH <- 0.30

LEFT_PADDING <- 0.10

RIGHT_PADDING <- 0.18

BOX_LINEWIDTH <- 0.35

OUTPUT_WIDTH <- 8
OUTPUT_HEIGHT <- 4.7

read_asv_data <- function(file_path,
                          oral_group,
                          fecal_group,
                          asv_order = NULL) {

  raw <- read_excel(
    file_path,
    col_names = FALSE
  )

  sample_groups <- as.character(
    unlist(raw[1, -1])
  )

  asv_names <- as.character(
    unlist(raw[-1, 1])
  )

  abundance_matrix <- raw[-1, -1]

  abundance_matrix[] <- lapply(
    abundance_matrix,
    as.numeric
  )

  long_list <- vector(
    "list",
    ncol(abundance_matrix)
  )

  for (i in seq_len(ncol(abundance_matrix))) {

    long_list[[i]] <- data.frame(

      ASV = asv_names,

      Group = sample_groups[i],

      Relative_abundance =
        as.numeric(abundance_matrix[[i]]),

      Sample = paste0(
        "Sample_",
        i
      ),

      stringsAsFactors = FALSE
    )
  }

  df <- bind_rows(long_list)

  df <- df %>%
    filter(
      Group %in% c(
        oral_group,
        fecal_group
      ),
      !is.na(Relative_abundance),
      !is.na(ASV)
    )

  if (is.null(asv_order)) {
    asv_order <- asv_names
  }

  df$ASV <- factor(
    df$ASV,
    levels = asv_order
  )

  df$Group <- factor(
    df$Group,
    levels = c(
      oral_group,
      fecal_group
    )
  )

  return(df)
}

draw_boxplot <- function(df,
                         species,
                         oral_group,
                         fecal_group,
                         output_dir,
                         asv_spacing = ASV_SPACING,
                         group_offset = GROUP_OFFSET,
                         box_width = BOX_WIDTH,
                         left_padding = LEFT_PADDING,
                         right_padding = RIGHT_PADDING,
                         box_linewidth = BOX_LINEWIDTH) {

  asv_levels <- levels(df$ASV)

  n_asv <- length(asv_levels)

  asv_centers <- 1 +
    (seq_len(n_asv) - 1) * asv_spacing

  names(asv_centers) <- asv_levels

  df <- df %>%
    mutate(

      ASV_center =
        asv_centers[
          as.character(ASV)
        ],

      X_position = case_when(

        Group == oral_group ~
          ASV_center - group_offset,

        Group == fecal_group ~
          ASV_center + group_offset,

        TRUE ~ ASV_center
      )
    )

  data_max <- max(
    df$Relative_abundance,
    na.rm = TRUE
  )

  y_axis_break <- switch(
    species,
    "GP" = c(40, 60),
    "RP" = c(60, 80),
    "MD" = c(5, 20),
    "RB" = c(5, 20),
    NULL
  )

  lower_y_breaks <- switch(
    species,
    "GP" = c(0, 20, 40),
    "RP" = c(0, 20, 40, 60),
    "MD" = 0:5,
    "RB" = 0:5,
    NULL
  )

  upper_y_breaks <- switch(
    species,
    "GP" = c(60, 80),
    "RP" = c(80, 100),
    "MD" = c(20, 30),
    "RB" = c(20, 30),
    NULL
  )

  y_axis_ticks <- sort(unique(c(
    lower_y_breaks,
    upper_y_breaks
  )))

  break_endpoints <- y_axis_break

  manual_tick_values <- y_axis_ticks[
    !(y_axis_ticks %in% break_endpoints)
  ]

  manual_tick_length <- 0.045

  y_break_scale <- if (species %in% c("MD", "RB")) {
    0.35
  } else {
    "fixed"
  }

  y_max <- max(
    data_max * 1.02,
    max(upper_y_breaks) * 1.005,
    1
  )

  x_axis_gap <- y_axis_break[1] * 0.03
  x_axis_y <- -x_axis_gap

  break_low  <- y_axis_break[1]
  break_high <- y_axis_break[2]

  box_stats <- df %>%
    group_by(ASV, Group, X_position) %>%
    summarise(
      n = sum(!is.na(Relative_abundance)),
      q1 = quantile(Relative_abundance, 0.25, na.rm = TRUE, names = FALSE, type = 7),
      median = median(Relative_abundance, na.rm = TRUE),
      q3 = quantile(Relative_abundance, 0.75, na.rm = TRUE, names = FALSE, type = 7),
      iqr = IQR(Relative_abundance, na.rm = TRUE, type = 7),
      values = list(Relative_abundance[!is.na(Relative_abundance)]),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      lower_fence = q1 - 1.5 * iqr,
      upper_fence = q3 + 1.5 * iqr,
      ymin = min(values[values >= lower_fence], na.rm = TRUE),
      ymax = max(values[values <= upper_fence], na.rm = TRUE),
      outliers = list(values[values < ymin | values > ymax]),
      box_xmin = X_position - box_width / 2,
      box_xmax = X_position + box_width / 2,

      box_hits_break =
        pmin(q1, q3) < break_high &
        pmax(q1, q3) > break_low,

      lower_whisker_hits_break =
        pmin(ymin, q1) < break_high &
        pmax(ymin, q1) > break_low,

      upper_whisker_hits_break =
        pmin(q3, ymax) < break_high &
        pmax(q3, ymax) > break_low,

      affected_by_break =
        box_hits_break |
        lower_whisker_hits_break |
        upper_whisker_hits_break
    ) %>%
    ungroup()

  unaffected_keys <- box_stats %>%
    filter(!affected_by_break) %>%
    select(ASV, Group)

  df_unaffected <- df %>%
    semi_join(
      unaffected_keys,
      by = c("ASV", "Group")
    )

  affected_stats <- box_stats %>%
    filter(affected_by_break)

  box_fill_lower <- affected_stats %>%
    filter(q1 < break_low) %>%
    transmute(
      ASV, Group,
      xmin = box_xmin, xmax = box_xmax,
      ymin = q1, ymax = pmin(q3, break_low)
    ) %>%
    filter(ymax > ymin)

  box_fill_upper <- affected_stats %>%
    filter(q3 > break_high) %>%
    transmute(
      ASV, Group,
      xmin = box_xmin, xmax = box_xmax,
      ymin = pmax(q1, break_high), ymax = q3
    ) %>%
    filter(ymax > ymin)

  box_fill_parts <- bind_rows(
    box_fill_lower,
    box_fill_upper
  )

  split_only_crossing_segments <- function(seg_df, low, high) {

    if (nrow(seg_df) == 0) {
      return(seg_df)
    }

    pieces <- lapply(seq_len(nrow(seg_df)), function(i) {

      one <- seg_df[i, , drop = FALSE]

      y_start <- one$y
      y_end   <- one$yend

      y1 <- min(y_start, y_end)
      y2 <- max(y_start, y_end)

      hits_break <- (y1 < high) && (y2 > low)

      if (!hits_break) {
        return(one)
      }

      out <- list()

      if (y1 < low) {
        low_piece <- one
        low_piece$y <- y1
        low_piece$yend <- min(y2, low)

        if (low_piece$yend > low_piece$y) {
          out[[length(out) + 1]] <- low_piece
        }
      }

      if (y2 > high) {
        high_piece <- one
        high_piece$y <- max(y1, high)
        high_piece$yend <- y2

        if (high_piece$yend > high_piece$y) {
          out[[length(out) + 1]] <- high_piece
        }
      }

      if (length(out) == 0) {
        one[0, , drop = FALSE]
      } else {
        bind_rows(out)
      }
    })

    result <- bind_rows(pieces)

    if (nrow(result) == 0) {
      return(seg_df[0, , drop = FALSE])
    }

    result
  }

  raw_box_side_segments <- bind_rows(
    affected_stats %>%
      transmute(ASV, Group, x = box_xmin, y = q1, yend = q3),
    affected_stats %>%
      transmute(ASV, Group, x = box_xmax, y = q1, yend = q3)
  )

  box_vertical_segments <- split_only_crossing_segments(
    raw_box_side_segments,
    break_low,
    break_high
  )

  raw_whisker_segments <- bind_rows(
    affected_stats %>%
      transmute(ASV, Group, x = X_position, y = ymin, yend = q1),
    affected_stats %>%
      transmute(ASV, Group, x = X_position, y = q3, yend = ymax)
  )

  box_whisker_segments <- split_only_crossing_segments(
    raw_whisker_segments,
    break_low,
    break_high
  )

  box_horizontal_edges <- bind_rows(
    affected_stats %>%
      transmute(ASV, Group, x = box_xmin, xend = box_xmax, y = q1),
    affected_stats %>%
      transmute(ASV, Group, x = box_xmin, xend = box_xmax, y = q3)
  ) %>%
    filter(y <= break_low | y >= break_high)

  box_median_segments <- affected_stats %>%
    transmute(ASV, Group, x = box_xmin, xend = box_xmax, y = median) %>%
    filter(y <= break_low | y >= break_high)

  if (nrow(affected_stats) == 0 ||
      sum(lengths(affected_stats$outliers)) == 0) {

    box_outliers <- tibble(
      ASV = factor(levels = levels(df$ASV)),
      Group = factor(levels = levels(df$Group)),
      X_position = numeric(0),
      y = numeric(0)
    )

  } else {

    box_outliers <- affected_stats %>%
      select(ASV, Group, X_position, outliers) %>%
      tidyr::unnest_longer(outliers, values_to = "y") %>%
      filter(!is.na(y))
  }

  first_center <- min(asv_centers)

  last_center <- max(asv_centers)

  x_left <- first_center -
    group_offset -
    box_width / 2 -
    left_padding

  x_right <- last_center +
    group_offset +
    box_width / 2 +
    right_padding

  break_slash_dx <- 0.055

  x_panel_left <- x_left - break_slash_dx - 0.005

  lower_range <- y_axis_break[1] - x_axis_y
  break_slash_dy_lower <- lower_range * 0.04

  upper_range <- y_max - y_axis_break[2]

  upper_height_ratio <- if (is.numeric(y_break_scale)) {
    y_break_scale
  } else {
    upper_range / lower_range
  }

  upper_units_per_y_ratio <-
    upper_height_ratio * lower_range / upper_range

  break_slash_dy_upper <-
    break_slash_dy_lower / upper_units_per_y_ratio

  break_slash_dy_upper <- min(
    break_slash_dy_upper,
    upper_range * 0.18
  )

  lower_axis_end <-
    y_axis_break[1] - break_slash_dy_lower / 2

  upper_axis_start <-
    y_axis_break[2] + break_slash_dy_upper / 2

  p <- ggplot(
    df,
    aes(
      x = X_position,
      y = Relative_abundance,
      fill = Group,
      group = interaction(
        ASV,
        Group
      )
    )
  ) +

  geom_segment(
    data = data.frame(
      y_tick = manual_tick_values
    ),
    aes(
      x = x_left,
      xend = x_left - manual_tick_length,
      y = y_tick,
      yend = y_tick
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.45,
    lineend = "butt",
    show.legend = FALSE
  ) +

  geom_segment(
    data = data.frame(
      x = c(x_left, x_left),
      xend = c(x_left, x_left),
      y = c(x_axis_y, upper_axis_start),
      yend = c(lower_axis_end, y_max)
    ),
    aes(
      x = x, xend = xend,
      y = y, yend = yend
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.5,
    lineend = "butt",
    show.legend = FALSE
  ) +

  geom_segment(
    data = data.frame(
      x = x_left,
      xend = x_right,
      y = x_axis_y,
      yend = x_axis_y
    ),
    aes(
      x = x, xend = xend,
      y = y, yend = yend
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.5,
    lineend = "butt",
    show.legend = FALSE
  ) +

  geom_segment(
    data = data.frame(
      x = c(
        x_left - break_slash_dx,
        x_left - break_slash_dx
      ),
      xend = c(
        x_left + break_slash_dx,
        x_left + break_slash_dx
      ),
      y = c(
        y_axis_break[1],
        y_axis_break[2] + break_slash_dy_upper
      ),
      yend = c(
        y_axis_break[1] - break_slash_dy_lower,
        y_axis_break[2]
      )
    ),
    aes(
      x = x, xend = xend,
      y = y, yend = yend
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.5,
    lineend = "butt",
    show.legend = FALSE
  ) +

  geom_boxplot(
    data = df_unaffected,
    width = box_width,
    colour = "black",
    linewidth = box_linewidth,
    alpha = 1,
    outlier.shape = 16,
    outlier.size = 1.3,
    outlier.colour = "black",
    show.legend = TRUE
  ) +

    geom_rect(
      data = box_fill_parts,
      aes(
        xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax,
        fill = Group
      ),
      inherit.aes = FALSE,
      colour = NA,
      alpha = 1,
      show.legend = TRUE
    ) +

    geom_segment(
      data = box_vertical_segments,
      aes(x = x, xend = x, y = y, yend = yend),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = box_linewidth,
      lineend = "butt",
      show.legend = FALSE
    ) +

    geom_segment(
      data = box_whisker_segments,
      aes(x = x, xend = x, y = y, yend = yend),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = box_linewidth,
      lineend = "butt",
      show.legend = FALSE
    ) +

    geom_segment(
      data = box_horizontal_edges,
      aes(x = x, xend = xend, y = y, yend = y),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = box_linewidth,
      lineend = "butt",
      show.legend = FALSE
    ) +

    geom_segment(
      data = box_median_segments,
      aes(x = x, xend = xend, y = y, yend = y),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = box_linewidth,
      lineend = "butt",
      show.legend = FALSE
    ) +

    geom_point(
      data = box_outliers,
      aes(x = X_position, y = y),
      inherit.aes = FALSE,
      shape = 16,
      size = 1.3,
      colour = "black",
      show.legend = FALSE
    ) +

  scale_fill_manual(

    values = setNames(

      c(
        "#C7E3F6",
        "#F1AB78"
      ),

      c(
        oral_group,
        fecal_group
      )
    ),

    breaks = c(
      oral_group,
      fecal_group
    )
  ) +

  scale_x_continuous(

    breaks = asv_centers,

    labels = asv_levels,

    limits = c(
      x_panel_left,
      x_right
    ),

    expand = c(
      0,
      0
    )
  ) +

  scale_y_continuous(

    limits = c(
      x_axis_y,
      y_max
    ),

    breaks = y_axis_ticks,

    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  ggbreak::scale_y_break(
    breaks = y_axis_break,

    scales = y_break_scale,

    space = 0.02,

    symbol = NULL
  ) +

  labs(

    x = NULL,

    y = "Relative abundance (%)",

    fill = NULL
  ) +

  theme_classic(
    base_size = 11
  ) +

    theme(

      plot.background = element_blank(),
      panel.background = element_blank(),

      panel.grid.major = element_blank(),

      panel.grid.minor = element_blank(),

      axis.text.x = element_text(

        colour = "black",

        size = 10,

        angle = 0,

        hjust = 0.5,

        vjust = 0.5,

        margin = margin(
          t = 5
        )
      ),

      axis.text.y = element_text(
        colour = "black",
        size = 10,
        margin = margin(r = -1)
      ),

      axis.title.y = element_text(

        colour = "black",

        size = 11,

        margin = margin(
          r = 7
        )
      ),

      axis.line.x = element_blank(),

      axis.line.y.left = element_blank(),

      axis.line.y.right = element_blank(),
      axis.text.y.right = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.title.y.right = element_blank(),

      axis.ticks.x = element_line(
        colour = "black",
        linewidth = 0.45
      ),

      axis.ticks.y.left = element_blank(),

      axis.ticks.length = unit(
        0.15,
        "cm"
      ),

      legend.position = "right",

      legend.title = element_blank(),

      legend.background = element_blank(),

      legend.box.background = element_blank(),

      legend.key = element_blank(),

      legend.text = element_text(

        colour = "black",

        size = 10
      ),

      legend.key.width = unit(
        0.48,
        "cm"
      ),

      legend.key.height = unit(
        0.48,
        "cm"
      ),

      legend.spacing.x = unit(
        0.08,
        "cm"
      ),

      legend.spacing.y = unit(
        0.08,
        "cm"
      ),

      plot.margin = margin(
        t = 8,
        r = 8,
        b = 6,
        l = 6
      )
    ) +

  guides(

    fill = guide_legend(

      override.aes = list(

        colour = NA,

        linewidth = 0,

        alpha = 1
      )
    )
  )

  print(p)

  ggsave(

    filename = file.path(
      output_dir,
      paste0(
        species,
        "_RA_boxplot.png"
      )
    ),

    plot = p,

    width = OUTPUT_WIDTH,

    height = OUTPUT_HEIGHT,

    units = "in",

    dpi = 600,

    bg = "white"
  )

  pdf_file <- file.path(
    output_dir,
    paste0(species, "_RA_boxplot.pdf")
  )

  single_page_grob <- grid::grid.grabExpr(
    print(p),
    wrap = TRUE,
    width = OUTPUT_WIDTH,
    height = OUTPUT_HEIGHT
  )

  grDevices::pdf(
    file = pdf_file,
    width = OUTPUT_WIDTH,
    height = OUTPUT_HEIGHT,
    onefile = TRUE,
    bg = "transparent",
    useDingbats = FALSE
  )

  grid::grid.newpage()

  grid::grid.rect(
    x = 0.5,
    y = 0.5,
    width = 1,
    height = 1,
    gp = grid::gpar(
      fill = "white",
      col = NA
    )
  )

  grid::grid.draw(single_page_grob)

  grDevices::dev.off()

  return(p)
}

GP_data <- read_asv_data(

  file_path = file_GP,

  oral_group = "GPO",

  fecal_group = "GPF",

  asv_order = c(
    "ASV1",
    "ASV2",
    "ASV4",
    "ASV50",
    "ASV415"
  )
)

GP_plot <- draw_boxplot(

  df = GP_data,

  species = "GP",

  oral_group = "GPO",

  fecal_group = "GPF",

  output_dir = output_dir
)

RP_data <- read_asv_data(

  file_path = file_RP,

  oral_group = "RPO",

  fecal_group = "RPF",

  asv_order = c(
    "ASV1",
    "ASV2",
    "ASV4",
    "ASV5",
    "ASV52"
  )
)

RP_plot <- draw_boxplot(

  df = RP_data,

  species = "RP",

  oral_group = "RPO",

  fecal_group = "RPF",

  output_dir = output_dir
)

MD_data <- read_asv_data(

  file_path = file_MD,

  oral_group = "MDO",

  fecal_group = "MDF",

  asv_order = c(
    "ASV44",
    "ASV51",
    "ASV81",
    "ASV84",
    "ASV90"
  )
)

MD_plot <- draw_boxplot(

  df = MD_data,

  species = "MD",

  oral_group = "MDO",

  fecal_group = "MDF",

  output_dir = output_dir
)

RB_data <- read_asv_data(

  file_path = file_RB,

  oral_group = "RBO",

  fecal_group = "RBF",

  asv_order = c(
    "ASV28",
    "ASV30",
    "ASV67",
    "ASV128",
    "ASV3"
  )
)

RB_plot <- draw_boxplot(

  df = RB_data,

  species = "RB",

  oral_group = "RBO",

  fecal_group = "RBF",

  output_dir = output_dir
)
