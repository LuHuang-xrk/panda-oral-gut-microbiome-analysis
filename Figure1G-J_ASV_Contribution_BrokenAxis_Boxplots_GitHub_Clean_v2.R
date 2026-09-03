required_packages <- c(
  "readxl",
  "ggplot2",
  "dplyr",
  "tidyr"
)

for (pkg in required_packages) {

  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }

}

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)

task_list <- list(

  list(
    file_in = "GPO contribution.xlsx",
    dir_out = ".",
    name = "GP",
    ylab = "GPO contribution (%)",
    break_range = c(30, 70),
    lower_ticks = c(0, 10, 20, 30),
    upper_ticks = c(70, 80)
  ),

  list(
    file_in = "RPO contribution.xlsx",
    dir_out = ".",
    name = "RP",
    ylab = "RPO contribution (%)",
    break_range = c(60, 80),
    lower_ticks = c(0, 20, 40, 60),
    upper_ticks = c(80, 90, 100)
  ),

  list(
    file_in = "MDO contribution.xlsx",
    dir_out = ".",
    name = "MD",
    ylab = "MDO contribution (%)",
    break_range = c(3, 15),
    lower_ticks = c(0, 1, 2, 3),
    upper_ticks = c(16, 20)
  ),

  list(
    file_in = "RBO contribution.xlsx",
    dir_out = ".",
    name = "RB",
    ylab = "RBO contribution (%)",
    break_range = c(3, 20),
    lower_ticks = c(0, 1, 2, 3),
    upper_ticks = c(20, 24)
  )

)

main_output_dir <- "."

if (!dir.exists(main_output_dir)) {

  dir.create(
    main_output_dir,
    recursive = TRUE
  )

}

global_palette_27 <- c(

  "#C7E3F6",
  "#CBCFE9",
  "#F3A8A8",
  "#C884B6",
  "#3D73B9",
  "#F1AB78",
  "#80C7A7",
  "#E1C4A6",
  "#BCCEA8",
  "#5EB0A2",
  "#8FA9D1",
  "#A9B6D8",
  "#D6A6C6",
  "#B58AAE",
  "#8E6FA8",
  "#6F86B7",
  "#79A6C9",
  "#74B7B2",
  "#93C6B3",
  "#A8C08D",
  "#C7B97F",
  "#D8B07C",
  "#D69A7E",
  "#C98F8F",
  "#B89A88",
  "#8FA39B",
  "#9C9FB7"

)

priority_taxon_colors <- c(

  "Escherichia–Shigella" = "#C7E3F6",

  "Streptococcus" = "#CBCFE9",

  "Leuconostoc" = "#F3A8A8",

  "f_Muribaculaceae" = "#C884B6",

  "Bacteroides" = "#3D73B9"

)

remaining_palette_22 <- c(

  "#F1AB78",
  "#74B7B2",
  "#D6A6C6",
  "#A8C08D",
  "#79A6C9",
  "#D69A7E",
  "#93C6B3",
  "#8E6FA8",
  "#D8B07C",
  "#8FA39B",
  "#8FA9D1",
  "#C98F8F",
  "#BCCEA8",
  "#B58AAE",
  "#E1C4A6",
  "#5EB0A2",
  "#A9B6D8",
  "#C7B97F",
  "#80C7A7",
  "#B89A88",
  "#6F86B7",
  "#9C9FB7"

)

normalize_taxon <- function(x) {

  x <- trimws(
    as.character(x)
  )

  x <- gsub(
    "\\s+",
    " ",
    x
  )

  x <- gsub(
    "Escherichia[-–—−]Shigella",
    "Escherichia–Shigella",
    x
  )

  x <- gsub(
    "^f__Muribaculaceae$",
    "f_Muribaculaceae",
    x
  )

  return(x)

}

read_wide_matrix_to_long <- function(file_path) {

  if (!file.exists(file_path)) {

    stop(
      sprintf(
        "文件不存在！路径：%s",
        file_path
      )
    )

  }

  raw <- read_excel(
    file_path,
    col_names = FALSE
  )

  n_row <- nrow(raw)

  cat(
    sprintf(
      "读取原始表格，总行数：%d\n",
      n_row
    )
  )

  taxon_row <- raw[
    n_row,
    -1
  ]

  asv_header <- raw[
    1,
    -1
  ]

  sample_data <- raw[
    2:(n_row - 1),
  ]

  colnames(sample_data)[1] <- "SampleID"

  colnames(sample_data)[-1] <- as.character(
    asv_header
  )

  taxon_map <- tibble(

    ASV = as.character(
      asv_header
    ),

    Taxon = normalize_taxon(
      as.character(
        taxon_row
      )
    )

  )

  orig_asv_order <- as.character(
    asv_header
  )

  long_df <- sample_data %>%

    pivot_longer(

      cols = -SampleID,

      names_to = "ASV",

      values_to = "Value"

    ) %>%

    mutate(

      Value = suppressWarnings(
        as.numeric(Value)
      )

    ) %>%

    left_join(

      taxon_map,

      by = "ASV"

    ) %>%

    filter(

      !is.na(Value),

      !is.na(Taxon),

      Taxon != "",

      ASV != ""

    )

  cat(
    sprintf(
      "转换后长表行数：%d；样本数：%d；ASV数目：%d\n",

      nrow(long_df),

      n_distinct(
        long_df$SampleID
      ),

      n_distinct(
        long_df$ASV
      )
    )
  )

  return(

    list(

      df = long_df,

      asv_order = orig_asv_order,

      taxon_map = taxon_map

    )

  )

}

data_list <- list()

for (t in task_list) {

  cat(
    "\n=====准备读取：",
    t$file_in,
    "=====\n"
  )

  if (!file.exists(t$file_in)) {

    warning(
      paste(
        "跳过！找不到文件：",
        t$file_in
      )
    )

    next

  }

  res <- read_wide_matrix_to_long(
    t$file_in
  )

  data_list[[t$name]] <- res

}

if (length(data_list) == 0) {

  stop(
    "没有成功读取任何表格，请核对文件路径！"
  )

}

all_taxa <- unique(

  unlist(

    lapply(

      data_list,

      function(x) {
        x$df$Taxon
      }

    )

  )

)

cat(
  "\n全部Taxon数量：",
  length(all_taxa),
  "\n"
)

priority_present <- intersect(

  names(
    priority_taxon_colors
  ),

  all_taxa

)

priority_map <- priority_taxon_colors[
  priority_present
]

remaining_taxa <- setdiff(

  all_taxa,

  priority_present

)

unused_priority_colors <- setdiff(

  unname(
    priority_taxon_colors
  ),

  unname(
    priority_map
  )

)

available_colors <- unique(

  c(
    remaining_palette_22,
    unused_priority_colors
  )

)

available_colors <- available_colors[

  !available_colors %in%
    unname(priority_map)

]

if (
  length(remaining_taxa) >
  length(available_colors)
) {

  stop(
    paste0(
      "Taxon数量超过当前颜色数量！",
      "剩余Taxon：",
      length(remaining_taxa),
      "；剩余颜色：",
      length(available_colors)
    )
  )

}

remaining_map <- setNames(

  if (length(remaining_taxa) > 0) {

    available_colors[
      seq_along(
        remaining_taxa
      )
    ]

  } else {

    character(0)

  },

  remaining_taxa

)

global_taxon_colors <- c(

  priority_map,

  remaining_map

)

global_taxon_colors <- global_taxon_colors[
  all_taxa
]

if ("UCG-005" %in% names(global_taxon_colors)) {

  global_taxon_colors["UCG-005"] <- "#686C91"

}

color_map_df <- data.frame(

  Taxon = names(
    global_taxon_colors
  ),

  Color = unname(
    global_taxon_colors
  ),

  stringsAsFactors = FALSE

)

write.csv(

  color_map_df,

  file.path(

    main_output_dir,

    "global_taxon_color.csv"

  ),

  row.names = FALSE,

  fileEncoding = "UTF-8"

)

make_manual_y_break <- function(
    df,
    break_range
) {

  y <- df$Value
  y <- y[is.finite(y)]

  if (length(y) < 2) {
    return(NULL)
  }

  y_max <- max(y, na.rm = TRUE)

  if (
    length(break_range) != 2 ||
    any(!is.finite(break_range))
  ) {
    warning("break_range 必须是两个有效数值，例如 c(35, 70)。")
    return(NULL)
  }

  break_low <- min(break_range)
  break_high <- max(break_range)

  if (
    break_low <= 0 ||
    break_high <= break_low ||
    break_low >= y_max
  ) {
    warning(
      sprintf(
        "断轴区间 %.3f ~ %.3f 与当前数据范围不匹配，将保留连续Y轴。",
        break_low,
        break_high
      )
    )
    return(NULL)
  }

  if (break_high >= y_max) {
    warning(
      sprintf(
        paste0(
          "断轴上限 %.3f >= 数据最大值 %.3f，",
          "无法形成有效上段Y轴，将保留连续Y轴。"
        ),
        break_high,
        y_max
      )
    )
    return(NULL)
  }

  box_stats <- df %>%
    filter(is.finite(Value)) %>%
    group_by(ASV) %>%
    summarise(
      Q1 = as.numeric(
        quantile(
          Value,
          probs = 0.25,
          na.rm = TRUE,
          type = 7
        )
      ),
      Q3 = as.numeric(
        quantile(
          Value,
          probs = 0.75,
          na.rm = TRUE,
          type = 7
        )
      ),
      .groups = "drop"
    )

  cut_boxes <- box_stats %>%
    filter(
      Q1 < break_high,
      Q3 > break_low
    )

  if (nrow(cut_boxes) > 0) {
    warning(
      sprintf(
        paste0(
          "当前断轴 %.3f ~ %.3f 会穿过 %d 个箱体：%s。",
          "建议重新检查断轴区间。"
        ),
        break_low,
        break_high,
        nrow(cut_boxes),
        paste(as.character(cut_boxes$ASV), collapse = ", ")
      )
    )
  }

  n_cut_points <- sum(
    y > break_low & y < break_high,
    na.rm = TRUE
  )

  upper_scale <- 0.30

  list(
    lower = break_low,
    upper = break_high,
    y_max = y_max,
    n_cut_points = n_cut_points,
    cut_boxes = cut_boxes,
    upper_scale = upper_scale
  )
}

draw_boxplot <- function(
    df,
    orig_asv_order,
    taxon_map,
    dir_out,
    plotname,
    ylab_text,
    break_range,
    lower_ticks,
    upper_ticks,
    global_color_map
) {

  if (!dir.exists(dir_out)) {

    dir.create(
      dir_out,
      recursive = TRUE
    )

  }

  valid_asv <- intersect(

    orig_asv_order,

    unique(
      df$ASV
    )

  )

  df$ASV <- factor(

    df$ASV,

    levels = valid_asv

  )

  legend_order_df <- taxon_map %>%

    filter(
      ASV %in% valid_asv
    ) %>%

    mutate(

      ASV = factor(
        ASV,
        levels = valid_asv
      )

    ) %>%

    arrange(
      ASV
    )

  legend_taxa <- unique(
    legend_order_df$Taxon
  )

  legend_taxa <- legend_taxa[
    !is.na(legend_taxa) &
      legend_taxa != ""
  ]

  df$Taxon <- factor(

    df$Taxon,

    levels = legend_taxa

  )

  panel_pal <- global_color_map[
    legend_taxa
  ]

  cat(
    "\n[",
    plotname,
    "] ASV顺序：\n",
    paste(
      valid_asv,
      collapse = " -> "
    ),
    "\n",
    sep = ""
  )

  cat(
    "[",
    plotname,
    "] 图例Taxon顺序：\n",
    paste(
      legend_taxa,
      collapse = " -> "
    ),
    "\n",
    sep = ""
  )

  y_data_max <- max(
    df$Value,
    na.rm = TRUE
  )

  y_max_auto <- y_data_max * 1.15

  if (is.numeric(upper_ticks) && length(upper_ticks) > 0 &&
      all(is.finite(upper_ticks))) {
    y_max_auto <- max(
      y_max_auto,
      max(upper_ticks) * 1.02
    )
  }

  if (
    !is.finite(y_max_auto) ||
    y_max_auto <= 0
  ) {
    y_max_auto <- 1
  }

  break_info <- make_manual_y_break(
    df,
    break_range
  )

  if (!is.numeric(lower_ticks) || length(lower_ticks) < 1 ||
      any(!is.finite(lower_ticks))) {
    stop(
      sprintf(
        "[%s] lower_ticks 必须是至少包含1个有效数值的数值向量。",
        plotname
      )
    )
  }

  if (!is.numeric(upper_ticks) || length(upper_ticks) < 1 ||
      any(!is.finite(upper_ticks))) {
    stop(
      sprintf(
        "[%s] upper_ticks 必须是至少包含1个有效数值的数值向量。",
        plotname
      )
    )
  }

  if (!is.null(break_info)) {

    lower_breaks <- sort(unique(lower_ticks))
    upper_breaks <- sort(unique(upper_ticks))

    lower_breaks <- lower_breaks[
      lower_breaks >= 0 &
        lower_breaks <= break_info$lower
    ]

    upper_breaks <- upper_breaks[
      upper_breaks >= break_info$upper
    ]

    y_max_auto <- max(
      y_max_auto,
      max(upper_breaks) * 1.02
    )

    y_axis_breaks <- sort(
      unique(
        c(
          lower_breaks,
          upper_breaks
        )
      )
    )

  } else {

    y_axis_breaks <- sort(
      unique(
        c(lower_ticks, upper_ticks)
      )
    )

  }

  cat(
    sprintf(
      "[%s] 数据实际最大值 = %.3f；Y轴上限(含顶部留白) = %.3f\n",
      plotname,
      y_data_max,
      y_max_auto
    )
  )

  if (!is.null(break_info)) {

    cat(
      sprintf(
        "[%s] 下半段刻度 = %s；上半段刻度 = %s
",
        plotname,
        paste(lower_breaks, collapse = ", "),
        paste(upper_breaks, collapse = ", ")
      )
    )

    cat(
      sprintf(
        paste0(
          "[%s] 固定断轴区间 = %.3f ~ %.3f；",
          "该区间内共有 %d 个真实数据点不会显示；",
          "下半段刻度 = %s；",
          "上段比例 = %.2f\n"
        ),
        plotname,
        break_info$lower,
        break_info$upper,
        break_info$n_cut_points,
        paste(lower_breaks, collapse = ", "),
        break_info$upper_scale
      )
    )

    if (nrow(break_info$cut_boxes) > 0) {
      cat(
        sprintf(
          "[%s] 警告：断轴穿过箱体：%s\n",
          plotname,
          paste(
            as.character(break_info$cut_boxes$ASV),
            collapse = ", "
          )
        )
      )
    }

  } else {

    cat(
      sprintf(
        "[%s] 当前指定区间无法形成有效断轴，保留连续Y轴。\n",
        plotname
      )
    )

  }

  if (is.null(break_info)) {
    stop(
      sprintf(
        "[%s] 当前固定断轴区间无法使用，请检查 break_range。",
        plotname
      )
    )
  }

  break_low <- break_info$lower
  break_high <- break_info$upper
  upper_scale <- break_info$upper_scale

  if (length(lower_breaks) >= 2) {

    lower_tick_step <- min(
      diff(lower_breaks)
    )

  } else {

    lower_tick_step <- max(
      break_low,
      1
    )

  }

  gap_height <- max(
    lower_tick_step * 0.25,
    0.15
  )

  x_axis_bottom_gap <- lower_tick_step * 0.18

  y_plot_min <- -x_axis_bottom_gap

  map_broken_y <- function(y) {

    y <- as.numeric(y)

    out <- rep(
      NA_real_,
      length(y)
    )

    lower_id <- is.finite(y) &
      y <= break_low

    upper_id <- is.finite(y) &
      y >= break_high

    out[lower_id] <- y[lower_id]

    out[upper_id] <-
      break_low +
      gap_height +
      (y[upper_id] - break_high) * upper_scale

    out

  }

  box_rows <- vector(
    "list",
    length(valid_asv)
  )

  outlier_rows <- list()

  for (i in seq_along(valid_asv)) {

    current_asv <- valid_asv[i]

    current_id <- which(
      as.character(df$ASV) == current_asv &
        is.finite(df$Value)
    )

    current_values <- df$Value[
      current_id
    ]

    current_taxon <- as.character(
      df$Taxon[
        current_id[1]
      ]
    )

    q <- as.numeric(
      quantile(
        current_values,
        probs = c(
          0.25,
          0.50,
          0.75
        ),
        na.rm = TRUE,
        type = 7
      )
    )

    q1 <- q[1]
    med <- q[2]
    q3 <- q[3]

    iqr_value <- q3 - q1

    lower_fence <- q1 -
      1.5 * iqr_value

    upper_fence <- q3 +
      1.5 * iqr_value

    inside_values <- current_values[
      current_values >= lower_fence &
        current_values <= upper_fence
    ]

    if (length(inside_values) == 0) {

      whisker_low <- q1
      whisker_high <- q3

    } else {

      whisker_low <- min(
        inside_values,
        na.rm = TRUE
      )

      whisker_high <- max(
        inside_values,
        na.rm = TRUE
      )

    }

    current_outliers <- current_values[
      current_values < whisker_low |
        current_values > whisker_high
    ]

    box_rows[[i]] <- data.frame(
      x = i,
      ASV = current_asv,
      Taxon = current_taxon,
      ymin = whisker_low,
      lower = q1,
      middle = med,
      upper = q3,
      ymax = whisker_high,
      stringsAsFactors = FALSE
    )

    if (length(current_outliers) > 0) {

      outlier_rows[[
        length(outlier_rows) + 1
      ]] <- data.frame(
        x = i,
        ASV = current_asv,
        Taxon = current_taxon,
        Value = current_outliers,
        stringsAsFactors = FALSE
      )

    }

  }

  box_df <- bind_rows(
    box_rows
  )

  box_df$Taxon <- factor(
    box_df$Taxon,
    levels = legend_taxa
  )

  if (length(outlier_rows) > 0) {

    outlier_df <- bind_rows(
      outlier_rows
    )

  } else {

    outlier_df <- data.frame(
      x = numeric(0),
      ASV = character(0),
      Taxon = character(0),
      Value = numeric(0),
      stringsAsFactors = FALSE
    )

  }

  cut_box_body <- box_df[
    box_df$lower < break_high &
      box_df$upper > break_low,
    ,
    drop = FALSE
  ]

  if (nrow(cut_box_body) > 0) {

    stop(
      sprintf(
        paste0(
          "[%s] 固定断轴 %.3f~%.3f 会穿过箱体主体：%s。",
          "为避免生成错误图形，已停止绘图。"
        ),
        plotname,
        break_low,
        break_high,
        paste(
          cut_box_body$ASV,
          collapse = ", "
        )
      )
    )

  }

  box_df$lower_plot <- map_broken_y(
    box_df$lower
  )

  box_df$middle_plot <- map_broken_y(
    box_df$middle
  )

  box_df$upper_plot <- map_broken_y(
    box_df$upper
  )

  make_visible_vertical_segment <- function(
    x_value,
    y1,
    y2
  ) {

    y_min <- min(
      y1,
      y2
    )

    y_max <- max(
      y1,
      y2
    )

    result <- list()

    if (y_min <= break_low) {

      lower_end <- min(
        y_max,
        break_low
      )

      if (lower_end >= y_min) {

        result[[
          length(result) + 1
        ]] <- data.frame(
          x = x_value,
          xend = x_value,
          y = map_broken_y(y_min),
          yend = map_broken_y(lower_end)
        )

      }

    }

    if (y_max >= break_high) {

      upper_start <- max(
        y_min,
        break_high
      )

      if (y_max >= upper_start) {

        result[[
          length(result) + 1
        ]] <- data.frame(
          x = x_value,
          xend = x_value,
          y = map_broken_y(upper_start),
          yend = map_broken_y(y_max)
        )

      }

    }

    if (length(result) == 0) {

      return(
        NULL
      )

    }

    bind_rows(
      result
    )

  }

  whisker_rows <- list()

  for (i in seq_len(nrow(box_df))) {

    lower_whisker <- make_visible_vertical_segment(
      box_df$x[i],
      box_df$ymin[i],
      box_df$lower[i]
    )

    if (!is.null(lower_whisker)) {

      whisker_rows[[
        length(whisker_rows) + 1
      ]] <- lower_whisker

    }

    upper_whisker <- make_visible_vertical_segment(
      box_df$x[i],
      box_df$upper[i],
      box_df$ymax[i]
    )

    if (!is.null(upper_whisker)) {

      whisker_rows[[
        length(whisker_rows) + 1
      ]] <- upper_whisker

    }

  }

  whisker_df <- bind_rows(
    whisker_rows
  )

  if (nrow(outlier_df) > 0) {

    outlier_df <- outlier_df[
      outlier_df$Value <= break_low |
        outlier_df$Value >= break_high,
      ,
      drop = FALSE
    ]

    outlier_df$y_plot <- map_broken_y(
      outlier_df$Value
    )

  }

  y_plot_max <- map_broken_y(
    y_max_auto
  )

  lower_tick_plot <- map_broken_y(
    lower_breaks
  )

  upper_tick_plot <- map_broken_y(
    upper_breaks
  )

  y_break_positions <- c(
    lower_tick_plot,
    upper_tick_plot
  )

  y_break_labels <- c(
    as.character(lower_breaks),
    as.character(upper_breaks)
  )

  x_axis_left <- 1 - 0.62

  x_axis_right <- length(valid_asv) +
    0.62

  box_half_width <- 0.82 / 2

  box_df$xmin <- box_df$x -
    box_half_width

  box_df$xmax <- box_df$x +
    box_half_width

  slash_half_width <- 0.075

  slash_half_height <- gap_height * 0.16

  slash_centers <- c(
    break_low,
    break_low + gap_height
  )

  axis_break_mark_df <- data.frame(
    x = rep(
      x_axis_left - slash_half_width,
      2
    ),
    xend = rep(
      x_axis_left + slash_half_width,
      2
    ),
    y = slash_centers + slash_half_height,
    yend = slash_centers - slash_half_height
  )

  y_axis_segment_df <- data.frame(
    x = c(
      x_axis_left,
      x_axis_left
    ),
    xend = c(
      x_axis_left,
      x_axis_left
    ),
    y = c(
      y_plot_min,
      break_low + gap_height
    ),
    yend = c(
      break_low,
      y_plot_max
    )
  )

  y_tick_values <- c(
    lower_breaks[lower_breaks < break_low],
    upper_breaks[upper_breaks > break_high]
  )

  y_tick_positions <- map_broken_y(
    y_tick_values
  )

  y_tick_length <- 0.08

  y_tick_df <- data.frame(
    x = rep(
      x_axis_left - y_tick_length,
      length(y_tick_positions)
    ),
    xend = rep(
      x_axis_left,
      length(y_tick_positions)
    ),
    y = y_tick_positions,
    yend = y_tick_positions
  )

  p <- ggplot() +

  geom_rect(
    data = box_df,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = lower_plot,
      ymax = upper_plot,
      fill = Taxon
    ),
    colour = "#404040",
    linewidth = 0.3,
    alpha = 1,
    show.legend = TRUE
  ) +

  geom_segment(
    data = box_df,
    aes(
      x = xmin,
      xend = xmax,
      y = middle_plot,
      yend = middle_plot
    ),
    inherit.aes = FALSE,
    colour = "#404040",
    linewidth = 0.3
  ) +

  geom_segment(
    data = whisker_df,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    colour = "#404040",
    linewidth = 0.3
  ) +

  geom_point(
    data = outlier_df,
    aes(
      x = x,
      y = y_plot
    ),
    inherit.aes = FALSE,
    shape = 16,
    colour = "black",
    size = 1.8,
    alpha = 1
  ) +

  geom_segment(
    data = y_axis_segment_df,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.55
  ) +

  geom_segment(
    data = y_tick_df,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.45,
    lineend = "butt"
  ) +

  geom_segment(
    data = axis_break_mark_df,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.55,
    lineend = "butt"
  ) +

  scale_fill_manual(
    values = panel_pal,
    breaks = legend_taxa,
    limits = legend_taxa,
    drop = FALSE
  ) +

  scale_x_continuous(
    limits = c(
      x_axis_left,
      x_axis_right
    ),
    oob = scales::oob_keep,
    breaks = seq_along(
      valid_asv
    ),
    labels = valid_asv,
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  scale_y_continuous(
    limits = c(
      y_plot_min,
      y_plot_max
    ),
    breaks = y_break_positions,
    labels = y_break_labels,
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  labs(
    x = NULL,
    y = ylab_text,
    fill = NULL
  ) +

  guides(
    fill = guide_legend(
      keywidth = unit(
        0.8,
        "cm"
      ),
      keyheight = unit(
        0.9,
        "cm"
      ),
      byrow = TRUE
    )
  ) +

  theme_classic() +

    theme(

      axis.title.y = element_text(
        size = 16,
        margin = margin(
          r = 7
        )
      ),

      axis.text.x = element_text(
        size = 12,
        angle = 45,
        hjust = 1,
        vjust = 1
      ),

      axis.text.y = element_text(
        size = 13
      ),

      axis.line.y.left = element_blank(),

      axis.ticks.y.left = element_blank(),

      axis.title.y.right = element_blank(),
      axis.text.y.right = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.line.y.right = element_blank(),

      legend.position = "right",
      legend.title = element_blank(),

      legend.text = element_text(
        size = 12
      ),

      legend.key = element_rect(
        colour = "black",
        linewidth = 0.3
      ),

      legend.spacing.x = unit(
        0.05,
        "cm"
      ),

      legend.spacing.y = unit(
        0.02,
        "cm"
      ),

      plot.margin = margin(
        10,
        5,
        10,
        5,
        unit = "mm"
      )

    ) +

    coord_cartesian(
      clip = "off"
    )

  pdf_file <- file.path(
    dir_out,
    paste0(
      plotname,
      "_box.pdf"
    )
  )

  ggsave(
    filename = pdf_file,
    plot = p,
    device = grDevices::pdf,
    width = 8.5,
    height = 4.8,
    units = "in",
    useDingbats = FALSE,
    bg = "white"
  )

  png_file <- file.path(

    dir_out,

    paste0(
      plotname,
      "_box.png"
    )

  )

  ggsave(

    filename = png_file,

    plot = p,

    width = 8.5,

    height = 4.8,

    units = "in",

    dpi = 600,

    bg = "white"

  )

  print(p)

  cat(
    sprintf(
      "[%s] 已保存：\nPDF：%s\nPNG：%s\n",
      plotname,
      pdf_file,
      png_file
    )
  )

  invisible(p)

}

plot_out <- list()

for (t in task_list) {

  if (
    is.null(
      data_list[[t$name]]
    )
  ) {

    cat(
      "跳过 ",
      t$name,
      "：读取失败\n"
    )

    next

  }

  res <- data_list[[t$name]]

  plot_out[[t$name]] <- draw_boxplot(

    df = res$df,

    orig_asv_order = res$asv_order,

    taxon_map = res$taxon_map,

    dir_out = t$dir_out,

    plotname = t$name,

    ylab_text = t$ylab,

    break_range = t$break_range,

    lower_ticks = t$lower_ticks,

    upper_ticks = t$upper_ticks,

    global_color_map = global_taxon_colors

  )

}
