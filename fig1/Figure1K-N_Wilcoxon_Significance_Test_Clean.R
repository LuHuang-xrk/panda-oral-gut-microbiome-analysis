library(readxl)
library(cellranger)

run_test <- function(file_path, oral_prefix, fecal_prefix, species, panel, asv_order) {
  df <- read_excel(file_path)

  oral_cols <- grep(paste0("^", oral_prefix), colnames(df), value = FALSE)
  fecal_cols <- grep(paste0("^", fecal_prefix), colnames(df), value = FALSE)

  if (length(oral_cols) == 0) stop(paste0("未识别到", oral_prefix, "列"))
  if (length(fecal_cols) == 0) stop(paste0("未识别到", fecal_prefix, "列"))

  asv_ids <- df[[1]]
  keep <- asv_ids %in% asv_order
  df <- df[keep, , drop = FALSE]
  asv_ids <- df[[1]]

  p_values <- sapply(seq_len(nrow(df)), function(i) {
    oral <- as.numeric(unlist(df[i, oral_cols]))
    fecal <- as.numeric(unlist(df[i, fecal_cols]))

    oral <- oral[!is.na(oral)]
    fecal <- fecal[!is.na(fecal)]

    if (length(oral) == 0 || length(fecal) == 0) return(NA)

    wilcox.test(oral, fecal, exact = FALSE)$p.value
  })

  fdr_values <- p.adjust(p_values, method = "BH")

  result <- data.frame(
    Panel = panel,
    Species = species,
    ASV_ID = asv_ids,
    p_value = p_values,
    FDR = fdr_values
  )

  result$significance <- cut(
    result$FDR,
    breaks = c(-Inf, 1e-4, 1e-3, 1e-2, 5e-2, Inf),
    labels = c("****", "***", "**", "*", "ns")
  )

  result$ASV_ID <- factor(result$ASV_ID, levels = asv_order)
  result <- result[order(result$ASV_ID), ]
  result$ASV_ID <- as.character(result$ASV_ID)

  result
}

task_list <- list(
  list(
    file = "GP relative abundance.xlsx",
    oral = "GPO",
    fecal = "GPF",
    species = "GP",
    panel = "K",
    asv_order = c("ASV1", "ASV2", "ASV4", "ASV50", "ASV415")
  ),
  list(
    file = "RP relative abundance.xlsx",
    oral = "RPO",
    fecal = "RPF",
    species = "RP",
    panel = "L",
    asv_order = c("ASV1", "ASV2", "ASV4", "ASV5", "ASV52")
  ),
  list(
    file = "MD relative abundance.xlsx",
    oral = "MDO",
    fecal = "MDF",
    species = "MD",
    panel = "M",
    asv_order = c("ASV44", "ASV51", "ASV81", "ASV84", "ASV90")
  ),
  list(
    file = "RB relative abundance.xlsx",
    oral = "RBO",
    fecal = "RBF",
    species = "RB",
    panel = "N",
    asv_order = c("ASV28", "ASV30", "ASV67", "ASV128", "ASV3")
  )
)

result_all <- do.call(
  rbind,
  lapply(task_list, function(x) {
    run_test(
      file_path = x$file,
      oral_prefix = x$oral,
      fecal_prefix = x$fecal,
      species = x$species,
      panel = x$panel,
      asv_order = x$asv_order
    )
  })
)

write.csv(
  result_all,
  "Figure1K-N_Wilcoxon_FDR_significance.csv",
  row.names = FALSE
)
