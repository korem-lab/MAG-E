suppressPackageStartupMessages({
  library(arrow)
  library(tidyverse)
  library(lmerTest)   # also attaches lme4
  library(emmeans)
})
emm_options(pval.digits = 6)

# ---- args ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript script.R <input.parquet> <out_dir>", call. = FALSE)
}
in_file <- args[1]
out_dir <- args[2]

# ---- data ----
factor_cols <- c("assembler", "coverage", "binner", "binning_mode", "target", "genome", "metric")

df <- read_parquet(in_file) |>
  mutate(value = as.numeric(value), across(all_of(factor_cols), factor))


# ---- analysis ----
metric_names <- c(cov_fs = "fscore", cov_pr = "precision", cov_rc = "recall")
pipeline_factors <- c("binner", "assembler", "binning_mode")

write_out <- function(x, prefix, label) {
  write_csv(x, file.path(out_dir, paste0(prefix, "_", label, ".csv")))
}

analyze_metric <- function(df, met) {
  label  <- metric_names[[met]]
  df_run <- filter(df, metric == met)
  if (met == "cov_pr") df_run$value <- 1 - df_run$value

  m <- lmer(value ~ binner * assembler * binning_mode + (1 | target) + (1 | genome),
            data = df_run)

  # marginal means for each pipeline component
  emms <- map(set_names(pipeline_factors), \(f) emmeans(m, as.formula(paste("~", f))))
  walk(emms, \(e) { print(e); print(pairs(e)) })

  categories <- imap_dfr(emms, \(e, f)
    as.data.frame(e) |> rename(cat = all_of(f)) |> mutate(config = f))
  pairwise <- imap_dfr(emms, \(e, f)
    as.data.frame(pairs(e)) |> mutate(config = f))

  write_out(categories, "factor_emmeans", label)
  write_out(pairwise,   "factor_contrasts",      label)

  # full pipeline combinations
  emm_all <- emmeans(m, ~ binner * assembler * binning_mode)
  write_out(as.data.frame(emm_all),        "pipeline_emmeans", label)
  write_out(as.data.frame(pairs(emm_all)), "pipeline_contrasts",     label)
}

walk(names(metric_names), \(met) analyze_metric(df, met))