# Load packages
library(lme4)
library(tidyverse)
library(lmerTest)
library(emmeans)
library(dplyr)

# read in df
args <- commandArgs(trailingOnly = True)
usage <- "Usage: Rscript <exe> <path> <exclude_refiners ? true|false>"
if (length(args) < 2) {
  stop(paste("Not enough arguments.\n", usage), call. = FALSE)
}
path <- args[1]
exclude_refiners <- tolower(args[2]) == "true"
df <- read.csv(paste(path, "per_genome_metrics.csv", sep="/")

# prep data
df <- mutate(df,
             value = as.numeric(value),
             binner = factor(binner),
             sample = factor(sample),
             is_refiner = factor(is_refiner),
             assembler = factor(assembler),
             genome = factor(genome),
             binning_mode = factor(binning_mode),
             metric = factor(metric)
)
# drop refiners from analysis
if (exclude_refiners) {
    df <- df[!df$is_refiner,]
}

# core analysis function
analyze_data = function(df, metric) {
  df_run = df[df$metric == metric,]
  if (metric == 'cov_pr') {
    df_run$value = 1-df_run$value
  }
  
  # make model
  m <- lmer(value ~ binner * assembler * binning_mode + (1|sample) + (1|genome) , data=df_run)
  emm_options(pval.digits=6)
  if (metric == 'cov_fs') {
    mout = 'fscore'
  }
  else if (metric == 'cov_pr') {
    mout = 'precision'
  }
  else {
    mout = 'recall'
  }

  # add the binner
  emm <- emmeans(m, ~ binner)
  print(emm)
  print(pairs(emm))
  res_df <- rename(as.data.frame(emm), cat=binner)
  res_df$config = rep('binner', nrow(res_df))

  # add the assembler
  emm <- emmeans(m, ~ assembler)
  print(pairs(emm))
  res_df_tmp <- rename(as.data.frame(emm), cat=assembler)
  res_df_tmp$config = rep('assembler', nrow(res_df_tmp))
  res_df = bind_rows(res_df, res_df_tmp)

  # add the binning_mode
  emm <- emmeans(m, ~ binning_mode)
  print(emm)
  print(pairs(emm))
  res_df_tmp <- rename(as.data.frame(emm), cat=binning_mode)
  res_df_tmp$config = rep('binning_mode', nrow(res_df_tmp))
  res_df = bind_rows(res_df, res_df_tmp)
  
  # write out
  write.csv(res_df, paste0(path, '/optimal_pipeline_cat_', mout, '.csv'), quote=FALSE)

  # Testing the best pipeline overall
  emm <- emmeans(m, ~ binner * assembler * binning_mode)
  emm_df <- as.data.frame(emm)
  if (metric == 'cov_pr') {
    best <- emm_df[which.min(emm_df$emmean),]
  }
  else {
    best <- emm_df[which.max(emm_df$emmean),]
  }
  pairs_df = as.data.frame(pairs(emm))
  best_label = with(best, paste(binner, assembler, binning_mode, sep= " "))
  against_best = subset(pairs_df, grepl(best_label, contrast))
  write.csv(against_best, paste0(path, '/optimal_pipeline_vs_best_', mout, '.csv'), quote=FALSE)
  write.csv(emm_df, paste0(path, '/all_pipelines_', mout, '.csv'), quote = FALSE)
  write.csv(pairs_df, paste0(path, '/all_pairs_', mout, '.csv'), quote=FALSE)
}

analyze_data(df, 'cov_fs')
analyze_data(df, 'cov_pr')
analyze_data(df, 'cov_rc')
