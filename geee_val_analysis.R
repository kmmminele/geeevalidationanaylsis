# Copyright (c) 2026 Katharina Mminele
# R script accompanying the publication:
# Longitudinal Validation of the Generic Rating Scale for Previous Treatment Experience, 
# Treatment Expectations, and Treatment Effects (G-EEE) in Endometriosis Patients
#
# All rights reserved.

# Reproducibility note:
# This code reflects the analyses used to generate the reported results.
# The original placement of random seeds, the ordering of variable pairs,
# and the sequence of bootstrap analyses were retained to ensure exact
# reproduction of the reported bootstrap confidence intervals. Some
# computations are therefore presented separately rather than being
# consolidated into a single function.

# The provided code covers all statistical analyses.
# Code used solely for formatting results for tables and figures is not included.

# Data are available from the corresponding author upon reasonable request.

# ============================================================
# 0. Packages
# ============================================================

library(tidyverse)
library(psych)
library(flextable)
library(boot)
library(lme4)
library(correlation)
library(tidybayes)
library(brms)
library(bayestestR)
library(performance)

# ============================================================
# 1. Multidimensionality of Improvement, worsening, and Side effects
# Residual correlations, clustered multilevel correlations, between-person associations
# ============================================================

# 1. Data preparation
# smaller data frame with expectation / effect variables only

mlm_data_exp <- geee_val_data %>% 
  select(VPN_Code_corr,
         contains("Geetreat_exp_")) %>% 
  pivot_longer(cols = -VPN_Code_corr,
               names_to = c("item", "time"),
               values_to = "value",
               names_pattern = "Geetreat_exp_(.*)_rec_T(.*)") %>% 
  mutate(item = factor(item, levels = c("01", "02", "03"), labels = c("improvement", "worsening", "side_effect")),
         time = factor(time, levels = c("0", "1", "2", "3", "4", "5", "6"))) %>% 
  pivot_wider(names_from = "item",
              values_from = "value")

mlm_data_eff <- geee_val_data %>% 
  select(VPN_Code_corr,
         contains("Gee_aktuell_")) %>% 
  pivot_longer(cols = -VPN_Code_corr,
               names_to = c("item", "time"),
               values_to = "value",
               names_pattern = "Gee_aktuell_(.*)_T(.*)") %>% 
  
  filter(time != "0") %>% 
  mutate(item = factor(item, levels = c("01", "02", "03"), labels = c("improvement", "worsening", "side_effect")),
         time = factor(time, levels = c("1", "2", "3", "4", "5", "6"))) %>% 
  pivot_wider(names_from = "item",
              values_from = "value")

# Due to the forced-choice response format, improvement, worsening, and
# side-effect ratings are always missing jointly. 

mlm_data_exp_complete <- mlm_data_exp %>% 
  filter(!is.na(improvement))
mlm_data_eff_complete <- mlm_data_eff %>% 
  filter(!is.na(improvement))


# 2. Extraction of residuals from mixed-effects models

m_imp_exp <- lmer(improvement ~ time + (1|VPN_Code_corr), data = mlm_data_exp_complete)
mlm_data_exp_complete$imp_res <- resid(m_imp_exp)

m_worse_exp <- lmer(worsening ~ time + (1|VPN_Code_corr), data = mlm_data_exp_complete)
mlm_data_exp_complete$wor_res <- resid(m_worse_exp)

m_side_exp <- lmer(side_effect ~ time + (1|VPN_Code_corr), data = mlm_data_exp_complete)
mlm_data_exp_complete$side_res <- resid(m_side_exp)

m_imp_eff <- lmer(improvement ~ time + (1|VPN_Code_corr), data = mlm_data_eff_complete)
mlm_data_eff_complete$imp_res <- resid(m_imp_eff)

m_worse_eff <- lmer(worsening ~ time + (1|VPN_Code_corr), data = mlm_data_eff_complete)
mlm_data_eff_complete$wor_res <- resid(m_worse_eff)

m_side_eff <- lmer(side_effect ~ time + (1|VPN_Code_corr), data = mlm_data_eff_complete)
mlm_data_eff_complete$side_res <- resid(m_side_eff)


# 3. Time-specific residual correlations and bootstrap CIs

boot_ci_rho <- function(df, x, y, R = 2000, method = "spearman") {
  
  df_xy <- df %>% select(all_of(c(x, y))) %>% drop_na()
  n_pair <- nrow(df_xy)

  # too little data -> NA
   if (n_pair < 5) {
     return(tibble(n = n_pair, rho = NA_real_, CI_low = NA_real_, CI_high = NA_real_))
   }
  
  b <- boot(
    data = df_xy,
    statistic = function(data, idx) {
      cor(data[idx, 1], data[idx, 2], method = method)
    },
    R = R
  )
  
  # Percentile CI
  ci <- suppressWarnings(boot.ci(b, type = "perc"))
  
  tibble(
    n = n_pair,
    rho = as.numeric(b$t0),
    CI_low = ci$percent[4],
    CI_high = ci$percent[5]
  )
}

set.seed(1649)

exp_cor_tbl <- list()

for(i in unique(mlm_data_exp_complete$time)){
  df <- mlm_data_exp_complete %>% filter(time == i)
  
  cor_tbl <- as.data.frame(correlation(data = df, select=c("imp_res", "wor_res", "side_res"), method= "spearman"))
  cor_tbl$time = i
  cor_tbl$n = nrow(df)
  
  tmp <- boot_ci_rho(df, cor_tbl[1,1], cor_tbl[1,2])
  cor_tbl[1, c("rho_boot","CI_low_boot","CI_high_boot","n_boot")] <- tmp[1, c("rho","CI_low","CI_high","n")]
  tmp <- boot_ci_rho(df, cor_tbl[2,1], cor_tbl[2,2])
  cor_tbl[2, c("rho_boot","CI_low_boot","CI_high_boot","n_boot")] <- tmp[1, c("rho","CI_low","CI_high","n")]
  tmp <- boot_ci_rho(df, cor_tbl[3,1], cor_tbl[3,2])
  cor_tbl[3, c("rho_boot","CI_low_boot","CI_high_boot","n_boot")] <- tmp[1, c("rho","CI_low","CI_high","n")]
  
  exp_cor_tbl[[i]] <- cor_tbl
  
}

set.seed(1650)

eff_cor_tbl <- list()

for(i in unique(mlm_data_eff_complete$time)){
  df <- mlm_data_eff_complete %>% filter(time == i)
  
  cor_tbl <- as.data.frame(correlation(data = df, select=c("imp_res", "wor_res", "side_res"), method= "spearman"))
  cor_tbl$time = i
  cor_tbl$n = nrow(df)
  
  tmp <- boot_ci_rho(df, cor_tbl[1,1], cor_tbl[1,2])
  cor_tbl[1, c("rho_boot","CI_low_boot","CI_high_boot","n_boot")] <- tmp[1, c("rho","CI_low","CI_high","n")]
  tmp <- boot_ci_rho(df, cor_tbl[2,1], cor_tbl[2,2])
  cor_tbl[2, c("rho_boot","CI_low_boot","CI_high_boot","n_boot")] <- tmp[1, c("rho","CI_low","CI_high","n")]
  tmp <- boot_ci_rho(df, cor_tbl[3,1], cor_tbl[3,2])
  cor_tbl[3, c("rho_boot","CI_low_boot","CI_high_boot","n_boot")] <- tmp[1, c("rho","CI_low","CI_high","n")]
  
  eff_cor_tbl[[i]] <- cor_tbl
  
}

# 4. Between-person correlations based on person-level means

# previous treatment experience, assessed only in subsample of participants
# with treatment experience, assessed only at T0
personmeans_experience <- geee_val_data %>% 
  filter(GeeTreat_experience_dicho == 1) %>% 
  select(VPN_Code_corr, contains("Geetreat_experience") & !contains("01") & !contains("dicho")) %>% 
  mutate(improvement = Geetreat_experience_02_T0,
         worsening = Geetreat_experience_03_T0,
         side_effect = Geetreat_experience_04_T0) %>% 
  select(VPN_Code_corr, improvement, worsening, side_effect)
personmeans_exp <- aggregate(cbind(improvement, worsening, side_effect) ~ VPN_Code_corr, data = mlm_data_exp_complete, mean)
personmeans_eff <- aggregate(cbind(improvement, worsening, side_effect) ~ VPN_Code_corr, data = mlm_data_eff_complete, mean)


boot_spearman <- function(df, x, y, R = 2000, seed = 1){
  set.seed(seed)
  d <- df[, c("VPN_Code_corr", x, y)]
  d <- d[complete.cases(d), ]
  ids <- unique(d$VPN_Code_corr)
  
  r_boot <- replicate(R, {
    ids_b <- sample(ids, length(ids), replace = TRUE)
    d_b <- do.call(rbind, lapply(ids_b, function(id){
      d[d$VPN_Code_corr == id, ]
    }))
    suppressWarnings(cor(d_b[[x]], d_b[[y]], method = "spearman"))
  })
  
  r_hat <- suppressWarnings(cor(d[[x]], d[[y]], method = "spearman"))
  ci <- quantile(r_boot, c(.025, .975), na.rm = TRUE)
  list(r = r_hat, ci = ci, boot = r_boot)
}

make_between_cor_tbl <- function(personmeans_df,
                                 subscale,
                                 vars = c("improvement", "worsening", "side_effect"),
                                 id_col = "VPN_Code_corr",
                                 R = 2000,
                                 seed = 1103){
  
  df <- personmeans_df[, c(id_col, vars)]
  names(df)[1] <- "VPN_Code_corr"
  
  # pairs
  pair_mat <- combn(vars, 2)
  n_pairs <- ncol(pair_mat)
  
  out <- vector("list", n_pairs)
  
  for (k in seq_len(n_pairs)) {
    x <- pair_mat[1, k]
    y <- pair_mat[2, k]
    
    res <- boot_spearman(df, x, y, R = R, seed = seed + k)
    out[[k]] <- data.frame(
      subscale   = subscale,
      Parameter1 = x,
      Parameter2 = y,
      rho        = unname(res$r),
      low        = unname(res$ci[1]),
      high       = unname(res$ci[2])
    )
  }
  
  do.call(rbind, out)
}

between_experience <- make_between_cor_tbl(personmeans_experience, subscale = "experience",
                                           vars = c("improvement","worsening","side_effect"),
                                           R = 2000, seed = 1411)

between_exp <- make_between_cor_tbl(personmeans_exp, subscale = "expectations",
                                    vars = c("improvement","worsening","side_effect"),
                                    R = 2000, seed = 1103)

between_eff <- make_between_cor_tbl(personmeans_eff, subscale = "effects",
                                    vars = c("improvement","worsening","side_effect"),
                                    R = 2000, seed = 2103)

# 5. clustered multilevel correlations
mlm_data_exp_complete_cluster <- mlm_data_exp_complete %>%
  group_by(VPN_Code_corr) %>%
  mutate(
    imp_mean  = mean(improvement,  na.rm = TRUE),
    wor_mean  = mean(worsening,  na.rm = TRUE),
    side_mean = mean(side_effect, na.rm = TRUE)
  )

mlm_data_eff_complete_cluster <- mlm_data_eff_complete %>%
  group_by(VPN_Code_corr) %>%
  mutate(
    imp_mean  = mean(improvement,  na.rm = TRUE),
    wor_mean  = mean(worsening,  na.rm = TRUE),
    side_mean = mean(side_effect, na.rm = TRUE)
  )

# look at residual correlations
within_cluster_exp <- misty::multilevel.cor(
  mlm_data_exp_complete_cluster,
  imp_res, wor_res, side_res,
  imp_mean, wor_mean, side_mean,
  cluster = "VPN_Code_corr"
)

within_cluster_eff <- misty::multilevel.cor(
  mlm_data_eff_complete_cluster,
  imp_res, wor_res, side_res,
  imp_mean, wor_mean, side_mean,
  cluster = "VPN_Code_corr"
)


# ============================================================
# 2. Cluster bootstrap for convergent validity
# ============================================================

# 1. Time-specific correlations and bootstrap confidence intervals

# 1.1 Initialize output table
expectations_cv_ci <- data.frame("time" = c("T0", "T1", "T2", "T3", "T4", "T5", "T6"),
                                 "G_improve_NRS_improve_n" = rep(NA, times = 7),
                                 "G_improve_NRS_improve_r" = rep(NA, times = 7),
                                 "G_improve_NRS_improve_lower" = rep(NA, times = 7),
                                 "G_improve_NRS_improve_upper" = rep(NA, times = 7),
                                 "G_improve_NRS_worse_n" = rep(NA, times = 7),
                                 "G_improve_NRS_worse_r" = rep(NA, times = 7),
                                 "G_improve_NRS_worse_lower" = rep(NA, times = 7),
                                 "G_improve_NRS_worse_upper" = rep(NA, times = 7),
                                 "G_improve_PDI_sum_n" = rep(NA, times = 7),
                                 "G_improve_PDI_sum_r" = rep(NA, times = 7),
                                 "G_improve_PDI_sum_lower" = rep(NA, times = 7),
                                 "G_improve_PDI_sum_upper" = rep(NA, times = 7),
                                 "G_worse_NRS_improve_n" = rep(NA, times = 7),
                                 "G_worse_NRS_improve_r" = rep(NA, times = 7),
                                 "G_worse_NRS_improve_lower" = rep(NA, times = 7),
                                 "G_worse_NRS_improve_upper" = rep(NA, times = 7),
                                 "G_worse_NRS_worse_n" = rep(NA, times = 7),
                                 "G_worse_NRS_worse_r" = rep(NA, times = 7),
                                 "G_worse_NRS_worse_lower" = rep(NA, times = 7),
                                 "G_worse_NRS_worse_upper" = rep(NA, times = 7),
                                 "G_worse_PDI_sum_n" = rep(NA, times = 7),
                                 "G_worse_PDI_sum_r" = rep(NA, times = 7),
                                 "G_worse_PDI_sum_lower" = rep(NA, times = 7),
                                 "G_worse_PDI_sum_upper" = rep(NA, times = 7),
                                 "G_side_NRS_improve_n" = rep(NA, times = 7),
                                 "G_side_NRS_improve_r" = rep(NA, times = 7),
                                 "G_side_NRS_improve_lower" = rep(NA, times = 7),
                                 "G_side_NRS_improve_upper" = rep(NA, times = 7),
                                 "G_side_NRS_worse_n" = rep(NA, times = 7),
                                 "G_side_NRS_worse_r" = rep(NA, times = 7),
                                 "G_side_NRS_worse_lower" = rep(NA, times = 7),
                                 "G_side_NRS_worse_upper" = rep(NA, times = 7),
                                 "G_side_PDI_sum_n" = rep(NA, times = 7),
                                 "G_side_PDI_sum_r" = rep(NA, times = 7),
                                 "G_side_PDI_sum_lower" = rep(NA, times = 7),
                                 "G_side_PDI_sum_upper" = rep(NA, times = 7)
)

effects_cv_ci <- data.frame("time" = c("T1", "T2", "T3", "T4", "T5", "T6"),
                            "G_improve_NRS_n" = rep(NA, times = 6),
                            "G_improve_NRS_r" = rep(NA, times = 6),
                            "G_improve_NRS_lower" = rep(NA, times = 6),
                            "G_improve_NRS_upper" = rep(NA, times = 6),
                            "G_improve_PDI_sum_n" = rep(NA, times = 6),
                            "G_improve_PDI_sum_r" = rep(NA, times = 6),
                            "G_improve_PDI_sum_lower" = rep(NA, times = 6),
                            "G_improve_PDI_sum_upper" = rep(NA, times = 6),
                            "G_worse_NRS_n" = rep(NA, times = 6),
                            "G_worse_NRS_r" = rep(NA, times = 6),
                            "G_worse_NRS_lower" = rep(NA, times = 6),
                            "G_worse_NRS_upper" = rep(NA, times = 6),
                            "G_worse_PDI_sum_n" = rep(NA, times = 6),
                            "G_worse_PDI_sum_r" = rep(NA, times = 6),
                            "G_worse_PDI_sum_lower" = rep(NA, times = 6),
                            "G_worse_PDI_sum_upper" = rep(NA, times = 6),
                            "G_side_NRS_n" = rep(NA, times = 6),
                            "G_side_NRS_r" = rep(NA, times = 6),
                            "G_side_NRS_lower" = rep(NA, times = 6),
                            "G_side_NRS_upper" = rep(NA, times = 6),
                            "G_side_PDI_sum_n" = rep(NA, times = 6),
                            "G_side_PDI_sum_r" = rep(NA, times = 6),
                            "G_side_PDI_sum_lower" = rep(NA, times = 6),
                            "G_side_PDI_sum_upper" = rep(NA, times = 6)
)

# 1.2 Calculate correlations at T0/1-T6
for(i in expectations_cv_ci$time){
  G_improve <- paste0("Geetreat_exp_01_rec_", i)
  G_worse <- paste0("Geetreat_exp_02_rec_", i)
  G_side <- paste0("Geetreat_exp_03_rec_", i)
  NRS_improve <- paste0("NRS_exp_verbesserung_mean_", i) 
  NRS_worse <- paste0("NRS_exp_verschlechterung_mean_", i)
  PDI_sum <- paste0("PDI_exp_sum_", i)
  
  #G_improve - NRS_improve
  d <- geee_val_data[, c(G_improve, NRS_improve)]
  
  d <- d[complete.cases(d), ]
  n_complete <- nrow(d)
  
  set.seed(659)
  
  boot_correlation <- boot(
    data = d,
    statistic = function(data, i) {
      cor(data[i, 1], data[i, 2], use = "complete.obs", method = "spearman")
    },
    R = 2000
  )
  
  ci <- boot.ci(boot_correlation, type = "perc")
  
  t0   <- round(boot_correlation$t0, 2)
  low  <- round(ci$percent[4], 2)
  high <- round(ci$percent[5], 2)
  
  expectations_cv_ci[expectations_cv_ci$time == i, 2:5] <- c(n_complete, t0, low, high)
  
  #G_improve_NRS_worse
  d <- geee_val_data[, c(G_improve, NRS_worse)]
  
  d <- d[complete.cases(d), ]
  n_complete <- nrow(d)
  
  set.seed(710)
  
  boot_correlation <- boot(
    data = d,
    statistic = function(data, i) {
      cor(data[i, 1], data[i, 2], use = "complete.obs", method = "spearman")
    },
    R = 2000
  )
  
  ci <- boot.ci(boot_correlation, type = "perc")
  
  t0   <- round(boot_correlation$t0, 2)
  low  <- round(ci$percent[4], 2)
  high <- round(ci$percent[5], 2)
  
  expectations_cv_ci[expectations_cv_ci$time == i, 6:9] <- c(n_complete, t0, low, high)
  
  #G_improve_PDI_sum
  d <- geee_val_data[, c(G_improve, PDI_sum)]
  
  d <- d[complete.cases(d), ]
  n_complete <- nrow(d)
  
  set.seed(711)
  
  boot_correlation <- boot(
    data = d,
    statistic = function(data, i) {
      cor(data[i, 1], data[i, 2], use = "complete.obs", method = "spearman")
    },
    R = 2000
  )
  
  ci <- boot.ci(boot_correlation, type = "perc")
  
  t0   <- round(boot_correlation$t0, 2)
  low  <- round(ci$percent[4], 2)
  high <- round(ci$percent[5], 2)
  
  expectations_cv_ci[expectations_cv_ci$time == i, 10:13] <- c(n_complete, t0, low, high)
  
  
  #G_worse - NRS_improve
  d <- geee_val_data[, c(G_worse, NRS_improve)]
  
  d <- d[complete.cases(d), ]
  n_complete <- nrow(d)
  
  set.seed(712)
  
  boot_correlation <- boot(
    data = d,
    statistic = function(data, i) {
      cor(data[i, 1], data[i, 2], use = "complete.obs", method = "spearman")
    },
    R = 2000
  )
  
  ci <- boot.ci(boot_correlation, type = "perc")
  
  t0   <- round(boot_correlation$t0, 2)
  low  <- round(ci$percent[4], 2)
  high <- round(ci$percent[5], 2)
  
  expectations_cv_ci[expectations_cv_ci$time == i, 14:17] <- c(n_complete, t0, low, high)
  
  #G_worse_NRS_worse
  d <- geee_val_data[, c(G_worse, NRS_worse)]
  
  d <- d[complete.cases(d), ]
  n_complete <- nrow(d)
  
  set.seed(713)
  
  boot_correlation <- boot(
    data = d,
    statistic = function(data, i) {
      cor(data[i, 1], data[i, 2], use = "complete.obs", method = "spearman")
    },
    R = 2000
  )
  
  ci <- boot.ci(boot_correlation, type = "perc")
  
  t0   <- round(boot_correlation$t0, 2)
  low  <- round(ci$percent[4], 2)
  high <- round(ci$percent[5], 2)
  
  expectations_cv_ci[expectations_cv_ci$time == i, 18:21] <- c(n_complete, t0, low, high)
  
  #G_worse_PDI_sum
  d <- geee_val_data[, c(G_worse, PDI_sum)]
  
  d <- d[complete.cases(d), ]
  n_complete <- nrow(d)
  
  set.seed(714)
  
  boot_correlation <- boot(
    data = d,
    statistic = function(data, i) {
      cor(data[i, 1], data[i, 2], use = "complete.obs", method = "spearman")
    },
    R = 2000
  )
  
  ci <- boot.ci(boot_correlation, type = "perc")
  
  t0   <- round(boot_correlation$t0, 2)
  low  <- round(ci$percent[4], 2)
  high <- round(ci$percent[5], 2)
  
  expectations_cv_ci[expectations_cv_ci$time == i, 22:25] <- c(n_complete, t0, low, high)
  
  # #G_side_NRS_improve
  d <- geee_val_data[, c(G_side, NRS_improve)]
  
  d <- d[complete.cases(d), ]
  n_complete <- nrow(d)
  
  set.seed(715)
  
  boot_correlation <- boot(
    data = d,
    statistic = function(data, i) {
      cor(data[i, 1], data[i, 2], use = "complete.obs", method = "spearman")
    },
    R = 2000
  )
  
  ci <- boot.ci(boot_correlation, type = "perc")
  
  t0   <- round(boot_correlation$t0, 2)
  low  <- round(ci$percent[4], 2)
  high <- round(ci$percent[5], 2)
  
  expectations_cv_ci[expectations_cv_ci$time == i, 26:29] <- c(n_complete, t0, low, high)
  
  #G_side_NRS_worse
  d <- geee_val_data[, c(G_side, NRS_worse)]
  
  d <- d[complete.cases(d), ]
  n_complete <- nrow(d)
  
  set.seed(716)
  
  boot_correlation <- boot(
    data = d,
    statistic = function(data, i) {
      cor(data[i, 1], data[i, 2], use = "complete.obs", method = "spearman")
    },
    R = 2000
  )
  
  ci <- boot.ci(boot_correlation, type = "perc")
  
  t0   <- round(boot_correlation$t0, 2)
  low  <- round(ci$percent[4], 2)
  high <- round(ci$percent[5], 2)
  
  expectations_cv_ci[expectations_cv_ci$time == i, 30:33] <- c(n_complete, t0, low, high)
  
  #G_side_PDI_sum
  d <- geee_val_data[, c(G_side, PDI_sum)]
  
  d <- d[complete.cases(d), ]
  n_complete <- nrow(d)
  
  set.seed(717)
  
  boot_correlation <- boot(
    data = d,
    statistic = function(data, i) {
      cor(data[i, 1], data[i, 2], use = "complete.obs", method = "spearman")
    },
    R = 2000
  )
  
  ci <- boot.ci(boot_correlation, type = "perc")
  
  t0   <- round(boot_correlation$t0, 2)
  low  <- round(ci$percent[4], 2)
  high <- round(ci$percent[5], 2)
  
  expectations_cv_ci[expectations_cv_ci$time == i, 34:37] <- c(n_complete, t0, low, high)
}

for(i in effects_cv_ci$time){
  G_improve <- paste0("Gee_aktuell_01_", i)
  G_worse <- paste0("Gee_aktuell_02_", i)
  G_side <- paste0("Gee_aktuell_03_", i)
  NRS <- paste0("NRS_mean_", i) 
  PDI_sum <- paste0("PDI_sum_", i)
  
  #G_improve - NRS
  d <- geee_val_data[, c(G_improve, NRS)]
  
  d <- d[complete.cases(d), ]
  n_complete <- nrow(d)
  
  set.seed(1335)
  
  boot_correlation <- boot(
    data = d,
    statistic = function(data, i) {
      cor(data[i, 1], data[i, 2], use = "complete.obs", method = "spearman")
    },
    R = 2000
  )
  
  ci <- boot.ci(boot_correlation, type = "perc")
  
  t0   <- round(boot_correlation$t0, 2)
  low  <- round(ci$percent[4], 2)
  high <- round(ci$percent[5], 2)
  
  effects_cv_ci[effects_cv_ci$time == i, 2:5] <- c(n_complete, t0, low, high)
  
  #G_improve_PDI_sum
  d <- geee_val_data[, c(G_improve, PDI_sum)]
  
  d <- d[complete.cases(d), ]
  n_complete <- nrow(d)
  
  set.seed(1338)
  
  boot_correlation <- boot(
    data = d,
    statistic = function(data, i) {
      cor(data[i, 1], data[i, 2], use = "complete.obs", method = "spearman")
    },
    R = 2000
  )
  
  ci <- boot.ci(boot_correlation, type = "perc")
  
  t0   <- round(boot_correlation$t0, 2)
  low  <- round(ci$percent[4], 2)
  high <- round(ci$percent[5], 2)
  
  effects_cv_ci[effects_cv_ci$time == i, 6:9] <- c(n_complete, t0, low, high)
  
  #G_worse - NRS
  d <- geee_val_data[, c(G_worse, NRS)]
  
  d <- d[complete.cases(d), ]
  n_complete <- nrow(d)
  
  set.seed(1339)
  
  boot_correlation <- boot(
    data = d,
    statistic = function(data, i) {
      cor(data[i, 1], data[i, 2], use = "complete.obs", method = "spearman")
    },
    R = 2000
  )
  
  ci <- boot.ci(boot_correlation, type = "perc")
  
  t0   <- round(boot_correlation$t0, 2)
  low  <- round(ci$percent[4], 2)
  high <- round(ci$percent[5], 2)
  
  effects_cv_ci[effects_cv_ci$time == i, 10:13] <- c(n_complete, t0, low, high)
  
  #G_worse_PDI_sum
  d <- geee_val_data[, c(G_worse, PDI_sum)]
  
  d <- d[complete.cases(d), ]
  n_complete <- nrow(d)
  
  set.seed(1340)
  
  boot_correlation <- boot(
    data = d,
    statistic = function(data, i) {
      cor(data[i, 1], data[i, 2], use = "complete.obs", method = "spearman")
    },
    R = 2000
  )
  
  ci <- boot.ci(boot_correlation, type = "perc")
  
  t0   <- round(boot_correlation$t0, 2)
  low  <- round(ci$percent[4], 2)
  high <- round(ci$percent[5], 2)
  
  effects_cv_ci[effects_cv_ci$time == i, 14:17] <- c(n_complete, t0, low, high)
  
  #G_side - NRS
  d <- geee_val_data[, c(G_side, NRS)]
  
  d <- d[complete.cases(d), ]
  n_complete <- nrow(d)
  
  set.seed(1341)
  
  boot_correlation <- boot(
    data = d,
    statistic = function(data, i) {
      cor(data[i, 1], data[i, 2], use = "complete.obs", method = "spearman")
    },
    R = 2000
  )
  
  ci <- boot.ci(boot_correlation, type = "perc")
  
  t0   <- round(boot_correlation$t0, 2)
  low  <- round(ci$percent[4], 2)
  high <- round(ci$percent[5], 2)
  
  effects_cv_ci[effects_cv_ci$time == i, 18:21] <- c(n_complete, t0, low, high)
  
  #G_side_PDI_sum
  d <- geee_val_data[, c(G_side, PDI_sum)]
  
  d <- d[complete.cases(d), ]
  n_complete <- nrow(d)
  
  set.seed(1335)
  
  boot_correlation <- boot(
    data = d,
    statistic = function(data, i) {
      cor(data[i, 1], data[i, 2], use = "complete.obs", method = "spearman")
    },
    R = 2000
  )
  
  ci <- boot.ci(boot_correlation, type = "perc")
  
  t0   <- round(boot_correlation$t0, 2)
  low  <- round(ci$percent[4], 2)
  high <- round(ci$percent[5], 2)
  
  effects_cv_ci[effects_cv_ci$time == i, 22:25] <- c(n_complete, t0, low, high)
  
}

# 1.3 Reshape time-specific results
expectations_cv_ci <- expectations_cv_ci %>% 
  pivot_longer(cols = -time,
               names_to = c("cor", ".value"),
               names_pattern = "([^_]+_[^_]+_[^_]+_[^_]+)_(n|r|lower|upper)")

effects_cv_ci <- effects_cv_ci %>% 
  pivot_longer(cols = -time,
               names_to = c("cor", ".value"),
               names_pattern = "(.+?)_(n|r|lower|upper)") 

# 2. Correlations between person-level means

# 2.1 Prepare longitudinal data
data_exp_GEEE <- geee_val_data %>% 
  select(VPN_Code_corr,
         contains("Geetreat_exp_")) %>% 
  pivot_longer(cols = -VPN_Code_corr,
               names_to = c("item", "time"),
               values_to = "value",
               names_pattern = "Geetreat_exp_(.*)_rec_T(.*)") %>% 
  mutate(item = factor(item, levels = c("01", "02", "03"), labels = c("improvement", "worsening", "side_effect")),
         time = factor(time, levels = c("0", "1", "2", "3", "4", "5", "6"))) %>% 
  pivot_wider(names_from = "item",
              values_from = "value")

data_exp_NRS <- geee_val_data %>% 
  select(VPN_Code_corr,
         contains("NRS_exp_") & !contains(c("max", "sum", "rec"))) %>% 
  pivot_longer(cols = -VPN_Code_corr,
               names_to = c("item", "time"),
               values_to = "value",
               names_pattern = "NRS_exp_(.*)_mean_T(.*)") %>% 
  mutate(item = factor(item, levels = c("verbesserung", "verschlechterung"), labels = c("NRS_exp_imp", "NRS_exp_worse")),
         time = factor(time, levels = c("0", "1", "2", "3", "4", "5", "6"))) %>% 
  pivot_wider(names_from = "item",
              values_from = "value")

data_exp_PDI_sum <- geee_val_data %>% 
  select(VPN_Code_corr,
         contains("PDI_exp_sum")) %>% 
  pivot_longer(cols = -VPN_Code_corr,
               names_to = c("time"),
               values_to = "PDI_exp",
               names_pattern = "PDI_exp_sum_T(.*)") %>% 
  mutate(time = factor(time, levels = c("0", "1", "2", "3", "4", "5", "6"))) 

data_exp <- data_exp_GEEE %>% 
  left_join(data_exp_NRS, by = c("VPN_Code_corr", "time")) %>% 
  left_join(data_exp_PDI_sum, by = c("VPN_Code_corr", "time"))

# Due to the forced-choice response format, improvement, worsening, and
# side-effect ratings are always missing jointly. 
data_exp <- data_exp %>% 
  filter(!is.na(improvement))

data_eff_GEEE <- geee_val_data %>% 
  select(VPN_Code_corr,
         contains("Gee_aktuell_")) %>% 
  pivot_longer(cols = -VPN_Code_corr,
               names_to = c("item", "time"),
               values_to = "value",
               names_pattern = "Gee_aktuell_(.*)_T(.*)") %>% 
  
  filter(time != "0") %>% 
  mutate(item = factor(item, levels = c("01", "02", "03"), labels = c("improvement", "worsening", "side_effect")),
         time = factor(time, levels = c("1", "2", "3", "4", "5", "6"))) %>% 
  pivot_wider(names_from = "item",
              values_from = "value")

data_eff_NRS <- geee_val_data %>% 
  select(VPN_Code_corr,
         contains("NRS_mean_")) %>% 
  pivot_longer(cols = -VPN_Code_corr,
               names_to = c("time"),
               values_to = "NRS",
               names_pattern = "NRS_mean_T(.*)") %>% 
  filter(time != "0") %>% 
  mutate(time = factor(time, levels = c("1", "2", "3", "4", "5", "6")))

data_eff_PDI_sum <- geee_val_data %>% 
  select(VPN_Code_corr,
         contains("PDI_sum_")) %>% 
  pivot_longer(cols = -VPN_Code_corr,
               names_to = c("time"),
               values_to = "PDI",
               names_pattern = "PDI_sum_T(.*)") %>% 
  filter(time != "0") %>% 
  mutate(time = factor(time, levels = c("1", "2", "3", "4", "5", "6")))

data_eff <- data_eff_GEEE %>% 
  left_join(data_eff_NRS, by = c("VPN_Code_corr", "time")) %>% 
  left_join(data_eff_PDI_sum, by = c("VPN_Code_corr", "time"))

# Due to the forced-choice response format, improvement, worsening, and
# side-effect ratings are always missing jointly. 
data_eff <- data_eff %>% 
  filter(!is.na(improvement))

# 2.2 Calculate person-level means
personmeans_exp <- aggregate(cbind(improvement, worsening, side_effect, NRS_exp_imp, NRS_exp_worse, PDI_exp)~VPN_Code_corr, data=data_exp, mean)

personmeans_eff <- aggregate(cbind(improvement, worsening, side_effect, NRS, PDI)~VPN_Code_corr, data=data_eff, mean)

# 2.3 Bootstrap correlations between person-level means
exp_cv <- data.frame(GEEE_var = rep(c("improvement", "worsening", "side_effect"), times = 3),
                     Ext_var = rep(c("NRS_exp_imp", "NRS_exp_worse", "PDI_exp"), each = 3),
                     r = NA,
                     lower = NA,
                     upper = NA)

set.seed(1753)

for(i in 1:nrow(exp_cv)){
  d <- personmeans_exp
  x <- exp_cv[i, "GEEE_var"]
  y <- exp_cv[i, "Ext_var"]
  
  boot_correlation <- boot(
    data = d,
    statistic = function(data, i) {
      cor(data[i, x], data[i, y], use = "complete.obs", method = "spearman")
    },
    R = 2000
  )
  
  ci <- boot.ci(boot_correlation, type = "perc")
  
  t0   <- round(boot_correlation$t0, 2)
  low  <- round(ci$percent[4], 2)
  high <- round(ci$percent[5], 2)
  
  exp_cv[i, "r"] <- t0
  exp_cv[i, "lower"] <- low
  exp_cv[i, "upper"] <- high
}

eff_cv <- data.frame(GEEE_var = rep(c("improvement", "worsening", "side_effect"), times = 2),
                     Ext_var = rep(c("NRS", "PDI"), each = 3),
                     r = NA,
                     lower = NA,
                     upper = NA)

set.seed(1656)

for(i in 1:nrow(eff_cv)){
  d <- personmeans_eff
  x <- eff_cv[i, "GEEE_var"]
  y <- eff_cv[i, "Ext_var"]
  
  boot_correlation <- boot(
    data = d,
    statistic = function(data, i) {
      cor(data[i, x], data[i, y], use = "complete.obs", method = "spearman")
    },
    R = 2000
  )
  
  ci <- boot.ci(boot_correlation, type = "perc")
  
  t0   <- round(boot_correlation$t0, 2)
  low  <- round(ci$percent[4], 2)
  high <- round(ci$percent[5], 2)
  
  eff_cv[i, "r"] <- t0
  eff_cv[i, "lower"] <- low
  eff_cv[i, "upper"] <- high
}

# Sensitivity analysis

# The analyses above were repeated using alternative
# operationalizations of the external validity measures:
#
# Primary analysis             Sensitivity analysis
# NRS mean                     NRS maximum
# PDI sum score                PDI global rating
#
# The analytical procedure, including Spearman correlations,
# 2,000 bootstrap resamples, and percentile confidence
# intervals, was otherwise identical.

# ============================================================
# 3. Hurdle models
# ============================================================

# 0. Data preparation, smaller data frames with relevant vars only
# Note: Assessment times were calculated precisely in months for the analyses 
# but are reported as monthly time points in the manuscript for easier 
# interpretation.

exp_imp_df <- geee_val_data %>%
  dplyr::select(VPN_Code_corr, contains("geetreat_exp_01"))
exp_imp_df <- exp_imp_df %>%
  pivot_longer(cols = -VPN_Code_corr,
               names_to = "time",
               values_to = "improvement",
               names_pattern = "Geetreat_exp_01_rec_T(.*)") %>%
  mutate(time = as.numeric(time),
         time_months = case_when(time == 0 ~ 0,
                                 time == 1 ~ 1.384615,
                                 time == 2 ~ 2.307692,
                                 time == 3 ~ 3.230769, 
                                 time == 4 ~ 4.153846, 
                                 time == 5 ~ 5.076923,
                                 time == 6 ~ 6))

exp_worse_df <- geee_val_data %>%
  dplyr::select(VPN_Code_corr, contains("geetreat_exp_02"))
exp_worse_df <- exp_worse_df %>%
  pivot_longer(cols = -VPN_Code_corr,
               names_to = "time",
               values_to = "worsening",
               names_pattern = "Geetreat_exp_02_rec_T(.*)") %>%
  mutate(time = as.numeric(time),
         time_months = case_when(time == 0 ~ 0,
                                 time == 1 ~ 1.384615,
                                 time == 2 ~ 2.307692,
                                 time == 3 ~ 3.230769, 
                                 time == 4 ~ 4.153846, 
                                 time == 5 ~ 5.076923,
                                 time == 6 ~ 6))

exp_side_df <- geee_val_data %>% 
  dplyr::select(VPN_Code_corr, contains("geetreat_exp_03"))
exp_side_df <- exp_side_df %>% 
  pivot_longer(cols = -VPN_Code_corr,
               names_to = "time",
               values_to = "side",
               names_pattern = "Geetreat_exp_03_rec_T(.*)") %>% 
  mutate(time = as.numeric(time),
         time_months = case_when(time == 0 ~ 0,
                                 time == 1 ~ 1.384615,
                                 time == 2 ~ 2.307692,
                                 time == 3 ~ 3.230769, 
                                 time == 4 ~ 4.153846, 
                                 time == 5 ~ 5.076923,
                                 time == 6 ~ 6))

eff_imp_df <- geee_val_data %>%
  dplyr::select(VPN_Code_corr, contains("Gee_aktuell_01"))
eff_imp_df <- eff_imp_df %>%
  pivot_longer(cols = -VPN_Code_corr,
               names_to = "time",
               values_to = "improvement",
               names_pattern = "Gee_aktuell_01_T(.*)") %>%
  mutate(time = as.numeric(time),
         time_months = case_when(time == 0 ~ 0,
                                 time == 1 ~ 1.384615,
                                 time == 2 ~ 2.307692,
                                 time == 3 ~ 3.230769, 
                                 time == 4 ~ 4.153846, 
                                 time == 5 ~ 5.076923,
                                 time == 6 ~ 6))

eff_worse_df <- geee_val_data %>% 
  dplyr::select(VPN_Code_corr, contains("Gee_aktuell_02"))
eff_worse_df <- eff_worse_df %>% 
  pivot_longer(cols = -VPN_Code_corr,
               names_to = "time",
               values_to = "worsening",
               names_pattern = "Gee_aktuell_02_T(.*)") %>% 
  mutate(time = as.numeric(time),
         time_months = case_when(time == 0 ~ 0,
                                 time == 1 ~ 1.384615,
                                 time == 2 ~ 2.307692,
                                 time == 3 ~ 3.230769, 
                                 time == 4 ~ 4.153846, 
                                 time == 5 ~ 5.076923,
                                 time == 6 ~ 6))

eff_side_df <- geee_val_data %>% 
  dplyr::select(VPN_Code_corr, contains("Gee_aktuell_03"))
eff_side_df <- eff_side_df %>% 
  pivot_longer(cols = -VPN_Code_corr,
               names_to = "time",
               values_to = "side",
               names_pattern = "Gee_aktuell_03_T(.*)") %>% 
  mutate(time = as.numeric(time),
         time_months = case_when(time == 0 ~ 0,
                                 time == 1 ~ 1.384615,
                                 time == 2 ~ 2.307692,
                                 time == 3 ~ 3.230769, 
                                 time == 4 ~ 4.153846, 
                                 time == 5 ~ 5.076923,
                                 time == 6 ~ 6))

# 1. Model Specification
# 1.1 Prior Definition
pr_hurdle_exp_imp <- c(
  # hurdle component
  set_prior("normal(0.0,0.75)", class = "b", dpar = "hu"),
  set_prior("normal(-3, 1.0)", class = "Intercept", dpar = "hu"),
  
  # lognormal component
  set_prior("normal(0.0,0.3)", class = "b", dpar = "mu"),
  set_prior("normal(1.9, 0.2)", class = "Intercept", dpar = "mu"),
  set_prior("normal(0, 0.5)", class = "sigma", lb = 0),
  
  # random effects
  set_prior("exponential(2)", class = "sd", dpar = "hu"),
  set_prior("exponential(2)", class = "sd", dpar = "mu")
)

pr_hurdle_exp_worse <- c(
  # hurdle component
  set_prior("normal(0.0,0.75)", class = "b", dpar = "hu"),
  set_prior("normal(0.0, 1.0)", class = "Intercept", dpar = "hu"),
  
  # lognormal component
  set_prior("normal(0.0,0.3)", class = "b", dpar = "mu"),
  set_prior("normal(.9,.4)", class = "Intercept", dpar = "mu"),
  set_prior("normal(0, 0.5)", class = "sigma", lb = 0),
  
  # random effects
  set_prior("exponential(2)", class = "sd", dpar = "hu"),
  set_prior("exponential(2)", class = "sd", dpar = "mu")
)

pr_hurdle_exp_side <- c(
  # hurdle component
  set_prior("normal(0.0,0.75)", class = "b", dpar = "hu"),
  set_prior("normal(-2, 1.0)", class = "Intercept", dpar = "hu"),
  
  # lognormal component
  set_prior("normal(0.0,0.3)", class = "b", dpar = "mu"),
  set_prior("normal(1.4, 0.5", class = "Intercept", dpar = "mu"),
  set_prior("normal(0, 0.5)", class = "sigma", lb = 0),
  
  # random effects
  set_prior("exponential(2)", class = "sd", dpar = "hu"),
  set_prior("exponential(2)", class = "sd", dpar = "mu")
)

pr_hurdle_eff_imp <- c(
  # hurdle component
  set_prior("normal(0.0,0.75)", class = "b", dpar = "hu"),
  set_prior("normal(-2, 1.0)", class = "Intercept", dpar = "hu"),
  
  # lognormal component
  set_prior("normal(0.0,0.3)", class = "b", dpar = "mu"),
  set_prior("normal(1.8, 0.3)", class = "Intercept", dpar = "mu"),
  set_prior("normal(0, 0.5)", class = "sigma", lb = 0),
  
  # random effects
  set_prior("exponential(2)", class = "sd", dpar = "hu"),
  set_prior("exponential(2)", class = "sd", dpar = "mu")
)

pr_hurdle_eff_worse <- c(
  # hurdle component
  set_prior("normal(0.0,0.75)", class = "b", dpar = "hu"),
  set_prior("normal(0.3, 1.0)", class = "Intercept", dpar = "hu"),
  
  # lognormal component
  set_prior("normal(0.0,0.3)", class = "b", dpar = "mu"),
  set_prior("normal(1.2, 0.5)", class = "Intercept", dpar = "mu"),
  set_prior("normal(0, 0.5)", class = "sigma", lb = 0),
  
  # random effects
  set_prior("exponential(2)", class = "sd", dpar = "hu"),
  set_prior("exponential(2)", class = "sd", dpar = "mu")
)

pr_hurdle_eff_side <- c(
  # hurdle component
  set_prior("normal(0.0,0.75)", class = "b", dpar = "hu"),
  set_prior("normal(-2, 1.0)", class = "Intercept", dpar = "hu"),
  
  # lognormal component
  set_prior("normal(0.0,0.3)", class = "b", dpar = "mu"),
  set_prior("normal(1.5, 0.4)", class = "Intercept", dpar = "mu"),
  set_prior("normal(0, 0.5)", class = "sigma", lb = 0),
  
  # random effects
  set_prior("exponential(2)", class = "sd", dpar = "hu"),
  set_prior("exponential(2)", class = "sd", dpar = "mu")
)

# 1.2 Models fitted only on prior distribution
prior_exp_imp <- brm(
  data = exp_imp_df,
  formula = bf(improvement ~ time_months + (1 | VPN_Code_corr),
               hu ~ time_months + (1 | VPN_Code_corr)),
  family = hurdle_lognormal(link = 'identity', link_sigma = 'log', link_hu = 'logit'), # default settings
  prior = pr_hurdle_exp_imp,
  sample_prior = "only",
  chains = 4, cores = 8, iter = 6000, warmup = 1500, seed = 333,
  control = list(adapt_delta = 0.97, max_treedepth = 15),
  backend = "cmdstanr"
)

prior_exp_worse <- brm(
  data = exp_worse_df,
  formula = bf(worsening ~ time_months + (1 | VPN_Code_corr),
               hu ~ time_months + (1 | VPN_Code_corr)),
  family = hurdle_lognormal(link = 'identity', link_sigma = 'log', link_hu = 'logit'), # default settings
  prior = pr_hurdle_exp_worse,
  sample_prior = "only",
  chains = 4, cores = 8, iter = 6000, warmup = 1500, seed = 333,
  control = list(adapt_delta = 0.97, max_treedepth = 15),
  backend = "cmdstanr"
)

prior_exp_side <- brm(
  data = exp_side_df,
  formula = bf(side ~ time_months + (1 | VPN_Code_corr),
               hu ~ time_months + (1 | VPN_Code_corr)),
  family = hurdle_lognormal(link = 'identity', link_sigma = 'log', link_hu = 'logit'), # default settings
  prior = pr_hurdle_exp_side,
  sample_prior = "only",
  chains = 4, cores = 8, iter = 6000, warmup = 1500, seed = 333,
  control = list(adapt_delta = 0.97, max_treedepth = 15),
  backend = "cmdstanr"
)

prior_eff_imp <- brm(
  data = eff_imp_df,
  formula = bf(improvement ~ time_months + (1 | VPN_Code_corr),
               hu ~ time_months + (1 | VPN_Code_corr)),
  family = hurdle_lognormal(link = 'identity', link_sigma = 'log', link_hu = 'logit'), # default settings
  prior = pr_hurdle_eff_imp,
  sample_prior = "only",
  chains = 4, cores = 8, iter = 6000, warmup = 1500, seed = 333,
  control = list(adapt_delta = 0.97, max_treedepth = 15),
  backend = "cmdstanr"
)

prior_eff_worse <- brm(
  data = eff_worse_df,
  formula = bf(worsening ~ time_months + (1 | VPN_Code_corr),
               hu ~ time_months + (1 | VPN_Code_corr)),
  family = hurdle_lognormal(link = 'identity', link_sigma = 'log', link_hu = 'logit'), # default settings
  prior = pr_hurdle_eff_worse,
  sample_prior = "only",
  chains = 4, cores = 8, iter = 6000, warmup = 1500, seed = 333,
  control = list(adapt_delta = 0.97, max_treedepth = 15),
  backend = "cmdstanr"
)

prior_eff_side <- brm(
  data = eff_side_df,
  formula = bf(side ~ time_months + (1 | VPN_Code_corr),
               hu ~ time_months + (1 | VPN_Code_corr)),
  family = hurdle_lognormal(link = 'identity', link_sigma = 'log', link_hu = 'logit'), # default settings
  prior = pr_hurdle_eff_side,
  sample_prior = "only",
  chains = 4, cores = 8, iter = 6000, warmup = 1500, seed = 333,
  control = list(adapt_delta = 0.97, max_treedepth = 15),
  backend = "cmdstanr"
)

# 1.3 Prior Checks
all_objects <- ls()
prior_models <- all_objects[grepl("prior_", all_objects)]

plots_prior_s1 <- list() # full distribution
plots_prior_s2 <-  list() # hurdle
tables_prior <- list()

time_labs <- c(
  "0" = "before laparoscopy",
  "1.384615" = "1 month after laparoscopy",
  "2.307692" = "2 months after laparoscopy",
  "3.230769" = "3 months after laparoscopy",
  "4.153846" = "4 months after laparoscopy",
  "5.076923" = "5 months after laparoscopy",
  "6" = "6 months after laparoscopy"
)

lower_bound <- 0
upper_bound <- 10

theme_prior <- theme_minimal() +
  theme(
    text = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    plot.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  )

for(i in prior_models){
  fit <- get(i)
  
  subscale <- ifelse(grepl("exp", i), "exp", "eff")
  item <- ifelse(grepl("imp", i), "imp", ifelse(grepl("worse", i), "worse", "side"))
  
  time <- fit$data$time_months
  
  # Full Distribution
  yrep <- posterior_predict(fit)
  
  df_prior_s1 <- as.data.frame(yrep) %>%
    mutate(draw = row_number()) %>%
    pivot_longer(-draw, names_to = "obs", values_to = "yrep") %>%
    mutate(time = time[as.numeric(gsub("V", "", obs))])
  
  title <- if(subscale == "exp"){
    if(item == "imp"){
      "A: Model 1"
    }else if(item == "worse"){
      "B: Model 2"
    }else{
      "C: Model 3"
    }
  }else{
    if(item == "imp"){
      "D: Model 4"
    }else if(item == "worse"){
      "E: Model 5"
    }else{
      "F: Model 6"
    }
  }
  
  outcome <- if(subscale == "exp"){
    if(item == "imp"){
      "expected improvement"
    }else if(item == "worse"){
      "expected worsening"
    }else{
      "expected side effects"
    }
  }else{
    if(item == "imp"){
      "perceived improvement"
    }else if(item == "worse"){
      "perceived worsening"
    }else{
      "perceived side effects"
    }
  }
  
  p_prior_s1 <- ggplot(df_prior_s1, aes(x = yrep)) +
    
    geom_histogram(
      breaks = 0:10,
      fill = "grey",
      color = "black",
      closed = "left"
    ) +
    
    facet_wrap(
      ~time,
      labeller = labeller(time = time_labs),
      nrow = 1
    ) +
    
    scale_x_continuous(breaks = 0:10) +
    
    labs(
      title = title, #Prior predictive distribution of outcome (including zeros)
      x = outcome,
    ) +
    
    theme_prior +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  plots_prior_s1[[i]] <- p_prior_s1
  
  # Hurdle Component
  hu_epred <- posterior_epred(fit, dpar = "hu")
  
  df_hu <- as.data.frame(hu_epred) %>%
    mutate(draw = row_number()) %>%
    pivot_longer(-draw, names_to = "obs", values_to = "p_zero") %>%
    mutate(time = time[as.numeric(gsub("V", "", obs))])
  
  df_hu_summary <- df_hu %>%
    group_by(time) %>%
    summarise(
      median = median(p_zero),
      lower = quantile(p_zero, 0.025),
      upper = quantile(p_zero, 0.975),
      .groups = "drop"
    )
  
  p_prior_s2 <- ggplot(df_hu_summary, aes(x = time)) +
    
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70", alpha = 0.5) +
    
    geom_line(aes(y = median), size = 1) +
    
    geom_point(aes(y = median), size = 2) +
    
    scale_x_continuous(breaks = 0:10, labels = 0:10) +
    
    labs(
      title = title, #"Prior predictive trajectory of P(Y = 0)"
      x = "months since laparoscopy",
      y = paste0("P(Y = 0) \n", outcome)) +
    
    
    theme_prior +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
    
    plots_prior_s2[[i]] <- p_prior_s2
  
  
  df_table <- df_prior_s1
  
  table_prior <- df_table %>%
    group_by(time) %>%
    summarise(
      
      model = i,
      
      prop_zero = round(mean(yrep == 0),digits = 2),
      prop_in_scale = round(mean(yrep >= lower_bound & yrep <= upper_bound), digits = 2),
      median_y = round(median(yrep), digits = 2),
      p90 = round(quantile(yrep, 0.90), digits = 2),
      p95 = round(quantile(yrep, 0.95), digits = 2),
      p99 = round(quantile(yrep, 0.99), digits = 2)
    )
  
  tables_prior[[i]] <- flextable(table_prior)
}

# 2. Fit Models

fit_exp_imp_lognorm_prior <- brm(
  data = exp_imp_df,
  formula = bf(improvement ~ time_months + (1 | VPN_Code_corr),
               hu ~ time_months + (1 | VPN_Code_corr)),
  family = hurdle_lognormal(link = 'identity', link_sigma = 'log', link_hu = 'logit'), # default settings
  prior = pr_hurdle_exp_imp,
  chains = 4, cores = 8, iter = 6000, warmup = 1500, seed = 333,
  control = list(adapt_delta = 0.97, max_treedepth = 15),
  backend = "cmdstanr"
)

fit_exp_worse_lognorm_prior <- brm(
  data = exp_worse_df,
  formula = bf(worsening ~ time_months + (1 | VPN_Code_corr),
               hu ~ time_months + (1 | VPN_Code_corr)),
  family = hurdle_lognormal(link = 'identity', link_sigma = 'log', link_hu = 'logit'),
  prior = pr_hurdle_exp_worse,
  chains = 4, cores = 8, iter = 6000, warmup = 1500, seed = 333,
  control = list(adapt_delta = 0.97, max_treedepth = 15),
  backend = "cmdstanr"
)

fit_exp_side_lognorm_prior <- brm(
  data = exp_side_df,
  formula = bf(side ~ time_months + (1 | VPN_Code_corr),
               hu ~ time_months + (1 | VPN_Code_corr)),
  family = hurdle_lognormal(link = 'identity', link_sigma = 'log', link_hu = 'logit'),
  prior = pr_hurdle_exp_side,
  chains = 4, cores = 8, iter = 6000, warmup = 1500, seed = 333,
  control = list(adapt_delta = 0.97, max_treedepth = 15),
  backend = "cmdstanr"
)

fit_eff_imp_lognorm_prior <- brm(
  data = eff_imp_df,
  formula = bf(improvement ~ time_months + (1 | VPN_Code_corr),
               hu ~ time_months + (1 | VPN_Code_corr)),
  family = hurdle_lognormal(link = 'identity', link_sigma = 'log', link_hu = 'logit'),
  prior = pr_hurdle_eff_imp,
  chains = 4, cores = 8, iter = 6000, warmup = 1500, seed = 333,
  control = list(adapt_delta = 0.97, max_treedepth = 15),
  backend = "cmdstanr"
)

fit_eff_worse_lognorm_prior <- brm(
  data = eff_worse_df,
  formula = bf(worsening ~ time_months + (1 | VPN_Code_corr),
               hu ~ time_months + (1 | VPN_Code_corr)),
  family = hurdle_lognormal(link = 'identity', link_sigma = 'log', link_hu = 'logit'),
  prior = pr_hurdle_eff_worse,
  chains = 4, cores = 8, iter = 6000, warmup = 1500, seed = 333,
  control = list(adapt_delta = 0.97, max_treedepth = 15),
  backend = "cmdstanr"
)

fit_eff_side_lognorm_prior <- brm(
  data = eff_side_df,
  formula = bf(side ~ time_months + (1 | VPN_Code_corr),
               hu ~ time_months + (1 | VPN_Code_corr)),
  family = hurdle_lognormal(link = 'identity', link_sigma = 'log', link_hu = 'logit'),
  prior = pr_hurdle_eff_side,
  chains = 4, cores = 8, iter = 6000, warmup = 1500, seed = 333,
  control = list(adapt_delta = 0.97, max_treedepth = 15),
  backend = "cmdstanr"
)

# 3. Posterior predictive checks
# 3.1. ECDF-overlay

post_exp_imp_ecdf <- 
  pp_check(fit_exp_imp_lognorm_prior, type = "ecdf_overlay_grouped", group = "time_months", ndraws = 50)+
  coord_cartesian(xlim = c(0, 11))+
  facet_wrap(
    ~group,
    labeller = labeller(group = time_labs))+
  theme_prior

post_exp_worse_ecdf <- 
  pp_check(fit_exp_worse_lognorm_prior, type = "ecdf_overlay_grouped", group = "time_months", ndraws = 50)+
  coord_cartesian(xlim = c(0, 11))+
  facet_wrap(
    ~group,
    labeller = labeller(group = time_labs))+
  theme_prior

post_exp_side_ecdf <- 
  pp_check(fit_exp_side_lognorm_prior, type = "ecdf_overlay_grouped", group = "time_months", ndraws = 50)+
  coord_cartesian(xlim = c(0, 11))+
  facet_wrap(
    ~group,
    labeller = labeller(group = time_labs))+
  theme_prior

post_eff_imp_ecdf <- 
  pp_check(fit_eff_imp_lognorm_prior, type = "ecdf_overlay_grouped", group = "time_months", ndraws = 50)+
  coord_cartesian(xlim = c(0, 11))+
  facet_wrap(
    ~group,
    labeller = labeller(group = time_labs))+
  theme_prior

post_eff_worse_ecdf <- 
  pp_check(fit_eff_worse_lognorm_prior, type = "ecdf_overlay_grouped", group = "time_months", ndraws = 50)+
  coord_cartesian(xlim = c(0, 11))+
  facet_wrap(
    ~group,
    labeller = labeller(group = time_labs))+
  theme_prior

post_eff_side_ecdf <- 
  pp_check(fit_eff_side_lognorm_prior, type = "ecdf_overlay_grouped", group = "time_months", ndraws = 50)+
  coord_cartesian(xlim = c(0, 11))+
  facet_wrap(
    ~group,
    labeller = labeller(group = time_labs))+
  theme_prior

# 3.2. Reproduction of Floor / Ceiling effects 

models <- list(
  exp_imp   = fit_exp_imp_lognorm_prior,
  exp_worse = fit_exp_worse_lognorm_prior,
  exp_side  = fit_exp_side_lognorm_prior,
  eff_imp   = fit_eff_imp_lognorm_prior,
  eff_worse = fit_eff_worse_lognorm_prior,
  eff_side  = fit_eff_side_lognorm_prior
)

# function for one model
get_ppc_summary <- function(model, ndraws = 1000) {
  
  # Posterior predictive draws
  yrep <- posterior_predict(model, draws = ndraws)
  
  # observed outcome values
  y <- model.response(model.frame(model))
  
  # statistic for each replicate
  zero_rep    <- apply(yrep, 1, function(x) mean(x == 0))
  ceiling_rep <- apply(yrep, 1, function(x) mean(x >= 10))
  out_rep     <- apply(yrep, 1, function(x) mean(x > 10))
  mean_rep    <- apply(yrep, 1, mean)
  sd_rep      <- apply(yrep, 1, sd)
  
  data.frame(
    Zero_observed = mean(y == 0),
    Zero_median   = median(zero_rep),
    Zero_lower    = quantile(zero_rep, .025),
    Zero_upper    = quantile(zero_rep, .975),
    
    Ceiling_observed = mean(y >= 10),
    Ceiling_median   = median(ceiling_rep),
    Ceiling_lower    = quantile(ceiling_rep, .025),
    Ceiling_upper    = quantile(ceiling_rep, .975),
    
    Mean_observed = mean(y),
    Mean_median   = median(mean_rep),
    Mean_lower    = quantile(mean_rep, .025),
    Mean_upper    = quantile(mean_rep, .975),
    
    SD_observed = sd(y),
    SD_median   = median(sd_rep),
    SD_lower    = quantile(sd_rep, .025),
    SD_upper    = quantile(sd_rep, .975)
  )
}

# all models
set.seed(2113)

ppc_summary <- do.call(
  rbind,
  lapply(models, get_ppc_summary)
)

# add model name to col
ppc_summary <- cbind(
  Model = rownames(ppc_summary),
  ppc_summary
)

rownames(ppc_summary) <- NULL

# round values
ppc_summary[, -1] <- round(ppc_summary[, -1], 3)

ppc_summary

# 3.3. Model results

models <- c(
  "fit_exp_imp_lognorm_prior",
  "fit_exp_worse_lognorm_prior",
  "fit_exp_side_lognorm_prior",
  "fit_eff_imp_lognorm_prior",
  "fit_eff_worse_lognorm_prior",
  "fit_eff_side_lognorm_prior"
)

hurdle_results <- data.frame(
  model = rep(models, each = 4),
  parameter = rep(c("b_hu_Intercept", "b_hu_time_months", "b_Intercept", "b_time_months"), times = 6),
  estimate = NA_real_,
  CrI = NA_character_,
  pd = NA_real_,
  OR_mult = NA_real_,
  stringsAsFactors = FALSE
)


for (i in models) {
  
  model <- get(i)
  
  model_results <- describe_posterior(
    model,
    effects = "fixed",
    component = "all",
    ci = 0.95
  ) %>%
    mutate(
      pd_std = ifelse(pd > 1, pd / 100, pd),
      estimate = Median,
      CrI = paste0("[", CI_low, ", ", CI_high, "]"),
      OR_mult = exp(Median)
    ) %>%
    select(Parameter, estimate, CrI, pd_std, OR_mult)
  
  rows_to_fill <- hurdle_results$model == i
  
  hurdle_results[rows_to_fill, "estimate"] <- model_results$estimate[
    match(hurdle_results$parameter[rows_to_fill], model_results$Parameter)
  ]
  
  hurdle_results[rows_to_fill, "CrI"] <- model_results$CrI[
    match(hurdle_results$parameter[rows_to_fill], model_results$Parameter)
  ]
  
  hurdle_results[rows_to_fill, "pd"] <- model_results$pd_std[
    match(hurdle_results$parameter[rows_to_fill], model_results$Parameter)
  ]
  
  hurdle_results[rows_to_fill, "OR_mult"] <- model_results$OR_mult[
    match(hurdle_results$parameter[rows_to_fill], model_results$Parameter)
  ]
}

hurdle_results

# 3.4. Descriptive statistics posterior distribution

posterior_ppc_tables <- list()

for(model in models){
  
  fit <- get(model)
  
  tab <- posterior_summary(
    fit,
    probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
    variable = c("^b_", "^sd_", "^sigma"),
    regex = TRUE
  )
  
  tab <- round(tab, 2)
  
  tab <- as.data.frame(tab) |>
    rownames_to_column(var = "Parameter")
  
  posterior_ppc_tables[[model]] <- flextable(tab)
}

# 3.5. Predicted trajectories

newdat <- tibble(
  time_months = seq(min(exp_imp_df$time_months),
                    max(exp_imp_df$time_months),
                    length.out = 100)
)

theme_hurdle <- 
  theme_classic()+
  theme(
    text = element_text(size = 10),          # default text
    axis.text = element_text(size = 10),     # axis numbers
    axis.title = element_text(size = 10),    # axis titles
    plot.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.position = "bottom",
    legend.justification = "center"
  )

# 3.5.1 Expectations

# 3.5.1.1 Stack predictions from multiple models

model_list_exp <- list(
  improvement = fit_exp_imp_lognorm_prior,
  worsening = fit_exp_worse_lognorm_prior,
  side_effect = fit_exp_side_lognorm_prior
)

all_linpred <- bind_rows(lapply(names(model_list_exp), function(mn) {
  newdat %>%
    add_linpred_draws(model_list_exp[[mn]], dpar = c("mu", "sigma", "hu"), re_formula = NA) %>%
    mutate(
      prob_zero = plogis(hu),
      prob_positive = 1 - prob_zero,
      mean_positive = exp(mu + 0.5 * sigma^2),
      item = mn
    )
}))

all_epred <- bind_rows(lapply(names(model_list_exp), function(mn) {
  newdat %>%
    add_epred_draws(model_list_exp[[mn]], re_formula = NA) %>%
    rename(overall_mean = .epred) %>%
    mutate(item = mn)
}))

point_times_exp <- tibble(time_months = 0:6)
points_hurdle_exp_df <- bind_rows(lapply(names(model_list_exp), function(mn) {
  
  point_times_exp %>%
    add_linpred_draws(
      model_list_exp[[mn]],
      dpar = c("mu","sigma","hu"),
      re_formula = NA
    ) %>%
    mutate(
      prob_zero = plogis(hu),
      prob_positive = 1 - prob_zero,
      item = mn
    ) %>%
    group_by(item, time_months) %>%
    summarise(
      prob_positive = mean(prob_positive),
      .groups = "drop"
    )
}))

points_conditional_exp_df <- bind_rows(lapply(names(model_list_exp), function(mn) {
  
  point_times_exp %>%
    add_linpred_draws(
      model_list_exp[[mn]],
      dpar = c("mu", "sigma"),
      re_formula = NA
    ) %>%
    mutate(
      mean_positive = exp(mu + 0.5 * sigma^2),
      item = mn
    ) %>%
    group_by(item, time_months) %>%
    summarise(
      mean_positive = mean(mean_positive),
      .groups = "drop"
    )
}))

points_overall_exp_df <- bind_rows(lapply(names(model_list_exp), function(mn) {
  
  point_times_exp %>%
    add_epred_draws(
      model_list_exp[[mn]],
      re_formula = NA
    ) %>%
    rename(overall_mean = .epred) %>%
    mutate(item = mn) %>%
    group_by(item, time_months) %>%
    summarise(
      overall_mean = mean(overall_mean),
      .groups = "drop"
    )
  
}))


# 3.5.1.2 Hurdle Component

p_hurdle_exp <- all_linpred %>%
  ggplot(aes(x = time_months, y = prob_positive * 100)) +  # multiply by 100 for %
  stat_lineribbon(aes(group = item),
                  fill = "grey50",
                  .width = 0.95, alpha = 0.2, show.legend = F) +
  geom_point(
    data = points_hurdle_exp_df,
    aes(y = prob_positive * 100, shape = item),
    size = 2
  ) +
  labs(y = "P(Y > 0) [%]", x = "months since laparoscopy", title = "(a) Hurdle component") +
  scale_shape_manual(
    values = c(
      "improvement" = 15,
      "worsening" = 16,
      "side_effect" = 17
    ),
    breaks = c("improvement", "worsening", "side_effect"),
    labels = c("Improvement", "Worsening", "Side Effects"),
    name = NULL,
    guide = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  scale_x_continuous(limits = c(0,6), breaks = seq(0, 6, 1)) +
  theme_hurdle

# 3.5.1.3 Conditional expectation (0–10)

p_conditional_exp <- all_linpred %>%
  ggplot(aes(x = time_months, y = mean_positive)) +
  stat_lineribbon(aes(group = item),
                  fill = "grey50",
                  .width = 0.95, alpha = 0.2, show.legend = F) +
  geom_point(
    data = points_conditional_exp_df,
    aes(y = mean_positive, shape = item),
    size = 2
  ) +
  scale_shape_manual(
    values = c(
      "improvement" = 15,
      "worsening" = 16,
      "side_effect" = 17
    ),
    breaks = c("improvement", "worsening", "side_effect"),
    labels = c("Improvement", "Worsening", "Side Effects"),
    name = NULL,
    guide = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  labs(y = "E[Y | Y > 0]", x = "months since laparoscopy", title = "(b) Positive outcome distribution") +
  scale_y_continuous(limits = c(0, 10), breaks = 0:10) +
  scale_x_continuous(limits = c(0,6), breaks = seq(0, 6, 1)) +
  theme_hurdle

# 3.5.1.4 Overall expectation (0–10)

p_overall_exp <- all_epred %>%
  ggplot(aes(x = time_months, y = overall_mean)) +
  stat_lineribbon(aes(group = item),
                  fill = "grey50",
                  .width = 0.95, alpha = 0.2, show.legend = F) +
  geom_point(
    data = points_overall_exp_df,
    aes(y = overall_mean, shape = item),
    size = 2
  ) +
  scale_shape_manual(
    values = c(
      "improvement" = 15,
      "worsening" = 16,
      "side_effect" = 17
    ),
    breaks = c("improvement", "worsening", "side_effect"),
    labels = c("Improvement", "Worsening", "Side Effects"),
    name = NULL,
    guide = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  labs(y = "Marginal Expected Y", x = "months since laparoscopy", title = "(c) Overall model prediction") +
  scale_y_continuous(limits = c(0, 10), breaks = 0:10) +
  scale_x_continuous(limits = c(0,6), breaks = seq(0, 6, 1)) +
  theme_hurdle

# 3.5.2 Effects

newdat <- tibble(
  time_months = seq(1,6,
                    length.out = 100)
)

model_list_eff <- list(
  improvement = fit_eff_imp_lognorm_prior,
  worsening = fit_eff_worse_lognorm_prior,
  side_effect = fit_eff_side_lognorm_prior
)

# 3.5.2.1 Stack predictions from multiple models

all_linpred <- bind_rows(lapply(names(model_list_eff), function(mn) {
  newdat %>%
    add_linpred_draws(model_list_eff[[mn]], dpar = c("mu", "sigma", "hu"), re_formula = NA) %>%
    mutate(
      prob_zero = plogis(hu),
      prob_positive = 1 - prob_zero,
      mean_positive = exp(mu + 0.5 * sigma^2),
      item = mn
    )
}))

all_epred <- bind_rows(lapply(names(model_list_eff), function(mn) {
  newdat %>%
    add_epred_draws(model_list_eff[[mn]], re_formula = NA) %>%
    rename(overall_mean = .epred) %>%
    mutate(item = mn)
}))

point_times_eff <- tibble(time_months = 1:6)
points_hurdle_eff_df <- bind_rows(lapply(names(model_list_eff), function(mn) {
  
  point_times_eff %>%
    add_linpred_draws(
      model_list_eff[[mn]],
      dpar = c("mu","sigma","hu"),
      re_formula = NA
    ) %>%
    mutate(
      prob_zero = plogis(hu),
      prob_positive = 1 - prob_zero,
      item = mn
    ) %>%
    group_by(item, time_months) %>%
    summarise(
      prob_positive = mean(prob_positive),
      .groups = "drop"
    )
}))

points_conditional_eff_df <- bind_rows(lapply(names(model_list_eff), function(mn) {
  
  point_times_eff %>%
    add_linpred_draws(
      model_list_eff[[mn]],
      dpar = c("mu", "sigma"),
      re_formula = NA
    ) %>%
    mutate(
      mean_positive = exp(mu + 0.5 * sigma^2),
      item = mn
    ) %>%
    group_by(item, time_months) %>%
    summarise(
      mean_positive = mean(mean_positive),
      .groups = "drop"
    )
}))

points_overall_eff_df <- bind_rows(lapply(names(model_list_eff), function(mn) {
  
  point_times_eff %>%
    add_epred_draws(
      model_list_eff[[mn]],
      re_formula = NA
    ) %>%
    rename(overall_mean = .epred) %>%
    mutate(item = mn) %>%
    group_by(item, time_months) %>%
    summarise(
      overall_mean = mean(overall_mean),
      .groups = "drop"
    )
  
}))


# 3.5.2.2 Hurdle component (%)

p_hurdle_eff <- all_linpred %>%
  ggplot(aes(x = time_months, y = prob_positive * 100)) +  # multiply by 100 for %
  stat_lineribbon(aes(group = item),
                  fill = "grey50",
                  .width = 0.95, alpha = 0.2, show.legend = F) +
  geom_point(
    data = points_hurdle_eff_df,
    aes(y = prob_positive * 100, shape = item),
    size = 2
  ) +
  scale_shape_manual(
    values = c(
      "improvement" = 15,
      "worsening" = 16,
      "side_effect" = 17
    ),
    breaks = c("improvement", "worsening", "side_effect"),
    labels = c("Improvement", "Worsening", "Side Effects"),
    name = NULL,
    guide = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  labs(y = "P(Y > 0) [%]", x = "months since laparoscopy", title = "(a) Hurdle component") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  scale_x_continuous(limits = c(0,6), breaks = seq(0, 6, 1)) +
  theme_hurdle

# 3.5.2.3 Conditional effects (0–10)

p_conditional_eff <- all_linpred %>%
  ggplot(aes(x = time_months, y = mean_positive)) +
  stat_lineribbon(aes(group = item),
                  fill = "grey50",
                  .width = 0.95, alpha = 0.2, show.legend = F) +
  geom_point(
    data = points_conditional_eff_df,
    aes(y = mean_positive, shape = item),
    size = 2
  ) +
  scale_shape_manual(
    values = c(
      "improvement" = 15,
      "worsening" = 16,
      "side_effect" = 17
    ),
    breaks = c("improvement", "worsening", "side_effect"),
    labels = c("Improvement", "Worsening", "Side Effects"),
    name = NULL,
    guide = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  labs(y = "E[Y | Y > 0]", x = "months since laparoscopy", title = "(b) Positive outcome distribution") +
  scale_y_continuous(limits = c(0, 10), breaks = 0:10) +
  scale_x_continuous(limits = c(0,6), breaks = seq(0, 6, 1)) +
  theme_hurdle


# 3.5.2.4 Overall expectation (0–10)

p_overall_eff <- all_epred %>%
  ggplot(aes(x = time_months, y = overall_mean)) +
  stat_lineribbon(aes(group = item),
                  fill = "grey50",
                  .width = 0.95, alpha = 0.2, show.legend = F) +
  geom_point(
    data = points_overall_eff_df,
    aes(y = overall_mean, shape = item),
    size = 2
  ) +
  scale_shape_manual(
    values = c(
      "improvement" = 15,
      "worsening" = 16,
      "side_effect" = 17
    ),
    breaks = c("improvement", "worsening", "side_effect"),
    labels = c("Improvement", "Worsening", "Side Effects"),
    name = NULL,
    guide = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  labs(y = "Marginal Expected Y", x = "months since laparoscopy", title = "(c) Overall model prediction") +
  scale_y_continuous(limits = c(0, 10), breaks = 0:10) +
  scale_x_continuous(limits = c(0,6), breaks = seq(0, 6, 1)) +
  theme_hurdle


# ============================================================
# 4. Further Analyses Conducted During Manuscript Revision
# ============================================================

# 1. Detailed assessment of participation across measurement occasions

df_mo <- data.frame(
  VP_Code = geee_val_data$VPN_Code_corr,
  T0 = rep(NA, 393),
  T1 = rep(NA, 393),
  T2 = rep(NA, 393),
  T3 = rep(NA, 393),
  T4 = rep(NA, 393),
  T5 = rep(NA, 393),
  T6 = rep(NA, 393)
)

match_ID <- match(df_mo$VP_Code, geee_val_data$VPN_Code_corr)

df_mo$T0 <- ifelse(
  is.na(geee_val_data$Geetreat_exp_01_rec_T0[match_ID]),
  0,
  1
)
df_mo$T1 <- ifelse(
  is.na(geee_val_data$Geetreat_exp_01_rec_T1[match_ID]),
  0,
  1
)
df_mo$T2 <- ifelse(
  is.na(geee_val_data$Geetreat_exp_01_rec_T2[match_ID]),
  0,
  1
)
df_mo$T3 <- ifelse(
  is.na(geee_val_data$Geetreat_exp_01_rec_T3[match_ID]),
  0,
  1
)
df_mo$T4 <- ifelse(
  is.na(geee_val_data$Geetreat_exp_01_rec_T4[match_ID]),
  0,
  1
)
df_mo$T5 <- ifelse(
  is.na(geee_val_data$Geetreat_exp_01_rec_T5[match_ID]),
  0,
  1
)
df_mo$T6 <- ifelse(
  is.na(geee_val_data$Geetreat_exp_01_rec_T6[match_ID]),
  0,
  1
)

timepoints <- c("T0", "T1", "T2", "T3", "T4", "T5", "T6")

df_mo$participation <- apply(df_mo[, timepoints], 1, function(x) {
  
  if (all(x == 1)) {
    return("Complete")
  }
    last_response <- max(which(x == 1), na.rm = TRUE)
  
  if (any(x[1:last_response] == 0)) {
    return("Intermittent")
  }
    return("Dropout")
})


# 2. log Regression, examining drop-out after T0 and predictor vars

df_mo <- df_mo %>% 
  mutate(MO = T0 + T1 + T2 + T3 + T4 + T5 + T6,
         drop_out = ifelse(MO == 1, 1, 0))

data_logreg <- geee_val_data %>%
  left_join(
    df_mo %>%
      select(VP_Code, drop_out),
    by = c("VPN_Code_corr" = "VP_Code")
  )

model_logreg <- glm(
  drop_out ~
    PDI_sum_T0 +
    PDI_exp_sum_T0 +
    Geetreat_exp_01_rec_T0 +
    Geetreat_exp_02_rec_T0 +
    Geetreat_exp_03_rec_T0 +
    SWEM__rec_sum_T0,
  data = data_logreg,
  family = binomial(link = "logit")
)

summary(model_logreg) #0 = no drop-out; 1 = drop-out

exp(cbind(
  OR = coef(model_logreg),
  confint(model_logreg)
))

check_collinearity(model_logreg)

p_hat <- predict(model_logreg, type = "response")
pred <- ifelse(p_hat >= .5, 1, 0)
y <- model.response(model.frame(model_logreg))
table(
  Observed = y,
  Predicted = pred
)

# 3. Sensitivity analysis for between-person correlations using 
#    time-point-centered values

# 3.1. Intercorrelations

# Center variables around the respective time-point mean

mlm_data_exp_adjusted <- mlm_data_exp_complete %>%
  group_by(time) %>%
  mutate(
    improvement_adj = improvement - mean(improvement, na.rm = TRUE),
    worsening_adj   = worsening - mean(worsening, na.rm = TRUE),
    side_effect_adj = side_effect - mean(side_effect, na.rm = TRUE)
  ) %>%
  ungroup()

mlm_data_eff_adjusted <- mlm_data_eff_complete %>%
  group_by(time) %>%
  mutate(
    improvement_adj = improvement - mean(improvement, na.rm = TRUE),
    worsening_adj   = worsening - mean(worsening, na.rm = TRUE),
    side_effect_adj = side_effect - mean(side_effect, na.rm = TRUE)
  ) %>%
  ungroup()

# Calculate person means based on time-point-centered values

personmeans_exp_adjusted <- mlm_data_exp_adjusted %>%
  group_by(VPN_Code_corr) %>%
  summarise(
    improvement = mean(improvement_adj, na.rm = TRUE),
    worsening   = mean(worsening_adj, na.rm = TRUE),
    side_effect = mean(side_effect_adj, na.rm = TRUE),
    .groups = "drop"
  )

personmeans_eff_adjusted <- mlm_data_eff_adjusted %>%
  group_by(VPN_Code_corr) %>%
  summarise(
    improvement = mean(improvement_adj, na.rm = TRUE),
    worsening   = mean(worsening_adj, na.rm = TRUE),
    side_effect = mean(side_effect_adj, na.rm = TRUE),
    .groups = "drop"
  )

# Bootstrap analyses were conducted identically to the primary analyses,
# using personmeans_exp_adjusted instead of the original person means.

# 3.2 Convergent Validity

# Time-point-centered values

data_exp_adjusted <- data_exp %>%
  group_by(time) %>%
  mutate(
    improvement_adj = improvement - mean(improvement, na.rm = TRUE),
    worsening_adj   = worsening - mean(worsening, na.rm = TRUE),
    side_effect_adj = side_effect - mean(side_effect, na.rm = TRUE),
    NRS_exp_imp_adj = NRS_exp_imp - mean(NRS_exp_imp, na.rm = TRUE),
    NRS_exp_worse_adj = NRS_exp_worse - mean(NRS_exp_worse, na.rm = TRUE),
    PDI_exp_adj = PDI_exp - mean(PDI_exp, na.rm = TRUE)
  ) %>%
  ungroup()

data_eff_adjusted <- data_eff %>%
  group_by(time) %>%
  mutate(
    improvement_adj = improvement - mean(improvement, na.rm = TRUE),
    worsening_adj   = worsening - mean(worsening, na.rm = TRUE),
    side_effect_adj = side_effect - mean(side_effect, na.rm = TRUE),
    NRS_adj = NRS - mean(NRS, na.rm = TRUE),
    PDI_adj = PDI - mean(PDI, na.rm = TRUE)
  ) %>%
  ungroup()

# Person means based on time-point-centered values

personmeans_exp_adjusted <- data_exp_adjusted %>%
  group_by(VPN_Code_corr) %>%
  summarise(
    improvement   = mean(improvement_adj, na.rm = TRUE),
    worsening     = mean(worsening_adj, na.rm = TRUE),
    side_effect   = mean(side_effect_adj, na.rm = TRUE),
    NRS_exp_imp   = mean(NRS_exp_imp_adj, na.rm = TRUE),
    NRS_exp_worse = mean(NRS_exp_worse_adj, na.rm = TRUE),
    PDI_exp       = mean(PDI_exp_adj, na.rm = TRUE),
    .groups = "drop"
  )

personmeans_eff_adjusted <- data_eff_adjusted %>%
  group_by(VPN_Code_corr) %>%
  summarise(
    improvement = mean(improvement_adj, na.rm = TRUE),
    worsening   = mean(worsening_adj, na.rm = TRUE),
    side_effect = mean(side_effect_adj, na.rm = TRUE),
    NRS         = mean(NRS_adj, na.rm = TRUE),
    PDI         = mean(PDI_adj, na.rm = TRUE),
    .groups = "drop"
  ) 

# Correlations and bootstrap confidence intervals were calculated
# identically to the primary analyses, using person means based on
# time-point-centered values.

# 4. Sensitivity analyses for Hurdle Models

# For expectations, hurdle models were refit using measurement occasions T1-T6
# to examine whether effects were mainly driven by differences in expectations
# before versus after treatment.
# For better interpretation of the intercept, time was centered at T1

# data
exp_imp_df_sens <- exp_imp_df %>%
  filter(time >= 1) %>%
  mutate(
    time_post = time_months - 1.384615
  )

exp_worse_df_sens <- exp_worse_df %>%
  filter(time >= 1) %>%
  mutate(
    time_post = time_months - 1.384615
  )

exp_side_df_sens <- exp_side_df %>%
  filter(time >= 1) %>%
  mutate(
    time_post = time_months - 1.384615
  )

# Priors
pr_hurdle_exp_imp_sens <- c(
  # hurdle component
  set_prior("normal(0.0,0.75)", class = "b", dpar = "hu"),
  set_prior("normal(-2, 1.0)", class = "Intercept", dpar = "hu"),
  
  # lognormal component
  set_prior("normal(0.0,0.3)", class = "b", dpar = "mu"),
  set_prior("normal(1.8, 0.3)", class = "Intercept", dpar = "mu"),
  set_prior("normal(0, 0.5)", class = "sigma", lb = 0),
  
  # random effects
  set_prior("exponential(2)", class = "sd", dpar = "hu"),
  set_prior("exponential(2)", class = "sd", dpar = "mu")
)

pr_hurdle_exp_worse_sens <- c(
  # hurdle component
  set_prior("normal(0.0,0.75)", class = "b", dpar = "hu"),
  set_prior("normal(1, 1.0)", class = "Intercept", dpar = "hu"),
  
  # lognormal component
  set_prior("normal(0.0,0.3)", class = "b", dpar = "mu"),
  set_prior("normal(.9,.4)", class = "Intercept", dpar = "mu"),
  set_prior("normal(0, 0.5)", class = "sigma", lb = 0),
  
  # random effects
  set_prior("exponential(2)", class = "sd", dpar = "hu"),
  set_prior("exponential(2)", class = "sd", dpar = "mu")
)

pr_hurdle_exp_side_sens <- c(
  # hurdle component
  set_prior("normal(0.0,0.75)", class = "b", dpar = "hu"),
  set_prior("normal(0, 1.0)", class = "Intercept", dpar = "hu"),
  
  # lognormal component
  set_prior("normal(0.0,0.3)", class = "b", dpar = "mu"),
  set_prior("normal(0, 1)", class = "Intercept", dpar = "mu"),
  set_prior("normal(0, 0.5)", class = "sigma", lb = 0),
  
  # random effects
  set_prior("exponential(2)", class = "sd", dpar = "hu"),
  set_prior("exponential(2)", class = "sd", dpar = "mu")
)

# Fitted models were identical to those in the primary analyses, as were prior
# and posterior predictive checks.