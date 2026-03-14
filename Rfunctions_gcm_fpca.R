library(RcppHungarian)
library(pracma)
library(ggplot2)
library(dplyr)
library(tidyr)

bin_time_matrix <- function(Y, time_old, time_new, check = FALSE) {
  
  stopifnot(length(time_old) == ncol(Y))
  
  time_new <- sort(time_new)
  
  # define bin breaks
  breaks <- c(time_new, Inf)
  
  # assign each old time to a bin index
  bin_index <- findInterval(time_old, breaks, rightmost.closed = FALSE)
  
  # map to actual bin labels
  bin_label <- ifelse(bin_index == 0, NA, time_new[bin_index])
  
  # optional diagnostic
  if (check) {
    valid <- !is.na(bin_label)
    example_label <- unique(bin_label[valid])[1]
    
    cat("Example binning check:\n")
    cat("time_new bin =", example_label, "\n")
    cat("time_old values mapped to this bin:\n")
    print(time_old[bin_label == example_label])
    cat("\n")
  }
  
  # aggregate within bins
  Y_binned <- sapply(time_new, function(t) {
    cols <- which(bin_label == t)
    
    if (length(cols) == 0) {
      return(rep(NA_real_, nrow(Y)))
    }
    
    m <- rowMeans(Y[, cols, drop = FALSE], na.rm = TRUE)
    m[is.nan(m)] <- NA_real_
    m
  })
  
  colnames(Y_binned) <- time_new
  rownames(Y_binned) <- rownames(Y)
  
  return(Y_binned)
}



train_test_split <- function(Y, prop_train = 0.5) {
  is_vector <- is.null(dim(Y))
  
  if (is_vector) {
    n_subj <- length(Y)
    idx <- sample(seq_len(n_subj))
    
    n_train <- floor(prop_train * n_subj)
    train_idx <- idx[1:n_train]
    test_idx  <- idx[(n_train + 1):n_subj]
    
    train <- Y[train_idx]
    test  <- Y[test_idx]
    
  } else {
    n_subj <- nrow(Y)
    idx <- sample(seq_len(n_subj))
    
    n_train <- floor(prop_train * n_subj)
    train_idx <- idx[1:n_train]
    test_idx  <- idx[(n_train + 1):n_subj]
    
    train <- Y[train_idx, , drop = FALSE]
    test  <- Y[test_idx,  , drop = FALSE]
  }
  
  list(
    train = train,
    test  = test,
    train_idx = train_idx,
    test_idx  = test_idx
  )
}


introduce_missingness <- function(Y, missingness = 0.1, full = TRUE) {
  # Y: matrix of size (n_subjects × n_timepoints) OR vector
  # missingness: proportion of time points to be set missing
  # full: if TRUE → induce missingness across entire dataset
  #       if FALSE → induce missingness per subject
  
  Y_missing <- Y
  is_vector <- is.null(nrow(Y))
  
  if (is_vector) {
    # vector case
    n_subj <- 1
    n_time <- length(Y)
  } else {
    n_subj <- nrow(Y)
    n_time <- ncol(Y)
  }
  
  if (full) {
    # ---- Full dataset missingness ----
    total_points <- n_subj * n_time
    n_missing <- ceiling(missingness * total_points)
    
    # Randomly select global positions
    miss_idx <- sample(total_points, n_missing, replace = FALSE)
    
    # Convert to matrix indices
    if (is_vector) {
      Y_missing[miss_idx] <- NA
    } else {
      row_idx <- ((miss_idx - 1) %% n_subj) + 1
      col_idx <- ((miss_idx - 1) %/% n_subj) + 1
      Y_missing[cbind(row_idx, col_idx)] <- NA
    }
    
  } else {
    # ---- Subject-by-subject missingness ----
    for (i in 1:n_subj) {
      n_missing <- ceiling(missingness * n_time)
      miss_idx <- sample(1:n_time, n_missing, replace = FALSE)
      
      if (is_vector) {
        Y_missing[miss_idx] <- NA
      } else {
        Y_missing[i, miss_idx] <- NA
      }
    }
  }
  
  return(Y_missing)
}

LyLt_to_Y <- function(Ly, Lt, time) {
  
  n_subj <- length(Ly)
  n_time <- length(time)
  
  Y_quad <- matrix(NA_real_, ncol = n_time, nrow = n_subj)
  
  for (i in seq_len(n_subj)) {
    idx <- match(Lt[[i]], time)
    Y_quad[i, idx] <- Ly[[i]]
  }
  
  Y_quad
}

run_mc <- function(iter, n, time_points, cor_ss, sigma_eps, f1_t, f2_t){

  set.seed(1000 + iter)

  # ----- Parameters -----
  n_subj <- n
  time <- seq(0, 1, length.out = time_points)
  n_time <- length(time)
  missingness <- .3

  # ----- Fixed effects -----
  beta_s <- 1
  beta_s2 <- 1

  # ----- Random effects -----
  sigma_s <- 1
  sigma_s2 <- sqrt(2)
  rho_ss2 <- cor_ss

  Sigma <- matrix(c(sigma_s^2, rho_ss2 * sigma_s * sigma_s2,
                    rho_ss2 * sigma_s * sigma_s2, sigma_s2^2), 2, 2)

  b_s <- MASS::mvrnorm(n_subj, mu = c(0, 0), Sigma = Sigma)
  s_i  <- b_s[,1]
  s2_i <- b_s[,2]

  # ----- Generate Data -----
  Y <- matrix(NA, nrow = n_subj, ncol = n_time)
  for(i in 1:n_subj){
    Y[i,] <- (beta_s + s_i[i]) * f1_t +
             (beta_s2 + s2_i[i]) * f2_t +
             rnorm(n_time, 0, sigma_eps)
  }

  Y_miss <- introduce_missingness(Y, missingness)
  Y_quad_split <- train_test_split(Y_miss)
  Y_quad <- Y_quad_split$train
  Y_test <- Y_quad_split$test

  # ----- FPCA -----
  Ly_quad <- lapply(seq_len(nrow(Y_quad)), function(i) Y_quad[i,])
  Lt_quad <- lapply(seq_len(nrow(Y_quad)), function(i) time)

  fpca_quad <- FPCA(Ly_quad, Lt_quad)

  fit_ng <- ConvertSupport(
    fromGrid = fpca_quad$workGrid,
    toGrid   = time,
    phi      = fpca_quad$phi
  )

  loading1 <- fit_ng[,1]
  loading2 <- fit_ng[,2]
  
  ise <- compare_efunctions_ise(true_phi = cbind(f1_t, f2_t), est_phi = fit_ng, T = time)$ise
  names(ise) <- c("ise_f1", "ise_f2")

  # ----- Lavaan Setup -----
  dat_lav <- as.data.frame(Y_test)
  colnames(dat_lav) <- paste0("y", 1:n_time)

  l1 <- paste0(loading1, "*y", 1:n_time, collapse = " + ")
  l2 <- paste0(loading2, "*y", 1:n_time, collapse = " + ")

  growth_mod <- paste0("
    eta1 =~ ", l1, "
    eta2 =~ ", l2, "

    eta1 ~~ var_eta1*eta1
    eta2 ~~ var_eta2*eta2
    eta1 ~~ cov_sl*eta2

    eta1 ~ mu1*1
    eta2 ~ mu2*1

    ", paste0(colnames(dat_lav), " ~~ vare*", colnames(dat_lav), collapse = "\n"), "
  ")

  fit_gcm <- try(
    growth(
      model   = growth_mod,
      data    = dat_lav,
      missing = "fiml",
      em.h1.iter.max = 10000
    ),
    silent = TRUE
  )

  if(inherits(fit_gcm, "try-error")){
    return(list(fit_est = rep(NA, 12), phi = c(loading1, loading2)))
  }

  f <- summary(fit_gcm, standardized = TRUE, fit.measures = TRUE)

  pe_vals <- f$pe %>%
    filter(label %in% c("var_eta1","var_eta2","cov_sl","mu1","mu2")) %>%
    pull(est) %>%
    abs()
  
  names(pe_vals) <- c("var_eta1","var_eta2","abs_cov_sl","mu1","mu2")

  vare_mean <- mean(f$pe[f$pe$label == "vare","est"])
  
  names(vare_mean) <- "vare"

  fit_vals <- f$fit[c("cfi","bic","rmsea","logl","npar","df")]

  list(fit_est = c(pe_vals, vare_mean, fit_vals, ise), phi = cbind(loading1, loading2))
}


run_mc_lb <- function(iter, n, time_points, cor_ss, sigma_eps, f1_t, f2_t){
  
  set.seed(1000 + iter)
  
  # ----- Parameters -----
  n_subj <- n
  time <- seq(0, 1, length.out = time_points)
  n_time <- length(time)
  missingness <- .3
  
  # ----- Fixed effects -----
  beta_s <- 1
  beta_s2 <- 1
  
  # ----- Random effects -----
  sigma_s <- 1
  sigma_s2 <- sqrt(2)
  rho_ss2 <- cor_ss
  
  Sigma <- matrix(c(sigma_s^2, rho_ss2 * sigma_s * sigma_s2,
                    rho_ss2 * sigma_s * sigma_s2, sigma_s2^2), 2, 2)
  
  b_s <- MASS::mvrnorm(n_subj, mu = c(0, 0), Sigma = Sigma)
  s_i  <- b_s[,1]
  s2_i <- b_s[,2]
  
  # ----- Generate Data -----
  Y <- matrix(NA, nrow = n_subj, ncol = n_time)
  for(i in 1:n_subj){
    Y[i,] <- (beta_s + s_i[i]) * f1_t +
      (beta_s2 + s2_i[i]) * f2_t +
      rnorm(n_time, 0, sigma_eps)
  }
  
  Y_miss <- introduce_missingness(Y, missingness)
  Y_quad_split <- train_test_split(Y_miss)
  Y_train <- Y_quad_split$train
  Y_test  <- Y_quad_split$test
  
  # LATENT BASIS TRAIN STEP 
  
  betas <- paste0("beta", 2:(n_time - 1))
  loading1 <- rep(1, n_time)
  loading2 <- c(0, betas, 1)
  
  dat_lav <- as.data.frame(Y_train)
  colnames(dat_lav) <- paste0("y", 1:n_time)
  
  l1 <- paste0(loading1, "*y", 1:n_time, collapse = " + ")
  l2 <- paste0(loading2, "*y", 1:n_time, collapse = " + ")
  
  growth_mod_train <- paste0("
    eta1 =~ ", l1, "
    eta2 =~ ", l2, "

    eta1 ~~ eta1
    eta2 ~~ eta2
    eta1 ~~ eta2

    ", paste0(colnames(dat_lav), " ~~ vare*", colnames(dat_lav), collapse = "\n"), "
  ")
  
  fit_train <- try(
    growth(model = growth_mod_train,
           data  = dat_lav,
           missing = "fiml",
           em.h1.iter.max = 10000),
    silent = TRUE
  )
  
  if(inherits(fit_train, "try-error")){
    return(list(fit_est = rep(NA, 12),
                phi = cbind(rep(NA,n_time), rep(NA,n_time))))
  }
  
  fit_train_sum <- summary(fit_train, standardized = TRUE)
  
  loading2_est <- c(
    0,
    fit_train_sum$pe[fit_train_sum$pe$label %in% betas, "est"],
    1
  )
  
  # TEST STEP (FIX LOADINGS) 
  
  dat_lav <- as.data.frame(Y_test)
  colnames(dat_lav) <- paste0("y", 1:n_time)
  
  l1 <- paste0(loading1, "*y", 1:n_time, collapse = " + ")
  l2 <- paste0(loading2_est, "*y", 1:n_time, collapse = " + ")
  
  growth_mod_test <- paste0("
    eta1 =~ ", l1, "
    eta2 =~ ", l2, "

    eta1 ~~ var_eta1*eta1
    eta2 ~~ var_eta2*eta2
    eta1 ~~ cov_sl*eta2

    eta1 ~ mu1*1
    eta2 ~ mu2*1

    ", paste0(colnames(dat_lav), " ~~ vare*", colnames(dat_lav), collapse = "\n"), "
  ")
  
  fit_test <- try(
    growth(model = growth_mod_test,
           data  = dat_lav,
           missing = "fiml",
           em.h1.iter.max = 10000),
    silent = TRUE
  )
  
  if(inherits(fit_test, "try-error")){
    return(list(fit_est = rep(NA, 12),
                phi = cbind(loading1, loading2_est)))
  }
  
  f <- summary(fit_test, standardized = TRUE, fit.measures = TRUE)
  
  # Extract parameters
  pe_vals <- f$pe %>%
    dplyr::filter(label %in% c("var_eta1","var_eta2","cov_sl","mu1","mu2")) %>%
    dplyr::pull(est) %>%
    abs()
  
  names(pe_vals) <- c("var_eta1","var_eta2","abs_cov_sl","mu1","mu2")
  
  vare_mean <- mean(f$pe[f$pe$label == "vare","est"])
  names(vare_mean) <- "vare"
  
  fit_vals <- f$fit[c("cfi","bic","rmsea","logl","npar","df")]
  
  # ISE CALCULATION

  ise <- compare_efunctions_ise(
    true_phi = cbind(f1_t, f2_t),
    est_phi  = cbind(loading1, loading2_est),
    T = time
  )$ise
  
  names(ise) <- c("ise_f1","ise_f2")
  
  list(
    fit_est = c(pe_vals, vare_mean, fit_vals, ise),
    phi     = cbind(loading1, loading2_est)
  )
}

summarize_condition <- function(cond_name, mc_results, metric_names) {

  # Extract all iterations for this condition
  vals <- lapply(mc_results[[cond_name]], function(x) {
    if (is.null(x) || all(is.na(x))) {
      rep(NA, length(metric_names))
    } else {
      x[metric_names]
    }
  })

  mat <- do.call(rbind, vals)

  means <- colMeans(mat, na.rm = TRUE)
  sds   <- apply(mat, 2, sd, na.rm = TRUE)

  out <- c(
    setNames(means, paste0(metric_names, "_mean")),
    setNames(sds,   paste0(metric_names, "_sd"))
  )

  return(out)
}

compare_efunctions_ise <- function(true_phi, est_phi, T = 1:125) {

  if (is.vector(true_phi)) true_phi <- matrix(true_phi, ncol = 1)
  if (is.vector(est_phi))  est_phi  <- matrix(est_phi,  ncol = 1)
  
  n_true <- ncol(true_phi)
  n_est  <- ncol(est_phi)
  
  # Step 1: Compute ISE matrix

  ise_matrix <- matrix(NA, nrow = n_true, ncol = n_est)
  sign_matrix <- matrix(1, nrow = n_true, ncol = n_est)
  
  for (i in 1:n_true) {
    for (j in 1:n_est) {
      
      diff_pos <- (est_phi[, j] - true_phi[, i])^2
      diff_neg <- (-est_phi[, j] - true_phi[, i])^2
      
      ise_pos <- trapz(T, diff_pos)
      ise_neg <- trapz(T, diff_neg)
      
      if (ise_neg < ise_pos) {
        ise_matrix[i, j]  <- ise_neg
        sign_matrix[i, j] <- -1
      } else {
        ise_matrix[i, j]  <- ise_pos
        sign_matrix[i, j] <- 1
      }
    }
  }
  
  # Step 2: Optimal assignment

  assignment <- HungarianSolver(ise_matrix)
  
  results <- lapply(1:n_true, function(i) {
    
    j <- assignment$pairs[i, 2]
    sgn <- sign_matrix[i, j]
    
    est_aligned <- sgn * est_phi[, j]
    true_vec    <- true_phi[, i]
    
    r2 <- summary(lm(est_aligned ~ true_vec))$r.squared
    
    data.frame(
      true_phi_column = i,
      best_match_est_phi_column = j,
      sign = sgn,
      ise = ise_matrix[i, j],
      r_squared = r2
    )
  })
  
  match_results <- do.call(rbind, results)
  
  return(match_results)
}

get_metric <- function(cond_name, summary_df, n_mc, metric) {
  
  mean_val <- summary_df[
    summary_df$metric == paste0(metric,"_mean"),
    cond_name
  ]
  
  sd_val <- summary_df[
    summary_df$metric == paste0(metric,"_sd"),
    cond_name
  ]
  
  ci <- 1.96 * as.numeric(sd_val) / sqrt(n_mc)
  
  list(
    mean = as.numeric(mean_val),
    lower = as.numeric(mean_val) - ci,
    upper = as.numeric(mean_val) + ci
  )
}

plot_simulation_metric <- function(summary_df,
                                   n_mc,
                                   metric,
                                   abline_value = NULL,
                                   main_title = NULL,
                                   fix_y = FALSE, 
                                   pad = .01) {
  
  get_metric_vals <- function(cond_name) {
    
    mean_val <- summary_df[
      summary_df$metric == paste0(metric, "_mean"),
      cond_name
    ]
    
    sd_val <- summary_df[
      summary_df$metric == paste0(metric, "_sd"),
      cond_name
    ]
    
    ci <- 1.96 * as.numeric(sd_val) / sqrt(n_mc)
    
    data.frame(
      condition = cond_name,
      mean  = as.numeric(mean_val),
      lower = as.numeric(mean_val) - ci,
      upper = as.numeric(mean_val) + ci
    )
  }
  
  panels <- list(
    n = list(
      conds = c(
        "n100_tp30_cor0_eps0.25",
        "n200_tp30_cor0_eps0.25",
        "n500_tp30_cor0_eps0.25"
      ),
      x_vals = c(100,200,500),
      xlab = "Sample Size (n)"
    ),
    
    tp = list(
      conds = c(
        "n500_tp12_cor0_eps0.25",
        "n500_tp30_cor0_eps0.25"
      ),
      x_vals = c(12,30),
      xlab = "Time Points"
    ),
    
    cor = list(
      conds = c(
        "n500_tp30_cor0_eps0.25",
        "n500_tp30_cor0.3_eps0.25",
        "n500_tp30_cor0.8_eps0.25"
      ),
      x_vals = c(0,0.3,0.8),
      xlab = "RI–S Covariance"
    ),
    
    eps = list(
      conds = c(
        "n500_tp30_cor0_eps0.25",
        "n500_tp30_cor0_eps0.5",
        "n500_tp30_cor0_eps1"
      ),
      x_vals = c(0.25,0.5,1),
      xlab = "Measurement Error"
    )
  )
  
  plot_df <- bind_rows(lapply(names(panels), function(p) {
    
    panel_info <- panels[[p]]
    
    df <- bind_rows(lapply(panel_info$conds, get_metric_vals))
    df$x <- panel_info$x_vals
    df$panel <- panel_info$xlab
    
    df
  }))
  
  palette_cols <- c(
    "#0072B2",
    "#D55E00",
    "#009E73",
    "#CC79A7"
  )
  
  p <- ggplot(plot_df,
              aes(x = x,
                  y = mean,
                  group = panel,
                  color = panel)) +
    
    geom_line(size = 1.1) +
    geom_point(size = 3) +
    
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  width = 0,
                  size = .8) +
    
    facet_wrap(~panel, scales = "free_x", ncol = 2) +
    
    scale_color_manual(values = palette_cols) +
    
    labs(
      x = NULL,
      y = toupper(metric),
      title = main_title
    ) +
    
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 12),
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  
  if (!is.null(abline_value)) {
    p <- p + geom_hline(
      yintercept = abline_value,
      linetype = "dashed",
      size = 1
    )
  }
  
  if (fix_y) {
    
    ymin <- min(plot_df$lower, na.rm = TRUE)
    ymax <- max(plot_df$upper, na.rm = TRUE)
    
    if (!is.null(abline_value)) {
      ymin <- min(ymin, abline_value)
      ymax <- max(ymax, abline_value)
    }
    
    padding <- pad*(ymax - ymin)
    
    p <- p + coord_cartesian(
      ylim = c(ymin - padding, ymax + padding)
    )
  }
  
  return(p)
}

plot_aligned_efunctions <- function(list_best_cond,
                                    f1_t,
                                    f2_t,
                                    n_mc,
                                    T_grid = 1:length(f1_t)) {
  
  true_phi <- cbind(f1_t, f2_t)
  
  all_estimates <- list()
  counter <- 1
  
  for (i in seq(2, 2*n_mc, by = 2)) {
    
    est_phi <- list_best_cond[i]$phi
    
    # check if NULL
    if (is.null(est_phi)) {
      counter <- counter + 1
      next
    }
    
    match_results <- compare_efunctions_ise(
      true_phi = true_phi,
      est_phi  = est_phi,
      T        = T_grid
    )
    
    aligned_est <- matrix(NA,
                          nrow = length(T_grid),
                          ncol = ncol(true_phi))
    
    for (k in 1:nrow(match_results)) {
      j   <- match_results$best_match_est_phi_column[k]
      sgn <- match_results$sign[k]
      aligned_est[, k] <- sgn * est_phi[, j]
    }
    
    all_estimates[[counter]] <- data.frame(
      T = rep(T_grid, 2),
      value = c(aligned_est[,1], aligned_est[,2]),
      efunc = rep(c("f1","f2"), each = length(T_grid)),
      replication = counter
    )
    
    counter <- counter + 1
  }
  
  df_est <- dplyr::bind_rows(all_estimates)
  
  df_true <- data.frame(
    T = rep(T_grid, 2),
    value = c(f1_t, f2_t),
    efunc = rep(c("f1","f2"), each = length(T_grid))
  )
  
  # ---- Plot f1 ----
  p_f1 <- ggplot() +
    geom_line(data = subset(df_est, efunc == "f1"),
              aes(x = T,
                  y = value,
                  group = replication),
              alpha = 0.2,
              linewidth = 0.6,
              color = "#0072B2") +
    
    geom_line(data = subset(df_true, efunc == "f1"),
              aes(x = T, y = value),
              linewidth = 1.6,
              color = "#0072B2") +
    
    theme_minimal(base_size = 14) +
    labs(
      x = "t",
      y = expression(phi[1](t)),
      title = "Eigenfunction 1 Across Replications"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  
  # ---- Plot f2 ----
  p_f2 <- ggplot() +
    geom_line(data = subset(df_est, efunc == "f2"),
              aes(x = T,
                  y = value,
                  group = replication),
              alpha = 0.2,
              linewidth = 0.6,
              color = "#D55E00") +
    
    geom_line(data = subset(df_true, efunc == "f2"),
              aes(x = T, y = value),
              linewidth = 1.6,
              color = "#D55E00") +
    
    theme_minimal(base_size = 14) +
    labs(
      x = "t",
      y = expression(phi[2](t)),
      title = "Eigenfunction 2 Across Replications"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  
  return(list(p_f1 = p_f1, p_f2 = p_f2))
}

plot_simulation_metric2 <- function(summary_df1,
                                   summary_df2,
                                   n_mc,
                                   metric,
                                   labels = c("Method 1", "Method 2"),
                                   abline_value = NULL,
                                   main_title = NULL) {
  
  get_metric_vals <- function(summary_df, cond_name, method_label) {
    
    mean_val <- summary_df[
      summary_df$metric == paste0(metric, "_mean"),
      cond_name
    ]
    
    sd_val <- summary_df[
      summary_df$metric == paste0(metric, "_sd"),
      cond_name
    ]
    
    ci <- 1.96 * as.numeric(sd_val) / sqrt(n_mc)
    
    data.frame(
      condition = cond_name,
      mean  = as.numeric(mean_val),
      lower = as.numeric(mean_val) - ci,
      upper = as.numeric(mean_val) + ci,
      method = method_label
    )
  }
  
  ### Define conditions
  panels <- list(
    n = list(
      conds = c(
        "n100_tp30_cor0_eps0.25",
        "n200_tp30_cor0_eps0.25",
        "n500_tp30_cor0_eps0.25"
      ),
      x_vals = c(100,200,500),
      xlab = "Sample Size (n)"
    ),
    
    tp = list(
      conds = c(
        "n500_tp12_cor0_eps0.25",
        "n500_tp30_cor0_eps0.25"
      ),
      x_vals = c(12,30),
      xlab = "Time Points"
    ),
    
    cor = list(
      conds = c(
        "n500_tp30_cor0_eps0.25",
        "n500_tp30_cor0.3_eps0.25",
        "n500_tp30_cor0.8_eps0.25"
      ),
      x_vals = c(0,0.3,0.8),
      xlab = "RI–S Covariance"
    ),
    
    eps = list(
      conds = c(
        "n500_tp30_cor0_eps0.25",
        "n500_tp30_cor0_eps0.5",
        "n500_tp30_cor0_eps1"
      ),
      x_vals = c(0.25,0.5,1),
      xlab = "Measurement Error"
    )
  )
  
  ### Build long dataframe for BOTH summaries
  plot_df <- bind_rows(lapply(names(panels), function(p) {
    
    panel_info <- panels[[p]]
    
    df1 <- bind_rows(lapply(panel_info$conds,
                            function(cond)
                              get_metric_vals(summary_df1, cond, labels[1])))
    
    df2 <- bind_rows(lapply(panel_info$conds,
                            function(cond)
                              get_metric_vals(summary_df2, cond, labels[2])))
    
    df <- bind_rows(df1, df2)
    
    df$x <- rep(panel_info$x_vals, times = 2)
    df$panel <- panel_info$xlab
    
    df
  }))
  
  ### Nice contrasting colors
  palette_cols <- c("#009E73", "#CC79A7")
  
  p <- ggplot(plot_df,
              aes(x = x,
                  y = mean,
                  group = method,
                  color = method)) +
    
    geom_line(size = 1.1) +
    geom_point(size = 3) +
    
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  width = 0,
                  size = .8) +
    
    facet_wrap(~panel, scales = "free_x", ncol = 2) +
    
    scale_color_manual(values = palette_cols) +
    
    labs(
      x = NULL,
      y = toupper(metric),
      color = NULL,
      title = main_title
    ) +
    
    theme_minimal(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold", size = 12),
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  
  if (!is.null(abline_value)) {
    p <- p + geom_hline(yintercept = abline_value,
                        linetype = "dashed",
                        size = 1)
  }
  
  return(p)
}
