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

run_mc <- function(iter, n, time_points, cor_ss, f1_t, f2_t){

  set.seed(1000 + iter)

  # ----- Parameters -----
  n_subj <- n
  time <- seq(0, 1, length.out = time_points)
  n_time <- length(time)
  missingness <- .3

  # ----- Fixed effects -----
  beta_s <- 0
  beta_s2 <- 0

  # ----- Random effects -----
  sigma_s <- 1
  sigma_s2 <- 1
  rho_ss2 <- cor_ss

  Sigma <- matrix(c(sigma_s^2, rho_ss2 * sigma_s * sigma_s2,
                    rho_ss2 * sigma_s * sigma_s2, sigma_s2^2), 2, 2)

  b_s <- MASS::mvrnorm(n_subj, mu = c(0, 0), Sigma = Sigma)
  s_i  <- b_s[,1]
  s2_i <- b_s[,2]

  sigma_eps <- 0.5

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
    pull(est)
  
  names(pe_vals) <- c("var_eta1","var_eta2","cov_sl","mu1","mu2")

  vare_mean <- mean(f$pe[f$pe$label == "vare","est"])
  
  names(vare_mean) <- "vare"

  fit_vals <- f$fit[c("cfi","bic","rmsea","logl","npar","df")]

  list(fit_est = c(pe_vals, vare_mean, fit_vals), phi = cbind(loading1, loading2))
}


run_mcmc_np <- function(iter, n, time_points, cor_ss, f1_t, f2_t){

  set.seed(1000 + iter)

  # ----- Setup -----
  n_subj <- n
  time <- seq(0, 1, length.out = time_points)
  n_time <- length(time)
  missingness <- .3

  beta_s  <- 0
  beta_s2 <- 0

  sigma_s  <- 1
  sigma_s2 <- 1
  rho_ss2  <- cor_ss

  Sigma <- matrix(
    c(sigma_s^2,
      rho_ss2 * sigma_s * sigma_s2,
      rho_ss2 * sigma_s * sigma_s2,
      sigma_s2^2),
    2, 2
  )

  b_s  <- MASS::mvrnorm(n_subj, mu = c(0, 0), Sigma = Sigma)
  s_i  <- b_s[,1]
  s2_i <- b_s[,2]

  sigma_eps <- 0.5

  # ----- Generate Data -----
  Y <- matrix(NA, nrow = n_subj, ncol = n_time)
  for(i in 1:n_subj){
    Y[i,] <- (beta_s + s_i[i])  * f1_t +
             (beta_s2 + s2_i[i]) * f2_t +
             rnorm(n_time, 0, sigma_eps)
  }

  Y_miss <- introduce_missingness(Y, missingness)
  Y_split <- train_test_split(Y_miss)

  Y_train <- Y_split$train
  Y_test  <- Y_split$test

  # =========================
  # ----- FPCA GCM ---------
  # =========================

  Ly <- lapply(seq_len(nrow(Y_train)), function(i) Y_train[i,])
  Lt <- lapply(seq_len(nrow(Y_train)), function(i) time)

  fpca_fit <- FPCA(Ly, Lt)

  fit_ng <- ConvertSupport(
    fromGrid = fpca_fit$workGrid,
    toGrid   = time,
    phi      = fpca_fit$phi
  )

  loading1 <- fit_ng[,1]
  loading2 <- fit_ng[,2]

  dat_lav <- as.data.frame(Y_test)
  colnames(dat_lav) <- paste0("y", 1:n_time)

  l1 <- paste0(loading1, "*y", 1:n_time, collapse = " + ")
  l2 <- paste0(loading2, "*y", 1:n_time, collapse = " + ")

  growth_mod_np <- paste0("
    eta1 =~ ", l1, "
    eta2 =~ ", l2, "

    eta1 ~~ eta1
    eta2 ~~ eta2
    eta1 ~~ eta2
    
    eta1 ~ mu1*1
    eta2 ~ mu2*1

    ", paste0(colnames(dat_lav), " ~~ vare*", colnames(dat_lav), collapse = "\n"), "
  ")

  fit_np <- try(
    growth(
      model = growth_mod_np,
      data  = dat_lav,
      missing = "fiml",
      em.h1.iter.max = 10000
    ),
    silent = TRUE
  )

  if(inherits(fit_np, "try-error")){
    return(list(
      fit_est       = rep(NA, 12),
      fit_est_param = rep(NA, 12),
      phi           = cbind(loading1, loading2)
    ))
  }

  sum_np <- summary(fit_np, standardized = TRUE, fit.measures = TRUE)

  fit_est <- c(
    var_eta1 = sum_np$pe$est[sum_np$pe$lhs=="eta1" & sum_np$pe$rhs=="eta1"],
    var_eta2 = sum_np$pe$est[sum_np$pe$lhs=="eta2" & sum_np$pe$rhs=="eta2"],
    cov_sl   = sum_np$pe$est[sum_np$pe$lhs=="eta1" & sum_np$pe$rhs=="eta2"],
    mu1      = sum_np$pe$est[sum_np$pe$lhs=="eta1" & sum_np$pe$op=="~" & sum_np$pe$rhs=="1"],
    mu2      = sum_np$pe$est[sum_np$pe$lhs=="eta2" & sum_np$pe$op=="~" & sum_np$pe$rhs=="1"],
    vare     = mean(sum_np$pe$est[sum_np$pe$label=="vare"]),
    sum_np$fit[c("cfi","bic","rmsea","logl","npar","df")]
  )

  # =========================
  # ----- Parametric GCM ---
  # =========================

  loading1_p <- rep(1, n_time)
  loading2_p <- seq(0, 1, length.out = n_time)

  l1p <- paste0(loading1_p, "*y", 1:n_time, collapse = " + ")
  l2p <- paste0(loading2_p, "*y", 1:n_time, collapse = " + ")

  growth_mod_param <- paste0("
    eta1 =~ ", l1p, "
    eta2 =~ ", l2p, "

    eta1 ~~ eta1
    eta2 ~~ eta2
    eta1 ~~ eta2
    
    eta1 ~ mu1*1
    eta2 ~ mu2*1

    ", paste0(colnames(dat_lav), " ~~ vare*", colnames(dat_lav), collapse = "\n"), "
  ")

  fit_param <- try(
    growth(
      model = growth_mod_param,
      data  = dat_lav,
      missing = "fiml",
      em.h1.iter.max = 10000
    ),
    silent = TRUE
  )

  if(inherits(fit_param, "try-error")){
    fit_est_param <- rep(NA, length(fit_est))
  } else {

    sum_param <- summary(fit_param, standardized = TRUE, fit.measures = TRUE)

    fit_est_param <- c(
      var_eta1 = sum_param$pe$est[sum_param$pe$lhs=="eta1" & sum_param$pe$rhs=="eta1"],
      var_eta2 = sum_param$pe$est[sum_param$pe$lhs=="eta2" & sum_param$pe$rhs=="eta2"],
      cov_sl   = sum_param$pe$est[sum_param$pe$lhs=="eta1" & sum_param$pe$rhs=="eta2"],
      mu1      = sum_param$pe$est[sum_param$pe$lhs=="eta1" & sum_param$pe$op=="~" & sum_param$pe$rhs=="1"],
      mu2      = sum_param$pe$est[sum_param$pe$lhs=="eta2" & sum_param$pe$op=="~" & sum_param$pe$rhs=="1"],
      vare     = mean(sum_param$pe$est[sum_param$pe$label=="vare"]),
      sum_param$fit[c("cfi","bic","rmsea","logl","npar","df")]
    )
  }

  # ----- Return -----
  list(
    fit_est       = fit_est,
    fit_est_param = fit_est_param,
    phi           = cbind(loading1, loading2)
  )
}


