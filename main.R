rm(list = ls())  
library(dplyr) 
library(ggplot2)
library(tidyverse)
library(ggplot2)
library(mvtnorm)
library(MASS) 
# Hierachical model
simu_hier = function(alpha1, alpha2, tau1, tau2, sigma1, sigma2, n1, n2, ncells){
  #simulates hierarchical data with y_ij ~ N (alpha_i, sigma^2) where alpha_i ~ N(alpha, tau^2) where y_ij ~ f_i and E[f_i] ~ N(alpha, sigma^2 + tau^2)
  #sigma, tau will represent variance, NOT the sd
  #n1, n2 = # of individuals in each group respectively
  #ncells = # of cells 
  Y1 = matrix(data = NA, nrow = n1, ncol = ncells)
  Y2 = matrix(data = NA, nrow = n2, ncol = ncells)
  for (i in 1:n1) { #GROUP 1
    alpha_1i = rnorm(1, alpha1, sqrt(tau1))
    Y1[i,] = rnorm(ncells, mean = alpha_1i, sd = sqrt(sigma1))
  }
  for (i in 1:n2) { #GROUP 2 
    alpha_2i = rnorm(1, alpha2, sqrt(tau2))
    Y2[i,] = rnorm(ncells, mean = alpha_2i, sd = sqrt(sigma2))
  }
  return(list(Y1 = Y1, Y2 = Y2))
} 
# Compute the moments 
compute_moments = function(y, p) {
  return(sapply(1:p, function(k) y^k))
} 
# Compute the Covariance matrix for Moments 
mom_cov = function(Y, p) {
  n = nrow(Y)
  total_samples <- 0
  overall_sum <- rep(0, p)  # Sum for overall mean
  within_var_sum <- matrix(0, nrow = p, ncol = p)  # Within-individual variance
  between_var_sum <- matrix(0, nrow = p, ncol = p)  # Between-individual variance 
  
  mean_list = matrix(0, nrow = nrow(Y), ncol = p)
  m_vector <- numeric(n)  # Store the number of samples for each individual
  for (i in 1:nrow(Y)){
    samples_i = Y[i, ]
    m_i = length(samples_i) 
    m_vector[i] = m_i
    moments_i = t(sapply(samples_i, compute_moments, p = p)) 
    # Compute the mean of the moments for this individual
    mean_t_i = colMeans(moments_i) 
    mean_list[i, ] = mean_t_i
    # Update overall sum for computing the overall mean
    overall_sum = overall_sum + m_i * mean_t_i
    total_samples = total_samples + m_i 
    # Compute within-individual variance
    diffs = t(moments_i) - mean_t_i
    within_var_sum = within_var_sum + diffs%*%t(diffs) 
    # for (j in 1:m_i) {
    #   diff <- moments_i[j, ] - mean_t_i
    #   within_var_sum <- within_var_sum + (diff %*% t(diff))
    # } 
  }
  # Compute the overall mean of the moments
  overall_mean <- overall_sum / total_samples
  for (i in 1:nrow(Y)){
    m_i = m_vector[i]
    mean_t_i = mean_list[i, ]
    diff <- mean_t_i - overall_mean
    between_var_sum <- between_var_sum + m_i^2 * (diff %*% t(diff)) # Should be m_i^2 not m_i
  }
  Cov = within_var_sum/total_samples^2 + between_var_sum/total_samples^2
  return(Cov)
}

# Function to get bin counts for a whole matrix 
get_bin_counts = function(Y, bin_edges, K) {
  # Flatten the matrix to use cut and table efficiently
  Y_flat = as.vector(t(Y))  # t(Y) or Y? 
  # Use cut() to classify each element of Y into bins defined by bin_edges
  bin_indices <- cut(Y_flat, breaks = bin_edges, include.lowest = TRUE, labels = FALSE)
  # Count the occurrences in each bin per row
  counts_matrix <- matrix(0, nrow = nrow(Y), ncol = K) 
  for (i in 1:nrow(Y)) {
    row_indices = bin_indices[ ((i - 1) * ncol(Y) + 1):(i * ncol(Y)) ]
    counts_matrix[i, ] = tabulate(row_indices, nbins = K) 
  }
  return(counts_matrix)
}
## Jackknife estimator 
jackknife_se = function(Y1, Y2, est_grid, K, p, ncells) {
  n1 <- nrow(Y1)
  n2 <- nrow(Y2)
  n = n1 + n2 
  # Store all Jackknife estimates
  jackknife_estimates <- matrix(0, nrow = n, ncol = p + 1) 
  # p + 1 for Intercept and coefficients 
  for (i in 1:n) {
    if (i <= n1){
      # Remove the ith individual from group 1
      Y1_jack <- Y1[-i, ]
      Y2_jack <- Y2
    } else {
      # Remove the (i - n1)th individual from group 2
      Y1_jack <- Y1
      Y2_jack <- Y2[-(i - n1), ]
    }
    # Calculate Smat1 and Smat2 using the helper function `get_bin_counts`
    Smat1_jack <- get_bin_counts(Y1_jack, bin_edges = est_grid, K)
    Smat2_jack <- get_bin_counts(Y2_jack, bin_edges = est_grid, K)
    Smat_jack <- rbind(Smat1_jack, Smat2_jack) 
    tk <- (est_grid[-1] + est_grid[-length(est_grid)]) / 2 
    X <- rep(1, length(tk))
    for (dim in 1:p) {
      X <- cbind(X, tk^dim)
    }
    varb <- paste0("X_", 1:p)
    colnames(X) <- c("Intercept", varb)
    
    
    S1sum_jack <- colSums(Smat1_jack)
    S2sum_jack <- colSums(Smat2_jack)
    carrier_jack <- density(x = c(as.vector(Y1_jack), as.vector(Y2_jack)), kernel = "gaussian", n = K, from = min(c(Y1_jack, Y2_jack)), to = max(c(Y1_jack, Y2_jack))) 
    scale_factor <- (ncells * nrow(Smat_jack)) / sum(carrier_jack$y) 
    carrier_est_jack <- carrier_jack$y * scale_factor 
    
    df1_jack <- as.data.frame(cbind(S1sum_jack, carrier_est_jack, X))
    df2_jack <- as.data.frame(cbind(S2sum_jack, carrier_est_jack, X))
    colnames(df1_jack)[1] <- "sum_cts"
    colnames(df2_jack)[1] <- "sum_cts"
    formula <- as.formula(paste0('sum_cts ~ offset(log(carrier_est_jack)) + ', paste(varb, collapse = '+')))
    # Fit GLM model for the Jackknife sample
    fit_gp1_jack <- tryCatch(glm(formula, family = poisson(link = "log"), data = df1_jack), error = function(e) NULL)
    fit_gp2_jack <- tryCatch(glm(formula, family = poisson(link = "log"), data = df2_jack), error = function(e) NULL)
    if (!is.null(fit_gp1_jack) && !is.null(fit_gp2_jack)) {
      beta1_jack <- as.vector(fit_gp1_jack$coefficients)
      beta2_jack <- as.vector(fit_gp2_jack$coefficients)
      # Store the difference in estimates
      jackknife_estimates[i, ] <- beta1_jack - beta2_jack
    } else {
      jackknife_estimates[i, ] <- NA
    }
  }
  mean_estimates = colMeans(jackknife_estimates, na.rm = TRUE) 
  jackknife_var = (n - 1) / n * rowSums((t(jackknife_estimates) - mean_estimates)^2, na.rm = TRUE)
  jackknife_se = sqrt(jackknife_var) 
  return(jackknife_se)
}

# Specify the parameters  
alpha1 = 0; alpha2 = 0
#tau1 = 1; tau2 = 1
tau1 = 0.1; tau2 = 0.1
sigma1 = 1; sigma2 = 1 
n1 = 200; n2 = 200
ncells = 1000 ; p = 3 
K = 200 # Number of bins 
maxIter = 30


repID = 2 # Make sure the code is reproducible 
nCPUS = 5 


#logfile = paste0("./logs/singlecell-", alpha1, "-", alpha2, "-ncell-", ncell, "-n1-", n1, "-n2-", n2, ".txt") 
#cl = makeCluster(nCPUS, type = 'SOCK', outfile = logfile ) 
#registerDoSNOW(cl)
# You can use different combine functions 
#combine = function(x, ...){
#  mapply(rbind, x, ..., SIMPLIFY = FALSE)
#}
pvMat_mom = NULL
pvMat_beta = NULL
pvMat_jack = NULL
for (iter in 1:maxIter) {
  simData = simu_hier(alpha1 = alpha1, alpha2 = alpha2, tau1 = tau1, tau2 = tau2, sigma1 = sigma1, sigma2 = sigma2, n1 = n1, n2 = n2, ncells = ncells)
  Y1 = simData$Y1
  Y2 = simData$Y2
  y1 = as.vector(c(Y1)) 
  y2 = as.vector(c(Y2))
  y_agg = c(y1, y2)
  # Carrier density 
  l = min(y_agg) #lower bound
  u = max(y_agg) #upper bound
  carrier = density(x = y_agg, kernel = c("gaussian"), n = K, from = l, to = u ) 
  est_grid = seq(l, u, length.out = K + 1) #range is not 0,1 like in the original example
  est_midpoints  = (est_grid[-1] + est_grid[-length(est_grid)])/2  
  Smat1 = matrix(NA, nrow = n1, ncol = K) #group 1 COUNTS matrix for one epoch
  Smat2 = matrix(NA, nrow = n2, ncol = K) #group 2 COUNTS for one epoch
  # Counts in each bin 
  Smat1 = get_bin_counts(Y1, bin_edges = est_grid, K)
  Smat2 = get_bin_counts(Y2, bin_edges = est_grid, K)
  # for (i in 1:n1) {
  #   hist_res = hist(Y1[i, ], breaks = est_grid, plot=FALSE)
  #   Smat1[i, ] = hist_res$counts
  # }
  # for (i in 1:n2) {
  #   hist_res = hist(Y2[i, ], breaks = est_grid, plot=FALSE)
  #   Smat2[i, ] = hist_res$counts
  # } 
  Smat = rbind(Smat1, Smat2) 
  tk = est_midpoints
  X = rep(1, length(est_midpoints))
  for (dim in 1:p){
    X = cbind(X, tk^dim)
  } 
  varb = paste0("X_", 1:p)
  colnames(X) = c("Intercept", varb) 
  # For now we consider sum first to avoid numerical issue
  S1sum = colSums(Smat1)
  S2sum = colSums(Smat2)
  carrier = density(x = y_agg, kernel = c("gaussian"), n = K, from = l, to = u )
  scale_factor = (ncells*nrow(Smat))/sum(carrier$y) 
  carrier_est = carrier$y*scale_factor # Better to sum  
  df1 = as.data.frame(cbind( S1sum, carrier_est, X))
  df2 = as.data.frame(cbind( S2sum, carrier_est, X))
  colnames(df1)[1] = "sum_cts"
  colnames(df2)[1] = "sum_cts" 
  formula = as.formula( paste0('sum_cts~offset(log(carrier_est))+', paste(varb, collapse = '+'))) 
  
  # Fit glm model for each group
  sef_gp1_sum = tryCatch(glm(
    formula,
    family = poisson(link="log"), 
    data = df1
  ))  
  sef_gp2_sum = tryCatch(glm(
    formula,
    family = poisson(link="log"), 
    data = df2
  ))   
  
  beta_sum_est1 = as.vector(sef_gp1_sum$coefficients)
  beta_sum_est2 = as.vector(sef_gp2_sum$coefficients)
  sef_sum_df1 = as.vector( carrier_est * exp(X %*% beta_sum_est1) )
  sef_sum_df2 = as.vector( carrier_est * exp(X %*% beta_sum_est2) )  
  
  # Moment estimator 
  Cov1 = mom_cov(Y1, p)
  Cov2 = mom_cov(Y2, p) 
  
  T1 = y1 
  for (dim in 2:p){
    T1 = cbind(T1, y1^dim)
  }  
  T2 = y2 
  for (dim in 2:p){
    T2 = cbind(T2, y2^dim)
  } 
  T1_bar = colMeans(T1) 
  T2_bar = colMeans(T2)
  
  Zscore_mom = (T1_bar - T2_bar)/sqrt( diag( Cov1 + Cov2 ) ) 
  pv_mom = 2*pnorm(-abs(Zscore_mom))  
  # Asymptotic Covariance for beta 
  X_dec = t( t(X) - colMeans(X))[ ,-1] 
  G_1 = t(X_dec) %*% ( sef_sum_df1 * X_dec )  
  G_2 = t(X_dec) %*% ( sef_sum_df2 * X_dec ) 
  Sand1 = G_1/(n1*ncells)
  Sand2 = G_2/(n2*ncells) 
  Sand1_inv = solve(Sand1)
  Sand2_inv = solve(Sand2) 
  Cov_beta1 = Sand1_inv%*%Cov1%*%Sand1_inv
  Cov_beta2 = Sand2_inv%*%Cov2%*%Sand2_inv 
  beta_sum_diff = (beta_sum_est1 - beta_sum_est2)[-1] 
  Zscore_beta = beta_sum_diff/sqrt(diag(Cov_beta1 + Cov_beta2) )
  pv_beta = 2*pnorm( -abs(Zscore_beta) ) 
  
  sd_jack = jackknife_se(Y1, Y2, est_grid=est_grid, K, p, ncells)
  Zscore_jack = beta_sum_diff/sd_jack[-1] 
  pv_jack = 2*pnorm(-abs(Zscore_jack) )
  
  # # Bootstrap method to evaluate the uncertainty of beta
  # B = 200
  # betaMat_diff_b = NULL
  # ZscoreMat_b = NULL
  # start.time <- Sys.time()
  # for (b in 1:B) {
  #   Y1b = NULL
  #   Y2b = NULL
  #   for (i in 1:n1 ) {
  #     Y1b = rbind(Y1b, Y1[i, ][ sample(ncol(Y1), replace = T) ] )
  #   }
  #   for (i in 1:n2 ) {
  #     Y2b = rbind(Y2b, Y2[i, ][ sample(ncol(Y2), replace = T) ] )
  #   }
  # 
  #   y1b = as.vector(c(Y1b))
  #   y2b = as.vector(c(Y2b))
  #   y_agg_b = c(y1b, y2b)
  #   Smat1_b = matrix(NA, nrow = n1, ncol = K) #group 1 COUNTS matrix for one epoch
  #   Smat2_b = matrix(NA, nrow = n2, ncol = K) #group 2 COUNTS for one epoch
  # 
  #   Smat1_b = get_bin_counts(Y1b, bin_edges = est_grid, K)
  #   Smat2_b = get_bin_counts(Y2b, bin_edges = est_grid, K)
  #   Smat_b = rbind(Smat1_b, Smat2_b)
  #   S1sum_b = colSums(Smat1_b)
  #   S2sum_b = colSums(Smat2_b)
  # 
  #   carrier_b = density(x = y_agg_b, kernel = c("gaussian"), n = K, from = l, to = u )
  #   scale_factor_b = (ncells*nrow(Smat_b))/sum(carrier_b$y)
  #   carrier_est_b = carrier_b$y*scale_factor_b # Better to sum
  #   df1_b = as.data.frame(cbind( S1sum_b, carrier_est_b, X))
  #   df2_b = as.data.frame(cbind( S2sum_b, carrier_est_b, X))
  #   colnames(df1_b)[1] = "sum_cts"
  #   colnames(df2_b)[1] = "sum_cts"
  #   formula = as.formula( paste0('sum_cts~offset(log(carrier_est))+', paste(varb, collapse = '+')))
  # 
  #   sef_gp1_sum_b = tryCatch(glm(
  #     formula,
  #     family = poisson(link="log"),
  #     data = df1_b
  #   ))
  # 
  #   sef_gp2_sum_b = tryCatch(glm(
  #     formula,
  #     family = poisson(link="log"),
  #     data = df2_b
  #   ))
  #   beta1_b = as.vector(sef_gp1_sum_b$coefficients)
  #   beta2_b = as.vector(sef_gp2_sum_b$coefficients)
  #   beta_diff_b = beta1_b - beta2_b
  #   betaMat_diff_b = rbind(betaMat_diff_b, beta_diff_b)
  # 
  #   sef_df1_b = as.vector( carrier_est_b * exp(X %*% beta1_b) )
  #   sef_df2_b = as.vector( carrier_est_b * exp(X %*% beta2_b) )
  # 
  #   X_dec = t( t(X) - colMeans(X))[ ,-1]
  #   Sand1_b = t(X_dec) %*% ( sef_df1_b * X_dec )/(n1*ncells)
  #   Sand2_b = t(X_dec) %*% ( sef_df2_b * X_dec )/(n2*ncells)
  #   Sand1_inv_b = solve(Sand1_b)
  #   Sand2_inv_b = solve(Sand2_b)
  #   Cov_beta1_b = Sand1_inv_b%*%Cov1%*%Sand1_inv_b
  #   Cov_beta2_b = Sand2_inv_b%*%Cov2%*%Sand2_inv_b
  #   Zscore_b = (beta_diff_b[-1] - beta_sum_diff)/sqrt(diag(Cov_beta1_b + Cov_beta2_b) )
  #   ZscoreMat_b = rbind(ZscoreMat_b, Zscore_b)
  #   
  # }
  # end.time <- Sys.time()
  # time.taken <- end.time - start.time
  # time.taken
  # # 1. Using the boostrap covariance matrix 
  # sd_boot = apply(betaMat_diff_b[ ,-1], 2,sd)
  # pvBoot_sd = 2*pnorm(-abs(beta_sum_diff/sd_boot ))
  # pvBoot = NULL
  # for (i in 1:ncol(ZscoreMat_b)) {
  #    pvBoot = c(pvBoot, mean( abs(ZscoreMat_b[ ,i ]) > abs(Zscore_beta[i]) ) )
  # }

  pvMat_mom = rbind(pvMat_mom, pv_mom)
  pvMat_beta = rbind(pvMat_beta, pv_beta)
  pvMat_jack = rbind(pvMat_jack, pv_jack)
}


  library(ggplot2)
gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    #geom_point(aes(expected, observed), shape = 1, size = 3) +
    geom_point(aes(expected, observed)) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po) + 
    theme_bw(base_size = 16)  + 
    theme(
      axis.ticks = element_line(size = 0.5),
      panel.grid = element_blank()
      # panel.grid = element_line(size = 0.5, color = "grey80")
    )
}


pvMat_mom
pvMat_beta
pvMat_jack

gg_qqplot(pvMat_mom[,1]) 
gg_qqplot(pvMat_mom[,2]) 
gg_qqplot(pvMat_mom[,3]) 

gg_qqplot(pvMat_beta[,1]) 
gg_qqplot(pvMat_beta[,2]) 
gg_qqplot(pvMat_beta[,3]) 

gg_qqplot(pvMat_jack[,1])
gg_qqplot(pvMat_jack[,2])
gg_qqplot(pvMat_jack[,3])



# Bootstrap method to evaluate the uncertainty of beta
B = 100 
# Function to perform bootstrap resampling and compute estimates
# For now, we only consider the percentile bootstrap and studentized Bootstrap
# Both of them are based on non-parametric bootstrap
bootstrap_analysis <- function(Y1, Y2, X, est_grid, K, l, u, B, ncells) {
  # Initialize matrices to store bootstrap results
  n1 = nrow(Y1)
  n2 = nrow(Y2)
  betaMat_diff_b <- matrix(NA, nrow = B, ncol = ncol(X) )
  ZscoreMat_b <- NULL  # This can be used if needed for z-scores
  X_dec = t( t(X) - colMeans(X))[ ,-1] # Centralized sufficient statistics 
  # Loop over the number of bootstrap iterations
  for (b in 1:B) {
    # Bootstrap resampling for each individual
    Y1b <- t(apply(Y1, 1, function(row) sample(row, replace = TRUE)))
    Y2b <- t(apply(Y2, 1, function(row) sample(row, replace = TRUE)))
    # Flatten the resampled data
    y1b <- as.vector(c(Y1b))
    y2b <- as.vector(c(Y2b))
    y_agg_b <- c(y1b, y2b)
    # Compute bin counts using the optimized method
    Smat1_b <- get_bin_counts(Y1b, bin_edges = est_grid, K)
    Smat2_b <- get_bin_counts(Y2b, bin_edges = est_grid, K)
    Smat_b <- rbind(Smat1_b, Smat2_b)
    # Compute the carrier density and scale factor
    carrier_b <- density(x = y_agg_b, kernel = "gaussian", n = K, from = l, to = u)
    scale_factor_b <- (ncells * nrow(Smat_b)) / sum(carrier_b$y)
    carrier_est_b <- carrier_b$y * scale_factor_b
    # Prepare data frames for Poisson regression
    df1_b <- data.frame(sum_cts = colSums(Smat1_b), carrier_est_b, X)
    df2_b <- data.frame(sum_cts = colSums(Smat2_b), carrier_est_b, X)
    # Construct the formula for GLM
    formula <- as.formula(paste0("sum_cts~offset(log(carrier_est_b))+", paste(colnames(X)[-1], collapse = "+")))
    # Fit Poisson GLM for both groups and compute coefficient differences
    beta1_b <- tryCatch(coef(glm(formula, family = poisson(link = "log"), data = df1_b)), error = function(e) rep(NA, ncol(X)))
    beta2_b <- tryCatch(coef(glm(formula, family = poisson(link = "log"), data = df2_b)), error = function(e) rep(NA, ncol(X)))
    beta_diff_b <- beta1_b - beta2_b
    # Store the results
    betaMat_diff_b[b, ] <- beta_diff_b
    
    sef_df1_b = as.vector( carrier_est_b * exp(X %*% beta1_b) )
    sef_df2_b = as.vector( carrier_est_b * exp(X %*% beta2_b) )  
    
    Cov1b = mom_cov(Y1b, p)
    Cov2b = mom_cov(Y2b, p) 
    Sand1_b = t(X_dec) %*% ( sef_df1_b * X_dec )/(n1*ncells)
    Sand2_b = t(X_dec) %*% ( sef_df2_b * X_dec )/(n2*ncells)
    Sand1_inv_b = solve(Sand1_b)
    Sand2_inv_b = solve(Sand2_b)
    Cov_beta1_b = Sand1_inv_b%*%Cov1b%*%Sand1_inv_b
    Cov_beta2_b = Sand2_inv_b%*%Cov2b%*%Sand2_inv_b
    Zscore_b = (beta_diff_b[-1] - beta_sum_diff)/sqrt(diag(Cov_beta1_b + Cov_beta2_b) )
    # Percentile t method 
    ZscoreMat_b = rbind(ZscoreMat_b, Zscore_b) 
  }
  
  return(list(betaMat_diff = betaMat_diff_b, ZscoreMat_b = ZscoreMat_b))
}

# Usage

start.time <- Sys.time()
result <- bootstrap_analysis(Y1, Y2, X, est_grid, K, l, u, B, ncells)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
# Percentile method  
betaMat_diff_percent = result$betaMat_diff[ ,-1]
betaMat_diff_percent = betaMat_diff_percent[!rowSums( is.na(betaMat_diff_percent)), ]
beta_sum_diff 
pv_percent = NULL 
for (i in 1:ncol(betaMat_diff_percent) ) {
  pv_percent = c(pv_percent, mean(abs(betaMat_diff_percent[ ,i] - beta_sum_diff[i]) > abs(beta_sum_diff[i]) ) )
}
pv_percent
