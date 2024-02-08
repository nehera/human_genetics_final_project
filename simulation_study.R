start_time <- Sys.time()

# Load packages
if (("pacman" %in% installed.packages()[,"Package"]) == FALSE) { install.packages("pacman") }
pacman::p_load(tidyverse, MASS, BIPnet, GraceAKO, glmnet, caret, parallel) 

# Define simulation parameters
M <- 30 # Number of simulation replicates
seeds <- 1:M
possible_n <- c(200, 400)
possible_p <- c(400, 600)
possible_gy1 <- c(4, 2)
experiments_df <- expand.grid(s=seeds, n= possible_n, p= possible_p, gy1= possible_gy1)

experiments_list <- split(experiments_df, seq(nrow(experiments_df)))
experiments_list <- lapply(experiments_list, as.numeric)

# Input: vector of s seed, n observations, p predictors, and gy1 groups important to y in view 1
# Output: data.frame agg_results, which includes experimental conditions, method, FDR, TPP, & n features selected

# experiments_list <- experiments_list[1:2]
# s_n_p_gy1 <- experiments_list[[1]]

# Define main function
simulate_estimate <- function(s_n_p_gy1) {
  
  s = s_n_p_gy1[1] 
  n = s_n_p_gy1[2] # Number of observations
  p = s_n_p_gy1[3] # Number of features in total
  gy1 = s_n_p_gy1[4]
  
  ## Simulate data
  set.seed(s)
  
  # Fix additional parameters
  p_m <- floor(p/ 2) # p_m features per view
  n_groups <- 10 # N groups
  group_size <- 10 # N features/ group
  n_primary <- n_groups # Number of primary features (1/ group)
  n_secondary_per_group <- group_size-1 # N secondary features/ group
  # Define group distribution across views
  groups_per_view <- floor(n_groups/ 2) 
  n_primary_per_view <- floor(n_primary/ 2)
  n_secondary_per_view <- groups_per_view*group_size - n_primary_per_view
  n_singletons_per_view <- p_m - n_primary_per_view - n_secondary_per_view
  n_important_groups <- 4 # Number of groups that impact Y
  gy2 <- n_important_groups - gy1
  # Define intra-group correlations
  corr_primary_secondary <- 0.4
  corr_secondary_secondary <- corr_primary_secondary^2
  mu <- matrix(0, nrow = group_size)
  group_Sigma <- matrix(corr_secondary_secondary, group_size, group_size) +
    diag(1-corr_secondary_secondary, group_size)
  group_Sigma[-1, 1] <- corr_primary_secondary
  group_Sigma[1, -1] <- corr_primary_secondary
  # Observe groups
  observe_group <- function(g, n, mu, group_Sigma) {
    # TODO Potentially anti-correlate the group's features
    # sgn_flip <- rbinom(1,1,.5)
    # if (sgn_flip==1) { group_Sigma <- -group_Sigma+2*diag(nrow(group_Sigma)) }
    return(mvrnorm(n, mu, group_Sigma))
  }
  group_results <- lapply(1:n_groups, observe_group, 
                          n, mu, group_Sigma)
  
  # Split groups evenly across views
  X1_primary <- do.call(cbind, group_results[1:groups_per_view])
  X2_primary <- do.call(cbind, group_results[(groups_per_view+1):n_groups])
  
  # Observe singletons
  X1_singletons <- matrix(rnorm(n*n_singletons_per_view), 
                          ncol = p_m - ncol(X1_primary), nrow = n)
  X2_singletons <- matrix(rnorm(n*n_singletons_per_view), 
                          ncol = p_m - ncol(X2_primary), nrow = n)
  
  # Combine groups and singletons
  X1 <- cbind(X1_primary, X1_singletons)
  X2 <- cbind(X2_primary, X2_singletons)
  X <- cbind(X1, X2)
  
  # Define true effects
  group_beta_vector <- c(5, rep(5/sqrt(group_size), n_secondary_per_group))

  # Split effects across views by direction and magnitude
  beta1 <- rep(group_beta_vector, gy1)
  n_unimportant_1 <- p_m - length(beta1) # Assume constant across views
  beta1 <- c(beta1, rep(0, n_unimportant_1))
  
  beta2 <- rep(group_beta_vector, gy2)
  n_unimportant_2 <- p_m - length(beta2) # Assume constant across views
  beta2 <- c(beta2, rep(0, n_unimportant_2))
  
  beta_true <- matrix(c(beta1, beta2), ncol = 1)
  active_features <- which(beta_true!=0)
  
  # Simulate outcome
  epsilon <- matrix(rnorm(n), ncol = 1)
  Y <- X %*% beta_true + epsilon
  
  ## Fit models
  
  ## GRACE-AKO and GRACE analysis
  
  # Initialize L and W
  L <- W <- matrix(0, nrow = p, ncol = p) 
  # Generate W
  block_vector <- rep(1, group_size)
  block_matrix <- kronecker(diag(groups_per_view), 
                            matrix(1, nrow = group_size, ncol = group_size))
  # Define X1 and X2's groups
  n_grouped_features_per_view <- (n_primary_per_view + n_secondary_per_view)
  W[1:n_grouped_features_per_view, 1:n_grouped_features_per_view] <- block_matrix
  W[(p_m+1):(p_m+n_grouped_features_per_view), 
    (p_m+1):(p_m+n_grouped_features_per_view)] <- block_matrix
  # Compute D
  D <- rowSums(W)
  # Compute L
  for (u in 1:p) {
    for (v in 1:p) {
      if (u==v & D[u]!=0) {
        L[u,v] <- 1-W[u,v]/D[u]
      } else if (W[u,v]!=0) {
        # Probability of connection is nonzero
        L[u,v] <- -W[u,v]/ sqrt(D[u]*D[v])
      }
      # 0 otherwise
    }
  }
  lambda.L = seq(0.1, 2, 0.5)
  lambda.1 =  seq(110, 150, 10)
  lambda.2 = seq(1, 10, 5)
  fdr = 0.1
  n_bootstraps = 25
  gamma = 0.1
  L_1 = build_L1(L)
  X = as.data.frame(X)
  std.X = X %>% mutate_all(~scale((.)%>% as.vector))
  Y = data.frame(Y)
  std.Y = Y %>% mutate_all(~scale((.)%>% as.vector))
  X = as.matrix(std.X)
  Y1 = as.vector(std.Y)$Y
  Y = as.vector(Y1)
  
  # Conduct Grace-AKO and Grace
  result <- Grace_multi_knockoff(X = X, Y = Y, L_1 = L_1, L = L, fdr = fdr, 
                                lambda.L = lambda.L, lambda.1 = lambda.1, 
                                lambda.2= lambda.2, n_bootstraps = n_bootstraps, gamma = gamma)
  
  # Extract summary statistics
  mfdr_tpp_GRACE = mfdr_tpp_base(q=fdr, beta_true, result$beta_est_baseline)
  n_selected_GRACE <- sum(result$beta_est_baseline!=0)
  mfdr_tpp_GRACE_AKO = mfdr_tpp_knockoff(q = fdr, result$w, result$t, beta_true, result$beta_est_knockoff)
  n_selected_GRACE_AKO <- sum(result$beta_est_knockoff!=0)
  
  summary_grace <- data.frame(method=c("GRACE", "GRACE_AKO"), 
                              FDR=c(mfdr_tpp_GRACE$mFDR_baseline, mfdr_tpp_GRACE_AKO$mFDR_withknockoff),
                              TPP=c(mfdr_tpp_GRACE$TPP_baseline, mfdr_tpp_GRACE_AKO$TPP_withknockoff),
                              n_features=c(n_selected_GRACE, n_selected_GRACE_AKO))

  ## BIPnet and BIP analysis
  Y = matrix(Y, ncol = 1) # Ensure Y in matrix form
  dataList=list(X1, X2, Y) 
  Path1 <- Path2 <- matrix(0, nrow = p_m, ncol = n_groups+1)
  # Define singleton grouping
  Path1[(n_grouped_features_per_view+1):p_m, 11] <- 1
  Path2[(n_grouped_features_per_view+1):p_m, 11] <- 1
  # Define correlated var grouping
  PathDesign <- kronecker(diag(n_groups), matrix(1, group_size))
  Path1[1:n_grouped_features_per_view, -11] <- PathDesign[1:n_grouped_features_per_view, ]
  Path2[1:n_grouped_features_per_view, -11] <- PathDesign[(n_grouped_features_per_view+1):(n_groups*group_size), ]
  PathList=list(Path1, Path2)
  BAnet = BIP(dataList=dataList, IndicVar=c(0,0,1), groupList=PathList, Method="BIPnet", 
         nbrcomp=4, sample=5000, burnin=1000)
  BA = BIP(dataList=dataList, IndicVar=c(0,0,1), Method="BIP", 
         nbrcomp=4, sample=5000, burnin=1000)
  

  summarize_BA <- function(BIP_fit, beta_true) {
    important_feature_index <- which(beta_true!=0)
    n_important <- length(important_feature_index)
    var_sel_global_post <- c(BIP_fit$VarSelMeanGlobal[[1]], BIP_fit$VarSelMeanGlobal[[2]])
    positives <- which(var_sel_global_post>0.5)
    n_TP <- sum(positives %in% important_feature_index)
    n_FP <- sum(!(positives %in% important_feature_index))
    FDR <- n_FP/ (n_TP+n_FP)
    TPP <- n_TP/ length(positives)
    n_selected <- length(positives)
    BIP_fit_summary <- c(FDR, TPP, n_selected)
    return(BIP_fit_summary)
  }
  
  summary_BIP <- data.frame(matrix(ncol = 4, nrow = 2))
  summary_BIP[, 1] <- c("BIP", "BIPnet")
  summary_BIP[1, -1] <- summarize_BA(BA, beta_true)
  summary_BIP[2, -1] <- summarize_BA(BAnet, beta_true)
  colnames(summary_BIP) <- c("method", "FDR", "TPP", "n_features")
  
  ## Elastic net analysis
  cv_10 = trainControl(method = "cv", number = 10)
  hit_elnet = train(
    Y~ ., data = data.frame(Y, X),
    method = "glmnet",
    trControl = cv_10
  )
  # Identify the best values by minimum RMSE
  best_id <- hit_elnet$results$RMSE %>% which.min
  best_alpha <- hit_elnet$results$alpha[best_id]

  # Fit with best_alpha
  cv_enet <- cv.glmnet(x = X, y = Y, alpha = best_alpha) 
  
  important_features <- which(beta_true!=0)
  beta <- coef(cv_enet)
  positives <- which(beta!=0)[-1] - 1 
  n_TP <- sum(positives %in% important_features)
  n_FP <- sum(!(positives %in% important_features))
  FDR <- n_FP/ (n_TP+n_FP)
  TPP <- n_TP/ length(positives)
  n_selected <- length(positives)
  
  summary_elnet <- data.frame(matrix(ncol = 4, nrow = 1))
  summary_elnet[, 1] <- "elnet"
  summary_elnet[1, -1] <- c(FDR, TPP, n_selected)
  colnames(summary_elnet) <- c("method", "FDR", "TPP", "n_features")
  
  ## Consolidate summary statistics
  results <- rbind(summary_grace, summary_BIP, summary_elnet)
  results$s <- s
  results$n <- n
  results$p <- p
  results$gy1 <- gy1
  
  return(results)
  
}

# Define how many cores to use
n_experiments <- length(experiments_list)
n_cores <- detectCores()
if (n_experiments <= n_cores) {
  job_cores <- n_experiments
} else { job_cores <- n_cores }

# Perform the simulation
agg_results_list <- mclapply(experiments_list, simulate_estimate, mc.cores = job_cores)

# Reshape results
agg_results <- do.call(rbind, agg_results_list)

# Save results as an .rds file
saveRDS(agg_results, file = "agg-results.rds")

end_time <- Sys.time()
simulation_duration <- end_time - start_time
print("Data simulation, model fitting, and summary statistic evaluation required:")
print(simulation_duration)