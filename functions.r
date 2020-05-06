
logit <- function(x){
  1 / (1 + exp(-x))
}

invLogit <- function(x){
  log(x / (1 - x))
}

ExtractCovariatesFromText <- function(covariates_text) {
  
  covariates <- sapply(strsplit(covariates_text, split = ",")[[1]], function(x){
    if(grepl("-",x)){
      as.numeric(seq(strsplit(x, split = "-")[[1]][1], strsplit(x, split = "-")[[1]][2]))
    } else {
      as.numeric(x)
    }
  })
  covariates <- as.vector(unlist(covariates))
  if(covariates[1] == 0) covariates <- c()
  
  covariates
}

# true or false if the covaraites is categorical or numerical
ExtractClassCovariates <- function(ncov, fac_covariates, column_covariate){
  classCovariates <- rep(T, ncov)
  classCovariates[match(fac_covariates, column_covariate)] <- F
  
  classCovariates
}

Extract_IndexesCovariates <- function(X, column_covariate, ncov, classCovariates) {
  
  indexes_covariates <- c()
  indexes_covariates[1] <- 1
  if(any(column_covariate != 0)){
    
    indexes_covariates <- c()
    indexes_covariates[1] <- 1
    k <- 2
    for (i in 1:ncov) {
      if(classCovariates[i]){
        indexes_covariates[k] <- i + 1
        k <- k + 1
      } else {
        num_levels <- length(unique(X[,i]))
        indexes_covariates[k + 0:(num_levels-2)] <- i + 1
        k <- k + num_levels - 1
      }
    }
    
  }
  
  indexes_covariates
}

logit <- function(x){
  1 / (1 + exp(-x))
}

invLogit <- function(x){
  log(x / (1 - x))
}

sample_z <- function(psi, p, k_s){
  
  z <- rep(NA, nrow(k_s))
  
  for (i in 1:length(z)) {
    
    if(k_s[i,3] == 1){
      
      z[i] <- 1
      
    } else {
      
      p_zsequal1 <- (psi[i] * prod(1 - p[i + seq_len(k_s[i,4])])) / (psi[i] * prod(1 - p[i + seq_len(k_s[i,4])]) + (1 - psi[i])) 
      
      z[i] <- rbinom(1, 1, p_zsequal1)
      
    }
    
  }
  
  z
}

update_psi <- function(beta_psi, X_psi, 
                       b_psi, B_psi, d_bar_psi, 
                       z, k_s, sites, Y, ncov_psi, X_y_index, indexes_covariates_psi,
                       a_s, c_s, sigma_a, alpha_dp, M_neal,
                       usingClustering){
  
  k <- z - .5
  n <- rep(1, length(k))
  
  list_beta_psi <- sampler_beta(beta_psi, a_s, X_psi, b_psi, B_psi,
                                n, k, Y, ncov_psi, X_y_index)
  
  beta_psi <- list_beta_psi$beta
  Omega <- list_beta_psi$Omega
  XtOmegaX <- list_beta_psi$XtOmegaX
  
  # sample clustering variables
  
  if(usingClustering){
    
    list_a_s <- sample_as_cpp_clustering(c_s, k_s$Site, sites, beta_psi, X_psi, 
                                         z, Omega, sigma_a, alpha_dp, M_neal)
    a_s <- list_a_s$a_s
    c_s <- list_a_s$c_s  
    
    # relabel clusters
    {
      c_s <- factor(c_s)
      levels(c_s) <- as.character(1:length(unique(c_s)))
      c_s <- as.numeric(c_s)
    }
      
  }
  
  psi <- as.vector(logit(X_psi %*% beta_psi + a_s))
  
  list("beta" = beta_psi, "psi" = psi,
       "a_s" = a_s, "c_s" = c_s)
}

update_psi_integrated <- function(beta_psi, X_psi, 
                       b_psi, B_psi, d_bar_psi, 
                       z, k_s, sites, Y, ncov_psi, X_y_index, indexes_covariates_psi,
                       a_s, c_s, sigma_a, alpha_dp, M_neal,
                       mu_clusters, sigma_clusters,
                       usingClustering){
  
  k <- z - .5
  n <- rep(1, length(k))
  
  list_beta_psi <- sampler_beta_integrated(beta_psi, a_s, mu_clusters, sigma_clusters,
                                           X_psi, b_psi, B_psi,
                                           n, k, Y, ncov_psi, X_y_index)
  
  beta_psi <- list_beta_psi$beta
  Omega <- list_beta_psi$Omega
  XtOmegaX <- list_beta_psi$XtOmegaX
  
  # sample clustering variables
  
  if(usingClustering){
    
    list_a_s <- sample_as_cpp_clustering_integrated(c_s, k_s$Site, sites, beta_psi, X_psi, 
                                                    z, Omega, sigma_a, alpha_dp, M_neal)
    a_s <- list_a_s$a_s
    c_s <- list_a_s$c_s  
    mu_clusters <- list_a_s$mu_clusters  
    sigma_clusters <- list_a_s$sigma_clusters  
    
    mu_clusters <- a_s#mu_clusters[c_s[k_s$Site]]
    sigma_clusters <- sigma_clusters[c_s[k_s$Site]]
    
    # relabel clusters
    {
      c_s <- factor(c_s)
      levels(c_s) <- as.character(1:length(unique(c_s)))
      c_s <- as.numeric(c_s)
    }
    
  }
  
  psi <- as.vector(logit(X_psi %*% beta_psi + a_s))
  
  list("beta" = beta_psi, "psi" = psi,
       "a_s" = a_s, "c_s" = c_s,
       "mu_clusters" = mu_clusters, "sigma_clusters" = sigma_clusters)
}

update_l <- function(l_gp, a_l, b_l, beta_psi, sigma_gp, sd_l, Y, B_psi, l_beta_proposal, index_l, l_all){
  
  l_gp <- sample_l_adaptive(l_gp, a_l, b_l, beta_psi[2:(Y+1)], 
                            sigma_gp, sd_l, l_beta_proposal, index_l, l_all)
  B_psi[2:(Y+1), 2:(Y+1)] <- K(1:Y, 1:Y, sigma_gp, l_gp)
  
  list("l_gp" = l_gp, "B_psi" = B_psi)
}

update_p <- function(beta_p, y, z, X_p, indexes_covariates_p, b_p, B_p){
  
  k <- y - .5
  n <- rep(1, length(k))
  
  k <- k[z==1]
  n <- n[z==1]
  X_p_present <- X_p[z == 1,,drop = FALSE]
  
  list_beta_p <- sample_beta_omega_cpp(beta_p, X_p_present, b_p, B_p, n, k)
  beta_p <- list_beta_p$beta
  p <- as.vector(logit(X_p %*% beta_p))
  
  list("p" = p, "beta" = beta_p)
}

sample_l_simple <- function(l, a_l, b_l, b_t, sigma_gp, sd_l){
  
  Y <- length(b_t)
  l_star <- rnorm(1, l, sd_l)
  
  if(l_star > 0){
    logPrior <- dgamma(l, a_l, b_l, log = T) 
    logPrior_star <- dgamma(l_star, a_l, b_l, log = T) 
    
    loglikelihood <- dmvnorm(b_t, sigma =  K(1:Y, 1:Y, sigma_gp, l), log = T)
    loglikelihood_star <- dmvnorm(b_t, sigma =  K(1:Y, 1:Y, sigma_gp, l_star), log = T)
    
    logposterior <- logPrior + loglikelihood
    logposterior_star <- logPrior_star + loglikelihood_star
    
    if(runif(1) < exp(logposterior_star - logposterior)){
      l <- l_star
    }
  }
  
  l
}

sample_l_adaptive <- function(l, a_l, b_l, b_t, sigma_gp, sd_l, l_beta_proposal, index_l, l_all){
  
  Y <- length(b_t)
  
  if(index_l < 200){
    proposal_l <- sd_l  
  } else {
    Sigma_l <- (2.38) * sd(l_all[1:(index_l - 1)])
    proposal_l <- (1 - l_beta_proposal) * Sigma_l + l_beta_proposal * sd_l
  }
  
  l_star <- rnorm(1, l, proposal_l)
  
  if(l_star > 0){
    logPrior <- dgamma(l, a_l, b_l, log = T) 
    logPrior_star <- dgamma(l_star, a_l, b_l, log = T) 
    
    loglikelihood <- dmvnorm(b_t, sigma =  K(1:Y, 1:Y, sigma_gp, l), log = T)
    loglikelihood_star <- dmvnorm(b_t, sigma =  K(1:Y, 1:Y, sigma_gp, l_star), log = T)
    
    logposterior <- logPrior + loglikelihood
    logposterior_star <- logPrior_star + loglikelihood_star
    
    if(runif(1) < exp(logposterior_star - logposterior)){
      l <- l_star
    }
  }
  
  l
}

computeClusteringStartingPoints <- function(z, beta_psi, k_s, ncov_psi, X_psi, X_y_index, b_psi, B_psi, Y, S,
                                            sigma_a, numOfClusters){
  
  k <- z - .5
  n <- rep(1, length(k))
  
  list_beta_psi <- sampler_beta(beta_psi, rep(0, nrow(k_s)), X_psi, b_psi, B_psi,
                                n, k, Y, ncov_psi, X_y_index)
  
  Omega <- list_beta_psi$Omega
  
  a_s_clusters <- sample_cluster_locations(1:S, k_s$Site, sort(unique(k_s$Site)),
                                  beta_psi, X_psi,
                                  z, Omega, sigma_a)
  
  quantiles_as <- quantile(a_s_clusters, probs = seq(0, 1, length.out = numOfClusters + 1))
  a_s_groups <- (quantiles_as[-1] + quantiles_as[-length(quantiles_as)]) / 2
  
  quantiles_as[1] <- quantiles_as[1] - 1
  
  c_s <- cut(a_s_clusters, as.numeric(quantiles_as), labels = 1:numOfClusters)
  
  a_s <- a_s_groups[c_s[k_s$Site]]
  
  list("a_s" = a_s,
       "c_s" = as.numeric(c_s))
}

# final model functions

cleanData <- function(data, index_year, index_site, index_occ, 
                      num_covariates_psi_text, fac_covariates_psi_text,
                      num_covariates_p_text, fac_covariates_p_text){
  
  colnames(data)[c(index_year, index_site, index_occ)] <- c("Year","Site","Occ")
  
  data <- data %>% arrange(Site, Year)
  
  S <- length(unique(data$Site))
  Y <- length(unique(data$Year))
  years <- sort(unique(data$Year))
  
  # define data for occupancies
  
  {
    
    k_s <- data %>% group_by(Year, Site) %>% summarise(Occupied = as.numeric(sum(Occ) > 0),
                                                       Visits = n()) %>% arrange(Site, Year) 
    k_s <- as.data.frame(k_s)
    
    originalsites <- k_s$Site
    k_s$Site <- as.numeric(as.factor(k_s$Site))  
    
    # match row in data_p to row in k_s
    indexes_occobs <- rep(1:nrow(k_s), k_s$Visits)
    
    Occs <- data$Occ
  }
  
  # define data for site covariates
  
  {
    # covariate psi
    {
      num_covariates_psi <- ExtractCovariatesFromText(num_covariates_psi_text)
      fac_covariates_psi <- ExtractCovariatesFromText(fac_covariates_psi_text)
      
      column_covariate_psi <- sort(c(num_covariates_psi,fac_covariates_psi))
    }
    
    data_psi <- data[,c(which("Site" == colnames(data)),
                        which("Year" == colnames(data)),
                        column_covariate_psi)]
    data_psi <- data_psi[!duplicated(data_psi[,c("Site","Year")]),]
    
    data_psi <- data_psi %>% arrange(Site, Year)
    
    index_numericalCovariates_psi <- match(colnames(data)[num_covariates_psi],colnames(data_psi))
    index_categoricalCovariates_psi <- match(colnames(data)[fac_covariates_psi],colnames(data_psi))
    
  }
  
  # covariates for psi
  
  {
    
    X_psi <- data_psi
    X_psi$Year <- as.factor(X_psi$Year)
    
    nameVariables_psi <- colnames(X_psi)[-c(1,2)]
    ncov_psi <- ncol(X_psi) - 2
    
    indexes_covariates_psi <- c()
    indexes_covariates_psi[1] <- 1
    k <- 2
    indexes_covariates_psi[k + 0:(Y-1)] <- 2
    k <- k + Y 
    
    if(ncov_psi > 0){
      numOfLevels_psi <- rep(NA, ncov_psi)
      for (i in seq_len(ncov_psi)) {
        if((i + 2) %in% index_numericalCovariates_psi){
          indexes_covariates_psi[k] <- i + 2
          k <- k + 1
          X_psi[,i+2] <- as.numeric(X_psi[,i+2])
          numOfLevels_psi[i] <- 1
        } else {
          num_levels <- length(unique(X_psi[,i+2]))
          indexes_covariates_psi[k + 0:(num_levels-2)] <- i + 2
          k <- k + num_levels - 1
          X_psi[,i+2] <- as.factor(X_psi[,i+2])
          numOfLevels_psi[i] <- length(unique(X_psi[,i+2]))
        }
      }
    } else {
      numOfLevels_psi <- NULL
    }
    
    X_psi_year <- model.matrix( ~ . - 1, data = X_psi[,c("Year"),drop = F])
    
    if(ncov_psi > 0){
      X_psi_cov <- model.matrix(~ ., data = X_psi[,-c(1,2)])[,-1]
      X_psi <- cbind(1,X_psi_year, X_psi_cov)
    } else {
      X_psi <- cbind(1,X_psi_year)
    }
    
    X_y <- X_psi[,2 + 0:(Y-1)]
    X_y_index <- apply(X_y, 1, function(x) {which(x != 0)})
    X_y_index <- unlist(X_y_index)
    
  }
  
  # covariates for p
  
  {
    data_p <- data
    
    num_covariates_p <- ExtractCovariatesFromText(num_covariates_p_text)
    fac_covariates_p <- ExtractCovariatesFromText(fac_covariates_p_text)
    
    column_covariate_p <- sort(c(num_covariates_p,fac_covariates_p))
    
    X_p <- data_p[,column_covariate_p, drop = F]
    index_numericalCovariates_p <- match(colnames(data)[num_covariates_p],colnames(X_p))
    index_categoricalCovariates_p <- match(colnames(data)[fac_covariates_p],colnames(X_p))
    
    ncov_p <- ncol(X_p)
    
    if(ncov_p > 0){
      nameVariables_p <- colnames(data_p)[column_covariate_p]  
    } else {
      nameVariables_p <- c()
    }
    
    indexes_covariates_p <- c()
    indexes_covariates_p[1] <- 1
    k <- 2
    if(ncov_p > 0){
      numOfLevels_p <- rep(NA, ncov_p)
      for (i in seq_len(ncov_p)) {
        if(i %in% index_numericalCovariates_p){
          indexes_covariates_p[k] <- i + 1
          k <- k + 1
          X_p[,i] <- as.numeric(X_p[,i])
          numOfLevels_p[i] <- 1
        } else {
          num_levels <- length(unique(X_p[,i]))
          indexes_covariates_p[k + 0:(num_levels-2)] <- i + 1
          k <- k + num_levels - 1
          X_p[,i] <- as.factor(X_p[,i])
          numOfLevels_p[i] <- length(unique(X_p[,i]))
        }
      }
      
    } else {
      numOfLevels_p <- NULL
    }
    
    
    if(ncov_p == 0){
      X_p <- matrix(1, nrow = nrow(data_p), ncol = 1)
    } else {
      X_p <- model.matrix(~ ., X_p)  
    }
    
  }
  
  list("Occs" = Occs, "k_s" = k_s, "years" = years,
       "X_psi" = X_psi, "X_y_index" = X_y_index, "Y" = Y, "S" = S, "indexes_occobs" = indexes_occobs,
       "indexes_covariates_psi" = indexes_covariates_psi, "numOfLevels_psi" = numOfLevels_psi, "ncov_psi" = ncov_psi,
       "nameVariables_psi" = nameVariables_psi, "nameVariables_p" = nameVariables_p,
       "index_numericalCovariates_psi" = index_numericalCovariates_psi, "index_numericalCovariates_p" = index_numericalCovariates_p,
       "X_p" = X_p, "indexes_covariates_p" = indexes_covariates_p, "numOfLevels_p" = numOfLevels_p, "ncov_p" = ncov_p) 
}

createPrior <- function(prior_psi, sigma_psi, phi_psi,
                        sigma_gp, a_l_gp, b_l_gp,
                        prior_p, sigma_p, phi_p,
                        Y, X_psi, ncov_psi, index_numericalCovariates_psi, numOfLevels_psi,
                        X_p, ncov_p, index_numericalCovariates_p, numOfLevels_p){
  
  l_gp <- a_l_gp / b_l_gp
  
  mu_psi <- invLogit(prior_psi)
  mu_p <- invLogit(prior_p)
  
  # priors on psi
  
  {
    b_psi <- c(mu_psi, rep(0, ncol(X_psi) - 1))
    B_psi <- matrix(0, nrow = ncol(X_psi), ncol = ncol(X_psi))
    B_psi[1, 1] <- sigma_psi
    
    if(ncov_psi > 0){
      
      C <- matrix(0, nrow = ncol(X_psi) - 1, ncol = ncol(X_psi) - 1)
      C[1:Y, 1:Y] <- K(1:Y,1:Y, sigma_gp, l_gp)
      
      l <- Y 
      for (i in seq_len(ncov_psi)) {
        if((i + 2) %in% index_numericalCovariates_psi){
          C[l + 1, l + 1] <- 1
          l <- l + 1
        } else {
          num_levels <- numOfLevels_psi[i]
          C[l + 1:(num_levels-1), l + 1:(num_levels-1)] <- .5 * (diag(1, nrow = (num_levels-1)) + 
                                                                   matrix(1, nrow = (num_levels-1), ncol = (num_levels-1)))
          l <- l + num_levels - 1
        }
      }
      
      B_psi[2:(ncol(X_psi)), 2:(ncol(X_psi))] <- C
      B_psi[(Y + 2):ncol(X_psi), (Y + 2):ncol(X_psi)] <- B_psi[(Y + 2):ncol(X_psi), (Y + 2):ncol(X_psi)] * phi_psi
      
    } else {
      
      C <- matrix(0, nrow = ncol(X_psi) - 1, ncol = ncol(X_psi) - 1)
      C[1:Y, 1:Y] <- K(1:Y,1:Y, sigma_gp, l_gp)
      
      B_psi[2:(ncol(X_psi)), 2:(ncol(X_psi))] <- phi_psi * C
      
    }  
  }
  
  # prior on p
  
  {
    b_p <- c(mu_p, rep(0, ncol(X_p) - 1))
    B_p <- matrix(0, nrow = ncol(X_p), ncol = ncol(X_p))
    B_p[1, 1] <- sigma_p
    
    if(ncov_p > 0){
      
      C <- matrix(0, nrow = ncol(X_p) - 1, ncol = ncol(X_p) - 1)
      
      l <- 0
      for (i in seq_len(ncov_p)) {
        if((i) %in% index_numericalCovariates_p){
          C[l + 1, l + 1] <- 1
          l <- l + 1
        } else {
          num_levels <- numOfLevels_p[i]
          C[l + 1:(num_levels-1), l + 1:(num_levels-1)] <- .5 * (diag(1, nrow = (num_levels-1)) + 
                                                                   matrix(1, nrow = (num_levels-1), ncol = (num_levels-1)))
          l <- l + num_levels - 1
        }
      }
      
      B_p[2:(ncol(X_p)), 2:(ncol(X_p))] <- phi_p * C
      
    }   
  }
 
  list("b_psi" = b_psi, "B_psi" = B_psi,
       "b_p" = b_p, "B_p" = B_p) 
}

runMCMC <- function(Occs, X_psi, indexes_covariates_psi, ncov_psi, X_y_index,
                    S, Y, k_s, indexes_occobs,
                    X_p, indexes_covariates_p, ncov_p,
                    b_psi, B_psi, 
                    b_p, B_p,
                    sigma_gp, sd_l, l_beta_proposal, a_l_gp, b_l_gp, 
                    usingClustering, sigma_a, alpha_dp, M_neal,
                    nchain, nburn, niter){
  
  # initialize output
  {
    # z_output <- array(0, dim = c(nchain, niter, S, Y))
    
    # psi_output <- array(NA,  dim = c(nchain, niter, S, Y))
    beta_psi_output <- array(NA, dim = c(nchain , niter, length(indexes_covariates_psi)),
                             dimnames = list(c(), c(), colnames(X_psi)))
    
    l_gp_psi_output <- matrix(NA, nrow = nchain, ncol = niter)
    
    maxClusters <- 100
    n_s_output <- array(NA, dim = c(nchain, niter, maxClusters))
    # a_s_output <- array(NA, dim = c(nchain, niter, nrow(k_s)))
    
    # p_output <- array(NA,  dim = c(nchain, niter, Y, S, V))
    beta_p_output <- array(NA, dim = c(nchain , niter, length(indexes_covariates_p)),
                           dimnames = list(c(), c(), colnames(X_p)))
    
  }
  
  for (chain in 1:nchain) {
    
    # initialize parameters
    {
      
      print("Initializing the parameters")
      
      # covariates for psi
      {
        beta_psi <- b_psi
        psi <- as.vector(logit(X_psi %*% beta_psi))
        
      }
      
      # covariates for p
      {
        beta_p <- b_p
        p <- as.vector(logit(X_p %*% beta_p))
      }
      
      z <- rep(NA, nrow(k_s))
      for (i in 1:nrow(k_s)) {
        if(k_s[i,1] == 1){
          z[i] <- 1
        } else {
          z[i] <- rbinom(1, 1, psi[i])
        }
      }
      # z but same length of occs
      z_all <- z[indexes_occobs]
      
      # clustering parameters
      {
        if(usingClustering){
          list_cs <- computeClusteringStartingPoints(z, beta_psi, k_s, ncov_psi, X_psi, X_y_index, b_psi, B_psi, Y, S,
                                                     sigma_a, numOfClusters = 30)
          c_s <- list_cs$c_s
          a_s <- list_cs$a_s  
          # mu_clusters <- rep(0, S)
          # sigma_clusters <- rep(.01, S)
        } else {
          c_s <- rep(1, S)
          a_s <- rep(0, nrow(k_s))
        }
        
      }
      
      l_gp <- a_l_gp / b_l_gp
      
      l_all <- rep(NA, nburn + niter)
      index_l <- 1
      
    }
    
    for (iter in seq_len(nburn + niter)) {
      
      if(iter <= nburn){
        print(paste0("Chain = ",chain," - Burn-in Iteration = ",iter))  
        # print(paste0("Number of clusters = ",length(unique(c_s))))
      } else {
        print(paste0("Chain = ",chain," - Iteration = ",iter - nburn))  
        # print(paste0("Number of clusters = ",length(unique(c_s))))
      }
      
      # sample z ----
      
      z <- sample_z_cpp(psi, p, as.matrix(k_s[,3:4]))
      z_all <- z[indexes_occobs]
      
      # sample psi ----
      
      list_psi <- update_psi(beta_psi, X_psi, 
                             b_psi, B_psi, d_bar_psi, 
                             z, k_s, unique(k_s$Site), Y, ncov_psi, X_y_index, indexes_covariates_psi,
                             a_s, c_s, sigma_a, alpha_dp, M_neal, usingClustering) 
      psi <- list_psi$psi
      beta_psi <- list_psi$beta
      a_s <- list_psi$a_s
      c_s <- list_psi$c_s
      
      # sample l ---------
      
      list_l <- update_l(l_gp, a_l_gp, b_l_gp, beta_psi, sigma_gp, sd_l, Y,
                         B_psi, l_beta_proposal, index_l, l_all)
      l_gp <- list_l$l_gp
      B_psi <- list_l$B_psi
      l_all[index_l] <- l_gp
      index_l <- index_l + 1
      
      # sample p ----------------
      
      if(sum(z) > 0){
        list_p <- update_p(beta_p, Occs, z_all, X_p, indexes_covariates_p, b_p, B_p)
        p <- list_p$p
        beta_p <- list_p$beta
      }
      
      # write parameters  -------------------------------------------------------
      
      if(iter > nburn){
        # z_output[chain,iter - nburn,,] <- z
        
        # psi_output[chain, iter - nburn,,] <- psi
        beta_psi_output[chain,iter - nburn,] <- beta_psi
        # gamma_psi_output[chain,iter - nburn,] <- gamma_psi[-c(1,2)]
        
        l_gp_psi_output[chain, iter - nburn] <- l_gp
        
        if(usingClustering){
          n_s_output[chain, iter - nburn,1:maxClusters] <- sort(table(c_s))[1:maxClusters]  
        }
        
        # a_s_output[chain, iter - nburn,] <- a_s
        
        # p_output[chain, iter - nburn,,,] <- p
        beta_p_output[chain,iter - nburn,] <- beta_p
        # gamma_p_output[chain,iter - nburn,] <- gamma_p[-c(1,2)]
      }
      
    }
    
  }
  
  if(usingClustering){
    
    # resize n_s
    {
      maxNumOfClusters <- max(apply(n_s_output, c(1,2), function(x){
        sum(!is.na(x))
      }))
      maxNumOfClusters <- 100
      n_s_output_new <- n_s_output[,,seq_len(maxNumOfClusters),drop = F]
    }
      
  }
  
  list("beta_psi_output" = beta_psi_output,
       "l_gp_psi_output" = l_gp_psi_output,
       "n_s_output" = n_s_output_new,
       # "a_s_output" = a_s_output,
       "beta_p_output" = beta_p_output)
}

runModel <- function(data, index_year, index_site, index_occ, 
                     num_covariates_psi_text, fac_covariates_psi_text,
                     num_covariates_p_text, fac_covariates_p_text,
                     prior_psi, sigma_psi, phi_psi,
                     sigma_gp, a_l_gp, b_l_gp,
                     prior_p, sigma_p, phi_p,
                     usingClustering, sigma_a, alpha_dp,
                     nchain, nburn, niter){
  
  print("Analyzing the data..")
  
  {
    listData <- cleanData(data, index_year, index_site, index_occ, 
                          num_covariates_psi_text, fac_covariates_psi_text,
                          num_covariates_p_text, fac_covariates_p_text)
    Occs <- listData$Occs
    k_s <- listData$k_s
    X_psi <- listData$X_psi
    X_y_index <- listData$X_y_index
    indexes_covariates_psi <- listData$indexes_covariates_psi
    X_p <- listData$X_p
    indexes_covariates_p <- listData$indexes_covariates_p
    numOfLevels_psi <- listData$numOfLevels_psi
    numOfLevels_p <- listData$numOfLevels_p
    ncov_psi <- listData$ncov_psi
    ncov_p <- listData$ncov_p
    Y <- listData$Y
    years <- listData$years
    S <- listData$S
    indexes_occobs <- listData$indexes_occobs
    index_numericalCovariates_psi <- listData$index_numericalCovariates_psi
    index_numericalCovariates_p <- listData$index_numericalCovariates_p
    nameVariables_psi <- listData$nameVariables_psi
    nameVariables_p <- listData$nameVariables_p  
  }
  
  dataCharacteristics <- list("Years" = years,
                              "nameVariables_psi" = nameVariables_psi,
                              "nameVariables_p" = nameVariables_p) 
  
  # assign the prior
  {
    listPrior <- createPrior(prior_psi, sigma_psi, phi_psi,
                             sigma_gp, a_l_gp, b_l_gp,
                             prior_p, sigma_p, phi_p,
                             Y, X_psi, ncov_psi, index_numericalCovariates_psi, numOfLevels_psi,
                             X_p, ncov_p, index_numericalCovariates_p, numOfLevels_p)
    b_psi <- listPrior$b_psi
    B_psi <- listPrior$B_psi
    b_p <- listPrior$b_p
    B_p <- listPrior$B_p  
  }
  
  # run MCMCM
  listMCMC <- runMCMC(Occs, X_psi, indexes_covariates_psi, ncov_psi, X_y_index,
                      S, Y, k_s, indexes_occobs,
                      X_p, indexes_covariates_p, ncov_p,
                      b_psi, B_psi, 
                      b_p, B_p,
                      sigma_gp, sd_l = .025, l_beta_proposal = .05, a_l_gp, b_l_gp,
                      usingClustering, sigma_a, alpha_dp, M_neal,
                      nchain, nburn, niter)
  beta_psi_output <- listMCMC$beta_psi_output
  l_gp_psi_output <- listMCMC$l_gp_psi_output
  beta_p_output <- listMCMC$beta_p_output
  n_s_output <- listMCMC$n_s_output
  # a_s_output <- listMCMC$a_s_output
  
  list("dataCharacteristics" = dataCharacteristics,
       "beta_psi_output" = beta_psi_output,
       "l_gp_psi_output" = l_gp_psi_output,
       "n_s_output" = n_s_output,
       # "a_s_output" = a_s_output,
       "beta_p_output" = beta_p_output)
}

# simulating data

simulateData <- function(Y, S, V,
                         mu_psi, beta_psi, b_t, 
                         mu_p, beta_p){
  
  # simulate occupancies
  {
    ncov_psi <- length(beta_psi)
    X_psi_cov <- matrix(rnorm(S * Y * ncov_psi, sd = 1), nrow = S * Y, ncol = ncov_psi)
    
    index_year <- rep(1:Y, each = S)
    X_psi_year <- data.frame(Year = factor(index_year))
    X_psi_year <- model.matrix(~ . - 1, X_psi_year)
    
    X_psi <- cbind(1, X_psi_year, X_psi_cov)
    
    beta_psi_all <- c(mu_psi, b_t, beta_psi)
    
    psi <- logit(X_psi %*% beta_psi_all)
    
    # occupancies
    z <- rep(NA, Y * S)
    for (i in 1:(Y * S)) {
      z[i] <- rbinom(1, 1, prob = psi[i])
    } 
    
    z_all <- rep(z, each = V)
  }
  
  # simulate detections
  {
    ncov_p <- length(beta_p)
    X_p <- matrix(rnorm(ncov_p * S * Y * V), nrow = S * Y * V, ncol = ncov_p)
    
    probsDetection <- logit(mu_p + X_p %*% beta_p)
    
    y_ysv <- rep(NA, nrow = S * Y * V)
    l <- 1
    for (y in 1:Y) {
      for (s in 1:S) {
        for (v in 1:V) {
          if(z_all[l] == 0){
            y_ysv[l] <- 0
          } else {
            y_ysv[l] <- rbinom(1, 1, probsDetection[l])
          }
          l <- l + 1
        }
      }
    }
    
  }
  
  # put data in data.frame
  {
    data <- data.frame(Year = rep(index_year, each = V),
                       Site = rep(rep(1:S, each = V), times = Y),
                       X_psi_cov[rep(1:nrow(X_psi_cov), each = V),],
                       X_p,
                       Occ = y_ysv)
    
  }
  
  return(data)
}

# plotting

plotYearseffect <- function(modelResults){
  
  years <- modelResults$dataCharacteristics$Years
  Y <- length(years)
  
  beta_psi_output <- modelResults$beta_psi_output
  
  beta_psi_output <- apply(beta_psi_output, 3, c)
  
  CI_yearseffect <- apply(beta_psi_output[,2 + 0:(Y-1)], 2, function(x){
    quantile(x, probs = c(0.025,0.5,0.975))
  })
  
  ggplot(data = NULL, aes(x = years,
                          y = CI_yearseffect[2,],
                          ymin = CI_yearseffect[1,],
                          ymax = CI_yearseffect[3,])) + geom_point() + geom_errorbar() + geom_line() + 
    ylab("Years Effect") + scale_x_continuous(name = "Year", breaks = years) +
    # geom_line(data = NULL, aes(x = years, y = b_t), color = "red") +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 13, face = "bold", angle = 90),
          panel.grid.major = element_line(colour="grey", size=0.015),
          panel.background = element_rect(fill = "white", color = "black"))
  
  
}

plotYearsOccProbs <- function(modelResults){
  
  years <- modelResults$dataCharacteristics$Years
  Y <- length(years)
  
  beta_psi_output <- modelResults$beta_psi_output
  
  beta_psi_output <- apply(beta_psi_output, 3, c)
  
  CI_yearseffect <- apply(beta_psi_output[,2 + 0:(Y-1)], 2, function(x){
    # quantile(x, probs = c(0.025,0.5,0.975))
    quantile(logit(beta_psi_output[,1] + x))
  })
  
  ggplot(data = NULL, aes(x = years,
                          y = CI_yearseffect[2,],
                          ymin = CI_yearseffect[1,],
                          ymax = CI_yearseffect[3,])) + geom_point() + geom_errorbar() + geom_line() + 
    ylab("Occupancy probability") + scale_x_continuous(name = "Year", breaks = years) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 13, face = "bold", angle = 90),
          panel.grid.major = element_line(colour="grey", size=0.015),
          panel.background = element_rect(fill = "white", color = "black"))
  
  
}

plotCovariatesPsiEffect <- function(modelResults){
  
  namesCovariates <- modelResults$dataCharacteristics$nameVariables_psi
  ncov_psi <- length(namesCovariates)
  
  years <- modelResults$dataCharacteristics$Years
  Y <- length(years)
  
  beta_psi_output <- modelResults$beta_psi_output
  
  beta_psi_output <- apply(beta_psi_output, 3, c)
  
  betaeffect <- apply(beta_psi_output[,(Y + 1) + 1:ncov_psi], 2, function(x){
    quantile(x, probs = c(0.025,0.5,0.975))
  })
  
  ggplot(data = NULL, aes(x = namesCovariates,
                          y = betaeffect[2,],
                          ymin = betaeffect[1,],
                          ymax = betaeffect[3,])) + geom_point() + geom_errorbar() +
    ylab("Effect") + scale_x_discrete(name = "Covariates") +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 16, face = "bold"),
          panel.grid.major = element_line(colour="grey", size=0.015),
          panel.background = element_rect(fill = "white", color = "black"))
  
  
}

plotBaseLineCaptureProb <- function(modelResults){
  
  beta_p_output <- modelResults$beta_p_output
  
  beta_p_output <- apply(beta_p_output, 3, c)
  
  peffect <- quantile(logit(beta_p_output[,1]), probs = c(0.025,0.5,0.975))

  ggplot(data = NULL, aes(x = "",
                          y = peffect[2],
                          ymin = peffect[1],
                          ymax = peffect[3])) + geom_point() + geom_errorbar() +
    ylab("Value") + scale_x_discrete(name = "Capture probability") +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 16, face = "bold"),
          panel.grid.major = element_line(colour="grey", size=0.015),
          panel.background = element_rect(fill = "white", color = "black"))
  
}

plotCovariatesPEffect <- function(modelResults){
  
  namesCovariates <- modelResults$dataCharacteristics$nameVariables_p
  ncov_p <- length(namesCovariates)
  
  beta_p_output <- modelResults$beta_p_output
  
  beta_p_output <- apply(beta_p_output, 3, c)
  
  betaeffect <- apply(beta_p_output[,-1,drop = F], 2, function(x){
    quantile(x, probs = c(0.025,0.5,0.975))
  })
  
  ggplot(data = NULL, aes(x = namesCovariates,
                          y = betaeffect[2,],
                          ymin = betaeffect[1,],
                          ymax = betaeffect[3,])) + geom_point() + geom_errorbar() +
    ylab("Effect") + scale_x_discrete(name = "Covariates") +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 16, face = "bold"),
          panel.grid.major = element_line(colour="grey", size=0.015),
          panel.background = element_rect(fill = "white", color = "black"))
  
  
}

plotNumberOfClusters <- function(modelResults, minElements){
  
  n_s_output <- modelResults$n_s_output
  
  n_s_output <- apply(n_s_output, 3, c)
  
  num_of_clusters <- apply(n_s_output, 1, function(x){
    sum(x >= minElements, na.rm = T)
  })
  
  ggplot(data = NULL, aes(x = num_of_clusters, y = ..density..)) + 
    geom_histogram(fill = "cornsilk", color = "black") +
    ylab("Density") + scale_x_continuous(breaks = min(num_of_clusters):max(num_of_clusters),
                                         name = "Number of Clusters") +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 16, face = "bold"),
          panel.grid.major = element_line(colour="grey", size=0.015),
          panel.background = element_rect(fill = "white", color = "black"))
}

# diagnostics

tracePlot_psi_intercept <- function(modelResults){
  
  beta_psi_output <- modelResults$beta_psi_output 
  
  nchain <- nrow(beta_psi_output)
  niter <- ncol(beta_psi_output)
  
  beta_psi_output <- matrix(beta_psi_output[,,1], nrow = nchain, ncol = niter)
  
  beta_psi_output_long <- melt(beta_psi_output)
  
  diagnosticsPlot <- ggplot() + geom_line(data = beta_psi_output_long, aes(x = Var2, y = value, group = Var1, 
                                                                           color = factor(Var1))) + 
    xlab("Iterations") + ylab("Value") +
    theme(plot.title = element_text(hjust = 0.5, size = 17),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 11, face = "bold", hjust = 1),
          axis.line = element_line(colour="black", size=0.15),
          # panel.grid.minor = element_line(colour="grey", size=0.15),
          panel.grid.major = element_line(colour="grey", size=0.15),
          panel.background = element_rect(fill = "white", color = "black"),
          legend.position = "none") 

  plotTitle <- createPlotTitle(beta_psi_output, nchain)
  
  diagnosticsPlot <- diagnosticsPlot + ggtitle(plotTitle)
  
  diagnosticsPlot
}

tracePlot_psi_beta <- function(modelResults, index_cov){
  
  beta_psi_output <- modelResults$beta_psi_output 
  
  Y <- length(modelResults$dataCharacteristics$Years)
  
  if((Y + 1 + index_cov) <= dim(beta_psi_output)[3]){
    
    nchain <- nrow(beta_psi_output)
    niter <- ncol(beta_psi_output)
    
    beta_psi_output <- matrix(beta_psi_output[,,Y + 1 + index_cov], nrow = nchain, ncol = niter)
    
    beta_psi_output_long <- melt(beta_psi_output)
    
    diagnosticsPlot <- ggplot(data = beta_psi_output_long, aes(x = Var2, y = value, group = Var1, 
                                                               color = factor(Var1))) + geom_line() + 
      xlab("Iterations") + ylab("Value") +
      theme(plot.title = element_text(hjust = 0.5, size = 17),
            axis.title = element_text(size = 16, face = "bold"),
            axis.text.y = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(size = 11, face = "bold", hjust = 1),
            axis.line = element_line(colour="black", size=0.15),
            # panel.grid.minor = element_line(colour="grey", size=0.15),
            panel.grid.major = element_line(colour="grey", size=0.15),
            panel.background = element_rect(fill = "white", color = "black"),
            legend.position = "none") 
    
    plotTitle <- createPlotTitle(beta_psi_output, nchain)
    
    diagnosticsPlot <- diagnosticsPlot + ggtitle(plotTitle)
    
    return(diagnosticsPlot)
    
  } else {
    
    print("Index above the number of covariates")
    
  }
  
  
}

tracePlot_bt <- function(modelResults, index_year){
  
  beta_psi_output <- modelResults$beta_psi_output
  
  nchain <- nrow(beta_psi_output)
  niter <- ncol(beta_psi_output)
  
  beta_psi_output <- matrix(beta_psi_output[,,index_year + 1], nrow = nchain, ncol = niter)
  
  beta_psi_output_long <- melt(beta_psi_output)
  
  diagnosticsPlot <- ggplot(data = beta_psi_output_long, aes(x = Var2, y = value, group = Var1, 
                                          color = factor(Var1))) + geom_line() + 
    xlab("Iterations") + ylab("Value") +
    theme(plot.title = element_text(hjust = 0.5, size = 17),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 11, face = "bold", hjust = 1),
          axis.line = element_line(colour="black", size=0.15),
          # panel.grid.minor = element_line(colour="grey", size=0.15),
          panel.grid.major = element_line(colour="grey", size=0.15),
          panel.background = element_rect(fill = "white", color = "black"),
          legend.position = "none") 
  
  
  plotTitle <- createPlotTitle(beta_psi_output, nchain)
  
  diagnosticsPlot <- diagnosticsPlot + ggtitle(plotTitle)
  
  diagnosticsPlot
}

tracePlot_p_intercept <- function(modelResults){
  
  beta_p_output <- modelResults$beta_p_output 
  
  nchain <- nrow(beta_p_output)
  niter <- ncol(beta_p_output)
  
  beta_p_output <- matrix(beta_p_output[,,1], nrow = nchain, ncol = niter)
  
  beta_psi_output_long <- melt(beta_p_output)
  
  diagnosticsPlot <- ggplot(data = beta_psi_output_long, aes(x = Var2, y = value, group = Var1, 
                                          color = factor(Var1))) + geom_line() + 
    xlab("Iterations") + ylab("Value") +
    theme(plot.title = element_text(hjust = 0.5, size = 17),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 11, face = "bold", hjust = 1),
          axis.line = element_line(colour="black", size=0.15),
          # panel.grid.minor = element_line(colour="grey", size=0.15),
          panel.grid.major = element_line(colour="grey", size=0.15),
          panel.background = element_rect(fill = "white", color = "black"),
          legend.position = "none") 

  plotTitle <- createPlotTitle(beta_p_output, nchain)
  
  diagnosticsPlot <- diagnosticsPlot + ggtitle(plotTitle)
  
  diagnosticsPlot
}

tracePlot_p_beta <- function(modelResults, index_cov){
  
  beta_p_output <- modelResults$beta_p_output 
  
  if((1 + index_cov) <= dim(beta_p_output)[3]){
    
    nchain <- nrow(beta_p_output)
    niter <- ncol(beta_p_output)
    
    beta_p_output <- matrix(beta_p_output[,,1 + index_cov], nrow = nchain, ncol = niter)
    
    beta_psi_output_long <- melt(beta_p_output)
    
    diagnosticsPlot <- ggplot(data = beta_psi_output_long, aes(x = Var2, y = value, group = Var1, 
                                                               color = factor(Var1))) + geom_line() + 
      xlab("Iterations") + ylab("Value") +
      theme(plot.title = element_text(hjust = 0.5, size = 17),
            axis.title = element_text(size = 16, face = "bold"),
            axis.text.y = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(size = 11, face = "bold", hjust = 1),
            axis.line = element_line(colour="black", size=0.15),
            # panel.grid.minor = element_line(colour="grey", size=0.15),
            panel.grid.major = element_line(colour="grey", size=0.15),
            panel.background = element_rect(fill = "white", color = "black"),
            legend.position = "none") 
    
    plotTitle <- createPlotTitle(beta_p_output, nchain)
    
    diagnosticsPlot <- diagnosticsPlot + ggtitle(plotTitle)
    
    return(diagnosticsPlot)
    
  } else {
    
    print("Index above the number of covariates")
    
  }
  
}

tracePlot_l <- function(modelResults){
  
  l_gp_psi_output <- modelResults$l_gp_psi_output
  
  nchain <- nrow(l_gp_psi_output)
  
  l_gp_psi_output_long <- melt(l_gp_psi_output)
  
  diagnosticsPlot <- ggplot(data = l_gp_psi_output_long, aes(x = Var2, y = value, group = Var1, 
                                          color = factor(Var1))) + geom_line() + 
    xlab("Iterations") + ylab("Value") +
    theme(plot.title = element_text(hjust = 0.5, size = 17),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 11, face = "bold", hjust = 1),
          axis.line = element_line(colour="black", size=0.15),
          # panel.grid.minor = element_line(colour="grey", size=0.15),
          panel.grid.major = element_line(colour="grey", size=0.15),
          panel.background = element_rect(fill = "white", color = "black"),
          legend.position = "none")
  
  plotTitle <- createPlotTitle(l_gp_psi_output, nchain)
  
  diagnosticsPlot <- diagnosticsPlot + ggtitle(plotTitle)
  
  diagnosticsPlot

}

tracePlot_numOfClusters <- function(modelResults){

  n_s_output <- modelResults$n_s_output
  numOfClusters_output <- apply(n_s_output, c(1,2), function(x){
    sum(!is.na(x))
  })
  
  nchain <- nrow(numOfClusters_output)
  
  numOfClusters_output_long <- melt(numOfClusters_output)
  
  diagnosticsPlot <- ggplot(data = numOfClusters_output_long, aes(x = Var2, y = value, group = Var1, 
                                          color = factor(Var1))) + geom_line() + 
    xlab("Iterations") + ylab("Value") +
    theme(plot.title = element_text(hjust = 0.5, size = 17),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 11, face = "bold", hjust = 1),
          axis.line = element_line(colour="black", size=0.15),
          # panel.grid.minor = element_line(colour="grey", size=0.15),
          panel.grid.major = element_line(colour="grey", size=0.15),
          panel.background = element_rect(fill = "white", color = "black"),
          legend.position = "none")
  
  plotTitle <- createPlotTitle(numOfClusters_output, nchain)
  
  diagnosticsPlot <- diagnosticsPlot + ggtitle(plotTitle)
  
  diagnosticsPlot

}

createPlotTitle <- function(mcmc_output, nchain){
  
  eff_samplesize <- ess(mcmc_output)
  
  plotTitle <- paste0("Effective sample size = ",round(eff_samplesize, 1))
  
  if(nchain > 1){
    Rhat <- compute_rhat(mcmc_output)
    
    plotTitle <- paste0(plotTitle,paste0(" / R.hat = ",round(Rhat, 3)))
  }
  
  
}

compute_rhat <- function(mcmc_output){
  
  mcmc_output_list <- lapply(1:nrow(mcmc_output), function(i){
    mcmc(mcmc_output[i,])
  })
  mcmc_output_list_2 <- as.mcmc.list(mcmc_output_list)

  Rhat <- gelman.diag(mcmc_output_list_2)  
  
  Rhat$psrf[2]
}

ess <- function(mcmc_output){
    
  mcmc_output_list <- lapply(1:nrow(mcmc_output), function(i){
    mcmc(mcmc_output[i,])
  })
  mcmc_output_list_2 <- as.mcmc.list(mcmc_output_list)
  
  effectiveSize(mcmc_output_list_2)
  
}

compute_rhat_alex <- function(mcmc_output){

  nchain <- nrow(mcmc_output)
  niter <- ncol(mcmc_output)

  chains_mean <- apply(mcmc_output, 1, mean)
  B_rhat <- niter * var(chains_mean)

  s_m_sq <- apply(mcmc_output, 1, var)
  W_rhat <- mean(s_m_sq)

  post_variance <- ((niter - 1) / niter) * W_rhat + ((B_rhat) / niter)

  sqrt(post_variance / W_rhat)
}



