library(MASS); library(Rcpp); library(RcppArmadillo); library(dplyr);  library(mvtnorm)
library(reshape2); library(ggplot2)
# setwd("C:/Users/ad603/Dropbox/R Folder/PhD/Byron Occupancy")
setwd("C:/Users/alexd/Dropbox/R Folder/PhD/Byron Occupancy")
sourceCpp("codecpp.cpp")
source('functions.r', echo=F)

# simulated data --------------------------------------------------------------------

# fixed parameters
{
  Y <- 20 # years
  S <- 100 # sites
  V <- 10 # mean number of visits
  
  # number of covariates
  ncov_psi <- 0
  ncov_p <- 0
}

# model parameters
{
  prior_psi <- .5
  
  l_gp <- 1
  sigma_gp <- .1
  
  prior_p <- .2
}

# simulate occupancies
{
  mu_psi <- invLogit(prior_psi)
  
  b_t <- mvrnorm(1, rep(0, Y), K(1:Y, 1:Y, sigma_gp, l_gp))
  index_year <- rep(1:Y, each = S)
  b_t_true <- b_t
  
  X_psi_cov <- matrix(rnorm(S * Y * ncov_psi, sd = 1), nrow = S * Y, ncol = ncov_psi)
  beta_psi_cov <- rnorm(ncov_psi, 0, sd = 1)
  
  X_psi_year <- data.frame(Year = factor(index_year))
  X_psi_year <- model.matrix(~ . - 1, X_psi_year)
  
  X_psi <- cbind(1, X_psi_year, X_psi_cov)
  
  beta_psi <- c(mu_psi, b_t, beta_psi_cov)
  
  psi <- logit(X_psi %*% beta_psi)
  
  # occupancies
  z <- rep(NA, Y * S)
  for (i in 1:(Y * S)) {
    z[i] <- rbinom(1, 1, prob = psi[i])
  } 
  
  z_all <- rep(z, each = V)
}

# simulate detections
{
  X_p <- matrix(rnorm(ncov_p * S * Y * V), nrow = S * Y * V, ncol = ncov_p)

  mu_p <- invLogit(prior_p)
 
  beta_p <- rnorm(ncov_p, 0, sd = 1) 
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
  data <- data.frame(Year = 1950 + rep(index_year, each = V),
                     Site = rep(rep(1:S, each = V), times = Y),
                     X_psi_cov[rep(1:nrow(X_psi_cov), each = V),],
                     X_p,
                     Occ = y_ysv)
  
  data <- data[sample.int(nrow(data)),]
}

write.csv(data, file = "simulatedData.csv", row.names = F)

# real data ---------------------------------------------------------------

# input
{
  # data_psi <- read.csv(file = "data_psi_1995.csv", stringsAsFactors = F)
  data <- read.csv(file = "Ringlet_BNM_1970_2014_processed.csv", stringsAsFactors = F)
  # data <- read.csv(file = "simulatedData.csv", stringsAsFactors = F)
  
  # index_year <- 1
  # index_site <- 2
  # index_occ <- 3
  index_year <- 9
  index_site <- 3
  index_occ <- 7
  
  # num_covariates_p_text <- "0"
  # fac_covariates_p_text <- "0"
  num_covariates_p_text <- "8"
  fac_covariates_p_text <- "0"
  
  # num_covariates_psi_text <- "0"
  # fac_covariates_psi_text <- "0"
  num_covariates_psi_text <- "4-5"
  fac_covariates_psi_text <- "0"
  
  usingCov <- rep(T, 2)
  usingGamma <- rep(F, 2)
  
}

# cleaning data
{
  
  colnames(data)[c(index_year, index_site, index_occ)] <- c("Year","Site","Occ")

  data <- data %>% arrange(Site, Year)
  
  S <- length(unique(data$Site))
  Y <- length(unique(data[,index_year]))
  
  # define data for occupancies
  
  {
    
    k_s <- data %>% group_by(Year, Site) %>% summarise(Occupied = as.numeric(sum(Occ) > 0),
                                                       Visits = n()) %>% arrange(Site, Year) 
    k_s <- as.data.frame(k_s)
    
    originalsites <- k_s$Site
    k_s$Site <- as.numeric(as.factor(k_s$Site))  
    
    # match row in data_p to row in k_s
    indexes_occobs <- rep(1:nrow(k_s), k_s$Visits)
    
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
    
    if(usingCov[1]){
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
    }
    
    X_psi_year <- model.matrix( ~ . - 1, data = X_psi[,c("Year"),drop = F])
    
    if(usingCov[1]){
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
    
    if(usingCov[2]){
      ncov_p <- ncol(X_p)
      nameVariables_p <- colnames(data_p)[column_covariate_p]  
    } else {
      ncov_p <- 0
      nameVariables_p <- c()
    }
    
    indexes_covariates_p <- c()
    indexes_covariates_p[1] <- 1
    k <- 2
    if(usingCov[2]){
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
      
    }
    
    
    if(!usingCov[2]){
      X_p <- matrix(1, nrow = nrow(data_p), ncol = 1)
    } else {
      X_p <- model.matrix(~ ., X_p)  
    }
    
  }
  
}

# priors -----------

# input priors
{
  prior_psi <- .2
  sigma_psi <- 2
  phi_psi <- 2
  
  sigma_gp <- .1
  a_l_gp <- 2
  b_l_gp <- 1
  sd_l <- .05
  ggplot(data = NULL, aes(x = c(0,3))) + stat_function(fun = dgamma, args = list(shape = a_l_gp,
                                                                                 rate = b_l_gp))
  sigma_a <- .25
  alpha_dp <- 1
  M_neal <- 10
  
  prior_p <- .5
  sigma_p <- 2
  phi_p <- 2
  
}

l_gp <- a_l_gp / b_l_gp

d_bar_psi <- floor(ncov_psi / 3)
d_bar_p <- floor(ncov_p / 3)

mu_psi <- invLogit(prior_psi)
mu_p <- invLogit(prior_p)

if(usingCov[1]){
  
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
  
  b_psi <- c(mu_psi, rep(0, ncol(X_psi) - 1))
  
  B_psi <- matrix(0, nrow = ncol(X_psi), ncol = ncol(X_psi))
  B_psi[1, 1] <- sigma_psi
  B_psi[2:(ncol(X_psi)), 2:(ncol(X_psi))] <- C
  B_psi[(Y + 2):ncol(X_psi), (Y + 2):ncol(X_psi)] <- B_psi[(Y + 2):ncol(X_psi), (Y + 2):ncol(X_psi)] * phi_psi

} else {
  
  C <- matrix(0, nrow = ncol(X_psi) - 1, ncol = ncol(X_psi) - 1)
  l_gp <- a_l_gp / b_l_gp
  C[1:Y, 1:Y] <- K(1:Y,1:Y, sigma_gp, l_gp)
  
  b_psi <- c(mu_psi, rep(0, ncol(X_psi) - 1))
  B_psi <- matrix(0, nrow = ncol(X_psi), ncol = ncol(X_psi))
  B_psi[1, 1] <- sigma_psi
  B_psi[2:(ncol(X_psi)), 2:(ncol(X_psi))] <- phi_psi * C
    
}

if(usingCov[2]){
  
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
  
  b_p <- c(mu_p, rep(0, ncol(X_p) - 1))
  B_p <- matrix(0, nrow = ncol(X_p), ncol = ncol(X_p))
  B_p[1, 1] <- sigma_p
  B_p[2:(ncol(X_p)), 2:(ncol(X_p))] <- phi_p * C

} else {

  b_p <- mu_p
  B_p <- matrix(sigma_p, nrow = ncol(X_p), ncol = ncol(X_p))
    
}

# mcmc --------

nchain <- 1
nburn <- 1000
niter <- 1000

{
  # z_output <- array(0, dim = c(nchain, niter, S, Y))
  
  # psi_output <- array(NA,  dim = c(nchain, niter, S, Y))
  beta_psi_output <- array(NA, dim = c(nchain , niter, length(indexes_covariates_psi)),
                           dimnames = list(c(), c(), colnames(X_psi)))
  gamma_psi_output <- array(NA, dim = c(nchain , niter, ncov_psi),
                            dimnames = list(c(), c(), nameVariables_psi))
  
  l_gp_psi_output <- matrix(NA, nrow = nchain, ncol = niter)
  acceptances_l <- rep(NA, niter)
  
  # p_output <- array(NA,  dim = c(nchain, niter, Y, S, V))
  beta_p_output <- array(NA, dim = c(nchain , niter, length(indexes_covariates_p)),
                            dimnames = list(c(), c(), colnames(X_p)))
  gamma_p_output <- array(NA, dim = c(nchain , niter, ncov_p),
                            dimnames = list(c(), c(), nameVariables_p))
  
}

for (chain in 1:nchain) {

  # initialize parameters
  {
    
    # covariates for psi
    {
      beta_psi <- rep(0, length(indexes_covariates_psi))
      if(usingGamma[1]){
        gamma_psi <- c(1, 1, rbinom(ncov_psi, 1, d_bar_psi / ncov_psi))  
      } else {
        gamma_psi <- rep(1, ncov_psi + 2)
      }
      
      if(usingCov[1]){
        beta_psi[1] <- mu_psi
        index_present <- indexes_covariates_psi %in% which(gamma_psi == 1)
        X_psi_gamma <- X_psi[,index_present,drop = F]
        psi <- as.vector(logit(X_psi_gamma %*% beta_psi[index_present,drop = F]))
      } else {
        beta_psi[1] <- mu_psi
        psi <- as.vector(logit(X_psi %*% beta_psi))
      }  
      
      # a_s variables
      {
        c_s <- sample(1:5, S, replace = T)
        a_s <- rep(0, nrow(k_s))
      }
      
    }
    
    # covariates for p
    {
      beta_p <- rep(0, length(indexes_covariates_p))
      if(usingGamma[2]){
        gamma_p <- c(1, 1, rbinom(ncov_p, 1, d_bar_p / ncov_p))  
      } else {
        gamma_p <- rep(1, ncov_p + 2)
      }
      
      if(usingCov[2]){
        beta_p[1] <- mu_p
        index_present <- indexes_covariates_p %in% which(gamma_p == 1)
        X_p_gamma <- X_p[,index_present,drop = F]
        p <- as.vector(logit(X_p_gamma %*% beta_p[index_present,drop = F]))
      } else {
        beta_p[1] <- mu_p
        p <- as.vector(logit(X_p %*% beta_p))
      }  
    }
    
    z <- rep(NA, nrow(k_s))
    for (i in 1:nrow(k_s)) {
      if(k_s[i,1] == 1){
        z[i] <- 1
      } else {
        z[i] <- rbinom(1, 1, psi[i])
      }
    }
    # z but same length of data_p
    z_all <- z[indexes_occobs]
    
    l_gp <- a_l_gp / b_l_gp
    
  }

  for (iter in seq_len(nburn + niter)) {

    print(paste0("Iteration = ",iter,
                 " - l = ",l_gp,
                 " - Number of clusters = ",length(table(c_s))))
    
    # sample z ----
    
    z <- sample_z_cpp(psi, p, as.matrix(k_s[,3:4]))
    z_all <- z[indexes_occobs]
    
    # sample psi ----
    
    list_psi <- update_psi(beta_psi, gamma_psi,  X_psi, 
                           b_psi, B_psi, d_bar_psi, 
                           z, k_s, unique(k_s$Site), Y, ncov_psi, X_y_index, indexes_covariates_psi,
                           usingCov[1], usingGamma[1], 
                           a_s, c_s, sigma_a, alpha_dp, M_neal,
                           l_gp, sigma_gp, sd_l = .1, a_l_gp, b_l_gp,
                           T, T, T) 
    psi <- list_psi$psi
    beta_psi <- list_psi$beta
    if(!isTRUE(all.equal(l_gp,list_psi$l_gp))) acceptances_l[iter - nburn] <- 1
    l_gp <- list_psi$l_gp
    a_s <- list_psi$a_s
    c_s <- list_psi$c_s
    B_psi <- list_psi$B_psi
    if(usingGamma[1]){
      gamma_psi <- list_psi$gamma
    }

    # sample p ----------------
    
    if(sum(z) > 0){
      list_p <- update_p(beta_p, gamma_p, data_p$Occ, z_all, X_p, indexes_covariates_p, b_p,
                           B_p, usingCov[2], usingGamma[2], d_bar_p)
      p <- list_p$p
      beta_p <- list_p$beta
      if(usingGamma[2]){
        gamma_p <- list_p$gamma
      }
    }

    # write parameters  -------------------------------------------------------
  
    if(iter > nburn){
      # z_output[chain,iter - nburn,,] <- z
      
      # psi_output[chain, iter - nburn,,] <- psi
      beta_psi_output[chain,iter - nburn,] <- beta_psi
      # gamma_psi_output[chain,iter - nburn,] <- gamma_psi[-c(1,2)]
      
      l_gp_psi_output[chain, iter - nburn] <- l_gp
      
      # p_output[chain, iter - nburn,,,] <- p
      beta_p_output[chain,iter - nburn,] <- beta_p
      # gamma_p_output[chain,iter - nburn,] <- gamma_p[-c(1,2)]
    }

  }

}

mcmc_results <- list("z_output" = z_output,
                     "psi_output" = psi_output,
                     "beta_psi_output" = beta_psi_output,
                     "gamma_psi_output" = gamma_psi_output,
                     "p_output" = p_output,
                     "beta_p_output" = beta_p_output,
                     "gamma_p_output" = gamma_p_output,
                     "indexes_covariates_psi" = indexes_covariates_psi,
                     "indexes_covariates_p" = indexes_covariates_p)


# YEARS EFFECT ------------------------------------------------------------

CI_yearseffect <- apply(beta_psi_output[1,,2 + 0:(Y-1)], 2, function(x){
  quantile(x, probs = c(0.025,0.5,0.975))
  # mean(logit(beta_psi_output[1,,1] + x))
})

qplot(1:Y, CI_yearseffect)

ggplot(data = NULL, aes(x = sort(unique(k_s$Year)),
                        y = CI_yearseffect[2,],
                        ymin = CI_yearseffect[1,],
                        ymax = CI_yearseffect[3,])) + geom_point() + geom_errorbar() + geom_line() + 
# geom_point(data = NULL, aes(x = sort(unique(k_s$Year)), y = b_t_true), color = "red") +
    ylab("Years Effect") + scale_x_continuous(name = "Year", breaks = seq(1970, 2014, by = 5))

# PLOTS -----------------------------------------------------------------------

years_effect <- apply(beta_psi_output[1,,2 + 0:(Y-1)], 2, mean)
qplot(sort(unique(k_s$Year)), years_effect, geom = "line") + geom_point() +
  theme_bw() + ylab("Years Effect") + scale_x_continuous(name = "Year", breaks = seq(1970, 2014, by = 5))

ggplot(data = NULL, aes(x = beta_psi_output[1,,47], y = ..density..)) + 
  geom_histogram(fill = "cornsilk", color = "black", binwidth =  .003) +
  theme_bw() + 
  geom_vline(aes(xintercept = quantile(beta_psi_output[1,,47], probs = c(0.025,0.975))), color = "red", size = 1) +
  xlab("Effect") + ylab("Density")

ggplot(data = NULL, aes(x = beta_psi_output[1,,48], y = ..density..)) + 
  geom_histogram(fill = "cornsilk", color = "black", binwidth =  .004) +
  theme_bw() + 
  geom_vline(aes(xintercept = quantile(beta_psi_output[1,,48], probs = c(0.025,0.975))), color = "red", size = 1) +
  xlab("Effect") + ylab("Density")

ggplot(data = NULL, aes(x = logit(beta_psi_output[1,,1]), y = ..density..)) + 
  geom_histogram(fill = "cornsilk", color = "black", binwidth =  .0125) +
  theme_bw() + geom_vline(aes(xintercept = mean(logit(beta_psi_output[1,,1]))), color = "red", size = 1) +
  # geom_vline(aes(xintercept = quantile(beta_psi_output[1,,1], probs = c(0.025,0.975))), color = "red", size = 1) +
  xlab("Probability") + ylab("Density")

# PLOTS -------------------------------------------------------------------

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

{
  
  # beta_psi - gamma_psi
  if(usingCov[1]){
    
    {
      beta_psi_output2 <- matrix(NA, nrow = niter * nchain, ncol = dim(beta_psi_output)[3])
      for (chain in 1:nchain) {
        beta_psi_output2[(chain - 1)*niter + 1:niter,] <- beta_psi_output[chain,,]
      }
      gamma_psi_output2 <- matrix(NA, nrow = niter * nchain, ncol = dim(gamma_psi_output)[3])
      for (chain in 1:nchain) {
        gamma_psi_output2[(chain - 1)*niter + 1:niter,] <- gamma_psi_output[chain,,]
      }
      
      beta_psi_output2[,1] <- logit(beta_psi_output2[,1])
      colnames(beta_psi_output2) <- dimnames(beta_psi_output)[[3]]
      
      CICoefficients_psi  <- sapply(1:dim(beta_psi_output2)[2], function(i){
        c(quantile(beta_psi_output2[,i], probs = c(0.025,0.975)),
          mean(beta_psi_output2[,i]))
      })
      
      beta0_psi_plot <- ggplot(data = NULL, aes(x = beta_psi_output2[,1], y = ..density..)) + 
        geom_histogram(fill = "cornsilk", color = "black") + ylab("") + xlab("Probability") + 
        theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,5,0), size = 14, face = "bold"),
              panel.background = element_rect(fill = "white"), 
              panel.border = element_rect(fill = NA, colour = "grey20"),
              panel.grid.major = element_line(colour = "grey92"), 
              panel.grid.minor = element_line(colour = "grey92", size = 0.25), 
              strip.background = element_rect(fill = "grey85", colour = "grey20"), 
              legend.key = element_rect(fill = "white", colour = NA)) + 
        ggtitle("Baseline occupancy probability")
      
      beta_psi_years_plot <- ggplot() +
        geom_errorbar(data = NULL, aes(x = colnames(beta_psi_output2)[2:Y], ymax=CICoefficients_psi[1,2:Y], 
                                       ymin=CICoefficients_psi[2,2:Y]),
                      width=0.2, size=1, color="black") +
        geom_point(data = NULL, aes(x = colnames(beta_psi_output2)[2:Y], y=CICoefficients_psi[3,2:Y]), size=4, shape=21, fill="white") +
        theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,5,0), size = 14, face = "bold"),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(fill = NA, colour = "grey20"),
              panel.grid.major = element_line(colour = "grey92"),
              panel.grid.minor = element_line(colour = "grey92", size = 0.25),
              strip.background = element_rect(fill = "grey85", colour = "grey20"),
              legend.key = element_rect(fill = "white", colour = NA)) +
        xlab("")+ ylab("")+ ggtitle("Covariates coefficient")
      
      beta_psi_plot <- ggplot() +
        geom_errorbar(data = NULL, aes(x = colnames(beta_psi_output2)[-seq(1,Y)], ymax=CICoefficients_psi[1,-seq(1,Y)], 
                                       ymin=CICoefficients_psi[2,-seq(1,Y)]),
                      width=0.2, size=1, color="black") +
        geom_point(data = NULL, aes(x = colnames(beta_psi_output2)[-seq(1,Y)], y=CICoefficients_psi[3,-seq(1,Y)]), size=4, shape=21, fill="white") +
        theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,5,0), size = 14, face = "bold"),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(fill = NA, colour = "grey20"),
              panel.grid.major = element_line(colour = "grey92"),
              panel.grid.minor = element_line(colour = "grey92", size = 0.25),
              strip.background = element_rect(fill = "grey85", colour = "grey20"),
              legend.key = element_rect(fill = "white", colour = NA)) +
        xlab("")+ ylab("")+ ggtitle("Covariates coefficient")
      
      # data_plot <- data.frame(name = dimnames(gamma_psi_output)[[3]],
      #                         prob = apply(gamma_psi_output2, 2, mean))
      # 
      # gamma_psi_plot <- ggplot(data_plot, aes(x=prob, y=reorder(name, prob))) +
      #   geom_point(size=3) + # Use a larger dot
      #   theme_bw() +
      #   xlab("Variables") + ylab("Probability") + 
      #   theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,5,0), size = 14, face = "bold"),
      #         panel.grid.major.x = element_blank(),
      #         panel.grid.minor.x = element_blank(),
      #         panel.grid.major.y = element_line(colour="grey60", linetype="dashed")) + 
      #   ggtitle("Covariates Importance")
    }
    
  } else {
    
    beta_psi_output2 <- matrix(NA, nrow = niter * nchain, ncol = dim(beta_psi_output)[3])
    for (chain in 1:nchain) {
      beta_psi_output2[(chain - 1)*niter + 1:niter,] <- beta_psi_output[chain,,]
    }
   
    beta_psi_output2[,1] <- logit(beta_psi_output2[,1])
    colnames(beta_psi_output2) <- dimnames(beta_psi_output)[[3]]
    
    CICoefficients_psi  <- sapply(1:dim(beta_psi_output2)[2], function(i){
      c(quantile(beta_psi_output2[,i], probs = c(0.025,0.975)),
        mean(beta_psi_output2[,i]))
    })
    
    beta0_psi_plot <- ggplot(data = NULL, aes(x = beta_psi_output2[,1], y = ..density..)) + 
      geom_histogram(fill = "cornsilk", color = "black") + ylab("") + xlab("Probability") + 
      theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,5,0), size = 14, face = "bold"),
            panel.background = element_rect(fill = "white"), 
            panel.border = element_rect(fill = NA, colour = "grey20"),
            panel.grid.major = element_line(colour = "grey92"), 
            panel.grid.minor = element_line(colour = "grey92", size = 0.25), 
            strip.background = element_rect(fill = "grey85", colour = "grey20"), 
            legend.key = element_rect(fill = "white", colour = NA)) + 
      ggtitle("Baseline occupancy probability")
    
    beta_psi_years_plot <- ggplot() +
      geom_errorbar(data = NULL, aes(x = colnames(beta_psi_output2)[2:Y], ymax=CICoefficients_psi[1,2:Y], 
                                     ymin=CICoefficients_psi[2,2:Y]),
                    width=0.2, size=1, color="black") +
      geom_point(data = NULL, aes(x = colnames(beta_psi_output2)[2:Y], y=CICoefficients_psi[3,2:Y]), size=4, shape=21, fill="white") +
      theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,5,0), size = 14, face = "bold"),
            panel.background = element_rect(fill = "white"),
            panel.border = element_rect(fill = NA, colour = "grey20"),
            panel.grid.major = element_line(colour = "grey92"),
            panel.grid.minor = element_line(colour = "grey92", size = 0.25),
            strip.background = element_rect(fill = "grey85", colour = "grey20"),
            legend.key = element_rect(fill = "white", colour = NA)) +
      xlab("")+ ylab("")+ ggtitle("Covariates coefficient")
    
    beta_psi_plot <- NULL
  }
  
  # beta_p - gamma_p
  if(usingCov[2]){
    
    {
      beta_p_output2 <- matrix(NA, nrow = niter * nchain, ncol = dim(beta_p_output)[3])
      for (chain in 1:nchain) {
        beta_p_output2[(chain - 1)*niter + 1:niter,] <- beta_p_output[chain,,]
      }
      gamma_p_output2 <- matrix(NA, nrow = niter * nchain, ncol = dim(gamma_p_output)[3])
      for (chain in 1:nchain) {
        gamma_p_output2[(chain - 1)*niter + 1:niter,] <- gamma_p_output[chain,,]
      }
      
      beta_p_output2[,1] <- logit(beta_p_output2[,1])
      colnames(beta_p_output2) <- dimnames(beta_p_output)[[3]]
      
      CICoefficients  <- sapply(1:dim(beta_p_output2)[2], function(i){
        c(quantile(beta_p_output2[,i], probs = c(0.025,0.975)),
          mean(beta_p_output2[,i]))
      })
      
      beta0_p_plot <- ggplot(data = NULL, aes(x = beta_p_output2[,1], y = ..density..)) + 
        geom_histogram(fill = "cornsilk", color = "black") + ylab("") + xlab("Probability") + 
        theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,5,0), size = 14, face = "bold"),
              panel.background = element_rect(fill = "white"), 
              panel.border = element_rect(fill = NA, colour = "grey20"),
              panel.grid.major = element_line(colour = "grey92"), 
              panel.grid.minor = element_line(colour = "grey92", size = 0.25), 
              strip.background = element_rect(fill = "grey85", colour = "grey20"), 
              legend.key = element_rect(fill = "white", colour = NA)) + 
        ggtitle("Baseline detection probability")
      
      beta_p_plot <- ggplot() +
        geom_errorbar(data = NULL, aes(x = colnames(beta_p_output2)[-1], ymax=CICoefficients[1,-1], 
                                       ymin=CICoefficients[2,-1]),
                      width=0.2, size=1, color="black") +
        geom_point(data = NULL, aes(x = colnames(beta_p_output2)[-1], y=CICoefficients[3,-1]), size=4, shape=21, fill="white") +
        theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,5,0), size = 14, face = "bold"),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(fill = NA, colour = "grey20"),
              panel.grid.major = element_line(colour = "grey92"),
              panel.grid.minor = element_line(colour = "grey92", size = 0.25),
              strip.background = element_rect(fill = "grey85", colour = "grey20"),
              legend.key = element_rect(fill = "white", colour = NA)) +
        xlab("")+ ylab("")+ ggtitle("Covariates coefficient")
      
      # data_plot <- data.frame(name = dimnames(gamma_psi_output)[[3]],
      #                         prob = apply(gamma_psi_output2, 2, mean))
      # 
      # gamma_psi_plot <- ggplot(data_plot, aes(x=prob, y=reorder(name, prob))) +
      #   geom_point(size=3) + # Use a larger dot
      #   theme_bw() +
      #   xlab("Variables") + ylab("Probability") + 
      #   theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,5,0), size = 14, face = "bold"),
      #         panel.grid.major.x = element_blank(),
      #         panel.grid.minor.x = element_blank(),
      #         panel.grid.major.y = element_line(colour="grey60", linetype="dashed")) + 
      #   ggtitle("Covariates Importance")
    }
    
  } else {
   
    {
      beta_p_output2 <- matrix(NA, nrow = niter * nchain, ncol = dim(beta_p_output)[3])
      for (chain in 1:nchain) {
        beta_p_output2[(chain - 1)*niter + 1:niter,] <- beta_p_output[chain,,]
      }
      
      beta_p_output2[,1] <- logit(beta_p_output2[,1])
      colnames(beta_p_output2) <- dimnames(beta_p_output)[[3]]
      
      CICoefficients  <- sapply(1:dim(beta_p_output2)[2], function(i){
        c(quantile(beta_p_output2[,i], probs = c(0.025,0.975)),
          mean(beta_p_output2[,i]))
      })
      
      beta0_p_plot <- ggplot(data = NULL, aes(x = beta_p_output2[,1], y = ..density..)) + 
        geom_histogram(fill = "cornsilk", color = "black") + ylab("") + xlab("Probability") + 
        theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,5,0), size = 14, face = "bold"),
              panel.background = element_rect(fill = "white"), 
              panel.border = element_rect(fill = NA, colour = "grey20"),
              panel.grid.major = element_line(colour = "grey92"), 
              panel.grid.minor = element_line(colour = "grey92", size = 0.25), 
              strip.background = element_rect(fill = "grey85", colour = "grey20"), 
              legend.key = element_rect(fill = "white", colour = NA)) + 
        ggtitle("Baseline detection probability")
      
      beta_p_plot <- NULL
      # data_plot <- data.frame(name = dimnames(gamma_psi_output)[[3]],
      #                         prob = apply(gamma_psi_output2, 2, mean))
      # 
      # gamma_psi_plot <- ggplot(data_plot, aes(x=prob, y=reorder(name, prob))) +
      #   geom_point(size=3) + # Use a larger dot
      #   theme_bw() +
      #   xlab("Variables") + ylab("Probability") + 
      #   theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,5,0), size = 14, face = "bold"),
      #         panel.grid.major.x = element_blank(),
      #         panel.grid.minor.x = element_blank(),
      #         panel.grid.major.y = element_line(colour="grey60", linetype="dashed")) + 
      #   ggtitle("Covariates Importance")
    }
    
  }
  
 
}

# OUTPUTS ------------------------------------------------------------------

# psi
{
  list_results <- mcmc_results
  psi_output <- list_results$psi_output
  
  nchain <- dim(psi_output)[1]
  niter <- dim(psi_output)[2]
  S <- dim(psi_output)[3]
  
  psi_output2 <- matrix(NA, nrow = niter * nchain, ncol = S)
  for (chain in 1:nchain) {
    psi_output2[(chain - 1)*niter + 1:niter,] <- psi_output[chain,,]
  }
  
  psi_means <- apply(psi_output2, 2, mean)
  ci95_psi <- apply(psi_output2, 2, function(x){
    quantile(x, probs = .975)
  })
  ci25_psi <- apply(psi_output2, 2, function(x){
    quantile(x, probs = .025)
  })
  
  df <- data.frame("Site" = 1:S, 
                   "Occupancy probability" = psi_means)
  
  column_presence <- input$presence_column
  if(column_presence != 0){
    k_s <- as.vector(presenceInput())
  } else {
    k_s <- rep(0, S)
  }
  
  df$Occupancy.probability[as.logical(k_s)] <- 1
  df$`97.5Credible Interval` <- ci95_psi
  df$`2.5Credible Interval` <- ci25_psi
  df$`97.5Credible Interval`[as.logical(k_s)] <- 1
  df$`2.5Credible Interval`[as.logical(k_s)] <- 1
  
}

# beta psi
{
  list_results <- mcmc_results
  
  beta_psi_output <- list_results$beta_psi_output
  gamma_psi_output <- list_results$gamma_psi_output
  indexes_covariates <- list_results$indexes_covariates
  
  nchain <- dim(beta_psi_output)[1]
  niter <- dim(beta_psi_output)[2]
  ncov <- dim(gamma_psi_output)[3] 
  ncov_all <- dim(beta_psi_output)[3]
  
  beta_psi_output2 <- matrix(NA, nrow = niter * nchain, ncol = ncov_all)
  for (chain in 1:nchain) {
    beta_psi_output2[(chain - 1)*niter + 1:niter,] <- beta_psi_output[chain,,]
  }
  gamma_psi_output2 <- matrix(NA, nrow = niter * nchain, ncol = ncov)
  for (chain in 1:nchain) {
    gamma_psi_output2[(chain - 1)*niter + 1:niter,] <- gamma_psi_output[chain,,]
  }
  
  # beta_psi_output2[,1] <- logit(beta_psi_output2[,1])
  colnames(beta_psi_output2) <- dimnames(beta_psi_output)[[3]]

  beta_psi_mean <- sapply(1:ncol(beta_psi_output2), function(i){
    if(i == 1){
      mean(beta_psi_output2[,i])
    } else {
      mean(beta_psi_output2[gamma_psi_output2[,indexes_covariates[i]-1]!=0,i])  
    }
  })
  ci95_psi <- sapply(1:ncol(beta_psi_output2), function(i){
    if(i == 1){
      quantile(beta_psi_output2[,i], probs = .975)  
    } else {
      quantile(beta_psi_output2[gamma_psi_output2[,indexes_covariates[i]-1]!=0,i], probs = .975)  
    }
  })
  ci25_psi <- sapply(1:ncol(beta_psi_output2), function(i){
    if(i == 1){
      quantile(beta_psi_output2[,i], probs = .025)
    } else {
      quantile(beta_psi_output2[gamma_psi_output2[,indexes_covariates[i] - 1]!=0,i], probs = .025)  
    }
  })
  
  df <- data.frame("Means" = beta_psi_mean)
  
  df$`97.5Credible Interval` <- ci95_psi
  df$`2.5Credible Interval` <- ci25_psi
  
  row.names(df) <- colnames(beta_psi_output2)
}

# theta11
{
  list_results <- mcmc_results
  theta11_output <- list_results$theta11_output
  
  nchain <- dim(theta11_output)[1]
  niter <- dim(theta11_output)[2]
  S <- dim(theta11_output)[3]
  
  theta11_output2 <- matrix(NA, nrow = niter * nchain, ncol = S)
  for (chain in 1:nchain) {
    theta11_output2[(chain - 1)*niter + 1:niter,] <- theta11_output[chain,,]
  }
  
  theta11_means <- apply(theta11_output2, 2, mean)
  ci95_theta11 <- apply(theta11_output2, 2, function(x){
    quantile(x, probs = .975)
  })
  ci25_theta11 <- apply(theta11_output2, 2, function(x){
    quantile(x, probs = .025)
  })
  
  df <- data.frame("Site" = 1:S, 
                   "Probability" = theta11_means)
  
  
  df$`97.5Credible Interval` <- ci95_theta11
  df$`2.5Credible Interval` <- ci25_theta11
  
  
}

# beta theta11
{
  list_results <- mcmc_results
  
  beta_theta11_output <- list_results$beta_theta11_output
  gamma_theta11_output <- list_results$gamma_theta11_output
  indexes_covariates <- list_results$indexes_covariates
  
  nchain <- dim(beta_theta11_output)[1]
  niter <- dim(beta_theta11_output)[2]
  ncov <- dim(gamma_theta11_output)[3] 
  ncov_all <- dim(beta_theta11_output)[3]
  
  beta_theta11_output2 <- matrix(NA, nrow = niter * nchain, ncol = ncov_all)
  for (chain in 1:nchain) {
    beta_theta11_output2[(chain - 1)*niter + 1:niter,] <- beta_theta11_output[chain,,]
  }
  gamma_theta11_output2 <- matrix(NA, nrow = niter * nchain, ncol = ncov)
  for (chain in 1:nchain) {
    gamma_theta11_output2[(chain - 1)*niter + 1:niter,] <- gamma_theta11_output[chain,,]
  }
  
  # beta_theta11_output2[,1] <- logit(beta_theta11_output2[,1])
  colnames(beta_theta11_output2) <- dimnames(beta_theta11_output)[[3]]
  
  beta_theta11_mean <- sapply(1:ncol(beta_theta11_output2), function(i){
    if(i == 1){
      mean(beta_theta11_output2[,i])
    } else {
      mean(beta_theta11_output2[gamma_theta11_output2[,indexes_covariates[i]-1]!=0,i])  
    }
  })
  ci95_theta11 <- sapply(1:ncol(beta_theta11_output2), function(i){
    if(i == 1){
      quantile(beta_theta11_output2[,i], probs = .975)  
    } else {
      quantile(beta_theta11_output2[gamma_theta11_output2[,indexes_covariates[i]-1]!=0,i], probs = .975)  
    }
  })
  ci25_theta11 <- sapply(1:ncol(beta_theta11_output2), function(i){
    if(i == 1){
      quantile(beta_theta11_output2[,i], probs = .025)
    } else {
      quantile(beta_theta11_output2[gamma_theta11_output2[,indexes_covariates[i] - 1]!=0,i], probs = .025)  
    }
  })
  
  df <- data.frame("Means" = beta_theta11_mean)
  
  df$`97.5Credible Interval` <- ci95_theta11
  df$`2.5Credible Interval` <- ci25_theta11
  
  row.names(df) <- colnames(beta_theta11_output2)
}

# theta10
{
  list_results <- mcmc_results
  theta10_output <- list_results$theta10_output
  
  nchain <- dim(theta10_output)[1]
  niter <- dim(theta10_output)[2]
  S <- dim(theta10_output)[3]
  
  theta10_output2 <- matrix(NA, nrow = niter * nchain, ncol = S)
  for (chain in 1:nchain) {
    theta10_output2[(chain - 1)*niter + 1:niter,] <- theta10_output[chain,,]
  }
  
  theta10_means <- apply(theta10_output2, 2, mean)
  ci95_theta10 <- apply(theta10_output2, 2, function(x){
    quantile(x, probs = .975)
  })
  ci25_theta10 <- apply(theta10_output2, 2, function(x){
    quantile(x, probs = .025)
  })
  
  df <- data.frame("Site" = 1:S, 
                   "Probability" = theta10_means)
  
  
  df$`97.5Credible Interval` <- ci95_theta10
  df$`2.5Credible Interval` <- ci25_theta10
  
  
}

# beta theta10
{
  list_results <- mcmc_results
  
  beta_theta10_output <- list_results$beta_theta10_output
  gamma_theta10_output <- list_results$gamma_theta10_output
  indexes_covariates <- list_results$indexes_covariates
  
  nchain <- dim(beta_theta10_output)[1]
  niter <- dim(beta_theta10_output)[2]
  ncov <- dim(gamma_theta10_output)[3] 
  ncov_all <- dim(beta_theta10_output)[3]
  
  beta_theta10_output2 <- matrix(NA, nrow = niter * nchain, ncol = ncov_all)
  for (chain in 1:nchain) {
    beta_theta10_output2[(chain - 1)*niter + 1:niter,] <- beta_theta10_output[chain,,]
  }
  gamma_theta10_output2 <- matrix(NA, nrow = niter * nchain, ncol = ncov)
  for (chain in 1:nchain) {
    gamma_theta10_output2[(chain - 1)*niter + 1:niter,] <- gamma_theta10_output[chain,,]
  }
  
  # beta_theta10_output2[,1] <- logit(beta_theta10_output2[,1])
  colnames(beta_theta10_output2) <- dimnames(beta_theta10_output)[[3]]
  
  beta_theta10_mean <- sapply(1:ncol(beta_theta10_output2), function(i){
    if(i == 1){
      mean(beta_theta10_output2[,i])
    } else {
      mean(beta_theta10_output2[gamma_theta10_output2[,indexes_covariates[i]-1]!=0,i])  
    }
  })
  ci95_theta10 <- sapply(1:ncol(beta_theta10_output2), function(i){
    if(i == 1){
      quantile(beta_theta10_output2[,i], probs = .975)  
    } else {
      quantile(beta_theta10_output2[gamma_theta10_output2[,indexes_covariates[i]-1]!=0,i], probs = .975)  
    }
  })
  ci25_theta10 <- sapply(1:ncol(beta_theta10_output2), function(i){
    if(i == 1){
      quantile(beta_theta10_output2[,i], probs = .025)
    } else {
      quantile(beta_theta10_output2[gamma_theta10_output2[,indexes_covariates[i] - 1]!=0,i], probs = .025)  
    }
  })
  
  df <- data.frame("Means" = beta_theta10_mean)
  
  df$`97.5Credible Interval` <- ci95_theta10
  df$`2.5Credible Interval` <- ci25_theta10
  
  row.names(df) <- colnames(beta_theta10_output2)
}

# p11
{
  list_results <- mcmc_results
  p11_output <- list_results$p11_output
  
  nchain <- dim(p11_output)[1]
  niter <- dim(p11_output)[2]
  S <- dim(p11_output)[3]
  
  p11_output2 <- matrix(NA, nrow = niter * nchain, ncol = S)
  for (chain in 1:nchain) {
    p11_output2[(chain - 1)*niter + 1:niter,] <- p11_output[chain,,]
  }
  
  p11_means <- apply(p11_output2, 2, mean)
  ci95_p11 <- apply(p11_output2, 2, function(x){
    quantile(x, probs = .975)
  })
  ci25_p11 <- apply(p11_output2, 2, function(x){
    quantile(x, probs = .025)
  })
  
  df <- data.frame("Site" = 1:S, 
                   "Probability" = p11_means)
  
  
  df$`97.5Credible Interval` <- ci95_p11
  df$`2.5Credible Interval` <- ci25_p11
  
  
}

# beta p11
{
  list_results <- mcmc_results
  
  beta_p11_output <- list_results$beta_p11_output
  gamma_p11_output <- list_results$gamma_p11_output
  indexes_covariates <- list_results$indexes_covariates
  
  nchain <- dim(beta_p11_output)[1]
  niter <- dim(beta_p11_output)[2]
  ncov <- dim(gamma_p11_output)[3] 
  ncov_all <- dim(beta_p11_output)[3]
  
  beta_p11_output2 <- matrix(NA, nrow = niter * nchain, ncol = ncov_all)
  for (chain in 1:nchain) {
    beta_p11_output2[(chain - 1)*niter + 1:niter,] <- beta_p11_output[chain,,]
  }
  gamma_p11_output2 <- matrix(NA, nrow = niter * nchain, ncol = ncov)
  for (chain in 1:nchain) {
    gamma_p11_output2[(chain - 1)*niter + 1:niter,] <- gamma_p11_output[chain,,]
  }
  
  # beta_p11_output2[,1] <- logit(beta_p11_output2[,1])
  colnames(beta_p11_output2) <- dimnames(beta_p11_output)[[3]]
  
  beta_p11_mean <- sapply(1:ncol(beta_p11_output2), function(i){
    if(i == 1){
      mean(beta_p11_output2[,i])
    } else {
      mean(beta_p11_output2[gamma_p11_output2[,indexes_covariates[i]-1]!=0,i])  
    }
  })
  ci95_p11 <- sapply(1:ncol(beta_p11_output2), function(i){
    if(i == 1){
      quantile(beta_p11_output2[,i], probs = .975)  
    } else {
      quantile(beta_p11_output2[gamma_p11_output2[,indexes_covariates[i]-1]!=0,i], probs = .975)  
    }
  })
  ci25_p11 <- sapply(1:ncol(beta_p11_output2), function(i){
    if(i == 1){
      quantile(beta_p11_output2[,i], probs = .025)
    } else {
      quantile(beta_p11_output2[gamma_p11_output2[,indexes_covariates[i] - 1]!=0,i], probs = .025)  
    }
  })
  
  df <- data.frame("Means" = beta_p11_mean)
  
  df$`97.5Credible Interval` <- ci95_p11
  df$`2.5Credible Interval` <- ci25_p11
  
  row.names(df) <- colnames(beta_p11_output2)
}

# p10
{
  list_results <- mcmc_results
  p10_output <- list_results$p10_output
  
  nchain <- dim(p10_output)[1]
  niter <- dim(p10_output)[2]
  S <- dim(p10_output)[3]
  
  p10_output2 <- matrix(NA, nrow = niter * nchain, ncol = S)
  for (chain in 1:nchain) {
    p10_output2[(chain - 1)*niter + 1:niter,] <- p10_output[chain,,]
  }
  
  p10_means <- apply(p10_output2, 2, mean)
  ci95_p10 <- apply(p10_output2, 2, function(x){
    quantile(x, probs = .975)
  })
  ci25_p10 <- apply(p10_output2, 2, function(x){
    quantile(x, probs = .025)
  })
  
  df <- data.frame("Site" = 1:S, 
                   "Probability" = p10_means)
  
  
  df$`97.5Credible Interval` <- ci95_p10
  df$`2.5Credible Interval` <- ci25_p10
  
}

# beta p10
{
  list_results <- mcmc_results
  
  beta_p10_output <- list_results$beta_p10_output
  gamma_p10_output <- list_results$gamma_p10_output
  indexes_covariates <- list_results$indexes_covariates
  
  nchain <- dim(beta_p10_output)[1]
  niter <- dim(beta_p10_output)[2]
  ncov <- dim(gamma_p10_output)[3] 
  ncov_all <- dim(beta_p10_output)[3]
  
  beta_p10_output2 <- matrix(NA, nrow = niter * nchain, ncol = ncov_all)
  for (chain in 1:nchain) {
    beta_p10_output2[(chain - 1)*niter + 1:niter,] <- beta_p10_output[chain,,]
  }
  gamma_p10_output2 <- matrix(NA, nrow = niter * nchain, ncol = ncov)
  for (chain in 1:nchain) {
    gamma_p10_output2[(chain - 1)*niter + 1:niter,] <- gamma_p10_output[chain,,]
  }
  
  # beta_p10_output2[,1] <- logit(beta_p10_output2[,1])
  colnames(beta_p10_output2) <- dimnames(beta_p10_output)[[3]]
  
  beta_p10_mean <- sapply(1:ncol(beta_p10_output2), function(i){
    if(i == 1){
      mean(beta_p10_output2[,i])
    } else {
      mean(beta_p10_output2[gamma_p10_output2[,indexes_covariates[i]-1]!=0,i])  
    }
  })
  ci95_p10 <- sapply(1:ncol(beta_p10_output2), function(i){
    if(i == 1){
      quantile(beta_p10_output2[,i], probs = .975)  
    } else {
      quantile(beta_p10_output2[gamma_p10_output2[,indexes_covariates[i]-1]!=0,i], probs = .975)  
    }
  })
  ci25_p10 <- sapply(1:ncol(beta_p10_output2), function(i){
    if(i == 1){
      quantile(beta_p10_output2[,i], probs = .025)
    } else {
      quantile(beta_p10_output2[gamma_p10_output2[,indexes_covariates[i] - 1]!=0,i], probs = .025)  
    }
  })
  
  df <- data.frame("Means" = beta_p10_mean)
  
  df$`97.5Credible Interval` <- ci95_p10
  df$`2.5Credible Interval` <- ci25_p10
  
  row.names(df) <- colnames(beta_p10_output2)
}

# plots - Sites-specific occupancy probabilities ----------------------------------------

library(ggplot2)

psi_output2 <- matrix(NA, nrow = niter * nchain, ncol = S)
for (chain in 1:nchain) {
  psi_output2[(chain - 1)*niter + 1:niter,] <- psi_output[chain,,]
}

psi_output_long <- melt(psi_output2)

if(usingCov[1]){
  ggplot(data = psi_output_long, aes(x = as.factor(Var2), y = value)) + geom_boxplot(fill = "cornsilk") +
    # theme(panel.background = element_rect(fill = "white")) +
    theme_bw() + scale_y_continuous(breaks = seq(0, 1, .1), name = "Probability") +
    xlab("Sites")
} else {
  ggplot(data = NULL, aes(x = "1", y = psi_output2[,1])) + geom_boxplot() +
    # theme(panel.background = element_rect(fill = "white")) +
    theme_bw() + scale_y_continuous(breaks = seq(0, 1, .1), name = "Probability") +
    xlab("")
}

psi_means <- apply(psi_output2, 2, mean)

df <- data.frame("Site" = 1:S, "Occupancy probability" = psi_means)

df$Occupancy.probability[as.logical(k_s)] <- 1

write.csv(df, file = "output.csv", row.names = F)

# plots - Baseline occupancy probabilities and covariates effect ----------------------------------------

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  # plots <- c(list(...), plotlist)

  # alex addition
  plots <- list()
  for(i in 1:length(plotlist)){
    for(j in 1:length(plotlist[[1]])){
      plots[[j + (i - 1) * length(plotlist[[1]])]] <- plotlist[[i]][[j]]
    }
  }
  #

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

library(ggplot2)

beta_psi_output2 <- matrix(NA, nrow = niter * nchain, ncol = dim(beta_psi_output)[3])
for (chain in 1:nchain) {
  beta_psi_output2[(chain - 1)*niter + 1:niter,] <- beta_psi_output[chain,,]
}
gamma_psi_output2 <- matrix(NA, nrow = niter * nchain, ncol = dim(gamma_psi_output)[3])
for (chain in 1:nchain) {
  gamma_psi_output2[(chain - 1)*niter + 1:niter,] <- gamma_psi_output[chain,,]
}

beta_psi_output2[,1] <- logit(beta_psi_output2[,1])
colnames(beta_psi_output2) <- dimnames(beta_psi_output)[[3]]

# beta_psi_output_long <- melt(beta_psi_output2)

ggplot(data = NULL, aes(x = beta_psi_output2[,1], y = ..density..)) +
  geom_histogram(fill = "cornsilk", color = "black") + ylab("") + xlab("Probability") +
  theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,5,0), size = 14, face = "bold"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "grey20"),
        panel.grid.major = element_line(colour = "grey92"),
        panel.grid.minor = element_line(colour = "grey92", size = 0.25),
        strip.background = element_rect(fill = "grey85", colour = "grey20"),
        legend.key = element_rect(fill = "white", colour = NA), complete = TRUE) +
  ggtitle("Baseline occupancy probability")

# mean_coefs <- apply(beta_psi_output2[,-1,drop = FALSE], 2, mean)
# upper95_coefs <- apply(beta_psi_output2[,-1,drop = FALSE], 2, function(x) {
#   quantile(x, probs = c(.975))
# })
# lower95_coefs <- apply(beta_psi_output2[,-1,drop = FALSE], 2, function(x) {
#   quantile(x, probs = c(.025))
# })

CICoefficients  <- sapply(1:dim(beta_psi_output2)[2], function(i){
  if(i == 1){
    c(quantile(beta_psi_output2[,1], probs = c(0.025,0.975)),
      mean(beta_psi_output2[,1]))
  } else {
    c(quantile(beta_psi_output2[gamma_psi_output2[,indexes_covariates[i]-1]!= 0,i], probs = c(0.025,0.975)),
      mean(beta_psi_output2[gamma_psi_output2[,indexes_covariates[i]-1]!= 0,i]))
  }
})

data_plot <- data.frame(name = dimnames(gamma_psi_output)[[3]],
                        prob = apply(gamma_psi_output2, 2, mean))

ggplot(data_plot, aes(x=prob, y=reorder(name, prob))) +
  geom_point(size=3) + # Use a larger dot
  theme_bw() +
  xlab("Variables") + ylab("Probability") + 
  theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,5,0), size = 14, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour="grey60", linetype="dashed")) + 
  ggtitle("Covariates Importance")

ggplot() +
  geom_errorbar(data = NULL, aes(x = colnames(beta_psi_output2)[-1], ymax=CICoefficients[1,-1], 
                                 ymin=CICoefficients[2,-1]),
                width=0.2, size=1, color="black") +
  geom_point(data = NULL, aes(x = colnames(beta_psi_output2)[-1], y=CICoefficients[3,-1]), size=4, shape=21, fill="white") +
  theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,5,0), size = 14, face = "bold"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "grey20"),
        panel.grid.major = element_line(colour = "grey92"),
        panel.grid.minor = element_line(colour = "grey92", size = 0.25),
        strip.background = element_rect(fill = "grey85", colour = "grey20"),
        legend.key = element_rect(fill = "white", colour = NA), complete = TRUE) +
  xlab("")+ ylab("")+ ggtitle("Covariates coefficient")


beta_psi_mean <- apply(beta_psi_output2, 2, mean)

beta_psi_mean[1] <- logit(beta_psi_mean[1])

write.csv(as.data.frame(t(beta_psi_mean)), file = "output.csv", row.names = F)

# plots - Sites-specific eDNA presence probabilities (given species presence) -------------------

library(ggplot2)

library(reshape2)

theta11_output2 <- matrix(NA, nrow = niter * nchain, ncol = S)
for (chain in 1:nchain) {
  theta11_output2[(chain - 1)*niter + 1:niter,] <- theta11_output[chain,,]
}

theta11_output_long <- melt(theta11_output2)

if(usingCov[2]){
  ggplot(data = theta11_output_long, aes(x = as.factor(Var2), y = value)) +
    geom_boxplot(fill = "cornsilk") +
    # theme(panel.background = element_rect(fill = "white")) +
    theme_bw() + scale_y_continuous(breaks = seq(0, 1, .1), name = "Probability") +
    xlab("Sites")
} else {
  ggplot(data = NULL, aes(x = "Site", y = theta11_output2[,1])) + geom_boxplot() +
    # theme(panel.background = element_rect(fill = "white")) +
    theme_bw() + scale_y_continuous(breaks = seq(0, 1, .1), name = "Probability") +
    xlab("")
}

theta11_means <- apply(theta11_output2, 2, mean)

df <- data.frame("Site" = 1:S, "Probability" = theta11_means)

write.csv(df, file = "output.csv", row.names = F)

# plots - Sites-specific eDNA presence probabilities (given species absence) -------------------

library(ggplot2)

library(reshape2)

theta10_output2 <- matrix(NA, nrow = niter * nchain, ncol = S)
for (chain in 1:nchain) {
  theta10_output2[(chain - 1)*niter + 1:niter,] <- theta10_output[chain,,]
}

theta10_output_long <- melt(theta10_output2)

if(usingCov[2]){
  ggplot(data = theta10_output_long, aes(x = as.factor(Var2), y = value)) +
    geom_boxplot(fill = "cornsilk") +
    # theme(panel.background = element_rect(fill = "white")) +
    theme_bw() + scale_y_continuous(breaks = seq(0, 1, .1), name = "Probability") +
    xlab("Sites")
} else {
  ggplot(data = NULL, aes(x = "Site", y = theta10_output2[,1])) + geom_boxplot() +
    # theme(panel.background = element_rect(fill = "white")) +
    theme_bw() + scale_y_continuous(breaks = seq(0, 1, .1), name = "Probability") +
    xlab("")
}

theta10_means <- apply(theta10_output2, 2, mean)

df <- data.frame("Site" = 1:S, "Probability" = theta10_means)

write.csv(df, file = "output.csv", row.names = F)
