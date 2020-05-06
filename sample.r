library(MASS); library(dplyr); 
library(mvtnorm); library(coda)
library(reshape2); library(ggplot2)

library(FastOccupancy)

# simulate data -----------------------------------------------------------

Y <- 5 # years
S <- 1000 # sites
V <- 50 # number of visits

mu_psi <- -1
beta_psi <- c(-1,1)
b_t <- mvrnorm(1, mu = rep(0, Y), Sigma = diag(0.1, nrow = Y))

mu_p <- -1
beta_p <- c(1,0)

simulatedData <- simulateData(Y, S, V,
                  mu_psi, beta_psi, b_t, 
                  mu_p, beta_p)

# data parameters
{
  index_year <- 1
  index_site <- 2
  index_occ <- 7
  num_covariates_psi_text <- "3-4"
  fac_covariates_psi_text <- "0"
  num_covariates_p_text <- "5-6" 
  fac_covariates_p_text <- "0"  
}

# prior parameters
{
  prior_psi <- .2
  sigma_psi <- 2
  phi_psi <- 2
  
  sigma_gp <- .1
  a_l_gp <- 2
  b_l_gp <- 1
  sd_l <- .05
  l_beta_proposal <- .05
  
  sigma_a <- .25
  alpha_dp <- 1
  M_neal <- 10
  
  prior_p <- .5
  sigma_p <- 2
  phi_p <- 2
  
}

# mcmc parameter
{
  nchain <- 2
  nburn <- 500
  niter <- 500
}

modelResults <- runModel(simulatedData, # data
                         index_year, # index of the column with year
                         index_site, # index of the column with site
                         index_occ, # index of the column with captures
                         num_covariates_psi_text, # index of numerical covariates for psi
                         fac_covariates_psi_text, # index of categorical covariates for psi
                         num_covariates_p_text, # index of numerical covariates for p
                         fac_covariates_p_text, # index of categorical covariates for p
                         prior_psi, # prior mean of occupancy probability
                         sigma_psi, # just as in the notes
                         phi_psi, # just as in the notes
                         sigma_gp,  # sigma_b of the notes
                         a_l_gp, # a_l in the notes
                         b_l_gp, # b_l in the notes
                         prior_p, # prior mean of capture probability
                         usingClustering = T,
                         sigma_p, # just as in the notes
                         phi_p, # just as in the notes
                         sigma_a, # just as in the notes
                         alpha_dp, # just as in the notes
                         nchain, # number of chains
                         nburn, # burn-in iterations
                         niter) # number of iterations

plotYearseffect(modelResults)
plotYearsOccProbs(modelResults)
plotCovariatesPsiEffect(modelResults)
plotBaseLineCaptureProb(modelResults)
plotCovariatesPEffect(modelResults)
plotNumberOfClusters(modelResults, 1)

tracePlot_psi_intercept(modelResults)
tracePlot_psi_beta(modelResults, 1)
tracePlot_bt(modelResults, 1)
tracePlot_p_intercept(modelResults)
tracePlot_p_beta(modelResults, 1)
tracePlot_l(modelResults)
tracePlot_numOfClusters(modelResults)

# true data ----------

# input
{
  setwd("C:/Users/alexd/OneDrive/Desktop/byron files")
  data <- read.csv(file = "Ringlet_BNM_1970_2014_processed.csv", stringsAsFactors = F)
  
  index_year <- 9
  index_site <- 3
  index_occ <- 7
  
  num_covariates_p_text <- "8"
  fac_covariates_p_text <- "0"
  
  num_covariates_psi_text <- "4-5"
  fac_covariates_psi_text <- "0"
  
}

# prior parameters
{
  prior_psi <- .2
  sigma_psi <- 2
  phi_psi <- 2
  
  sigma_gp <- .1
  a_l_gp <- 1
  b_l_gp <- 1
  sd_l <- .05
  
  sigma_a <- .25
  alpha_dp <- 1
  M_neal <- 10
  
  prior_p <- .5
  sigma_p <- 2
  phi_p <- 2
  
}

# mcmc parameter
{
  nchain <- 2
  nburn <- 10000
  niter <- 10000
}

modelResults <- runModel(data, 
                         index_year, # index of the column with year
                         index_site, # index of the column with site
                         index_occ, # index of the column with captures
                         num_covariates_psi_text, # index of numerical covariates for psi
                         fac_covariates_psi_text, # index of categorical covariates for psi
                         num_covariates_p_text, # index of numerical covariates for p
                         fac_covariates_p_text, # index of categorical covariates for p
                         prior_psi, # prior mean of occupancy probability
                         sigma_psi, # just as in the notes
                         phi_psi, # just as in the notes
                         sigma_gp,  # sigma_b of the notes
                         a_l_gp, # a_l in the notes
                         b_l_gp, # b_l in the notes
                         prior_p, # prior mean of capture probability
                         usingClustering = F,
                         sigma_p, # just as in the notes
                         phi_p, # just as in the notes
                         sigma_a, # just as in the notes
                         alpha_dp, # just as in the notes
                         nchain, # number of chains
                         nburn, # burn-in iterations
                         niter) # number of iterations

plotYearseffect(modelResults)
plotYearsOccProbs(modelResults)
plotCovariatesPsiEffect(modelResults)
plotBaseLineCaptureProb(modelResults)
plotCovariatesPEffect(modelResults)
plotNumberOfClusters(modelResults, 5)

tracePlot_psi_intercept(modelResults)
tracePlot_psi_beta(modelResults, 2)
tracePlot_bt(modelResults, 3)
tracePlot_p_intercept(modelResults)
tracePlot_p_beta(modelResults, 1)
tracePlot_l(modelResults)
tracePlot_numOfClusters(modelResults)
