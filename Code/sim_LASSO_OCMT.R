## MC simulation

# load required packages
library(tidyverse)
library(bayesSurv)
library(glmnet)
library(parallel)
library(MASS)
library(rootSolve)
library(bivariate)
library(pracma)
library(mvtnorm)
library(rmutil)
library(genlasso)
library(HDPenReg)





####################################################################
##                                                                ##
##   SigmaMat_generator: Covariance matrix for X                  ##
##                                                                ##
##        1. p: number of variables                               ##
##        2. base: speed of decay which default is 0.5            ##
##                                                                ##
##        RETURN: A covariance matrix                             ##
##                                                                ##
####################################################################

SigmaMat_generator <- function(p, base = 0.5) {
  
  # create matrix: pxp
  resMat <- matrix(nrow = p, ncol = p)
  
  # apply two for loops to generate the matrix
  for (i in 1:p) {
    for (j in 1:p) {
      resMat[i, j] = base^(abs(i-j))
    }
  }
  
  resMat[1, 1] <- 3
  
  return(resMat)
}


##########################################################################
##                                                                      ##
##   simulate_DGP: simulate DGP according to the design                 ##
##                                                                      ##
##        1. Design: which DGP u wanna use for simulation               ##
##        2. n: number of observations                                  ##
##        3. p: number of variables                                     ##
##        4. alpha: true value for parameter of interest                ##
##        5. Rsq_d: R-square in second equation                         ##
##        6. Rsq_y: R-square in first equation                          ##  
##                                                                      ##
##        RETURN: a list containing the elements of DGP                 ##
##                                                                      ##
##########################################################################

simulate_DGP <- function(Design, n, p, alpha, Rsq_d, Rsq_y) {
  
  ####################### Simulating the data ##########################
  
  # list being returned when the function is called
  res <- list()
  
  # epsilon: error term in 1st equation
  e <- rnorm(n, 0, 1) %>% matrix(nrow = n, ncol = 1)
  
  # v: error term in 2nd equation
  v <- rnorm(n, 0, 1) %>% matrix(nrow = n, ncol = 1)
  
  
  ##################### Local setting for different DESIGN  ##################### 
  
  # Design 1: homoscedastic
  # Design 2: heteroscedastic
  # Design 3: combination of deterministic and random coefficients
  
  if (Design == 1) {
    
    # beta
    beta_0 = ((1/(1:p))^2 %>% matrix(nrow = p, ncol = 1))
    
    # Sigma: covariance matrix in X
    SigmaMat <- SigmaMat_generator(p, base = 0.5)
    
    # X: design matrix
    meanMat <- matrix(0, nrow = p, ncol = 1) # zero-mean vector
    X <- rMVNorm(n, mean = meanMat, Sigma = SigmaMat) %>% matrix(nrow = n, ncol = p)
    
    
    ###################### decided constant for y and d ######################
    
    ## decide constant for d and constant for y (c_d & c_y)
    
    # c_d: closed form solution
    c_d <- sqrt(((1/(1-Rsq_d)) - 1) / (t(beta_0)%*%SigmaMat%*%beta_0)) %>% as.numeric()
    
    # c_y: solve second-degree polynomial equation
    a <- t(beta_0)%*%SigmaMat%*%beta_0
    b <- 2*alpha*c_d*t(beta_0)%*%SigmaMat%*%beta_0
    c <- -((alpha^2+1)/(1-Rsq_y^2)-alpha^2-1)
    
    func <- function(x) a*x^2 + b*x + c # define second order polynomial equation
    solver <- multiroot(func, start = 0) # using Newton method 
    c_y <- solver$root
    
    ############################################################################
    
    
    # regression coefficient: theta = c * beta 
    # constant: control the R_square for each reduced form
    theta_g = c_y*beta_0
    theta_m = c_d*beta_0

    
    # treatment variabe: d
    sigma_d = 1
    d <- X%*%theta_m + sigma_d*v  # simulate d 
    
    # outcome variable: y
    sigma_y = 1
    y <- alpha*d + X%*%theta_g + sigma_y*e  # simulate y 
    
  } else if (Design == 2) {
    
    # beta
    beta_0 = (1/(1:p))^2 %>% matrix(nrow = p, ncol = 1)
    
    # Sigma: covariance matrix in X
    SigmaMat <- SigmaMat_generator(p, base = 0.5)
    
    # X: design matrix
    meanMat <- matrix(0, nrow = p, ncol = 1) # zero-mean vector
    X <- rMVNorm(n, mean = meanMat, Sigma = SigmaMat) %>% matrix(nrow = n, ncol = p)

    
    ########################## decided constant for y and d #########################
    
    ## decide constant for d and constant for y (c_d & c_y)
    
    # c_d: closed form solution
    c_d <- sqrt(((1/(1-Rsq_d)) - 1) / (t(beta_0)%*%SigmaMat%*%beta_0)) %>% as.numeric()
    
    # c_y: solve second-degree polynomial equation
    a <- alpha^2*(1-Rsq_y)*t(beta_0)%*%SigmaMat%*%beta_0
    b <- c_d*(1-Rsq_y)*(t(beta_0)%*%beta_0)
    c <- (c_d^2*t(beta_0)%*%SigmaMat%*%beta_0+1.25)*(1-Rsq_y)-1.25
    
    func <- function(x) a*x^2 + b*x + c # define second order polynomial equation
    solver <- multiroot(func, start = 0) # using Newton method 
    c_y <- solver$root
  
    
    ##############################################################################
    
    # regression coefficient: theta = c * beta 
    theta_g = c_y*beta_0
    theta_m = c_d*beta_0
    
    # treatment variable: d
    sigma_d <- sqrt((1+X%*%beta_0)^2 / (mean((1+X%*%beta_0)^2)))
    
    # outcome variable: y
    sigma_y <- sqrt((1+alpha*d+X%*%beta_0)^2 / (mean((1+alpha*d+X%*%beta_0)^2)))
    y <- alpha*d + X%*%theta_g  + sigma_y*e 
    
  } else if (Design == 3) {
    
    # beta
    beta_0 = (1/(1:p))^2 %>% matrix(nrow = p, ncol = 1) 
    
    # Sigma
    SigmaMat <- SigmaMat_generator(p, base = 0.5)
    
    # X: design matrix
    meanMat <- matrix(0, nrow = p, ncol = 1) # zero-mean vector 
    X <- rMVNorm(n, mean = meanMat, Sigma = SigmaMat) %>% matrix(nrow = n, ncol = p)
    
    # Generating thetas for Design 3
    theta_g <- matrix(nrow = p)
    theta_m <- matrix(nrow = p)
    
    ###################### decided constant for y and d ######################
    
    # decide c_d & c_y
    c_d <- sqrt(((1/(1-Rsq_d)) - 1) / (t(beta_0)%*%SigmaMat%*%beta_0)) %>% as.numeric()
    
    a <- alpha^2*(1-Rsq_y)*t(beta_0)%*%SigmaMat%*%beta_0
    b <- c_d*(1-Rsq_y)*(t(beta_0)%*%beta_0)
    c <- (c_d^2*t(beta_0)%*%SigmaMat%*%beta_0+1.25)*(1-Rsq_y)-1.25

    func <- function(x) a*x^2 + b*x + c
    solver <- multiroot(func, start = 0) # using Newton method 
    c_y <- solver$root
    
    
    ###############################################################################
    
    for (j in 1:200) {
      if (j <= 5) {
        theta_g[j, 1] <- c_y * ((1/j)^2) 
        theta_m[j, 1] <- c_d * ((1/j)^2)
      } else {
        theta_g[j, 1] <- rnorm(1, mean = 0, sd = sqrt(1/p))
        theta_m[j, 1] <- rnorm(1, mean = 0, sd = sqrt(1/p))
      }
    }
    
    # treatment variable: d
    sigma_d <- sqrt((1+X%*%beta_0)^2 / (mean((1+X%*%beta_0)^2)))
    d <- X%*%theta_m + sigma_d*v
    
    # outcome variable: y
    sigma_y <- sqrt((1+alpha*d+X%*%beta_0)^2 / (mean((1+alpha*d+X%*%beta_0)^2)))
    y <- alpha*d + X%*%theta_g + sigma_y*e 
    
  } else {
    
    print("No such design :(")
    return(1)
    
  } 
  
  res[["X"]] <- X
  res[["d"]] <- d
  res[["y"]] <- y
  res[["c_y"]] <- c_y
  res[["c_d"]] <- c_d
  res[["theta_g"]] <- theta_g
  res[["theta_m"]] <- theta_m
  res[["cov_mat"]] <- SigmaMat 
  
  return(res)
  
}


##################################################################################
##                                                                              ##
##   postDouble_lasso: return estimation of alpha by post-double LASSO method   ##
##                                                                              ##
##        1. n: number of observations                                          ##
##        2. p: number of variables                                             ##
##        1. lst: list containing simulated data                                ##
##                                                                              ##
##        RETURN: estimation of alpha by postDouble method                      ##
##                                                                              ##
##################################################################################

postDouble_LASSO <- function(n, p, lst) {
  
  # storing the index of remaining variables after LASSO 
  idx_fit1 <- vector(mode = "integer")
  idx_fit2 <- vector(mode = "integer")
  
  X <- lst$X
  y <- lst$y
  d <- lst$d
  
  ##################### Estimating procedure #######################
  
  # lambda proposed by (2.14)
  lambda <- 2*c*sqrt(n)*qnorm(1-gamma/(2*p))
  
  # penalty loading proposed by APPENDIX A # problem occurs here
  loading1 <- iterated_loading(n, p, X, y, nu = 0, K = 20)
  
  # fit the LASSO using lambda and loading defined above # regress X on y
  lasso1 <- 
    glmnet(X, y, family = "gaussian", alpha = 1, lambda = lambda/(10*n), 
           penalty.factor = loading1, intercept = F)
  
  # which variables still remain in the model?
  beta_lasso1 <- as.vector(lasso1$beta)
  idx_fit1 <- which(beta_lasso1 != 0)
  
  # penalty loading proposed by APPENDIX A
  loading2 <- iterated_loading(n, p, X, d, nu = 0, K = 20)
  
  # fit the LASSO using lambda and loading defined above # regress X on d
  lasso2.cv <- 
    glmnet(X, d, family = "gaussian", alpha = 1, lambda = lambda/(10*n), 
              penalty.factor = loading2, intercept = F)
  
  # which variables still remain in the model?
  beta_lasso2 <- as.vector(lasso2$beta)
  idx_fit2 <- which(beta_lasso2 != 0)
  
  # union of idx_fit1 and idx_fit2
  union_idx <- sort(union(idx_fit1, idx_fit2))
  
  # s_hat: number of elements in the union set
  s_hat <- length(union_idx)
  
  # define new design matrix used in post-OLS stage
  X_new <- cbind(d, X[, union_idx])
  
  # estimate parameter of interest by OLS # using gernalized inverse
  para_hat <- ginv(t(X_new)%*%X_new)%*%t(X_new)%*%y
  alpha_hat <- para_hat[1, 1]
  
  return(alpha_hat)
}


#######################################################################
##                                                                   ##
##   post_LASSO: return estimation of alpha by post-LASSO method     ##
##                                                                   ##
##        1. lst: list containing simuated data                      ##
##                                                                   ##
##        RETURN: estimation of alpha by post LASSO method           ##
##                                                                   ##
#######################################################################

post_LASSO <- function(n, p, lst) {
  
  # storting the index of remaining variables after LASSO which excluding alpha in penalized term.
  idx_fit <- vector(mode = "integer")
  
  X <- lst$X
  y <- lst$y
  d <- lst$d
  

  # define new design matrix used in single-LASSO
  X_new <- cbind(d, X)
  
  ######################### Estimating procedure ##########################
  
  # lambda proposed by (2.14)
  lambda <- 2*c*sqrt(n)*qnorm(1-gamma/(2*p))
  
  # penalty loading proposed by APPENDIX A
  loading <- iterated_loading(n, p, X, y, 0, 20)
  loading <- c(0, loading)
  
  # refit LASSO using adquent lambda
  lasso.naive <- 
    glmnet(X_new, y, family = "gaussian", alpha = 1, lambda = lambda/(10*n), intercept = F, 
           penalty.factor = loading)
  
  # which variables still remain in the model after model selection (contains d for sure)
  beta_lasso <- as.vector(lasso.naive$beta)
  idx_fit <- which(beta_lasso != 0)
  
  # define new design matrix used in post-staged OLS
  X_new <- subset(X_new, select = idx_fit)
  
  # estimate parameter by OLS
  para_hat <- ginv(t(X_new)%*%X_new)%*%t(X_new)%*%y
  alpha_hat <- para_hat[1, 1]
  
  return(alpha_hat)
}

##############################################################################
##                                                                          ##
##   iterated_loading: choose panelty loading by itereated method           ##
##                                                                          ##
##        1. n: number of observations                                      ##
##        2. p: number of variables                                         ##
##        3. X: design matrix containing covariates                         ##
##        4. Y: vector of respose variable                                  ##
##        4. nu: tolerance level                                            ##
##        5. K: maximum number of iteration                                 ##
##                                                                          ##
##        RETURN: a p-dimensional vector                                    ##
##                                                                          ##
##############################################################################

iterated_loading <- function(n, p, X, y, nu, K) {
  
  # constant
  c <- 1.1
  
  # gamma (significant level)
  gamma <- 0.05
  
  # lambda proposed by (2.14)
  init_lambda <- 2*c*sqrt(n)*qnorm(1-gamma/(2*p))
  
  # least squares estimator 
  beta_bar <- ginv(t(X)%*%X)%*%t(X)%*%y
  
  
  ############## estimating the penalty loadings #################
  
  # penalty loading in zero stage 
  init_loading <- vector(length = p, mode = "numeric")
  
  residuals <- as.vector(y - X%*%beta_bar)
  
  for (j in 1:p) {
    x_ij <- as.vector(subset(X, select = j))
    init_loading[j] = sqrt(mean(x_ij^2*residuals^2))
  }

  
  # number of stage
  k <- 0
  
  # tolerance level
  nu <- nu
  
  # upper bound on the number of iterations
  K <- K
  
  # define penalty loading for k stage
  current_loading <- vector(length = p, mode = "numeric")
  
  # define penalty loading for k+1 stage
  next_loading <- vector(length = p, mode = "numeric")
  
  # initial setting in order to make the iteration start # not important 
  next_loading <- init_loading
  
  while(1) {
    
    # increment k in each iteration
    k = k+1
    print(paste("This is", k, "stage"))
    
    # reassign l_k_next_hat to l_k_hat from previous iteration
    current_loading <- next_loading
    
    ### (1) compute post-LASSO estimator based on loading current_loading
    # D_mat <- diag(current_loading)
    # lasso_obj <- genlasso(y, X = X, D = D_mat, svd = T, approx = T, minlam = lambda)
    # beta_lasso <- coef(lasso_obj, lambda = init_lambda)$beta
    
    lasso_k <- 
      glmnet(X, y, family = "gaussian", alpha = 1, lambda = init_lambda/(10*n),
             intercept = F, penalty.factor = current_loading)
    
    beta_lasso <- as.vector(lasso_k$beta)
    
    # index set
    idx_fit <- which(beta_lasso != 0)
    
    # estimate post-LASSO estimator 
    X_temp <- subset(X, select = idx_fit)
    beta_tilda <- ginv(t(X_temp)%*%X_temp)%*%t(X_temp)%*%y
    
    
    ### (2) calculate next loading
    s_hat <- length(idx_fit)
    print(paste("the variable remaining after LASSO in", k, "stage:", s_hat))
    
    residuals <- y - X_temp%*%beta_tilda
    
    for (j in 1:p) {
      x_ij <- subset(X, select = j)
      next_loading[j] =  sqrt(mean(x_ij^2*residuals^2)) * sqrt(n / (n - s_hat))
    }
    
    ### (3): stopping criteria
    criteria <- max(abs(current_loading - next_loading))
    
    if (criteria <= nu | k > K) {
      print(paste("complete at", k, "stage!"))
      print(paste("criteria:", criteria))
      print(next_loading)
      print(paste("S_hat: ", s_hat))
      return(next_loading)
    }
  }
}


# another implementation for iterated loading proposed by Sparse models and methods for ....

iterated_loading2 <- function(n, p, X, y, nu, K) {
  
  # constant
  c <- 1.1
  
  # gamma (significant level)
  gamma <- 0.1/log(p)
  
  # lambda proposed by (2.14)
  init_lambda <- 2*c*sqrt(n)*qnorm(1-gamma/(2*p))

  
  ############## estimating the penalty loadings #################
  
  # penalty loading in zero stage 
  init_loading <- vector(length = p, mode = "numeric")
  
  deviations <- as.vector(y - mean(y))

  for (j in 1:p) {
    x_ij <- as.vector(subset(X, select = j))
    init_loading[j] = sqrt(mean(x_ij^2*deviations^2))
  }
  
  
  # number of stage
  k <- 0
  
  # tolerance level
  nu <-  0
  
  # upper bound on the number of iterations
  K <- 15
  
  # define penalty loading for k stage
  current_loading <- vector(length = p, mode = "numeric")
  
  # define penalty loading for k+1 stage
  next_loading <- vector(length = p, mode = "numeric")
  
  # initial setting in order to make the iteration start # not important 
  next_loading <- init_loading
  
  while(1) {
    
    # increment k in each iteration
    k = k+1
    print(paste("This is", k, "stage"))
    
    # reassign l_k_next_hat to l_k_hat from previous iteration
    current_loading <- next_loading
    
    ### (1) compute post-LASSO estimator based on loading current_loading
    
    
    lasso_k <- 
      glmnet(X, y, family = "gaussian", alpha = 1, lambda = init_lambda/(10*n),
             intercept = F, penalty.factor = current_loading)
    
    beta_lasso <- as.vector(lasso_k$beta)
    
    # index set
    idx_fit <- which(beta_lasso != 0)
    
    # estimate post-LASSO estimator 
    X_temp <- subset(X, select = idx_fit)
    beta_tilda <- ginv(t(X_temp)%*%X_temp)%*%t(X_temp)%*%y
    
    
    ### (2) calculate next loading
    s_hat <- length(idx_fit)
    print(paste("the variable remaining after LASSO in", k, "stage:", s_hat))
    
    residuals <- as.vector(y - X_temp%*%beta_tilda)
    
    for (j in 1:p) {
      x_ij <- as.vector(subset(X, select = j))
      next_loading[j] = sqrt(mean(x_ij^2*residuals^2))
    }
    
    ### (3): stopping criteria
    criteria <- max(abs(current_loading - next_loading))
    
    if (criteria <= nu | k > K) {
      print(paste("complete at", k, "stage!"))
      print(paste("criteria:", criteria))
      print(next_loading)
      print(paste("S_hat: ", s_hat))
      return(next_loading)
    }
  }
}








################################################################################
##                                                                            ##
##   Calc_t: calculate t-ratios in first stage by the closed form from (13)   ##
##                                                                            ##
##        1. n: number of observations                                        ##
##        2. p: number of variables                                           ##
##        3. X: design matrix containing covariates                           ##
##        4. Y: vector of respose variable                                    ##
##                                                                            ##
##        RETURN: a p-dimensional vector                                      ##
##                                                                            ##
################################################################################

Calc_t <- function(n, p, X, y) {
  
  # define tau
  tau <- matrix(rep(1, n), nrow = n, ncol = 1)
  
  # define M_tau
  M_tau <- diag(1, nrow = n, ncol = n) - (tau%*%t(tau)) / n
  
  # define vector used for storing sigma hat
  sigma_hat <- matrix(nrow = p, ncol = 1)
  
  # consider n bivariate regressions of y on x_it for i = 1, ..., p
  for (i in 1:p) {
    
    # define a design matrix containing only one variable
    X_new <- matrix(X[,i], nrow = n, ncol = 1)
    
    # OLS estimator 
    OLS <- ginv(t(X_new)%*%X_new)%*%t(X_new)%*%y
    
    # fitted value 
    y_hat <- X_new%*%OLS
    
    # residuals 
    e <- y - y_hat
    
    # sigma_hat
    sigma_hat[i] <- (t(e)%*%e)/n
  }
  
  # define vector used for storing t-ratio
  t_ratio <- matrix(nrow = p, ncol = 1)
  
  # calcuate t-ratio 
  for (i in 1:p) {
    nume <- t(X[, i]) %*% M_tau %*% y
    deno <- sigma_hat[i, ] * sqrt(t(X[,i] %*% M_tau %*% X[,i]))
    t_ratio[i, ] <- nume / deno
  }
  
  return(t_ratio)
}

################################################################################
##                                                                            ##
##   Calc_t_subsequentStage: calculate t-ratios in subsequent stage by (16)   ##
##                                                                            ##
##        1. n: number of observations                                        ##
##        2. p: number of variables                                           ##
##        3. X: design matrix containing covariates                           ##
##        4. Y: vector of respose variable                                    ##
##                                                                            ##
##        RETURN: a p-dimensional vector                                      ##
##                                                                            ##
################################################################################


Calc_t_subsequentStage <- function(n, p, X, y, idx_prev) {
  
  # define X_{j-1}
  X_prev <- subset(X, select = idx_prev)
  
  # define the dimension of column space of design matrix X
  k <- dim(X_prev)[2]
  
  # define M_{j-1}
  M_prev <- diag(1, nrow = n, ncol = n) - X_prev%*%ginv(t(X_prev)%*%X_prev)%*%t(X_prev)
  
  # create a vector used for storing sigma hat in subsequent stage 
  sigma_hat_subsequent <- matrix(nrow = p, ncol = 1)
  
  # create an index set excluding the index of variable selected in previous stages
  idx_for_reg <- seq(1, p)[-idx_prev]
    
  for (i in idx_for_reg) {
    
    # define a design matrix containing the variables selected in previous stage and x_it in idx_for_reg, one at a time.
    X_new <- subset(X, select = c(idx_prev, i))
    
    # OLS estimator 
    OLS <- ginv(t(X_new)%*%X_new)%*%t(X_new)%*%y
    
    # fitted value 
    y_hat <- X_new%*%OLS
    
    # residuals 
    e <- y - y_hat
    
    # sigma_hat
    sigma_hat_subsequent[i] <- (t(e)%*%e)/n
  }
  
  ############### t-ratios #################
  
  # create a vector used for storing t-ratio
  t_ratio_subsequent <- matrix(nrow = p, ncol = 1)
  
  # calcuate t-ratio 
  for (i in idx_for_reg) {
    nume <- t(X[, i]) %*% M_prev %*% y
    deno <- sigma_hat_subsequent[i, ] * sqrt(t(X[, i] %*% M_prev %*% X[,i]))
    t_ratio_subsequent[i, ] <- as.vector(nume) / as.vector(deno)
  }
  
  return(t_ratio_subsequent)
}



#######################################################################
##                                                                   ##
##   Calc_criticalValue: calculate critical value                    ##
##                                                                   ##
##        1. n: number of observations                               ##
##        2. gamma: significance level                               ##
##        3. t-ratio: t-ratios w.r.t. the corresponding stage        ##
##        4. c: constant                                             ##
##        5. delta: critical value exponent                          ##
##                                                                   ##
##        RETURN: a p-dimensional vector                             ##
##                                                                   ##
#######################################################################


Calc_criticalValue <- function(n, gamma, t_ratio, c, delta) {
  temp <- qnorm(1 - (gamma/(2*c*n^delta)))
  return(temp)
}


#######################################################################
##                                                                   ##
##   OCMT: use OCMT as a model selection approach                    ##
##                                                                   ##
##        1. n: number of observations                               ##
##        3. p: number of variables                                  ##
##        4. X: desgin matrix                                        ##
##        5. y: response variable                                    ##
##                                                                   ##
##        RETURN: a p-dimensional vector                             ##
##                                                                   ##
#######################################################################

OCMT <- function(n, p, X, y, gamma, const, delta, delta_star) {
  
  # index set for first stage and subsequent stage
  idx_set_firstStage <- integer() # initial
  idx_set_cur <- integer() # k
  idx_set_prev <- integer() # k-1
  
  
  ################## first stage #################
  
  # calculate t-ratios for first reduced equation 
  t_ratios <- as.vector(Calc_t(n, p, X, y))
  
  # critical value for first stage
  cV_firstStage <- Calc_criticalValue(n, gamma, t_ratios, const, delta) # delta = 1 i.e. Bonferroni
  
  # first-stage OCMT selection indicator 
  OCMT_selector <- abs(t_ratios) > cV_firstStage
  
  # which variable are selected in this stage
  idx_set_firstStage <- which(OCMT_selector)
  idx_set_cur <- idx_set_firstStage
  idx_set_prev <- idx_set_cur
  
  
  ############ subsequent stage #############
  
  # number of stage
  k <- 2
  
  while(1) {
    
    # once entering the new iteration, relpace the number of current idx to previous idx 
    idx_set_prev <- c(idx_set_cur, idx_set_prev)
    idx_set_prev <- unique(idx_set_prev)
    
    # calculate t-ratios in current stage
    t_ratios_cur <- Calc_t_subsequentStage(n, p, X, y, idx_prev = idx_set_prev)
    
    # critical value for first stage
    cV_cur <- Calc_criticalValue(n, gamma, t_ratios, const, delta_star) # delta* = 2
    
    # subsequent-stage OCMT selection indicator 
    OCMT_selector <- abs(t_ratios_cur) > cV_cur
    
    # which variables are selected in the current stage
    idx_set_cur <- which(OCMT_selector)
    
    # define a criteria used for stopping the algorthim
    criteria <- sum(OCMT_selector, na.rm = T)
    
    # if the criteria is satistified then break the while loop.
    if (criteria == 0) {
      print(paste("OCMT stopped in", k, "stage."))
      break()
    }
    
    k = k + 1
  }
  
  return(idx_set_prev)
} 


##########################################################################
##                                                                      ##
##   postDouble_OCMT: replace OCMT with LASSO in postDouble literature  ##
##                                                                      ##
##        1. n: number of observations                                  ##
##        3. p: number of variables                                     ##
##        4. X: desgin matrix                                           ##
##        5. y: response variable                                       ##
##                                                                      ##
##        RETURN: a p-dimensional vector                                ##
##                                                                      ##
##########################################################################

postDouble_OCMT <- function(n, p, lst, gamma, const, delta, delta_star) {
  
  # storing the index of remaining variables after LASSO 
  idx_fit1 <- vector(mode = "integer")
  idx_fit2 <- vector(mode = "integer")
  
  X <- lst$X
  y <- lst$y
  d <- lst$d
  
  ##################### Estimating procedure #######################
  
  # apply OCMT to first reduced equation
  idx_fit1 <- OCMT(n, p, X, y, gamma, const, delta, delta_star)

  # apply OCMT to second reducred equation
  idx_fit2 <- OCMT(n, p, X, d, gamma, const, delta, delta_star)
  
  # union of idx_fit1 and idx_fit2
  union_idx <- sort(union(idx_fit1, idx_fit2))
  
  # s_hat: number of elements in the union set
  s_hat <- length(union_idx)
  
  # define new design matrix used in post-OLS stage
  X_new <- cbind(d, X[, union_idx])
  
  # estimate parameter of interest by OLS # using gernalized inverse
  para_hat <- ginv(t(X_new)%*%X_new)%*%t(X_new)%*%y
  alpha_hat <- para_hat[1, 1]
  
  return(alpha_hat)
}


##############################################################################
##                                                                          ##
##   main: main function to make estimation based on different method       ##
##                                                                          ##
##        1. seed: control randomness                                       ##
##        2. n: number of observations                                      ##
##        3. p: number of covariates                                        ##
##        4, alpha: true value of parameter of interest                     ##
##        5. Design: which DGP u wanna use for simulation                   ##
##        6. Rsq_d: desired R-square for equation (2)                       ##
##        7. Rsq_y: desired R-square for equation (1)                       ##
##                                                                          ##
##############################################################################

main <- function(seed, Design, n, p, alpha, Rsq_d, Rsq_y) {
  
  # return object
  est_result <- list()
  
  # control the randomness
  set.seed(seed)
  
  # simulate the data
  sim <- simulate_DGP(Design, n, p, alpha, Rsq_d, Rsq_y)
  
  # estimate the alpha from postDouble LASSO method
  postDouble_LASSO_est <- postDouble_LASSO(n, p, sim)
  est_result[["postDouble_LASSO_est"]] <- postDouble_LASSO_est
  
  # estimate the alpha from post LASSO method
  post_LASSO_est <- post_LASSO(n, p, sim)
  est_result[["post_LASSO_est"]] <- post_LASSO_est
  
  # estimate the alpha from postDouble OCMT method
  postDouble_OCMT_est <- postDouble_OCMT(n, p, sim, 0.05, 1, 1, 2)
  est_result[["post_OCMT_est"]] <- postDouble_OCMT_est
  
  return(est_result)
}


##################### parallel computing ######################


################  Design 1: homoscedastic  ####################

numSim <- 1000
numCores <- 3
seeds <- sample(1:numSim, size = numSim, replace = F)

double_design1_loading1 <- 
  mclapply(seeds, main, Design = 1, n = 100, p = 200, alpha = 0.5, Rsq_d = 0.8, Rsq_y = 0.8, mc.cores = numCores) 

res <- data.frame()

for (i in 1:numSim) {
  temp <- data.frame(double_design1[[i]])
  res <- rbind(res, temp)
}


res %>%
  gather(key = "approach", value = "alpha_hat") %>%
  ggplot() +
  geom_histogram(aes(x = alpha_hat, fill = approach, color = approach), alpha = 0.3, bins = 80, position = "identity") +
  geom_vline(xintercept = 0.5, color = "black", alpha = 0.6, linetype="dashed") +
  scale_x_continuous(limits=c(-0.5, 1.5))

apply(res-0.5, 2, function(x){sqrt(mean(x))})




