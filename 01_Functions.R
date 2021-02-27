##################
#Functions for 
#
#Author: Anthony Weishampel
#Date: 2/27/2021
#
##################





library(splines)
library(fields)
library(class)
library(gamm4)
library(mvtnorm)
library(MLmetrics)
library(parallel)
library(naivebayes)
library(xtable)
library(fda)
library(refund)
library(Matrix)
library(MASS)
library(fields)
library(arm)
library(kernlab)
library(svd)

inv.mat <- function(mat){
  
  x2 = svd2(mat)
  x2$d=ifelse(x2$d>0, x2$d, 0)
  mat=round(x2$u%*%diag(x2$d)%*%t(x2$v), 7)
  
  chol_mat = chol(mat)
  chol_mat_inv = forwardsolve(t(chol_mat), diag(dim(mat)[1]))
  chol_mat_inv = backsolve(chol_mat, chol_mat_inv)
  
  return(chol_mat_inv)
  
}


####
#Function to generate data in the various scenarios
#
#Inputs:
# Scenario: Which scenario to generate data for: 1=A, 2=B, 3=C In paper
# grid: which set of points within [0,1] will data be observed on
# N: Total number of subjects to generate data for
# p: vector showing distirbution of the two groups 
# sigma: variance in the latent curves
# binary: Reports binary or latent curves
#
# Outputs: Binary or latent curves
#####
generate_data <- function(scenario = 1 , grid=seq(from = 0,to = 1, length.out = 48),
                          N, p = rep(0.5, 2), sigma = 0.1, binary = T){
  
  #For Scenario A in the paper
  if(scenario == 1){
    
    #number of outputs
    D = length(grid)
    
    #Number of subjects per groups
    N1 = round(p[1]*N)
    N2 = N-N1
    
    #set variance struture or scneario A
    theta1<- theta0 <-1/(c(1:40)^2)
    #set mean functions
    mu0<-  c(0, -0.5, 1, -0.5,1,-0.5)
    mu1<- c(0, -0.75, 0.75, -0.15,1.4,0.1)
    
    #generate latent functions for both groups given the mean and variance
    X0 <- generate_latent_process(theta=theta0, mu_coeff =mu0, tt=seq(0,1, len=D), n=N1)
    X1 <- generate_latent_process(theta=theta1, mu_coeff = mu1, tt=seq(0,1, len=D), n=N2)
    
    #if binary curves wanted
    if(binary){
      #binary_values
      Y0 <- generate_binary_fns(X0) 
      Y1 <- generate_binary_fns(X1) 
      Curves_binary = t(cbind(Y0, Y1))
    }else{
      #non_binary
      epsilon = matrix(rnorm(N*D, 0, sigma), ncol=D, nrow=N)
      Curves_binary = t(cbind(X0, X1))+epsilon
    }
    
    return(Curves_binary)
    
  }
  #For Scenario B
  if(scenario == 2){
    
    # length of response
    D = length(grid)
    
    # number of curves for the two groups
    N1 = round(p[1]*N)
    N2 = N-N1
    
    #mean functions
    beta0_true = rep(0, D)
    beta1_true = rep(0, D)

    #variance for each group
    k=51
    Ks= 1:k
    lambda_k0 = exp(-Ks/3)
    lambda_k1 = exp(-Ks/2)
    
    #matrix to hold function values
    psis = matrix(rep(1, D), ncol = 1)
    for(i in 1:((k-1)/2)){
      p1 = (i)*2
      psi_1 = sqrt(2)*cos(p1*pi*grid)
      psi_2 = sqrt(2)*sin(p1*pi*grid)
      psis = cbind(psis , psi_1, psi_2)
    }
    
    for(i in 1:k){
      #set variance
      sigma0 = lambda_k0[i]
      sigma1 = lambda_k1[i]
      #get scores
      c_11_g1=rnorm(N1, 0, sigma0)
      c_11_g2=rnorm(N2, 0, sigma1)
      c_11=c(c_11_g1, c_11_g2)
      
      #set matrix values
      if(i==1){
        c1 = c_11
      }else{
        c1 = cbind(c1, c_11)
      }
      
    }
    
    #set mean functions matrixes
    beta0_true_mat = kronecker(matrix(1, nrow=N1), matrix(beta0_true, ncol=D))
    beta1_true_mat = kronecker(matrix(1, nrow=N2), matrix(beta1_true, ncol=D))
    beta0_true_mat = rbind (beta0_true_mat, beta1_true_mat)
    
    #effects of functions
    c1_mat = c1%*%t(psis)
    
    # estimate variance
    epsilon = matrix(rnorm(N*D, 0, sigma), ncol=D, nrow=N*J)
    
    Curves = beta0_true_mat + c1_mat + epsilon
    
    if(binary){
      Curves_binary = t(matrix(rbinom(D*N,1,  p=c(t(invlogit(Curves)))), nrow = D))
    }else{
      Curves_binary = t(matrix(c(t(Curves)), nrow = D))
    }
    
    return(Curves_binary)
    
  }
  #Scenario C
  if(scenario == 3){
    
    # get length of responses
    D = length(grid)
    
    # number of curves for the two groups
    N1 = round(p[1]*N)
    N2 = N-N1
    
    n1 = N1
    n2 = N2 
    
    #ndf is global variable for running scenario C
    is=1:ndf
    # estimate the scores for each component
    scores1 = sapply(is, function(x) get.scores(x, n1, n2))
    
    #mean function globally defined
    mu_hat = mu_t_hat
    eigen_funcs1 = eigen_funcs_true
    #generate latent curves 
    latent_curves = t(matrix(rep(mu_hat, n1+n2), ncol = n1+n2)) +
      t(eigen_funcs1%*%t(scores1))
    
    #get_indexes of where the lantent curves need to be evaluated
    grid2 = seq(from = 0,to = 1, length.out = dim(latent_curves)[2])
    grid_index = sapply(grid, function(x) which.min(abs(grid2-x)))
    Curves = latent_curves[,grid_index]
    
    if(binary){
      Curves_binary = t(matrix(rbinom(D*N,1,  p=c(t(invlogit(Curves)))), nrow = D))
    }else{
      Curves_binary = t(matrix(c(t(Curves)), nrow = D))
    }
    
    return(Curves_binary)
    
  }
  
}


####
#Function that makes a matrix positive semi definite
#
#Inputs: A square-symmetric matrix
#Outputs: A positive matrix
####
make_pos_semi_def = function(x){
  x2 = svd2(x)
  #remove negative vals
  x2$d=ifelse(x2$d>0, x2$d, 0)
  #remake matrix
  x2=x2$u%*%diag(x2$d)%*%t(x2$v)
  return(x2)
}

####
#Function: Get number of functions based on the pvs and eigenvalues
#
#Inputs: 
# PVE: Value [0,1] to determine number of eigenfunctions
# vec: vector of eigenvalues
# set_max_number: If you want a maximum number of values 
#
#Outputs: 
# The number of eigenfunctions and eigenvalues for KL-approximation
#
####
get_length_pve = function(pve, vec, set_max_number = NA){
  
  #total sum of eigenvals
  s=sum(vec)
  vec2 = cumsum(vec)
  vec=vec2/s
  #return(which(vec>=pve)[1])
  if(!is.na(set_max_number)){
    if(which(vec>=pve)[1] < set_max_number){
      return(which(vec>=pve)[1])
    }else{
      return(set_max_number)
    }
  }
  return(which(vec>=pve)[1])
}


#####
#Fucnction to estimate eiginfunctions
#
#Inputs: 
# K_b: estimate for covariance matrix
# pve: Proportion of variance explained
# fix_num_of_functions: set the number of eigenfuctions to be returned
#
#Outputs: 
# 
#
#####
estimate_eigenfunctions = function(K_b, pve=0.98, fix_num_of_functions=0){
  
  #SVD
  svd_kb = svd2(K_b)
  
  #Only get the pve length if the number of functions is not given
  if(fix_num_of_functions==0){
    pb = get_length_pve(pve, svd_kb$d)
  }else{
    pb = fix_num_of_functions
  }
  pb = ifelse(pb>1, pb, 2)
  eigen_vals1 = svd_kb$d[1:pb]
  eigen_funcs1 = svd_kb$v[,1:pb]
  
  return(list(eigen_vals1=eigen_vals1, eigen_funcs1=eigen_funcs1))
  
}

fpca_estimate_functions <-  function(Curves, pve=0.98, fix_number_of_functions = 0){
  
  K_s = var(Curves)
  
  D = dim(K_s)[1]
  Ys = cbind(as.vector(K_s), rep(1:D, D), as.vector(t(kronecker(matrix(1, ncol=D) , 1:D) )))
  #remove diagonals to remove excess error
  Ys = Ys[which(Ys[,2]!=Ys[,3]), ]
  K_s2 = smooth.2d(Ys[,1], ind = Ys[, 2:3] , nrow = D, ncol= D, theta = 4)
  K_s = K_s2$z
  K_s = make_pos_semi_def(K_s)
  
  fpca_results = estimate_eigenfunctions(K_s, pve, fix_number_of_functions)
  fpca_results$eigen_funcs1 = fpca_results$eigen_funcs1*sqrt(D)
  fpca_results$eigen_vals1 = fpca_results$eigen_vals1/D
  sigma = mean(diag(var(Curves)-K_s))
  
  # fit = smooth.spline(as.vector(colMeans(Curves)))
  # mu_hat = fit$y
  
  z1 = as.vector(colMeans(Curves))
  gam1 <- gam(z1~s(tt, bs = "cr", m=2, k = 10), 
              family="gaussian", method = "REML")
  mu_hat = gam1$fitted.values
  
  
  return(list(eigen_vals=fpca_results$eigen_vals1,
              eigen_funcs=fpca_results$eigen_funcs1,
              sigma = sigma, mu_hat = mu_hat))
  
}


#Function to return the logit
logit <- function(x){
  return(log(x/(1-x)))
}

#Function to return the inverse of the logit
invlogit <- function(x){
  return(1/(1+exp(-x)))
}

#function to return the derivative of the logit
d_logit <- function(x){
  return( logit(x)*(1-logit(x)))
}


#Method for extracting scores (training an testing sets) based on Hall function 2008
GFPCA_estimate_para_with_test <- function(Curves_binary, Curves_binary_test){
  
  D = dim(Curves_binary)[2]
  N = dim(Curves_binary)[1]
  N_test = dim(Curves_binary_test)[1]
  tt = seq(0,1, len=D)
  
  #estimated smooth mean function
  z1 = as.vector(colMeans(Curves_binary))
  gam1 <- gam(z1~s(tt, bs = "cr", m=2, k = 10),
              family="gaussian", method = "REML")
  #because probability needs to ensure values are between 0 and 1
  alpha_t = ifelse(gam1$fitted.values<0.001, 0.001, gam1$fitted.values)
  alpha_t = ifelse(alpha_t>0.999, 0.999, alpha_t)
  
  # fit = smooth.spline(as.vector(colMeans(Curves_binary)))#, lambda = 0.01)
  # #make sure probabilities aren't 0 or 1 after smoothing
  # alpha_t = ifelse(fit$y<0.001, 0.001, fit$y) 
  # alpha_t = ifelse(alpha_t>0.999, 0.999, alpha_t)
  
  v_t_hat = as.vector(invlogit(alpha_t))
  alpha_t_hat = as.vector(alpha_t)
  
  beta_ts  = t(Curves_binary)%*%Curves_binary/N
  
  Bs = cbind(as.vector(beta_ts), rep(1:D, D), as.vector(t(kronecker(matrix(1, ncol=D) , 1:D) )))
  #remove diagonals to remove excess error
  Bs = Bs[which(Bs[,2]!=Bs[,3]), ]
  #smooth the covariance matrix
  beta_ts2 = smooth.2d(Bs[,1], ind = Bs[, 2:3] , nrow = D, ncol= D, theta = 4)
  beta_ts_hat = beta_ts2$z
  tau_ts = (beta_ts_hat - matrix(alpha_t_hat, ncol = 1) %*% t(alpha_t_hat))/
    (matrix(d_logit(v_t_hat), ncol = 1)%*%t(d_logit(v_t_hat)))
  
  #make sym and semi-post def
  tau_ts_hat = round(make_pos_semi_def(tau_ts), 5)
  mu_t_hat = v_t_hat
  
  #decomposition of resulting covariance matrix
  tau_eigens  =  estimate_eigenfunctions(tau_ts_hat, pve = 0.98)
  
  #estimate from Hall 2008
  d = t((Curves_binary - matrix(rep(logit(mu_t_hat), each = N), ncol = D)) /
          matrix(rep(d_logit(mu_t_hat), each = N), ncol = D))
  
  #it's small enough to use solve
  #sigma_inv = solve(tau_ts_hat+diag(0.01, D))
  sigma_inv = solve(tau_ts_hat)
  
  d_test = t((Curves_binary_test - matrix(rep(logit(mu_t_hat), each = N_test), ncol = D)) /
               matrix(rep( d_logit(mu_t_hat), each = N_test), ncol = D))
  
  #estimate scores via approximation presented in Hall 2008
  score_coef = t(t(matrix(rep(tau_eigens$eigen_vals1, each = N), nrow = N))*
                   (t(tau_eigens$eigen_funcs1)%*%sigma_inv%*%d))
  
  #estimate scores in testing set via approximation presented in Hall 2008
  score_coef_test = t(t(matrix(rep(tau_eigens$eigen_vals1, each = N_test), nrow = N_test))*
                        (t(tau_eigens$eigen_funcs1)%*%sigma_inv%*%d_test))
  
  return(list(score_coef, score_coef_test))
  
}

# Here I define various scenarios:

# Scenarios inspired by Delaigle and Hall(2012)

# Define latent process X_i using Fourier basis
# X1 is length_timepoints x n
generate_latent_process <- function(theta, mu_coeff, tt, n=20){
  # tt<- seq(0,1, len=101)
  K <- length(theta); Khalf <- round(K/2)
  Kseqeven <- 2*(1:Khalf); Kseqodd<- Kseqeven-1
  
  Phitt1 <- sapply(c(1:Khalf), function(j) sqrt(2) * sin(2*pi*j*tt))
  Phitt2 <- sapply(c(1:Khalf), function(j) sqrt(2) * cos(2*pi*j*tt))
  Phitt <- matrix(0, ncol=2*Khalf, nrow=length(tt))
  Phitt[,Kseqeven] <- Phitt1
  Phitt[,Kseqodd] <- Phitt2
  Phitt <-cbind(rep(1, length(tt)), Phitt)[, 1:K]
  
  Z <- matrix(rnorm (K*n), ncol=n)
  Lambda_half <- diag(sqrt(theta))
  
  mu<- rep(0, K)
  mu[1:length(mu_coeff)] <- mu_coeff
  X <- Phitt%*% (Lambda_half%*% Z + mu)
  X
}

# Create the binary valued curves from latent curves
generate_binary_fns <- function(X){
  length_tt <- nrow(X)
  probs<-  1/(1+exp(-X))
  out <- matrix(rbinom(length(X), 1,  c(probs)), nrow=length_tt)
  out
}


fpca_estimate_functions_dai <-  function(Curves, Classes, pve=0.98, fix_number_of_functions = 0){
  
  Ks1 = var(Curves[Classes==1,])
  Ks2 = var(Curves[Classes==2,])
  Classes_train = Classes
  priors =  c(table(Classes_train)/length(Classes_train))
  K_s = Ks1*priors[1]+Ks2*priors[2]
  D = dim(K_s)[1]
  #Ys = cbind(as.vector(K_s), rep(1:D, D), as.vector(t(kronecker(matrix(1, ncol=D) , 1:D) )))
  #remove diagonals to remove excess error
  #Ys = Ys[which(Ys[,2]!=Ys[,3]), ]
  #K_s2 = smooth.2d(Ys[,1], ind = Ys[, 2:3] , nrow = D, ncol= D, theta = 4)
  #K_s = K_s2$z
  K_s = make_pos_semi_def(K_s)
  fpca_results = estimate_eigenfunctions(K_s, pve, fix_number_of_functions)
  fpca_results$eigen_funcs1 = fpca_results$eigen_funcs1*sqrt(D)
  fpca_results$eigen_vals1 = fpca_results$eigen_vals1/D
  sigma = mean(diag(var(Curves)-K_s))
  
  tt = seq(0,1, len=D)
  
  z1 = as.vector(colMeans(Curves))
  gam1 <- gam(z1~s(tt, bs = "cr", m=2, k = 10), 
              family="gaussian", method = "REML")
  mu_hat = gam1$fitted.values
  #fit = smooth.spline(as.vector(colMeans(Curves)))
  #mu_hat = fit$y
  
  return(list(eigen_vals=fpca_results$eigen_vals1,
              eigen_funcs=fpca_results$eigen_funcs1,
              sigma = sigma, mu_hat = mu_hat))
  
}


qda_updated2 = function(scores, classes, prior_g, scores_test){
  
  train_model = qda(scores, classes, prior = prior_g)
  test_dat = predict(train_model, scores_test)
  guess = apply(test_dat$posterior, 1, which.max)
  return(guess)
  
  
}

nb_updated2 = function(scores, classes, prior_g, scores_test){
  
  naive_bayes_dat = rbind.data.frame(scores, scores_test)
  naive_bayes_dat = cbind.data.frame(classes = c(classes, rep(NA, dim(scores_test)[1])),
                                     naive_bayes_dat)
  naive_bayes_dat$classes = factor(naive_bayes_dat$classes)
  
  N_train = length(classes)
  
  naive_bayes_dat_train = naive_bayes_dat[1:N_train,]
  naive_bayes_dat_test = naive_bayes_dat[-(1:N_train),-1]
  rf1 = naive_bayes(classes~., data = naive_bayes_dat_train, usekernel = T, bw = "nrd0")
  guess = predict(rf1, naive_bayes_dat_test, type = "class")
  
  return(guess)
  
}

get_cdf_den = function(densities, scores, classes, i){
  
  score_cur = scores[i,]
  class_cur = classes[i]
  
  dens_cur = densities[[class_cur]]
  cdf.vals = rep(NA, length(dens_cur))
  for(k in 1:length(dens_cur)){
    
    dens_cur_k = dens_cur[[k]]
    cdf.vals[k] = cumsum(dens_cur_k$y)[which.min(abs(dens_cur_k$x-score_cur[k]))]
    cdf.vals[k] = cdf.vals[k]/sum(dens_cur_k$y)
    
  }
  return(cdf.vals) 
}

get_pdf_den = function(densities, scores, classes, i){
  
  score_cur = scores[i,]
  class_cur = classes[i]
  
  dens_cur = densities[[class_cur]]
  cdf.vals = rep(NA, length(dens_cur))
  for(k in 1:length(dens_cur)){
    
    dens_cur_k = dens_cur[[k]]
    cdf.vals[k] = dens_cur_k$y[which.min(abs(dens_cur_k$x-score_cur[k]))]*(max(dens_cur_k$x)-min(dens_cur_k$x))/length(dens_cur_k$x)
    
  }
  return(cdf.vals) 
}

get_pdf_den2 = function(densities, scores, classes, i){
  
  score_cur = scores[i,]
  class_cur = classes[i]
  dens_cur = densities[[class_cur]]
  cdf.vals = rep(NA, length(dens_cur))
  
  for(k in 1:length(dens_cur)){
    
    dens_cur_k = dens_cur[[k]]
    
    approx_fun_den = approxfun(dens_cur_k)
    xnew = score_cur[k]
    if(xnew<=min(dens_cur_k$x)){
      xnew=min(dens_cur_k$x)
    }
    if(xnew>=max(dens_cur_k$x)){
      xnew=max(dens_cur_k$x)
    }
    
    cdf.vals[k]  = approx_fun_den(xnew)
    
  }
  return(cdf.vals) 
}

#function to estimate the predicted values when h is known
nb_updated = function(scores, classes, prior_g, scores_test, h = 1.06){
  
  #grid search to find h
  
  if(length(scores)==length(classes)){
    scores = matrix(scores, ncol = 1) 
    scores_test = matrix(scores_test, ncol = 1) 
  }
  
  nd = dim(scores)[2]
  densities = list()
  for(i in 1:2){
    densities[[i]]=list()
    for(k in 1:nd){
      densities[[i]][[k]] = density(scores[classes==i,k],
                                    kernel = "gaussian", 
                                    bw = h*sd(scores[classes==i,k]))
    }
  }
  
  n_test = dim(scores_test)[1]
  p.mat = matrix(NA, nrow = n_test, ncol = length(prior_g))
  vec = matrix(1:n_test, ncol = 1)
  for(i in 1:2){
    
    pdf.vals.test.1 = t(apply(vec, 1,
                              function(x)  get_pdf_den2(densities, scores_test, rep(i, n_test), x)))
    pdf_vals = apply(pdf.vals.test.1, 1, prod)
    p.mat[,i] = prior_g[i]* pdf_vals
    
  }
  guess = apply(p.mat, 1, which.max)
  
  return(guess)
  
}

#grid search to estimate predicted values and estimate the values of h
nb_updated_grid = function(scores, classes, prior_g, scores_test,
                           grid_search = T, min.h = 0.01, max.h = 1,
                           CV = 10, n_grid = 25, return_h = F, return_prob = F){
  
  vec.cv = matrix(1:CV, ncol = 1)
  #grid search to find h
  
  #add extra incase unequal distribution of groups 
  cvgroups = rep(1:CV, (length(classes)/CV+1))
  #remove unneeded values
  cvgroups = cvgroups[1:length(classes)]
  #randomly assign groups
  cvgroups = sample(cvgroups, size = length(classes), replace = F)
  
  if(length(scores)==length(classes)){
    scores = matrix(scores, ncol = 1) 
    scores_test = matrix(scores_test, ncol = 1) 
  }
  
  get.cv.h = function(h.val){
    groups.probs = 
      apply(vec.cv, 1, function(x) 
        mean(#get guess
          nb_updated(scores[cvgroups!=x,], classes[cvgroups!=x], 
                     c(table(classes[cvgroups!=x])/length(classes[cvgroups!=x])) , #define the new prior probs 
                     scores[cvgroups==x,], h = h.val) ==  classes[cvgroups==x]))
    #get cv.acc
    mean(groups.probs)
  }
  
  grid.vals.h = matrix(seq(min.h, max.h, length.out = n_grid), ncol = 1)
  h.accs = apply(grid.vals.h, 1, function(h) get.cv.h(h))
  
  h = grid.vals.h[which.max(h.accs)]
  #h = grid.vals.h[max(which(h.accs==max(h.accs)))]
  
  
  nd = dim(scores)[2]
  ng = length(unique(classes))
  densities = list()
  for(i in 1:ng){
    densities[[i]]=list()
    for(k in 1:nd){
      densities[[i]][[k]] = density(scores[classes==i,k],
                                    kernel = "gaussian", 
                                    bw = h*sd(scores[classes==i,k]))
    }
  }
  
  n_test = dim(scores_test)[1]
  p.mat = matrix(NA, nrow = n_test, ncol = length(prior_g))
  vec = matrix(1:n_test, ncol = 1)
  for(i in 1:ng){
    
    pdf.vals.test.1 = t(apply(vec, 1,
                              function(x)  get_pdf_den2(densities, scores_test, rep(i, n_test), x)))
    pdf_vals = apply(pdf.vals.test.1, 1, prod)
    p.mat[,i] = prior_g[i]* pdf_vals
    
  }
  
  if(return_prob){
    return(p.mat/rowSums(p.mat))
  }
  
  guess = apply(p.mat, 1, which.max)
  
  #print(h)
  
  if(return_h){
    return(h)
  }
  return(guess)
  
}



regression_f = function(z, dta, lm_structure){
  lm2 <- lm(lm_structure, data = subset(dta, dta$id==z))
  return(lm2$coefficients)
}

regression_g = function(z, Curves, tt, k=10, method="REML"){
  z1 = Curves[z,]
  gam1 <- gam(z1~s(tt, bs = "cr", m=2, k = k), 
              family="binomial", method = method)
  return(gam1$fitted.values)
}

regression_bf = function(z, dta, lm_structure, prior_scales_test){
  bayesglm1 = bayesglm(lm_structure,
                       family = gaussian,
                       data = subset(dta, dta$id==z),
                       prior.scale = sqrt(prior_scales_test),
                       prior.df = Inf, scaled = F)
  return(bayesglm1$coefficients)
}

regression_bf2 = function(z, dta, lm_structure, prior_scales_test){
  bayesglm1 = bayesglm(lm_structure,
                       family = binomial(),
                       data = subset(dta, dta$id==z),
                       prior.scale = sqrt(prior_scales_test),
                       prior.df = Inf, scaled = F)
  return(bayesglm1$coefficients)
}


eigen_funcs1 = fpca.cur$efunctions
eigen_funcs_true = eigen_funcs1
mu_t_hat = fpca.cur$mu
eigen_vals = fpca.cur$evalues
eigen_vals_true = fpca.cur$evalues
mu_t_hat_true = mu_t_hat
scores_train_true = scores_train

#h.val = nb_updated_grid(scores_train, Classes_train, prior_g, scores_test, return_h = T )
#from previous analaysis
h.val=0.09
dim(scores_train)
data.train = data.frame(scores_train, Classes_train)
#data.test = data.frame(scores_test, Classes_test)
dens.dat = list()
ndf = dim(scores_train)[2]
for(i in 1:ndf){
  
  psi_class1 = data.train[,i][data.train$Classes_train==1]
  psi_class2 = data.train[,i][data.train$Classes_train==2]
  
  dens1 = density(psi_class1, kernel = "gaussian", from = min(psi_class1, psi_class2)-0.25*sd(psi_class2), 
                  to = max(psi_class1, psi_class2)+0.5*sd(psi_class2), n = 500, 
                  bw = h.val*sd(psi_class1))
  dens2 = density(psi_class2, kernel = "gaussian", from = min(psi_class1, psi_class2)-0.25*sd(psi_class2), 
                  to = max(psi_class1, psi_class2)+0.5*sd(psi_class2), n = 500, 
                  bw = h.val*sd(psi_class2))
  
  dens.dat[[i]] = data.frame(score1 = dens1$x, d1.not = dens1$y, d1.bot = dens2$y)
}

get.scores = function(i, n1,n2){
  
  cdf.vals1 = runif(n1, 0, 1)
  cdf.vals2 = runif(n2, 0, 1)
  dens.dat.cur = dens.dat[[i]]
  
  dens.dat.cur$cdf1 = cumsum(dens.dat.cur$d1.bot)/sum(dens.dat.cur$d1.bot)
  dens.dat.cur$cdf2 = cumsum(dens.dat.cur$d1.not)/sum(dens.dat.cur$d1.not)
  
  scores.n1 = sapply(cdf.vals1, function(x) dens.dat.cur$score1[which.min(abs(x-dens.dat.cur$cdf1))])
  scores.n2 = sapply(cdf.vals2, function(x) dens.dat.cur$score1[which.min(abs(x-dens.dat.cur$cdf2))])
  
  return(c(scores.n1, scores.n2))
  
}

