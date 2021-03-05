##################
#Functions for "Classification of Social Media Users Using a Generalized Functional Analysis"
#
#Author: Anthony Weishampel
#Date: 2/27/2021
#
##################

#dependencies to run all of the code
library(parallel)
library(xtable)
library(fda)
library(refund)
library(Matrix)
library(MASS)
library(arm)
library(mgcv)
library(readr)

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
    epsilon = matrix(rnorm(N*D, 0, sigma), ncol=D, nrow=N)
    
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
    # dens.dat is a global variable which estimates the 
    scores1 = sapply(is, function(x) get.scores(dens.dat, x, n1, n2))
    
    #mean function globally defined
    mu_hat = mu_t_hat_true
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
#Function to estimate eiginfunctions
#
#Inputs: 
# K_b: estimate for covariance matrix
# pve: Proportion of variance explained
# fix_num_of_functions: set the number of eigenfunctions to be returned
#
#Outputs: 
# eigenvalues and eigenvectors
#
#####
estimate_eigenfunctions = function(K_b, pve=0.98, fix_num_of_functions=0){
  
  #SVD on matrix
  svd_kb = svd2(K_b)
  
  #Only get the pve length if the number of functions is not given
  if(fix_num_of_functions==0){
    #get number of components
    pb = get_length_pve(pve, svd_kb$d)
  }else{
    pb = fix_num_of_functions
  }
  pb = ifelse(pb>1, pb, 2)
  eigen_vals1 = svd_kb$d[1:pb]
  eigen_funcs1 = svd_kb$v[,1:pb]
  
  return(list(eigen_vals1=eigen_vals1, eigen_funcs1=eigen_funcs1))
  
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


#####
#Function: To get the density values of new scores for an individual in a given class
#
#Inputs: 
# densities: List of density lists 
# scores: mtrix of scores for the N individual 
# classes: matrix of classes (these are not the true classes but the value to be evaluated)
# i: which individual to evaluate
#
#Output: 
# Density values of the new scores for individual i
#
#####
get_pdf_den2 = function(densities, scores, classes, i){
  
  #get score for current user
  score_cur = scores[i,]
  #get class for current user
  class_cur = classes[i]
  #get density for current class
  dens_cur = densities[[class_cur]]
  
  pdf.vals = rep(NA, length(dens_cur))
  #for each component
  for(k in 1:length(dens_cur)){
    
    #get K comp density
    dens_cur_k = dens_cur[[k]]
    
    #approximate the function
    approx_fun_den = approxfun(dens_cur_k)
    
    #get new score
    xnew = score_cur[k]
    
    #if new value is outside of defined range
    if(xnew<=min(dens_cur_k$x)){
      xnew=min(dens_cur_k$x)
    }
    #if new vlaue is outside of defined range
    if(xnew>=max(dens_cur_k$x)){
      xnew=max(dens_cur_k$x)
    }
    
    #get value in the desnity
    pdf.vals[k]  = approx_fun_den(xnew)
    
  }
  return(pdf.vals) 
}

#####
#Function: Grid search to estimate predicted values and estimate the values of h
#
#Inputs: 
# scores: N x K matrix of scores in the training set
# classes: Group labels vector of length N
# prior_g: vector of prior probability of being in each group, sums up to 1
# scores_test: N_test x K matrix of scores in the testing set
# h: multiplier for kernel based density function
#
#Output: 
# predicted classes for the Bayes Classifier
#
#####
nb_updated = function(scores, classes, prior_g, scores_test, h = 1.06){
  
  #in case scores are N x 1 matrix that is read as a vector 
  if(length(scores)==length(classes)){
    scores = matrix(scores, ncol = 1) 
    scores_test = matrix(scores_test, ncol = 1) 
  }
  
  #get K 
  nd = dim(scores)[2]
  #get number of groups
  ng = length(unique(classes))
  #initialize list
  densities = list()
  #estimate density at each K for each group
  for(i in 1:ng){
    densities[[i]]=list()
    for(k in 1:nd){
      densities[[i]][[k]] = density(scores[classes==i,k],
                                    kernel = "gaussian", 
                                    bw = h*sd(scores[classes==i,k]))
    }
  }
  
  # get number of test functions
  n_test = dim(scores_test)[1]
  p.mat = matrix(NA, nrow = n_test, ncol = length(prior_g))
  vec = matrix(1:n_test, ncol = 1)
  #For each group get the Bayes classifier value of probability of being in that group 
  for(i in 1:ng){
    #apply for each user get the probability for each component in that group
    pdf.vals.test.1 = t(apply(vec, 1,
                              function(x)  get_pdf_den2(densities, scores_test, rep(i, n_test), x)))
    #apply for each user to get the product of those K probabilities
    pdf_vals = apply(pdf.vals.test.1, 1, prod)
    #multiply by prior probability
    p.mat[,i] = prior_g[i]* pdf_vals #p.mat is now matrix of Bayes probabilities
  }
  
  #get guess based on which is the max
  guess = apply(p.mat, 1, which.max)
  
  return(guess)
  
}


#####
#Function: Grid search to estimate predicted values and estimate the values of h
#
#Inputs: 
# scores: N x K matrix of scores in the training set
# classes: Group labels vector of length N
# prior_g: vector of prior probability of being in each group, sums up to 1
# scores_test: N_test x K matrix of scores in the testing set
# min.h: min possible value for the multiplier
# max.h: maximum possible value for the multiplier
# n_grid: number of vlaues between min.h and max.h to search over
# CV: Number of folds for cross validation
# return_h: T/F to return the value of the multiplier
# return_prob: T/F to return group Bayes classifier probability for each individual in testing set
#
#Output: 
# predictions from the grid search
#
#####
nb_updated_grid = function(scores, classes, prior_g, scores_test,
                           min.h = 0.01, max.h = 1,
                           CV = 10, n_grid = 25, return_h = F, return_prob = F){
  
  #create matrix for apply functions
  vec.cv = matrix(1:CV, ncol = 1)
  

  #create vector of possible CV groups to each account 
  #add extra incase unequal distribution of groups
  cvgroups = rep(1:CV, (length(classes)/CV+1))
  #remove the unneeded values
  cvgroups = cvgroups[1:length(classes)]
  #randomly assign CV group to each account
  cvgroups = sample(cvgroups, size = length(classes), replace = F)
  
  #in case scores are N x 1 matrix that is read as a vector 
  if(length(scores)==length(classes)){
    #If k == 1 change to vector
    scores = matrix(scores, ncol = 1) 
    scores_test = matrix(scores_test, ncol = 1) 
  }
  
  # define function here to use the cvgroups this function 
  # will be used in the following apply statement
  get.cv.h = function(h.val){
    groups.probs = 
      apply(vec.cv, 1, function(x) 
        mean(#get guess using the nb_updated function
          nb_updated(scores[cvgroups!=x,], classes[cvgroups!=x], 
                     c(table(classes[cvgroups!=x])/length(classes[cvgroups!=x])) , #define the new prior probs 
                     scores[cvgroups==x,], h = h.val) ==  classes[cvgroups==x]))
    #return the accuracy of the prediction
    mean(groups.probs)
  }
  
  #initialize matrix for apply which contains all of the possible grid values
  grid.vals.h = matrix(seq(min.h, max.h, length.out = n_grid), ncol = 1)
  #apply the previously defined functions to get the CV accuracies at each h value
  h.accs = apply(grid.vals.h, 1, function(h) get.cv.h(h))
  
  # assign h value based on the one with the largest CV accuracy
  h = grid.vals.h[which.max(h.accs)]
  #h = grid.vals.h[max(which(h.accs==max(h.accs)))]
  
  #get K 
  nd = dim(scores)[2]
  #get number of groups
  ng = length(unique(classes))
  #initialize list
  densities = list()
  #estimate density at each K for each group
  for(i in 1:ng){
    densities[[i]]=list()
    for(k in 1:nd){
      densities[[i]][[k]] = density(scores[classes==i,k],
                                    kernel = "gaussian", 
                                    bw = h*sd(scores[classes==i,k]))
    }
  }
  
  # get number of test functions
  n_test = dim(scores_test)[1]
  p.mat = matrix(NA, nrow = n_test, ncol = length(prior_g))
  vec = matrix(1:n_test, ncol = 1)
  #For each group get the Bayes classifier value of probability of being in that group 
  for(i in 1:ng){
    #apply for each user get the probability for each component in that group
    pdf.vals.test.1 = t(apply(vec, 1,
                              function(x)  get_pdf_den2(densities, scores_test, rep(i, n_test), x)))
    #apply for each user to get the product of those K probabilities
    pdf_vals = apply(pdf.vals.test.1, 1, prod)
    #multiply by prior probability
    p.mat[,i] = prior_g[i]* pdf_vals #p.mat is now matrix of Bayes probabilities
  }
  
  #returns matrix of probabilies for each group
  if(return_prob){
    return(p.mat/rowSums(p.mat))
  }
  
  #group prediction is based on maximum posterior probability
  guess = apply(p.mat, 1, which.max)
  
  #print(h)

  if(return_h){
    return(h)
  }
  
  #return guesses
  return(guess)
  
}


#####
#Function: wrapper function for gam() which outputs the fitted values
#
#Inputs: 
# z : index z = 1,...,N 
# Curves : N x D matrix of observed binary series
# tt : grid of timepoints going from 0 to 1 with D observations
# k : number of basis functions
# method: method used to evaluate the gam
#
#Output: 
# Fitted values from the game function for subject z 
#
#####
regression_g = function(z, Curves, tt, k=10, method="REML"){
  z1 = Curves[z,]
  gam1 <- gam(z1~s(tt, bs = "cr", m=2, k = k), 
              family="binomial", method = method)
  return(gam1$fitted.values)
}


#####
#Function: wrapper function for bayesglm() which outputs expected coefficients for the gaussian responses
#
#Inputs: 
# z : index z = 1,...,N 
# dta : data frame contain the subject id, mean function, eigenfunctions and observe binary values
# lm_structure: formula displaying the bayesglm function
# eigen_vals1: eigen_values for the psi coefficients
#
#Output: 
# Fitted values from the game function for subject z 
#
#Notes: This function is used in simulation scenario Spline+FPCA
#####
regression_bf = function(z, dta, lm_structure, eigen_vals1){
  bayesglm1 = bayesglm(lm_structure,
                       family = gaussian,
                       data = subset(dta, dta$id==z),
                       prior.scale = sqrt(eigen_vals1),
                       prior.df = Inf, scaled = F)
  return(bayesglm1$coefficients)
}

#####
#Function: wrapper function for bayesglm() which outputs expected coefficients for the gaussian responses
#
#Inputs: 
# z : index z = 1,...,N 
# dta : data frame contain the subject id, mean function, eigenfunctions and observe binary values
# lm_structure: formula displaying the bayesglm function
# eigen_vals1: eigen_values for the psi coefficients
#
#Output: 
# Fitted values from the game function for subject z 
#
#####
regression_bf2 = function(z, dta, lm_structure, eigen_vals1){
  bayesglm1 = bayesglm(lm_structure,
                       family = binomial(link = "logit"),
                       data = subset(dta, dta$id==z),
                       prior.scale = sqrt(eigen_vals1), #set scales
                       prior.df = Inf, #normal priors
                       scaled = F ) #Do not adjust the scales
  return(bayesglm1$coefficients)
}


#####
#Function: Randomly draw new scores from the two density distributions 
#
#Inputs: 
# dens.dat: list containing the following information about the two densities
#       dens.dat[[k]] = list containing the information on the kth component
#       For each k: Data.frame defined as
#             dens.dat[[k]]$score: vector containing defining the x values two groups
#             dens.dat[[k]]$d1.not: vector of same length containing the density values for not bots
#             dens.dat[[k]]$d1.bot: vector of same length containing the density values for bots
# k: which component to draw from k = 1,...,K
# n1: number of bots to draw from
# n2: number of genuine accounts to draw from
#
#Outputs: 
# vector of scores for the kth component (bot scores, genuine scores)
#
#####
get.scores = function(dens.dat, k, n1, n2){
  
  #randomly draw from uniform distributions
  q.vals1 = runif(n1, 0, 1)
  q.vals2 = runif(n2, 0, 1)
  # get the estimated density 
  # dens.dat is a global variable
  dens.dat.cur = dens.dat[[k]]
  
  #estimate empirical cdfs
  dens.dat.cur$cdf1 = cumsum(dens.dat.cur$d1.bot)/sum(dens.dat.cur$d1.bot)
  dens.dat.cur$cdf2 = cumsum(dens.dat.cur$d1.not)/sum(dens.dat.cur$d1.not)
  
  #Recover scores for each group 
  #Find the closest value of the random uniform dist to the empirical cdf 
  #   and the score corresponding to that value 
  scores.n1 = sapply(q.vals1, function(x) dens.dat.cur$score[which.min(abs(x-dens.dat.cur$cdf1))])
  scores.n2 = sapply(q.vals2, function(x) dens.dat.cur$score[which.min(abs(x-dens.dat.cur$cdf2))])
  
  return(c(scores.n1, scores.n2))
  
}



#####
#Function: Method for extracting scores (training an testing sets) 
# based on methodology presented in Hall 2008
#
#Inputs: 
# Curves_binary: N x M matrix of binary observations for the training set accounts
# Curves_binary_test: N_test x M matrix of binary observations for the testing set accounts
#
#Outputs: 
# scores for both training and testing set accounts
#
#Notes: This is only used at the mFPCA comparable
#####
GFPCA_estimate_para_with_test <- function(Curves_binary, Curves_binary_test){
  
  D = dim(Curves_binary)[2]
  N = dim(Curves_binary)[1]
  N_test = dim(Curves_binary_test)[1]
  tt = seq(0,1, len=D)
  
  #estimate smooth mean function
  z1 = as.vector(colMeans(Curves_binary))
  gam1 <- gam(z1~s(tt, bs = "cr", m=2, k = 10),
              family="gaussian", method = "REML")
  #because probability needs to ensure values are between 0 and 1
  alpha_t = ifelse(gam1$fitted.values<0.001, 0.001, gam1$fitted.values)
  alpha_t = ifelse(alpha_t>0.999, 0.999, alpha_t)
  
  #from Hall 2008 paper 
  v_t_hat = as.vector(invlogit(alpha_t))
  alpha_t_hat = as.vector(alpha_t)
  beta_ts  = t(Curves_binary)%*%Curves_binary/N
  
  #Format for smooth function 
  Bs = cbind(as.vector(beta_ts), rep(1:D, D), as.vector(t(kronecker(matrix(1, ncol=D) , 1:D) )))
  #Bs = Bs[which(Bs[,2]!=Bs[,3]), ]
  #remove diagonal values
  Bs[which(Bs[,2]==Bs[,3]), 1] = NA
  Bs[,2:3] = Bs[,2:3]/D
  Bs = data.frame(val = Bs[,1], t1 = Bs[,2], t2 = Bs[,3])
  #smooth the covariance matrix, need large theta because we have lots of data points
  beta_ts2 = gam(val~te(t1, t2, 
                        bs = "cr", m=2, k=5), 
                 method = "REML", family="gaussian", data = Bs)
  newdat = data.frame(t1 = (1:D)/D, t2 = (1:D)/D)
  diag.vals = predict(beta_ts2, newdata = newdat)
  return.mat = matrix(NA, nrow = D, ncol = D)
  diag(return.mat) = diag.vals
  return.mat[which(Bs[,2]!=Bs[,3])] = beta_ts2$fitted.value
  beta_ts_hat = return.mat  
  
  tau_ts = (beta_ts_hat - matrix(alpha_t_hat, ncol = 1) %*% t(alpha_t_hat))/
    (matrix(d_logit(v_t_hat), ncol = 1)%*%t(d_logit(v_t_hat)))
  
  #make sym and semi-post def
  tau_ts_hat = round(make_pos_semi_def(tau_ts), 5)
  mu_t_hat = v_t_hat
  
  #decomposition of resulting covariance matrix
  #tau_eigens  =  estimate_eigenfunctions(tau_ts_hat, pve = 0.98)
  tau_eigens = estimate_eigenfunctions(tau_ts_hat, pve=0.98)
  tau_eigens$eigen_funcs1 = tau_eigens$eigen_funcs1*sqrt(D)
  tau_eigens$eigen_vals1 = tau_eigens$eigen_vals1/D
  
  
  #estimate from Hall 2008
  d = t((Curves_binary - matrix(rep(logit(mu_t_hat), each = N), ncol = D)) /
          matrix(rep(d_logit(mu_t_hat), each = N), ncol = D))
  
  #it's small enough to use solve in simulation scenarios
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


#####
#Function: To apply the methods form Dai et al 2017
#
#Inputs: 
# Curves: N x M observed for N individuals at M time points. This are Gaussian observations not binary
# Classes: Vector of classes for each individual
# pve: Proportion of variance explained 
# fix_number_of_functions: If user wanted to fix number of eigenfunctions and not use PVE
#
#Output: 
# Parsimonious representation for the functions,
# i.e. mean function, eigenfunction, variance, eigenvalues
#
# Notes: This is used for the Latent Curves classifier in the simulation study
#####
fpca_estimate_functions_dai <-  function(Curves, Classes, pve=0.98, 
                                         fix_number_of_functions = 0){
  
  #Get sample  estimates of the covariance matrices for each group as per paper
  Ks1 = cov(Curves[Classes==1,])
  Ks2 = cov(Curves[Classes==2,])
  Classes_train = Classes
  #get prior distributions
  priors =  c(table(Classes_train)/length(Classes_train))
  #Estimate pooled covariance
  K_s = Ks1*priors[1]+Ks2*priors[2]
  D = dim(K_s)[1]
  
  #ensure estimated matrix is semi-positive definite
  K_s = make_pos_semi_def(K_s)
  #get eigenfunctions for these matrix
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








