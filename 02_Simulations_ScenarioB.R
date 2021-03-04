##################
#Simulation for Scenario 2 of 
#"Classification of Social Media Users Using a Generalized Functional Analysis"
#
#Author: Anthony Weishampel
#Date: 2/27/2021
#
# Make sure to run 01_functions file first
##################


#Needs to be changed for the three values of D
D = 150

#define the grid
grid = seq(0, 1, length = D)

#Number of users in the testing sets
Ns = c(50, 300, 1000)

#Number of cores able to run code on
numCores_to_run = 16

#for parallel processing
numCores <- detectCores() # get the number of cores available


#number of iterations
ITER = 500

#number of curves in the testing set
N_test = 500

#set seed
set.seed(1234)

###define the lists and matrices for each scenario

#Records the miclassification rates
misclass_list = list()
misclass_mat = matrix(NA, nrow = 5,  ncol = ITER)

#Records the times
time_mat = matrix(NA, nrow = 5,  ncol = ITER)
time_list = list()

#measures the precision
precision_mat = matrix(NA, nrow = 5,  ncol = ITER)
precision_list = list()

#Records the sensitivity
sens_mat = matrix(NA, nrow = 5,  ncol = ITER)
sens_list = list()

#Records the specificity
spec_mat = matrix(NA, nrow = 5,  ncol = ITER)
spec_list = list()

####
#Function: Main Simulation function to pass through to different cores
#
#Input: 
# i: iteration to run simulation one
#
#Outputs: 
# misclassification, time, precision, sensitivity, and specificity for each tested classifier on iteration i
####
core_function = function(i){
  
  #Need to seed for each different core to ensure same results 
  set.seed(i)
  
  #matrices to return one row for each tested classifier
  misclass_mat = matrix(NA, nrow = 5,  ncol = 1)
  time_mat = matrix(NA, nrow = 5,  ncol = 1)
  precision_mat = matrix(NA, nrow = 5,  ncol = 1)
  sens_mat = matrix(NA, nrow = 5,  ncol = 1)
  spec_mat = matrix(NA, nrow = 5,  ncol = 1)

  #output to know how much is left
  if(i%%25 == 0 || i==1){
      print(paste("On Iteration ", i , " out of ", ITER))
  }
  
  #N2 is the global variable for number of 
  N=N2
  
  #Needed for generating data ensures that we get 
  # same number of functions from each group. we 
  # assign the
  p=c(0.5,0.5)
  
  #define grid of timepoints which curves are observed on
  grid = seq(0,1, len=D)
  tt=seq(0,1, len=D)

  #generate N*2 (N for each group) curves for the given scenario, we are generating one curve from each class per individual and then assigning class after
  Curves = generate_data(scenario = 2, grid = grid,  N=2*N,  p = p, 
                           sigma = 0, binary = T)

  #organize generated curves by class
  Y0 = t(Curves[1:(dim(Curves)[1]/2),])
  Y1 = t(Curves[-(1:(dim(Curves)[1]/2)),])
  
  #generate N_test*2 (N_test for each group) curves for the given scenario  
  Curves = generate_data( scenario = 2, grid = grid,  N=2*N_test,  p = p,
                                    sigma = 0, binary = T)
  
  #organize generated testing curves by class
  Y0_test = t(Curves[1:(dim(Curves)[1]/2),])
  Y1_test = t(Curves[-(1:(dim(Curves)[1]/2)),])
  
  #temp assignment curves for training and testing 
  # currently listed as first N as group 0 and second N as group 1
  Curves_train = t(cbind(Y0, Y1))
  Curves_test = t(cbind(Y0_test, Y1_test))
  
  #randomly assign N/N_test individuals into the two classes 
  # each class has a 0.50 probability of being observed
  # +1 is added for formatting of curves in the next step 
  Classes_train = rbinom(N, 1, 0.5)+1
  Classes_test = rbinom(N_test, 1, 0.5)+1
  
  #Get the curves based on the classes 
  # for individual i: 
  #     if from group 0 get the ith row from of the curves matrix
  #     if from group 1 get the 2*ith row from the curves matrix
  Curves_train = Curves_train[(1:N)+N*(Classes_train-1),]
  Curves_test = Curves_test[(1:(N_test))+(N_test)*(Classes_test-1),]
  
  
  N = dim(Curves_train)[1]/2
  
  #treat the categorical classes as factors and not numerical values  
  classes_train = factor(Classes_train)
  classes_test = factor(Classes_test)
  
  #N_train = N but just easier for following from now on  
  N_train = length(Classes_train)
  
  ####
  #Code for the proposed Method
  # The following steps are for the proposed method as outlined in the paper
  ###
  
  #get start time
  st = Sys.time()
  
  tt=seq(0,1, len=D)
  
  ##
  #Step 1 of the proposed method 
  ##
  vec = matrix(1:(N_train), ncol = 1)
  smoothed_x = logit(t(apply(vec, 1, function(x) regression_g(x, Curves_train, tt))))
  
  ##
  #Step 2 of the proposed method
  ##
  fpca.cur2 = fpca.face(smoothed_x, pve = 0.98, p=3, m=2, knots = 6) #lambda selected via grid search optim, #p=degree of splines
  #correct fpca.cur2 eigenval values b/c too large
  #get multiplier is function of D in fact its 1/D
  #get_multiplier = sum((fpca.cur2$efunctions[,1])^2/length(fpca.cur2$efunctions[,1]))
  get_multiplier = 1/D
  fpca.cur = fpca.cur2
  #correct eigenfunctions
  fpca.cur$efunctions = fpca.cur2$efunctions/sqrt(get_multiplier) 
  #correct eigenvalues
  fpca.cur$evalues = fpca.cur2$evalues*get_multiplier
  #correct scores
  fpca.cur$scores = fpca.cur2$scores*sqrt(get_multiplier)
  
  ##
  #STEP 3: 
  # Set up and apply Bayesglm framework
  ##
  fit = list(mu = fpca.cur$mu,
             evalues = fpca.cur$evalues,
             efunctions = fpca.cur$efunctions)
  
  mu_t_hat = fit$mu
  eigen_vals1 = fit$evalues
  eigen_funcs1 = fit$efunctions
  
  #data frame used in bayesglm
  dta = data.frame(index = rep(tt, N_train),
                   value = c(t(Curves_train)),
                   id = rep(1:N_train, each = D))
  
  npc = length(eigen_vals1)
  if(npc>1){
    for (z in 1:npc) {
      dta <- cbind(dta, rep(eigen_funcs1[,z], N_train))
    }
  }else{
    dta = cbind(dta, matrix(eigen_funcs1, ncol =1))
  }
  
  #assign names to data frame
  names(dta)[4:(4 + npc - 1)] <- c(paste0("psi", 1:npc))
  #repeat mean function in data frame once per user
  dta$mu = rep(mu_t_hat , N_train)
  
  #get formula for glm
  glm_structure = paste(paste0("psi", 1:npc), collapse = "+")
  glm_structure = paste("value ~ -1 + offset(mu) +" , glm_structure , sep="")
  #set scale for the glm 
  prior_scales_test = eigen_vals1
  
  #Estimate the Scores for the training set
  vec = matrix(1:N_train, ncol = 1)
  scores_train = t(apply(vec, 1, function(x) regression_bf2(x, dta, glm_structure, prior_scales_test)))
  
  #Step 3 for Testing Data
  #Get the socres for the testing set
  
  #just like before define data frame 
  dta = data.frame(index = rep(tt, N_test),
                   value = c(t(Curves_test)),
                   id = rep(1:N_test, each = D))
  
  npc = length(eigen_vals1)
  
  if(npc>1){
    for (z in 1:npc) {
      dta <- cbind(dta, rep(eigen_funcs1[,z], N_test))
    }
  }else{
    dta = cbind(dta, matrix(eigen_funcs1, ncol =1))
  }
  names(dta)[4:(4 + npc - 1)] <- c(paste0("psi", 1:npc))
  dta$mu = rep(mu_t_hat , N_test)
  
  glm_structure = paste(paste0("psi", 1:npc), collapse = "+")
  glm_structure = paste("value ~ -1 + offset(mu) +" , glm_structure , sep="")
  
  vec = matrix(1:N_test, ncol = 1)
  scores_test = t(apply(vec, 1, function(x) regression_bf2(x, dta, glm_structure, prior_scales_test)))
  
  #step 4 
  #get propability of being in each group
  prior_g = c(table(Classes_train)/length(Classes_train))
  #run non parametric bayes classifier
  pred_classes = nb_updated_grid(scores_train, Classes_train, prior_g, scores_test, min.h = 0.5, max.h = 1.5)

  #clock the time
  et = Sys.time()
  time_diff = et - st
  
  #Determine how well the classifier performed
  t1 = table(pred_classes, Classes_test)
  #record results for proposed method
  misclass_mat[1, 1] = 1- sum(pred_classes==Classes_test)/length(Classes_test)
  sens_mat[1,1] = sum(pred_classes==1 & Classes_test == 1)/(sum(pred_classes==2 & Classes_test == 1)+sum(pred_classes==1 & Classes_test == 1))
  spec_mat[1,1] = sum(pred_classes==2 & Classes_test == 2)/(sum(pred_classes==2 & Classes_test == 2)+sum(pred_classes==1 & Classes_test == 2))
  precision_mat[1,1] = sum(pred_classes==1 & Classes_test == 1)/(sum(pred_classes= 1 & Classes_test == 1)+sum(pred_classes==1 & Classes_test == 2))
  time_mat[1, 1] = as.numeric(time_diff)
  

  ##
  #Method 2: Marginal FPCA based on Hall 2008
  ##
  #get start time
  st = Sys.time()

  #get probability of being in each group
  prior_g = c(table(Classes_train)/length(Classes_train))

  #Get scores for both the training and testing set curves
  gfpca_scores = GFPCA_estimate_para_with_test(Curves_train, Curves_test)

  #Training score estimates
  scores_train = gfpca_scores[[1]]
  #Testing score estimates
  scores_test = gfpca_scores[[2]]

  #Fit the non-parametric Bayes classifier
  pred_classes = nb_updated_grid(scores_train, Classes_train, prior_g, scores_test, min.h = 0.5, max.h = 1.5)

  et = Sys.time()
  time_diff = et - st

  #record results
  misclass_mat[2, 1] = 1- sum(pred_classes==Classes_test)/length(Classes_test)
  sens_mat[2,1] = sum(pred_classes==1 & Classes_test == 1)/(sum(pred_classes==2 & Classes_test == 1)+sum(pred_classes==1 & Classes_test == 1))
  spec_mat[2, 1] = sum(pred_classes==2 & Classes_test == 2)/(sum(pred_classes==2 & Classes_test == 2)+sum(pred_classes==1 & Classes_test == 2))
  precision_mat[2, 1] = sum(pred_classes==1 & Classes_test == 1)/(sum(pred_classes= 1 & Classes_test == 1)+sum(pred_classes==1 & Classes_test == 2))
  time_mat[2, 1] = as.numeric(time_diff)


  ##
  #Method 3: Spline+mFPCA
  ##

  #set start time
  st = Sys.time()

  tt=seq(0,1, len=D)

  #apply step one from proposed classifier to both the training and testing sets
  # since will are building a classifier based on these estimates
  vec = matrix(1:(N_train), ncol = 1)
  smoothed_x = logit(t(apply(vec, 1, function(x) regression_g(x, Curves_train, tt))))

  vec = matrix(1:(N_test), ncol = 1)
  smoothed_x_test = logit(t(apply(vec, 1, function(x) regression_g(x, Curves_test, tt))))

  ##Step 2 for mFPCA
  fpca.cur2 = fpca.face(smoothed_x, pve = 0.98, p=3, m=2, knots = 6) #lambda selected via grid search optim, #p=degree of splines
  #correct fpca.cur2 eigenval values b/c too large
  get_multiplier = 1/D
  fpca.cur = fpca.cur2
  fpca.cur$efunctions = fpca.cur2$efunctions/sqrt(get_multiplier) 
  fpca.cur$evalues = fpca.cur2$evalues*get_multiplier
  fpca.cur$scores = fpca.cur2$scores*sqrt(get_multiplier)
  
  #testing set
  fpca.cur2 = fpca.face(smoothed_x, pve = 0.98, p=3, m=2, knots = 6, 
                        Y.pred = smoothed_x_test) #lambda selected via grid search optim, #p=degree of splines
  #correct fpca.cur2 eigenval values b/c too large
  get_multiplier = 1/D
  fpca.cur.test = fpca.cur2
  fpca.cur.test$efunctions = fpca.cur2$efunctions/sqrt(get_multiplier) 
  fpca.cur.test$evalues = fpca.cur2$evalues*get_multiplier
  fpca.cur.test$scores = fpca.cur2$scores*sqrt(get_multiplier)
  
  
  #Step Trian classifier on scores estimated by FPCA.face function
  scores_train = fpca.cur$scores
  scores_test = fpca.cur.test$scores
  prior_g = c(table(Classes_train)/length(Classes_train))

  pred_classes = nb_updated_grid(scores_train, Classes_train, prior_g, scores_test, min.h = 0.5, max.h = 1.5)
  et = Sys.time()
  time_diff = et - st

  #write results for second method
  misclass_mat[3, 1] = 1- sum(pred_classes==Classes_test)/length(Classes_test)
  sens_mat[3, 1] = sum(pred_classes==1 & Classes_test == 1)/(sum(pred_classes==2 & Classes_test == 1)+sum(pred_classes==1 & Classes_test == 1))
  spec_mat[3, 1] = sum(pred_classes==2 & Classes_test == 2)/(sum(pred_classes==2 & Classes_test == 2)+sum(pred_classes==1 & Classes_test == 2))
  precision_mat[3, 1] = sum(pred_classes==1 & Classes_test == 1)/(sum(pred_classes= 1 & Classes_test == 1)+sum(pred_classes==1 & Classes_test == 2))
  time_mat[3, 1] = as.numeric(time_diff)


  ###
  #Oracle number 1: Binary but real eigenfunctions
  ###
  
  #start time
  st = Sys.time()

  #set the true mean functions
  mu_true0 = rep(0, length(grid))
  mu_true1 = rep(0, length(grid))
  mu_true = mu_true0
  
  #Use the true  eigenfunctions
  k=9
  Ks= 1:k
  lambda_k = exp(-Ks/3)
  psis = matrix(rep(1, D), ncol = 1)
  
  for(z4 in 1:((k-1)/2)){
    p1 = (z4)*2
    psi_1 = sqrt(2)*cos(p1*pi*grid)
    psi_2 = sqrt(2)*sin(p1*pi*grid)
    psis = cbind(psis , psi_1, psi_2)
    
  }
  
  #set up the data frame for the bayesglm classifier
  fit = list(mu = mu_true,
             evalues = lambda_k,
             efunctions = psis)
  mu_t_hat = fit$mu
  eigen_vals1 = fit$evalues
  eigen_funcs1 = fit$efunctions
  
  dta = data.frame(index = rep(grid, N_train),
                   value = c(t(Curves_train)),
                   id = rep(1:N_train, each = D))
  
  npc = length(eigen_vals1)
  if(npc>1){
    for (z in 1:npc) {
      dta <- cbind(dta, rep(eigen_funcs1[,z], N_train))
    }
  }else{
    dta = cbind(dta, matrix(eigen_funcs1, ncol =1))
  }
  names(dta)[4:(4 + npc - 1)] <- c(paste0("psi", 1:npc))
  #dta$mu = rep(mu_t_hat , N_train)
  dta$mu = c(rep(mu_true0, N), rep(mu_true1, N))
  
  glm.structure = paste(paste0("psi", 1:npc), collapse = "+")
  glm_structure = paste("value ~ -1 + offset(mu) +" , glm.structure , sep="")
  
  #for class 1 change the variance
  lambda_k = exp(-Ks/3)
  prior_scales_test = sqrt(lambda_k)
  vec = matrix(which(Classes_train==1), ncol = 1)
  scores_train1 = t(apply(vec, 1, function(x) regression_bf2(x, dta, glm_structure, prior_scales_test)))
  
  #for class 2 change the variance
  lambda_k = exp(-Ks/2)
  prior_scales_test = sqrt(lambda_k)
  vec = matrix(which(Classes_train==2), ncol = 1)
  scores_train2 = t(apply(vec, 1, function(x) regression_bf2(x, dta, glm_structure, prior_scales_test)))
  
  #combine the scores...notice the change in the order from before
  scores_train = rbind(scores_train1, scores_train2)
  
  #set up data frame for bayesglm for testing set same as before
  dta = data.frame(index = rep(grid, N_test),
                   value = c(t(Curves_test)),
                   id = rep(1:N_test, each = D))
  
  npc = length(eigen_vals1)
  
  if(npc>1){
    for (z in 1:npc) {
      dta <- cbind(dta, rep(eigen_funcs1[,z], N_test))
    }
  }else{
    dta = cbind(dta, matrix(eigen_funcs1, ncol =1))
  }
  names(dta)[4:(4 + npc - 1)] <- c(paste0("psi", 1:npc))
  
  dta$mu = c(rep(mu_true0, N_test/2), rep(mu_true1, N_test/2))
  glm.structure = paste(paste0("psi", 1:npc), collapse = "+")
  glm_structure = paste("value ~ -1 + offset(mu) +" , glm.structure , sep="")
  
  #set true eigenvalues
  lambda_k = exp(-Ks/3)
  #prior_scales_test = sqrt(lambda_k)
  prior_scales_test = lambda_k
  
  #vec = matrix(1:(N_test/2), ncol = 1)
  vec = matrix(which(Classes_test==1), ncol = 1)
  scores_test1 = t(apply(vec, 1, function(x) regression_bf2(x, dta, glm_structure, prior_scales_test)))
  
  lambda_k = exp(-Ks/2)
  #prior_scales_test = sqrt(lambda_k)
  prior_scales_test = lambda_k
  
  #Get scores 
  vec = matrix(which(Classes_test==2), ncol = 1)
  scores_test2 = t(apply(vec, 1, function(x) regression_bf2(x, dta, glm_structure, prior_scales_test)))
  
  #again notice that scores were rearranged
  scores_test = rbind(scores_test1, scores_test2)
  
  #prior probability
  prior_g = c(table(Classes_train)/length(Classes_train))
  
  #data was rearranged so we have to change order of scores
  Classes_train2 = c(rep(1, sum(Classes_train==1)), rep(2, sum(Classes_train==2)))
  Classes_test2 = c(rep(1, sum(Classes_test==1)), rep(2, sum(Classes_test==2)))
  
  pred_classes = nb_updated_grid(scores_train, Classes_train2, prior_g, scores_test, min.h = 0.5, max.h = 1.5)
  
  et = Sys.time()
  time_diff = et - st
  
  misclass_mat[4, 1] = 1- sum(pred_classes==Classes_test2)/length(Classes_test2)
  sens_mat[4,1] = sum(pred_classes==1 & Classes_test2 == 1)/(sum(pred_classes==2 & Classes_test2 == 1)+sum(pred_classes==1 & Classes_test2 == 1))
  spec_mat[4,1] = sum(pred_classes==2 & Classes_test2 == 2)/(sum(pred_classes==2 & Classes_test2 == 2)+sum(pred_classes==1 & Classes_test2 == 2))
  precision_mat[4,1] = sum(pred_classes==1 & Classes_test2 == 1)/(sum(pred_classes= 1 & Classes_test2 == 1)+sum(pred_classes==1 & Classes_test2 == 2))
  time_mat[4, 1] = as.numeric(time_diff)
  
  
  
  ##
  #Oracle Latent_curves 1
  ##

  N=N2
  p=c(0.5,0.5)
  grid = seq(0,1, len=D)
  tt=seq(0,1, len=D)
  
  #estimate non-binary functions. Just like before generating 2*N curves (N for each group)
  Curves = generate_data(scenario = 2, grid = grid,  N=2*N,  p = p,
                         sigma = 0.01, binary = F)

  Y0 = t(Curves[1:(dim(Curves)[1]/2),])
  Y1 = t(Curves[-(1:(dim(Curves)[1]/2)),])

  #estimate non-binary functions, 2*N_test curves (N for each group)
  Curves = generate_data(scenario = 2, grid = grid,  N=2*N_test,  p = p,
                         sigma = 0.01, binary = F)
  Y0_test = t(Curves[1:(dim(Curves)[1]/2),])
  Y1_test = t(Curves[-(1:(dim(Curves)[1]/2)),])

  Curves_train = t(cbind(Y0, Y1))
  Curves_test = t(cbind(Y0_test, Y1_test))

  #get the Curves based on the randomly assigned classes
  Curves_train = Curves_train[(1:N)+N*(Classes_train-1),]
  Curves_test = Curves_test[(1:(N_test))+(N_test)*(Classes_test-1),]

  N_train = length(Classes_train)

  st = Sys.time()

  #get eigenvals and eigenfunctions based on estimates from Dai et al. 2017
  kb_e = fpca_estimate_functions_dai(Curves_train, Classes_train)

  #Function to approximate integral to get scores. 
  #Curves are dense enough for this approximation to provide accurate estimates
  get_scores = function(x){
    apply(kb_e$eigen_funcs, 2,  function(y) sum(y*x)/D)
  }

  #get Groups for training and testing sets
  scores_train2 = t(apply(t(t(Curves_train)-kb_e$mu_hat), 1, function(x) get_scores(x)))
  scores_test2 = t(apply(t(t(Curves_test)-kb_e$mu_hat), 1, function(x) get_scores(x)))

  #Fit Bayes classifier
  pred_classes = nb_updated_grid(scores_train2, Classes_train, prior_g, scores_test2, min.h = 0.5, max.h = 1.5)

  et = Sys.time()
  time_diff = et - st

  #Record results
  misclass_mat[5, 1] = 1- sum(pred_classes==Classes_test)/length(Classes_test)
  sens_mat[5,1] = sum(pred_classes==1 & Classes_test == 1)/(sum(pred_classes==2 & Classes_test == 1)+sum(pred_classes==1 & Classes_test == 1))
  spec_mat[5,1] = sum(pred_classes==2 & Classes_test == 2)/(sum(pred_classes==2 & Classes_test == 2)+sum(pred_classes==1 & Classes_test == 2))
  precision_mat[5,1] = sum(pred_classes==1 & Classes_test == 1)/(sum(pred_classes= 1 & Classes_test == 1)+sum(pred_classes==1 & Classes_test == 2))
  time_mat[5, 1] = as.numeric(time_diff)

  return(c(misclass_mat, sens_mat, spec_mat, precision_mat, time_mat))
  
}

#For loop over N values of interest
for(i2 in 1:length(Ns)){

  N2 = Ns[i2]
  
  #apply core function to each core to parallelize the functions
  results = mclapply(1:ITER, function(x) core_function(x), mc.cores = numCores_to_run)
  
  #uncomment if you want to see results 
  #print(results)
  
  #format results based on outputs
  results = matrix(unlist(results), ncol = ITER)
  #print(results)
  misclass_mat = results[1:5,]
  time_mat = results[21:25,]
  sens_mat = results[6:10,]
  spec_mat = results[11:15,]
  precision_mat = results[16:20,]
  
  #save results in list
  misclass_list[[i2]] = misclass_mat
  time_list[[i2]] = time_mat
  sens_list[[i2]] = sens_mat
  spec_list[[i2]] = spec_mat
  precision_list[[i2]] = precision_mat
  
  #display results at each N
  print("MisClassificaiton")
  print(rowMeans(misclass_mat, na.rm = T))
  print(apply(misclass_mat, 1, function(y) sd(y, na.rm=T)))
  print("Time")
  print(rowMeans(time_mat, na.rm = T))
  print(apply(time_mat, 1, function(y) sd(y, na.rm=T)))
  print("Sensitivity")
  print(rowMeans(sens_mat, na.rm = T))
  print(apply(sens_mat, 1, function(y) sd(y, na.rm=T)))
  print("Specificity")
  print(rowMeans(spec_mat, na.rm = T))
  print(apply(spec_mat, 1, function(y) sd(y, na.rm=T)))
  print("Precision")
  print(rowMeans(precision_mat, na.rm = T))
  print(apply(precision_mat, 1, function(y) sd(y, na.rm=T)))

}

#display results at end 

#Misclassification Results
#mean
lapply(misclass_list, function(x) rowMeans(x, na.rm=T))
#sd
lapply(misclass_list, function(x) apply(x, 1, function(y) sd(y, na.rm=T)))

#Timing Results
#mean
lapply(time_list, function(x) rowMeans(x, na.rm=T))
#sd
lapply(time_list, function(x) apply(x, 1, function(y) sd(y, na.rm=T)))

#Sensitivity
#mean
lapply(sens_list, function(x) rowMeans(x, na.rm=T))
#sd
lapply(sens_list, function(x) apply(x, 1, function(y) sd(y, na.rm=T)))

#Specificity
#mean
lapply(spec_list, function(x) rowMeans(x, na.rm=T))
#sd
lapply(spec_list, function(x) apply(x, 1, function(y) sd(y, na.rm=T)))

#Precision
#mean
lapply(precision_list, function(x) rowMeans(x, na.rm=T))
#sd
lapply(precision_list, function(x) apply(x, 1, function(y) sd(y, na.rm=T)))



####
# Format Results for latex tables
####

##
#Misclassification
##

m1 = lapply(misclass_list, function(x) rowMeans(x, na.rm=T))
s1 = lapply(misclass_list, function(x) apply(x, 1, function(y) sd(y, na.rm=T)))

xtable(t(matrix(paste(round(unlist(m1), 2), " (", round(unlist(s1), 2), ")", sep = ""), nrow = 5)))

##
#Times
##

m1 = lapply(time_list, function(x) rowMeans(x, na.rm=T))
s1 = lapply(time_list, function(x) apply(x, 1, function(y) sd(y, na.rm=T)))

xtable(t(matrix(paste(round(unlist(m1), 2), " (", round(unlist(s1), 2), ")", sep = ""), nrow = 5)))

##
#Sensitivity
##

m1 = lapply(sens_list, function(x) rowMeans(x, na.rm=T))
s1 = lapply(sens_list, function(x) apply(x, 1, function(y) sd(y, na.rm=T)))

xtable(t(matrix(paste(round(unlist(m1), 2), " (", round(unlist(s1), 2), ")", sep = ""), nrow = 5)))

##
# Specificity
##

m1 = lapply(spec_list, function(x) rowMeans(x, na.rm=T))
s1 = lapply(spec_list, function(x) apply(x, 1, function(y) sd(y, na.rm=T)))

xtable(t(matrix(paste(round(unlist(m1), 2), " (", round(unlist(s1), 2), ")", sep = ""), nrow = 5)))

##
#Precision results
##

m1 = lapply(precision_list, function(x) rowMeans(x, na.rm=T))
s1 = lapply(precision_list, function(x) apply(x, 1, function(y) sd(y, na.rm=T)))

xtable(t(matrix(paste(round(unlist(m1), 2), " (", round(unlist(s1), 2), ")", sep = ""), nrow = 5)))

