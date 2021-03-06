##################
#Data Analysis for "Classification of Social Media Users Using a Generalized Functional Analysis"
#
#Author:
#Date: 2/27/2021
#
#
#Notes: Make sure to compile the 01_functions.R file first
##################

library(readr)
library(randomForest)

#load each dataset...make sure to change the directory to the location of the data
timelines_genuine <- read_csv("/dir/timelines_genuine.csv")
timelines_bots <- read_csv("/dir/timelines_bots.csv")
timelines_fs1<- read_csv("/dir/timelines_fs1.csv")
timelines_fs2 <- read_csv("/dir/timelines_fs2.csv")

#Number of days to analyze
J_for_analysis = 14

#number of minutes per break (multiple of 5)
num_mins = 30

#random_days
#variable to select random starting day (used in days sensitivity analysis)
Random_start_day = F

#Number of cores able to run code on
numCores_to_run = 16

#Get Number for iterations
ITER = 50

#Format all data 
# this is not the number of days wanted to analyze
J = 28

#Set D to the end of the 28 days
# number of days * number of 5 minute long intervals in a day
D = J*24*12+1

#Format data for each group
timelines_genuine = timelines_genuine[,2:D]
timelines_bots = timelines_bots[,2:D]
timelines_fs1 = timelines_fs1[,1:(D-1)]
timelines_fs2 = timelines_fs2[,1:(D-1)]

#format to matrix
timelines_genuine = data.matrix(timelines_genuine)
timelines_bots = data.matrix(timelines_bots)
timelines_fs1 = data.matrix(timelines_fs1)
timelines_fs2 = data.matrix(timelines_fs2)


#Function to bin binary series 
#nrow = 2 is 10 minutes 
#nrow = 3 is 15 minutes
#nrow = 4 is 20 minutes
#nrow = 6 is 30 minutes
#nrow = 9 is 45 minutes
#nrow = 12 is 60 minutes
#nrow = 24 is 120 minutes
#nrow = 72 is 360 minutes
format_into_x_min = function(timeline){
  
  sum_vec = t(matrix(timeline, nrow = num_mins/5))
  return_vec = rowSums(sum_vec)
  return_vec = ifelse(return_vec==0, 0, 1)
  return(return_vec)
  
}

#format timelines to appropriate number of minutes per interval
timelines_genuine = t(apply(timelines_genuine, 1, function(x) format_into_x_min(x)))
timelines_bots = t(apply(timelines_bots, 1, function(x) format_into_x_min(x)))
timelines_fs1 = t(apply(timelines_fs1, 1, function(x) format_into_x_min(x)))
timelines_fs2 = t(apply(timelines_fs2, 1, function(x) format_into_x_min(x)))

set.seed(1234)

core_function = function(i){
  
  #set seed for each iteration needed because code is run on multiple cores
  set.seed(i)
  
  #matrix to record accuracies
  acc_mat = matrix(NA, nrow = 2,  ncol = 1)
  
  if(i%%25 == 0 || i==1){
    print(paste("On Iteration ", i , " out of ", ITER))
  }
  
  #only select the first 14 days
  J = J_for_analysis
  #J times number of observations in a day
  D = J*(60*24)/num_mins
  nd = dim(timelines_genuine)[2]-D
  
  if(Random_start_day){
    start_time = sample((1:nd),1)
  }else{
    start_time = 1
  }
  
  timelines_genuine2 = timelines_genuine[,start_time:(start_time+D-1)]
  timelines_bots2 = timelines_bots[,start_time:(start_time+D-1)]
  timelines_fs1_2 = timelines_fs1[,start_time:(start_time+D-1)]
  timelines_fs2_2 = timelines_fs2[,start_time:(start_time+D-1)]
  
  #only include accounts which were active during those J days
  # can remove if you want all accounts but run into issues with GLM and GAM. 
  min_num_tweets = 1
  
  timelines_genuine2 = timelines_genuine2[rowSums(timelines_genuine2)>min_num_tweets,]
  timelines_bots2 = timelines_bots2[rowSums(timelines_bots2)>min_num_tweets,]
  timelines_fs1_2 = timelines_fs1_2[rowSums(timelines_fs1_2)>min_num_tweets,]
  timelines_fs2_2 = timelines_fs2_2[rowSums(timelines_fs2_2)>min_num_tweets,]
  
  #select training and testing sets for genuine accounts
  N_genuine = dim(timelines_genuine2)[1]
  p_train = 0.8
  N_gtrain = round(N_genuine * p_train)
  train_genuine = sample(1:N_genuine, N_gtrain)
  Curves_gtrain = timelines_genuine2[train_genuine, ]
  Curves_gtest = timelines_genuine2[-train_genuine, ]
  
  #select training and testing sets for bots
  N_bots = dim(timelines_bots2)[1]
  N_btrain = round(N_bots * p_train)
  train_bots = sample(1:N_bots, N_btrain)
  Curves_btrain = timelines_bots2[train_bots, ]
  Curves_btest = timelines_bots2[-train_bots, ]
  
  #select training and testing sets for fs2
  N_r = dim(timelines_fs2_2)[1]
  N_rtrain = round(N_r * p_train)
  train_r = sample(1:N_r, N_rtrain)
  Curves_rtrain = timelines_fs2_2[train_r, ]
  Curves_rtest = timelines_fs2_2[-train_r, ]
  
  #select training and testing sets for fs1
  N_i = dim(timelines_fs1_2)[1]
  N_itrain = round(N_i * p_train)
  train_i = sample(1:N_i, N_itrain)
  Curves_itrain = timelines_fs1_2[train_i, ]
  Curves_itest = timelines_fs1_2[-train_i, ]
  
  #Group all data together into one training set matrix and one testing set matrix
  Curves_train = rbind(Curves_gtrain, Curves_btrain)
  Curves_test = rbind(Curves_gtest, Curves_btest)
  D = dim(Curves_train)[2]

  #Two Levels Levels Can easily change to 4 levels by changing the levels
  Classes_train = c(rep(1, N_gtrain), rep(2, N_btrain))
  Classes_train = as.factor(Classes_train)
  Classes_test = c(rep(1, N_genuine-N_gtrain), rep(2, N_bots-N_btrain))
  Classes_test = as.factor(Classes_test)

  individuals = 1:length(Classes_train)
  N_train = length(Classes_train)
  N_test = length(Classes_test)
  
  tt=seq(0,1, len=D)
  
  #Step 1 of proposed method
  vec = matrix(1:(N_train), ncol = 1)
  smoothed_x = logit(t(apply(vec, 1, function(x) regression_g(x, Curves_train, tt, k=10))))
  
  #Step 2 of proposed method
  fpca.cur2 = fpca.face(smoothed_x, pve = 0.99, p=3, m=2, knots = 6) #lambda selected via grid search optim, #p=degree of splines
  #correct fpca.cur2 eigenval values b/c too large and functions aren't normalized
  #get multiplier is function of D
  get_multiplier = 1/D
  fpca.cur = fpca.cur2
  fpca.cur$efunctions = fpca.cur2$efunctions/sqrt(get_multiplier)
  fpca.cur$evalues = fpca.cur2$evalues*get_multiplier
  fpca.cur$scores = fpca.cur2$scores*sqrt(get_multiplier)
  
  #Step 3 of proposed method to estimate scores
  fit = list(mu = fpca.cur$mu,
             evalues = fpca.cur$evalues,
             efunctions = fpca.cur$efunctions)
  
  mu_t_hat = fit$mu
  eigen_vals1 = fit$evalues
  eigen_funcs1 = fit$efunctions
  
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
  names(dta)[4:(4 + npc - 1)] <- c(paste0("psi", 1:npc))
  dta$mu = rep(mu_t_hat , N_train)
  
  glm_structure = paste(paste0("psi", 1:npc), collapse = "+")
  glm_structure = paste("value ~ -1 + offset(mu) +" , glm.structure , sep="")
  #set scale for the glm 
  prior_scales_test = eigen_vals1
  
  vec = matrix(1:N_train, ncol = 1)
  scores_train = t(apply(vec, 1, function(x) regression_bf2(x, dta, glm_structure, prior_scales_test)))
  
  ###
  #Get scores for testing set
  ###
  
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
  
  #step 4 of proposed method
  # train the bayes classifier
  prior_g = c(table(Classes_train)/length(Classes_train))
  guess = nb_updated_grid(scores_train, Classes_train, prior_g, scores_test, CV = 5)
  
  t1 = table(guess, Classes_test)
  
  pred_classes = factor(guess, levels = levels(as.factor(Classes_test)))
  t1 = table(pred_classes, Classes_test)
  print(t1)
  acc_mat[1,1] = sum(diag(t1))/sum(t1)
  t3 = t1
  
  
  #Run random forest classifier
  rf1 = randomForest(y = Classes_train, x = Curves_train, ntree = 100)
  guess = apply(Curves_test, 1, function(x) predict(rf1, x))
  t2 = table(guess, Classes_test)
  acc_mat[2,1] = sum(diag(t2))/sum(t2)
  
  # return list
  return_list = list()
  return_list[[1]] = c(acc_mat) #accuracy
  return_list[[2]] = t3 #confusion matrix for proposed method
  return_list[[3]] = t2 #confusion matrix for random forest
	
  return(return_list)
  
  
  
}

#apply core function to each core to parallelize the analysis
results = mclapply(1:ITER, function(x) core_function(x), mc.cores = numCores_to_run)

#get results from the previous apply function
results1 = list()
results2 = list()

for(i in 1:ITER){
  results1[[i]] = results[[i]][[1]]
  results2[[i]] = results[[i]][[2]]
}

#format results & print them
results = matrix(unlist(results1), ncol = ITER)
print(results)

results2 = matrix(unlist(results2), ncol = ITER)
print(results2)

#Print mean and sd
print(apply(results, 1, function(x) mean(x, na.rm = T)))
print(apply(results, 1, function(x) sd(x, na.rm = T)))

#Print mean and sd
print(apply(results2, 1, function(x) mean(x, na.rm = T)))
print(apply(results2, 1, function(x) sd(x, na.rm = T)))

#remove data before saving RData
remove(timelines_genuine, timelines_bots, 
         timelines_fs1, timelines_fs2)

#save results
#save.image("/dir/bot_genuine.RData")
