##################
#Debot Analysis for "Classification of Social Media Users Using a Generalized Functional Analysis"
#
#Author: Anthony Weishampel
#Date: 2/27/2021
#
#
#Notes: Make sure to compile the 01_functions.R file first
##################

library(dtw)
library(parallel)

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



####
#Function to get the k closest neighbors guess based on dtw
#
#Inputs:
# k: number of neighbors
# x: current binary series in testing set
# binary_train: N x D binary series matrix of training set
# classes_train: vector of training set labels
#
# Outputs: K-NN classifier results
#####
dtw_KNN = function(k=5, x, binary_train, classes_train){
  
  dtws = dtwDist(binary_train, matrix(x, nrow= 1))
  k_val = sort(dtws)[k]
  dtws_under_k = which(dtws <= k_val)
  k_classes = classes_train[dtws_under_k]
  guess = names(which.max(table(k_classes)))
  return(guess)
  
}

####
#Function to get the k closest neighbors based on dtw
#
#Inputs:
# k: number of neighbors
# x: current binary series in testing set
# binary_train: N x D binary series matrix of training set
# classes_train: vector of training set labels
#
# Outputs: K-NN classifier all neighbors labels not just guess
#####
dtw_KNN_total = function(k=5, x, binary_train, classes_train){
  dtws = dtwDist(binary_train, matrix(x, nrow= 1))
  k_val = order(dtws)[1:k]
  k_classes = classes_train[k_val]
  return(k_classes)
}

####
#Function to run on multiple cores
#
#Inputs:
# i: for the set.seed on each core
#
# Outputs: results of DTW_KNN 
#####
core_function = function(i, k_val=5){
  
  set.seed(i)
  
  acc_mat = matrix(NA, nrow = 1,  ncol = 1)
  precision_mat = matrix(NA, nrow = 1,  ncol = 1)
  sens_mat = matrix(NA, nrow = 1,  ncol = 1)
  spec_mat = matrix(NA, nrow = 1,  ncol = 1)
  
  if(i%%25 == 0 || i==1){
    print(paste("On Iteration ", i , " out of ", ITER))
  }
  
  ###
  #DTW Method
  ###
  
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
  Curves_train = rbind(Curves_gtrain, Curves_btrain, Curves_itrain, Curves_rtrain)
  Curves_test = rbind(Curves_gtest, Curves_btest, Curves_itest, Curves_rtest)
  D = dim(Curves_train)[2]
  
  #Two Levels Levels Can easily change to 4 levels by changing the levels
  Classes_train = c(rep(1, N_gtrain), rep(2, N_btrain), rep(3, N_itrain), rep(4, N_rtrain))
  Classes_train = as.factor(Classes_train)
  Classes_test = c(rep(1, N_genuine-N_gtrain), rep(2, N_bots-N_btrain), rep(3, N_i-N_itrain), rep(4, N_r-N_rtrain))
  Classes_test = as.factor(Classes_test)
  
  individuals = 1:length(Classes_train)
  N_train = length(Classes_train)
  N_test = length(Classes_test)
  
  #set number of accounts in training set 
  #Needs to be smaller because of very high computation time. 
  N_train = 100
  training_accounts = sample(1:length(Classes_train), N_train)
  Classes_train = Classes_train[training_accounts]
  Curves_train = Curves_train[training_accounts,]
  
  pred_classes = apply(Curves_test, 1, function(x) dtw_KNN(k=k_val, x, Curves_train, Classes_train))
  #t1 = table(guess, Classes_test)
  #acc_mat[1,z_counter] = sum(diag(t1))/sum(t1)

  acc_mat[1, 1] = sum(pred_classes==Classes_test)/length(Classes_test)
  sens_mat[1,1] = sum(pred_classes==1 & Classes_test == 1)/(sum(pred_classes==2 & Classes_test == 1)+sum(pred_classes==1 & Classes_test == 1))
  spec_mat[1,1] = sum(pred_classes==2 & Classes_test == 2)/(sum(pred_classes==2 & Classes_test == 2)+sum(pred_classes==1 & Classes_test == 2))
  precision_mat[1,1] = sum(pred_classes==1 & Classes_test == 1)/(sum(pred_classes= 1 & Classes_test == 1)+sum(pred_classes==1 & Classes_test == 2))

  print(c(acc_mat, sens_mat, spec_mat, precision_mat))
  pred_classes = factor(pred_classes, levels = levels(as.factor(Classes_test)))
  t1 = table(pred_classes,Classes_test)
  print(t1)
  
  return_list = list()
  return_list[[1]] = c(acc_mat, sens_mat, spec_mat, precision_mat)
  return_list[[2]] = t1
  
  return(return_list)
}



####
#Function to get the best values of k based on grid search
#
#Inputs:
# k_min: lowest possible value of k
# k_max: largest possible value of k
#
# Outputs: k value used in classifier
#####
get_k = function(k_min = 5, k_max = 20){
  
  set.seed(12345)
  
  ###
  #DTW Method
  ###
  
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
  Curves_train = rbind(Curves_gtrain, Curves_btrain, Curves_itrain, Curves_rtrain)
  Curves_test = rbind(Curves_gtest, Curves_btest, Curves_itest, Curves_rtest)
  D = dim(Curves_train)[2]
  
  #Two Levels Levels Can easily change to 4 levels by changing the levels
  Classes_train = c(rep(1, N_gtrain), rep(2, N_btrain), rep(3, N_itrain), rep(4, N_rtrain))
  Classes_train = as.factor(Classes_train)
  Classes_test = c(rep(1, N_genuine-N_gtrain), rep(2, N_bots-N_btrain), rep(3, N_i-N_itrain), rep(4, N_r-N_rtrain))
  Classes_test = as.factor(Classes_test)
  
  individuals = 1:length(Classes_train)
  N_train = length(Classes_train)
  N_test = length(Classes_test)

  
  #set number of accounts in training set 
  #Needs to be smaller because of large computation time. 
  N_train = 100
  training_accounts = sample(1:length(Classes_train), N_train)
  Classes_train = Classes_train[training_accounts]
  Curves_train = Curves_train[training_accounts,]

  
  #smaller test size to find the kl
  N_test = 100
  training_accounts = sample(1:length(Classes_test), N_test)
  Classes_test = Classes_test[training_accounts]
  Curves_test = Curves_test[training_accounts,]
  
  #dtw_KNN_first = function(Curves_test, Curves_train, Classes_train, Classes_test, k_min=3, k_max=20)
  ks = k_min:k_max
  accs = rep(NA, length(ks))
  
  pred_classes = apply(Curves_test, 1, function(x) dtw_KNN_total(k=k_max, x, Curves_train, Classes_train))
  pred_classes_mat = matrix(unlist(pred_classes), ncol = length(Classes_test))
  
  for(j in 1:length(ks)){
    
    k = ks[j]
    guesses = apply(pred_classes_mat[1:k,], 2, function(x) names(which.max(table(x))))
    accs[j] = sum(guesses==Classes_test)/length(Classes_test)
    
  }
  
  #get the k
  k=ks[which.max(accs)]
  
  #print(k)
  return(k)

}



####
# Get Initial K value
####
k_val = get_k()
#print(k_val)

#apply function for each cluster
results = mclapply(1:ITER, function(x) core_function(x, k_val), mc.cores = numCores_to_run)
print(results)

results1 = list()
results2 = list()

for(i in 1:ITER){
  results1[[i]] = results[[i]][[1]]
  results2[[i]] = results[[i]][[2]]
}

results = matrix(unlist(results1), ncol = ITER)
#print(results)
results2 = matrix(unlist(results2), ncol = ITER)
#print(results2)

##
#Get mean results for accuracy, sensitivity, specificity, and precision
##
print(apply(results, 1, function(x) mean(x, na.rm = T)))

##
#Get sd results for accuracy, sensitivity, specificity, and precision
##
print(apply(results, 1, function(x) sd(x, na.rm = T)))

##
#Get average numbers in confusion matrix 
##
print(apply(results2, 1, function(x) mean(x, na.rm = T)))

##
#Get sd numbers in confusion matrix 
##
print(apply(results2, 1, function(x) sd(x, na.rm = T)))



