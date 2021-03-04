README.txt

#######
#Project: "Classification of Social Media Users Using a Generalized Functional Analysis"
#Authors: Anthony Weishampel
#Date: 03/04/2021
#######

The binary series for the data are provided and analyzed by the accompanying R code. This repository contains all of the code required to replicate the analyses in "Classification of Social Media Users Using a Generalized Functional Analysis". 

The follow R dependencies are required for compilation: 
	parallel, xtable, fda, refund, Matrix, MASS, arm, mgcv, randomforest, dtw


We did have access to cluster and ran much of the code across 16 different cores. We highly recommend access to multiple cores/cloud computing for compiling this code. The code was run using (7) DELL R7425 Dual Processor AMD Epyc 32 core 2.2 GHz machines with 512GB RAM each running 64Bit Ubuntu Linux Version 18.04. 
Instructions for Use

Data: 
	The formatted binary series are found in the data fold. In the provided binary time series, there are four csv files, one for each type of account. Each row in these csv files correspond to the binary series of one account. The binary series are summarized using five minute windows, where a 1 is recorded if the account posted during those five minutes and a 0 is presented otherwise. We provided code to convert the provided binary series to different window sizes by binning the data. 

	The Twitter data for the genuine and automated accounts is located at http://mib.projects.iit.cnr.it/dataset.html. Twitter data for the foreign state linked accounts are provided in the Data Archive of Twitter Civic Integrity https://about.twitter.com/en/our-priorities/civic-integrity. All of the data are available for academic research and provided as csv files. In these files the User ID field is the unique identifier for the accounts. The data contain information about the Tweets and the user profile information. 

Reproducibility:

	All of the tables for the simulation study (Tables 1, A1, A2) can be reproduced using the provided code. To reproduce these tables: 

		1) Compile the 01_functions.R code. 
		2) In each of the 02_Simulations_Scenario files, you can select the number of cores to run the code on by changing the numCores_to_run line in the code. 
		3) You also need to select the value for the number of observations per curves denoted by m in the paper and D in the code. This is accomplished by changing the value of D to be either 150, 250, or 500.
		4)  Additionally, for Scenario C the scenario3_data.Rdata needs to be download and it location needs to be included by replacing the “/dir/” to the appropriate directory. 
		5)  Compile the resulting 02_Simulation_Scenario files. Each of these files will produce the latex syntax for the three sets of three lines (N = 50, 300, 1000) for the corresponding lines (Scenario and m) of table 1, A1, A2. Thus in total each one of the three 02_Simulation_Scenario R code will be compiled three times, one per value of D. 

		The simulation code also produces lines for tables to record computation time and precision. 

	To reproduce the results for the data analyses, there are four different files: one for each grouping structure. You need to download the 4 .csv data files in the data folder and change the dir location in the code to the location of the files. The provided R files are set to 30 minute intervals and 14 days of analyses. To run the sensitivity analyses you have to change the values of J_for_analysis and num_min to the desired values. These files only include the proposed method and random forest classifier. You also need to change the /dir/ where the output will be saved.  The DeBot classifier is included in a separate file and similar changes are needed. 

	The accuracies of Table 2 and the confusion matrices of table 3 can be reproduced by:
		1) Downloading all four .csv data files and changing the /dir/ in each of the files to the appropriate directory.  
		2) Compile the 01_functions.R code.
		3) In each of the 03 data analysis files, you can select the number of cores to run the code on by changing the numCores_to_run line in the code. The 03 data analysis files provide the accuracies and the results of the confusion matrices at each iteration. 
		4) Make the appropriate /dir/ changes in the Debot file and compile the code. 
		5) Sensitivity and Specificity can be obtained by analyzing the saved output data from the confusion matrices. 
