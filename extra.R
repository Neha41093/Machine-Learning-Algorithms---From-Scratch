#Extra Question -- Part 1 -- (Initialization using K-means++ algorithm)

rm(list=ls())
#set according to your working directory
old.dir <- getwd()
setwd("C:/Users/Neha Rawat/Desktop/IU-Data Science/AML/Assignment 2")

#Ionosphere dataset
ionosphere_data_with_labels <- read.csv("ionosphere_data.csv",header=FALSE)
ionosphere_data <- as.matrix(ionosphere_data_with_labels[,1:34])
#getting true centers
good_mean_indices <- which(ionosphere_data_with_labels[,35] == "g")
bad_mean_indices <- which(ionosphere_data_with_labels[,35] == "b")
true_mean_good_data <- as.matrix(ionosphere_data[good_mean_indices,])
true_mean_bad_data <- as.matrix(ionosphere_data[bad_mean_indices,])
true_mean_good <- apply(true_mean_good_data,2,mean)
true_mean_bad <- apply(true_mean_bad_data,2,mean)

#****************************EM using normal initialization method*************
library(mvtnorm)

#Initializations
k_values1 <- c(2:5)
max_iter1 <- 50
convergence_check1 <- 0.0000005
max_seed_iter1 <- 20
sum_error2 <- 0
gb_errors2 <- matrix(NA, nrow=max(k_values1)-1, ncol=max_seed_iter1)
mark <- 0
s_prob_m <- matrix(NA, nrow=max(k_values1)-1, ncol=max_seed_iter1)


#Expectation (E-Step) Functions
BayesNumerator <- function(all.matrix,cl.mean,cl.covar,cl.prior) {
  b_num <- as.matrix(outer(cl.prior,as.matrix(dmvnorm(all.matrix,cl.mean,cl.covar)))) 
  return(b_num)
}

BayesWeightComplete <- function(k.value,all.matrix,k.mean,k.covar,k.prior){
  k_vec <- as.matrix(c(1:k.value))
  by_wt <- t(as.matrix(apply(k_vec, 1 ,function(x) BayesNumerator(all.matrix,k.mean[x,],k.covar[[x]],k.prior[x,]))))
  by_prop <- as.matrix(apply(by_wt,2,function(x) x/sum(x)))
  return(by_prop)
}

#Maximization (M-Step) Functions
ParametersUpdate <- function(e.matrix,all.matrix){
  mean.total <- matrix(data=0,nrow=nrow(e.matrix),ncol=ncol(all.matrix))
  cov.total <- list()
  #mean part
  mean.num <- as.matrix(e.matrix %*% all.matrix)
  for(i in c(1:nrow(e.matrix))){
    if (sum(e.matrix[i,]) == 0){ p <- 1
    }else
    {p <-sum(e.matrix[i,]) 
    }
    mean.total[i,] <- mean.num[i,]/p
  }
  #priors part
  prior.total <- matrix(apply(e.matrix,1,mean),nrow=nrow(e.matrix))
  #covariance part
  for(i in c(1:nrow(e.matrix))){
    singular <- 0
    if (sum(e.matrix[i,]) == 0){ p <- 1
    }else
    {p <-sum(e.matrix[i,]) 
    }
    cov.total1 <- matrix(0,nrow=ncol(all.matrix),ncol=ncol(all.matrix))
    cov.diff <- as.matrix(t(t(all.matrix)-mean.total[i,]))
    for(j in c(1:nrow(all.matrix))){
      cov.step1 <- matrix(outer(e.matrix[i,j],tcrossprod(cov.diff[j,])),ncol=ncol(all.matrix),nrow=ncol(all.matrix))
      cov.total1 <- cov.total1 + cov.step1
    }
    cov.total[[i]]<- matrix(cov.total1/p,ncol=ncol(all.matrix),nrow=ncol(all.matrix))
    
    #dealing with singular covariance matrices by reseeding
    if (det(cov.total[[i]]) == 0){
      cov.total[[i]] <- diag(ncol(all.matrix))
      mean.total[i,] <- all.matrix[sample(nrow(all.matrix),1,replace=FALSE), ]
      prior.total[i] <- 1/nrow(e.matrix)
    }
  }
  #return all 3 for all clusters
  parameters <- list(meanupd=mean.total,covupd=cov.total,priorupd=prior.total)  
  return(parameters)
}

#Helper function to calculate Euclidean distance between a vector and matrix
EuclideanDistanceMini <- function(one.vec, all.matrix) {
  diff.matrix <- t(t(all.matrix) - one.vec)
  this.row <- apply(diff.matrix, 1, function(x) sum(x**2)**0.5)
  return(this.row)
}

#Function to calculate Euclidean distance between two full matrices
EuclideanDistanceMain <- function(matrix_one, matrix_two) {
  dist.matrix <- apply(matrix_one, 1, function(x) EuclideanDistanceMini(x, matrix_two))
  dist.matrix <- matrix(dist.matrix, nrow=nrow(matrix_one), byrow=TRUE)
  return(dist.matrix)
}


#Function to calculate the error proportion i.e bad/(good+bad)
GoodBadError <- function(k,assigned_clusters,cluster.centers,all.matrix) {
  good_bad_ratio <- list()
  if (k==2){
    for (k.value in c(1:k)) {
      good_bad_ratio[[k.value]] <- length(which(as.matrix(all.matrix[which(assigned_clusters == k.value),35]) == "b"))/length(which(assigned_clusters == k.value))
    }
  }else {
    dist_mat_good <- EuclideanDistanceMini(true_mean_good,cluster.centers)
    dist_mat_bad <- EuclideanDistanceMini(true_mean_bad,cluster.centers)
    dist_mat <- rbind(dist_mat_good,dist_mat_bad)
    dist_assign <- apply(dist_mat,2,which.min)
    assigned_clust_temp <- assigned_clusters
    for(i in 1:length(assigned_clust_temp)){
      assigned_clust_temp[i] <- dist_assign[assigned_clust_temp[i]]
    }
    for (k.value in c(1:2)) {
      if (length(which(assigned_clust_temp == k.value)) != 0){
        good_bad_ratio[[k.value]] <- length(which(as.matrix(all.matrix[which(assigned_clust_temp == k.value),35]) == "b"))/length(which(assigned_clust_temp == k.value))
      } else {
        good_bad_ratio[[k.value]] <- 0
      }
    }
  } 
  return(good_bad_ratio)
}

start.time <- Sys.time()
#Start main code for the algo
for (k1 in k_values1) {
  print(c("K-Value: ", k1))

  for (seed_iter1 in 1:max_seed_iter1) {
    print(c("Initialization Number: ", seed_iter1))
    
    #randomly initializing means
    cluster_means_initial <- as.matrix(ionosphere_data[sample(nrow(ionosphere_data),k1,replace=FALSE), ])
    cluster_cov_initial <- rep(list(diag(ncol(ionosphere_data))),k1)
    cluster_prior_initial <- matrix(rep(1/k1),nrow=k1,ncol=1)
    
    converged1 = FALSE
    
    while (mark != 1)
    {
      for (iter1 in 1:max_iter1) {
        
        #Running E-Step
        E_step_mat <- BayesWeightComplete(k1,ionosphere_data,cluster_means_initial,cluster_cov_initial,cluster_prior_initial)
        
        #Running M-Step
        M_step_list <- ParametersUpdate(E_step_mat,ionosphere_data)
        cluster_means_upd <- M_step_list$meanupd
        cluster_cov_upd <- M_step_list$covupd
        cluster_prior_upd <- M_step_list$priorupd
        
        #Calculate error/movement of means
        error1 <- EuclideanDistanceMain(cluster_means_initial, cluster_means_upd)
        for(f1 in c(1:k1)) {
          sum_error2 <- sum_error2 + error1[f1,f1]}
        sum_error21 <- sum_error2
        sum_error2 <- 0
        
        #Check convergence criteria
        if(sum_error21 <= convergence_check1 ) {
          converged1 <- TRUE
          break
        }
        else {
          cluster_means_initial <- cluster_means_upd
          cluster_cov_initial <- cluster_cov_upd
          cluster_prior_initial <- cluster_prior_upd
        }
      }
      #getting cluster assignments
      E_step_mat <- BayesWeightComplete(k1,ionosphere_data,cluster_means_upd,cluster_cov_upd,cluster_prior_upd) 
      cluster_assignment1 <- apply(E_step_mat,2,which.max)
      #dealing with empty clusters - reseeding
      k_list <- as.matrix(c(1:k1))
      cluster_count <- as.matrix(apply(k_list,1,function(x) length(which(cluster_assignment1==x))))
      #Tagging empty clusters
      check1 <- which(cluster_count == 0) 
      check1 <- as.matrix(check1)
      if (nrow(check1) == 0){
        print("No empty clusters")
        mark <- 1
        break
      } else {
        print("Empty cluster present")
        cluster_means_initial <- cluster_means_upd
        cluster_cov_initial <- cluster_cov_upd
        cluster_prior_initial <- cluster_prior_upd
        for (b in 1:nrow(check1)){
          cluster_means_initial[check1[b,],] <- ionosphere_data[sample(nrow(ionosphere_data),1,replace=FALSE), ]
          cluster_cov_initial[[check1[b,]]] <- diag(ncol(ionosphere_data))
          cluster_prior_initial[check1[b,],] <- 1/k1
        }
      }
    }
    
    mark <- 0

    gb_errors2[k1-1, seed_iter1] <- sum(do.call(rbind,GoodBadError(k1,cluster_assignment1,cluster_means_upd,ionosphere_data_with_labels)))
    
    #sse calculation (maximizing log-likelihood equivalents instead of minimizing euclidean distance)
    cluster_id2 <- as.matrix(c(1:k1))
    s_prob_mat <- t(as.matrix(apply(cluster_id2,1,function(x) BayesNumerator(ionosphere_data,cluster_means_upd[x,],cluster_cov_upd[[x]],cluster_prior_upd[x,]))))
    s_prob <- as.matrix(apply(as.matrix(apply(s_prob_mat,2,sum)),2,function(x) log(x)))
    s_prob_m[k1-1,seed_iter1] <- sum((s_prob))
  }
  
}
end.time <- Sys.time()
print(end.time - start.time)

#Plotting bad("b" label) proportions for k and iterations
plot(2:5, apply(gb_errors2,1,mean),
     type="b", pch = 19, frame = FALSE, 
     main = "Ionosphere data",
     xlab="Number of clusters K",
     ylab="Total good-bad errors")

# log-likelihood plot
plot(2:5, apply(s_prob_m,1,sum),
     type="b", pch = 19, frame = FALSE, 
     main = "Ionosphere data",
     xlab="Number of clusters K",
     ylab="Total sum of log-likelihoods")

#***********************EM using kmeans++ initialization method**********************
#Initializations
k_values2 <- c(2:5)
max_iter2 <- 50
convergence_check2 <- 0.0000005
max_seed_iter2 <- 20
sum_error3 <- 0
gb_errors3 <- matrix(NA, nrow=max(k_values1)-1, ncol=max_seed_iter1)
mark1 <- 0
s_prob_m1 <- matrix(NA, nrow=max(k_values1)-1, ncol=max_seed_iter1)


#Expectation (E-Step) Functions
BayesNumerator <- function(all.matrix,cl.mean,cl.covar,cl.prior) {
  b_num <- as.matrix(outer(cl.prior,as.matrix(dmvnorm(all.matrix,cl.mean,cl.covar)))) 
  return(b_num)
}

BayesWeightComplete <- function(k.value,all.matrix,k.mean,k.covar,k.prior){
  k_vec <- as.matrix(c(1:k.value))
  by_wt <- t(as.matrix(apply(k_vec, 1 ,function(x) BayesNumerator(all.matrix,k.mean[x,],k.covar[[x]],k.prior[x,]))))
  by_prop <- as.matrix(apply(by_wt,2,function(x) x/sum(x)))
  return(by_prop)
}

#Maximization (M-Step) Functions
ParametersUpdate <- function(e.matrix,all.matrix){
  mean.total <- matrix(data=0,nrow=nrow(e.matrix),ncol=ncol(all.matrix))
  cov.total <- list()
  #mean part
  mean.num <- as.matrix(e.matrix %*% all.matrix)
  for(i in c(1:nrow(e.matrix))){
    if (sum(e.matrix[i,]) == 0){ p <- 1
    }else
    {p <-sum(e.matrix[i,]) 
    }
    mean.total[i,] <- mean.num[i,]/p
  }
  #priors part
  prior.total <- matrix(apply(e.matrix,1,mean),nrow=nrow(e.matrix))
  #covariance part
  for(i in c(1:nrow(e.matrix))){
    singular <- 0
    if (sum(e.matrix[i,]) == 0){ p <- 1
    }else
    {p <-sum(e.matrix[i,]) 
    }
    cov.total1 <- matrix(0,nrow=ncol(all.matrix),ncol=ncol(all.matrix))
    cov.diff <- as.matrix(t(t(all.matrix)-mean.total[i,]))
    for(j in c(1:nrow(all.matrix))){
      cov.step1 <- matrix(outer(e.matrix[i,j],tcrossprod(cov.diff[j,])),ncol=ncol(all.matrix),nrow=ncol(all.matrix))
      cov.total1 <- cov.total1 + cov.step1
    }
    cov.total[[i]]<- matrix(cov.total1/p,ncol=ncol(all.matrix),nrow=ncol(all.matrix))
    
    #dealing with singular covariance matrices by reseeding
    if (det(cov.total[[i]]) == 0){
      cov.total[[i]] <- diag(ncol(all.matrix))
      mean.total[i,] <- all.matrix[sample(nrow(all.matrix),1,replace=FALSE), ]
      prior.total[i] <- 1/nrow(e.matrix)
    }
  }
  #return all 3 for all clusters
  parameters <- list(meanupd=mean.total,covupd=cov.total,priorupd=prior.total)  
  return(parameters)
}

#Helper function to calculate Euclidean distance between a vector and matrix
EuclideanDistanceMini <- function(one.vec, all.matrix) {
  diff.matrix <- t(t(all.matrix) - one.vec)
  this.row <- apply(diff.matrix, 1, function(x) sum(x**2)**0.5)
  return(this.row)
}

#Function to calculate Euclidean distance between two full matrices
EuclideanDistanceMain <- function(matrix_one, matrix_two) {
  dist.matrix <- apply(matrix_one, 1, function(x) EuclideanDistanceMini(x, matrix_two))
  dist.matrix <- matrix(dist.matrix, nrow=nrow(matrix_one), byrow=TRUE)
  return(dist.matrix)
}


#Function to calculate the error proportion i.e bad/(good+bad)
GoodBadError <- function(k,assigned_clusters,cluster.centers,all.matrix) {
  good_bad_ratio <- list()
  if (k==2){
    for (k.value in c(1:k)) {
      good_bad_ratio[[k.value]] <- length(which(as.matrix(all.matrix[which(assigned_clusters == k.value),35]) == "b"))/length(which(assigned_clusters == k.value))
    }
  }else {
    dist_mat_good <- EuclideanDistanceMini(true_mean_good,cluster.centers)
    dist_mat_bad <- EuclideanDistanceMini(true_mean_bad,cluster.centers)
    dist_mat <- rbind(dist_mat_good,dist_mat_bad)
    dist_assign <- apply(dist_mat,2,which.min)
    assigned_clust_temp <- assigned_clusters
    for(i in 1:length(assigned_clust_temp)){
      assigned_clust_temp[i] <- dist_assign[assigned_clust_temp[i]]
    }
    for (k.value in c(1:2)) {
      if (length(which(assigned_clust_temp == k.value)) != 0){
        good_bad_ratio[[k.value]] <- length(which(as.matrix(all.matrix[which(assigned_clust_temp == k.value),35]) == "b"))/length(which(assigned_clust_temp == k.value))
      } else {
        good_bad_ratio[[k.value]] <- 0
      }
    }
  } 
  return(good_bad_ratio)
}

start.time <- Sys.time()
#Start main code for the algo
for (k2 in k_values2) {
  print(c("K-Value: ", k2))
 
  for (seed_iter2 in 1:max_seed_iter2) {
    print(c("Initialization Number: ", seed_iter2))

    #randomly initializing means using ++ method
    cluster_mean_others <- matrix(NA,ncol=ncol(ionosphere_data),nrow=k2-1)
    cluster_mean_index1 <- sample(nrow(ionosphere_data),1,replace=FALSE)
    cluster_means_initial1 <- ionosphere_data[cluster_mean_index1, ]
    dist_cent_vec <- t(EuclideanDistanceMini(cluster_means_initial1,ionosphere_data[-cluster_mean_index1,]))
    dist_cent_prob <- t(apply(dist_cent_vec,1,function(x) x/sum(x)))
    cluster_sort <- sort(dist_cent_prob,decreasing=TRUE)
    for (h in 1:(k2-1)){
     cluster_mean_others[h,] <- ionosphere_data[min(which(dist_cent_prob == cluster_sort[h])), ]
    }
    cluster_means_initial2 <- rbind(cluster_means_initial1,cluster_mean_others)

    cluster_cov_initial2 <- rep(list(diag(ncol(ionosphere_data))),k2)
    cluster_prior_initial2 <- matrix(rep(1/k2),nrow=k1,ncol=1)
    
    converged2 = FALSE
    
    while (mark1 != 1)
    {
      for (iter2 in 1:max_iter2) {
        
        #Running E-Step
        E_step_mat1 <- BayesWeightComplete(k2,ionosphere_data,cluster_means_initial2,cluster_cov_initial2,cluster_prior_initial2)
        
        #Running M-Step
        M_step_list1 <- ParametersUpdate(E_step_mat1,ionosphere_data)
        cluster_means_upd2 <- M_step_list1$meanupd
        cluster_cov_upd2 <- M_step_list1$covupd
        cluster_prior_upd2 <- M_step_list1$priorupd
        
        #Calculate error/movement of means
        error2 <- EuclideanDistanceMain(cluster_means_initial2, cluster_means_upd2)
        for(f2 in c(1:k2)) {
          sum_error3 <- sum_error3 + error2[f2,f2]}
        sum_error22 <- sum_error3
        sum_error3 <- 0
        
        #Check convergence criteria
        if(sum_error22 <= convergence_check2 ) {
          converged2 <- TRUE
          break
        }
        else {
          cluster_means_initial2 <- cluster_means_upd2
          cluster_cov_initial2 <- cluster_cov_upd2
          cluster_prior_initial2 <- cluster_prior_upd2
        }
      }
      #getting cluster assignments
      E_step_mat1 <- BayesWeightComplete(k2,ionosphere_data,cluster_means_upd2,cluster_cov_upd2,cluster_prior_upd2) 
      cluster_assignment2 <- apply(E_step_mat1,2,which.max)
      #dealing with empty clusters - reseeding
      k_list1 <- as.matrix(c(1:k2))
      cluster_count1 <- as.matrix(apply(k_list1,1,function(x) length(which(cluster_assignment2==x))))
      #Tagging empty clusters
      check2 <- which(cluster_count1 == 0) 
      check2 <- as.matrix(check2)
      if (nrow(check2) == 0){
        print("No empty clusters")
        mark1 <- 1
        break
      } else {
        print("Empty cluster present")
        cluster_means_initial2 <- cluster_means_upd2
        cluster_cov_initial2 <- cluster_cov_upd2
        cluster_prior_initial2 <- cluster_prior_upd2
        for (b1 in 1:nrow(check2)){
          cluster_means_initial2[check2[b1,],] <- ionosphere_data[sample(nrow(ionosphere_data),1,replace=FALSE), ]
          cluster_cov_initial2[[check2[b1,]]] <- diag(ncol(ionosphere_data))
          cluster_prior_initial2[check2[b1,],] <- 1/k2
        }
      }
    }
    
    mark1 <- 0

    gb_errors3[k2-1, seed_iter2] <- sum(do.call(rbind,GoodBadError(k2,cluster_assignment2,cluster_means_upd2,ionosphere_data_with_labels)))
    
    #sse calculation (maximizing log-likelihood equivalents instead of minimizing euclidean distance)
    cluster_id3 <- as.matrix(c(1:k2))
    s_prob_mat1 <- t(as.matrix(apply(cluster_id3,1,function(x) BayesNumerator(ionosphere_data,cluster_means_upd2[x,],cluster_cov_upd2[[x]],cluster_prior_upd2[x,]))))
    s_prob1 <- as.matrix(apply(as.matrix(apply(s_prob_mat1,2,sum)),2,function(x) log(x)))
    s_prob_m1[k2-1,seed_iter2] <- sum((s_prob1))
  }
  
}
end.time <- Sys.time()
print(end.time - start.time)

#Plotting bad("b" label) proportions for k and iterations
plot(2:5, apply(gb_errors3,1,mean),
     type="b", pch = 19, frame = FALSE, 
     main = "Ionosphere data",
     xlab="Number of clusters K",
     ylab="Total good-bad errors")

# log-likelihood plot
plot(2:5, apply(s_prob_m1,1,sum),
     type="b", pch = 19, frame = FALSE, 
     main = "Ionosphere data",
     xlab="Number of clusters K",
     ylab="Total sum of log-likelihoods")


#Boxplot - EM normal vs EM++
library(reshape2)
library(ggplot2)

em_old_errors <- melt(gb_errors2)
em_errors <- melt(gb_errors3)

box_data <- data.frame(label=matrix(rep(c("EM normal","EM++"),c(80,80)),ncol=1), rbind(em_old_errors,em_errors))


g <- ggplot(data = box_data, aes(x=Var2, y=value)) + geom_boxplot(aes(fill=label))
g + facet_wrap( ~ (Var1+1), scales="free")

#************************************************************************************************

#Part 2 -- EM using Breast Cancer Wisconsin Dataset

rm(list=ls())
#set according to your working directory
old.dir <- getwd()
setwd("C:/Users/Neha Rawat/Desktop/IU-Data Science/AML/Assignment 2")

#Breast Cancer Wisconsin dataset
wisconsin_data_with_labels <- read.csv("wisconsin_data.csv",header=FALSE)
na_value_l <- length(which(wisconsin_data_with_labels[,7] == "?"))
na_values <- which(wisconsin_data_with_labels[,7] == "?")
data_rows <- c(1:699)
accepted_rows <- as.matrix(which(!(data_rows %in% na_values)))
wisconsin_data <- t(simplify2array(as.matrix(apply(accepted_rows, 1, function(x) strtoi(wisconsin_data_with_labels[x,c(2:10)]))),higher= FALSE))
wisconsin_data_with_labels <- t(simplify2array(as.matrix(apply(accepted_rows, 1, function(x) strtoi(wisconsin_data_with_labels[x,]))),higher= FALSE))

#getting true centers
good_mean_indices <- which(wisconsin_data_with_labels[,11] == 2)
bad_mean_indices <- which(wisconsin_data_with_labels[,11] == 4)
true_mean_good_data <- as.matrix(wisconsin_data[good_mean_indices,])
true_mean_bad_data <- as.matrix(wisconsin_data[bad_mean_indices,])
true_mean_good <- apply(true_mean_good_data,2,mean)
true_mean_bad <- apply(true_mean_bad_data,2,mean)

#****************************EM using normal initialization method*************
library(mvtnorm)

#Initializations
k_values1 <- c(2:5)
max_iter1 <- 50
convergence_check1 <- 0.0000005
max_seed_iter1 <- 20
sum_error2 <- 0
gb_errors2 <- matrix(NA, nrow=max(k_values1)-1, ncol=max_seed_iter1)
mark <- 0
s_prob_m <- matrix(NA, nrow=max(k_values1)-1, ncol=max_seed_iter1)


#Expectation (E-Step) Functions
BayesNumerator <- function(all.matrix,cl.mean,cl.covar,cl.prior) {
  b_num <- as.matrix(outer(cl.prior,as.matrix(dmvnorm(all.matrix,cl.mean,cl.covar)))) 
  return(b_num)
}

BayesWeightComplete <- function(k.value,all.matrix,k.mean,k.covar,k.prior){
  k_vec <- as.matrix(c(1:k.value))
  by_wt <- t(as.matrix(apply(k_vec, 1 ,function(x) BayesNumerator(all.matrix,k.mean[x,],k.covar[[x]],k.prior[x,]))))
  by_prop <- as.matrix(apply(by_wt,2,function(x) x/sum(x)))
  return(by_prop)
}

#Maximization (M-Step) Functions
ParametersUpdate <- function(e.matrix,all.matrix){
  mean.total <- matrix(data=0,nrow=nrow(e.matrix),ncol=ncol(all.matrix))
  cov.total <- list()
  #mean part
  mean.num <- as.matrix(e.matrix %*% all.matrix)
  for(i in c(1:nrow(e.matrix))){
    if (sum(e.matrix[i,]) == 0){ p <- 1
    }else
    {p <-sum(e.matrix[i,]) 
    }
    mean.total[i,] <- mean.num[i,]/p
  }
  #priors part
  prior.total <- matrix(apply(e.matrix,1,mean),nrow=nrow(e.matrix))
  #covariance part
  for(i in c(1:nrow(e.matrix))){
    singular <- 0
    if (sum(e.matrix[i,]) == 0){ p <- 1
    }else
    {p <-sum(e.matrix[i,]) 
    }
    cov.total1 <- matrix(0,nrow=ncol(all.matrix),ncol=ncol(all.matrix))
    cov.diff <- as.matrix(t(t(all.matrix)-mean.total[i,]))
    for(j in c(1:nrow(all.matrix))){
      cov.step1 <- matrix(outer(e.matrix[i,j],tcrossprod(cov.diff[j,])),ncol=ncol(all.matrix),nrow=ncol(all.matrix))
      cov.total1 <- cov.total1 + cov.step1
    }
    cov.total[[i]]<- matrix(cov.total1/p,ncol=ncol(all.matrix),nrow=ncol(all.matrix))
    
    #dealing with singular covariance matrices by reseeding
    if (det(cov.total[[i]]) == 0){
      cov.total[[i]] <- diag(ncol(all.matrix))
      mean.total[i,] <- all.matrix[sample(nrow(all.matrix),1,replace=FALSE), ]
      prior.total[i] <- 1/nrow(e.matrix)
    }
  }
  #return all 3 for all clusters
  parameters <- list(meanupd=mean.total,covupd=cov.total,priorupd=prior.total)  
  return(parameters)
}

#Helper function to calculate Euclidean distance between a vector and matrix
EuclideanDistanceMini <- function(one.vec, all.matrix) {
  diff.matrix <- t(t(all.matrix) - one.vec)
  this.row <- apply(diff.matrix, 1, function(x) sum(x**2)**0.5)
  return(this.row)
}

#Function to calculate Euclidean distance between two full matrices
EuclideanDistanceMain <- function(matrix_one, matrix_two) {
  dist.matrix <- apply(matrix_one, 1, function(x) EuclideanDistanceMini(x, matrix_two))
  dist.matrix <- matrix(dist.matrix, nrow=nrow(matrix_one), byrow=TRUE)
  return(dist.matrix)
}


#Function to calculate the error proportion i.e bad/(good+bad)
GoodBadError <- function(k,assigned_clusters,cluster.centers,all.matrix) {
  good_bad_ratio <- list()
  if (k==2){
    for (k.value in c(1:k)) {
      good_bad_ratio[[k.value]] <- length(which(as.matrix(all.matrix[which(assigned_clusters == k.value),11]) == 4))/length(which(assigned_clusters == k.value))
    }
  }else {
    dist_mat_good <- EuclideanDistanceMini(true_mean_good,cluster.centers)
    dist_mat_bad <- EuclideanDistanceMini(true_mean_bad,cluster.centers)
    dist_mat <- rbind(dist_mat_good,dist_mat_bad)
    dist_assign <- apply(dist_mat,2,which.min)
    assigned_clust_temp <- assigned_clusters
    for(i in 1:length(assigned_clust_temp)){
      assigned_clust_temp[i] <- dist_assign[assigned_clust_temp[i]]
    }
    for (k.value in c(1:2)) {
      if (length(which(assigned_clust_temp == k.value)) != 0){
        good_bad_ratio[[k.value]] <- length(which(as.matrix(all.matrix[which(assigned_clust_temp == k.value),11]) == 4))/length(which(assigned_clust_temp == k.value))
      } else {
        good_bad_ratio[[k.value]] <- 0
      }
    }
  } 
  return(good_bad_ratio)
}

start.time <- Sys.time()
#Start main code for the algo
for (k1 in k_values1) {
  print(c("K-Value: ", k1))
  for (seed_iter1 in 1:max_seed_iter1) {
    print(c("Initialization Number: ", seed_iter1))
    
    #randomly initializing means
    cluster_means_initial <- as.matrix(wisconsin_data[sample(nrow(wisconsin_data),k1,replace=FALSE), ])
    cluster_cov_initial <- rep(list(diag(ncol(wisconsin_data))),k1)
    cluster_prior_initial <- matrix(rep(1/k1),nrow=k1,ncol=1)
    
    converged1 = FALSE
    
    while (mark != 1)
    {
      for (iter1 in 1:max_iter1) {
        
        #Running E-Step
        E_step_mat <- BayesWeightComplete(k1,wisconsin_data,cluster_means_initial,cluster_cov_initial,cluster_prior_initial)
        
        #Running M-Step
        M_step_list <- ParametersUpdate(E_step_mat,wisconsin_data)
        cluster_means_upd <- M_step_list$meanupd
        cluster_cov_upd <- M_step_list$covupd
        cluster_prior_upd <- M_step_list$priorupd
        
        #Calculate error/movement of means
        error1 <- EuclideanDistanceMain(cluster_means_initial, cluster_means_upd)
        for(f1 in c(1:k1)) {
          sum_error2 <- sum_error2 + error1[f1,f1]}
        sum_error21 <- sum_error2
        sum_error2 <- 0
        
        #Check convergence criteria
        if(sum_error21 <= convergence_check1 ) {
          converged1 <- TRUE
          break
        }
        else {
          cluster_means_initial <- cluster_means_upd
          cluster_cov_initial <- cluster_cov_upd
          cluster_prior_initial <- cluster_prior_upd
        }
      }
      #getting cluster assignments
      E_step_mat <- BayesWeightComplete(k1,wisconsin_data,cluster_means_upd,cluster_cov_upd,cluster_prior_upd) 
      cluster_assignment1 <- apply(E_step_mat,2,which.max)
      #dealing with empty clusters - reseeding
      k_list <- as.matrix(c(1:k1))
      cluster_count <- as.matrix(apply(k_list,1,function(x) length(which(cluster_assignment1==x))))
      #Tagging empty clusters
      check1 <- which(cluster_count == 0) 
      check1 <- as.matrix(check1)
      if (nrow(check1) == 0){
        print("No empty clusters")
        mark <- 1
        break
      } else {
        print("Empty cluster present")
        cluster_means_initial <- cluster_means_upd
        cluster_cov_initial <- cluster_cov_upd
        cluster_prior_initial <- cluster_prior_upd
        for (b in 1:nrow(check1)){
          cluster_means_initial[check1[b,],] <- wisconsin_data[sample(nrow(wisconsin_data),1,replace=FALSE), ]
          cluster_cov_initial[[check1[b,]]] <- diag(ncol(wisconsin_data))
          cluster_prior_initial[check1[b,],] <- 1/k1
        }
      }
    }
    
    mark <- 0
    gb_errors2[k1-1, seed_iter1] <- sum(do.call(rbind,GoodBadError(k1,cluster_assignment1,cluster_means_upd,wisconsin_data_with_labels)))
    
    #sse calculation (maximizing log-likelihood equivalents instead of minimizing euclidean distance)
    cluster_id2 <- as.matrix(c(1:k1))
    s_prob_mat <- t(as.matrix(apply(cluster_id2,1,function(x) BayesNumerator(wisconsin_data,cluster_means_upd[x,],cluster_cov_upd[[x]],cluster_prior_upd[x,]))))
    s_prob <- as.matrix(apply(as.matrix(apply(s_prob_mat,2,sum)),2,function(x) log(x)))
    s_prob_m[k1-1,seed_iter1] <- sum((s_prob))
  }
  
}
end.time <- Sys.time()
print(end.time - start.time)

#Plotting bad(4 label) proportions for k and iterations
plot(2:5, apply(gb_errors2,1,mean),
     type="b", pch = 19, frame = FALSE, 
     main = "Breast Cancer Wisconsin data",
     xlab="Number of clusters K",
     ylab="Total good-bad errors")

# log-likelihood plot
plot(2:5, apply(s_prob_m,1,sum),
     type="b", pch = 19, frame = FALSE, 
     main = "Breast Cancer Wisconsin data",
     xlab="Number of clusters K",
     ylab="Total sum of log-likelihoods")

