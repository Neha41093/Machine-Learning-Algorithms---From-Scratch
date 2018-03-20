#Question 3.2 -- KMeans update and EM comparison (no cov and prior update)

#*******************Ringnorm dataset****************************************
#*******************Using K-means algorithm********************************

rm(list=ls())
#set according to your working directory
old.dir <- getwd()
setwd("C:/Users/Neha Rawat/Desktop/IU-Data Science/AML/Assignment 2")

#for parallel computing to reduce time
library(parallel)

#Ringnorm dataset
ringnorm_data_with_labels <- read.csv("ringnorm_data.csv",header=FALSE)
ringnorm_data <- as.matrix(ringnorm_data_with_labels[,2:21])
#getting true centers
good_mean_indices <- which(ringnorm_data_with_labels[,1] == 1)
bad_mean_indices <- which(ringnorm_data_with_labels[,1] == -1)
true_mean_good_data <- as.matrix(ringnorm_data[good_mean_indices,])
true_mean_bad_data <- as.matrix(ringnorm_data[bad_mean_indices,])
true_mean_good <- apply(true_mean_good_data,2,mean)
true_mean_bad <- apply(true_mean_bad_data,2,mean)

#constant and variable initializations
k_values <- c(2:5) #list of k-values for k-means
max_seed_iter <- 20 # Number of random initializations to try
max_iter <- 300 #maximum iterations to stop searching for k-means convergence
convergence_check <- 0.000000000005 #convergence threshold
gb_errors <- matrix(NA, nrow=max(k_values)-1, ncol=max_seed_iter)
cluster_centers_store <- list()
cluster_centers_store1 <- list()
sum_error <- 0
s_distance_m1 <- matrix(NA, nrow=max(k_values)-1, ncol=max_seed_iter)

#defining functions for the process

#Helper function to calculate Euclidean distance between a vector and matrix
EuclideanDistanceMini <- function(one.vec, all.matrix) {
  diff.matrix <- t(t(all.matrix) - one.vec)
  this.row <- apply(diff.matrix, 1, function(x) sum(x**2)**0.5)
  return(this.row)
}

#Function to calculate Euclidean distance between two full matrices
EuclideanDistance <- function(cl,matrix_one, matrix_two) {
  clusterExport(cl, list("matrix_one", "EuclideanDistanceMini", "matrix_two"), envir=environment())
  dist.matrix <- parRapply(cl, matrix_one,function(x) EuclideanDistanceMini(x, matrix_two))
  dist.matrix <- matrix(dist.matrix, nrow=nrow(matrix_one), byrow=TRUE)
  return(dist.matrix)
}

#Helper function to compute centroid of a particular k value
UpdateCentroidsMini <- function(k.value, assigned_clusters, all.matrix) {
  n_dimensions = dim(all.matrix)[2]
  cluster_indices <- which(assigned_clusters == k.value)
  #dealing with empty clusters i.e. maintaining clusters by reseeding
  if (length(cluster_indices) == 0) {
    centroid <- matrix(all.matrix[sample(nrow(all.matrix),1,replace=FALSE), ], ncol=n_dimensions)
    return(centroid)
  }
  centroid <- matrix(apply(matrix(all.matrix[cluster_indices,], ncol=n_dimensions), 2, sum) * 1.0 / length(cluster_indices), ncol=n_dimensions)
  return(centroid)
}


#Function to compute centroids for all the k-values
UpdateCentroids <- function(cl, k, assigned_clusters, all.matrix) {
  cluster_ids = as.matrix(c(1:k))
  clusterExport(cl, list("all.matrix", "UpdateCentroidsMini", "cluster_ids", "assigned_clusters"), envir=environment())
  cluster_centers1 <- t(simplify2array(parApply(cl, cluster_ids, MARGIN = 1, function(x) UpdateCentroidsMini(x, assigned_clusters, all.matrix)), higher=FALSE))
  return(as.matrix(cluster_centers1))
}

#Function to calculate the error proportion i.e bad/(good+bad)
GoodBadError <- function(k,assigned_clusters,cluster.centers,all.matrix) {
  good_bad_ratio <- list()
  if (k==2){
    for (k.value in c(1:k)) {
      good_bad_ratio[[k.value]] <- length(which(as.matrix(all.matrix[which(assigned_clusters == k.value),1]) == -1))/length(which(assigned_clusters == k.value))
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
        good_bad_ratio[[k.value]] <- length(which(as.matrix(all.matrix[which(assigned_clust_temp == k.value),1]) == -1))/length(which(assigned_clust_temp == k.value))
      } else {
        good_bad_ratio[[k.value]] <- 0
      }
    }
  } 
  return(good_bad_ratio)
}

#Clusters for parallel computing
cl = makeCluster(detectCores())

start.time <- Sys.time()
#Start main code for the algo
for (k in k_values) {
  print(c("K-Value: ", k))
  
  for (seed_iter in 1:max_seed_iter) {
    print(c("Initialization Number: ", seed_iter))
    
    #randomly initializing centroids
    cluster_centers_initial <- as.matrix(ringnorm_data[sample(nrow(ringnorm_data),k,replace=FALSE), ])
    cluster_centers_store1[[seed_iter]] <- cluster_centers_initial
    converged = FALSE
    
    #Start K-means computation for the max number of iterations and convergence criteria
    for (iter in 1:max_iter) {
      
      # Calculate distance from earlier centroids and assigning clusters
      dist.matrix <- EuclideanDistance(cl,cluster_centers_initial, ringnorm_data)
      cluster_assignment <- apply(dist.matrix, 2, which.min)
      
      # Computing new centroids
      cluster_centers1 <- UpdateCentroids(cl, k, cluster_assignment, ringnorm_data)
      
      #Calculate error/movement
      error <- EuclideanDistance(cl,cluster_centers_initial, cluster_centers1)
      for(f in c(1:k)) {
        sum_error <- sum_error + error[f,f]}
      sum_error1 <- sum_error
      sum_error <- 0
      
      #Check convergence criteria
      if( sum_error1 <= convergence_check ) {
        converged = TRUE
        break
      }
      else {
        cluster_centers_initial <- cluster_centers1
      }
    }
    if (converged == FALSE) {
      dist.matrix <- EuclideanDistance(cl, cluster_centers1, ringnorm_data)
      cluster_assignment <- apply(dist.matrix, 2, which.min)
    }
    gb_errors[k-1, seed_iter] <- sum(do.call(rbind,GoodBadError(k,cluster_assignment,cluster_centers1,ringnorm_data_with_labels)))
    
    #sse
    cluster_id1 <- as.matrix(c(1:k))
    s_distance <- as.matrix(apply(cluster_id1,1,function(x) sum(dist.matrix[x,which(cluster_assignment == x)])))
    s_distance_m1[k-1, seed_iter] <- sum(s_distance)
  }
  cluster_centers_store[[k-1]] <- do.call(rbind,cluster_centers_store1)
}
end.time <- Sys.time()
print(end.time - start.time)

#Plotting bad(-1 label) proportions for k and iterations
plot(2:5, apply(gb_errors,1,mean),
     type="b", pch = 19, frame = FALSE, 
     main = "Ringnorm data",
     xlab="Number of clusters K",
     ylab="Total good-bad errors")

#Plotting SSE
plot(2:5, apply(s_distance_m1,1,sum),
     type="b", pch = 19, frame = FALSE, 
     main = "Ringnorm data",
     xlab="Number of clusters K",
     ylab="Total sum of squared errors")

#*******************Using EM algorithm********************************
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
      good_bad_ratio[[k.value]] <- length(which(as.matrix(all.matrix[which(assigned_clusters == k.value),1]) == -1))/length(which(assigned_clusters == k.value))
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
        good_bad_ratio[[k.value]] <- length(which(as.matrix(all.matrix[which(assigned_clust_temp == k.value),1]) == -1))/length(which(assigned_clust_temp == k.value))
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
  cursor <- 1  
  for (seed_iter1 in 1:max_seed_iter1) {
    print(c("Initialization Number: ", seed_iter1))
    
    #randomly initializing means
    cluster_means_initial <- t(as.matrix(apply(as.matrix(c(cursor:(cursor+k1-1))),1,function(x) cluster_centers_store[[k1-1]][x,])))
    cluster_cov_initial <- rep(list(diag(ncol(ringnorm_data))),k1)
    cluster_prior_initial <- matrix(rep(1/k1),nrow=k1,ncol=1)
    
    converged1 = FALSE
    
    while (mark != 1)
    {
      for (iter1 in 1:max_iter1) {
        
        #Running E-Step
        E_step_mat <- BayesWeightComplete(k1,ringnorm_data,cluster_means_initial,cluster_cov_initial,cluster_prior_initial)
        
        #Running M-Step
        M_step_list <- ParametersUpdate(E_step_mat,ringnorm_data)
        cluster_means_upd <- M_step_list$meanupd
        cluster_cov_upd <- cluster_cov_initial
        cluster_prior_upd <- cluster_prior_initial
        
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
      E_step_mat <- BayesWeightComplete(k1,ringnorm_data,cluster_means_upd,cluster_cov_upd,cluster_prior_upd) 
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
          cluster_means_initial[check1[b,],] <- ringnorm_data[sample(nrow(ringnorm_data),1,replace=FALSE), ]
          cluster_cov_initial[[check1[b,]]] <- diag(ncol(ringnorm_data))
          cluster_prior_initial[check1[b,],] <- 1/k1
        }
      }
    }
    
    mark <- 0
    cursor <- cursor+k1
    gb_errors2[k1-1, seed_iter1] <- sum(do.call(rbind,GoodBadError(k1,cluster_assignment1, cluster_means_upd,ringnorm_data_with_labels)))
    
    #sse calculation (maximizing log-likelihood equivalents instead of minimizing euclidean distance)
    cluster_id2 <- as.matrix(c(1:k1))
    s_prob_mat <- t(as.matrix(apply(cluster_id2,1,function(x) BayesNumerator(ringnorm_data,cluster_means_upd[x,],cluster_cov_upd[[x]],cluster_prior_upd[x,]))))
    s_prob <- as.matrix(apply(as.matrix(apply(s_prob_mat,2,sum)),2,function(x) log(x)))
    s_prob_m[k1-1,seed_iter1] <- sum((s_prob))
  }
  
}
end.time <- Sys.time()
print(end.time - start.time)

#Plotting bad(-1 label) proportions for k and iterations
plot(2:5, apply(gb_errors2,1,mean),
     type="b", pch = 19, frame = FALSE, 
     main = "Ringnorm data",
     xlab="Number of clusters K",
     ylab="Total good-bad errors")

# log-likelihood plot
plot(2:5, apply(s_prob_m,1,sum),
     type="b", pch = 19, frame = FALSE, 
     main = "Ringnorm data",
     xlab="Number of clusters K",
     ylab="Total sum of log-likelihoods")

#Boxplot - kmeans vs EM
library(reshape2)
library(ggplot2)

kmeans_errors <- melt(gb_errors)
em_errors <- melt(gb_errors2)

box_data <- data.frame(label=matrix(rep(c("kmeans","EM"),c(80,80)),ncol=1), rbind(kmeans_errors,em_errors))


g <- ggplot(data = box_data, aes(x=Var2, y=value)) + geom_boxplot(aes(fill=label))
g + facet_wrap( ~ (Var1+1), scales="free")

##***************************************************************************************************##

#*******************Ionosphere dataset****************************************
#*******************Using K-means algorithm********************************

rm(list=ls())
#set according to your working directory
old.dir <- getwd()
setwd("C:/Users/Neha Rawat/Desktop/IU-Data Science/AML/Assignment 2")

#for parallel computing to reduce time
library(parallel)

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

#constant and variable initializations
k_values <- c(2:5) #list of k-values for k-means
max_seed_iter <- 20 # Number of random initializations to try
max_iter <- 300 #maximum iterations to stop searching for k-means convergence
convergence_check <- 0.000000000005 #convergence threshold
gb_errors <- matrix(NA, nrow=max(k_values)-1, ncol=max_seed_iter)
cluster_centers_store <- list()
cluster_centers_store1 <- list()
sum_error <- 0
s_distance_m1 <- matrix(NA, nrow=max(k_values)-1, ncol=max_seed_iter)

#defining functions for the process

#Helper function to calculate Euclidean distance between a vector and matrix
EuclideanDistanceMini <- function(one.vec, all.matrix) {
  diff.matrix <- t(t(all.matrix) - one.vec)
  this.row <- apply(diff.matrix, 1, function(x) sum(x**2)**0.5)
  return(this.row)
}

#Function to calculate Euclidean distance between two full matrices
EuclideanDistance <- function(cl,matrix_one, matrix_two) {
  clusterExport(cl, list("matrix_one", "EuclideanDistanceMini", "matrix_two"), envir=environment())
  dist.matrix <- parRapply(cl, matrix_one,function(x) EuclideanDistanceMini(x, matrix_two))
  dist.matrix <- matrix(dist.matrix, nrow=nrow(matrix_one), byrow=TRUE)
  return(dist.matrix)
}

#Helper function to compute centroid of a particular k value
UpdateCentroidsMini <- function(k.value, assigned_clusters, all.matrix) {
  n_dimensions = dim(all.matrix)[2]
  cluster_indices <- which(assigned_clusters == k.value)
  #dealing with empty clusters i.e. maintaining clusters by reseeding
  if (length(cluster_indices) == 0) {
    centroid <- matrix(all.matrix[sample(nrow(all.matrix),1,replace=FALSE), ], ncol=n_dimensions)
    return(centroid)
  }
  centroid <- matrix(apply(matrix(all.matrix[cluster_indices,], ncol=n_dimensions), 2, sum) * 1.0 / length(cluster_indices), ncol=n_dimensions)
  return(centroid)
}


#Function to compute centroids for all the k-values
UpdateCentroids <- function(cl, k, assigned_clusters, all.matrix) {
  cluster_ids = as.matrix(c(1:k))
  clusterExport(cl, list("all.matrix", "UpdateCentroidsMini", "cluster_ids", "assigned_clusters"), envir=environment())
  cluster_centers1 <- t(simplify2array(parApply(cl, cluster_ids, MARGIN = 1, function(x) UpdateCentroidsMini(x, assigned_clusters, all.matrix)), higher=FALSE))
  return(as.matrix(cluster_centers1))
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

#Clusters for parallel computing
cl = makeCluster(detectCores())

start.time <- Sys.time()
#Start main code for the algo
for (k in k_values) {
  print(c("K-Value: ", k))
  
  for (seed_iter in 1:max_seed_iter) {
    print(c("Initialization Number: ", seed_iter))
    
    #randomly initializing centroids
    cluster_centers_initial <- as.matrix(ionosphere_data[sample(nrow(ionosphere_data),k,replace=FALSE), ])
    cluster_centers_store1[[seed_iter]] <- cluster_centers_initial
    converged = FALSE
    
    #Start K-means computation for the max number of iterations and convergence criteria
    for (iter in 1:max_iter) {
      
      # Calculate distance from earlier centroids and assigning clusters
      dist.matrix <- EuclideanDistance(cl,cluster_centers_initial, ionosphere_data)
      cluster_assignment <- apply(dist.matrix, 2, which.min)
      
      # Computing new centroids
      cluster_centers1 <- UpdateCentroids(cl, k, cluster_assignment, ionosphere_data)
      
      #Calculate error/movement
      error <- EuclideanDistance(cl,cluster_centers_initial, cluster_centers1)
      for(f in c(1:k)) {
        sum_error <- sum_error + error[f,f]}
      sum_error1 <- sum_error
      sum_error <- 0
      
      #Check convergence criteria
      if( sum_error1 <= convergence_check ) {
        converged = TRUE
        break
      }
      else {
        cluster_centers_initial <- cluster_centers1
      }
    }
    if (converged == FALSE) {
      dist.matrix <- EuclideanDistance(cl, cluster_centers1, ionosphere_data)
      cluster_assignment <- apply(dist.matrix, 2, which.min)
    }
    gb_errors[k-1, seed_iter] <- sum(do.call(rbind,GoodBadError(k,cluster_assignment,cluster_centers1,ionosphere_data_with_labels)))
    
    #sse
    cluster_id1 <- as.matrix(c(1:k))
    s_distance <- as.matrix(apply(cluster_id1,1,function(x) sum(dist.matrix[x,which(cluster_assignment == x)])))
    s_distance_m1[k-1, seed_iter] <- sum(s_distance)
  }
  cluster_centers_store[[k-1]] <- do.call(rbind,cluster_centers_store1)
}
end.time <- Sys.time()
print(end.time - start.time)

#Plotting bad("b" label) proportions for k and iterations
plot(2:5, apply(gb_errors,1,mean),
     type="b", pch = 19, frame = FALSE, 
     main = "Ionosphere data",
     xlab="Number of clusters K",
     ylab="Total good-bad errors")

#Plotting SSE
plot(2:5, apply(s_distance_m1,1,sum),
     type="b", pch = 19, frame = FALSE, 
     main = "Ionosphere data",
     xlab="Number of clusters K",
     ylab="Total sum of squared errors")

#*******************Using EM algorithm********************************
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
  cursor <- 1  
  for (seed_iter1 in 1:max_seed_iter1) {
    print(c("Initialization Number: ", seed_iter1))
    
    #randomly initializing means
    cluster_means_initial <- t(as.matrix(apply(as.matrix(c(cursor:(cursor+k1-1))),1,function(x) cluster_centers_store[[k1-1]][x,])))
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
        cluster_cov_upd <- cluster_cov_initial
        cluster_prior_upd <- cluster_prior_initial
        
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
    cursor <- cursor+k1
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

#Boxplot - kmeans vs EM
library(reshape2)
library(ggplot2)

kmeans_errors <- melt(gb_errors)
em_errors <- melt(gb_errors2)

box_data <- data.frame(label=matrix(rep(c("kmeans","EM"),c(80,80)),ncol=1), rbind(kmeans_errors,em_errors))


g <- ggplot(data = box_data, aes(x=Var2, y=value)) + geom_boxplot(aes(fill=label))
g + facet_wrap( ~ (Var1+1), scales="free")
