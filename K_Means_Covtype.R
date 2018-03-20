#!/usr/bin/env Rscript
#Question 3a- KMeans on subset of Covtype Dataset

#**************************************KMeans Covtype Implementation***************************
rm(list=ls())
#set according to your working directory
#old.dir <- getwd()
#setwd("C:/Users/Neha Rawat/Desktop/IU-Data Science/Data Mining/Assignment 2")

#for parallel computing to reduce time
library(parallel)

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


#Function to check purity of clusters
PurityCheck <- function(k,assigned_clusters, all.matrix) {
  purity.check <- numeric()
  for (k.value in c(1:k)) {
    label1 <- length(which(as.matrix(all.matrix[which(assigned_clusters == k.value),55]) == 1))/length(which(assigned_clusters == k.value))
    label2 <- length(which(as.matrix(all.matrix[which(assigned_clusters == k.value),55]) == 2))/length(which(assigned_clusters == k.value))
    label3 <- length(which(as.matrix(all.matrix[which(assigned_clusters == k.value),55]) == 3))/length(which(assigned_clusters == k.value))
    label4 <- length(which(as.matrix(all.matrix[which(assigned_clusters == k.value),55]) == 4))/length(which(assigned_clusters == k.value))
    label5 <- length(which(as.matrix(all.matrix[which(assigned_clusters == k.value),55]) == 5))/length(which(assigned_clusters == k.value))
    label6 <- length(which(as.matrix(all.matrix[which(assigned_clusters == k.value),55]) == 6))/length(which(assigned_clusters == k.value))
    label7 <- length(which(as.matrix(all.matrix[which(assigned_clusters == k.value),55]) == 7))/length(which(assigned_clusters == k.value))
    purity.check[k.value] <- max(label1,label2,label3,label4,label5,label6,label7)
    purity.check[is.na(purity.check)] <- 0
  }
  return(purity.check)
}

#Main Function
#Clusters for parallel computing
cl = makeCluster(detectCores())

#cov dataset
cov_data_w_labels <- read.csv("covtype.csv",header=FALSE)
set.seed(12)
cov_data_w_labels <- as.matrix(cov_data_w_labels[sample(nrow(cov_data_w_labels),1500,replace=FALSE), ])
cov_data <- cov_data_w_labels[,-55]

#Uncomment when using command line
#args <- commandArgs(trailingOnly = TRUE)
#k_one <- args[1]
#k_two <- args[2]
#k_three <- args[3]
#k_values <- c(k_one,k_two,k_three)


#Uncomment when using R Studio
#Edit this according to choice
k_values <- c(3,20,100)

#constant and variable initializations
max_seed_iter <- 2 # Number of random initializations to try
max_iter <- 200 #maximum iterations to stop searching for k-means convergence
convergence_check <- 0.000000000005 #convergence threshold
s_sse <- 0
s_sse_init <- 0
s_distance_m1 <- matrix(NA, nrow=length(k_values), ncol=max_seed_iter)
accuracy <- matrix(NA, nrow=length(k_values), ncol=max_seed_iter)
iterations <- matrix(NA, nrow=length(k_values), ncol=max_seed_iter)
seed_array <- c(100:110)
a <- 0
ctr <- 0

start.time <- Sys.time()
#Start main code for the algo
for (k in k_values) {
  print(c("K-Value: ", k))
  a <- a+1
  
  for (seed_iter in 1:max_seed_iter) {
    print(c("Initialization Number: ", seed_iter))
    
    #randomly initializing centroids
    set.seed(seed_array[seed_iter])
    cluster_centers_initial <- as.matrix(cov_data[sample(nrow(cov_data),k,replace=FALSE), ])
    converged = FALSE
    
    #Start K-means computation for the max number of iterations and convergence criteria
    for (iter in 1:max_iter) {
      
      # Calculate distance from earlier centroids and assigning clusters
      dist.matrix <- EuclideanDistance(cl,cluster_centers_initial, cov_data)
      ctr <- ctr + nrow(cluster_centers_initial)*nrow(cov_data)
      cluster_assignment <- apply(dist.matrix, 2, which.min)
      
      #calculate sse
      cluster_id1 <- as.matrix(c(1:k))
      s_distance <- as.matrix(apply(cluster_id1,1,function(x) sum(dist.matrix[x,which(cluster_assignment == x)])))
      s_sse_init <- sum(s_distance)
      
      # Computing new centroids
      cluster_centers1 <- UpdateCentroids(cl, k, cluster_assignment, cov_data)
      
      #Calculate cluster and new sse
      dist.matrix <- EuclideanDistance(cl,cluster_centers1, cov_data)
      ctr <- ctr + nrow(cluster_centers1)*nrow(cov_data)
      cluster_assignment <- apply(dist.matrix, 2, which.min)
      cluster_id1 <- as.matrix(c(1:k))
      s_distance <- as.matrix(apply(cluster_id1,1,function(x) sum(dist.matrix[x,which(cluster_assignment == x)])))
      s_sse <- sum(s_distance)
      
      #Check convergence criteria
      if(abs(s_sse - s_sse_init) <= convergence_check ) {
        converged = TRUE
        break
      }
      else {
        cluster_centers_initial <- cluster_centers1
      }
    }
    
    #Number of Iterations
    iterations[a, seed_iter] <- iter
    cat("Iterations: ", iterations[a, seed_iter])
    
    #Cluster Assignment
    print("Cluster Assignment:")
    print(cluster_assignment)
    
    accuracy[a, seed_iter] <- mean(PurityCheck(k,cluster_assignment,cov_data_w_labels))
    cat("Average Purity/Accuracy: ", accuracy[a, seed_iter])
    
    #Final sse for the iteration
    s_distance_m1[a, seed_iter] <- s_sse
  }
}
end.time <- Sys.time()
print(end.time - start.time)

#Distance Computations
print(ctr)
#16152000

#Plotting SSE
plot(k_values, apply(s_distance_m1,1,mean),
     type="b", pch = 19, frame = FALSE, 
     main = "covtype data",
     xlab="Number of clusters K",
     ylab="Average sum of squared errors")

#Plotting Average Accuracy
plot(k_values, apply(accuracy,1,mean),
     type="b", pch = 19, frame = FALSE, 
     main = "covtype data",
     xlab="Number of clusters K",
     ylab="Average Accuracy")

#Plotting Average Number of Iterations among all restarts
plot(k_values, apply(iterations,1,mean),
     type="b", pch = 19, frame = FALSE, 
     main = "covtype data",
     xlab="Number of clusters K",
     ylab="Average number of iterations among all restarts")
