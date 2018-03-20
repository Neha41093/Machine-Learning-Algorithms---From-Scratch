#!/usr/bin/env Rscript
#Question 3c - KMeans on Iris Dataset using different distance functions

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


#Helper function to calculate Cityblock distance between a vector and matrix
CityblockDistanceMini <- function(one.vec, all.matrix) {
  diff.matrix <- t(t(all.matrix) - one.vec)
  this.row <- apply(diff.matrix, 1, function(x) abs(sum(x)))
  return(this.row)
}

#Function to calculate Cityblock distance between two full matrices
CityblockDistance <- function(cl,matrix_one, matrix_two) {
  clusterExport(cl, list("matrix_one", "CityblockDistanceMini", "matrix_two"), envir=environment())
  dist.matrix <- parRapply(cl, matrix_one,function(x) CityblockDistanceMini(x, matrix_two))
  dist.matrix <- matrix(dist.matrix, nrow=nrow(matrix_one), byrow=TRUE)
  return(dist.matrix)
}

#Cosine Distance Function
CosineDistance <- function(cl,matrix_one,matrix_two){
  distance <- as.matrix(acos(matrix_one%*%t(matrix_two)/(sqrt(rowSums(matrix_one^2) %*% t(rowSums(matrix_two^2)))))/3.14) 
  distance[is.nan(distance)] <- 0
  return(distance)
}

#Helper function to calculate first Custom distance between a vector and matrix
C1DistanceMini <- function(one.vec, all.matrix) {
  diff.matrix <- t(t(all.matrix) - one.vec)
  this.row <- apply(diff.matrix, 1, function(x) ((sum(x[which(x > 0)])**2) + (abs(sum(x[which(x < 0)]))**2))**0.5)
  return(this.row)
}

#Function to calculate first Custom distance between two full matrices
C1Distance <- function(cl,matrix_one, matrix_two) {
  clusterExport(cl, list("matrix_one", "C1DistanceMini", "matrix_two"), envir=environment())
  dist.matrix <- parRapply(cl, matrix_one,function(x) C1DistanceMini(x, matrix_two))
  dist.matrix <- matrix(dist.matrix, nrow=nrow(matrix_one), byrow=TRUE)
  return(dist.matrix)
}

#Helper function to calculate second Custom distance between a vector and matrix
C2DistanceMini <- function(one.vec, all.matrix) {
  temp <- numeric()
  diff.matrix <- t(t(all.matrix) - one.vec)
  this.row1 <- apply(diff.matrix, 1, function(x) ((sum(x[which(x > 0)])**2) + (abs(sum(x[which(x < 0)]))**2))**0.5)
  for(i in 1:nrow(diff.matrix)){
    temp[i] <- sum(pmax(diff.matrix[i,],all.matrix[i,],one.vec))
    }
  this.row <- this.row1 / temp
  return(this.row)
}

#Function to calculate second Custom distance between two full matrices
C2Distance <- function(cl,matrix_one, matrix_two) {
  clusterExport(cl, list("matrix_one", "C2DistanceMini", "matrix_two"), envir=environment())
  dist.matrix <- parRapply(cl, matrix_one,function(x) C2DistanceMini(x, matrix_two))
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


#Main Function
#Clusters for parallel computing
cl = makeCluster(detectCores())

#Iris dataset
iris_data_w_labels <- read.csv("iris.csv",header=FALSE)
iris_data <- iris_data_w_labels[,-5]
iris_data <- data.matrix(iris_data)

#Uncomment when using command line
#args <- commandArgs(trailingOnly = TRUE)
#k_one <- args[1]
#k_two <- args[2]
#k_three <- args[3]
#k_four <- args[4]
#distance_meth <- args[5]
#k_values <- c(k_one,k_two,k_three,k_four)


#Uncomment when using R Studio
#Edit this according to choice
k_values <- c(2,3,20,50)
distance_meth <- "C2"

#constant and variable initializations
max_seed_iter <- 5 # Number of random initializations to try
max_iter <- 200 #maximum iterations to stop searching for k-means convergence
convergence_check <- 0.000000000005 #convergence threshold
s_sse <- 0
s_sse_init <- 0
s_distance_m1 <- matrix(NA, nrow=length(k_values), ncol=max_seed_iter)
iterations <- matrix(NA, nrow=length(k_values), ncol=max_seed_iter)
seed_array <- c(100:120)
a <- 0

start.time <- Sys.time()
#Start main code for the algo
for (k in k_values) {
  print(c("K-Value: ", k))
  a <- a+1
 
  for (seed_iter in 1:max_seed_iter) {
    print(c("Initialization Number: ", seed_iter))
    
    #randomly initializing centroids
    set.seed(seed_array[seed_iter])
    cluster_centers_initial <- as.matrix(iris_data[sample(nrow(iris_data),k,replace=FALSE), ])
    converged = FALSE
    
    #Start K-means computation for the max number of iterations and convergence criteria
    for (iter in 1:max_iter) {
      
      # Calculate distance from earlier centroids and assigning clusters
      dist.matrix <- match.fun(paste(distance_meth,"Distance",sep=""))(cl,cluster_centers_initial, iris_data)
      cluster_assignment <- apply(dist.matrix, 2, which.min)
      
      #calculate sse
      cluster_id1 <- as.matrix(c(1:k))
      s_distance <- as.matrix(apply(cluster_id1,1,function(x) sum(dist.matrix[x,which(cluster_assignment == x)])))
      s_sse_init <- sum(s_distance)
      
      # Computing new centroids
      cluster_centers1 <- UpdateCentroids(cl, k, cluster_assignment, iris_data)
      
      #Calculate cluster and new sse
      dist.matrix <- match.fun(paste(distance_meth,"Distance",sep=""))(cl,cluster_centers1, iris_data)
      cluster_assignment <- apply(dist.matrix, 2, which.min)
      cluster_id1 <- as.matrix(c(1:k))
      s_distance <- as.matrix(apply(cluster_id1,1,function(x) sum(dist.matrix[x,which(cluster_assignment == x)])))
      s_sse <- sum(s_distance)
      
      #Check convergence criteria
      if( abs(s_sse - s_sse_init) <= convergence_check ) {
        converged = TRUE
        break
      }
      else {
        cluster_centers_initial <- cluster_centers1
      }
    }

    #Cluster Assignment
    print("Cluster Assignment:")
    print(cluster_assignment)
    
    
    #Final sse for the iteration
    s_distance_m1[a, seed_iter] <- s_sse
  }
}
end.time <- Sys.time()
print(end.time - start.time)

#Plotting SSE
plot(k_values, apply(s_distance_m1,1,mean),
     type="b", pch = 19, frame = FALSE, 
     main = "Iris data",
     xlab="Number of clusters K",
     ylab="Average sum of squared errors")


