library(mvtnorm)
library(LaplacesDemon)
library(MCMCpack)

alpha = c(0.5,0.5)
mu0 = mu =  matrix(c(2,55,4.4,80),nrow=2,byrow=F)
sigma1 = matrix(c(0.8,7,7,70),nrow = 2)
tao1 = solve(sigma1)
sigma2 = matrix(c(0.8,7,7,70),nrow = 2)
tao2 = solve(sigma2)
psi = sigma0 = sigma <- array(c(sigma1,sigma2),dim = c(2,2,2))
tao0 = tao <- array(c(tao1,tao2),dim = c(2,2,2))
mj = 1
a1_a2 = c(1,1)
X <- faithful
Y <- sample(1:2,nrow(X),replace = T)
dmvnorm(X,mu[,1],sigma[,,1])
n <- as.vector(table(Y))

mu_result = array(c(2,55,4.4,80),dim = c(2,2))
sigma_result = array(c(sigma1,sigma2),dim = c(2,2,2))
alpha_result = c(0.5,0.5)

for(i in 1:n_iter){
  p_x <- matrix(c(dmvnorm(X,mu[,1],sigma[,,1]),dmvnorm(X,mu[,2],sigma[,,2])),nrow = 2,byrow = T)
  p_y <- t(t(alpha*p_x)/colSums(alpha*p_x))
  u <- runif(length(X))
  Y[u <= p_y[1,]] <- 1
  Y[!(u <= p_y[1,])] <- 2
  n <- as.vector(table(Y))
  
  mu_bar_1 = (det(tao0[,,1])*mu0[,1]+c(n[1]*mean(X[,1][Y==1]),n[1]*mean(X[,2][Y==1])))/(n[1]+det(tao0[,,1]))
  mu_bar_2 = (det(tao0[,,2])*mu0[,2]+c(n[2]*mean(X[,1][Y==2]),n[2]*mean(X[,2][Y==2])))/(n[2]+det(tao0[,,2]))
  mu_bar = matrix(c(mu_bar_1,mu_bar_2),nrow = 2,byrow = F)
  sigma_bar_1 = sigma[,,1]/(det(tao0[,,1])+n[1])
  sigma_bar_2 = sigma[,,2]/(det(tao0[,,2])+n[2])
  sigma_bar = array(c(sigma_bar_1,sigma_bar_2),dim = c(2,2,2))
  mu =  matrix(c(mvrnorm(1,mu_bar[,1],sigma_bar[,,1]),mvrnorm(1,mu_bar[,2],sigma_bar[,,2])),nrow = 2,byrow = F)
  
  matrix_1 = matrix(c(X[,1][Y==1]-mean(X[,1][Y==1]),X[,2][Y==1]-mean(X[,2][Y==1])),ncol = 2)
  matrix_2 = matrix(c(X[,1][Y==2]-mean(X[,1][Y==2]),X[,2][Y==2]-mean(X[,2][Y==2])),ncol = 2)
  lambda_1 = matrix(c(0,0,0,0),nrow = 2)
  for(i in length(matrix_1[,1])){
    lambda_1 = lambda_1 + matrix_1[i,]%*%t(matrix_1[i,])
  }
  lambda_2 = matrix(c(0,0,0,0),nrow = 2)
  for(i in length(matrix_2[,1])){
    lambda_2 = lambda_2 + matrix_2[i,]%*%t(matrix_2[i,])
  }
  x_bar_1 = c(mean(X[,1][Y==1]),mean(X[,2][Y==1]))
  x_bar_2 = c(mean(X[,1][Y==2]),mean(X[,2][Y==2]))
  rest_1 = (n[1]*det(tao0[,,1])/(n[1]+det(tao0[,,1])))*((x_bar_1-mu_bar[,1])%*%t(x_bar_1-mu_bar[,1]))
  rest_2 = (n[2]*det(tao0[,,2])/(n[2]+det(tao0[,,2])))*((x_bar_2-mu_bar[,2])%*%t(x_bar_2-mu_bar[,2]))
  sigma_1 = rinvwishart(n[1]+1,psi[,,1]+lambda_1+rest_1)
  sigma_2 = rinvwishart(n[2]+1,psi[,,2]+lambda_2+rest_2)
  sigma = array(c(sigma_1,sigma_2),dim = c(2,2,2))
  
  alpha <- as.vector(rdirichlet(1,alpha = a1_a2+n))
  
  mu_result = array(c(mu_result,mu),dim = c(2,2,i+1))
  sigma_result = array(c(sigma_result,sigma),dim = c(2,2,2,i+1))
  alpha_result = matrix(c(alpha_result,alpha),ncol = 2,byrow = T)
}

rWishart(1,180,psi[,,2]+lambda_2+rest_2)
