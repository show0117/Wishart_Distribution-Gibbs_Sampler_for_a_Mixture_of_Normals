library(mvtnorm)

alpha = c(0.5,0.5)
mu =  matrix(c(2,55,4.4,80),nrow=2,byrow=F)
sigma1 = matrix(c(0.8,7,7,70),nrow = 2)
sigma2 = matrix(c(0.8,7,7,70),nrow = 2)
sigma <- array(c(sigma1,sigma2),dim = c(2,2,2))
X <- faithful
Y <- sample(1:2,nrow(X),replace = T)

data_EM <- data.frame('mu1.1'=mu[1,1],'mu1.2'=mu[2,1],
                   'mu2.1'=mu[1,2],'mu2.2'=mu[2,2],
                   'sigma1.1.1'=sigma[1,1,1],
                   'sigma1.1.2'=sigma[1,2,1],
                   'sigma1.2.1'=sigma[2,1,1],
                   'sigma1.2.2'=sigma[2,2,1],
                   'sigma2.1.1'=sigma[1,1,2],
                   'sigma2.1.2'=sigma[1,2,2],
                   'sigma2.2.1'=sigma[2,1,2],
                   'sigma2.2.2'=sigma[2,2,2],
                   'alpha1'=alpha[1],'alpha2'=alpha[2])

n_iter = 100
for (i in 1:n_iter){
  p_x <- c()
  for(j in 1:(length(alpha))){   
    p_x <- rbind(p_x,dmvnorm(X,mu[,j],sigma[,,j]))
  }
  p_x <- alpha*p_x
  p_y <- t(t(p_x)/colSums(p_x))
  u <- runif(nrow(X))
  judge <- apply(p_y,MARGIN = 2,cumsum)
  
  index <- c()
  for(j in 1:length(alpha)){
    new_index <- setdiff(which(u < judge[j,]),index)
    index <- c(index,new_index)
    Y[new_index] <- j
    sigma[,,j] <- (as.matrix(t(X)-mu[,j]) %*% (as.matrix(t((t(X)-mu[,j])))*(p_y[j,])))/sum(p_y[j,])
    mu[,j] <- colSums(X*(p_y[j,]))/sum(p_y[j,])
  }
  
  alpha <- rowSums(p_y)/nrow(X)
  
  new_data = c(mu,sigma,alpha)
  data_EM = rbind(data_EM,new_data)
}


