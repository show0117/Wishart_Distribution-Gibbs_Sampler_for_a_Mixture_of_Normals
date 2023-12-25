library(mvtnorm)
library(LaplacesDemon)
library(MCMCpack)

m1 = m2 = 2
phi1 = phi2 = matrix(c(2,1,1,2),nrow=2,byrow=F)
mu01 = mu02 = matrix(c(5,100,8,150,10,200),nrow=3,byrow=T)

tao01 = tao02 = 1
a1_a2 = c(2.5,2.5)
alpha = c(0.5,0.5)

mu0 = mu =  matrix(c(2,55,4.4,80),nrow=2,byrow=F)
sigma1 = sigma2 = matrix(c(0.8,7,7,70),nrow = 2)
sigma0 = sigma <- array(c(sigma1,sigma2),dim = c(2,2,2))

X <- faithful
Y <- sample(1:2,nrow(X),replace = T)
n <- as.vector(table(Y))
n_iter = 1000

for (e in 1:3){
  data <- data.frame('mu1.1'=mu[1,1],'mu1.2'=mu[2,1],
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
  for(i in 1:n_iter){
    p_x <- matrix(c(dmvnorm(X,mu[,1],sigma[,,1]),dmvnorm(X,mu[,2],sigma[,,2])),nrow = 2,byrow = T)
    p_y <- t(t(alpha*p_x)/colSums(alpha*p_x))
    u <- runif(length(X))
    Y[u <= p_y[1,]] <- 1
    Y[!(u <= p_y[1,])] <- 2
    n <- as.vector(table(Y))
    
    mu_bar_1 = (tao01*mu01[e,]+c(n[1]*mean(X[,1][Y==1]),n[1]*mean(X[,2][Y==1])))/(n[1]+tao01)
    mu_bar_2 = (tao02*mu02[e,]+c(n[2]*mean(X[,1][Y==2]),n[2]*mean(X[,2][Y==2])))/(n[2]+tao02)
    mu_bar = matrix(c(mu_bar_1,mu_bar_2),nrow = 2,byrow = F)
    sigma_bar_1 = sigma[,,1]/(tao01+n[1])
    sigma_bar_2 = sigma[,,2]/(tao02+n[2])
    sigma_bar = array(c(sigma_bar_1,sigma_bar_2),dim = c(2,2,2))
    mu =  matrix(c(mvrnorm(1,mu_bar[,1],sigma_bar[,,1]),mvrnorm(1,mu_bar[,2],sigma_bar[,,2])),nrow = 2,byrow = F)
    
    matrix_1 = matrix(c(X[,1][Y==1]-mean(X[,1][Y==1]),X[,2][Y==1]-mean(X[,2][Y==1])),ncol = 2)
    matrix_2 = matrix(c(X[,1][Y==2]-mean(X[,1][Y==2]),X[,2][Y==2]-mean(X[,2][Y==2])),ncol = 2)
    lambda_1 = matrix(c(0,0,0,0),nrow = 2)
    for(j in 1:length(matrix_1[,1])){
      lambda_1 = lambda_1 + matrix_1[j,]%*%t(matrix_1[j,])
    }
    lambda_2 = matrix(c(0,0,0,0),nrow = 2)
    for(k in 1:length(matrix_2[,1])){
      lambda_2 = lambda_2 + matrix_2[k,]%*%t(matrix_2[k,])
    }
    x_bar_1 = c(mean(X[,1][Y==1]),mean(X[,2][Y==1]))
    x_bar_2 = c(mean(X[,1][Y==2]),mean(X[,2][Y==2]))
    rest_1 = (n[1]*tao01/(n[1]+tao01))*((x_bar_1-mu_bar[,1])%*%t(x_bar_1-mu_bar[,1]))
    rest_2 = (n[2]*tao02/(n[2]+tao02))*((x_bar_2-mu_bar[,2])%*%t(x_bar_2-mu_bar[,2]))
    sigma_1 = rinvwishart(n[1]+m1,phi1+lambda_1+rest_1)
    sigma_2 = rinvwishart(n[2]+m2,phi2+lambda_2+rest_2)
    sigma = array(c(sigma_1,sigma_2),dim = c(2,2,2))
    
    alpha <- as.vector(rdirichlet(1,alpha = a1_a2+n))
    new_data = c(mu_bar_1,mu_bar_2,sigma_1,sigma_2,alpha)
    data = rbind(data,new_data)
    print(i)
  }
  if (e==1){
    data_Gibbs_1 = data
  }else if (e==2){
    data_Gibbs_2 = data
  }else{
    data_Gibbs_3 = data
  }
}

#for (i in colnames(data_Gibbs_2)){
  #dt = data_Gibbs_2[i][52:101,]
  #print(i)
  #print(mean(dt))
  #print(var(dt))
#}

#£g12
x =  c(1:100)
y = data_EM$mu1.2[1:100]
y2 = data_Gibbs$mu1.2[1:100]
y3 = data_Gibbs_1$mu1.2[1:100]
y4 = data_Gibbs_2$mu1.2[1:100]
y5 = data_Gibbs_3$mu1.2[1:100]
plot(x,y,xlab = 'iteration',ylab = '£g',type = "l",ylim = c(54, 57),main = '£g21')
lines(x,y2,col = 'blue')
lines(x,y3,col = 'red')
lines(x,y4,col = 'gray')
lines(x,y5,col = 'green')
legend(65, 57, legend=c("EM","£g = [3,70]","£g = [5,100]",
                         "£g = [8,150]","£g = [10,200]"),
       col=c( 'black',"blue","red",'gray','green'), lty=1, cex=0.8)

#£g21
x =  c(1:100)
y = data_EM$mu2.1[1:100]
y2 = data_Gibbs$mu2.1[1:100]
y3 = data_Gibbs_1$mu2.1[1:100]
y4 = data_Gibbs_2$mu2.1[1:100]
y5 = data_Gibbs_3$mu2.1[1:100]
plot(x,y,xlab = 'iteration',ylab = '£g',type = "l",ylim = c(4.25, 4.4),main = '£g21')
lines(x,y2,col = 'blue')
lines(x,y3,col = 'red')
lines(x,y4,col = 'gray')
lines(x,y5,col = 'green')
legend(65, 4.4, legend=c("EM","£g = [3,70]","£g = [5,100]",
                        "£g = [8,150]","£g = [10,200]"),
       col=c( 'black',"blue","red",'gray','green'), lty=1, cex=0.8)

#£g22
x =  c(1:100)
y = data_EM$mu2.2[1:100]
y2 = data_Gibbs$mu2.2[1:100]
y3 = data_Gibbs_1$mu2.2[1:100]
y4 = data_Gibbs_2$mu2.2[1:100]
y5 = data_Gibbs_3$mu2.2[1:100]
plot(x,y,xlab = 'iteration',ylab = '£g',type = "l",ylim = c(79.5, 81.5),main = '£g22')
lines(x,y2,col = 'blue')
lines(x,y3,col = 'red')
lines(x,y4,col = 'gray')
lines(x,y5,col = 'green')
legend(65, 81.5, legend=c("EM","£g = [3,70]","£g = [5,100]",
                        "£g = [8,150]","£g = [10,200]"),
       col=c( 'black',"blue","red",'gray','green'), lty=1, cex=0.8)
