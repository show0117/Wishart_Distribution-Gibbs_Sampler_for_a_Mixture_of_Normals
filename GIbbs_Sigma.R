library(mvtnorm)
library(LaplacesDemon)
library(MCMCpack)

m1 = m2 = 2
phi1 = phi2 = matrix(c(2,1,1,2,5,10,10,100,25,50,50,500,
                       50,100,100,1000),ncol=2,byrow=T)
mu01 = mu02 = c(3,70)
tao01 = tao02 = 1
a1_a2 = c(2.5,2.5)
alpha = c(0.5,0.5)

mu0 = mu =  matrix(c(2,55,4.4,80),nrow=2,byrow=F)
sigma1 = sigma2 = matrix(c(0.8,7,7,70),nrow = 2)
sigma0 = sigma <- array(c(sigma1,sigma2),dim = c(2,2,2))

X <- faithful
Y <- sample(1:2,nrow(X),replace = T)
n <- as.vector(table(Y))

n_iter = 100

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
  
  mu_bar_1 = (tao01*mu01+c(n[1]*mean(X[,1][Y==1]),n[1]*mean(X[,2][Y==1])))/(n[1]+tao01)
  mu_bar_2 = (tao02*mu02+c(n[2]*mean(X[,1][Y==2]),n[2]*mean(X[,2][Y==2])))/(n[2]+tao02)
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
  sigma_1 = rinvwishart(n[1]+m1,phi1[(1+2*(e-1)):(2+2*(e-1)),]+lambda_1+rest_1)
  sigma_2 = rinvwishart(n[2]+m2,phi2[(1+2*(e-1)):(2+2*(e-1)),]+lambda_2+rest_2)
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

#£U111
x =  c(1:100)
y = data_EM$sigma1.1.1[1:100]
y2 = data_Gibbs$sigma1.1.1[1:100]
y3 = data_Gibbs_1$sigma1.1.1[1:100]
y4 = data_Gibbs_2$sigma1.1.1[1:100]
y5 = data_Gibbs_3$sigma1.1.1[1:100]
plot(x,y,xlab = 'iteration',ylab = '£U',type = "l",main = '£U111')
lines(x,y2,col = 'blue')
lines(x,y3,col = 'red')
lines(x,y4,col = 'gray')
lines(x,y5,col = 'green')
legend(80, 0.8, legend=c("EM","£U1",'£U2','£U3','£U4'),
       col=c( 'black',"blue","red",'gray','green'), lty=1, cex=0.8)

#£U112
x =  c(1:100)
y = data_EM$sigma1.2.2[1:100]
y2 = data_Gibbs$sigma1.2.2[1:100]
y3 = data_Gibbs_1$sigma1.2.2[1:100]
y4 = data_Gibbs_2$sigma1.2.2[1:100]
y5 = data_Gibbs_3$sigma1.2.2[1:100]
plot(x,y,xlab = 'iteration',ylim = c(20,70),ylab = '£U',type = "l",main = '£U122')
lines(x,y2,col = 'blue')
lines(x,y3,col = 'red')
lines(x,y4,col = 'gray')
lines(x,y5,col = 'green')
legend(80, 70, legend=c("EM","£U1",'£U2','£U3','£U4'),
       col=c( 'black',"blue","red",'gray','green'), lty=1, cex=0.8)
plot(faithful$eruptions,faithful$waiting)
F