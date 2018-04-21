
#Small problem involving 5 potential predictors

#Problem1

#generate data
set.seed(8416)

data <- data.frame(matrix(data = c(rep(1,60),rnorm(60*5, mean = 0,sd = 1)), nrow = 60, ncol = 6))
X_vars <- c('x0','X1','X2','X3','X4','X5')
names(data) <- X_vars

data['y'] <- mapply(FUN = function(a,b) a+1.2*b+rnorm(1,mean = 0, sd = 2.5),data[,'X4'],data[,'X5'])

#fit linear models
fitP1 <- lm('y ~ X1 + X2 + X3 + X4 + X5',data)

#discuss results in report
t <- summary(fitP1)


#Functions
set_a <- function(g,ci){
  #given vegtor of gammas retins vector of a
  l <- length(g)
  a = rep(1,l)
  for (i in (1:l)){
    if (g[i] == 1){
      a[i] <- ci
    }
  }
  return(a)
}
get_inv_D <- function (g,ci){
  # returns D^-1 given vector gamma
  a <- set_a(g,ci)
  d <- 1/(tau*a) #diagonal elements
  return(diag(d))
}
get_A <- function(sigma_jm1,XtX,inv_D_jm1,R_jm1){
  A <- solve(sigma_jm1^-2*XtX + solve(inv_D_jm1) %*% solve(R_jm1) %*% solve(inv_D_jm1))
  return(A)
}

library(mvtnorm)
library(invgamma)
library(Rlab)

get_beta_j <- function(XtX,beta_LS_est,gamma_jm1,sigma_jm1,ci,R_jm1) { #(13)
  A_jm1 <- get_A(sigma_jm1,XtX,get_inv_D(gamma_jm1,ci),R_jm1)
  return(rmvnorm(1, mean = (A_jm1*(sigma_jm1)^-2)%*%XtX%*%beta_LS_est,sigma = A_jm1)) 
}
get_sigma_j <- function(n,u_gamma_jm1,lambda_gamma_jm1,s_res) rinvgamma(1,(n+u_gamma_jm1)/2,(s_res+u_gamma_jm1*lambda_gamma_jm1)/2) #(14)
get_gamma_j <- function(n,ci,tau,pi) {
  gamma_j <- gamma_jm1
  for (i in 1:n){
    a = rnorm(1,ci*tau)*pi #(2) with gamma i = 1
    b = rnorm(1,tau)*(1-pi) #(2) with gamma i = 0 ##ask!
    p_bern = a/(a+b)
    if (p_bern>1){  #does this make sense??? ##ask!
      p_bern <- 1
    } else if (p_bern<0){
      p_bern <- 0
    }
    gamma_j[i] = rbern(1 , p_bern) #(16)
  }
  return(gamma_j)
}


#Parameters
#For this problem
pi <- 1/2 # fgamma follows priopr (7) and is equal to (1/2)^5
tau <- 1/3
ci <- 10
R <- diag(rep(1,6)) #Identity
u_gamma <- 0
lambda_gamma <- 1 #anything

#always
beta_0 <- summary(fitP1)$coefficients[,1]
sigma_0 <- summary(fitP1)$sigma
gamma_0 <- rep(1,length(beta_0))

#functions and precomputations
X <- as.matrix(data[X_vars])
y <- as.matrix(data$y)
XtX <- t(X) %*% X #precompute now to avoid computing it every time
beta_LS_est <-beta_0 #cte?  ##ask!


#step 0
gamma_jm1 <- gamma_0
sigma_jm1 <- sigma_0
results <- data.frame(combination=paste(gamma_0, collapse = ' '),stringsAsFactors = FALSE)

m <- 5000

for (i in 1:m){
  #step
  beta_j <- get_beta_j(XtX,beta_LS_est,gamma_jm1,sigma_jm1,ci,R)
  s_res <- mean((y -X %*% t(beta_j))^2) # CHECK THIS!!!! ##ask!
  sigma_j <- get_sigma_j(length(beta_j),u_gamma,lambda_gamma,s_res)
  gamma_j <- get_gamma_j(length(beta_j),ci,tau,pi)
  results[nrow(results) + 1,] = c(paste(gamma_j, collapse = ' '))
  
  #
  sigma_jm1 <- sigma_j
  gamma_jm1 <- gamma_j
}

results$combination <- as.factor(results$combination)
table(results)


#Problem 2

#data2 <- data
#data2$X3 <- sapply(data2$X5,function(a) a+0.15*rnorm(1,mean = 0,sd = 1)) #replace X3 to introduce extreme colinearity

#cor(data2$X3,data2$X5)

#data2['y'] <- mapply(FUN = function(a,b) a+1.2*b+rnorm(1,mean = 0, sd = 2.5),data2[,'X4'],data2[,'X5'])

#fitP2 <- lm('y ~ X1 + X2 + X3 + X4 + X5',data2)

#summary(fitP2)



