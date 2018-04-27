

set_a <- function(g,c_i){
  #given vegtor of gammas retins vector of a
  l <- length(g)
  a = rep(1,l)
  for (i in (1:l)){
    if (g[i] == 1){
      a[i] <- c_i
    }
  }
  return(a)
}
get_D <- function(g,c_i,tau){    #NEW
  a <- set_a(g,c_i)
  d <- tau*a #diagonal elements
  return(diag(d))
}

get_inv_D <- function (g,c_i,tau){
  # returns D^-1 given vector gamma
  a <- set_a(g,c_i)
  d <- 1/(tau*a) #diagonal elements
  return(diag(d))
}
get_A <- function(sigma_jm1,XtX,inv_D_jm1,R_jm1){
  A <- solve(sigma_jm1^-2*XtX + solve(inv_D_jm1) %*% solve(R_jm1) %*% solve(inv_D_jm1))
  return(A)
}


get_beta_j <- function(XtX,beta_LS_est,gamma_jm1,sigma_jm1,c_i,tau,R_jm1) { #(13)
  A_jm1 <- get_A(sigma_jm1,XtX,get_inv_D(gamma_jm1,c_i,tau),R_jm1)
  return(rmvnorm(1, mean = (A_jm1*(sigma_jm1)^-2)%*%XtX%*%beta_LS_est,sigma = A_jm1))
}
get_sigma_j <- function(n,u_gamma_jm1,lambda_gamma_jm1,s_res) rinvgamma(1,(n+u_gamma_jm1)/2,(s_res+u_gamma_jm1*lambda_gamma_jm1)/2) #(14)


# get_gamma_j <- function(beta_j,gamma_jm1,c_i,tau,p_i) {
#   #cat('\ngamma_jm1   ',gamma_jm1)
#   gamma_j <- gamma_jm1
#   #precompute to avoid computing every loop
#   dnorm_ctau <- dnorm(beta_j,mean = 0, sd = c_i*tau)
#   dnorm_tau <- dnorm(beta_j,mean = 0, sd = tau)
#   
#   for (i in sample(1:length(gamma_j))){
#     gamma_j[i]=1
#     #cat('\ngamma_j i=1 ', gamma_j)
#     a = prod((1-gamma_j)*dnorm_ctau+gamma_j*dnorm_tau)
#     #cat('\na ',a)
#     gamma_j[i]=0
#     #cat('\ngamma_j i=0 ', gamma_j)
#     b = prod((1-gamma_j)*dnorm_ctau+gamma_j*dnorm_tau)
#     #cat('\nb ',b)
#     p_bern = a/(a+b)
#     #cat('\np_bern ', p_bern)
#     gamma_j[i] = rbern(1 , p_bern) #(16)
#     #cat('\ngamma_j     ',gamma_j)
#   }
#   #print(gamma_j)
#   return(gamma_j)
# }

get_gamma_j <- function(beta_j,gamma_jm1,c_i,tau,p_i,R){
  gamma_j <- gamma_jm1
  for (i in sample(1:length(gamma_j))){
    gamma_j[i]=1
    D <- get_D(gamma_j,c_i,tau)
    a = dmvnorm(beta_j,sigma = D%*%R%*%D)
    gamma_j[i]=0
    D <- get_D(gamma_j,c_i,tau)
    b = dmvnorm(beta_j,sigma = D%*%R%*%D)
    p_bern = a/(a+b)
    gamma_j[i] = rbern(1 , p_bern) #(16)
  }
  return(gamma_j)
}

gibbsSamplingVarSelection <- function (data,m = 5000,p_i,tau,c_i,R,u_gamma,lambda_gamma){
  #data = dataframe. Target variable must be called y
  #m number of iterations (default 5000)
  #output table of results
  
  library(mvtnorm)
  library(invgamma)
  library(Rlab)
  
  
  #fit linear model
  fitP1 <- lm('y ~.',data)
  
  #always
  beta_0 <- summary(fitP1)$coefficients[,1]
  sigma_0 <- summary(fitP1)$sigma
  gamma_0 <- rep(1,length(beta_0))
  
  #functions and precomputations
  X <- as.matrix(subset(data,select = -y))
  y <- as.matrix(data$y)
  XtX <- t(X) %*% X #precompute now to avoid computing it every time
  beta_LS_est <-beta_0 
  
  
  #step 0
  gamma_jm1 <- gamma_0
  sigma_jm1 <- sigma_0
  results <- data.frame(combination=paste(gamma_0, collapse = ' '),stringsAsFactors = FALSE)
  
  for (i in 1:m){
    #step
    beta_j <- get_beta_j(XtX,beta_LS_est,gamma_jm1,sigma_jm1,c_i,tau,R)
    s_res <- sum((y -X %*% t(beta_j))^2) # CHECK THIS!!!! ##ask!
    sigma_j <- get_sigma_j(nrow(X),u_gamma,lambda_gamma,s_res)
    gamma_j <- get_gamma_j(beta_j,gamma_jm1,c_i,tau,p_i,R)
    results[nrow(results) + 1,] = c(paste(gamma_j, collapse = ' '))
    #cat('\nbeta_j=',beta_j,' sigma_j=',sigma_j,' gamma_j=',gamma_j)
    sigma_jm1 <- sigma_j
    gamma_jm1 <- gamma_j
  }
  
  results$combination <- as.factor(results$combination)
  return(sort(table(results)))
  
}