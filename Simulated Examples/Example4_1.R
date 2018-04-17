
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
summary(fitP1)

#Parameters
#For this problem
f_gamma <- function(g) (.5)^5
tau <- 1/3
ci <- 10
R <- diag(rep(1,6)) #Identity
u_gamma <- 0

#always
beta_0 <- summary(fitP1)$coefficients[,1]
sigma_0 <- summary(fitP1)$coefficients[,2]
gamma_0 <- rep(1,6)

#functions and precomputations
X = as.matrix(data[X_vars])
XtX = t(X) %*% X #precompute now to avoid computing it every time

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

get_A <- function(sigma_j,XtX,inv_D,R){
  A <- solve((sigma_j^-2)*XtX + inv_D %*% solve(R) %*% inv_D)
  return(A)
}

get_A(sigma_0,XtX,get_inv_D(gamma_0,ci),R)


#Problem 2

#data2 <- data
#data2$X3 <- sapply(data2$X5,function(a) a+0.15*rnorm(1,mean = 0,sd = 1)) #replace X3 to introduce extreme colinearity

#cor(data2$X3,data2$X5)

#data2['y'] <- mapply(FUN = function(a,b) a+1.2*b+rnorm(1,mean = 0, sd = 2.5),data2[,'X4'],data2[,'X5'])

#fitP2 <- lm('y ~ X1 + X2 + X3 + X4 + X5',data2)

#summary(fitP2)



