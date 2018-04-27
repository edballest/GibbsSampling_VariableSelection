#Example 4.1

#import custom functions
source("../gibbsSamplingVarSelection.R") 

#generate data
set.seed(8416)
data <- data.frame(matrix(data = c(rep(1,60),rnorm(60*5, mean = 0,sd = 1)), nrow = 60, ncol = 6))
X_vars <- c('x0','X1','X2','X3','X4','X5')
names(data) <- X_vars
data['y'] <- mapply(FUN = function(a,b) a+1.2*b+rnorm(1,mean = 0, sd = 2.5),data[,'X4'],data[,'X5'])

#use custom function
gibbsSamplingVarSelection(data,m = 1000,p_i = 1/2,tau = .33,c_i = 10,R = diag(rep(1,6)),u_gamma = 0,lambda_gamma = 1)
