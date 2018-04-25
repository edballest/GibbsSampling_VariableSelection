#Example 5.1

#import custom functions
source("../gibbsSamplingVarSelection.R") 

#import data
library(wle)
data(hald)
HaldData <- data.frame(y = hald[,1], X0=1, X1=hald[,2],X2=hald[,3],X3=hald[,4],X4=hald[,5])


#use custom function
gibbsSamplingVarSelection(HaldData,m = 100,p_i = 1/2,tau = 1/3,c_i = 10,R = diag(rep(1,5)),u_gamma = 0,lambda_gamma = 1)
