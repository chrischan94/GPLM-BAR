#This script DOES NOT perform variable selection on the alphas
#p=300,n=600
#Our model doesn't include an intercept for the simulation of the response variable

rm(list=ls())

#Load necessary libraries
library(BrokenAdaptiveRidge)
library(mvtnorm)
library(splines2)
library(glmnet) #for Lasso and Adaptive Lasso
library(latex2exp)

#Global Variables
random.seed = 123
rho = 0.25  
p = 300 # number of high-dimensional covariates
n = 600 # sample size
pz = 5 # number of non-zero high-dimensional parameters
pnonz = p - pz 
qW = 5 # number of categorical variables
qZ = 4 # number of non-linear variables
B = 20 # number of replications
prob = 0.5 
beta <- c(1,-1,rep(0,pnonz),-1,0.75,0.75) #vector of high-dimensional "genetic" covariates
alpha <- c(1,-0.5,-0.5,0.75,-1) #vector of alpha
nzpos <- c(1,2,(p-2):p) #index position of non-zero signals in beta
zpos <- (1:p)[-nzpos] #index position of zero signals in beta
b = 3 #number of basis functions (m_j)

######
#l1,l2,l3,r1,r2,r3,adj1,adj2 are control sequences for the plotting of our non-linear functions
l1 = 0.85 #left
r1 = 5.15 #right
l2 = -0.05
r2 = 1.05
l3 = -3.15
r3 = 1.15
adj1 = 0.2
adj2 = 0.05
#######

int = TRUE #returns Bernstein Polynomials with intercept 
reg.estBARAIC <- matrix(0,B,p+qW) #stores the regression estimates (beta + alpha) of BAR
reg.estLasso <- matrix(0,B,p+qW) #stores the regression estimates (beta + alpha) of Lasso penalty
reg.estALasso <- matrix(0,B,p+qW) #stores the regression estimates (beta + alpha) of Alasso penalty
reg.estBARBIC <- matrix(0,B,p+qW) #stores the regression estimates (beta + alpha) of SCAD
reg.estMLE <- matrix(0,B,pz+qW) #stores regression estimates of MLE(Oracle)

#control sequences of Z1-Z4
Z1.c <- seq(l1+adj1,r1-adj1,by=0.001) 
Z2.c <- seq(l2+adj2,r2-adj2,by=0.001) 
Z3.c <- seq(l2+adj2,r2-adj2,by=0.001)
Z4.c <- seq(l3+adj1,r3-adj1,by=0.001)

nonlin.estBARAIC1 <- matrix(0,B,length(Z1.c)) #stores the gammas of Bernstein polynomials of $\psi_1(Z_1)$ using AIC penalty
nonlin.estBARAIC2 <- matrix(0,B,length(Z2.c)) 
nonlin.estBARAIC3 <- matrix(0,B,length(Z3.c)) 
nonlin.estBARAIC4 <- matrix(0,B,length(Z4.c)) 
nonlin.estBARBIC1 <- matrix(0,B,length(Z1.c)) #stores the gammas of Bernstein polynomials of $\psi_1(Z_1)$ using BIC penalty
nonlin.estBARBIC2 <- matrix(0,B,length(Z2.c)) 
nonlin.estBARBIC3 <- matrix(0,B,length(Z3.c)) 
nonlin.estBARBIC4 <- matrix(0,B,length(Z4.c)) 

#stores the median mean squared error (MMSE)
MSE.BARAIC <- numeric(B) 
MSE.Lasso <- numeric(B) 
MSE.ALasso <- numeric(B)
MSE.BARBIC <- numeric(B)
MSE.GLMMLE <- numeric(B)

set.seed(random.seed)
a=1

while(a <= B){
  
  mu_X <- rep(0,p) #mean vector
  Sigma_X <- matrix(0, nrow = p, ncol = p) #variance-covariance matrix
  
  for(ii in 1:p){
    for(jj in 1:p){
      Sigma_X[ii,jj] = rho^(abs(ii-jj))
    }
  }
  
  X <- rmvnorm(n, mean=mu_X, sigma=Sigma_X)
  W <- replicate(qW, rbinom(n, size=1, prob=prob))
  Z1 <- runif(n, l1, r1)
  Z <- replicate(qZ-1, runif(n,l2,r2))
  Z2 <- Z[,1]
  Z3 <- Z[,2]
  Z4 <- runif(n, l3,r3)
  
  #transform to non-linear effects
  Z1_t <- sapply(1:n, FUN=function(i) 0.1*(Z1[i]-3)^2 )
  Z2_t <- sapply(1:n, FUN=function(i) 0.2*(cos(2*pi*Z2[i])+1) )
  Z3_t <- sapply(1:n, FUN=function(i) 0.2*sin(2*pi*Z3[i]) )
  Z4_t <- sapply(1:n, FUN=function(i) 0.2*(Z4[i]+1)^3 )
  
  Y <- rbinom(n, 1, 1/(1+exp(-X%*%beta - W%*%alpha - Z1_t - Z2_t - Z3_t - Z4_t)))
  
  #Data estimation using Bernstein polynomials
  Z1_b <- bernsteinPoly(Z1, degree=b, intercept=int)
  Z2_b <- bernsteinPoly(Z2, degree=b, intercept=int)
  Z3_b <- bernsteinPoly(Z3, degree=b, intercept=int)
  Z4_b <- bernsteinPoly(Z4, degree=b, intercept=int)
  
  Z1_bt <- Z1_b - matrix(rep(predict(Z1_b,newx=3),n),nrow=n,byrow=TRUE)
  Z2_bt <- Z2_b - matrix(rep(predict(Z2_b,newx=0.5),n),nrow=n,byrow=TRUE) 
  Z3_bt <- Z3_b - matrix(rep(predict(Z3_b,newx=0.5),n),nrow=n,byrow=TRUE)
  Z4_bt <- Z4_b - matrix(rep(predict(Z4_b,newx=-1),n),nrow=n,byrow=TRUE)
  
  Xmat <- cbind(X,W,Z1_bt,Z2_bt,Z3_bt,Z4_bt)
  cData <- createCyclopsData(Y ~ Xmat, modelType = "lr")
  PrAIC <- createBarPrior(penalty=2, exclude = c(1,(p+2):(p+qW+qZ*b+qZ*as.integer(int)+1)), initialRidgeVariance = 1)
  PrBIC <- createBarPrior(penalty=log(n), exclude = c(1,(p+2):(p+qW+qZ*b+qZ*as.integer(int)+1)), initialRidgeVariance = 1)
  FitBARAIC <- fitCyclopsModel(cData, prior=PrAIC) 
  FitBARBIC <- fitCyclopsModel(cData, prior=PrBIC)
  
  cvRidge <- cv.glmnet(x=Xmat, y=Y, family="binomial", type.measure="deviance", alpha=0)
  cvLasso <- cv.glmnet(x=Xmat, y=Y, family="binomial", type.measure="deviance", alpha=1, penalty.factor=c(rep(1,p),rep(0,qW+qZ*(b+as.integer(int)))))
  FitLasso <- glmnet(x=Xmat, y=Y, family="binomial", alpha=1, lambda=cvLasso$lambda.1se, penalty.factor=c(rep(1,p),rep(0,qW+qZ*(b+as.integer(int)))))
  coefRidge <- coef(cvRidge, s=cvRidge$lambda.min)[-1]
  cvAlasso <- cv.glmnet(x=Xmat, y=Y, family="binomial", type.measure="deviance", alpha = 1, penalty.factor = 1/abs(coefRidge))
  Alassofit <- glmnet(x=Xmat, y = Y, family = "binomial", alpha = 1, lambda = cvAlasso$lambda.1se,penalty.factor = c(1/abs(coefRidge[1:p]), rep(0,qW+qZ*(b+as.integer(int)))))
  Xmat_red <- Xmat[,-zpos] #reduces the data matrix to the "known" betas with signals
  GLMest <- coef(glm(Y ~ Xmat_red, family = "binomial"))[(2:(pz+qW+1))]
  
  reg.estBARAIC[a,] <- coef(FitBARAIC)[-c(1,(p+qW+2):(p+qW+qZ*b+qZ*as.integer(int)+1))] 
  reg.estLasso[a,] <- as.matrix(t(FitLasso$beta)[1:(p+qW)])
  reg.estALasso[a,] <- as.matrix(t(Alassofit$beta)[1:(p+qW)])
  reg.estBARBIC[a,] <- coef(FitBARBIC)[-c(1,(p+qW+2):(p+qW+qZ*b+qZ*as.integer(int)+1))] 
  reg.estMLE[a,] <- GLMest
  
  p1 <- predict(object=Z1_b,newx=Z1.c)-matrix(rep(predict(Z1_b,newx=3),length(Z1.c)), nrow=length(Z1.c), byrow=TRUE)
  p2 <- predict(object=Z2_b,newx=Z2.c)-matrix(rep(predict(Z2_b,newx=0.5), length(Z2.c)), nrow=length(Z2.c), byrow=TRUE)
  p3 <- predict(object=Z3_b,newx=Z3.c)-matrix(rep(predict(Z3_b,newx=0.5), length(Z3.c)), nrow=length(Z3.c), byrow=TRUE)
  p4 <- predict(object=Z4_b,newx=Z4.c)-matrix(rep(predict(Z4_b,newx=-1),length(Z4.c)), nrow=length(Z4.c), byrow=TRUE)

  #######
  #Manually change the chunk of code below if you want to have a different number of non-linear functions.               
  gammaAIC1_hat <- coef(FitBARAIC)[(p+qW+2):(p+qW+b+as.integer(int)+1)]
  gammaAIC2_hat <- coef(FitBARAIC)[(p+qW+b+as.integer(int)+2):(p+qW+2*b+2*as.integer(int)+1)]
  gammaAIC3_hat <- coef(FitBARAIC)[(p+qW+2*b+2*as.integer(int)+2):(p+qW+3*b+3*as.integer(int)+1)]
  gammaAIC4_hat <- coef(FitBARAIC)[(p+qW+3*b+3*as.integer(int)+2):(p+qW+qZ*b+qZ*as.integer(int)+1)]
  gammaBIC1_hat <- coef(FitBARBIC)[(p+qW+2):(p+qW+b+as.integer(int)+1)]
  gammaBIC2_hat <- coef(FitBARBIC)[(p+qW+b+as.integer(int)+2):(p+qW+2*b+2*as.integer(int)+1)]
  gammaBIC3_hat <- coef(FitBARBIC)[(p+qW+2*b+2*as.integer(int)+2):(p+qW+3*b+3*as.integer(int)+1)]
  gammaBIC4_hat <- coef(FitBARBIC)[(p+qW+3*b+3*as.integer(int)+2):(p+qW+qZ*b+qZ*as.integer(int)+1)]
  #######
  
  nonlin.estBARAIC1[a,] <- p1%*%gammaAIC1_hat
  nonlin.estBARAIC2[a,] <- p2%*%gammaAIC2_hat
  nonlin.estBARAIC3[a,] <- p3%*%gammaAIC3_hat
  nonlin.estBARAIC4[a,] <- p4%*%gammaAIC4_hat
  nonlin.estBARBIC1[a,] <- p1%*%gammaBIC1_hat
  nonlin.estBARBIC2[a,] <- p2%*%gammaBIC2_hat
  nonlin.estBARBIC3[a,] <- p3%*%gammaBIC3_hat
  nonlin.estBARBIC4[a,] <- p4%*%gammaBIC4_hat

  ########   
                 
  MSE.BARAIC[a] <- t(reg.estBARAIC[a,1:p]-beta)%*%Sigma_X%*%(reg.estBARAIC[a,1:p]-beta)
  MSE.Lasso[a] <- t(reg.estLasso[a,1:p]-beta)%*%Sigma_X%*%(reg.estLasso[a,1:p]-beta)
  MSE.ALasso[a] <- t(reg.estALasso[a,1:p]-beta)%*%Sigma_X%*%(reg.estALasso[a,1:p]-beta)
  MSE.BARBIC[a] <-t(reg.estBARBIC[a,1:p]-beta)%*%Sigma_X%*%(reg.estBARBIC[a,1:p]-beta)
  MSE.GLMMLE[a] <-t(reg.estMLE[a,1:pz] - beta[nzpos])%*%Sigma_X[nzpos,nzpos]%*%(reg.estMLE[a,1:pz] - beta[nzpos])
  
  a=a+1
  cat("Iteration", a-1, "is done.", "\n")
}

### Summarize Results 
true<- c(beta)
TP1 <- numeric(B) #True positives: the number of non-zero regression parameters estimated as non-zero
TN1 <- numeric(B) #True negatives: the number of zero regression parameters estimated as zero
FP1 <- numeric(B) #False positives: the number of zero regression parameters estimated as non-zero
TM1 <- numeric(B)
TP2 <- numeric(B) 
TN2 <- numeric(B) 
FP2 <- numeric(B) #False positives: the number of zero regression parameters estimated as non-zero
TM2 <- numeric(B)
TP3 <- numeric(B) #True positives: the number of non-zero regression parameters estimated as non-zero
TN3 <- numeric(B) #True negatives: the number of zero regression parameters estimated as zero
FP3 <- numeric(B) #False positives: the number of zero regression parameters estimated as non-zero
TM3 <- numeric(B)
TP4 <- numeric(B) #True positives: the number of non-zero regression parameters estimated as non-zero
TN4 <- numeric(B) #True negatives: the number of zero regression parameters estimated as zero
FP4 <- numeric(B) #False positives: the number of zero regression parameters estimated as non-zero
TM4 <- numeric(B)

for(i in 1:B){
  TP1[i] <- sum(abs(sign(reg.estBARAIC[i,nzpos])) == abs(sign(true[nzpos])))
  TN1[i] <- sum(sign(reg.estBARAIC[i,zpos]) == true[zpos])
  FP1[i] <- sum(sign(reg.estBARAIC[i,zpos]) != true[zpos])
  TP2[i] <- sum(abs(sign(reg.estLasso[i,nzpos])) == abs(sign(true[nzpos])))
  TN2[i] <- sum(sign(reg.estLasso[i,zpos]) == true[zpos])
  FP2[i] <- sum(sign(reg.estLasso[i,zpos]) != true[zpos])
  TP3[i] <- sum(abs(sign(reg.estALasso[i,nzpos])) == abs(sign(true[nzpos])))
  TN3[i] <- sum(sign(reg.estALasso[i,zpos]) == true[zpos])
  FP3[i] <- sum(sign(reg.estALasso[i,zpos]) != true[zpos])
  TP4[i] <- sum(abs(sign(reg.estBARBIC[i,nzpos])) == abs(sign(true[nzpos])))
  TN4[i] <- sum(sign(reg.estBARBIC[i,zpos]) == true[zpos])
  FP4[i] <- sum(sign(reg.estBARBIC[i,zpos]) != true[zpos])
  TM1[i] <- ifelse(TP1[i] + TN1[i] == p, 1, 0)
  TM2[i] <- ifelse(TP2[i] + TN2[i] == p, 1, 0)
  TM3[i] <- ifelse(TP3[i] + TN3[i] == p, 1, 0)
  TM4[i] <- ifelse(TP4[i] + TN4[i] == p, 1, 0)
}

m_TP1 <- mean(TP1); m_TN1 <- mean(TN1); m_FP1 <- mean(FP1); m_TM1 <- mean(TM1)
m_TP2 <- mean(TP2); m_TN2 <- mean(TN2); m_FP2 <- mean(FP2); m_TM2 <- mean(TM2)
m_TP3 <- mean(TP3); m_TN3 <- mean(TN3); m_FP3 <- mean(FP3); m_TM3 <- mean(TM3)
m_TP4 <- mean(TP4); m_TN4 <- mean(TN4); m_FP4 <- mean(FP4); m_TM4 <- mean(TM4)

MMSE.BARAIC <- median(MSE.BARAIC); MMSE.Lasso <- median(MSE.Lasso); MMSE.ALasso <- median(MSE.ALasso); MMSE.BARBIC <- median(MSE.BARBIC)
MMSE.GLMMLE <- median(MSE.GLMMLE)
SD.BARAIC <- sd(MSE.BARAIC); SD.Lasso <- sd(MSE.Lasso); SD.ALasso <- sd(MSE.ALasso); SD.BARBIC <- sd(MSE.BARBIC); SD.GLMMLE <- sd(MSE.GLMMLE)

#Estimate
estBARAIC <- apply(reg.estBARAIC[,nzpos],2,mean) 
estBARBIC <- apply(reg.estBARBIC[,nzpos],2,mean) 
estLasso <- apply(reg.estLasso[,nzpos],2,mean) 
estALasso <- apply(reg.estALasso[,nzpos],2,mean)
estOracle <- apply(reg.estMLE[,1:5], 2, mean)

b_BARAIC <- estBARAIC - beta[nzpos]
b_BARBIC <- estBARBIC- beta[nzpos]
b_Lasso <- estLasso - beta[nzpos]
b_ALasso <- estALasso - beta[nzpos]
b_Oracle <- estOracle - beta[nzpos]

sdBARAIC <- apply(reg.estBARAIC[,nzpos],2,sd)
sdBARBIC <- apply(reg.estBARBIC[,nzpos],2,sd)
sdLasso <- apply(reg.estLasso[,nzpos],2,sd)
sdALasso <- apply(reg.estALasso[,nzpos],2,sd)
sdOracle <- apply(reg.estMLE[,1:pz],2,sd)

##
estqWBARAIC = colMeans(reg.estBARAIC[,(p+1):(p+qW)])
estqWBARBIC = colMeans(reg.estBARBIC[,(p+1):(p+qW)])
estqWLasso= colMeans(reg.estLasso[,(p+1):(p+qW)])
estqWALasso = colMeans(reg.estALasso[,(p+1):(p+qW)])
estqWOracle = colMeans(reg.estMLE[,(pz+1):(pz+qW)])

alphab_BARAIC <- estqWBARAIC - alpha
alphab_BARBIC <- estqWBARBIC - alpha
alphab_Lasso <- estqWLasso - alpha
alphab_ALasso <- estqWALasso - alpha
alphab_Oracle <- estqWOracle - alpha

sdqwBARAIC = apply(reg.estBARAIC[,(p+1):(p+qW)],2,sd)
sdqwBARBIC = apply(reg.estBARBIC[,(p+1):(p+qW)],2,sd)
sdqwLasso = apply(reg.estLasso[,(p+1):(p+qW)],2,sd)
sdqwALasso = apply(reg.estALasso[,(p+1):(p+qW)],2,sd)
sdqWOracle = apply(reg.estMLE[,(pz+1):(pz+qW)],2,sd)

fileConn = "VariableSelectionResultsScenario1.txt"
write(x=paste("sample size is", n, sep=""),fileConn,append = FALSE)
write(x=paste("p is", p, sep=""),fileConn,append=TRUE)
write(x=paste("rho is", rho, sep=""),fileConn,append=TRUE)
write(x=paste("TP of BAR-AIC is :", round(m_TP1,2), sep = ""),fileConn, append = TRUE)
write(x=paste("TP of Lasso is :", round(m_TP2,2), sep = ""), fileConn, append = TRUE)
write(x=paste("TP of ALasso is :", round(m_TP3,2), sep = ""), fileConn, append = TRUE)
write(x=paste("TP of BAR-BIC is :", round(m_TP4,2), sep = ""), fileConn, append = TRUE)
write(x=paste("FP of BAR-AIC is :", round(m_FP1,2), sep = ""), fileConn, append = TRUE)
write(x=paste("FP of Lasso is :", round(m_FP2,2), sep = ""), fileConn, append = TRUE)
write(x=paste("FP of ALasso is :", round(m_FP3,2), sep = ""), fileConn, append = TRUE)
write(x=paste("FP of BAR-BIC is :", round(m_FP4,2), sep = ""), fileConn, append = TRUE)
write(x=paste("TM of BAR-AIC is :", 100*m_TM1, sep = ""), fileConn, append = TRUE)
write(x=paste("TM of Lasso is :", 100*m_TM2, sep = ""), fileConn, append = TRUE)
write(x=paste("TM of ALasso is :", 100*m_TM3, sep = ""), fileConn, append = TRUE)
write(x=paste("TM of BAR-BIC is :", 100*m_TM4, sep = ""),fileConn, append = TRUE )
write(x=paste("MMSE of BAR-AIC is :", round(MMSE.BARAIC,3), sep = ""),fileConn, append = TRUE)
write(x=paste("MMSE of Lasso is :", round(MMSE.Lasso,3), sep = ""),fileConn, append = TRUE)
write(x=paste("MMSE of ALasso is :", round(MMSE.ALasso,3), sep = ""),fileConn, append = TRUE)
write(x=paste("MMSE of BAR-BIC is :", round(MMSE.BARBIC,3), sep = ""),fileConn, append = TRUE)
write(x=paste("MMSE of Oracle is :", round(MMSE.GLMMLE,3), sep = ""),fileConn, append = TRUE)
write(x=paste("SD of BAR-AIC is :", round(SD.BARAIC,3), sep = ""),fileConn, append = TRUE)
write(x=paste("SD of Lasso is :", round(SD.Lasso,3), sep = ""),fileConn, append = TRUE)
write(x=paste("SD of ALasso is :", round(SD.ALasso,3), sep = ""),fileConn, append = TRUE)
write(x=paste("SD of SCAD is :", round(SD.BARBIC,3), sep = ""),fileConn, append = TRUE)
write(x=paste("SD of Oracle is :", round(SD.GLMMLE,3), sep = ""),fileConn, append = TRUE)

fileConn2 = "EstimationResultsScenario1.txt"
write(x=paste("Bias of Betas for BAR-AIC:", round(b_BARAIC,2), sep=""), fileConn2, append=FALSE)
write(x=paste("Bias of Betas for BAR-BIC:", round(b_BARBIC,2), sep=""), fileConn2, append=T)
write(x=paste("Bias of Betas for Lasso:", round(b_Lasso,2), sep=""), fileConn2, append=T)
write(x=paste("Bias of Betas for ALasso:", round(b_ALasso,2), sep=""), fileConn2, append=T)
write(x=paste("Bias of Betas for Oracle:", round(b_Oracle,2), sep=""), fileConn2, append=T)
write(x=paste("Bias of Alphas for BAR-AIC:", round(alphab_BARAIC,2), sep=""), fileConn2, append=T)
write(x=paste("Bias of Alphas for BAR-BIC:", round(alphab_BARBIC,2), sep=""), fileConn2, append=T)
write(x=paste("Bias of Alphas for Lasso:", round(alphab_Lasso,2), sep=""), fileConn2, append=T)
write(x=paste("Bias of Alphas for ALasso:", round(alphab_ALasso,2), sep=""), fileConn2, append=T)
write(x=paste("Bias of Alphas for Oracle:", round(alphab_Oracle,2), sep=""), fileConn2, append=T)
write(x=paste("SD of Betas for BAR-AIC:", round(sdBARAIC,2), sep=""), fileConn2, append=T)
write(x=paste("SD of Betas for BAR-BIC:", round(sdBARBIC,2), sep=""), fileConn2, append=T)
write(x=paste("SD of Betas for Lasso:", round(sdLasso,2), sep=""), fileConn2, append=T)
write(x=paste("SD of Betas for ALasso:", round(sdALasso,2), sep=""), fileConn2, append=T)
write(x=paste("SD of Betas for Oracle:", round(sdOracle,2), sep=""), fileConn2, append=T)
write(x=paste("SD of Alphas for BAR-AIC:", round(sdqwBARAIC,2), sep=""), fileConn2, append=T)
write(x=paste("SD of Alphas for BAR-BIC:", round(sdqwBARBIC,2), sep=""), fileConn2, append=T)
write(x=paste("SD of Alphas for Lasso:", round(sdqwLasso,2), sep=""), fileConn2, append=T)
write(x=paste("SD of Alphas for ALasso:", round(sdqwALasso,2), sep=""), fileConn2, append=T)
write(x=paste("SD of Alphas for Oracle:", round(sdqWOracle,2), sep=""), fileConn2, append=T)

psiZ1true <- 0.1*(Z1.c-3)^2
psiz2true <- 0.2*(cos(2*pi*Z2.c)+1)
psiz3true <- 0.2*sin(2*pi*Z3.c)
psiz4true <- 0.2*(Z4.c+1)^3 
m.predZ1AIC <- apply(nonlin.estBARAIC1,2,mean)
m.predZ2AIC <- apply(nonlin.estBARAIC2,2,mean)
m.predZ3AIC <- apply(nonlin.estBARAIC3,2,mean)
m.predZ4AIC <- apply(nonlin.estBARAIC4,2,mean) 
m.predZ1BIC <- apply(nonlin.estBARBIC1,2,mean)
m.predZ2BIC <- apply(nonlin.estBARBIC2,2,mean)
m.predZ3BIC <- apply(nonlin.estBARBIC3,2,mean)
m.predZ4BIC <- apply(nonlin.estBARBIC4,2,mean) 

##Save the R plot 
pdf("Scenario1.pdf")
par(mfrow=c(2,2))
plot(Z1.c,psiZ1true,type="l",xlab =TeX(r'($Z_1$)'),ylab=TeX(r'($\psi_1(Z_1)$)'), lwd =2.35)
lines(x=Z1.c,y=m.predZ1AIC, col="#FFCC00", lwd = 1.5)
lines(x=Z1.c,y=m.predZ1BIC, col="#0000FF", lwd = 1.5)

plot(Z2.c, psiz2true, type = "l", xlab = TeX(r'($Z_2$)'), ylab=TeX(r'($\psi_2(Z_2)$)'), lwd =2.35)
lines(x=Z2.c,y=m.predZ2AIC, col="#FFCC00", lwd = 1.5)
lines(x=Z2.c,y=m.predZ2BIC, col="#0000FF", lwd = 1.5)

plot(Z3.c,psiz3true,type="l",xlab=TeX(r'($Z_3$)'),ylab=TeX(r'($\psi_3(Z_3)$)'), lwd=2.35)
lines(x=Z3.c,y=m.predZ3AIC, col="#FFCC00", lwd = 1.5)
lines(x=Z3.c,y=m.predZ3BIC, col="#0000FF", lwd = 1.5)

plot(Z4.c,psiz4true,type="l",xlab=TeX(r'($Z_4$)'),ylab=TeX(r'($\psi_4(Z_4)$)'), lwd=2.7 )
lines(x=Z4.c,y=m.predZ4AIC, col="#FFCC00", lwd = 1.5)
lines(x=Z4.c,y=m.predZ4BIC, col="#0000FF", lwd = 1.5)

dev.off()

cat("The program ends at", date(), "\n")
