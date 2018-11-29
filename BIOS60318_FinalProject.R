##########################################################
##########################################################
################ BIOS 60318 Final Project ################
##########################################################
##########################################################
rm(list=ls()) # Clear workspace
if(!is.null(dev.list())) dev.off() # Clear plots
par(mfrow=c(1,1)) # Setup plot parameters
par(ps = 12, font.lab = 1)
set.seed(1) # set random seed generator for reproducibility

# Let's load relevant libraries
#library(MASS)
#library(olsrr)
#library(rcompanion)
#library(gridExtra)
#source('~/FHplot.R', encoding = 'UTF-8')
#source('~/symplot.R', encoding = 'UTF-8')
library(ggplot2)
library(ggpmisc)
library(CorReg)
source('~/extractPVal.R', encoding = 'UTF-8')
source('~/linearNLL.R', encoding = 'UTF-8')
source('~/customNll.R', encoding = 'UTF-8')
cat("\014") # Clear console

# Let's load the data
sugarData = read.table("sugar.csv", header=TRUE, sep=",")
antibioticsData = read.table("antibiotics.csv", header=TRUE, sep=",")

##########################################################
##########################################################
##################### Part I: ANOVA ######################
##########################################################
##########################################################
# Relevel the data to have the control as the reference
releveledTrt = relevel(antibioticsData$trt, ref=4)

# Define the nll function with 3 betas set as binary factors
nll_aov = function(p,x,y){
  # Define the parameters for the betas and residual standard error
  b0=p[1]
  b1=p[2]
  b2=p[3]
  b3=p[4]
  sigma=exp(p[5])
  
  # Define separate functions as binary responses
  expected1 = b0+(x[5:8]-1)*b1
  expected2 = b0+(x[9:12]-2)*b2
  expected3 = b0+(x[13:16]-3)*b3

  # Compute nll for given Y values
  nll1=-sum(dnorm(x=y[5:8],mean=expected1,sd=sigma,log=TRUE))
  nll2=-sum(dnorm(x=y[9:12],mean=expected2,sd=sigma,log=TRUE))
  nll3=-sum(dnorm(x=y[13:16],mean=expected3,sd=sigma,log=TRUE))
  nll = nll1+nll2+nll3
  return(nll)
}

# Define the initial guesses to be the means of the data set as a good starting point
initialGuess_aov=c(mean(antibioticsData$growth[1:4]),
               -mean(antibioticsData$growth[1:4])+mean(antibioticsData$growth[5:8]),
               -mean(antibioticsData$growth[1:4])+mean(antibioticsData$growth[9:12]),
               -mean(antibioticsData$growth[1:4])+mean(antibioticsData$growth[13:16]),
               1)
# Optimize the nll function
nllAnova=optim(par=initialGuess_aov, fn=nll_aov, x=as.numeric(releveledTrt), y = antibioticsData$growth)

# Create the simplified model
simpleMod_aov = function(p,x,y){
  # Define the parameters of the simple model
  b0=p[1]
  sigma=exp(p[2])
  # Define the function of the simple model
  expected=b0
  # Compute the nll of the simple model
  nll=-sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
  return(nll)
}

# Initial guess for the simple model
simpleGuess_aov = c(mean(antibioticsData$growth[1:4]),1)

# Optimize the simplified model
nllSimple_aov = optim(par=simpleGuess_aov, fn=simpleMod_aov, x=as.numeric(releveledTrt), y = antibioticsData$growth)

# Compute the LR statistic
teststat_aov = 2*(nllSimple_aov$value-nllAnova$value)

# Compute the degrees of freedom
df_aov=length(nllAnova$par)-length(nllSimple_aov$par)

# Compute the P value from the nll functions
nllPVal_aov = 1-pchisq(teststat_aov,df_aov)

# Fit using aov() that passes each group through lm()
aov.fit = aov(growth ~ releveledTrt, data=antibioticsData)
sum.aov = summary(aov.fit)

# PValue from aov() using F statistic
aovPVal=sum.aov[[1]][["Pr(>F)"]][1]

# Parameters found from nll
nllParams_aov = c(nllAnova$par[1], nllAnova$par[2], nllAnova$par[3], nllAnova$par[4], exp(nllAnova$par[5]), nllPVal_aov)

# Parameters found from lm() and aov()
aovParams = c(aov.fit$coefficients, sigma(aov.fit), aovPVal)

# Comparison of nll versus lm()
comparison_aov = cbind(nllParams_aov, aovParams)
dimnames(comparison_aov)[[1]]=c("Control", "Treatment 1", "Treatment 2", "Treatment 3", "Residual Std Error", "P Value")
comparison_aov
write.table(comparison_aov, "clipboard", sep="\t", row.names=FALSE)

# 95% Confidence interval of the parameters 
CI_AOV = confint(aov.fit, level = 0.95)
dimnames(CI_AOV)[[1]]=c("Control", "Treatment 1", "Treatment 2", "Treatment 3")
CI_AOV

# Boxplot of the data comparing control versus treatment with 95% CI of the means for each group in red
BoxPlot(antibioticsData$growth, 
        antibioticsData$trt, 
        AnoVa = TRUE, ylab="Growth of Bacteria", 
        names=c("Treatment 1", "Treatment 2", "Treatment 3", "Control"))

##########################################################
##########################################################
############### Part II: Linear Regression ###############
##########################################################
##########################################################

# Linear log likelihood function
nllLinear<-function(p,x,y){
  B0=p[1] 
  B1=p[2] 
  sigma=exp(p[3])
  expected=B0+x*B1
  nll=-sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
  return(nll) 
}

# Initial guess for optimization
initialGuess_lin=c(1,1,1)

# Optimize the nll function
linearFit=optim(par=initialGuess_lin,fn=nllLinear,x=sugarData$sugar,y=sugarData$growthSugar)

# Simple model for nll function
nllSimple_lin<-function(p,x,y){
  B0=p[1]
  sigma=exp(p[2])
  expected=B0
  nll=-sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
  return(nll) 
}

# Initial guess for simple model
simpleGuess_lin=c(1,1)

# Optimize the simple model
simpleFit=optim(par=simpleGuess_lin,fn=nllSimple_lin,x=sugarData$sugar,y=sugarData$growthSugar)

# Compute the LR statistic
teststat_lin=2*(simpleFit$value-linearFit$value)

# Compute the degrees of freedom
df_lin=length(linearFit$par)-length(simpleFit$par)

# P value from nll function
nllPVal_lin=1-pchisq(teststat_lin,df_lin)

# Combine all parameters for nll
nllParams_lin = c(linearFit$par[1], linearFit$par[2], exp(linearFit$par[3]), nllPVal_lin)

# Linear regression via lm()
linear.mod = lm(growthSugar ~ sugar, data=sugarData)

# All parameters from lm()
lmParams = c(linear.mod$coefficients, sigma(linear.mod), extractPVal(linear.mod))

# Comparison of the parameters
comparison_lin = cbind(nllParams_lin, lmParams)
dimnames(comparison_lin)[[1]]=c("Beta 0", "Beta 1", "Residual Std Error", "P Value")
comparison_lin

# 95% confidence interval of the parameters b0 and b1
CI_LM = confint(linear.mod)
dimnames(CI_LM)[[1]]=c("Beta 0", "Beta 1")
CI_LM

# Plot the data and fit
ggplot(sugarData,aes(sugarData$sugar,sugarData$growthSugar))+
  geom_point()+
  geom_abline(slope=linear.mod$coefficients[2],intercept=linear.mod$coefficients[1])+
  geom_abline(slope=0,intercept=simpleFit$par[1],color="red")+
  theme_classic()+
  xlab("Sugar Concentration")+
  ylab("Growth of E.Coli")+
  geom_smooth(method="lm")

##########################################################
##########################################################
############### Part III: Power Analysis #################
##########################################################
##########################################################

Nsim=100 # Number of simulations to run
beta0=10 # The intercept of the line
beta1=0.4 # The slope of the line
sigmaVals = c(1,2,4,6,8,12,16,24) # Standard Deviations for the error term

N=24 # Number of experimental units
maxX = 50 # Maximum value of experimental data
minX = 0 # Minimum value of experimental data
nLevelsVals = c(2,4,8) # Number of levels in ANOVA

paramComparison = NULL
paramComparison_nll = NULL

# Loop over all sigma values
for (nLevels in nLevelsVals){
  for (sigma in sigmaVals){

    # Linear Regression
    X=seq(minX,maxX,(maxX+1)/N) # Sequenced X values

    # Declare Linear Regression Matrices
    Y=matrix(rep(N*Nsim),nrow=N,ncol=Nsim) # Create our Y matrix
    coeff_matrix_lm=matrix(rep(2*Nsim),nrow=2,ncol=Nsim) # Create our coefficent matrix
    sigma_matrix_lm=matrix(rep(Nsim),nrow=1,ncol=Nsim) # Create our sigma matrix
    p_matrix_lm=matrix(rep(Nsim),nrow=1,ncol=Nsim) # Create our p value matrix

    # Declare Anova matrices
    coeff_matrix_aov=matrix(rep(nLevels*Nsim),nrow=nLevels,ncol=Nsim) # Create our coefficent matrix for aov
    sigma_matrix_aov=matrix(rep(Nsim),nrow=1,ncol=Nsim) # Create our sigma matrix for aov
    p_matrix_aov=matrix(rep(Nsim),nrow=1,ncol=Nsim) # Create our p value matrix for aov

    # running the simulations using lm() and aov()
    for (i in 1:Nsim){
      epsilon=rnorm(N,mean=0,sd=sigma) # Normalized error
      Y[,i]=beta0+beta1*X+epsilon  # Generate data points
      
      mod=lm(Y[,i]~X) # Fit the data using lm()
      coeff_matrix_lm[,i]=mod$coefficients # Extract coefficients
      sigma_matrix_lm[i]=summary(mod)$sigma # Extract sigmas
      p_matrix_lm[i]=extractPVal(mod) # Extract P Value
  
      anovaResults=data.frame(levels=rep(seq(0,nLevels-1,1),each=N/nLevels),YVal=Y[,i])
      aovMod = aov(anovaResults$YVal~as.factor(anovaResults$levels))
      coeff_matrix_aov[,i]=aovMod$coefficients
      sigma_matrix_aov[i]=sigma(aovMod)
      p_matrix_aov[i]=summary(aovMod)[[1]][["Pr(>F)"]][1]
    }

    # Calculate the averages of the 100 simulations
    averageb0_reg = mean(coeff_matrix_lm[1])
    averageb1_reg = mean(coeff_matrix_lm[2])
    averagep_reg = mean(p_matrix_lm)
    significantP_reg = length(which(p_matrix_lm<0.05))

    averageb0_aov = mean(coeff_matrix_aov[1])
    averageb1_aov = mean(coeff_matrix_aov[2])
    averagep_aov = mean(p_matrix_aov)
    significantP_aov = length(which(p_matrix_aov<0.05))
    
    # Compare the two parameters and P Values
    linearParameters = c(averageb0_reg, averageb1_reg,averagep_reg,significantP_reg)
    anovaParameters = c(averageb0_aov, averageb1_aov,averagep_aov,significantP_aov)
    comparison_III = cbind(linearParameters,anovaParameters)
    dimnames(comparison_III)[[1]]=c("Beta 0", "Beta 1", "P Value","Num Significant P")
    paramComparison = cbind(paramComparison, sigma, comparison_III)
    
    ##########################################################################################################
    
    # Declare Linear Regression Matrices
    Y_nll=matrix(rep(N*Nsim),nrow=N,ncol=Nsim) # Create our Y matrix
    coeff_matrix_lm_nll=matrix(rep(2*Nsim),nrow=2,ncol=Nsim) # Create our coefficent matrix
    sigma_matrix_lm_nll=matrix(rep(Nsim),nrow=1,ncol=Nsim) # Create our sigma matrix
    p_matrix_lm_nll=matrix(rep(Nsim),nrow=1,ncol=Nsim) # Create our p value matrix
    
    # Declare Anova matrices
    coeff_matrix_aov_nll=matrix(rep(nLevels*Nsim),nrow=nLevels,ncol=Nsim) # Create our coefficent matrix for aov
    sigma_matrix_aov_nll=matrix(rep(Nsim),nrow=1,ncol=Nsim) # Create our sigma matrix for aov
    p_matrix_aov_nll=matrix(rep(Nsim),nrow=1,ncol=Nsim) # Create our p value matrix for aov
    
    # running the simulations using custom nll() functions
    for (i in 1:Nsim){
      epsilon=rnorm(N,mean=0,sd=sigma) # Normalized error
      Y_nll[,i]=beta0+beta1*X+epsilon  # Generate data points
  
      mod_nll=linearNLL(X,Y[,i]) # Use nll to obtain parameter estimates
      coeff_matrix_lm_nll[,i]=mod_nll$coefficients # Extract coefficients
      sigma_matrix_lm_nll[i]=mod_nll$sigma # Extract sigmas
      p_matrix_lm_nll[i]=mod_nll$pValue # Extract P Value
    
      anovaResults=data.frame(levels=rep(seq(0,nLevels-1,1),each=N/nLevels),YVal=Y[,i])
      aovMod_nll = customNll(anovaResults$levels, anovaResults$YVal, nLevels)
      coeff_matrix_aov_nll[,i]=aovMod_nll$coefficients
      sigma_matrix_aov_nll[i]=aovMod_nll$sigma
      p_matrix_aov_nll[i]=aovMod_nll$pValue
    }
    
    # Calculate the averages of the 100 simulations
    averageb0_reg_nll = mean(coeff_matrix_lm_nll[1])
    averageb1_reg_nll = mean(coeff_matrix_lm_nll[2])
    averagep_reg_nll = mean(p_matrix_lm_nll)
    significantP_reg_nll = length(which(p_matrix_lm_nll<0.05))
    
    averageb0_aov_nll = mean(coeff_matrix_aov_nll[1])
    averageb1_aov_nll = mean(coeff_matrix_aov_nll[2])
    averagep_aov_nll = mean(p_matrix_aov_nll)
    significantP_aov_nll = length(which(p_matrix_aov_nll<0.05))
    
    # Compare the two parameters and P Values
    linearParameters_nll = c(averageb0_reg_nll, averageb1_reg_nll,averagep_reg_nll,significantP_reg_nll)
    anovaParameters_nll = c(averageb0_aov_nll, averageb1_aov_nll,averagep_aov_nll,significantP_aov_nll)
    comparison_III_nll = cbind(linearParameters_nll,anovaParameters_nll)
    dimnames(comparison_III_nll)[[1]]=c("Beta 0", "Beta 1", "P Value","Num Significant P")
    paramComparison_nll = cbind(paramComparison_nll, sigma, comparison_III_nll)
  }
}

twoLevelANOVA = cbind(paramComparison[,1:24])
fourLevelANOVA = cbind(paramComparison[,25:48])
eightLevelANOA = cbind(paramComparison[,49:72])

twoLevelANOVA_nll = cbind(paramComparison_nll[,1:24])
fourLevelANOVA_nll = cbind(paramComparison_nll[,25:48])
eightLevelANOA_nll = cbind(paramComparison_nll[,49:72])

write.table(twoLevelANOVA_nll, "clipboard", sep="\t", row.names=TRUE)

