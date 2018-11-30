##########################################################
##########################################################
################ BIOS 60318 Final Project ################
##########################################################
##########################################################
rm(list=ls()) # Clear workspace
if(!is.null(dev.list())) dev.off() # Clear plots
par(mfrow=c(1,1)) # Setup plot parameters
par(ps = 12, font.lab = 1) # Plot parameters
set.seed(1) # set random seed generator for reproducibility

# If you don't have the packages necessary, let's make sure you do
list.of.packages <- c("ggplot2", "CorReg")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Let's load relevant libraries
library(ggplot2) # Allow plotting capabilities
library(CorReg) # Allows use of modified BoxPlot()
source('~/extractPVal.R', encoding = 'UTF-8') # Custom function to extract a Pvalue from an lm() model
source('~/superNll.R', encoding = 'UTF-8') # Custom function to run MLLE on a data set in either regression or ANOVA format
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
antibioticsData$trt = as.numeric(relevel(antibioticsData$trt, ref=4))

# Use a custom nll function to obtain parameters and pValue for log likelihood ratio test
antibioticNll = superNll(antibioticsData$trt, antibioticsData$growth, 4, length(antibioticsData$growth), anova = TRUE)

# Fit using aov() that passes each group through lm() as a comparison
aov.fit = aov(growth ~ as.factor(trt), data=antibioticsData)

# PValue from aov() using F statistic
aovPVal=summary(aov.fit)[[1]][["Pr(>F)"]][1]

# Parameters found from nll
antibioticParam_nll = c(antibioticNll$coefficients, antibioticNll$sigma, antibioticNll$pValue)

# Parameters found from lm() and aov()
antibioticParam_aov = c(aov.fit$coefficients, sigma(aov.fit), aovPVal)

# Comparison of nll versus lm()
comparison_antibiotics = cbind(antibioticParam_nll, antibioticParam_aov)
dimnames(comparison_antibiotics)[[1]]=c("Control", "Treatment 1", "Treatment 2", "Treatment 3", "Residual Std Error", "P Value")
comparison_antibiotics
#write.table(comparison_antibiotics, "clipboard", sep="\t", row.names=FALSE) # Transport comparison to excel

# 95% Confidence interval of the parameters 
antibiotic_CI = confint(aov.fit, level = 0.95)
dimnames(antibiotic_CI)[[1]]=c("Control", "Treatment 1", "Treatment 2", "Treatment 3")
antibiotic_CI

# Boxplot of the data comparing control versus treatment with 95% CI of the means for each group in red
BoxPlot(antibioticsData$growth, 
        as.factor(antibioticsData$trt), 
        AnoVa = TRUE, ylab="Growth of Bacteria", 
        names=c("Control", "Treatment 1", "Treatment 2", "Treatment 3"),
        verbose=FALSE)

##########################################################
##########################################################
############### Part II: Linear Regression ###############
##########################################################
##########################################################

sugarNll = superNll(sugarData$sugar,sugarData$growthSugar,2,length(sugarData$sugar),anova=FALSE)

# Combine all parameters for nll
sugarParam_nll = c(sugarNll$coefficients, sugarNll$sigma, sugarNll$pValue)

# Linear regression via lm()
linear.mod = lm(growthSugar ~ sugar, data=sugarData)

# All parameters from lm()
sugarParam_lm = c(linear.mod$coefficients, sigma(linear.mod), extractPVal(linear.mod))

# Comparison of the parameters
comparison_sugar = cbind(sugarParam_nll, sugarParam_lm)
dimnames(comparison_sugar)[[1]]=c("Beta 0", "Beta 1", "Residual Std Error", "P Value")
comparison_sugar

# 95% confidence interval of the parameters b0 and b1
CI_LM = confint(linear.mod)
dimnames(CI_LM)[[1]]=c("Beta 0", "Beta 1")
CI_LM

# Plot the data and fit
ggplot(sugarData,aes(sugarData$sugar,sugarData$growthSugar))+
  geom_point()+
  geom_abline(slope=linear.mod$coefficients[2],intercept=linear.mod$coefficients[1])+
  geom_abline(slope=0,intercept=3.33,color="red")+
  theme_classic()+
  xlab("Sugar Concentration")+
  ylab("Growth of E.Coli")+
  geom_smooth(method="lm")

##########################################################
##########################################################
############### Part III: Power Analysis #################
##########################################################
##########################################################

Nsim=10000 # Number of simulations to run
beta0=10 # The intercept of the line
beta1=0.4 # The slope of the line
sigmaVals = c(1,2,4,6,8,12,16,24) # Standard Deviations for the error term

N=24 # Number of experimental units
maxX = 50 # Maximum value of experimental data
minX = 0 # Minimum value of experimental data
nLevelsVals = c(2,4,8) # Number of levels in ANOVA

# Define NULL vectors to hold data later in the code
paramComparison = NULL
paramComparison_nll = NULL
pMatrixAll = NULL
pMatrixRun = matrix(0, Nsim, 6)

# Loop over all levels and all sigmas
for (nLevels in nLevelsVals){
  for (sigma in sigmaVals){

    # Linear Regression
    X=seq(minX,maxX,(maxX+1)/N) # Sequenced X values

    # Declare Linear Regression Matrices for lm()
    Y=matrix(rep(N*Nsim),nrow=N,ncol=Nsim) # Create our Y matrix
    coeff_matrix_lm=matrix(rep(2*Nsim),nrow=2,ncol=Nsim) # Create our coefficent matrix
    sigma_matrix_lm=matrix(rep(Nsim),nrow=1,ncol=Nsim) # Create our sigma matrix
    p_matrix_lm=matrix(rep(Nsim),nrow=1,ncol=Nsim) # Create our p value matrix

    # Declare Anova matrices for aov()
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

    # Calculate the averages of the simulations
    averageb0_reg = mean(coeff_matrix_lm[1])
    averageb1_reg = mean(coeff_matrix_lm[2])
    averagep_reg = mean(p_matrix_lm)
    significantP_reg = length(which(p_matrix_lm<0.05))

    averageb0_aov = mean(coeff_matrix_aov[1])
    averageb1_aov = mean(coeff_matrix_aov[2])
    averagep_aov = mean(p_matrix_aov)
    significantP_aov = length(which(p_matrix_aov<0.05))
    
    # Compare the two parameters and P Values between lm() and aov()
    linearParameters = c(averageb0_reg, averageb1_reg,averagep_reg,significantP_reg)
    anovaParameters = c(averageb0_aov, averageb1_aov,averagep_aov,significantP_aov)
    comparison_III = cbind(linearParameters,anovaParameters)
    dimnames(comparison_III)[[1]]=c("Beta 0", "Beta 1", "P Value","Num Significant P")
    paramComparison = cbind(paramComparison, sigma, comparison_III)
    
    ##########################################################################################################
    
    # Declare Linear Regression Matrices for custom nll function
    Y_nll=matrix(rep(N*Nsim),nrow=N,ncol=Nsim) # Create our Y matrix
    coeff_matrix_lm_nll=matrix(rep(2*Nsim),nrow=2,ncol=Nsim) # Create our coefficent matrix
    sigma_matrix_lm_nll=matrix(rep(Nsim),nrow=1,ncol=Nsim) # Create our sigma matrix
    p_matrix_lm_nll=matrix(rep(Nsim),nrow=1,ncol=Nsim) # Create our p value matrix
    
    # Declare Anova matrices for custom nll function
    coeff_matrix_aov_nll=matrix(rep(nLevels*Nsim),nrow=nLevels,ncol=Nsim) # Create our coefficent matrix for aov
    sigma_matrix_aov_nll=matrix(rep(Nsim),nrow=1,ncol=Nsim) # Create our sigma matrix for aov
    p_matrix_aov_nll=matrix(rep(Nsim),nrow=1,ncol=Nsim) # Create our p value matrix for aov
    
    # running the simulations using custom nll() functions
    for (i in 1:Nsim){
      epsilon=rnorm(N,mean=0,sd=sigma) # Normalized error
      Y_nll[,i]=beta0+beta1*X+epsilon  # Generate data points
  
      mod_nll=superNll(X,Y[,i],nLevels = 2,nExpUnits = N,anova = FALSE) # Use nll to obtain parameter estimates for linear regression
      coeff_matrix_lm_nll[,i]=mod_nll$coefficients # Extract coefficients
      sigma_matrix_lm_nll[i]=mod_nll$sigma # Extract sigmas
      p_matrix_lm_nll[i]=mod_nll$pValue # Extract P Value
    
      anovaResults=data.frame(levels=rep(seq(0,nLevels-1,1),each=N/nLevels),YVal=Y[,i]) # set up anova data frame and bin data
      aovMod_nll = superNll(anovaResults$levels, anovaResults$YVal, nLevels, N,anova = TRUE) # use nll to obtain parameter estimates from anova regression
      coeff_matrix_aov_nll[,i]=aovMod_nll$coefficients # Extract coefficients
      sigma_matrix_aov_nll[i]=aovMod_nll$sigma # Extract sigmas
      p_matrix_aov_nll[i]=aovMod_nll$pValue # Extract p values
    }
    
    # Calculate the averages of the simulations
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
    
    # Combine all P Values for each MLLE and LSE run
    pMatrixRun[,1] = p_matrix_lm
    pMatrixRun[,2] = p_matrix_aov
    pMatrixRun[,3] = p_matrix_lm_nll
    pMatrixRun[,4] = p_matrix_aov_nll
    pMatrixRun[,5] = sigma
    pMatrixRun[,6] = nLevels
    pMatrixAll = cbind(pMatrixAll, pMatrixRun)
  }
}

# Extract the two level comparison data for lm() vs. aov()
twoLevelANOVA = cbind(paramComparison[,1:24])
# Extract the four level comparison data for lm() vs. aov()
fourLevelANOVA = cbind(paramComparison[,25:48])
# Extract the eight level comparison data for lm() vs. aov()
eightLevelANOVA = cbind(paramComparison[,49:72])

# Extract the two level comparison data for linear regression nll vs. anova nll
twoLevelANOVA_nll = cbind(paramComparison_nll[,1:24])
# Extract the four level comparison data for linear regression nll vs. anova nll
fourLevelANOVA_nll = cbind(paramComparison_nll[,25:48])
# Extract the eight level comparison data for linear regression nll vs. anova nll
eightLevelANOVA_nll = cbind(paramComparison_nll[,49:72])

# Plot the histogram of the p-values for lm() vs. aov() for 8 levels and sigma = 24 (MSE)
# LSE For linear regression p distribution
ggplot(data.frame(pVal=t(p_matrix_lm)), aes(x=pVal))+
  geom_histogram(binwidth = 0.1, color="white",fill="blue")+
  theme_classic()+
  xlab("p-Values")+
  ylab("Frequency")+
  ggtitle(label="p-Value distribution of lm()")+
  theme(plot.title=element_text(hjust=0.5))+
  geom_vline(aes(xintercept=mean(pVal)),color="red",size=1.2)

# LSE For ANOVA p distribution
ggplot(data.frame(pVal=t(p_matrix_aov)), aes(x=pVal))+
  geom_histogram(binwidth = 0.1, color="white",fill="blue")+
  theme_classic()+
  xlab("p-Values")+
  ylab("Frequency")+
  ggtitle(label="p-Value distribution of aov()")+
  theme(plot.title=element_text(hjust=0.5))+
  geom_vline(aes(xintercept=mean(pVal)),color="red",size=1.2)

# Plot the histogram of the p-values for linear nll vs. anova nll for 8 levels and sigma = 24 (MLLE)
# MLLE For linear regression p distribution
ggplot(data.frame(pVal=t(p_matrix_lm_nll)), aes(x=pVal))+
  geom_histogram(binwidth = 0.1, color="white",fill="blue")+
  theme_classic()+
  xlab("p-Values")+
  ylab("Frequency")+
  ggtitle(label="p-Value distribution of linear regresssion nll()")+
  theme(plot.title=element_text(hjust=0.5))+
  geom_vline(aes(xintercept=mean(pVal)),color="red",size=1.2)

# MLLE For ANOVA p distribution
ggplot(data.frame(pVal=t(p_matrix_aov_nll)), aes(x=pVal))+
  geom_histogram(binwidth = 0.1, color="white",fill="blue")+
  theme_classic()+
  xlab("p-Values")+
  ylab("Frequency")+
  ggtitle(label="p-Value distribution of anova nll()")+
  theme(plot.title=element_text(hjust=0.5))+
  geom_vline(aes(xintercept=mean(pVal)),color="red",size=1.2)

# Plot the histogram of the p-values for lm() vs. aov() for 2 levels and sigma = 1 (LSE)
# LSE For linear regression p distribution
ggplot(data.frame(pVal=-log10(pMatrixAll[,1])), aes(x=pVal))+
  geom_histogram(binwidth = 0.1, color="white",fill="blue")+
  theme_classic()+
  xlab("-log10(p-Values)")+
  ylab("Frequency")+
  theme(plot.title=element_text(hjust=0.5))+
  geom_vline(aes(xintercept=mean(pVal)),color="red",size=1.2)+
  xlim(min(-log10(pMatrixAll[,1])),max(-log10(pMatrixAll[,1])))

# LSE For ANOVA p distribution
ggplot(data.frame(pVal=-log10(pMatrixAll[,2])), aes(x=pVal))+
  geom_histogram(binwidth = 0.05, color="white",fill="blue")+
  theme_classic()+
  xlab("-log10(p-Values)")+
  ylab("Frequency")+
  theme(plot.title=element_text(hjust=0.5))+
  geom_vline(aes(xintercept=mean(pVal)),color="red",size=1.2)+
  xlim(min(-log10(pMatrixAll[,2])),max(-log10(pMatrixAll[,2])))

# Plot the histogram of the p-values for custom nll for Linear regression and ANOVA for 2 levels and sigma = 1 (MLLE)
# MLLE For linear regression p distribution
# NOTE: 1e24 was added to all values due to the values being below the machine computable epsilon value
ggplot(data.frame(pVal=-log10(pMatrixAll[,3]+1e-24)), aes(x=pVal))+
  geom_histogram(binwidth = 0.1, color="white",fill="blue")+
  theme_classic()+
  xlab("-log10(p-Values)")+
  ylab("Frequency")+
  theme(plot.title=element_text(hjust=0.5))+
  geom_vline(aes(xintercept=mean(pVal)),color="red",size=1.2)+
  xlim(23,25)

# MLLE For ANOVA p distribution
ggplot(data.frame(pVal=-log10(pMatrixAll[,4])), aes(x=pVal))+
  geom_histogram(binwidth = 0.05, color="white",fill="blue")+
  theme_classic()+
  xlab("-log10(p-Values)")+
  ylab("Frequency")+
  theme(plot.title=element_text(hjust=0.5))+
  geom_vline(aes(xintercept=mean(pVal)),color="red",size=1.2)+
  xlim(min(-log10(pMatrixAll[,4])),max(-log10(pMatrixAll[,4])))

#write.table(eightLevelANOVA, "clipboard", sep="\t", row.names=TRUE) # Copy data to clipboard
#write.table(eightLevelANOVA_nll, "clipboard", sep="\t", row.names=TRUE) # copy data to clipboard
