# Custom function to optimize parameters using
# negative log likelihood for anova and linear regression
# x := The x data we wish to observe a relationship for
# y := The y data we wish to observe a relationship for
# nLevels := number of levels to be used; Default is 2 for simple linear regression
# nExpUnits := Number of data points in the data set
# anova := Will we be using linear regression or anvoa? Default is false

superNll = function(x,y,nLevels=2,nExpUnits, anova=FALSE){

# If we are using anova, use a more complex matrice solution for nll
if (anova == TRUE){
  # Define the function to optimize
  nll=function(p,x,y){
    # Declare empty matrices to hold data
    B = matrix(0, 1, nLevels)
    expected = 0
    XX = matrix(0, nExpUnits, nLevels)
    
    # Fill in a matrix of size [nExpUnits x nLevels]
    # This will hold binary values for whether or not
    # a specific beta parameter is used
    for (j in 1:nLevels){
      for (i in ((nExpUnits*j/nLevels)-((nExpUnits/nLevels)-1)):(nExpUnits*j/nLevels)){
        XX[i,j] = 1
      }
    }
    
    #Define the values for the beta parameters
    B[1] = p[1]
    for (i in 2:nLevels){
      B[i] = p[i]  
      expected = expected + B[i]*XX[,i]
    }
    expected = expected + B[1]
    
    # Define the value for the residual standard error parameters
    sigma=exp(p[nLevels+1])
    
    # Compute and return the negative log likelihood
    nll=-sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
    return(nll) 
  }  
  
  # Initial guess for optimization
  # Declare empty matrix to be filled in
  initialGuess = matrix(0, 1, nLevels+1)
  # First guess is the mean of the first subset of data
  initialGuess[1] = mean(y[((nExpUnits/nLevels)-((nExpUnits/nLevels)-1)):(nExpUnits/nLevels)])
  # Initial guess for residual standard error is 1
  initialGuess[nLevels+1] = 1
  # Subsequent guesses after the first are differences in
  # means of the current subset of the data compared to the first subset
  for (i in 2:nLevels){
    initialGuess[i] = mean(y[((nExpUnits*i/nLevels)-((nExpUnits/nLevels)-1)):(i*nExpUnits/nLevels)])-initialGuess[1]
  }
  
  # Optimize the nll function
  myNll=optim(par=initialGuess,fn=nll,x=x,y=y, hessian = TRUE)
  
  # Simple model for nll function
  nllSimple_lin<-function(p,x,y){
    B0=p[1]
    sigma=exp(p[2])
    expected=B0
    nll=-sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
    return(nll) 
  }
  
  # Initial guess for simple model
  simpleGuess=c(1,1)
  
  # Optimize the simple model
  simpleFit=optim(par=simpleGuess,fn=nllSimple_lin,x=x,y=y)
  
  # Compute the LR statistic
  teststat_lin=2*(simpleFit$value-myNll$value)
  
  # Compute the degrees of freedom
  df_lin=length(myNll$par)-length(simpleFit$par)
  
  # P value from nll function
  pValue=1-pchisq(teststat_lin,df_lin)
  
  # Coefficients from nll function
  coefficients = myNll$par[1:nLevels]
  
  # Sigma from nll function
  sigma = exp(myNll$par[nLevels+1])
  
  # Return the pValue, coefficients, and sigma
  nllReturns = NULL
  attributes(nllReturns)$pValue = pValue
  attributes(nllReturns)$coefficients = coefficients
  attributes(nllReturns)$sigma = sigma
  attributes(nllReturns)$hessian = myNll$hessian
  nllReturns = attributes(nllReturns)
  return(nllReturns)
}
  
# If we are not using anova, use simple linear regression for nll
if (anova == FALSE){
  # Linear log likelihood function
  nllLinear=function(p,x,y){
    B0=p[1] 
    B1=p[2] 
    sigma=exp(p[3])
    expected=B0+x*B1
    nll=-sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
    return(nll) 
  }  
  # Initial guess for optimization
  # First guess tries to approximate the intercept from slope
  # Second guess uses end points of the data to approximate the slope
  # Guess of the residual standard error is 1
  slope = (y[nExpUnits]-y[1]) / (x[nExpUnits]-x[1])
  firstGuess = y[1]-slope*x[1]
  initialGuess_lin=c(firstGuess,
                     slope,
                     1)
  # Optimize the nll function
  linearFit=optim(par=initialGuess_lin,fn=nllLinear,x=x,y=y,hessian = TRUE)
  
  # Simple model for nll function
  nllSimple_lin<-function(p,x,y){
    B0=p[1]
    sigma=exp(p[2])
    expected=B0
    nll=-sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
    return(nll) 
  }
  
  # Initial guess for simple model
  simpleGuess_lin=c(firstGuess,1)
  
  # Optimize the simple model
  simpleFit=optim(par=simpleGuess_lin,fn=nllSimple_lin,x=x,y=y, hessian = TRUE)
  
  # Compute the LR statistic
  teststat_lin=2*(simpleFit$value-linearFit$value)
  
  # Compute the degrees of freedom
  df_lin=length(linearFit$par)-length(simpleFit$par)
  
  # P value from nll function
  pValue=1-pchisq(teststat_lin,df_lin)
  
  # Coefficients from nll function
  coefficients = linearFit$par[1:2]
  
  # Sigma from nll function
  sigma = exp(linearFit$par[3])
  
  # Return the pValue, coefficients, and sigma
  nllReturns = NULL
  attributes(nllReturns)$pValue = pValue
  attributes(nllReturns)$coefficients = coefficients
  attributes(nllReturns)$sigma = sigma
  attributes(nllReturns)$hessian = linearFit$hessian
  nllReturns = attributes(nllReturns)
  return(nllReturns)
}  
  
}








