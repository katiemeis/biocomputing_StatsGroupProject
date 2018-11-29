customNll = function(x,y,nLevels){
  # custom log likelihood function
  
  nll=function(p,x,y){
    B = matrix(0, 1, nLevels)
    expected = 0
    XX = matrix(0, 24, nLevels)
    
    for (j in 1:nLevels){
      for (i in ((24*j/nLevels)-(nLevels+1)):(24*j/nLevels)){
        XX[i,j] = 1
      }
    }
    
    #B0 to BN
    for (i in 2:nLevels-1){
      B[i] = p[i]  
      expected = expected + B[i]*XX[,i+1]
    }
    expected = expected + B[1]
    
    sigma=exp(p[nLevels+1])
    
    nll=-sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
    return(nll) 
  }  
  
  # Initial guess for optimization
  initialGuess = matrix(0, 1, nLevels+1)
  initialGuess[1] = mean(y[((24/nLevels)-((24/nLevels)-1)):(24/nLevels)])
  initialGuess[nLevels+1] = 1
  for (i in 2:nLevels){
    initialGuess[i] = mean(y[((24*i/nLevels)-((24/nLevels)-1)):(i*24/nLevels)])-initialGuess[1]
  }
  
  # Optimize the nll function
  linearFit=optim(par=initialGuess,fn=nll,x=x,y=y)
  
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
  teststat_lin=2*(simpleFit$value-linearFit$value)
  
  # Compute the degrees of freedom
  df_lin=length(linearFit$par)-length(simpleFit$par)
  
  # P value from nll function
  pValue=1-pchisq(teststat_lin,df_lin)
  
  # Coefficients from nll function
  coefficients = linearFit$par[1:nLevels]
  
  # Sigma from nll function
  sigma = exp(linearFit$par[nLevels+1])
  
  # Return the pValue, coefficients, and sigma
  nllReturns = NULL
  attributes(nllReturns)$pValue = pValue
  attributes(nllReturns)$coefficients = coefficients
  attributes(nllReturns)$sigma = sigma
  nllReturns = attributes(nllReturns)
  return(nllReturns)
}








