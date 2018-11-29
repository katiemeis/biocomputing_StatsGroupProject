anovaNLL = function(x,y,nLevels){
  # log likelihood function
  
  nll=function(p,x,y){
    B = NULL
    expected = NULL
    for (i in 1:nLevels){
      B[i] = p[i]  
      expected = expected + B[i]*XX[i]
    }
    sigma=exp(p[nLevels+1])

    expected=B0+x*B1
    nll=-sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
    return(nll) 
  }  
  
  # Initial guess for optimization
  initialGuess_lin=c(y[which(x==0)],
                     (mean(y)/mean(x)),
                     1)
  
  # Optimize the nll function
  linearFit=optim(par=initialGuess_lin,fn=nllLinear,x=x,y=y)
  
  # Simple model for nll function
  nllSimple_lin<-function(p,x,y){
    B0=p[1]
    sigma=exp(p[2])
    expected=B0
    nll=-sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
    return(nll) 
  }
  
  # Initial guess for simple model
  simpleGuess_lin=c(y[which(x==0)],1)
  
  # Optimize the simple model
  simpleFit=optim(par=simpleGuess_lin,fn=nllSimple_lin,x=x,y=y)
  
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
  nllReturns = attributes(nllReturns)
  return(nllReturns)
}








