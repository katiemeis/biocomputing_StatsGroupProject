setwd("~/Documents/Biocomp/biocomputing_StatsGroupProject/")
rm(list=ls())

library(ggplot2)

Sugar<-read.csv("sugar.csv")

ggplot(Sugar, aes(x=sugar,y=growth))+
  geom_point()

nllLinear<-function(p,x,y){
  B0=p[1] 
  B1=p[2] 
  sigma=exp(p[3])
  expected=B0+x*B1
  nll=-sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
  return(nll) 
}

initialGuess=c(1,1,1)
linearFit=optim(par=initialGuess,fn=nllLinear,x=Sugar$sugar,y=Sugar$growth)
print(linearFit)

nllSimple<-function(p,x,y){
  B0=p[1]
  sigma=exp(p[2])
  expected=B0+x
  nll=-sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
  return(nll) 
}

firstGuess=c(1,1)
simpleFit=optim(par=firstGuess,fn=nllSimple,x=Sugar$sugar,y=Sugar$growth)
print(simpleFit)

ggplot(Sugar,aes(sugar,growth))+
  geom_point()+
  geom_abline(slope=1.73,intercept=-0.85)+
  geom_abline(slope=0,intercept=3.33,color="red")

teststat=2*(simpleFit$value-linearFit$value)
df=length(linearFit$par)-length(simpleFit$par)
1-pchisq(teststat,df)

ggplot(Sugar,aes(sugar,growth))+
  geom_point()+
  geom_abline(slope=1.73,intercept=-0.85)
