extractPVal <- function (linearFit) {
  fStat = summary(linearFit)$fstatistic
  pVal = pf(fStat[1],fStat[2],fStat[3],lower.tail=F)
  attributes(pVal) = NULL
  return(pVal)
}