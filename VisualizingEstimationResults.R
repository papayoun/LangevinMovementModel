rm(list = ls())
load(file = "Results_OzakiEulerComparison_Analytics.RData")

TimeSteps <- do.call(c, lapply(AllResults, function(x) diff(x$Data[, "Time"]))) 
barplot(table(round(TimeSteps,2)))
EulEstimates <- do.call(rbind.data.frame, lapply(AllResults, function(x){
  res = x$fitE$par
  data.frame(Beta1 = res[1], Beta2 = res[2], Alpha = res[3], method = "euler")
}))
OzaEstimates <- do.call(rbind.data.frame, lapply(AllResults, function(x){
  res = x$fitO$par
  data.frame(Beta1 = res[1], Beta2 = res[2], Alpha = res[3], method = "ozaki")
}))

Estimates <- rbind.data.frame(OzaEstimates, EulEstimates)
par(mfrow = c(1, 3))
boxplot(Estimates[, "Beta1"] ~ Estimates[, "method"])
abline(h = -1)
boxplot(Estimates[, "Beta2"] ~ Estimates[, "method"])
abline(h = 0.5)
boxplot(Estimates[, "Alpha"] ~ Estimates[, "method"])
abline(h = 0.05)
