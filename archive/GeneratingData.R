source("archive/SimulationFunctions.R")
source("archive/CovariateParameters.R")# contains ALpha, Beta1 and Beta2
set.seed(123)
A <- 0; B1 <- -1; B2 <- 0.5
CompleteSample <- SimulationProcess(Start = c(0, 0), Npoints = 5 * 10^5, 
                                    Delta = 0.001, A, B1, B2)
write.table(x = CompleteSample, file = "archive/SimulatedSample.txt", 
            row.names = F, col.names = T)
