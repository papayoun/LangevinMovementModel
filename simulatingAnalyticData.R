rm(list = ls())
source("archive/SimulationFunctions.R")
source("archive/CovariateParameters.R")
source("LinearAlgebraFunctions.R")
source("optionsParsingScript.R")
SAVEPATH <- "simulatedData/analyticCase/"
A <- 0.05; B1 <- -1; B2 <- 0.5
set.seed(seed)
CompleteSample <- SimulationProcess(Start = c(0, 0), Npoints = 5 * 10^5 + 1, 
                                    Delta = 0.01, A, B1, B2, 1, method = "Ozaki")
write.table(CompleteSample, file = paste0(SAVEPATH, "simAnalyticData_seed", seed, ".txt"), col.names = T, row.names = F, sep = ";")

CompleteSample <- SimulationProcess(Start = c(0, 0), Npoints = 5 * 10^5 + 1, 
                                    Delta = 0.01, A, B1, B2, 2, method = "Ozaki")