source("ExampleSimulationFunctionsAndParameters.R")
set.seed(123)
CompleteSample <- SimulationProcess(Start = c(0, 0), Npoints = 10^5, Delta = 0.005)
write.table(x = CompleteSample, file = "SimulatedSample.txt", 
            row.names = F, col.names = T)
