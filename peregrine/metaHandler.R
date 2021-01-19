


# Cleaning up, loading resources and setting a random seed
rm (list = ls( )) 

#setwd(dirname(parent.frame(2)$ofile))
#setwd("D:/Dropbox/Progetto/w/calibration study/Scripts")
setwd("./")

#require(parallel)
source("util.R")
source("script NI&PA.R")
source("geoAbm.R")

#seed  <-  sample.int(999999999, size = 1)
#seed  <-  123456789
#set.seed (seed)


# Setup parallel execution
#nCores <- detectCores()
#nCores <- detectCores() - 2 # <<<  This line is to be commented 
                            # when running simulation batteries.
#cluster <- makeCluster(nCores)

# ====== Calling parallel runs ==================================

#clusterExport(cluster) # Add global variables as arguments.

# Initiate parallel runs
#parLapply(cl, 
#          2:4, 
#          function(exponent) 
#          base^exponent
#)

# Terminate
#stopCluster(cluster)

# ==============================================================


# This makes r compile functions on first execution.
library(compiler)
enableJIT(1)

# ID of the simulation battery
batteryID <- sample.int(999999999, 1)

# We create a dataframe, where each row is a unique parameter configuration
# that we want to run
#parameterSpace <- expand.grid(c(9),c("uniform"),c(0.6, 0.9),c(1))
parameterSpace <- expand.grid(
  c(1:12),                          # wijk
  c("uniform", "beta", "groupBias"),# initialOpinionDistribution
  c(0.6, 0.9),                      # H
  c(1, 2, 3)                        # distance decay
)
parameterSpace <- cbind(
  parameterSpace,
  200,      # timeMax
  FALSE,   # printOpinionHistogram
  TRUE     # exportOutput
)
names(parameterSpace) <- c(
  "wijk",
  "initialOpinionDistribution",
  "H",
  "distanceDecay",
  "timeMax",
  "printOpinionHistogram",
  "exportOutput"
)


# We specify how many times we want to run a parameter configuration. Each run
# will have a different random seed.
nRunsPerCondition <- 1


# Next, we prepare an empty dataframe where to store the outcome of every
# simulation run. The variable "indexParameters" is the key to join simulation
# results and their parameters. Similarly, the variable "seed" will allow to
# find the world/agents dataframe associated with these outcomes.
outcomeVariables <- c(
  "seed",
  "indexParameters",
  "steps",
  "polarizationIndex",
  "meanOpinionGlobal",
  "absOpGlobal",
  "varOpinionGlobal",
  "meanOpinionG1",
  "absOpG1",
  "varOpinionG1",
  "meanOpinionG2",
  "absOpG2",
  "varOpinionG2",
  "opClustering1",
  "opClustering2",
  "opClustering3",
  "opClusteringA1",
  "opClusteringA2",
  "opClusteringA3",
  "opAlignment1",
  "opAlignment2",
  "opAlignment3"
)

# For every simulation run, we get as output a list containing [[1]] a vector
# with the results of the simulation, and [[2]] a dataframe with the agentset
# as it was when the simulation terminated.
# Here we create the data structure to store the [[1]] and [[2]] of every
# simulation run.
simResults <- data.frame(matrix(NA,ncol = length(outcomeVariables)))
names(simResults) <- outcomeVariables
simW <- list()


# And this is the function that starts the simulation battery and stores its
# outcomes into the results dataframe and end-of-simulation agentset list.
# Needless to say, this can take a very long time to run.
countSimRun <- 1
for(i in 1:nrow(parameterSpace)){
  for(r in 1:nRunsPerCondition){
    seed <- sample.int(999999999, size = 1)
    print(paste(
      "Starting simulation ",
      countSimRun,
      "of",
      nrow(parameterSpace) * nRunsPerCondition
    ))
    print(paste("Seed:", seed, "- Parameter index:", i))
    print(paste("Timestamp:", Sys.time()))
    outcome <- run(
      #calibrationMode = "none", ############################################
      seed = seed,
      indexParameters = i,
      wijk = parameterSpace$wijk[i],
      initialOpinionDistribution = parameterSpace$initialOpinionDistribution[i],
      H = parameterSpace$H[i],
      distanceDecay = parameterSpace$distanceDecay[i],
      timeMax = parameterSpace$timeMax[i],
      printOpinionHistogram = FALSE,
      exportOutput = TRUE
    )
    simResults[countSimRun,] <- outcome[[1]]
    simW[[countSimRun]] <- outcome[[2]]
    countSimRun <- countSimRun + 1
  }
}

# Saving the outcomes of the simulations.
save(
  simResults,
  simW,
  parameterSpace,
  file=paste0("./simOutput/sims_", batteryID, ".RData")
)

