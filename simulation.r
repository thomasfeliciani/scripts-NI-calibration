# This script defines the function that runs the ABM, run().
# An usage example can be found at the very bottom of the script.
#
#
# We start by loading the necessary resources, the first of which is the 
# calibration data file. The next block of code ensures that the data file 
# exists, and offers to download it if it is not found.
# 
cityDataURL <- "https://1drv.ms/u/s!AhmgAwgcjlrQhtFKT5Bv0y_f5TyYQQ?e=YMWlLn"

# If the calibration file cannot be found...
if (!file.exists("./cityData/geodata_Rotterdam.RData")) {
  
  #... ensure its directory exists ...
  if (!dir.exists("./cityData/")) dir.create("./cityData/")
  
  #... and ask if it should be downloaded.
  download_ <- askYesNo(paste(
    "Calibration data file not found. Do you wish to download it?\n",
    "By clicking 'yes', you will be directed to the download page.\n"
  ))
  
  # If the user does not wish to download the file right now, provide
  # instructions as to how to do it later; and then terminate the execution of
  # the script.
  if (!download_ | is.na(download_)){
    cat(
      "Results file is missing. You can download it from:",
      cityDataURL,
      "... and then place it in ./cityData/",
      "Or you can contact Thomas Feliciani at",
      paste0(
        c("cian","ail.","i@", "gm","com","tho","feli","mas.")
        [c(6,8,7,1,3,4,2,5)],
        collapse = ''), "",
      sep = "\n"
    )
    stop("")
  }
  
  # If the user wants to download the missing calibration file, then direct them
  # to the download page and tell them where the file should be placed.
  if (download_) {
    cat(
      "Please download the calibration data file 'geodata_Rotterdam.RData'",
      "and place it into the folder './cityData/'.",
      "Then try to run this script again.",
      "", sep = "\n"
    )
    browseURL(resultsURL)
    stop("The calibration data file is missing.")
  }
  rm(download_)
} 
rm(cityDataURL)

# Now we should have all we need. We can now load the calibration file,
# the library ggplot and and a utility script. Note that the order in which
# these resources are loaded matters.
load("./cityData/geodata_Rotterdam.RData", envir = globalenv())
library("ggplot2")
source("util.r")




##########################################################
################## Simulation procedure ##################
##########################################################
#
#
# This is the ABM: executing the function run() produces a simulation run,
# and returns a list containing the resulting agentset and global variables.
# The parameterization of "run" matters a great deal for the computing time:
# simulations run faster with fewer simulated time steps (i.e. lower timeMax);
# smaller districts (e.g. wijk=9); and lower size for the sample of angents on
# which we calculate the polarization index (i.e. lower polSampleSize).
#
run <- function (
  timeMax = 0,
  seed  =  NULL,
  wijk = 9, # District index from: citySummary$district
  initialOpinionDistribution = "groupBias", # "uniform", "beta" or "groupBias"
  H = 0.6,
  rateOpinionChange = 1,
  typeInteraction = "two-way",
  distanceDecay = 2, # 1 means "s=10"; 2 "s=100"; 3 "s=1000"
  polSampleSize = 50, # How many agents is the polarization index calculated on.
  exportTimeSeries = FALSE,
  spreadAgentsinCell = TRUE,
  printStatusMessages = TRUE, 
  printOpinionHistogram = TRUE
) {
  
  
  # Initialization______________________________________________________________
  
  if (printStatusMessages) print("Generating world dataframe.")
  #RNGversion("3.5") # For backward compatib. this can set an older RNG version
  if (is.null(seed)) seed <- sample.int(999999999, size = 1)
  set.seed(seed)
  
  # We initialize a time variable.
  steps <- 0
  
  # And we import from the calibration data file the appropriate district and
  # agent information.
  agents <- as.data.frame(worldList[[wijk]])
  populationSize <- nrow(agents)
  
  # We also import the interaction network (aka the proximity matrix between
  # square areas in the district). Because we have three distance decay
  # functions (1 is "s=10"; 2 is "s=100"; 3 is "s=1000"), we have three matrices
  # to choose from.
  if (distanceDecay == 1) proxmat <- proximityList1[[wijk]]
  if (distanceDecay == 2) proxmat <- proximityList2[[wijk]]
  if (distanceDecay == 3) proxmat <- proximityList3[[wijk]]
  
  
  # If necessary, here we adjust agents' coordinates to be randomly spread in
  # their cell. This is done for aesthetic purposes (e.g. for plotting the
  # agents on a district map) - it does not affect the ABM functioning:
  # the ABM always assumes that agents are located at the centroid of their 
  # square area.
  if (spreadAgentsinCell == TRUE) {
    for (i in 1:populationSize){
      agents$east[i] <- agents$east[i] + runif(n = 1, min = -50, max = 50)
      agents$north[i] <- agents$north[i] + runif(n = 1, min = -50, max = 50)
    }
    agents <- cbind(agents, rd2wgs84(agents$east, agents$north))
  }
  #agents <- subset(agents, select = -c(east, north))
  #agents <- subset(agents, select = -c(x_coor, y_coor))
  
  
  # Here we calculate the probability of interaction of agents in a given cell
  # with agents from other cells. This proability depends on the proximity
  # between the two cells (given by proxmat), times the number of residents in
  # the target cell.
  nCells <- length(table(agents$index))
  popDensity <- c()
  probmat <- matrix(NA, nrow = nCells, ncol = nCells)
  
  for (i in 1:nCells) {
    popDensity[i] <- nrow(subset(agents, agents$index == i))
    probmat[,i] <- proxmat[,i] * popDensity[i]
  }

  
  # Agents are divided into two groups:
  #   - non-western (g = 1, or G1 for short);
  #   - natives and western (g = -1, or G2 for short).
  # We save the indices of which agents belong to which group. G1 and G2 will
  # be useful later, when we'll need to select all agents from one group or the
  # other.
  G1 <- which(agents$group == 1, arr.ind = TRUE) # Non-western
  G2 <- which(agents$group != 1, arr.ind = TRUE) # Natives and western
  
  
  # Next, we give agents an opinion.
  # For this, we have three presents: 
  #   - uniform distibution (implemented as a beta distr. with alpha=beta=1)
  #   - bell-shaped distribution (beta distribution with alpha=beta=3)
  #   - bell-shaped distribution with different means for the two groups
  #   (beta distribution: alpha=3 and beta=3.5 for one group, 3.5 and 3 for the
  #    other group)
  if (initialOpinionDistribution == "beta"){
    agents$opinion <- rbeta(populationSize, 3, 3) * 2 - 1
  } else if (initialOpinionDistribution == "uniform") {
    agents$opinion <- rbeta(populationSize, 1, 1) * 2 - 1 # same as "runif()"
  } else if (initialOpinionDistribution == "groupBias") {
    o1 <- rbeta(length(G1), 3, 3.5) * 2 - 1
    o2 <- rbeta(length(G2), 3.5, 3) * 2 - 1
    agents$opinion[G1] <- o1
    agents$opinion[G2] <- o2
  }
  
  
  # And we define a function to plot the distribution of opinions across the
  # two groups. By updating this plot as the simulation progresses we can get a
  # sense of what is happening.
  histOpinion <- function(){
    ggplot(data = agents, aes(x = opinion)) +
      geom_histogram(binwidth = 0.05) +
      facet_grid(
        ~factor(as.character(group),
                levels = c("-1", "1"),
                labels = c("native and western", "non-western"))
      ) +
      scale_x_continuous(limits = c(-1.05,1)) +
      ggtitle(paste0(
        "District: ", citySummary$district[wijk],
        ". Time index: ", steps, " of ", timeMax))
  }
  
  
  # And we print the plot:
  if (printOpinionHistogram == TRUE) suppressWarnings(print(histOpinion()))
  
  
  # We also define a function that creates and populates the data structure
  # that we'll use to export the initialization parameters and some global
  # simulation statistics. Parameters and outcome stats will be saved in a 
  # vector (here implemented as a one-row data.frame).
  createoutput <- function (...){
    out <- data.frame(
      seed = seed,
      wijk = wijk,
      initialOpinionDistribution = initialOpinionDistribution,
      H = H,
      rateOpinionChange = rateOpinionChange,
      typeInteraction = typeInteraction,
      distanceDecay = distanceDecay,
      timeMax = timeMax,
      polSampleSize = polSampleSize,
      steps = steps,
      polarizationIndex = polarizationIndex,
      meanOpinionGlobal = meanOpinionGlobal,
      absOpGlobal = absOpGlobal,
      varOpinionGlobal = varOpinionGlobal,
      meanOpinionG1 = meanOpinionG1,
      absOpG1 = absOpG1,
      varOpinionG1 = varOpinionG1,
      meanOpinionG2 = meanOpinionG2,
      absOpG2 = absOpG2,
      varOpinionG2 = varOpinionG2,
      opAlignment1 = opAlignment1,
      opAlignment2 = opAlignment2,
      opAlignment3 = opAlignment3
    )
    return(out)
  }
  
  
  # We also define a function to calculate the spatial measures, i.e. agent- and
  # district-level alignment.
  computeSpatialMeasures <- function() {
    
    # Taking the district composition at the cell-level (aka "square unit").
    dat <- subset(cbs100_rot, cbs100_rot$WK_CODE == districtsList[wijk])
    dat <- dat[dat$nauto2014 != 0 & dat$nnwal2014 != 0,] 
    dat$dens <- dat$inw2014 / citySummary$n_pop[wijk]
    
    # For each cell, we take the average opinion of its residents:
    ops <- rep(NA, times = nrow(dat))
    for (cell in 1:length(dat)) {
      ops[cell] <- mean(agents[agents$location == dat$OBJECTID[cell],"opinion"])
    }
    dat$opinion <- ops
    
    # No we have all we need to calculate Moran's I:
    opAlignment1 <- moranI( # opinion-group alignment
      x = dat$pnwal2014,
      y = dat$opinion,
      proxmat = proximityList1[[wijk]],
      dens = dat$dens,
      N = citySummary$n_pop[wijk]
    )
    opAlignment2 <- moranI(
      x = dat$pnwal2014,
      y = dat$opinion,
      proxmat = proximityList2[[wijk]],
      dens = dat$dens,
      N = citySummary$n_pop[wijk]
    )
    opAlignment3 <- moranI(
      x = dat$pnwal2014,
      y = dat$opinion,
      proxmat = proximityList3[[wijk]],
      dens = dat$dens,
      N = citySummary$n_pop[wijk]
    )
    
    # We add the local I estimates to each agent, based on the cell in which
    # they reside.
    agents <- base::merge(
      x = agents,
      y = as.data.frame(cbind(
        dat$OBJECTID,
        opAlignment1 = opAlignment1$localI,
        opAlignment2 = opAlignment2$localI,
        opAlignment3 = opAlignment3$localI
      )),
      by.x="location",
      by.y = "V1"
    )
    opAlignment1 <- opAlignment1$globalI
    opAlignment2 <- opAlignment2$globalI
    opAlignment3 <- opAlignment3$globalI
    
    return(list(
      agents = agents,
      opAlignment1 = opAlignment1,
      opAlignment2 = opAlignment2,
      opAlignment3 = opAlignment3
    ))
  }
  
  
  # The initialization ends by declaring some variables: these are agent-level
  # attributes that we will update during the simulation.
  agents$timeFirstExtr <- agents$nIntFirstExtr <- NA
  agents$nInteractions <- agents$durationPol <- 0
  if (printStatusMessages) print(paste("Initialization completed:", Sys.time()))
  #if (printVideoFrames == TRUE){downloadBaseMap(zoom=16)}
  
  
  # We also prepare a data structure that allows us to export the time series:
  # a snapshot of the agentset and global statistics for each step of the
  # simulation.
  if (exportTimeSeries == TRUE){
    timeS <- list()
    timeS[[1]] <- list()
    timeS[[2]] <- list()
  }
  
  
  
  # Simulation__________________________________________________________________
  #
  # This is the main loop of the simulation. Each cycle "t" is a time step.
  for (t in 1:timeMax){
    
    # Setting the timeMax to zero means that the model simply initializes the
    # district, but the actual simulation of the interactions between agents do
    # not need to start. Thus:
    if (timeMax == 0) break()
    
    # At every time point, we ask all agents, taken one at a time and in
    # random order, to do the following.
    shuffledAgents <- sample(1:nrow(agents)) # shuffling the agentset.
    for (ego in shuffledAgents) {
      
      
      # For every agent "ego", we first select an interaction partner, "alter",
      # with a relative probability given by the probability matrix "probmat".
      # The first step to finding Alter is to choose Alter's cell, here called
      # "targetCell".
      egoCell <- agents$index[ego]
      targetCell <- sample(c(1:nCells), 1, prob = probmat[egoCell,])
      
      
      # Then, we find Alter as an agent from the target cell, chosen with
      # uniform probability. We also ensure that our Ego does not interact with
      # itself.
      repeat{ 
        alter <- sample(c(1:popDensity[targetCell]), 1)
        alter <- as.numeric(
          rownames(agents[which(agents$index == targetCell),][alter,])
        )
        if (ego != alter) break 
      }
      
      
      # Then, we simulate the interaction between Ego and Alter. At the end of
      # the interaction, Ego (and, if the interaction is two-way, also Alter)
      # will have updated their opinion.
      # The function NIcomputeOpinion() tells us what their opinion will be.
      newOp <- NIcomputeOpinion(
        agents, ego, alter, H, typeInteraction, rateOpinionChange
      )
      
      # Updating Ego's opinion accordingly:
      agents$opinion[ego] <- newOp$value[1]
      
      
      # We also keep track of how many interactions Ego has had throughout the
      # simulation. If Ego's new opinion is extreme (i.e. +1 or -1), then we
      # Write down the simulation time step in which Ego has become extremist,
      # and also how many interactions Ego has had at that point.
      agents$nInteractions[ego] <- agents$nInteractions[ego] + 1
      
      if (is.na(agents$timeFirstExtr[ego]) &
          # Notice how we test for equality using all.equal(): this is meant
          # to mitigate floating-point errors:
          isTRUE(all.equal(abs(newOp$value[1]), 1))) {
        agents$timeFirstExtr[ego] <- steps + 1
        agents$nIntFirstExtr[ego] <- agents$nInteractions[ego]
      }
      
      
      # If the interaction is two-way (the default), it means that both Ego and
      # Alter update their opinion at the end of their interaction. So, we
      # update Alter in the same way we did for Ego:
      if (typeInteraction == "two-way"){
        agents$opinion[alter] <- newOp$value[2]
        agents$nInteractions[alter] <- agents$nInteractions[alter] + 1
        if (is.na(agents$timeFirstExtr[alter]) &
            isTRUE(all.equal(abs(newOp$value[2]), 1))) {
          agents$timeFirstExtr[alter] <- steps + 1
          agents$nIntFirstExtr[alter] <- agents$nInteractions[alter]
        }
      }
    }
    
    
    # Now all agents have interacted as Ego once. We can measure some agent- and
    # district-level statistics.
    # Starting with the agent-level, we update a simple outcome variable: 
    # the amount of time each agent has held an extreme opinion.
    for(i in 1:nrow(agents)){
      if (isTRUE(all.equal(abs(agents$opinion[i]), 1))) {
        agents$durationPol[i] <- agents$durationPol[i] + 1
      }
    }
    
    
    # Next we move to the district-level outcome measures.
    # We first calculate the polarization index on a sample of agents.
    # Here we create the sample for the polarization index. Note that we
    # identify the sample agents at the start of the simulation, and we'll then
    # keep track of the agents in the sample throughout the simulation.
    if (steps == 0) {
      ifelse(
        populationSize > polSampleSize,
        polSampleIndices <- sample(c(1:populationSize), size = polSampleSize),
        polSampleIndices <- c(1:populationSize)
      )
    }
    polSample <- agents[polSampleIndices,]
    
    
    # The polarization index is the variance of opinion differences among the
    # (sample) agents. To calculate it, we start by calculating the opinion
    # difference between each pair of agents (in the sample).
    opinionDifferences <- c()
    for (i in 1:nrow(polSample)){
      for (j in 1:nrow(polSample)) {
        if(i != j) {
          oD <- abs(polSample$opinion[i] - polSample$opinion[j])
          opinionDifferences <- append (opinionDifferences, oD)
        }
      }
    }
    
    # And then take the variance:
    polarizationIndex <- var(opinionDifferences)
    
    
    # Here we measure some other outcome variables: some on the whole
    # population, and some on the two groups separately (G1 and G2).
    meanOpinionGlobal <- mean (agents$opinion)
    varOpinionGlobal <- var (agents$opinion)
    absOpGlobal <- mean(abs(agents$opinion))
    opinionG1 <- agents$opinion[G1] # G1 is non-western residents, or "g = 1"
    opinionG2 <- agents$opinion[G2] # G2 is natives and western, or "g = -1"
    meanOpinionG1 <- mean(opinionG1)
    meanOpinionG2 <- mean(opinionG2)
    absOpG1 <- mean(abs(opinionG1))
    absOpG2 <- mean(abs(opinionG2))
    varOpinionG1 <- var(opinionG1)
    varOpinionG2 <- var(opinionG2)
    
    
    # Updating the simulation time index:
    steps <- steps + 1
    
    
    # We run a naïve convergence test to get an idea of whether the simulation
    # has reached (or has nearly reached) an equilibrium. In equilibrium, no 
    # interactions can lead agents to change their opinion. We search for two 
    # equilibria: perfect consensus and perfect between-group polarization
    # (aka global alignment).
    if(
      # if consensus is reached...
      (isTRUE(all.equal(meanOpinionGlobal, meanOpinionG1)) &
       isTRUE(all.equal(meanOpinionGlobal, meanOpinionG2)) &
       isTRUE(all.equal(varOpinionGlobal, 0))) |
      
      # or if bipolarization is reached...
      (isTRUE(all.equal(varOpinionG1, 0)) &
       isTRUE(all.equal(varOpinionG2, 0)) &
       (isTRUE(all.equal(abs(meanOpinionG1), 1))) &
       (isTRUE(all.equal(abs(meanOpinionG2), 1)))
      )
    )
    
    # ... then the system has converged.
    {hasConverged <- TRUE} else {hasConverged <- FALSE}
    
    
    # We are concluding a simulation time step. Depending on the circumstances,
    # e.g. if we are saving all information for every time step
    # (exportTimeSeries == TRUE), or if the system has converged and thus this
    # is the last time step, we might have to calculate some more outcome 
    # variables and return/save them.
    if (exportTimeSeries == TRUE){
      spatialMeasures <- computeSpatialMeasures() # measuring alignment
      agents <- spatialMeasures$agents
      opAlignment1 <- spatialMeasures$opAlignment1
      opAlignment2 <- spatialMeasures$opAlignment2
      opAlignment3 <- spatialMeasures$opAlignment3
      
      timeS[[1]][[steps]] <- createoutput()
      timeS[[2]][[steps]] <- agents
      agents$opAlignment1 <- agents$opAlignment2 <- agents$opAlignment3 <- NULL
    }
    
    if (printStatusMessages) print(paste("Time:", steps, "of", timeMax))
    if (printOpinionHistogram == TRUE) {
      dev.off()
      suppressWarnings(print(histOpinion()))
    }
    
    
    # This is where we export/return the simulation results at the end of the
    # simulation run.
    if (hasConverged | steps == timeMax) {
      if (hasConverged & printStatusMessages){
        print (paste(
          "System converged to equilibrium. Polarization index:",
          polarizationIndex
        ))
      } else { # Instead, if steps == timeMax,
        if (printStatusMessages) {
          print(paste("System not converged within the simulated time.",
                      "Simulation terminated."))
        }
      }
      spatialMeasures <- computeSpatialMeasures()
      agents <- spatialMeasures$agents
      opAlignment1 <- spatialMeasures$opAlignment1
      opAlignment2 <- spatialMeasures$opAlignment2
      opAlignment3 <- spatialMeasures$opAlignment3
      break
    }
  }
  
  results <- list(
    districtStats = createoutput(),
    agents = agents
  )
  if (exportTimeSeries) results$timeSeries <- timeS
  return(results)
}





##########################################################
##################### usage example ######################
##########################################################
#
#
if (FALSE) { ################
  
  # This is how we can run a simulation of Pernis (wijk = 9), five time steps
  # long (timeMax = 5), and with our baseline parameter configuration.
  # To reproduce the exact same run twice, we also need to specify a seed
  # (e.g. seed = 12345). 
  test <- run(
    timeMax = 5,
    seed = NULL, 
    wijk = 9, # District index from: citySummary$district
    initialOpinionDistribution = "groupBias", # "uniform", "beta" or "groupBias"
    H = 0.6,
    distanceDecay = 2, # 1 means "s=10"; 2 "s=100"; 3 "s=1000"
    printOpinionHistogram = TRUE,
    exportTimeSeries = FALSE
  )
  
  # Parameter configuration and outcome variables for this run are saved here:
  c(test$districtStats)
  
  # While the agentset at the end of the simulation is saved here:
  View(test$agents)
  
  # We can visualize how opinions are distributed in the district:
  plotMap(test$agents)
  
} ################

