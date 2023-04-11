#


#rm (list = ls())
library(ggplot2)
#source("util.r")


computeWeight <- function(ego, alter){
  return(
    1 - (abs(agents$opinion[ego] - agents$opinion[alter]) * H +
           abs(agents$group[ego] - agents$group[alter]) * (1 - H)) / 1
  )
}

# Function that returns the new opinion of the interacting agent(s), resulting
# from the interaction between ego and alter.
#
NIcomputeOpinion <- function(ego, alter) {
  
  # We store the interaction weight and the opinion difference
  w <- computeWeight(ego, alter)
  opinionDiff <- abs (agents$opinion[alter] - agents$opinion[ego])
  
  # We run a naive convergence test
  hasConverged <<- TRUE
  if (w != 0 & opinionDiff != 0 & opinionDiff != 2 ) {hasConverged <<- FALSE}
  
  # We update the opinion of ego
  oEgo <- agents$opinion[ego] + rateOpinionChange * (agents$opinion[alter] - agents$opinion[ego]) * w / 2
  if (oEgo < -1) {oEgo <- -1}
  else if (oEgo > 1) {oEgo <- 1}
  
  # If interactions are two-way (i.e. if alter influences ego while at the same
  # time ego influences alter), then we also determine the new opinion of alter.
  if (typeInteraction == "two-way"){
    oAlter <- agents$opinion[alter] + rateOpinionChange * (agents$opinion[ego] - agents$opinion[alter]) * w / 2
    if (oAlter < -1) {oAlter <- -1}
    else if (oAlter > 1) {oAlter <- 1}
    return(c(oEgo, oAlter))
  } else {
    return(oEgo)
  }
}


# Calculation of outcome measures__________________________________________
# This function is run at the end of every simulation step, and produces
# basic statistics on the system.

computeSimpleMeasures <- function(w){
  
  # We first calculate the polarization index.
  # If the population of agents is too big, we select a sample of agents
  # on which to compute the polarization index.
  # The list of sampled agents is defined only once at time = 0, while the
  # information on the sample is updated everytime we compute the measures.
  if (steps == 0) {
    if (populationSize > polSampleSize){
      polSampleList <<- sample(c(1:populationSize), polSampleSize)
    } else {
      polSampleList <<- c(1:populationSize)
    }
  }
  polSample <<- subset (w, as.integer(rownames(w)) %in% polSampleList)
  
  # Here we compute the polarization index with a nested loop to measure
  # opininion differences in the population (or in a sample of it).
  opinionDifferences <- c()
  for (i in 1:nrow(polSample)){
    for (j in 1:nrow(polSample)) {
      if(i != j) {
        oD <- abs(polSample$opinion[i] - polSample$opinion[j])
        opinionDifferences <- append (opinionDifferences, oD)
      }
    }
  }
  
  # The polarization index is the variance in the opinion differences:
  polarizationIndex <<- var (opinionDifferences)
  
  # And we measure some other outcome variables: some on the whole
  # population, and some on the two groups separately (G1 and G2).
  meanOpinionGlobal <<- mean (w$opinion)
  varOpinionGlobal <<- var (w$opinion)
  absOpGlobal <<- mean(abs(w$opinion))
  opinionG1 <- opinionG2 <- c()
  for (i in G1){opinionG1 <- append (opinionG1, w[i,"opinion"])}
  for (i in G2){opinionG2 <- append (opinionG2, w[i,"opinion"])}
  meanOpinionG1 <<- mean (opinionG1)
  meanOpinionG2 <<- mean (opinionG2)
  absOpG1 <<- mean(abs(opinionG1))
  absOpG2 <<- mean(abs(opinionG2))
  ifelse(
    length(opinionG1)<=1,
    varOpinionG1 <<-0,
    varOpinionG1 <<- var (opinionG1)
  )
  ifelse(
    length(opinionG2)<=1,
    varOpinionG2 <<-0,
    varOpinionG2 <<- var (opinionG2)
  )
}

moranI <- function(x, y = NULL, proxmat, dens = NULL, N = length(x)) {
  # Adapted from Anselin (1995, eq. 7, 10, 11)
  # https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1538-4632.1995.tb00338.x
  
  #dens: the proportion of individuals in each cell over the district population
  #if individual level data dens is.null and N is simply length of input
  #if we have aggregate data then N should be total population size (or actually
  #just a large number)
  if(is.null(y)){y <- x}
  if(is.null(dens)){dens <- rep(1/N, times = N)}
  
  #correct scaling of opinions for densities
  v1dens_ind <- rep(x, times = (dens * N))
  v1dens <- (x - mean(v1dens_ind))/sd(v1dens_ind)
  v2dens_ind <- rep(y, times = (dens * N))
  v2dens <- (y - mean(v2dens_ind)) / sd(v2dens_ind)
  
  #(density) weighted proximity matrix
  w <- proxmat
  wdens <- t(dens * t(w))
  wdens <- wdens / rowSums(wdens)
  
  #density and proximity weighted locals
  localI <- (v1dens * wdens %*% v2dens) #formula 7
  
  #correct the normalization constants
  m2 <- sum(v1dens^2 * dens)
  S0 <- N #we know the weight matrix for the individual level should add up to N
  ydens <- S0 * m2
  globalI <- sum(localI * dens * N) / ydens # formula 10/11
  
  return(list(
    globalI = globalI,
    localI = as.numeric(localI)
  ))
}



computeExtraMeasures <- function(...){
  dat <- subset(cbs100_rot, cbs100_rot$WK_CODE == districtsList[wijk])
  dat <- dat[dat$nauto2014 != 0 & dat$nnwal2014!=0,]
  dat$dens <- dat$inw2014 / citySummary$n_pop[wijk]
  ops <- c(NA, length = length(dat))
  for (l in 1:length(dat)){ # for every cell
    cell <- subset(agents, agents$location == dat$OBJECTID[l])
    ops[l] <- mean(cell$opinion)
  }
  dat$opinion <- ops
  
  
  #opClustering1 <<- moranI( # opinion clustering
  #  x = dat$opinion,
  #  proxmat = proximityList1[[wijk]],
  #  dens = dat$dens,
  #  N = citySummary$n_pop[wijk]
  #)
  #opClustering2 <<- moranI(
  #  x = dat$opinion,
  #  proxmat = proximityList2[[wijk]],
  #  dens = dat$dens,
  #  N = citySummary$n_pop[wijk]
  #)
  #opClustering3 <<- moranI(
  #  x = dat$opinion,
  #  proxmat = proximityList3[[wijk]],
  #  dens = dat$dens,
  #  N = citySummary$n_pop[wijk]
  #)
  #agents <<- base::merge(
  #  x=agents,
  #  y=as.data.frame(cbind(
  #    dat$OBJECTID,
  #    opClustering1$localI,
  #    opClustering2$localI,
  #    opClustering3$localI
  #  )),
  #  by.x="location",
  #  by.y = "V1"
  #)
  #opClusteringA1 <<- mean(abs(opClustering1))
  #opClusteringA2 <<- mean(abs(opClustering2))
  #opClusteringA3 <<- mean(abs(opClustering3))
  #opClustering1 <<- opClustering1$globalI
  #opClustering2 <<- opClustering2$globalI
  #opClustering3 <<- opClustering3$globalI
  
  opAlignment1 <<- moranI( # opinion-group alignment
    x = dat$pnwal2014,
    y = dat$opinion,
    proxmat = proximityList1[[wijk]],
    dens = dat$dens,
    N = citySummary$n_pop[wijk]
  )
  opAlignment2 <<- moranI( # opinion-group alignment
    x = dat$pnwal2014,
    y = dat$opinion,
    proxmat = proximityList2[[wijk]],
    dens = dat$dens,
    N = citySummary$n_pop[wijk]
  )
  opAlignment3 <<- moranI( # opinion-group alignment
    x = dat$pnwal2014,
    y = dat$opinion,
    proxmat = proximityList3[[wijk]],
    dens = dat$dens,
    N = citySummary$n_pop[wijk]
  )
  agents <<- base::merge(
    x=agents,
    y=as.data.frame(cbind(
      dat$OBJECTID,
      opAlignment1$localI,
      opAlignment2$localI,
      opAlignment3$localI
    )),
    by.x="location",
    by.y = "V1"
  )
  opAlignment1 <<- opAlignment1$globalI
  opAlignment2 <<- opAlignment2$globalI
  opAlignment3 <<- opAlignment3$globalI
}



# 1) World import ______________________________________________________________
worldImport <- function(){
  if (printStatusMessages) {print("Generating world dataframe.")}
  
  # We define the main features of agents: their position on the map
  # and group identity If we are running an 'uncalibrated' model, these
  # features are assigned randomly. Conversely, with the 'calibrated' version
  # we import the precalculated data about any of the city districts.
  agents <<- data.frame()
  if (calibrationMode == "none"){
    x_coor <- rnorm(populationSize, 0, 5000)
    y_coor <- rnorm(populationSize, 0, 5000)
    location <- rep(1, populationSize)
    n <- round(populationSize * groupRatio)
    group <- rep(c(1, -1), times=c(n, populationSize - n))
    #group <- as.numeric(rbinom(populationSize, 1, groupRatio))
    #group <- replace(group, group == 0, -1)
    agents <<- data.frame(cbind(
      location,
      x_coor,
      y_coor,
      group
    ))
    
    # And we define an interaction network, by default a complete network.
    # Using the argument 'inputProbmat' from the run() function, it's possible
    # to specify an arbitrary interaction network instead.
    if(is.null(inputProbmat)){
      probmat <<- matrix(1, nrow=populationSize, ncol=populationSize)
    diag(probmat) <<- 0
    } else {
      probmat <<- inputProbmat
    }
    
    
  } else {
    
    # If we are running a calibrated model, the initial world features are
    # imported from the pre-calculated calibrated world.
    load("./cityData/geodata_Rotterdam.RData", envir = globalenv())
    #cbs100_rot <<- cbs100_rot
    agents <<- as.data.frame(worldList[[wijk]])
    
    # We also import the interaction network (aka the proximity matrix between
    # square areas in the district):
    if(is.null(inputProbmat)){
      if (distanceDecay==1){
        proxmat <<- as.matrix(proximityList1[[wijk]])
      }
      if (distanceDecay==2){
        proxmat <<- as.matrix(proximityList2[[wijk]])
      }
      if (distanceDecay==3){
        proxmat <<- as.matrix(proximityList3[[wijk]])
      }
    } else {
      proxmat <<- inputProbmat
    }
    
    # We tally the agentset
    populationSize <<- nrow(agents)
    
    # If necessary, here we adjust agents' coordinates to be randomly spread in
    # their cell. This is done for aesthetic purposes (e.g. for plotting the
    # agents on a district map) - it does not affect the ABM functioning:
    # the ABM always assumes that agents are located on the centroid of their 
    # square area.
    if (spreadAgentsinCell == TRUE){
      for (i in 1:populationSize){
        agents$east[i] <<- agents$east[i] + runif(n = 1, min = -50, max = 50)
        agents$north[i] <<- agents$north[i] + runif(n = 1, min = -50, max = 50)
      }
      agents <<- cbind(agents, rd2wgs84(agents$east, agents$north))
    }
    agents <<- subset(agents, select = -c(east, north))


    # Here we calculate the probability of interaction of agents in a given cell
    # with agents from other cells. This proability depends on the proximity
    # between the two cells (found in proxmat), times the number of residents in
    # the target cell.
    nCells <<- length(table(agents$index))
    popDensity <<- c()
    for (i in 1:nCells) {
      popDensity[i] <<- nrow(subset(agents, agents$index == i))
    }
    probmat <<- matrix(NA, nrow=nCells, ncol=nCells)
    for (i in 1:nCells){
      probmat[,i] <<- proxmat[,i] * popDensity[i]
    }
  }
  
  ############################################################
  #agents$group <<- 1
  #G2 <<- sample(1:nrow(agents), size=citySummary$n_nwa[wijk])
  #for (i in G2){agents$group[i] <<- -1}
  #for(i in as.numeric(row.names(agents[order(worldList[[3]]$north),]))[1:citySummary$n_nwa[wijk]]){agents$group[i] <<- -1} #
  ############################################################


  # We create lists of agents who belong to the two groups
  G1 <<- which(agents$group == 1, arr.ind = TRUE)
  G2 <<- which(agents$group != 1, arr.ind = TRUE)


  # And we initialize a time variable.
  steps <<- 0


  # Next, we give agents an opinion.
  # For this, we have three presents: 
  #  -uniform distibution (beta distribution with alpha=beta=1)
  #  -bell-shaped distribution (beta distribution with alpha=beta=3)
  #  -bell-shaped distribution with different means for the two groups
  #   (beta distribution: alpha=3 and beta=3.5 for one group, 3.5 and 3 for the
  #    other group)
  if (printStatusMessages) {print ("Generating initial opinions")}
  if (initialOpinionDistribution == "beta"){
    agents$opinion <<- rbeta(populationSize, 3, 3, ncp = 0)
    agents$opinion <<- agents$opinion * 2 - 1
  } else if (initialOpinionDistribution == "uniform") {
    #agents$opinion <<- runif(populationSize, -1, 1)
    agents$opinion <<- rbeta(populationSize, 1, 1, ncp = 0) * 2 - 1
  } else if (initialOpinionDistribution == "groupBias") {
    o1 <- rbeta(length(G1), 3, 3.5, ncp = 0)
    o2 <- rbeta(length(G2), 3.5, 3, ncp = 0)
    for (i in 1:populationSize){
      if (agents$group[i] == 1) {
        agents$opinion[i] <<- o1[1]
        o1 <- o1[-1]
      } else {
        agents$opinion[i] <<- o2[1]
        o2 <- o2[-1]
      }
    }
    
    # We scale the opinion to the correct range
    agents$opinion <<- agents$opinion * 2 - 1
  }
  
  opDiffBetwGroupst0 <<- 
    abs(mean(agents$opinion[G1]) - mean(agents$opinion[G2]))
  opSDt0 <<- sd(agents$opinion)
  
  if (printOpinionHistogram == TRUE) {histOpinion()}
  agents$timeFirstExtr <<- agents$nIntFirstExtr <<- NA
  agents$durationPol <<- 0
  
  # The setup is complete: the last thing we need to id is to initialize the
  # time counters and to calculate the outcome measures for time = 0.
  agents$nInteractions <<- 0
  computeSimpleMeasures(agents)
  if (printStatusMessages) {print(paste("Setup completed on", Sys.time()))}
}




# 2) Simulation procedure ______________________________________________________

# This is the core of the ABM: it defines what happens in one iteration
# of the model. This function is called by the last command of the script
# (see the very last line).
#
run <- function (
  timeMax = 0,
  resetWorld = TRUE,
  seed  =  sample.int(999999999, size = 1),
  
  initialOpinionDistribution = "groupBias", # "uniform", "beta" or "groupBias"
  calibrationMode = "Rotterdam",# "Rotterdam" or "none".
  #
  # If calibrationMode = "Rotterdam", then the district needs to be specified:
  wijk = 9,                     # This is the district index in districtList
  #
  # If calibrationMode = "None", then we need to specify the number of agents we
  # want (populationSize), and the relative share of group +1 in the population 
  # (groupRatio).
  populationSize = 10,
  groupRatio = 0.5,             # Only matters if calibrationMode = "none"
  
  inputProbmat = NULL,
  H = 0.6,
  #Hc = 3,
  #S = 4,
  Mechanism = "NI", #"PA" or "NI"
  rateOpinionChange = 1,
  typeInteraction = "two-way",
  distanceDecay = 2,
  
  # Secondary parameters:
  #
  # This is the sample size for the calculation of the polarization index:
  polSampleSize = 50,
  frequencyExactConvergenceTest = 100, #time steps
  indexParameters = 0,
  exportOutput = FALSE,
  exportTimeSeries = FALSE,
  spreadAgentsinCell = FALSE,
  printStatusMessages = TRUE,
  printOpinionHistogram = TRUE,
  printVideoFrames = FALSE,
  wb_palette = colorRampPalette(c("black", "orange"))
){
  
  # The run starts by setting the model parameters. These are passed on by the
  # arguments of the run() function, or taken from the global environment.
  timeMax <<- timeMax
  seed <<- seed
  initialOpinionDistribution <<- initialOpinionDistribution
  calibrationMode <<- calibrationMode
  wijk <<- wijk
  if(calibrationMode == "none"){populationSize <<- populationSize}
  groupRatio <<- groupRatio
  inputProbmat <<- inputProbmat
  H <<- H
  Mechanism <<- Mechanism
  distanceDecay <<- distanceDecay
  rateOpinionChange <<- rateOpinionChange
  typeInteraction <<- typeInteraction
  polSampleSize <<- polSampleSize
  frequencyExactConvergenceTest <<- frequencyExactConvergenceTest
  indexParameters <<- indexParameters
  exportOutput <<- exportOutput
  spreadAgentsinCell <<- spreadAgentsinCell
  printStatusMessages <<- printStatusMessages
  printOpinionHistogram <<- printOpinionHistogram
  printVideoFrames <<- printVideoFrames
  if (printVideoFrames == TRUE){spreadAgentsinCell <<-TRUE}
  if (calibrationMode == "none"){printVideoFrames <<- spreadAgentsinCell <<-FALSE}
  wb_palette <<- wb_palette
  
  # Here we set the random seed
  RNGversion("3.5")
  set.seed(seed)
  
  # By default, every time re run() the ABM, we initialize a new agentset. 
  # However, by setting the argument resetWorld=FALSE, we can have the model
  # continue running the previous simulation.
  if (resetWorld == TRUE){worldImport()}
  if(printVideoFrames == TRUE){downloadBaseMap(zoom=16)}
  
  if (exportTimeSeries == TRUE){
    timeS <<- list()
    timeS[[1]] <<- list()
    timeS[[2]] <<- list()
  }
  
  # This loop defines what happens in each simulation step.
  for (t in 1:timeMax){
    if (timeMax == 0){break()}
    
    # At every time point, we ask all agents, taken one at a time and in
    # random order, to do the following.
    #
    #
    # Persuasive arugment model:
    #
    if (Mechanism == "PA"){
      hasConverged <<- TRUE
      forEachShuffled(agents, function(ego){
        #print(findInteractionPartner(ego))
        
        # We construct the similarity vector of the agent who is about
        # to interact.
        # *** NOTE: the similarity vector defined here is not yet based
        # on distance - but it will as soon as I complete and wire in 
        # the findInteractionPartner function (which is now unplugged).
        sims <- sapply(1:nrow(agents), function (i){
          candidateAlter = agents[i,]
          sim <- ( 4 - abs(ego$group - candidateAlter$group) -
                     abs(ego$opinion - candidateAlter$opinion) * H 
          ) / (2 + 2 * H)
          if(sim < 0){
            #print("[ERROR]") # *To be fixed*
            sim <- 0
          } else {
            # Convergence test
            if (sim > 0) {hasConverged <<- FALSE}
          }
          return(sim ^ Hc)
        })
        
        # Lastly, we pick an interaction partner "alter" based on the
        # similarity vector, and update the opinion of ego accordingly.
        alter <- agents[sample(agents$ID, 1, prob=sims), ] # placeholder   
        agents[ego$ID, "opinion"] <<- PAcomputeOpinion(ego, alter)
      })
    } else {
      
      
      # Negative influence model
      #
      # We select agents in radom order:
      shuffledAgents <- sample(1:nrow(agents))
      for (ego in shuffledAgents) {
        
        # For every agent, we first select an interaction partner, with a
        # relative probability given by the probability matrix "probmat".
        if (calibrationMode == "none"){
          alter <- sample(1:populationSize, 1, prob = probmat[ego,])
          #print(alter)
          
        } else { 
          
          # In case we are running a calibrated model, an interaction partner is
          # selected with a probability based on the proximity between the
          # square area of agents ego and the square area of the possible
          # interaction partners.
          egoCell <- agents$index[ego]
          targetCell <- sample(c(1:nCells), 1, prob=probmat[egoCell,])
          #print(distmat[egoCell, targetCell])
          #distances1 <<- append(distances1, distmat[egoCell, targetCell])############################
          repeat{
            alter <- sample(c(1:popDensity[targetCell]),1)
            alter <- as.numeric(rownames(agents[which(agents$index==targetCell),][alter,]))
            if (ego != alter) {break}
          }
          #print(paste(egoCell,targetCell, alter))
        }
        
        # Then, we simulate the interaction, as a result of which ego updates
        # her opinion:
        #agents$opinion[ego] <<- NIcomputeOpinion(ego, alter)
        newOp <- NIcomputeOpinion(ego, alter)
        agents$opinion[ego] <<- newOp[1]
        agents$nInteractions[ego] <<- agents$nInteractions[ego] + 1
        
        if (is.na(agents$timeFirstExtr[ego]) &
            isTRUE(all.equal(abs(newOp[1]), 1))){
          agents$timeFirstExtr[ego] <<- steps + 1
          agents$nIntFirstExtr[ego] <<- agents$nInteractions[ego]
        }
        if (typeInteraction == "two-way"){
          agents$opinion[alter] <<- newOp[2]
          agents$nInteractions[alter] <<- agents$nInteractions[alter] + 1
          if (is.na(agents$timeFirstExtr[alter]) &
              isTRUE(all.equal(abs(newOp[2]), 1))){
            agents$timeFirstExtr[alter] <<- steps + 1
            agents$nIntFirstExtr[alter] <<- agents$nInteractions[alter]
          }
        }
      }
      
      # Lastly, we update a simple, agent-level outcome variable: the amount of time
      # the agent has held an extreme opinion
      for(i in 1:nrow(agents)){
        if (abs(agents$opinion[i]) == 1) {
          agents$durationPol[i] <<- agents$durationPol[i] + 1
        }
      }
    }
    
    # Once all agents have interacted, we update the time index:
    steps <<- steps + 1
    
    # Finally, we conclude the simulation cycle by calculating the outcome
    # measuers, and by stopping the simulation run if the system has converged
    # or if we have simulated all the time steps that we intended to.
    computeSimpleMeasures(agents)
    
    if(printVideoFrames == TRUE){
      printPlot(paste0("NI_", districtsNames[wijk]), 1180, 720)
    }
    
    # Convergence test (NI model)
    if (hasConverged == TRUE & Mechanism == "NI"){
      if(
        
        # if consensus is reached
        (isTRUE(all.equal(meanOpinionGlobal, meanOpinionG1)) &
         isTRUE(all.equal(meanOpinionGlobal, meanOpinionG2)) &
         isTRUE(all.equal(varOpinionGlobal, 0))) |
        
        # or if bipolarization is reached
        (isTRUE(all.equal(varOpinionG1, 0)) &
         isTRUE(all.equal(varOpinionG2, 0)) &
         (isTRUE(all.equal(meanOpinionG1, 1)) | isTRUE(all.equal(meanOpinionG1, -1))) &
         (isTRUE(all.equal(meanOpinionG2, 1)) | isTRUE(all.equal(meanOpinionG2, -1)))
        )
      )
      
        # Then the system has converged.
        {hasConverged <<- TRUE} else {hasConverged <<- FALSE}
    }

    #print(steps)
    if (exportTimeSeries == TRUE){
      computeExtraMeasures() ##################################################
      timeS[[1]][[steps]] <<- createoutput() # "out"
      timeS[[2]][[steps]] <<- agents
      agents$opClustering1 <<- agents$opClustering2 <<- agents$opClustering3 <<-
        agents$opAlignment1 <<- agents$opAlignment2 <<- agents$opAlignment3 <<- NULL
    }
    
    if (printStatusMessages) {
      print(paste("Time:", steps))
      #print (paste("Polarization index: ", polarizationIndex))
      }
    if (printOpinionHistogram == TRUE) {
      dev.off()
      histOpinion()
    }
    if (hasConverged == TRUE | steps == timeMax) { ####################################################
      if (hasConverged){
        print (paste(
          "System converged to equilibrium. Polarization index:",
          polarizationIndex
        ))
      } else { # Thus, if steps == timeMax
        if (printStatusMessages){
          print(paste("System not converged within the simulated time.",
                      "Simulation terminated."))
        }
      }
      computeExtraMeasures() ##################################################
      break
    }
  }
  if (exportOutput == TRUE){return(list(createoutput(),agents))}
  if (exportTimeSeries == TRUE){return(timeS)}
}


createoutput <- function (...){
  out <- as.data.frame(cbind(
    seed,
    indexParameters,
    steps,
    polarizationIndex,
    meanOpinionGlobal,
    absOpGlobal,
    varOpinionGlobal,
    meanOpinionG1,
    absOpG1,
    varOpinionG1,
    meanOpinionG2,
    absOpG2,
    varOpinionG2,
    #opClustering1,
    #opClustering2,
    #opClustering3,
    #opClusteringA1,
    #opClusteringA2,
    #opClusteringA3,
    opAlignment1,
    opAlignment2,
    opAlignment3,
    opDiffBetwGroupst0,
    opSDt0
  ))
  return(out)
}





if (FALSE){
  run(
    timeMax = 2,
    resetWorld = TRUE,
    seed = 158749486,
    initialOpinionDistribution = "groupBias",
    calibrationMode = "Rotterdam",
    wijk = 9,
    H = 0.6,
    Mechanism = "NI",
    typeInteraction = "two-way",
    distanceDecay = 2,
    polSampleSize = 50,
    frequencyExactConvergenceTest = 100,
    printOpinionHistogram = FALSE,
  )
}


