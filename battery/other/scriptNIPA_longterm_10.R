
# This script is divided into two parts:
# 1) World import - Sets up the model (agents generation, calibration)
# 2) Simulation procedure - the core of the ABM
#
# Parts (1) and (2) call to different functions that, for clarity, are
# not defined in this script but organized in the two auxiliary scripts,
# "util.r" and "geoAbm.r". Auxiliary functions are grouped based on 
# their purpose.


# If needed, we can set the script address as working directory in either of
# the following ways:
#
# setwd("<current script directory>")
#setwd(dirname(parent.frame(2)$ofile))
#setwd("D:/Dropbox/Progetto/w/calibration study/Scripts")

# Cleaning up environment and loading resources.
#rm (list = ls())
require(ggplot2)
source("util.R")
source("geoAbm.R")




# 1) World import ______________________________________________________________
worldImport <- function(){
  if (printStatusMessages == TRUE) {print("Generating world dataframe.")}
  
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
  if (printStatusMessages == TRUE) {print ("Generating initial opinions")}
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
  if (printOpinionHistogram == TRUE) {histOpinion()}
  agents$timeFirstPol <<- NA
  agents$durationPol <<- 0
  
  # The setup is complete: we finalize it by calculating the outcome
  # measures for time = 0.
  computeSimpleMeasures(agents)
  if (printStatusMessages==TRUE) {print(paste("Setup completed on", Sys.time()))}
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
  
  initialOpinionDistribution = "uniform", # "uniform", "beta" or "groupBias"
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
  spreadAgentsinCell = FALSE,
  printStatusMessages = TRUE,
  printOpinionHistogram = TRUE,
  printVideoFrames = FALSE,
  wb_palette = colorRampPalette(c("black", "orange"))
){
  
  # The run starts by setting the model parameters. These are passed on by the
  # arguments of the run() function, or taken from the global environment.
  timeMax <<- timeMax
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
  set.seed(seed)
  
  # By default, every time re run() the ABM, we initialize a new agentset. 
  # However, by setting the argument resetWorld=FALSE, we can have the model
  # continue running the previous simulation.
  if (resetWorld == TRUE){worldImport()}
  if(printVideoFrames == TRUE){downloadBaseMap(zoom=16)}
  
  agents$steps <<- 0
  agents$id <<- c(1:nrow(agents))
  computeExtraMeasures()
  out <- c(
    wijk,
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
    opClustering1,
    opClustering2,
    opClustering3,
    opClusteringA1,
    opClusteringA2,
    opClusteringA3,
    opAlignment1,
    opAlignment2,
    opAlignment3
  )
  longrun[1,] <<- out
  #print(names(agents))
  longrunW <<- agents
  # This loop defines what happens in each simulation step.
  for (t in 1:timeMax){
    agents$steps <<- t
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
        
        if (is.na(agents$timeFirstPol[ego]) &
            isTRUE(all.equal(abs(newOp[1]), 1))){
          agents$timeFirstPol[ego] <<- steps
        }
        if (typeInteraction == "two-way"){
          agents$opinion[alter] <<- newOp[2]
          if (is.na(agents$timeFirstPol[alter]) &
              isTRUE(all.equal(abs(newOp[2]), 1))){
            agents$timeFirstPol[alter] <<- steps
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

    
    if (printStatusMessages == TRUE) {
      print(paste("Time:", steps))
      #print (paste("Polarization index: ", polarizationIndex))
      }
    if (printOpinionHistogram == TRUE) {
      dev.off()
      histOpinion()
    }
    #if (hasConverged == TRUE) { ####################################################
    #  print (paste("System converged to equilibrium. Polarization index:", polarizationIndex))
    #  computeExtraMeasures() ##################################################
    #  break
    #}

    #updateGraphics()
    if (steps %% freqOutcomesG == 0){
      agents$opClustering1 <<- agents$opClustering2 <<- agents$opClustering3 <<- NULL
      agents$opAlignment1 <<- agents$opAlignment2 <<- agents$opAlignment3 <<- NULL
      computeExtraMeasures()
      out <- c(
        wijk,
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
        opClustering1,
        opClustering2,
        opClustering3,
        opClusteringA1,
        opClusteringA2,
        opClusteringA3,
        opAlignment1,
        opAlignment2,
        opAlignment3
      )
      #print(out)
      #return(list(out,agents))
      longrun[nrow(longrun) + 1,] <<- out
      
      if (steps %% freqOutcomesL == 0){
        #print(names(agents))
        longrunW<<-rbind(longrunW, agents)
        save(
          longrun,
          longrunW,
          file=paste0("./simOutput/longrun_w_", wijk, ".RData")
        )
      }
    }
    if (steps == timeMax) {
     # computeExtraMeasures() ##################################################
      if (printStatusMessages == TRUE){
        #print(paste("System not converged within the simulated time.",
        #      "Simulation terminated."))
      }
      break
    }
  }
}

# This function is what starts the model, producing as many iterations
# as defined by the parameter timeMax.
#run()
#run (wijk = 1, timeMax = 0)
#run(wijk=9,
#    timeMax=200,
#    distanceDecay = 1,
#    H=0.9,
#    initialOpinionDistribution = "uniform",
#    printVideoFrames=TRUE
#)

#system.time(run (wijk = 9, timeMax = 10))
#run (wijk = 9, timeMax = 3, exportOutput = TRUE)

#run(timeMax=50, populationSize = 10, calibrationMode = "none", resetWorld = TRUE)
#print(agents$timeFirstPol)


outcomeVariables <- c(
  "wijk",
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
freqOutcomesG <-10
freqOutcomesL <-100
longrun <<- data.frame(matrix(NA,ncol = length(outcomeVariables)))
names(longrun) <- outcomeVariables
run (
  timeMax = 5000,
  wijk = 10,
  initialOpinionDistribution = "beta",
  distanceDecay = 2,
  H = 0.6,
  exportOutput = TRUE
)



#save(
#  distances1, distances2, distances3,
#  file="./testProximityMatrices.RDATA"
#)
