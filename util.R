#  Utility functions

# This script contains general purpose, auxiliary functions
# for the main script: "script NI&PA".



# Opinion formation models__________________________________________________

# Function to calculate the interaction weight (or 'similarity') between
# two interacting agents
computeWeight <- function(ego, alter){
  w <- 1 - (abs(agents$opinion[ego] - agents$opinion[alter]) * H +
              abs(agents$group[ego] - agents$group[alter]) * (1 - H)) / 1
  return(w)
}


# Function that returns the new opinion of the interacting agent(s), resulting
# from the interaction between ego and alter.
#
NIcomputeOpinion <- function(ego, alter){
  
  # We store the interaction weight and the opinion difference
  w <- computeWeight(ego, alter)
  hasConverged <<- TRUE
  opinionDiff <- abs (agents$opinion[alter] - agents$opinion[ego])
  
  # We run a naive convergence test
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


# Functions that run the persuasive argument (PA) model
#
PAcomputeOpinion <- function(ego, alter){
  opinion <- ego$opinion + argument(ego$opinion, alter$opinion)
  if (opinion < -1) opinion <- -1
  else if (opinion > 1) opinion <- 1
  return(opinion)
}

# For the PA model, this function returns the effect of a "pseudo-argument":
# that is, the amount of influence on the opinion of ego that an interaction
# with alter would have produced, if alter had communicated an argument to ego.
argument <-function(opinion, j_opinion){
  if (rbinom(1, 1, ( j_opinion + 1 ) / 2) == 1){  
    # If j picks a pro argument...
    if (rbinom(1, 1, ( opinion + 1 ) / 2) == 1){ 
      # ...and i drops a pro argument, then a=0 (ineffective argument exchange)
      return(0)
    }else{                                                
      # ...and i drops a con argument, then i's opinion gets a positive push
      return(2 / S)
    }
  }else{                                                     
    # If j picks a con argument
    if (rbinom(1, 1, ( opinion + 1 ) / 2) == 1){
      # ...and i drops a pro argument, then i's opinion gets a negative push
      return(-2 / S)
    }else{                                                   
      # ...and i drops a con argument, then a=0 (ineffective argument exchange)
      return(0)
    }
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

moranI <- function(x, y=NULL, proxmat, dens=NULL, N=length(x)) {
  # Adapted from Anselin (1995, eq. 7, 10, 11)
  # https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1538-4632.1995.tb00338.x
  
  #dens: the proportion of individuals in each cell over the district population
  #if individual level data dens is.null and N is simply length of input
  #if we have aggregate data then N should be total population size (or actually
  #just a large number)
  if(is.null(y)){y <- x}
  if(is.null(dens)){dens <- rep(1/N, times=N)}
  
  #correct scaling of opinions for densities
  v1dens_ind <- rep(x, times=(dens*N))
  v1dens <- (x - mean(v1dens_ind))/sd(v1dens_ind)
  v2dens_ind <- rep(y, times=(dens*N))
  v2dens <- (y - mean(v2dens_ind))/sd(v2dens_ind)
  
  #(density) weighted proximity matrix
  w <- proxmat
  wdens <- t(dens * t(w))
  wdens <- wdens/rowSums(wdens)
  
  #density and proximity weighted locals
  localI <- (v1dens * wdens %*% v2dens) #formula 7
  
  #correct the normalization constants
  m2 <- sum(v1dens^2 * dens)
  S0 <- N #we know the weight matrix for the individual level should add up to N
  ydens <- S0 * m2
  globalI <- sum(localI * dens * N)/ydens # formula 10/11
  
  return(list(
    globalI = globalI,
    localI = as.numeric(localI)
  ))
}



computeExtraMeasures <- function(...){
  dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[wijk])
  dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
  dat$dens <- dat$inw2014 / citySummary$n_pop[wijk]
  ops <- c(NA, length = length(dat))
  for (l in 1:length(dat)){ # for every cell
    cell <- subset(agents, agents$location == dat$OBJECTID[l])
    ops[l] <- mean(cell$opinion)
  }
  dat$opinion <- ops
  

  opClustering1 <<- moranI( # opinion clustering
    x = dat$opinion,
    proxmat = proximityList1[[wijk]],
    dens = dat$dens,
    N = citySummary$n_pop[wijk]
  )
  opClustering2 <<- moranI(
    x = dat$opinion,
    proxmat = proximityList2[[wijk]],
    dens = dat$dens,
    N = citySummary$n_pop[wijk]
  )
  opClustering3 <<- moranI(
    x = dat$opinion,
    proxmat = proximityList3[[wijk]],
    dens = dat$dens,
    N = citySummary$n_pop[wijk]
  )
  agents <<- base::merge(
    x=agents,
    y=as.data.frame(cbind(
      dat$OBJECTID,
      opClustering1$localI,
      opClustering2$localI,
      opClustering3$localI
    )),
    by.x="location",
    by.y = "V1"
  )
  #opClusteringA1 <<- mean(abs(opClustering1))
  #opClusteringA2 <<- mean(abs(opClustering2))
  #opClusteringA3 <<- mean(abs(opClustering3))
  opClustering1 <<- opClustering1$globalI
  opClustering2 <<- opClustering2$globalI
  opClustering3 <<- opClustering3$globalI
  
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


# Graphics____________________________________________________________________
# This is to quickly visualize how opinions are distributed across space.
updateGraphics <- function(){
  agents$op_color <- wb_palette(10)[as.numeric(cut(agents$opinion, breaks = 10))]
  plot(
    agents$x_coor,
    agents$y_coor,
    main = 'World',
    xlab = 'Longitude',
    ylab = 'Latitude',
    axes = F,
    pch = 20,
    col = agents$op_color
  )
}


histOpinion <- function(){
  #par(mar=c(5,0,0,0))
  return(hist(
    agents$opinion,
    breaks = 20,
    main = NULL,
    xlim = c(-1,1),
    xlab = "Opinion distribution",
    yaxt = "n",
    ylab = NULL
  ))
}



# ABM general functions_________________________________________________________
# These functions are used a lot in agent-based modeling and I am likely
# to use them everytime I program an ABM in R: this is why I put them
# together here.

# This runs the function f for one random agent of the set m. In Netlogo,
# this corresponds to: "ask one-of m [f]"
forOneOf <- function (m, f){
  return (f(m[sample.int(nrow(m), 1), ]))
  # .. or, expanded:
  #i <- sample.int(nrow(m), 1)
  #target <- m[i, ]
  #return (f(target))
}

# This function runs the function f on every agent (row) of the set m, taken
# in random order. In Netlogo, this would be: "ask m [f]"
forEachShuffled <- function(m, f){
  return(sapply(sample(1:nrow(m)), function (i) f(m[i,])))
}

# This function allows to run functions f that take as argument two agents:
# i, and a randomly chosen interaction partner j, where i!=j.
# It runs the function for every agent i in m.
forEachWithPartner <- function(m, f){
  return(forEachShuffled(m, function(agent){
    id <- 0
    repeat{
      id <- sample.int(nrow(m), size = 1)
      if(id != agent$ID){
        break
      }
    }
    f(agent, m[id, ])
  }))
}



# Miscellanea___________________________________________________________________

writeToCSV <- function (variablesToExport){
  write.table(
    variablesToExport,
    file = "Raw output.csv",
    row.names=FALSE,
    na="",
    col.names=FALSE,
    sep=","
  )
}

# Old implementation - not in use:
# Function for finding an interaction partner for agent ego
if (FALSE) {
  findInteractionPartner <- function(ego){
    #alter <- 0
    if (Mechanism == "PA"){ 
      candidateAlter = agents[i,]
      sim <- ( 4 - abs(ego$group - candidateAlter$group) -
                 abs(ego$opinion - candidateAlter$opinion) * H ) / (2 + 2 * H)
      if(sim < 0){
        #print("[ERROR]") # *To be fixed*
        sim <- 0
      } else {
        if (sim > 0) {hasConverged <<- FALSE}
      }
      return(sim ^ Hc)
    } else {
      
      # For the Negative Influence model:
      
      # This implements interaction noise: with probability = interactionNoise
      # ego will get a random interaction partner regardless of their distance.
      if (rbinom(1, 1, interactionNoise) == 1 || length(nw$neighborhood[ego] == 0)){
        repeat{
          alter <- sample.int(populationSize, size = 1)
          if(alter != ego$ID){
            break
          }
        }
        
        # Standard selection for the NI model:
      } else {
        alter <- sample(nw$neighbors[[ego]], 1, prob = nw$nbWeight[[ego]])
      }
    }
    return(alter)
  }
  
  sample(agents$neighbors[[1]], 1, prob = agents$nbWeight[[1]])
  sample(agents$neighbors[[1]], 1, prob = agents$nbWeight[[1]])  
  tries <- sample(agents$neighbors[[1]], 10000, replace=TRUE, prob = agents$nbWeight[[1]])
  table (tries)
  cityData$inw2014[]
}


