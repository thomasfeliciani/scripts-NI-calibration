#  Utility functions

# This script contains general purpose, auxiliary functions
# for the main script, "simulation.r".
library("ggplot2")
library("ggpubr")
library("ggmap")
library("ggsn")


# Opinion formation models______________________________________________________
#
# This calculates the interaction weight (or 'similarity') between two
# interacting agents. "agents" refers to the dataframe containing the
# simulation agentset; "ego" and "alter" are two row indices that identify the
# pair of agents of which we wish to know the similarity; "H" is the model
# parameter.
computeWeight <- function(agents, ego, alter, H){
  return(
    1 - (abs(agents$opinion[ego] - agents$opinion[alter]) * H +
           abs(agents$group[ego] - agents$group[alter]) * (1 - H)) / 1
  )
}


# Function that returns the new opinion of the interacting agent(s), resulting
# from the interaction between ego and alter.
#
NIcomputeOpinion <- function(
  agents,
  ego,
  alter,
  H,
  typeInteraction = "two-way",
  rateOpinionChange = 1
) {
  
  # We store the interaction weight and the opinion difference
  w <- computeWeight(agents, ego, alter, H)
  opinionDiff <- abs (agents$opinion[alter] - agents$opinion[ego])
  
  # We run a convergence test
  hasConverged <- TRUE
  if (w != 0 & opinionDiff != 0 & opinionDiff != 2 ) hasConverged <- FALSE
  
  # We update the opinion of ego
  oEgo <- agents$opinion[ego] +
    rateOpinionChange * (agents$opinion[alter] - agents$opinion[ego]) * w / 2
  if (oEgo < -1) oEgo <- -1
  else if (oEgo > 1) oEgo <- 1
  
  # If interactions are two-way (i.e. if alter influences ego and at the same
  # time ego influences alter), then we also determine the new opinion of alter.
  if (typeInteraction == "two-way"){
    oAlter <- agents$opinion[alter] +
      rateOpinionChange * (agents$opinion[ego] - agents$opinion[alter]) * w / 2
    if (oAlter < -1) {oAlter <- -1}
    else if (oAlter > 1) {oAlter <- 1}
    return(list(value = c(oEgo, oAlter), hasConverged = hasConverged))
  } else {
    return(list(value = oEgo, hasConverged = hasConverged))
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
  if (rbinom(1, 1, (j_opinion + 1) / 2) == 1) {  
    # If j picks a pro argument...
    if (rbinom(1, 1, (opinion + 1) / 2) == 1) { 
      # ...and i drops a pro argument, then a=0 (ineffective argument exchange)
      return(0)
    } else {                                                
      # ...and i drops a con argument, then i's opinion gets a positive push
      return(2 / S)
    }
  } else {                                                     
    # If j picks a con argument
    if (rbinom(1, 1, (opinion + 1) / 2) == 1){
      # ...and i drops a pro argument, then i's opinion gets a negative push
      return(-2 / S)
    } else {                                                   
      # ...and i drops a con argument, then a=0 (ineffective argument exchange)
      return(0)
    }
  }
}


# Moran's I_____________________________________________________________________
moranI <- function(x, y = NULL, proxmat, dens = NULL, N = length(x)) {
  # Adapted from Anselin (1995, eq. 7, 10, 11)
  # https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1538-4632.1995.tb00338.x
  #
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
  
  # (density) weighted proximity matrix
  w <- proxmat
  wdens <- t(dens * t(w))
  wdens <- wdens / rowSums(wdens)
  
  # density and proximity weighted locals
  localI <- (v1dens * wdens %*% v2dens) #formula 7
  
  # correct the normalization constants
  m2 <- sum(v1dens^2 * dens)
  S0 <- N #we know the weight matrix for the individual level should add up to N
  ydens <- S0 * m2
  globalI <- sum(localI * dens * N) / ydens # formula 10/11
  
  return(list(
    globalI = globalI,
    localI = as.numeric(localI)
  ))
}



# Graphics______________________________________________________________________
# This function plots the agents on a map and shows how groups and opinions
# are distributed across space.
plotMap <- function(agents, zoom = 15) {
  myMap <- ggmap::get_stamenmap( # downloads map tiles by Stamen Design 2021
    bbox = c(
      left = min(agents$y_coor),
      bottom = min(agents$x_coor),
      right = max(agents$y_coor),
      top = max(agents$x_coor)),
    maptype = "toner-background",
    crop = FALSE, zoom = zoom
  )
  
  mapTheme <- ggplot2::theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )
  
  groupPlot <- ggmap(myMap, darken = c(0.5, "white")) + 
    geom_point( # group density (tile / square unit)
      aes(x = y_coor, y = x_coor, color = as.factor(group)),
      data = agents, size = 0.6, alpha = 0.4,
      position = position_jitter(width = 0.00065, height = 0.00045, seed = 1)
    ) +
    scale_color_manual(
      values = c("-1" = "darkorange", "1" = "darkmagenta"),
      labels = c("natives and western", "non-western")
    ) +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    ggtitle("group")
  
  opinionPlot <- ggmap(myMap, darken = c(0.5, "white")) + 
    geom_point( # group density (tile / square unit)
      aes(x = y_coor, y = x_coor, color = opinion),
      data = agents, size = 0.6, alpha = 0.3,
      position = position_jitter(width = 0.00065, height = 0.00045, seed = 1)
    ) +
    scale_colour_gradientn(colors = c("red", "gray", "blue")) +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    ggtitle("opinion")
  
  return(ggpubr::ggarrange(
    groupPlot + mapTheme, opinionPlot + mapTheme,
    ncol = 2, hjust = -1, common.legend = FALSE
  ))
}




# Miscellanea___________________________________________________________________
#
writeToCSV <- function(variablesToExport) {
  write.table(
    variablesToExport,
    file = "Raw output.csv",
    row.names = FALSE,
    na = "",
    col.names = FALSE,
    sep = ","
  )
}

