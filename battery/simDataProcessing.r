# ==============================================================================
# Loading output of simulation batteries from peregrine
#
# finding the results of the simulation batteries. There is one file per
# simulation battery. Each file contains:
#   - a dataframe with the parameter configuration for that battery
#   - a dataframe with the (global) outcome variables of each run in the battery
#   - a list containing the agentset ('world') dataframes as they were at
#     the end of the simulation run.
rm (list = ls( )) 

#source("geoAbm.R")
load("./cityData/geodata_Rotterdam.RData")
source("simulation.r")
source("util.r")
files <- list.files(path = "./simOutput/peregrine/")
library(compiler)
enableJIT(1)



#load(paste0("./simOutput/peregrine/", files[1]))
samples <- list()
for (d in 1:nrow(citySummary)){
  samples[[d]] <- sample(1:citySummary$n_pop[d], 500)
}

for (chunk in 1:9){
#for (chunk in 1:9){
  
  GsimResults <- data.frame()
  GsimAgents <- data.frame()
  
  #for (b in 1:2){ ################
  for (b in ((chunk * 100) - 99):(chunk * 100)) {
    #for (b in 1:length(files)){ # for every simulation battery (peregrine file)
    print(paste("Loading file", b, "of", length(files), "- Time:", Sys.time()))
    
    load(paste0("./simOutput/peregrine/", files[b]))
    
    simResults$meanAbsOpinion <- NA
    simResults$countExtremists <- NA
    simResults$countEverExtremists <- NA
    simResults$fileName <- NA
    
    #for (r in 1:2){ ##################
    for (r in 1:nrow(simResults)){ # for every simulation run of the battery
      #print(r)
      W <- simW[[r]]
      
      # Re-calculating moranI by first reconstructing the square-unit info:
      wijk <- parameterSpace$wijk[simResults$indexParameters[r]] # wijk number
      dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[wijk])
      dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
      dat$dens <- dat$inw2014 / citySummary$n_pop[wijk]
      ops <- c()
      for (l in 1:length(dat)){ # for every cell
        cell <- subset(W, W$location == dat$OBJECTID[l])
        ops[l] <- mean(cell$opinion)
      }
      dat$opinion <- ops
      
      #opClustering1 <- moranI( # opinion clustering
      #  x = dat$opinion,
      #  proxmat = proximityList1[[wijk]],
      #  dens = dat$dens,
      #  N = citySummary$n_pop[wijk]
      #)
      #opClustering2 <- moranI(
      #  x = dat$opinion,
      #  proxmat = proximityList2[[wijk]],
      #  dens = dat$dens,
      #  N = citySummary$n_pop[wijk]
      #)
      #opClustering3 <- moranI(
      #  x = dat$opinion,
      #  proxmat = proximityList3[[wijk]],
      #  dens = dat$dens,
      #  N = citySummary$n_pop[wijk]
      #)
      opAlignment1 <- moranI( # opinion-group alignment
        x = dat$pnwal2014,
        y = dat$opinion,
        proxmat = proximityList1[[wijk]],
        dens = dat$dens,
        N = citySummary$n_pop[wijk]
      )
      opAlignment2 <- moranI( # opinion-group alignment
        x = dat$pnwal2014,
        y = dat$opinion,
        proxmat = proximityList2[[wijk]],
        dens = dat$dens,
        N = citySummary$n_pop[wijk]
      )
      opAlignment3 <- moranI( # opinion-group alignment
        x = dat$pnwal2014,
        y = dat$opinion,
        proxmat = proximityList3[[wijk]],
        dens = dat$dens,
        N = citySummary$n_pop[wijk]
      )
      as <- base::merge(
        x=W[,1:2],
        y=as.data.frame(cbind(
          dat$OBJECTID,
          #    opClustering1 = opClustering1$localI,
          #    opClustering2 = opClustering2$localI,
          #    opClustering3 = opClustering3$localI,
          opAlignment1 = opAlignment1$localI,
          opAlignment2 = opAlignment2$localI,
          opAlignment3 = opAlignment3$localI
        )),
        by.x="location",
        by.y = "V1"
      )
      #W$opClustering1 <- as$opClustering1
      #W$opClustering2 <- as$opClustering2
      #W$opClustering3 <- as$opClustering3
      W$opAlignment1 <- as$opAlignment1
      W$opAlignment2 <- as$opAlignment2
      W$opAlignment3 <- as$opAlignment3
      #simResults$opClustering1[r] <- opClustering1$globalI
      #simResults$opClustering2[r] <- opClustering2$globalI
      #simResults$opClustering3[r] <- opClustering3$globalI
      simResults$opAlignment1[r] <- opAlignment1$globalI
      simResults$opAlignment2[r] <- opAlignment2$globalI
      simResults$opAlignment3[r] <- opAlignment3$globalI
      simResults$opClusteringA1 <- simResults$opClusteringA2 <-
        simResults$opClusteringA3 <- NULL
      
      # === end re-calculation ===
      
      # We can sample N agents from each run, and store them into the dataset
      # GsimAgents. These will be used for micro-level analyses.
      W$seed <- simResults$seed[r] # this is the key to join the two datasets,
      # GsimAgents and GsimResults
      sample <- data.frame()
      #sampleI <- sample(nrow(W), 200)
      #for (i in sampleI){
      for (i in samples[[parameterSpace$wijk[simResults$indexParameters[r]]]]) {
        #for (i in 1:nrow(W)){ ############################# turns OFF sampling
        sample <- rbind(sample, W[i,])
      }
      GsimAgents <- rbind(GsimAgents, sample)
      #GsimAgents <- rbind(GsimAgents, W)
      
      #simulatedWijk <- parameterSpace$wijk[simResults$indexParameters[r]]
      
      simResults$meanAbsOpinion[r] <- mean(abs(W$opinion))
      simResults$countExtremists[r] <- nrow(subset(W, W$opinion %in% c(1,-1)))
      simResults$countEverExtremists[r] <- nrow(subset(W, !is.na(W$timeFirstExtr)))
      simResults$fileName[r] <- files[b] #+ 569] #########################################
      
      
      #fileName <- append(fileName, files[b])
      #countExtremists <- append(countExtremists, nrow(subset(W, W$opinion%in%c(1,-1))))
      #countEverExtremists <- append(countEverExtremists, nrow(subset(W, !is.na(W$timeFirstPol))))
      #meanAbsOpinion <- append(meanAbsOpinion, mean(abs(W$opinion)))
      
      runResults <- cbind(
        parameterSpace[simResults$indexParameters[r],],
        simResults[r,]
      )
      runResults$initialOpinionDistribution <- 
        as.character(runResults$initialOpinionDistribution) ###
      
      if(b %in% ((1:9*100)-99) & r==1){
        GsimResults <- runResults
        #levels(GsimResults$initialOpinionDistribution) <- ###
        #  c("beta","uniform","groupBias")
      } else {
        GsimResults <- rbind(GsimResults,runResults)
        #GsimResults[(b - 1) * nrow(simResults) + r,] <- runResults
        #GsimResults$initialOpinionDistribution[(b - 1) * nrow(simResults) + r] <-
        #  as.character(runResults$initialOpinionDistribution)
      }
    }
    #rm(simW, simResults, parameterSpace)
  }
  
  print(paste("Exporting data chunk", chunk))
  save(
    GsimResults,
    GsimAgents,
    citySummary,
    file = paste0("./simOutput/completeDataset20191120_chunk_", chunk,".RDATA")
  )
}
enableJIT(0)

#
#rm (list = ls( ))
#load("./simOutput/completeDataset20191120_chunk_1.RDATA")
#nrow(GsimResults)
#nrow(GsimAgents)
#table(unique(GsimAgents$seed) %in% unique(GsimResults$seed))
#





rm (list = ls( ))
load("./cityData/geodata_Rotterdam.RData")
#load("./simOutput/completeDataset20191031.RDATA")

chunkNames <- list.files("./simOutput/")
chunkNames <- chunkNames[c(-1, -11, -12)]# Removing unnecessary items

# Re-binding the files
#d <- NULL
r <- ri <- NULL
for (f in 1:length(chunkNames)){
  print(paste("Loading and collating:", chunkNames[f]))
  load(paste0("./simOutput/", chunkNames[f]))
  if(f==1){
    r <- GsimResults
    ri <- GsimAgents
  } else {
    r <- rbind(r, GsimResults)
    ri <- rbind(ri, GsimAgents)
  }
}
GsimResults <- r
GsimAgents <- ri
rm(r, ri)



GsimResults$iniAli3 <- GsimResults$iniAli2 <- GsimResults$iniAli1 <- NA
source("script NI&PA.R")

for (i in 1: nrow(GsimResults)){
#for (i in 1:2){
  print(paste("Calc. initial alignment in run", i, "of", nrow(GsimResults)))
  run(
    timeMax =0,
    seed = GsimResults$seed[i],
    wijk = GsimResults$wijk[i],
    distanceDecay = GsimResults$distanceDecay[i],
    initialOpinionDistribution = GsimResults$initialOpinionDistribution[i],
    H = GsimResults$H[i],
    #exportOutput = FALSE,
    exportTimeSeries = FALSE
  )
  
  wijk <- GsimResults$wijk[i]
  dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[wijk])
  dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
  dat$dens <- dat$inw2014 / citySummary$n_pop[wijk]
  ops <- c()
  for (l in 1:length(dat)){ # for every cell
    cell <- subset(agents, agents$location == dat$OBJECTID[l])
    ops[l] <- mean(cell$opinion)
  }
  dat$opinion <- ops

  GsimResults$iniAli1[i] <- moranI( # opinion-group alignment
    x = dat$pnwal2014,
    y = dat$opinion,
    proxmat = proximityList1[[wijk]],
    dens = dat$dens,
    N = citySummary$n_pop[wijk]
  )$globalI
  GsimResults$iniAli2[i] <- moranI( # opinion-group alignment
    x = dat$pnwal2014,
    y = dat$opinion,
    proxmat = proximityList2[[wijk]],
    dens = dat$dens,
    N = citySummary$n_pop[wijk]
  )$globalI
  GsimResults$iniAli3[i] <- moranI( # opinion-group alignment
    x = dat$pnwal2014,
    y = dat$opinion,
    proxmat = proximityList3[[wijk]],
    dens = dat$dens,
    N = citySummary$n_pop[wijk]
  )$globalI
}
save(
  GsimResults,
  GsimAgents,
  citySummary,
  worldList,
  file = "./simOutput/completeDataset.RDATA"
)



# Now we use these values to normalize alignment and clustering scores in the
# results dataset. The max theoretical alignment is calculated in the 
# calibration script and can be accessed from the citySummary dataframe.
# For clustering, however, we do not have theoretical maxima and minima. We 
# shall therefore find the empirically observed maxima and minima.
# The citySummary dataframe is updated with the new info on min and max
# clustering.
for (wijk in 1:nrow(citySummary)){
  citySummary$minClustering1[wijk] <- min(GsimResults$opClustering1)
  citySummary$minClustering2[wijk] <- min(GsimResults$opClustering2)
  citySummary$minClustering3[wijk] <- min(GsimResults$opClustering3)
  citySummary$maxClustering1[wijk] <- max(GsimResults$opClustering1)
  citySummary$maxClustering2[wijk] <- max(GsimResults$opClustering2)
  citySummary$maxClustering3[wijk] <- max(GsimResults$opClustering3)
}
for (i in 1:nrow(GsimResults)){
  wijk <- GsimResults$wijk[i]
  GsimResults$opNetAli1[i] <-
    abs(GsimResults$opAlignment1[i]) / citySummary$maxAlignment1[wijk]
  GsimResults$opNetAli2[i] <-
    abs(GsimResults$opAlignment2[i]) / citySummary$maxAlignment2[wijk]
  GsimResults$opNetAli3[i] <-
    abs(GsimResults$opAlignment3[i]) / citySummary$maxAlignment3[wijk]
  
  GsimResults$opNetClus1[i] <-
    (GsimResults$opClustering1[i] + citySummary$minClustering1[wijk]) /
    citySummary$maxClustering1[wijk]
  GsimResults$opNetClus2[i] <-
    (GsimResults$opClustering2[i] + citySummary$minClustering2[wijk]) /
    citySummary$maxClustering2[wijk]
  GsimResults$opNetClus3[i] <-
    (GsimResults$opClustering3[i] + citySummary$minClustering3[wijk]) /
    citySummary$maxClustering3[wijk]
}
#enableJIT(0)


if (FALSE){ #################################################
# Recalculating and integrating the independent variables into
# the results dataframe.

  GsimResults$opNetAli1 <- abs(GsimResults$opAlignment1)
  GsimResults$opNetAli2 <- abs(GsimResults$opAlignment2)
  GsimResults$opNetAli3 <- abs(GsimResults$opAlignment3)
  
for (i in 1:length(worldList)){
#for (i in 1:2){
  print(paste("Adding independent variables for district", districtsNames[i]))
  #w <- worldList[[i]]
  dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[i])
  dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
  dat$grClustering1 <- MoranI( # group clustering
    x=dat$nnwal2014,
    proxmat=proximityList1[[i]],
    type = "local"
  )
  dat$grClustering2 <- MoranI( # group clustering
    x=dat$nnwal2014,
    proxmat=proximityList2[[i]],
    type = "local"
  )
  dat$grClustering3 <- MoranI( # group clustering
    x=dat$nnwal2014,
    proxmat=proximityList3[[i]],
    type = "local"
  )
  dat$grSegregation1 <- MoranI( # group segregation
    x=dat$inw2014 - dat$nnwal2014,
    y=dat$nnwal2014,
    proxmat=proximityList1[[i]],
    type = "local"
  )
  dat$grSegregation2 <- MoranI( # group segregation
    x=dat$inw2014 - dat$nnwal2014,
    y=dat$nnwal2014,
    proxmat=proximityList2[[i]],
    type = "local"
  )
  dat$grSegregation3 <- MoranI( # group segrecation
    x=dat$inw2014 - dat$nnwal2014,
    y=dat$nnwal2014,
    proxmat=proximityList3[[i]],
    type = "local"
  )
  citySummary$grClustering1[i] <- mean(dat$grClustering1)
  citySummary$grClustering2[i] <- mean(dat$grClustering2)
  citySummary$grClustering3[i] <- mean(dat$grClustering3)
  citySummary$grSegregation1[i] <- mean(dat$grSegregation1)
  citySummary$grSegregation2[i] <- mean(dat$grSegregation2)
  citySummary$grSegregation3[i] <- mean(dat$grSegregation3)
  
  citySummary$n_squ[i] <- nrow(dat)
  citySummary$n_pop[i] <- nrow(worldList[[i]])
  citySummary$n_nwa[i] <- nrow(subset(worldList[[i]], worldList[[i]]$group==1))
  citySummary$p_nwa[i] <- citySummary$n_nwa[i] / citySummary$n_pop[i]
  
  for(agent in 1:nrow(worldList[[i]])){
    locData <- dat[dat$OBJECTID==worldList[[i]]$location[agent],]
    worldList[[i]]$grClustering1[agent] <- locData$grClustering1
    worldList[[i]]$grClustering2[agent] <- locData$grClustering2
    worldList[[i]]$grClustering3[agent] <- locData$grClustering3
    
    worldList[[i]]$grSegregation1[agent] <- locData$grSegregation1
    worldList[[i]]$grSegregation2[agent] <- locData$grSegregation2
    worldList[[i]]$grSegregation3[agent] <- locData$grSegregation3
  }
}


# Here we do the same for individual agents' results, only this time we need
# the local measures. We stored the local measures in the list worldList, where
# each item is the list of agents in a given district.
# The complication (captured in the first two rows of the loop) comes from the 
# fact that individual results (GsimResults) only contain a random sample of 200
# agents per simulation run. For each of them we thus need to trace back their
# agent ID (rownames), and the wijk in which they were created. With this info,
# we can finally find which agent in the worldList corresponds to which agent
# in GsimAgents, and we can copy the independent variables.
GsimAgents$grClustering1 <- GsimAgents$grClustering2 <-
  GsimAgents$grClustering3 <- GsimAgents$grSegregation1 <-
  GsimAgents$grSegregation2 <- GsimAgents$grSegregation3 <- NA

# attaching indep cell-level stats
enableJIT(1) ##################################################################
progressBar <- round((nrow(GsimAgents) / 100) * c(1:100))
for(i in 1:nrow(GsimAgents)){
  if (i %in% progressBar){print(paste0("Progress: ",i/nrow(GsimAgents)*100,"%"))}
  if (is.na(GsimAgents$grSegregation3[i])){
    wijk <- GsimResults[GsimResults$seed == GsimAgents$seed[i],]$wijk
    w <- worldList[[wijk]][as.integer(rownames(GsimAgents[i,])),]
    GsimAgents$grClustering1[i] <- w$grClustering1
    GsimAgents$grClustering2[i] <- w$grClustering2
    GsimAgents$grClustering3[i] <- w$grClustering3
    GsimAgents$grSegregation1[i] <- w$grSegregation1
    GsimAgents$grSegregation2[i] <- w$grSegregation2
    GsimAgents$grSegregation3[i] <- w$grSegregation3
  }
}
enableJIT(0)
save(
  GsimResults,
  GsimAgents,
  citySummary,
  worldList,
  file = "./simOutput/completeDataset.RDATA"
)

require(ggplot2)
qplot(GsimAgents$grSegregation3,
      geom="histogram",
      binwidth = 0.1,  
      xlim=c(-1,1)
)
} ###############################################






# Joining citySummary dataframe (containing the indep. variables)
# and the GsimResults dataframe (containing the dep. variables)
#safe <- GsimResults
GsimResults$grClustering1 <- GsimResults$grClustering2 <- 
  GsimResults$grClustering3 <- GsimResults$grSegregation1 <-
  GsimResults$grSegregation2 <- GsimResults$grSegregation3 <- NULL

GsimResults$order <- c(1:nrow(GsimResults))
GsimResults$printOpinionHistogram <- GsimResults$exportOutput <- NULL
citySummary$WK_NR <- c(1:nrow(citySummary))
GsimResults <- base::merge(
  x = GsimResults,
  y = citySummary,
  by.x = "wijk",
  by.y = "WK_NR"
)
citySummary$WK_NR <- NULL
GsimResults <- GsimResults[order(GsimResults$order),]
rownames(GsimResults) <- c(1:nrow(GsimResults))
GsimResults$WK_NR <- GsimResults$order <- 
  GsimResults$district <- GsimResults$WK_CODE <- NULL
#table(GsimResults$seed == safe$seed)

GsimResults$intuitiveAlignment <-
  abs(GsimResults$meanOpinionG1 - GsimResults$meanOpinionG2)
hist(GsimResults$intuitiveAlignment)

r <- GsimResults
ri <- GsimAgents

ri <- base::merge(
  x = r[,c(1:4,6,8)],
  y = GsimAgents,
  by = "seed"
)


# Dealing with NaN which are in fact 0.
# These are the NaN that MoranI() returns when presented with
# perfect consensus. The only NaN that should present are thus
# the ones in the variables on opinion clustering and alignment.
for(c in 1:ncol(r)){
  if(any(is.nan(r[,c]))) {
    print(paste0("NaN found in col ", c, ": r$", names(r)[c]))
    r[,c][is.nan(r[,c])] <- 0
  }
}
for(c in 1:ncol(ri)){
  if(any(is.nan(ri[,c]))) {
    print(paste0("NaN found in col ", c, ": ri$", names(ri)[c]))
    ri[,c][is.nan(ri[,c])] <- 0
  }
}


for (i in 1:nrow(r)){
  r$expOutgr1[i] <- citySummary$expOutgr1[r$wijk[i]]
  r$expOutgr2[i] <- citySummary$expOutgr2[r$wijk[i]]
  r$expOutgr3[i] <- citySummary$expOutgr3[r$wijk[i]]
}

# We take the additive inverse of the raw scores for segregation
#r$grSegregation1 <- r$grSegregation1 * -1
#r$grSegregation2 <- r$grSegregation2 * -1
#r$grSegregation3 <- r$grSegregation3 * -1
#ri$grSegregation1 <- ri$grSegregation1 * -1
#ri$grSegregation2 <- ri$grSegregation2 * -1
#ri$grSegregation3 <- ri$grSegregation3 * -1

save(
  #GsimResults,
  #GsimAgents,
  r,
  ri,
  citySummary,
  worldList,
  file = "./simOutput/completeDataset.RDATA"
)



# ==============================================================================
# Functions to plot results

# Custom function to plot heatmaps
# This function takes two variables as argument, and prints a heatmap of their
# distribution on the  plane. It's equivalent to a scatterplot, except:
#  - it allows to plot massive numbers of observations while keeping the plot
#    readable;
#  - if the argument 'type' is different from the default "count", then the
#    columns of the resulting heatplot will be normalized (their values all sum
#    to one)
heatmap <- function(
  bins = 20,
  x,
  y,
  type = "count",
  colorLow = "white",
  colorHigh = "black",
  xlab = "",
  ylab = "",
  legend = TRUE,
  legendpositionX = 480, #480
  legendpositionY = 5,
  verbose = TRUE,
  ...
)
{
  if(verbose==TRUE) print("Printing heatmap")
  ma <- matrix(0,nrow=bins,ncol=bins)
  ds<-as.data.frame(cbind(x,y))
  for(c in 1:bins){
    if(verbose==TRUE){print(paste0(round(c/bins*100),"%"))}
    wbin <- abs(max(x) - min(x)) / bins
    hbin <- abs(max(y) - min(y)) / bins
    lowerBound <- (c - 1) * wbin + min(x)
    upperBound <- c * wbin + min(x)
    if(c == bins){
      strip <- subset(ds, ds$x >= lowerBound)
    } else {
      strip <- subset(ds, ds$x >= lowerBound & ds$x < upperBound)
    }
    
    if (nrow(strip)==0 | is.null(nrow(strip))){
      ma[,c] <- 0
    } else {
      ifelse(type == "count", denom<-1, denom<-nrow(strip))
      for (r in 1:bins){
        lowerBound <- (r - 1) * hbin + min(y)
        upperBound <- r* hbin + min(y)
        if (r==bins){
          z <- sum(strip$y >= lowerBound) / denom
        } else {
          z <- sum(strip$y >= lowerBound & strip$y < upperBound) / denom
        }
        ma[r,c] <- z
      }
    }
  }
  #ma<<-ma
  if(type=="count"){
    title <- "count of observations"
    legendLabels <- c(min(ma), round((max(ma)-min(ma))/2), max(ma))
  } else {
    title <- "proportion of observations\n(every vertical bin sums to 1)"
    legendLabels <- c(0, 0.5, 1)
  }
  colPalette <- colorRampPalette(
    c(colorLow, colorHigh),
    bias = 1,
    space = c("rgb", "Lab"),
    interpolate = c("linear", "spline"),
    alpha = FALSE
  )
  image(
    t(ma),
    col=colPalette(256),
    xlab=xlab,
    ylab=ylab,
    yaxt="n",
    main=title,
    ...
  )
  axis(
    2,
    at=seq(0,1,0.25),
    labels=round(c(
      min(y),
      (max(y) - min(y)) / 4,
      (max(y) - min(y)) / 4 * 2,
      (max(y) - min(y)) / 4 * 3,
      max(y)
    ), digits = 1)
  )
  if(legend == TRUE){
    legend(
      grconvertX(legendpositionX, "device"),
      grconvertY(legendpositionY, "device"),
      as.character(legendLabels),
      fill = colPalette(3),
      xpd = NA
    )
  }
}
if(FALSE){
  xTest <- rbeta(100000,3,9)
  yTest <- rbeta(100000,9,3)
  plot(xTest,yTest)
  heatmap(
    x = xTest,#rbeta(1000,3,9),
    y = yTest,#rbeta(1000,9,3),
    type = "count", #"count",
    bins = 10,
    #colorLow = "yellow",
    #colorHigh = "purple",
    xlab = "xTest",
    ylab = "yTest"
  )
}

# Function to plot violins.
districtViolins <- function(
  data=rr, depVar, indepVar, depVar2=NULL, depVarLabel, indepVarLabel
){
  labs_prep <- citySummary[order(citySummary[,indepVar]),]$district
  labs <- c()
  for (i in 1:length(labs_prep)){
    labs[i] <- paste0(
      labs_prep[i],
      " (",
      round(
        subset(citySummary, citySummary$district==labs_prep[i])[,indepVar],
        digits = 3
      ),
      ")"
    )
  }
  p <- ggplot(data, aes(factor(data[,indepVar]), data[,depVar]))+
    ylab(depVarLabel)
  if(!is.null(depVar2)){p <- p + geom_violin(
    aes(y = data[,depVar2]),
    fill = "white", color = "#ababab",
    scale = "width",
    draw_quantiles = 0.5
  )}
  p <- p +
    geom_violin(
      fill="gray",
      scale = "width",
      draw_quantiles = 0.5#,
      #aes(fill="gray")#, color="#ff8800")#, color="green")
    ) +  #stat_summary(fun.y=mean, geom="point", shape=4, size=2) + 
    ggtitle(paste("districts ordered by\n", indepVarLabel)) +
    theme(
      plot.margin=unit(c(0,0,0,40),"pt"),
      plot.title = element_text(hjust=0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    scale_x_discrete(labels=labs)
  return(p)
}

#testing function
#load(file= "oldRunsResults.RData")
#
#print(districtViolins(
#  depVar="opNetAli2",
#  indepVar="expOutgr2",
#  depVarLabel="Average local alignment (s=100)",
#  indepVarLabel="outgroup exposure (s=100)"))
violinPlot <- function(
  data = rri,
  depVar,
  indepVar,
  panels=NA,
  depVarLabel,
  indepVarLabel,
  panelLabel=NA, bins=10,
  fill="gray"
){
  temp <- data
  #temp[,indepVar] <- 
  #  as.numeric(cut(indepVar, breaks=(c(0:bins)/bins))) ####
  ggplot(data, 
         aes(
           cut(data[,indepVar], breaks=(c(0:bins)/bins)),
           #cut_width(indepVar, width=0.1),
           data[,depVar]#, fill=as.factor(rri$group)
         ))+
    ylab(depVarLabel) +
    xlab(indepVarLabel) +
    geom_violin(
      fill=fill,
      scale = "width",
      draw_quantiles = 0.5,
      bw="bcv"
    ) +  #stat_summary(fun.y=mean, geom="point", shape=4, size=2) +
    #facet_wrap(panels, labeller = as_labeller(panelLabel)) +
    ggplot(data, aes(data[,indepVar]), data[,depVar]) +
    geom_smooth() +
    
    theme(
      plot.margin=unit(c(0,0,0,40),"pt"),
      plot.title = element_text(hjust=0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    )
}


#
violinPlotOld <- function(
  data = rri,
  depVar,
  indepVar,
  panels=NA,
  depVarLabel,
  indepVarLabel,
  panelLabel=NA,
  bins=10,
  fill="gray"
){
  #temp <- data
  #temp[,indepVar] <- 
  #  as.numeric(cut(indepVar, breaks=(c(0:bins)/bins))) ####
  p <- ggplot(data, 
         aes(
           cut(indepVar, breaks=(c(0:bins)/bins)),
           #cut_width(indepVar, width=0.1),
           depVar#, fill=as.factor(rri$group)
         ))+
    ylab(depVarLabel) +
    xlab(indepVarLabel) +
    geom_violin(
      fill=fill,
      scale = "width",
      draw_quantiles = 0.5,
      bw="bcv"
    ) + 
    theme(
      plot.margin=unit(c(0,0,0,40),"pt"),
      plot.title = element_text(hjust=0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    )
  return(p)
}


save(
  #GsimResults,
  #GsimAgents,
  r,
  ri,
  citySummary,
  worldList,
  heatmap,
  districtViolins,
  violinPlot,
  violinPlotOld,
  file = "./simOutput/completeDataset.RDATA"
)

















































# ==============================================================================
# We still have a few missing metrics to calculate.
# These concern the opinion distribution at the start of the simulation (t=0).

rm (list = ls( ))
load("./cityData/geodata_Rotterdam.RData")
load("./simOutput/completeDataset.RDATA")
files <- list.files(path = "./simOutput/peregrine/")
source("simulation.r")
library("compiler")

rownames(r) <- 1:nrow(r)

r$iniPolarizationIndex <- r$iniIntuitiveAlignment <-
  r$iniSDopinions <- r$SDopinions <- NA
ri$iniOpAlignment3 <- ri$iniOpAlignment2 <- ri$iniOpAlignment1 <- NA
ri$maxOpAlignment3 <- ri$maxOpAlignment2 <- ri$maxOpAlignment1 <- NA

#samples <- list()
#for (d in 1:nrow(citySummary)){
#  samples[[d]] <- sample(1:citySummary$n_pop[d], 500)
#}

#for (chunk in 1:9){
#  #for (chunk in 1:9){
#  
#  GsimResults <- data.frame()
#  GsimAgents <- data.frame()
#  
#  #for (b in 1:2){ ################
#  for (b in ((chunk * 100) - 99):(chunk * 100)) {
#    #for (b in 1:length(files)){ # for every simulation battery (peregrine file)
#    print(paste("Loading file", b, "of", length(files), "- Time:", Sys.time()))
#    
#    load(paste0("./simOutput/peregrine/", files[b]))

enableJIT(1)
#for (b in 3:10) {
for (b in length(files):1) {
  print(paste("Loading file", b, "of", length(files), "- Time:", Sys.time()))
  load(paste0("./simOutput/peregrine/", files[b]))
  
  # Find which of these runs we need to plot - we'll calculate the missing
  # measures only for them.
  #c(r$H==0.6 & r$initialOpinionDistribution=="groupBias" & r$distanceDecay==2)
  
  #runs <- which(
  #  c(parameterSpace$initialOpinionDistribution=="groupBias" &
  #      parameterSpace$distanceDecay==2) |
  #    c(parameterSpace$H==0.6 & parameterSpace$distanceDecay==2) |
  #    c(parameterSpace$H==0.6 &
  #        parameterSpace$initialOpinionDistribution=="groupBias"))
  runs <- 1:nrow(parameterSpace)
  
  
  #if(length(runs) > 0){ for(q in runs){
  for(q in runs){
    print(paste("Processing simulation index", q))
    W <- simW[[q]]
    qq <- which(r$seed == simResults$seed[q])
    
    set.seed(r$seed[qq])
    
    # Updating / integrating variables at the end of the simulation run
    # (from W)
    #
    #as <- W[sample,]
    #opinionDifferences <- c()
    #for (i in 1:nrow(as)){
    #  for (j in 1:nrow(as)) {
    #    if(i != j) {
    #      oD <- abs(as$opinion[i] - as$opinion[j])
    #      opinionDifferences <- append (opinionDifferences, oD)
    #    }
    #  }
    #}
    #r$polarizationIndex[q] <- var (opinionDifferences)
    r$SDopinions[qq] <- sd(W$opinion)
    
    
    
    # Then we calculate it over the opinion distribution at the end of the
    # simulation.
    #suppressWarnings(run(
    #  timeMax = 0,
    #  seed = simResults$seed[q],
    #  wijk = parameterSpace$wijk[q],
    #  distanceDecay = parameterSpace$distanceDecay[q],
    #  initialOpinionDistribution = parameterSpace$initialOpinionDistribution[q],
    #  H = parameterSpace$H[q],
    #  exportOutput = FALSE,
    #  exportTimeSeries = FALSE,
    #  printStatusMessages = FALSE
    #))
    #sample <- sample(1:nrow(agents), size = 200, replace = FALSE)
    #as <- agents[sample,]
    #opinionDifferences <- c()
    #for (i in 1:nrow(as)){
    #  for (j in 1:nrow(as)) {
    #    if(i != j) {
    #      oD <- abs(as$opinion[i] - as$opinion[j])
    #      opinionDifferences <- append (opinionDifferences, oD)
    #    }
    #  }
    #}
    #r$iniPolarizationIndex[q] <- var (opinionDifferences)
    agents <- W
    G1 <- which(agents$group == 1, arr.ind = TRUE)
    G2 <- which(agents$group != 1, arr.ind = TRUE)
    if (parameterSpace$initialOpinionDistribution[q] == "beta"){
      agents$opinion <- rbeta(nrow(W), 3, 3, ncp = 0)
      agents$opinion <- agents$opinion * 2 - 1
    } else if (parameterSpace$initialOpinionDistribution[q] == "uniform") {
      agents$opinion <- rbeta(nrow(W), 1, 1, ncp = 0) * 2 - 1
    } else if (parameterSpace$initialOpinionDistribution[q] == "groupBias") {
      agents$opinion[G1] <- rbeta(length(G1), 3, 3.5, ncp = 0)
      agents$opinion[G2] <- rbeta(length(G2), 3.5, 3, ncp = 0)
      agents$opinion <- agents$opinion * 2 - 1
    }
    
    r$iniSDopinions[qq] <- sd(agents$opinion)
    r$iniIntuitiveAlignment[qq] <-
      abs(mean(agents$opinion[G1]) - mean(agents$opinion[G2]))
    
    if(FALSE){ ################################
    # And then we calculate what we need at the agent level.
    dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[r$wijk[qq]])
    dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
    dat$dens <- dat$inw2014 / citySummary$n_pop[r$wijk[qq]]
    ops <- c()
    for (l in 1:length(dat)){ # for every cell
      cell <- subset(agents, agents$location == dat$OBJECTID[l])
      ops[l] <- mean(cell$opinion)
    }
    dat$opinion <- ops
    
    iniAli1 <- moranI( # opinion-group alignment
      x = dat$pnwal2014, y = dat$opinion, proxmat =proximityList1[[r$wijk[qq]]],
      dens = dat$dens, N = citySummary$n_pop[r$wijk[qq]]
    )$localI
    
    iniAli2 <- moranI( # opinion-group alignment
      x = dat$pnwal2014, y = dat$opinion, proxmat =proximityList2[[r$wijk[qq]]],
      dens = dat$dens, N = citySummary$n_pop[r$wijk[qq]]
    )$localI
    
    iniAli3 <- moranI( # opinion-group alignment
      x = dat$pnwal2014, y = dat$opinion, proxmat =proximityList3[[r$wijk[qq]]],
      dens = dat$dens, N = citySummary$n_pop[r$wijk[qq]]
    )$localI
    
    
    for (location in 1:citySummary$n_squ[r$wijk[qq]]) {
      a <- which(ri$seed == r$seed[qq] & ri$index == location)
      ri$iniOpAlignment1[a] <- iniAli1[location]
      ri$iniOpAlignment2[a] <- iniAli2[location]
      ri$iniOpAlignment3[a] <- iniAli3[location]
    }
    } #############################
    #print("")
  }

}
enableJIT(0)



save(
  r,
  ri,
  citySummary,
  worldList,
  heatmap,
  districtViolins,
  violinPlot,
  violinPlotOld,
  file = "./simOutput/completeDataset.RDATA"
)


# ______________________________________________________________________________
# Initial local alignment

param <- expand.grid(
  wijk = 1:12,
  distr = unique(r$initialOpinionDistribution)
)

iniAli <- list(
  param = param,
  iniAli1 = list(),
  iniAli2 = list(),
  iniAli3 = list()
)


# We'll save the value for each of 100 repetitions in three matrices (one for
# each of the three alignment estimates). In each matrix, the row corresponds
# to the parameter configuration from iniAli, and each column is a repetition.
#squ1 <- squ2 <- squ3 <- matrix(NA, ncol=100, nrow=nrow(iniAli))

for (p in 1:nrow(param)){
  print(paste0("parameter config. ", p, " of ", nrow(param),". ", Sys.time()))
  wijk <- param$wijk[p]
  distr <- as.character(param$distr[p])
  agents <- worldList[[wijk]]
  
  ma1 <- ma2 <- ma3 <-  matrix(NA, ncol = 100, nrow = citySummary$n_squ[wijk])
  
  G1 <- which(agents$group == 1, arr.ind = TRUE)
  G2 <- which(agents$group != 1, arr.ind = TRUE)
  
  for (rep in 1:100) {
    print(paste("repetition", rep))
    if (distr == "beta"){
      agents$opinion <- rbeta(nrow(agents), 3, 3, ncp = 0)
    } else if (distr == "uniform") {
      agents$opinion <- rbeta(nrow(agents), 1, 1, ncp = 0)
    } else if (distr == "groupBias") {
      agents$opinion[G1] <- rbeta(length(G1), 3, 3.5, ncp = 0)
      agents$opinion[G2] <- rbeta(length(G2), 3.5, 3, ncp = 0)
    }
    agents$opinion <- agents$opinion * 2 - 1
    
    # And then we calculate what we need at the agent level.
    dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[wijk])
    dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
    dat$dens <- dat$inw2014 / citySummary$n_pop[wijk]
    ops <- c()
    for (l in 1:length(dat)){ # for every cell
      cell <- subset(agents, agents$location == dat$OBJECTID[l])
      ops[l] <- mean(cell$opinion)
    }
    dat$opinion <- ops
    
    ma1[,rep] <- moranI( # opinion-group alignment
      x = dat$pnwal2014, y = dat$opinion, proxmat =proximityList1[[wijk]],
      dens = dat$dens, N = citySummary$n_pop[wijk]
    )$localI
    
    ma2[,rep] <- moranI( # opinion-group alignment
      x = dat$pnwal2014, y = dat$opinion, proxmat =proximityList2[[wijk]],
      dens = dat$dens, N = citySummary$n_pop[wijk]
    )$localI
    
    ma3[,rep] <- moranI( # opinion-group alignment
      x = dat$pnwal2014, y = dat$opinion, proxmat =proximityList3[[wijk]],
      dens = dat$dens, N = citySummary$n_pop[wijk]
    )$localI
    
  }
  
  iniAli$iniAli1[[p]] <- ma1
  iniAli$iniAli2[[p]] <- ma2
  iniAli$iniAli3[[p]] <- ma3
}

save(iniAli, file = "./simOutput/iniAli.RDATA")

for (p in 1:nrow(iniAli$param)){
  print(paste0("Processing ", p, " of ", nrow(iniAli$param), ". ", Sys.time()))
  for(index in 1:nrow(iniAli$iniAli1[[p]])) {
    a <- which(
      ri$wijk == iniAli$param$wijk[p] & 
        ri$initialOpinionDistribution == as.character(iniAli$param$distr[p]) &
        ri$index == index
    )
    ri$iniOpAlignment1[a] <- iniAli$iniAli1[[p]][index,]
    ri$iniOpAlignment2[a] <- iniAli$iniAli2[[p]][index,]
    ri$iniOpAlignment3[a] <- iniAli$iniAli3[[p]][index,]
  }
}

save(iniAli, ri, file = "./simOutput/iniAli.RDATA")




# Re-calculating the max alignment in all districts.
# This updates the "maxAlignment" variables of citySummary (district level) and
# adds "maxAlignment" variables to the agent level simulation results.

#source("simulation.r")
source("util.r")

ri$maxOpAlignment3 <- ri$maxOpAlignment2 <- ri$maxOpAlignment1 <- NA
maxAlignment1 <- maxAlignment2 <- maxAlignment3 <- c()

for (wijk in 1:nrow(citySummary)){
  
  # Create wijk
  print(paste("Calculating max alignment in", districtsNames[wijk]))
  #run(timeMax=0, wijk=wijk)
  agents <- worldList[[wijk]]
  
  # Impose maximum alignment (give opinion=1 to all members of group 1, and
  # give opinion=-1 to the rest)
  for(i in 1:nrow(agents)){
    ifelse(
      agents$group[i] == 1,
      agents$opinion[i] <- -1,
      agents$opinion[i] <- 1
    ) 
  }
  
  # We take the square unit map of the district, and to each cell we give the
  # average opinion of its residents
  dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[wijk])
  dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
  dat$dens <- dat$inw2014 / citySummary$n_pop[wijk]
  ops <- c()
  for (l in 1:length(dat)){ # for every cell
    cell <- subset(agents, agents$location == dat$OBJECTID[l])
    ops[l] <- mean(cell$opinion)
  }
  dat$opinion <- ops
  
  
  # Now that he have the maximum theoretical degree of alignment, we can 
  # calculate and store how high global Moran's I scores can get.
  # We need to do this three times, one for each distance decay function.
  maxAlignment1 <- moranI( # opinion-group alignment
    x = dat$pnwal2014,
    y = dat$opinion,
    proxmat = proximityList1[[wijk]],
    dens = dat$dens,
    N = citySummary$n_pop[wijk]
  )
  maxAlignment2 <- moranI( # opinion-group alignment
    x = dat$pnwal2014,
    y = dat$opinion,
    proxmat = proximityList2[[wijk]],
    dens = dat$dens,
    N = citySummary$n_pop[wijk]
  )
  maxAlignment3 <- moranI( # opinion-group alignment
    x = dat$pnwal2014,
    y = dat$opinion,
    proxmat = proximityList3[[wijk]],
    dens = dat$dens,
    N = citySummary$n_pop[wijk]
  )
  
  # Adding max alignment to the global-level city summary
  citySummary$maxAlignment1[wijk] <- abs(maxAlignment1$globalI)
  citySummary$maxAlignment2[wijk] <- abs(maxAlignment2$globalI)
  citySummary$maxAlignment3[wijk] <- abs(maxAlignment3$globalI)
  
  # Adding max alignment to the agent-level simulation results
  for (v in 1:nrow(dat)) {
    ri$maxOpAlignment1[ri$wijk==wijk & ri$index==v] <- maxAlignment1$localI[v]
    ri$maxOpAlignment2[ri$wijk==wijk & ri$index==v] <- maxAlignment2$localI[v]
    ri$maxOpAlignment3[ri$wijk==wijk & ri$index==v] <- maxAlignment3$localI[v]
  }
}




# To citySummary we also add a column for the max SD of opinions: i.e. the
# SD measured when all western residents are given opinion 1, and non-western
# -1.
citySummary$maxSDopinions <- NA
for (d in 1:nrow(citySummary)){
  o <- rep(1, times = citySummary$n_pop[d])
  o[1:citySummary$n_nwa[d]] <- -1
  citySummary$maxSDopinions[d] <- sd(o)
}


save(
  r,
  ri,
  citySummary,
  worldList,
  heatmap,
  districtViolins,
  violinPlot,
  violinPlotOld,
  file = "./simOutput/completeDataset.RDATA"
)