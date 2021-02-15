
rm (list = ls( )) 

# ==============================================================================
# Looking at the effects of outgroup exposure on the polarization dynamics
# Loading resources
source("script NI&PA.R")
source("util.R")
source("geoAbm.R")




testSim <- run(
  timeMax = 10,
  seed  =  12345,
  initialOpinionDistribution = "groupBias", # "uniform", "beta" or "groupBias"
  calibrationMode = "Rotterdam",# "Rotterdam" or "none".
  wijk = 3,
  H = 0.6,
  distanceDecay = 2,
  exportOutput = TRUE,
  printOpinionHistogram = FALSE
)

a_mixed <- testSim[[2]]
a_segr <- testSim[[2]]

boxplot(log10(a_mixed$nIntFirstExtr), log10(a_segr$nIntFirstExtr))

boxplot(abs(a_mixed$opAlignment2), abs(a_segr$opAlignment2))









# ==============================================================================
# Inspecting ABM internal validity.
# We run the model in a stylized example, to see if it reproduces the known
# behavior of the NI model.

# Loading resources
source("script NI&PA.R")
source("util.R")
source("geoAbm.R")
require(network)


# Finding values of parameter H that make negative influence plausibly unlikely
testSim <- function(
  n = 1000,
  grp = c(1,-1),
  opi = 3
) {
  res <<- matrix(NA, nrow = n, ncol = 21)
  for (h in 0:20){
    H <<- 1/20 * h
    print(paste0("Testing for H=",H))
    for (i in 1:n){
      agents<<-as.data.frame(cbind(
        grp, # group
        c(rbeta(2, opi, opi, ncp = 0) *2 -1)
      ))
      names(agents) <<- c("group","opinion")
      w <- computeWeight(1,2)
      ifelse(w<0,res[i,h+1]<<-1,res[i,h+1]<<-0)
    }
  }
  r<-c()
  for (i in 1:21){r[i]<-sum(res[,i])/n}
  return(r)
}



plotTest <- function(a, b){
  plot(
    (c(1:21) -1) * 0.05,
    a,
    ylim = c(0,1),
    col = "black",
    #main=title,
    xlab="parameter H",
    ylab="probability of negative influence",
    cex.lab=2, cex.axis=2
  )
  mtext("Black = ingroup interaction", side=3, line=0.5, cex=2)
  mtext("Red = outgroup interaction", side=3, line=2, col="red",cex=2)
  lines(
    (c(1:21) -1) * 0.05,
    a
  )
  points(
    (c(1:21) -1) * 0.05,
    b,
    col = "red"
  )
  lines(
    (c(1:21) -1) * 0.05,
    b,
    col = "red",
    lty = 2
  )
}

for(i in c(1,3,5,10,20)){
  print(paste0("Printing plots for alpha=beta=",i))
  r1<-testSim(n=20000, grp=c(1,1), opi=i)
  r2<-testSim(n=20000, grp=c(1,-1), opi=i)
  png(paste0("outputGraphics/ab",i,".png"), 
      width = 1000, height = 500,
      units = "px", pointsize = 10, 
      res = NA)
  par(mfrow=c(1,2))
  hist(
    c(rbeta(20000, i, i, ncp = 0) *2 -1),
    breaks=20,
    xlim=c(-1,1),
    xlab="opinion",
    main=paste0("Beta distribution alpha=beta=", i),
    cex.lab=2, cex.axis=2, cex.main=2
  )
  plotTest(a=r1,b=r2)
  dev.off()
}

rm(r1,r2,res,agents,testSim,plotTest)



# Here we check that the ABM behaves as intended.
# We start by creating an interaction network with two components, one for 
# each demographic group.
# Because the two groups are disconnetted, we expect that it's very likely that
# the two groups develop perfect internal consensus, and that there is going to be
# a little opinion difference between the two groups.
n <- 10
testProbmat <- matrix(0, nrow=n, ncol=n)

# We create cliques between the first 5 (group 1) and last 5 agents (group -1)
for (i in 1:n){
  for (j in 1:n){
    if (i < (n / 2 + 1) & j < (n / 2 + 1)){ #
      testProbmat[i,j] <- 1
      testProbmat[j,i] <- 1
    }
    if (i > (n / 2) & j > (n / 2)){
      testProbmat[i,j] <- 1
      testProbmat[j,i] <- 1
    }
  }
}
diag(testProbmat) <- 0
set.seed(12345)
plot.network(
  as.network(testProbmat),
  vertex.col=agents$opinion,
  vertex.cex = 5#,
  #xlim=c(0,5),
  #ylim=c(0,10)
)

# Here's a quick function to run the ABM with our testProbmat.
testRun <- function(){
  simulation <-   run(
    inputProbmat = testProbmat,
    timeMax=100,
    populationSize = 10,
    H = 0.6,
    initialOpinionDistribution = "uniform",
    calibrationMode = "none",
    resetWorld = TRUE,
    printOpinionHistogram = FALSE,
    exportOutput =TRUE
  )
  return(simulation)
}

# We can see that we consistently obtain the expected results.
test <- testRun()
test[[2]]
print(paste("Polarization index =", test[[1]][4]))

# Now we do the same, but by connetting the two components, thereby allowing
# intergroup interactions.
testProbmat[5,6] <- 1
testProbmat[6,5] <- 1
plot.network(
  as.network(testProbmat),
  #vertex.col=agents$opinion,
  vertex.cex = 5
)


# This time we expect to observe between-group bi-polarization a lot more often,
# since agents 5 and 6, who come from different groups, can now interact,
# trigger mutual repulsion, and spread extreme attitudes in their groups.
#
# Again, the model behaves as expected:
test <- testRun()
test[[2]]
print(paste("Polarization index =", test[[1]][4]))


# Now we inspect the dynamic in the two connected cliques by plotting how the
# opinions in the network change over time. We save the output frames in a
# dedicated folder, ./outputGraphics/frames/
x <- plot.network(as.network(testProbmat)) # This saves the coordinates
                                           # of the nodes in the plot
x <- round(x, digits = 2)
for(t in 1:200){
  ifelse(t==1, init<-TRUE, init<-FALSE)
  run(
    seed = sample(1000,1),
    inputProbmat = testProbmat,
    timeMax=t,
    populationSize = 10,
    H = 0.6,
    initialOpinionDistribution = "uniform",
    calibrationMode = "none",
    resetWorld = init,
    printOpinionHistogram = FALSE,
    exportOutput =TRUE
  )
  colorizeOp <- colorRamp(c("red","blue"))
  cols <- colorizeOp((agents$opinion + 1)/2)
  agents$color <- rgb(cols, maxColorValue=256)
  #set.seed(1235) # Specifying a seed prevents nodes from being plotted in
                  # in different places across animation frames
  png(paste0("outputGraphics/frames/smallNW_", steps,  ".png"), 
      width = 500, height = 500,
      units = "px", pointsize = 10, 
      res = NA)
  plot.network(
    main=paste("time step:", steps),
    cex.main = 3,
    as.network(testProbmat),
    coord=x,
    vertex.col=agents$color,
    vertex.cex = 5
  )
  dev.off()
}

# The NI model can replicate known behavior, and the proximity matrix seems to
# function as intended.
rm(testRun, testProbmat, n, x)




if (FALSE){
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
n <- 2
testProbmat <- matrix(1, nrow=n, ncol=n)
diag(testProbmat) <- 0

testRun <- function(){
  simulation <-   run(
    inputProbmat = testProbmat,
    timeMax=100,
    populationSize = 2,
    H = 3,
    initialOpinionDistribution = "groupBias",
    calibrationMode = "none",
    resetWorld = TRUE,
    printOpinionHistogram = FALSE,
    exportOutput =TRUE
  )
  return(simulation)
}

# We can see that we consistently obtain the expected results.
test <- testRun()
results <- test[[1]]
names(results) <- outcomeVariables
print(results)
#print(paste("Polarization index =", test[[1]][4]))
}




# Looking for min and max values of Moran's I
rm (list = ls( )) 

# Loading resources
source("script NI&PA.R")
source("util.R")

maxAlignment1 <- maxAlignment2 <- maxAlignment3 <- c()

for (wijk in 1:12){
  
  # Create wijk
  run(timeMax=0, wijk=wijk)
  print(paste("Calculating max alignment in", districtsNames[wijk]))
  
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
  ops <- c(NA, length = length(dat))
  for (l in 1:length(dat)){ # for every cell
    cell <- subset(agents, agents$location == dat$OBJECTID[l])
    ops[l] <- mean(cell$opinion)
  }
  dat$opinion <- ops
  
  # Now that he have the maximum theoretical degree of alignment, we can 
  # calculate and store how high global Moran's I scores.
  # NOTE: we do this three times, one for each distance decay function.
  maxAlignment1[wijk] <- stats::weighted.mean(abs(moranI(
    x = dat$pnwal2014,
    y = dat$opinion,
    proxmat = proximityList1[[wijk]]
  )$localI) , w = )
  maxAlignment2[wijk] <- mean(abs(MoranI(
    x = dat$pnwal2014,
    y = dat$opinion,
    proxmat = proximityList2[[wijk]],
    type = "local"
  )))
  maxAlignment3[wijk] <- mean(abs(MoranI(
    x = dat$pnwal2014,
    y = dat$opinion,
    proxmat = proximityList3[[wijk]],
    type = "local"
  )))
}

# We add this info to the district summary stats and save.
citySummary <- as.data.frame(cbind(
  citySummary,
  maxAlignment1,
  maxAlignment2,
  maxAlignment3
))
save(
  districtsList,
  districtsNames,
  rd2wgs84,
  MoranI,
  MoranScatterplot,
  cbs100_rot,
  citySummary,
  proximityList1,
  proximityList2,
  proximityList3,
  worldList,
  file="./cityData/geodata_Rotterdam.RData"
)



# ==============================================================================
# Exploring long runs 
#
# Here we plot the change of the main outcome variable across time, for a few
# selected simulation runs that we ran for up to 5000 simulation steps.
rm (list = ls( )) 
load("./cityData/geodata_Rotterdam.RData")

# Here we select the simulation run to look at:
load("./simOutput/longrun/longrun_w_1.RData")

# This function plots a signal/time series
plt <- function(x, name){
  plot(longrun$steps, x, type="l", log="x", main=name, xlab="time", ylab="",
       ylim=c(0,1)#c(min(x),max(x))
  )
  abline(v=200, col="red")
  abline(h=x[1], lty=3)
}
#plt(longrun$polarizationIndex, "polarization index")

# This block of code produces a summary of the long term dynamics of the
# simulation run. It shows the time series for the main outcome variables.
par(mfrow=c(3,3))
plt(longrun$polarizationIndex, "polarization index")
plt(longrun$varOpinionGlobal, "opinion variance")
plt(longrun$absOpGlobal, "mean absolute opinion")
plt(longrun$opClustering1, "global op. opclustering (s=10)")
plt(longrun$opClustering2, "global op. clustering (s=100)")
plt(longrun$opClustering3, "global op. clustering (s=1000)")
plt(longrun$opAlignment1, "global op. alignment (s=10)")
plt(longrun$opAlignment2, "global op. alignment (s=100)")
plt(longrun$opAlignment3, "global op. alignment (s=1000)")
rm(plt)




# ==============================================================================
# Loading output of simulation batteries from peregrine
#
# finding the results of the simulation batteries. There should be one file per
# simulation battery. Each file contains:
#   - a dataframe with the parameter configuration for that battery
#   - a dataframe with the (global) outcome variables of each run in the battery
#   - a list containing the agentset ('world') dataframes as they were at
#     the end of the simulation run.
rm (list = ls( )) 

source("script NI&PA.R")
source("util.R")
source("geoAbm.R")
load("./cityData/geodata_Rotterdam.RData")
files <- list.files(path = "./simOutput/peregrine/")
#library(compiler)
#enableJIT(1)
GsimResults <- data.frame()
GsimAgents <- data.frame()


#load(paste0("./simOutput/peregrine/", files[1]))





#for (b in 570:length(files)){ ####################################
for (b in 1:length(files)){ # for every simulation battery (peregrine file)
  print(paste("Loading file", b, "of", length(files), "- Time:", Sys.time()))
  
  load(paste0("./simOutput/peregrine/", files[b]))
  
  simResults$meanAbsOpinion <- NA
  simResults$countExtremists <- NA
  simResults$countEverExtremists <- NA
  simResults$fileName <- NA
  #for (r in 1:2){
  for (r in 1:nrow(simResults)){ # for every simulation run of the battery
    #print(r)
    W <- simW[[r]]
    
    # We can sample 200 agents from each run, and store them into the dataset
    # GsimAgents. These will be used for micro-level analyses.
    W$seed <- simResults$seed[r] # this is the key to join the two datasets,
                                  # GsimAgents and GsimResults
    sampleI <- sample(nrow(W), 200)
    sample <- data.frame()
    for (i in sampleI){
      sample <- rbind(sample, W[i,])
    }
    GsimAgents <- rbind(GsimAgents, sample)
    #GsimAgents <- rbind(GsimAgents, W)

    #simulatedWijk <- parameterSpace$wijk[simResults$indexParameters[r]]
    
    #simResults$meanAbsOpinion[r] <- mean(abs(W$opinion))
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
    if(b==1 & r==1){
      GsimResults <- runResults
      levels(GsimResults$initialOpinionDistribution) <-
        c("beta","uniform","groupBias")
    } else {
      GsimResults[(b - 1) * nrow(simResults) + r,] <- runResults
      GsimResults$initialOpinionDistribution[(b - 1) * nrow(simResults) + r] <-
        as.character(runResults$initialOpinionDistribution)
    }
  }
  #rm(simW, simResults, parameterSpace)
}
enableJIT(0)


names(GsimAgents)
plot(GsimAgents$timeFirstExtr,GsimAgents$nIntFirstExtr)

table(GsimResults$initialOpinionDistribution)
#GsimResults_b <- GsimResults
#GsimAgents_b <- GsimAgents
#GsimResults_a <- GsimResults[GsimResults$fileName != "sims_671395284.RData",]
#GsimAgents_a <- subset(GsimAgents, GsimAgents$seed %in% unique(GsimResults_a$seed))

#GsimResults <- rbind(GsimResults_a, GsimResults_b)
#GsimAgents <- rbind(GsimAgents_a, GsimAgents_b)

#save(
#  GsimResults_a,
#  GsimAgents_a,
#  GsimResults_b,
#  GsimAgents_b,
#  citySummary,
#  file = "./simOutput/completeDataset.RDATA"
#)

save(
  GsimResults,
  GsimAgents,
  citySummary,
  file = "./simOutput/completeDataset.RDATA"
)



rm (list = ls( ))
load("./cityData/geodata_Rotterdam.RData")
load("./simOutput/completeDataset.RDATA")

# Here I tried to normalize the "global alignment", by dividing the value we
# measured at the end of each simulation run by what I thought was the theoretical
# maximum.
if (FALSE){
for (i in 1:nrow(GsimResults)){
  GsimResults$opNetAli1[i] <-
    GsimResults$opAlignment1[i] / citySummary$maxAlignment1[GsimResults$wijk[i]]
  GsimResults$opNetAli2[i] <-
    GsimResults$opAlignment2[i] / citySummary$maxAlignment2[GsimResults$wijk[i]]
  GsimResults$opNetAli3[i] <-
    GsimResults$opAlignment3[i] / citySummary$maxAlignment3[GsimResults$wijk[i]]
}

# Then I tried to see if there were alignment values higher than the theoretical
# maximum. If there were, the normalized score would be higher than 1
if (
  any(GsimResults$opNetAli1 > 1) |
  any(GsimResults$opNetAli2 > 1) |
  any(GsimResults$opNetAli3 > 1)
){print("problem")}
  
# Turns out that our theoretical maximum was flawed, as we had runs that generated
# stronger alignment than what we thought was the max.
# At this point I wanted to look into cases where the observed score was higher
# than thought possible, so that we can understand what is going on.
  test <- subset(
    GsimResults,
    GsimResults$opAlignment3 == max(GsimResults$opAlignment3, na.rm = TRUE)
  )
  load(paste0("./simOutput/peregrine/", test$fileName))
  W <- simW[[22]]
  hist(W$opinion)
  hist(W$opAlignment2)
  mean(abs(W$opAlignment2))
  citySummary$maxAlignment2[10]
  test1 <- subset(W, W$group == 1)
  test2 <- subset(W, W$group != 1)
  mean(test1$opinion)
  mean(test2$opinion)
  dat$op_color <- wb_palette(10)[as.numeric(cut(dat$opinion, breaks = 20))]
  plot(dat$l, dat$j, pch=21, bg = dat$op_color)
  dat$gr_color <- wb_palette(10)[as.numeric(cut(dat$pnwal2014, breaks = 20))]
  plot(dat$l, dat$j, pch=21, bg = dat$gr_color)
  
  dat$ali2 <- abs(MoranI(
  x = dat$pnwal2014,
  y = dat$opinion,
  proxmat = proximityList2[[wijk]],
  type = "local"
))
dat$ali_color <- wb_palette(10)[as.numeric(cut(dat$ali2, breaks = 20))]

for (i in 1:nrow(dat)){
  ifelse(
    dat$ali2[i] > citySummary$maxAlignment2[10],
    dat$ali_color[i] <- "black",
    dat$ali_color[i] <- "white"
  )
}
plot(dat$l, dat$j, pch=21, bg = dat$ali_color,
     main="Black=cells with out-of-bound alignment")

# Looking at these descriptives we learned that, when the simulation is close to
# convergence to a moderate opinion, tiny opinion difference can be distributed
# in such a way that the measure for alignment scores higher than it would score
# in case of max between-group bi-polarization.
# In essence, we were wrong about the max theoretically possible degree of
# alignment. For lack of a way to find the actual theoretical maximum, in the
# next lines of code we proceed with normalizing alignment (and clustering,
# which is based on the same index) using the *observed* maximum among
# simulations of the same district.
}
  

# Here we find the theoretical max alignment and clustering (and min clustering)
# for each district. We'll need this to normalize the alignment and clustering
# score and to be able to compare results across districts.
library(compiler)
enableJIT(1)
for (w in 1:length(districtsList)){
  districtResults <- subset(GsimResults, GsimResults$wijk == w)
  citySummary$maxAlignment1[w] <- max(districtResults$opAlignment1, na.rm=TRUE)
  citySummary$maxAlignment2[w] <- max(districtResults$opAlignment2, na.rm=TRUE)
  citySummary$maxAlignment3[w] <- max(districtResults$opAlignment3, na.rm=TRUE)
  citySummary$maxClustering1[w] <- max(districtResults$opClustering1, na.rm=TRUE)
  citySummary$maxClustering2[w] <- max(districtResults$opClustering2, na.rm=TRUE)
  citySummary$maxClustering3[w] <- max(districtResults$opClustering3, na.rm=TRUE)
  citySummary$minClustering1[w] <- min(districtResults$opClustering1, na.rm=TRUE)
  citySummary$minClustering2[w] <- min(districtResults$opClustering2, na.rm=TRUE)
  citySummary$minClustering3[w] <- min(districtResults$opClustering3, na.rm=TRUE)
}






rm (list = ls( ))
load("./cityData/geodata_Rotterdam.RData")
load("./simOutput/completeDataset.RDATA")


save(
  GsimResults,
  GsimAgents,
  citySummary,
  worldList,
  file = "./simOutput/completeDataset.RDATA"
)
save(
  districtsList,
  districtsNames,
  rd2wgs84,
  MoranI,
  MoranScatterplot,
  cbs100_rot,
  citySummary,
  proximityList1,
  proximityList2,
  proximityList3,
  worldList,
  file="./cityData/geodata_Rotterdam.RData"
)





# Now we use these values to normalize alignment scores in the results dataset:
for (i in 1:nrow(GsimResults)){
  wijk <- GsimResults$wijk[i]
  GsimResults$opNetAli1[i] <-
    GsimResults$opAlignment1[i] / citySummary$maxAlignment1[wijk]
  GsimResults$opNetAli2[i] <-
    GsimResults$opAlignment2[i] / citySummary$maxAlignment2[wijk]
  GsimResults$opNetAli3[i] <-
    GsimResults$opAlignment3[i] / citySummary$maxAlignment3[wijk]
  
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
enableJIT(0)


if (FALSE){
# Recalculating and integrating the independent variables into
# the results dataframe.
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
}

# Here we add the independent measures to the results dataset of simulation runs.
# The global measures that we need are stored in citySummary.
# This code is now removed because citySummary and GsimResults are joined later
#
#library(compiler)
#enableJIT(1)
#
#GsimResults$expIngroup <- GsimResults$expOutgroup <-
#  GsimResults$NWAexpIngroup <- GsimResults$NWAexpOutgroup <-
#  GsimResults$NATexpIngroup <- GsimResults$NATexpOutgroup <- NULL
#
#for (r in 1:nrow(GsimResults)){
#  wijk <- GsimResults$wijk[r]
#  GsimResults$expIngr1[r] <- citySummary$expIngr1[wijk]
#  GsimResults$expIngr2[r] <- citySummary$expIngr2[wijk]
#  GsimResults$expIngr3[r] <- citySummary$expIngr3[wijk]
#  
#  GsimResults$expOutgr1[r] <- citySummary$expOutgr1[wijk]
#  GsimResults$expOutgr2[r] <- citySummary$expOutgr2[wijk]
#  GsimResults$expOutgr3[r] <- citySummary$expOutgr3[wijk]
#  
#  GsimResults$NWAexpIngr1[r] <- citySummary$NWAexpIngr1[wijk]
#  GsimResults$NWAexpIngr2[r] <- citySummary$NWAexpIngr2[wijk]
#  GsimResults$NWAexpIngr3[r] <- citySummary$NWAexpIngr3[wijk]
#  
#  GsimResults$NWAexpOutgr1[r] <- citySummary$NWAexpOutgr1[wijk]
#  GsimResults$NWAexpOutgr2[r] <- citySummary$NWAexpOutgr2[wijk]
#  GsimResults$NWAexpOutgr3[r] <- citySummary$NWAexpOutgr3[wijk]
#
#  GsimResults$NATexpIngr1[r] <- citySummary$NATexpIngr1[wijk]
#  GsimResults$NATexpIngr2[r] <- citySummary$NATexpIngr2[wijk]
#  GsimResults$NATexpIngr3[r] <- citySummary$NATexpIngr3[wijk]
#  
#  GsimResults$NATexpOutgr1[r] <- citySummary$NATexpOutgr1[wijk]
#  GsimResults$NATexpOutgr2[r] <- citySummary$NATexpOutgr2[wijk]
#  GsimResults$NATexpOutgr3[r] <- citySummary$NATexpOutgr3[wijk]
#}

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




# (for backward compatibility) Dropping old or redundant variables
#names(citySummary)
#citySummary$globalI_propNWA <- citySummary$globalI_biv_prop <-
#  citySummary$globalI_countNWA <- citySummary$globalI_countNAT <-
#  citySummary$globalI_biv_count <- NULL


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

save(
  GsimResults,
  GsimAgents,
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


# ==============================================================================
# Testing agents' sampling procedure
rm (list = ls( ))
load("./cityData/geodata_Rotterdam.RData")
load("./simOutput/completeDataset.RDATA")
require(ggplot2)
r <- GsimResults
#ri <- GsimAgents

ri <- base::merge(
  x = r[,c(1:4,6)],
  y = GsimAgents,
  by = "seed"
)

rm(GsimResults, GsimAgents)

# We take the additive inverse of the raw scores for segregation
r$grSegregation1 <- r$grSegregation1 * -1
r$grSegregation2 <- r$grSegregation2 * -1
r$grSegregation3 <- r$grSegregation3 * -1
ri$grSegregation1 <- ri$grSegregation1 * -1
ri$grSegregation2 <- ri$grSegregation2 * -1
ri$grSegregation3 <- ri$grSegregation3 * -1

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

if (FALSE){
  # Making sure sampling works
  seeds <- sort(unique(r$seed))
  r$EmeanAbsOpinion <- r$EopClustering1 <- r$EopClustering2<- r$EopClustering3 <- NA
  r$EopAlignment1 <- r$EopAlignment2 <- r$EopAlignment3 <- NA
  #for (i in 1:5){
  for (i in 1:3000){
    #for (i in 1:nrow(r)){
    print(paste(i, "/", 5000))
    #s <- r$seed[i]
    s <- seeds[i]
    row <- which(r$seed == s)
    sample <- subset(ri, ri$seed == s)
    r$EmeanAbsOpinion[row] <- mean(abs(sample$opinion))
    r$EopClustering1[row] <- mean(sample$opClustering1)
    r$EopClustering2[row] <- mean(sample$opClustering2)
    r$EopClustering3[row] <- mean(sample$opClustering3)
    r$EopAlignment1[row] <- mean(abs(sample$opAlignment1))
    r$EopAlignment2[row] <- mean(abs(sample$opAlignment2))
    r$EopAlignment3[row] <- mean(abs(sample$opAlignment3))
  }
  for(d in c(1, 9, 10)){
    print(paste("plotting for district nr", d))
    t <- subset(r, r$wijk == d)
    par(mfrow=c(3,3))
    plot(t$opClustering1, t$EopClustering1, pch=20, cex=0.01)
    plot(t$opClustering2, t$EopClustering2, pch=20, cex=0.01)
    plot(t$opClustering3, t$EopClustering3, pch=20, cex=0.01)
    plot(t$opAlignment1, t$EopAlignment1, pch=20, cex=0.01)
    plot(t$opAlignment2, t$EopAlignment2, pch=20, cex=0.01)
    plot(t$opAlignment3, t$EopAlignment3, pch=20, cex=0.01)
    plot(t$meanAbsOpinion, t$EmeanAbsOpinion, pch=20, cex=0.01)
    mtext(districtsNames[d], side = 3, line = -2, outer = TRUE)
  }
  
  # Bootstrapping
  files <- list.files(path = "./simOutput/peregrine/")
  load(paste0("./simOutput/peregrine/", files[1]))
  w <- simW[[1]]
  
  estimates <- c()
  for (n in 1:200){
    sampleI <- sample(nrow(w), 400)
    sample <- data.frame()
    for (i in sampleI){
      sample <- rbind(sample, w[i,])
    }
    estimates[n] <- mean(abs(sample$opAlignment1))
  }
  hist(estimates, breaks=10, xlim=c(0,max(estimates)), main="op alignment (s=10, n=400)")
  abline(v=mean(abs(w$opAlignment1)), col="red")
  
  estimates <- c()
  for (n in 1:200){
    sampleI <- sample(nrow(w), 400)
    sample <- data.frame()
    for (i in sampleI){
      sample <- rbind(sample, w[i,])
    }
    estimates[n] <- mean(sample$opClustering3)
  }
  
  hist(estimates, breaks=10, xlim=c(0,max(estimates)), main="op clustering (s=1000, n=400)")
  abline(v=mean(w$opClustering3), col="red")
} ##




save(
  r, rr, ri, rri,
  file="./oldRunsResults.RData"
)

load(file="./oldRunsResults.RData")




# ==============================================================================
# Analysis of the simulation outcomes.

# loading results. We call r the macro-level, ri the micro.
rm (list = ls( ))
load("./cityData/geodata_Rotterdam.RData")
load("./simOutput/completeDataset.RDATA")
require(ggplot2)
r <- GsimResults
#ri <- GsimAgents

# For the micro we add the macro-level information about the simulation run that
# the agent was in.
ri <- base::merge(
  x = r[,c(1:4,6)],
  y = GsimAgents,
  by = "seed"
)
if(FALSE){
  rri$pol <- NA
  for (i in 1:nrow(rri)){
    rri$pol[i] <- r[r$seed==rri$seed[i],]$polarizationIndex
  }
}

# Free up some RAM
rm(GsimResults, GsimAgents)

# The segregation scores need to be fixed, and we do that now.
# We take the additive inverse of the raw scores for segregation
r$grSegregation1 <- r$grSegregation1 * -1
r$grSegregation2 <- r$grSegregation2 * -1
r$grSegregation3 <- r$grSegregation3 * -1
ri$grSegregation1 <- ri$grSegregation1 * -1
ri$grSegregation2 <- ri$grSegregation2 * -1
ri$grSegregation3 <- ri$grSegregation3 * -1

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


# Function to plot violins.
# version 1 
if(FALSE){
ggplot(rr, aes(factor(round(rr$NATexpOutgr2, digits=3)), rr$polarizationIndex))+
  ylab("Polarization index") +
   # geom_smooth() + ############################################
  geom_violin(
    scale = "width",
    draw_quantiles = 0.5,
    aes(fill=rr$wkName, color=rr$wkName)#, color="green")
  ) + xlab("Exposure to outgroup (s=100)") +
  #scale_x_discrete(limits = sort(unique(rr$NATexpOutgr2))) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    #axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )
}
# version 2 (preferred)
districtViolins <- function(
  data = rr,
  depVar, 
  indepVar, 
  depVarLabel, 
  indepVarLabel
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
  ggplot(data, aes(factor(data[,indepVar]), data[,depVar]))+
    ylab(depVarLabel) +
    #geom_smooth() + ###########################################################
    geom_violin(
      fill="gray",
      scale = "width",
      draw_quantiles = 0.5#,
      #aes(fill="gray")#, color="#ff8800")#, color="green")
    ) +  #stat_summary(fun.y=mean, geom="point", shape=4, size=2) + 
    ggtitle(paste("Districts ordered by\n", indepVarLabel)) +
    theme(
      plot.margin=unit(c(0,0,0,40),"pt"),
      plot.title = element_text(hjust=0.5),
      axis.text.x = element_text(angle = 30, hjust = 1),
      axis.title.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    scale_x_discrete(labels=labs)
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
      axis.text.x = element_text(angle = 30, hjust = 1),
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
  ggplot(data, 
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
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    )
}


test <- subset(
  ri,
  ( 
    ri$distanceDecay == 2 &#2 &
      ri$initialOpinionDistribution == "beta" &
      ri$H == 0.6
  )
)
par(mfrow=c(4,3))
for (d in 1:12){
  #for (d in c(1,2,11)){
  #temp <- subset(rri, rri$wijk == d)
  temp <- subset(test, test$wijk == d)
  p <-violinPlotOld(
    data = temp,
    depVar = temp$opAlignment2,
    #depVar = abs(temp$opAlignment2),
    indepVar = temp$exposureOutgroup,
    panels = temp$group,
    #bins=5,
    #depVarLabel = "Alignment (agent level, s=100)",
    depVarLabel = "Alignment (raw scores agent level, s=100)",
    indepVarLabel = "Outgroup exposure (s=100)"
  ) + facet_wrap(temp$group, labeller = as_labeller(c(
    "1"="Group 1: non-western",
    "-1" = "Group -1: western"
  ))) + ggtitle(districtsNames[d])
  print(p)
}

par(mfrow=c(2,1))
hist(subset(rri, rri$group==1)$exposureOutgroup, xlim=c(0,1))
hist(subset(rri, rri$group==-1)$exposureOutgroup, xlim=c(0,1))

dia <- subset(di, di$group == -1)
dib <- subset(di, di$group == 1)
par(mfrow=c(2,1))
ggplot(
  dia,
  aes(
    x=cut_interval(x=dia$expOutgr2, length=0.1),
    y=abs(dia$opAlignment2))
) + 
  geom_boxplot() + 
  xlab ("Exposure to outgroup (s=100)") +
  ylab ("Alignment (s=100)") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )
ggplot(
  dib,
  aes(
    x=cut_interval(x=dib$expOutgr2, length=0.1),
    y=abs(dib$opAlignment2))
) + 
  geom_boxplot() + 
  xlab ("Exposure to outgroup (s=100)") +
  ylab ("Alignment (s=100)") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )

violinPlotOld(
  data = di,
  depVar = abs(di$opAlignment2),
  indepVar = di$expOutgr2,
  panels = di$group,
  bins=7,
  depVarLabel = "Alignment (agent level, s=100)",
  indepVarLabel = "Outgroup exposure (s=100)"
) + facet_wrap(di$group, labeller = as_labeller(c(
  "1"="Group 1: non-western",
  "-1" = "Group -1: western"
)))





violinPlotOld(
  depVar = log10(rri$timeFirstPol),
  indepVar = rri$exposureOutgroup,
  panels = rri$wijk,
  depVarLabel = "Time of first extremization (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)",
  bins=5
) + facet_wrap(rri$wijk, labeller = as_labeller(c(
  "1" = "Stadscentrum",
  "2" = "Delfshaven",
  "3" = "Overschie",
  "4" = "Noord",
  "5" = "Hillegersberg-Schiebroek",
  "6" = "Kralingen-Crooswijk",
  "7" = "Feijenoord",
  "8" = "IJsselmonde",
  "9" = "Pernis",
  "10" = "Prins Alexander",
  "11" = "Charlois",
  "12" = "Hoogvliet"
)))



#head(  cut(di$expOutgr2, breaks=(c(0:10)/10))   )
#table(  as.numeric(cut(di$expOutgr2, breaks=(c(0:10)/10)))   )

# Exploring single runs from dataset ab
load(file="./ab.RData")

a[[1]][[200]]$polarizationIndex
b[[1]][[200]]$polarizationIndex
c[[1]][[50]]$polarizationIndex
d[[1]]$polarizationIndex
e[[1]]$polarizationIndex

abs(a[[1]][[200]]$meanOpinionG1 - a[[1]][[200]]$meanOpinionG2)
abs(b[[1]][[200]]$meanOpinionG1 - b[[1]][[200]]$meanOpinionG2)
abs(c[[1]][[50]]$meanOpinionG1 - c[[1]][[50]]$meanOpinionG2)
abs(d[[1]]$meanOpinionG1 - d[[1]]$meanOpinionG2)
abs(e[[1]]$meanOpinionG1 - e[[1]]$meanOpinionG2)
# polarization and strong global alignment in all three runs

ai <- a[[2]][[200]]
bi <- b[[2]][[200]]
ci <- c[[2]][[50]]
di <- d[[2]]
ei <- e[[2]]
plot(di$timeFirstExtr, di$nIntFirstExtr)

table(ai$nIntFirstExtr < ai$timeFirstExtr)
table(bi$nIntFirstExtr < bi$timeFirstExtr)
table(ci$nIntFirstExtr < ci$timeFirstExtr)
table(di$nIntFirstExtr < di$timeFirstExtr)
table(ei$nIntFirstExtr < ei$timeFirstExtr)

# time first extremization
violinPlot(
  data = ai,
  depVar = log10(ai$timeFirstExtr),
  indepVar = ai$expOutgr2,
  depVarLabel = "Time of first extremization (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)"#,
)
violinPlot(
  data = bi,
  depVar = log10(bi$timeFirstExtr),
  indepVar = bi$expOutgr2,
  depVarLabel = "Time of first extremization (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)"#,
)
violinPlot(
  data = ci,
  depVar = log10(ci$timeFirstExtr),
  indepVar = ci$expOutgr2,
  depVarLabel = "Time of first extremization (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)"#,
)
violinPlot(
  data = di,
  depVar = log10(di$timeFirstExtr),
  indepVar = di$expOutgr2,
  depVarLabel = "Time of first extremization (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)"#,
)
violinPlot(
  data = ei,
  depVar = log10(ei$timeFirstExtr),
  indepVar = ei$expOutgr2,
  depVarLabel = "Time of first extremization (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)"#,
)

# N interatctions before first extremization
violinPlot(
  data = ai,
  depVar = log10(ai$nIntFirstExtr),
  indepVar = ai$expOutgr2,
  depVarLabel = "N int. before first extr. (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)",
  fill = "purple"
)
violinPlot(
  data = bi,
  depVar = log10(bi$nIntFirstExtr),
  indepVar = bi$expOutgr2,
  depVarLabel = "N int. before first extr. (log_10))",
  indepVarLabel = "Outgroup exposure (s=100)",
  fill = "purple"
)
violinPlot(
  data = ci,
  depVar = log10(ci$nIntFirstExtr),
  indepVar = ci$expOutgr2,
  depVarLabel = "N int. before first extr. (log_10))",
  indepVarLabel = "Outgroup exposure (s=100)",
  fill = "purple"
)
violinPlot(
  data = di,
  depVar = log10(di$nIntFirstExtr),
  indepVar = di$expOutgr2,
  depVarLabel = "N int. before first extr. (log_10))",
  indepVarLabel = "Outgroup exposure (s=100)",
  fill = "purple"
)

#by group
violinPlot(
  data = ai,
  depVar = log10(ai$timeFirstExtr),
  indepVar = ai$expOutgr2,
  panels = ai$group,
  depVarLabel = "Time of first extremization (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)"
) + facet_wrap(ai$group, labeller = as_labeller(c(
  "1"="Group 1: non-western residents",
  "-1" = "Group -1: western residents"
)))
violinPlot(
  data = bi,
  depVar = log10(bi$timeFirstExtr),
  indepVar = bi$expOutgr2,
  panels = bi$group,
  depVarLabel = "Time of first extremization (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)"
) + facet_wrap(bi$group, labeller = as_labeller(c(
  "1"="Group 1: non-western residents",
  "-1" = "Group -1: western residents"
)))
violinPlot(
  data = ci,
  depVar = log10(ci$timeFirstExtr),
  indepVar = ci$expOutgr2,
  panels = ci$group,
  depVarLabel = "Time of first extremization (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)"
) + facet_wrap(ci$group, labeller = as_labeller(c(
  "1"="Group 1: non-western residents",
  "-1" = "Group -1: western residents"
)))
violinPlot(
  data = di,
  depVar = log10(di$timeFirstExtr),
  indepVar = di$expOutgr2,
  panels = di$group,
  depVarLabel = "Time of first extremization (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)"
) + facet_wrap(di$group, labeller = as_labeller(c(
  "1"="Group 1: non-western residents",
  "-1" = "Group -1: western residents"
)))
violinPlot(
  data = ai,
  depVar = log10(ai$nIntFirstExtr),
  indepVar = ai$expOutgr2,
  panels = ai$group,
  depVarLabel = "N int. before first extr. (log_10))",
  indepVarLabel = "Outgroup exposure (s=100)",
  fill = "purple"
) + facet_wrap(ai$group, labeller = as_labeller(c(
  "1"="Group 1: non-western residents",
  "-1" = "Group -1: western residents"
)))
violinPlot(
  data = bi,
  depVar = log10(bi$nIntFirstExtr),
  indepVar = bi$expOutgr2,
  panels = bi$group,
  depVarLabel = "N int. before first extr. (log_10))",
  indepVarLabel = "Outgroup exposure (s=100)",
  fill = "purple"
) + facet_wrap(bi$group, labeller = as_labeller(c(
  "1"="Group 1: non-western residents",
  "-1" = "Group -1: western residents"
)))
violinPlot(
  data = ci,
  depVar = log10(ci$nIntFirstExtr),
  indepVar = ci$expOutgr2,
  panels = ci$group,
  depVarLabel = "N int. before first extr. (log_10))",
  indepVarLabel = "Outgroup exposure (s=100)",
  fill = "purple"
) + facet_wrap(ci$group, labeller = as_labeller(c(
  "1"="Group 1: non-western residents",
  "-1" = "Group -1: western residents"
)))
violinPlot(
  data = di,
  depVar = log10(di$nIntFirstExtr),
  indepVar = di$expOutgr2,
  panels = di$group,
  depVarLabel = "N int. before first extr. (log_10))",
  indepVarLabel = "Outgroup exposure (s=100)",
  fill = "purple"
) + facet_wrap(di$group, labeller = as_labeller(c(
  "1"="Group 1: non-western residents",
  "-1" = "Group -1: western residents"
)))





#

require(plyr)
require(ggplot2)
#rri$plot <- mapvalues(rri$group, c(-1,1), c("group -1", "b"))
#panelLabel <- labeller(c(
#    "1" = "group 1",
#    "-1" = "group -1"
#  )
#)
violinPlot(
  depVar = log10(rri$timeFirstPol),
  indepVar = rri$exposureOutgroup,
  depVarLabel = "Time of first extremization (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)"#,
)

violinPlot(
  depVar = log10(rri$timeFirstPol),
  indepVar = rri$exposureOutgroup,
  panels = rri$group,
  depVarLabel = "Time of first extremization (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)"
) + facet_wrap(rri$group, labeller = as_labeller(c(
  "1"="Group 1: non-western residents",
  "-1" = "Group -1: western residents"
))) 


violinPlot(
  depVar = log10(rri$timeFirstPol),
  indepVar = rri$exposureOutgroup,
  panels = rri$wijk,
  depVarLabel = "Time of first extremization (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)",
  bins=5
) + facet_wrap(rri$wijk, labeller = as_labeller(c(
    "1" = "Stadscentrum",
    "2" = "Delfshaven",
    "3" = "Overschie",
    "4" = "Noord",
    "5" = "Hillegersberg-Schiebroek",
    "6" = "Kralingen-Crooswijk",
    "7" = "Feijenoord",
    "8" = "IJsselmonde",
    "9" = "Pernis",
    "10" = "Prins Alexander",
    "11" = "Charlois",
    "12" = "Hoogvliet"
)))


load(file= "oldRunsResults.RData")

print(districtViolins(
  depVar="opNetAli2",
  indepVar="expOutgr2",
  depVarLabel="Average local alignment (s=100)",
  indepVarLabel="outgroup exposure (s=100)"))

print(districtViolins(
  depVar="opAlignment2",
  indepVar="expOutgr2",
  depVarLabel="Average local alignment (s=100)",
  indepVarLabel="outgroup exposure (s=100)"))




# old measures of alignment

rr$oldAlignment <- rr$varOpinionGlobal - mean(c(rr$varOpinionG1, rr$varOpinionG2))
hist(rr$oldAlignment)
districtViolins(
  depVar="oldAlignment",
  indepVar="expOutgr2",
  depVarLabel="Alignment (difference between\nglobal- and group-variance",
  indepVarLabel="Outgroup exposure (s=100)")

rr$intuitiveAlignment <- abs(rr$varOpinionG1 - rr$varOpinionG2)
hist(rr$intuitiveAlignment)
districtViolins(
  depVar="intuitiveAlignment",
  indepVar="expOutgr2",
  depVarLabel="Alignment (abs. difference between\nthe average attitude of the two groups",
  indepVarLabel="Outgroup exposure (s=100)")


# We find and save separately the runs that were ran with the 'baseline'
# parameter configuration.
rr <- subset(
  r,
  ( 
    r$distanceDecay == 2 &
    r$initialOpinionDistribution == "groupBias" &
    r$H == 0.6
  )
)
for(i in 1:nrow(rr))(rr$wkName[i]<-districtsNames[rr$wijk[i]])


rri <- subset(
  ri,
  ( 
    ri$distanceDecay == 2 &#2 &
      ri$initialOpinionDistribution == "groupBias" &
      ri$H == 0.6
  )
)




# 1) segregation -> polarization

# Preliminary question: did we get some polarization at all?
# Here we use histograms to plot the distribution of the polarization index
# at the end of the simulation run. We look into the baseline configuration (rr)
# as well as all configurations together (r)
for(i in 1:nrow(r)){if(r$polarizationIndex[i]>1){r$polarizationIndex[i]<-1}}
par(mar=c(4,1,1,0.1))
hist(
  r$polarizationIndex, breaks=10,
  main="All parameter configurations",
  xlab="Polarization index", yaxt="n", ylab="", yaxt='n',
  col="black", border="white"
)
abline(v=0.3, col= "red")
title(ylab="Relative frequency", line=0)

for(i in 1:nrow(rr)){if(rr$polarizationIndex[i]>1){rr$polarizationIndex[i]<-1}}
par(mar=c(4,1,1,0.1))
hist(
  rr$polarizationIndex, breaks=10,
  main="Baseline parameter configuration",
  xlab="Polarization index", yaxt="n", ylab="", yaxt='n',
  col="black", border="white"
)
abline(v=0.3, col= "red")
title(ylab="Relative frequency", line=0)

for(i in 1:nrow(rr)){if(rr$polarizationIndex[i]>1){rr$polarizationIndex[i]<-1}}
par(mar=c(4,1,1,0.1))
hist(
  rri$opinion, breaks=10,
  main="",
  xlab="Attitude", yaxt="n", ylab="", yaxt='n',
  col="black", border="black"
)
title(ylab="Relative frequency", line=0)

par(mar=c(4,1,1,0.1))
hist(
  rri$timeFirstPol, breaks=10,
  main="",
  xlab="Time of first extremization", yaxt="n", ylab="", yaxt='n',
  col="black", border="black"
)
title(ylab="Relative frequency", line=0)

par(mar=c(4,1,1,0.1))
hist(
  rri$exposureOutgroup, breaks=10,
  main="",
  xlab="Exposure to outgroup (s=100)", yaxt="n", ylab="", yaxt='n',
  col="black", border="white"
)
title(ylab="Relative frequency", line=0)



# Plotting micro-level
#
# First we need to fetch the resources and polygon data.
load("./cityData/importedCBS.RData")

cit <- subset(nedb,nedb$WK_CODE %in% districtsList)
cit@data$WK_CODE <- droplevels(cit@data$WK_CODE)
cit <- spTransform(cit,CRS("+proj=longlat +ellps=WGS84"))

a <- "gray80"
b <- "gray60"
c <- "gray40"
d <- "gray20"
color <- c(a,b,c,d,b,c,b,a,c,d,d,a)
rm(a,b,c,d)

cit@data$COLOR <- NA
for (i in 1:length(cit@data$COLOR)){
  c <- citySummary
  c$color <- color
  l <- subset(c, c$WK_CODE == cit@data$WK_CODE[i])
  cit@data$COLOR[i] <- l$color
  #cit@polygons$COLOR[i] <- l$color
}
par(mar=c(0,0,0,0))
plot(cit, col=cit$COLOR, border=cit$COLOR)

ifelse(
  district == 0,
  dat <- subset(cbs100_rot, cbs100_rot$WK_CODE %in% districtsList),
  dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[district])
)
dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]


dat$tfp <- dat$absop <- #pol <- 
  dat$clus1 <- dat$clus2 <- dat$clus3 <-
  dat$li1 <- dat$ali2 <- dat$ali3 <-
  dat$exp1 <- dat$exp2 <- dat$exp3 <- NA
for (l in 1:length(dat)){ # for every cell
  cell <- subset(w, w$location == dat$OBJECTID[l])
  cell$timeFirstPol[is.na(cell$timeFirstPol)] <- 201
  dat$tfp[l] <- mean(cell$timeFirstPol, na.rm=TRUE)
  dat$absop[l] <- mean(abs(cell$opinion), na.rm=TRUE)
  #pol[l] <- mean(cell$pol, na.rm=TRUE)
  dat$clus1[l] <- mean(cell$opClustering1, na.rm=TRUE)
  dat$clus2[l] <- mean(cell$opClustering2, na.rm=TRUE)
  dat$clus3[l] <- mean(cell$opClustering3, na.rm=TRUE)
  dat$ali1[l] <- mean(abs(cell$opAlignment1), na.rm=TRUE)
  dat$ali2[l] <- mean(abs(cell$opAlignment2), na.rm=TRUE)
  dat$ali3[l] <- mean(abs(cell$opAlignment3), na.rm=TRUE)
  dat$exp2[l] <- mean(cell$exposureOutgroup, na.rm=TRUE)
}

cPalette <- colorRampPalette(c("blue", "red")) # low, high
cPaletteTrio <-colorRampPalette(c("blue", "white", "red"))
cPaletteReverse <- colorRampPalette(c("red", "blue")) # low, high
x <- log(dat$tfp +1, 10) 
dat$colo <- cPaletteReverse(20)[as.numeric(cut(x, breaks=20))]
par(mar=c(0,0,1,0))
plot(cit, col=cit$COLOR, border=cit$COLOR, main="Time of first polarization")
points(dat$l, dat$j, pch=16, col = dat$colo, cex=0.5)
par(mar=c(5,5,5,5))






# Time of first polarization
ggplot(
  rri,
  aes(
    x=cut_interval(x=rri$exposureOutgroup, length=0.1),
    y=log10(rri$timeFirstPol+1))
  ) + 
  geom_boxplot() + 
  xlab ("Exposure to outgroup (s=100)") +
  ylab ("Time of first extremization (log_10)") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )

require(dplyr)
require(tidyr)
ggplot(
  rri,
  aes(
    x=cut_interval(x=rri$exposureOutgroup, length=0.2),
    y=log10(rri$timeFirstPol+1))
) + 
  geom_boxplot() + 
  facet_wrap(rri$wijk) +
  xlab ("Exposure to outgroup (s=100)") +
  ylab ("Time of first extremization (log_10)") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    axis.line = element_line(colour = "black")
  )


# alignment

ggplot(
  rri,
  aes(
    x=cut_interval(x=rri$exposureOutgroup, length=0.1),
    y=abs(rri$opAlignment2))
) + 
  geom_boxplot() + 
  xlab ("Exposure to outgroup (s=100)") +
  ylab ("Alignment (s=100)") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )



# Absolute opinion
dat$colo <- cPalette(20)[as.numeric(cut(dat$absop, breaks=20))]
par(mar=c(0,0,1,0))
plot(cit, col=cit$COLOR, border=cit$COLOR, main="Absolute opinion")
points(dat$l, dat$j, pch=16, col = dat$colo, cex=0.5)
par(mar=c(5,5,5,5))

# Opinion clustering
x<-dat$clus2
for(i in 1:length(x)){
  if(is.na(x[i])){x[i] <- 0}
  if(x[i]>4){x[i]<-4}
}
dat$colo <- cPaletteTrio(20)[as.numeric(cut(x, breaks=seq(-4,4,8/20)))]
par(mar=c(0,0,1,0))
plot(cit, col=cit$COLOR, border=cit$COLOR, main="Opinion clustering (s=100)")
points(dat$l, dat$j, pch=16, col = dat$colo, cex=0.5)
par(mar=c(5,5,5,5))

heatmap(
  x = rri$exposureOutgroup,#rbeta(1000,3,9),
  y = x,#rbeta(1000,9,3),
  type = "proportion",
  bins = 11,
  xlab = "Exposure to outgroup (s=100)",
  ylab = "Opinion clustering (s=100)",
  legend = FALSE
)

# Opinion alignment
x<-dat$ali2
for(i in 1:length(x)){
  if(is.na(x[i])){x[i] <- 0}
  if(x[i]>4){x[i]<-4}
}
dat$colo <- cPalette(20)[as.numeric(cut(x, breaks=20))]
par(mar=c(0,0,1,0))
plot(cit, col=cit$COLOR, border=cit$COLOR, main="Group-attitude alignment (s=100)")
points(dat$l, dat$j, pch=16, col = dat$colo, cex=0.5)
par(mar=c(5,5,5,5))

heatmap(
  x = rri$exposureOutgroup,#rbeta(1000,3,9),
  y = x,#rbeta(1000,9,3),
  type = "count",
  bins = 11,
  xlab = "Exposure to outgroup (s=100)",
  ylab = "Group-attitude alignment (s=100)",
  legend = FALSE
)


# outgroup exposure
x<-dat$exp2
for(i in 1:length(x)){if(is.na(x[i])){x[i] <- 0}}
dat$colo <- cPalette(20)[as.numeric(cut(x, breaks=20))]
par(mar=c(0,0,1,0))
plot(cit, col=cit$COLOR, border=cit$COLOR, main="Outgroup exposure (s=100)")
points(dat$l, dat$j, pch=16, col = dat$colo, cex=0.5)
par(mar=c(5,5,5,5))

# proportion nwa
x<-dat$pnwal2014
#for(i in 1:length(x)){if(is.na(x[i])){x[i] <- 0}}
dat$colo <- cPalette(20)[as.numeric(cut(x, breaks=20))]
par(mar=c(0,0,1,0))
plot(cit, col=cit$COLOR, border=cit$COLOR, main="Proportion non-western residents")
points(dat$l, dat$j, pch=16, col = dat$colo, cex=0.5)
par(mar=c(5,5,5,5))





# clustering - alignment
x<-dat$clus2
for(i in 1:length(x)){
  if(is.na(x[i])){x[i] <- 0}
  if(x[i]>4){x[i]<-4}
}
y<-dat$ali2
for(i in 1:length(y)){
  if(is.na(y[i])){y[i] <- 0}
  if(y[i]>4){y[i]<-4}
}

heatmap(
  x = x,#rbeta(1000,3,9),
  y = y,#rbeta(1000,9,3),
  type = "prop",
  bins = 11,
  xlab = "Attitude clustering (s=100)",
  ylab = "Group-attitude alignment (s=100)",
  legend = FALSE
)



heatmap(
  x = dat$pnwal2014,#rbeta(1000,3,9),
  y = y,#rbeta(1000,9,3),
  type = "prop",
  bins = 11,
  xlab = "Attitude clustering (s=100)",
  ylab = "Group-attitude alignment (s=100)",
  legend = FALSE
)



print(districtViolins(
  depVar="opAlignment1",
  indepVar="expOutgr1",
  depVarLabel="Group-attitude alignment (s=100)",
  indepVarLabel="NWA outgroup exposure (s=100)"))




#plot(rri$exposureOutgroup, log(rri$timeFirstPol + 1,10))
y <- scale(log(rri$timeFirstPol + 1,10))
y <- y + min(y)


for(i in 1:nrow(rri)){if(is.na(rri$timeFirstPol[i])){rri$timeFirstPol[i]<-201}}
 
heatmap(
  x = rri$exposureOutgroup,#rbeta(1000,3,9),
  y = log(rri$timeFirstPol + 1,10),#rbeta(1000,9,3),
  type = "count",
  bins = 11,
  xlab = "Exposure to outgroup (s=100)",
  ylab = "Time of first polarization (log_10)",
  legend = FALSE
)

for(i in 1:nrow(rri)){if(is.na(rri$grSegregation2[i])){rri$grSegregation2[i]<-0}}
heatmap(
  x = rri$grSegregation2,#rbeta(1000,3,9),
  y = log(rri$timeFirstPol + 1,10),#rbeta(1000,9,3),
  type = "proportion",
  bins = 10,
  xlab = "Group segregation (s=100)",
  ylab = "Time of first polarization (log_10)",
  legend = FALSE
)


x <- log(dat$tfp +1, 10) 
dat$colo <- cPalette(20)[as.numeric(cut(dat$tfp, breaks=20))]
par(mar=c(0,0,1,0))
plot(
  cit, col=cit$COLOR, border=cit$COLOR,
  main="Time of first polarization"
)
points(dat$l, dat$j, pch=20, col = dat$colo, cex=1)


par(mar=c(5,5,5,5))


dat$colo <- cPalette(20)[as.numeric(cut(dat$tfp, breaks=20))]
plot(dat$l, dat$j, cex=2, pch=15, col=dat$colodistrict, main = 
       paste("time of first polarization \nin", placeName))
points(dat$l, dat$j, pch=20, col = dat$colo, cex=1)




dat$tfp <- tfp
dat$absop <- absop
dat$clus <- clus
dat$clus[dat$clus > 1.1] <- 1.1 #####
dat$clus[dat$clus < 0] <- 0############
#dat$pol <-pol
#dat$cclus <- cclus
#dat$aali <- aali
dat$ali <- ali
dat$ali[dat$ali > 1.1] <-1 #############
dat$exp <- exp


#dat$colodistrict <- cPaletteD(12)[as.numeric(cut(as.numeric(dat$WK_CODE), breaks=12))]







# Time first polarization
dat$colo <- cPalette(20)[as.numeric(cut(dat$tfp, breaks=20))]
plot(dat$l, dat$j, cex=2, pch=15, col=dat$colodistrict, main = 
       paste("time of first polarization \nin", placeName))
points(dat$l, dat$j, pch=20, col = dat$colo, cex=1)






# Plotting micro-level
#
# First we need to fetch the resources and polygon data.
load("./cityData/importedCBS.RData")
require(maptools)
require(rgeos)
require(tidyverse)
require(broom)
 





  
  



















district <- 0 # if 0, whole Rotterdam
placeName <- NA
ifelse(district==0, placeName<-"Rotterdam", placeName<-districtsNames[district])
cPalette <- colorRampPalette(c("blue", "red")) # low, high
cPaletteD <- colorRampPalette(c("beige","orange", "yellow","cyan", "green", "azure3", "black")) # low, high
#wijk <- rr[example,]$wijk
#load(paste0("./simOutput/peregrine/", rr[example,]$fileName))
#w <- worldList[[rr$indexParameters[example]]]
#w <- subset(rri, rri$wijk==example)
w <- rri
#w <- s
#safe <- w

polarizedSeeds <- subset(r, r$polarizationIndex > 0.7)$seed
w <- w[w$seed %in% polarizedSeeds,]


#w$cclus <- abs(w$opinion) * w$opClustering2
#w$aali <- abs(w$opinion) * w$opAlignment2
# plotting all districts together
#
ifelse(
  district == 0,
  dat <- subset(cbs100_rot, cbs100_rot$WK_CODE %in% districtsList),
  dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[district])
)
dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]


tfp <- absop <- pol <- clus <- ali <- exp <- c()
for (l in 1:length(dat)){ # for every cell
  cell <- subset(w, w$location == dat$OBJECTID[l])
  cell$timeFirstPol[is.na(cell$timeFirstPol)] <- 201
  tfp[l] <- mean(cell$timeFirstPol, na.rm=TRUE)
  absop[l] <- mean(abs(cell$opinion), na.rm=TRUE)
  #pol[l] <- mean(cell$pol, na.rm=TRUE)
  clus[l] <- mean(cell$opClustering2, na.rm=TRUE)
  ali[l] <- mean(abs(cell$opAlignment2), na.rm=TRUE)
  exp[l] <- mean(cell$exposureOutgroup, na.rm=TRUE)
}

dat$tfp <- tfp
dat$absop <- absop
dat$clus <- clus
dat$clus[dat$clus > 1.1] <- 1.1 #####
dat$clus[dat$clus < 0] <- 0############
#dat$pol <-pol
#dat$cclus <- cclus
#dat$aali <- aali
dat$ali <- ali
dat$ali[dat$ali > 1.1] <-1 #############
dat$exp <- exp


dat$colodistrict <- cPaletteD(12)[as.numeric(cut(as.numeric(dat$WK_CODE), breaks=12))]







# Time first polarization
dat$colo <- cPalette(20)[as.numeric(cut(dat$tfp, breaks=20))]
plot(dat$l, dat$j, cex=2, pch=15, col=dat$colodistrict, main = 
       paste("time of first polarization \nin", placeName))
points(dat$l, dat$j, pch=20, col = dat$colo, cex=1)


#plot(dat$l, dat$j, pch=20, cex=0.5)
       #paste("time of first polarization \n", districtsNames[example]))


# Abs opinion
dat$colo <- cPalette(10)[as.numeric(cut(dat$absop, breaks=10))]
plot(dat$l, dat$j, pch=20, col = dat$colo, cex=1, main = 
       paste("abs. opinion \nin", placeName))
       #paste("absolute opinion\n", districtsNames[example]))


# op clustering
dat$colo <- cPalette(10)[as.numeric(cut(dat$clus, breaks=10))]
plot(dat$l, dat$j, pch=15, col = dat$colodistrict, cex=2, main = 
       paste("op. clustering (s=100)\nin", placeName))

points(dat$l, dat$j, cex=1, pch=20, col=dat$colo)


# op alignment
dat$colo <- cPalette(20)[as.numeric(cut(dat$ali, breaks=20))]
plot(dat$l, dat$j, pch=15, col = dat$colodistrict, cex=2, main = 
       paste("op. alignment (s=100)\nin", placeName))
points(dat$l, dat$j, pch=20, col = dat$colo, cex=1)


# Abs opinion * clustering
dat$colo <- cPalette(10)[as.numeric(cut(dat$absop * dat$clus, breaks=10))]
plot(dat$l, dat$j, pch=20, col = dat$colo, cex=0.5, main = 
       paste("abs. op * clustering (s=100)\nin", placeName))
       #paste("alignment (s=100)\n", districtsNames[example]))


# Abs opinion * clustering * alignment
dat$colo <- cPalette(10)[as.numeric(cut(dat$absop * dat$clus * dat$clus, breaks=10))]
plot(dat$l, dat$j, pch=20, col = dat$colo, cex=0.5, main = 
       paste("abs. op * clust. * align. (s=100)\nin", placeName))
#paste("alignment (s=100)\n", districtsNames[example]))


# Outgroup exposure
dat$colo <- cPalette(10)[as.numeric(cut(dat$exp, breaks=10))]
plot(dat$l, dat$j, pch=20, col = dat$colo, cex=0.5, main = 
       paste("exposure outgroup\n", placeName))


# proportion non-western
dat$colo <- cPalette(10)[as.numeric(cut(dat$pnwal2014, breaks=10))]
plot(dat$l, dat$j, pch=20, col = dat$colo, cex=0.5, main = 
       paste("proportion non-western\nin", placeName))



















# differences between districts
#
# 1) polarization
#boxplot(rr$polarizationIndex ~ rr$wkName)
ggplot(rr, aes(factor(rr$wkName), rr$polarizationIndex)) +
  ylab("Polarization index") +
  geom_violin(scale = "width",draw_quantiles = 0.5, fill="#ff8800") +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )

# 2) Clustering
ggplot(rr, aes(factor(rr$wkName), rr$opNetClus2)) +
  ylab("Attitude clustering (s=100)") +
  geom_violin(scale = "width",draw_quantiles = 0.5, fill="#ff8800") +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )

# 3) alignment
ggplot(rr, aes(factor(rr$wkName), rr$opNetAli2)) +
  ylab("Alignment (s=100)") +
  geom_violin(scale = "width",draw_quantiles = 0.5, fill="#ff8800") +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )





rr <- subset(
  r,
  ( 
    r$distanceDecay == 2 &
      r$initialOpinionDistribution == "groupBias" &
      r$H == 0.6
  )
)
for(i in 1:nrow(rr))(rr$wkName[i]<-districtsNames[rr$wijk[i]])


print(districtViolins(
  depVar="polarizationIndex",
  indepVar="expOutgr1",
  depVarLabel="Polarization index",
  indepVarLabel="exposure to outgroup (s=10)"))

print(districtViolins(
  depVar="polarizationIndex",
  indepVar="expOutgr2",
  depVarLabel="Polarization index",
  indepVarLabel="exposure to outgroup (s=100)"))

print(districtViolins(
  depVar="polarizationIndex",
  indepVar="expOutgr3",
  depVarLabel="Polarization index",
  indepVarLabel="exposure to outgroup (s=1000)"))


print(districtViolins(
  depVar="polarizationIndex",
  indepVar="p_nwa",
  depVarLabel="Polarization index",
  indepVarLabel="proportion of non-western residents"))



print(districtViolins(
  depVar="opNetClus2",
  indepVar="expOutgr2",
  depVarLabel="Opinion clustering (s=100)",
  indepVarLabel="exposure to outgroup (s=100)"))

print(districtViolins(
  depVar="opNetAli2",
  indepVar="expOutgr2",
  depVarLabel="Average local alignment (s=100)",
  indepVarLabel="outgroup exposure (s=100)"))



# plotting micro stats










test<- subset(rr, rr$polarizationIndex > 0.7)

ggplot(test, aes(factor(test$wkName), test$opNetClus2)) +
  ylab("Attitude clustering (s=100)") +
  geom_violin(scale = "width",draw_quantiles = 0.5, fill="#ff8800") +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )





# Regression tables =======================================================

rr <- subset(r,r$distanceDecay==2)

hist(r$polarizationIndex - mean(r$polarizationIndex))
hist(rr$polarizationIndex - mean(rr$polarizationIndex))

hist(r$opNetClus2 - mean(r$opNetClus2))
hist(rr$opNetClus2 - mean(rr$opNetClus2))

hist(r$opNetAli2 - mean(r$opNetAli2))
hist(rr$opNetAli2 - mean(rr$opNetAli2))


# 1) segr -> polarization

summary(lm(formula = rr$polarizationIndex ~ rr$expOutgr2 + 
     as.factor(rr$initialOpinionDistribution) + rr$H ))

summary(lm(formula = rr$polarizationIndex ~ rr$p_nwa +
     rr$grSegregation2 + rr$initialOpinionDistribution + rr$H ))

summary(lm(formula = rr$polarizationIndex ~ rr$p_nwa + 
     rr$grSegregation2 + rr$expOutgr2 +  rr$initialOpinionDistribution + rr$H ))

summary(lm(formula = rr$polarizationIndex ~ as.factor(rr$wijk)))


# 2) segr -> op clustering

summary(lm(formula = rr$opNetClus2 ~ rr$expOutgr2 + 
             as.factor(rr$initialOpinionDistribution) + rr$H ))

summary(lm(formula = rr$opNetClus2 ~ rr$p_nwa +
             rr$grSegregation2 + rr$initialOpinionDistribution + rr$H ))

summary(lm(formula = rr$opNetClus2 ~ rr$p_nwa + 
             rr$grSegregation2 + rr$expOutgr2 +  rr$initialOpinionDistribution + rr$H ))

summary(lm(formula = rr$opNetClus2 ~ as.factor(rr$wijk)))

# 3) segr -> alignment

summary(lm(formula = rr$opNetAli2 ~ rr$expOutgr2 + 
             as.factor(rr$initialOpinionDistribution) + rr$H ))

summary(lm(formula = rr$opNetAli2 ~ rr$p_nwa +
             rr$grSegregation2 + rr$initialOpinionDistribution + rr$H ))

summary(lm(formula = rr$opNetAli2 ~ rr$p_nwa + 
             rr$grSegregation2 + rr$expOutgr2 +  rr$initialOpinionDistribution + rr$H ))

summary(lm(formula = rr$opNetAli2 ~ as.factor(rr$wijk)))

fit<-lm(formula = rr$polarizationIndex ~ rr$p_nwa + 
          rr$grSegregation2 + rr$expOutgr2 +  rr$initialOpinionDistribution + rr$H )
par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(fit)  # Plot the model information




#rr <- subset(GsimResults, GsimResults$initialOpinionDistribu tion!="groupBias")
rr <- subset(rr, rr$distanceDecay == 2)


hist(subset(rr,rr$initialOpinionDistribution=="beta")$polarizationIndex)
hist(subset(rr,rr$initialOpinionDistribution=="uniform")$polarizationIndex)

# 1)
summary(lm(rr$meanAbsOpinion~rr$p_nwa))
summary(lm(rr$meanAbsOpinion~rr$p_nwa + rr$expOutgroup))
summary(lm(
  rr$meanAbsOpinion ~ rr$p_nwa + rr$expOutgroup +
    rr$initialOpinionDistribution +
    rr$H +
    rr$wijk
))
summary(lm(
  rr$polarizationIndex ~ rr$p_nwa + rr$expOutgroup +
    rr$initialOpinionDistribution +
    rr$H +
    rr$wijk
))

# 2)
summary(lm(
  rr$opClusteringA2 ~ rr$p_nwa + rr$expOutgroup +
    rr$initialOpinionDistribution +
    rr$H +
    rr$wijk
))

# 3)
# measuring baseline alignment at t=0
source("./script NI&PA.R")
load("./cityData/geodata_Rotterdam.RData")
d <- ali2 <- c()
for (r in 0:9){
  print(paste0("replication #", r + 1 ," of 10"))
  for(i in 1:12){
    d[12*r + i] <- i
    w <- worldList[[i]]
    w$opinion <- rbeta(nrow(w), 3, 3, ncp = 0)
    w$opinion <- w$opinion * 2 - 1
    dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[i])
    dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
    ops <- c(NA, length = length(dat))
    for (l in 1:length(dat)){ # for every cell
      cell <- subset(w, w$location == dat$OBJECTID[l])
      ops[l] <- mean(cell$opinion)
    }
      dat$opinion <- ops
      ali2[12*r + i] <- mean(abs(MoranI( # global alignment
        x = dat$pnwal2014,
        y = dat$opinion,
        proxmat = proximityList2[[i]],
        type = "local"
      )))
  }
}
test <- as.data.frame(cbind(d,ali2))
for(i in 1:nrow(rr)){
  rr$relativeAlignment[i] <- rr$opAlignment2[i] - mean(subset(test, test$d==rr$wijk[i])$ali2)
}
rr$absoluteAlignment <- rr$opAlignment2
summary(lm(
  rr$relativeAlignment ~ rr$p_nwa + rr$expOutgroup +
    rr$initialOpinionDistribution +
    rr$H +
    rr$wijk
))
summary(lm(
  rr$absoluteAlignment ~ rr$p_nwa + rr$expOutgroup +
    rr$initialOpinionDistribution +
    rr$H +
    rr$wijk
))

r1 <- subset(rr, rr$H == 0.6 & rr$initialOpinionDistribution == "beta")
r2 <- subset(rr, rr$H == 0.6 & rr$initialOpinionDistribution != "beta")
r3 <- subset(rr, rr$H == 0.9 & rr$initialOpinionDistribution == "beta")
r4 <- subset(rr, rr$H == 0.9 & rr$initialOpinionDistribution != "beta")

#0)
png(
  filename = "./outputGraphics/0.png",
  width = 1000,
  height = 800,
  units = "px"
)
par(mfrow=c(4,4))#,oma=c(0,6,0,0))
hist(r1$polarizationIndex, xlab="", ylab="Frequency", cex.lab=1.5, main="H=.6, beta op.distr.",cex.main=1.5, xlim=c(0,1))
mtext("op. polarization", side=1, line=5, cex=1)
hist(r2$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="H=.6, uniform op.distr.", cex.main=1.5, xlim=c(0,1))
mtext("op. polarization", side=1, line=5, cex=1)
hist(r3$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="H=.9, beta op.distr.", cex.main=1.5, xlim=c(0,1))
mtext("op. polarization", side=1, line=5, cex=1)
hist(r4$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="H=.9, uniform op.distr.", cex.main=1.5, xlim=c(0,1))
mtext("op. polarization", side=1, line=5, cex=1)
#op clustering
hist(r1$opClustering2, xlab="", ylab="Frequency", cex.lab=1.5, main="", xlim=c(-0.2,0.2))
mtext("op. clustering 2", side=1, line=3, cex=1)
hist(r2$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", xlim=c(-0.2,0.2))
mtext("op. clustering 2", side=1, line=3, cex=1)
hist(r3$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", xlim=c(-0.2,0.2))
mtext("op. clustering 2", side=1, line=3, cex=1)
hist(r4$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", xlim=c(-0.2,0.2))
#op alignment absolute
hist(r1$absoluteAlignment, xlab="", ylab="Frequency", cex.lab=1.5, main="", xlim=c(-0.3,0.3))
mtext("op. alignment 2", side=1, line=3, cex=1)
hist(r2$absoluteAlignment, xlab="", ylab="", cex.lab=1.5, main="", xlim=c(-0.3,0.3))
mtext("op. alignment 2", side=1, line=3, cex=1)
hist(r3$absoluteAlignment, xlab="", ylab="", cex.lab=1.5, main="", xlim=c(-0.3,0.3))
mtext("op. alignment 2", side=1, line=3, cex=1)
hist(r4$absoluteAlignment, xlab="", ylab="", cex.lab=1.5, main="", xlim=c(-0.3,0.3))
mtext("op. alignment 2", side=1, line=3, cex=1)
#op alignment relative
hist(r1$relativeAlignment, xlab="", ylab="Frequency", cex.lab=1.5, main="", xlim=c(-0.25,0.25))
mtext("relative op. alignment 2", side=1, line=3, cex=1)
hist(r2$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", xlim=c(-0.25,0.25))
mtext("relative op. alignment 2", side=1, line=3, cex=1)
hist(r3$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", xlim=c(-0.25,0.25))
mtext("relative op. alignment 2", side=1, line=3, cex=1)
hist(r4$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", xlim=c(-0.25,0.25))
mtext("relative op. alignment 2", side=1, line=3, cex=1)
dev.off()

#1a)
png(
  filename = "./outputGraphics/1a.png",
  width = 1000,
  height = 1000,
  units = "px"
)
par(mfrow=c(5,4))#,oma=c(0,6,0,0))
# exposure to outgroup
plot(r1$expOutgroup, r1$polarizationIndex, xlab="", ylab="op. polarization", cex.lab=1.5, main="H=.6, beta op.distr.", pch=".")
abline(lm(r1$polarizationIndex~r1$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r2$expOutgroup, r2$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="H=.6, uniform op.distr.", pch=".")
abline(lm(r2$polarizationIndex~r2$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r3$expOutgroup, r3$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="H=.9, beta op.distr.", pch=".")
abline(lm(r3$polarizationIndex~r3$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r4$expOutgroup, r4$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="H=.9, uniform op.distr.", pch=".")
abline(lm(r4$polarizationIndex~r4$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
# proportion non-western
plot(r1$p_nwa, r1$polarizationIndex, xlab="", ylab="op. polarization", cex.lab=1.5, main="", pch=".")
abline(lm(r1$polarizationIndex~r1$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r2$p_nwa, r2$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$polarizationIndex~r2$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r3$p_nwa, r3$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$polarizationIndex~r3$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r4$p_nwa, r4$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$polarizationIndex~r4$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 1)
plot(r1$grSegregation1, r1$polarizationIndex, xlab="", ylab="op. polarization", cex.lab=1.5, main="", pch=".")
abline(lm(r1$polarizationIndex~r1$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r2$grSegregation1, r2$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$polarizationIndex~r2$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r3$grSegregation1, r3$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$polarizationIndex~r3$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r4$grSegregation1, r4$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$polarizationIndex~r4$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 2)
plot(r1$grSegregation2, r1$polarizationIndex, xlab="", ylab="op. polarization", cex.lab=1.5, main="", pch=".")
abline(lm(r1$polarizationIndex~r1$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r2$grSegregation2, r2$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$polarizationIndex~r2$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r3$grSegregation2, r3$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$polarizationIndex~r3$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r4$grSegregation2, r4$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$polarizationIndex~r4$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 3)
plot(r1$grSegregation3, r1$polarizationIndex, xlab="", ylab="op. polarization", cex.lab=1.5, main="", pch=".")
abline(lm(r1$polarizationIndex~r1$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r2$grSegregation3, r2$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$polarizationIndex~r2$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r3$grSegregation3, r3$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$polarizationIndex~r3$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r4$grSegregation3, r4$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$polarizationIndex~r4$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
dev.off()
#1b)
png(
  filename = "./outputGraphics/1b.png",
  width = 1000,
  height = 1000,
  units = "px"
)
par(mfrow=c(5,4))#,oma=c(0,6,0,0))
# exposure to outgroup
plot(r1$expOutgroup, r1$meanAbsOpinion, xlab="", ylab="mean abs. op.", cex.lab=1.5, main="H=.6, beta op.distr.", pch=".")
abline(lm(r1$meanAbsOpinion~r1$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r2$expOutgroup, r2$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="H=.6, uniform op.distr.", pch=".")
abline(lm(r2$meanAbsOpinion~r2$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r3$expOutgroup, r3$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="H=.9, beta op.distr.", pch=".")
abline(lm(r3$meanAbsOpinion~r3$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r4$expOutgroup, r4$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="H=.9, uniform op.distr.", pch=".")
abline(lm(r4$meanAbsOpinion~r4$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
# proportion non-western
plot(r1$p_nwa, r1$meanAbsOpinion, xlab="", ylab="mean abs. op.", cex.lab=1.5, main="", pch=".")
abline(lm(r1$meanAbsOpinion~r1$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r2$p_nwa, r2$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$meanAbsOpinion~r2$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r3$p_nwa, r3$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$meanAbsOpinion~r3$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r4$p_nwa, r4$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$meanAbsOpinion~r4$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 1)
plot(r1$grSegregation1, r1$meanAbsOpinion, xlab="", ylab="mean abs. op.", cex.lab=1.5, main="", pch=".")
abline(lm(r1$meanAbsOpinion~r1$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r2$grSegregation1, r2$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$meanAbsOpinion~r2$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r3$grSegregation1, r3$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$meanAbsOpinion~r3$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r4$grSegregation1, r4$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$meanAbsOpinion~r4$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 2)
plot(r1$grSegregation2, r1$meanAbsOpinion, xlab="", ylab="mean abs. op.", cex.lab=1.5, main="", pch=".")
abline(lm(r1$meanAbsOpinion~r1$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r2$grSegregation2, r2$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$meanAbsOpinion~r2$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r3$grSegregation2, r3$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$meanAbsOpinion~r3$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r4$grSegregation2, r4$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$meanAbsOpinion~r4$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 3)
plot(r1$grSegregation3, r1$meanAbsOpinion, xlab="", ylab="mean abs. op.", cex.lab=1.5, main="", pch=".")
abline(lm(r1$meanAbsOpinion~r1$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r2$grSegregation3, r2$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$meanAbsOpinion~r2$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r3$grSegregation3, r3$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$meanAbsOpinion~r3$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r4$grSegregation3, r4$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$meanAbsOpinion~r4$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
dev.off()

#1c)
png(
  filename = "./outputGraphics/1c.png",
  width = 1000,
  height = 1000,
  units = "px"
)
par(mfrow=c(5,4))#,oma=c(0,6,0,0))
# exposure to outgroup
plot(r1$expOutgroup, r1$meanOpinionGlobal, xlab="", ylab="mean op.", cex.lab=1.5, main="H=.6, beta op.distr.", pch=".")
abline(lm(r1$meanOpinionGlobal~r1$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r2$expOutgroup, r2$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="H=.6, uniform op.distr.", pch=".")
abline(lm(r2$meanOpinionGlobal~r2$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r3$expOutgroup, r3$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="H=.9, beta op.distr.", pch=".")
abline(lm(r3$meanOpinionGlobal~r3$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r4$expOutgroup, r4$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="H=.9, uniform op.distr.", pch=".")
abline(lm(r4$meanOpinionGlobal~r4$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
# proportion non-western
plot(r1$p_nwa, r1$meanOpinionGlobal, xlab="", ylab="mean op.", cex.lab=1.5, main="", pch=".")
abline(lm(r1$meanOpinionGlobal~r1$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r2$p_nwa, r2$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$meanOpinionGlobal~r2$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r3$p_nwa, r3$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$meanOpinionGlobal~r3$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r4$p_nwa, r4$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$meanOpinionGlobal~r4$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 1)
plot(r1$grSegregation1, r1$meanOpinionGlobal, xlab="", ylab="mean op.", cex.lab=1.5, main="", pch=".")
abline(lm(r1$meanOpinionGlobal~r1$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r2$grSegregation1, r2$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$meanOpinionGlobal~r2$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r3$grSegregation1, r3$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$meanOpinionGlobal~r3$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r4$grSegregation1, r4$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$meanOpinionGlobal~r4$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 2)
plot(r1$grSegregation2, r1$meanOpinionGlobal, xlab="", ylab="mean op.", cex.lab=1.5, main="", pch=".")
abline(lm(r1$meanOpinionGlobal~r1$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r2$grSegregation2, r2$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$meanOpinionGlobal~r2$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r3$grSegregation2, r3$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$meanOpinionGlobal~r3$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r4$grSegregation2, r4$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$meanOpinionGlobal~r4$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 3)
plot(r1$grSegregation3, r1$meanOpinionGlobal, xlab="", ylab="mean op.", cex.lab=1.5, main="", pch=".")
abline(lm(r1$meanOpinionGlobal~r1$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r2$grSegregation3, r2$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$meanOpinionGlobal~r2$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r3$grSegregation3, r3$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$meanOpinionGlobal~r3$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r4$grSegregation3, r4$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$meanOpinionGlobal~r4$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
dev.off()

#2)
png(
  filename = "./outputGraphics/2.png",
  width = 1000,
  height = 1000,
  units = "px"
)
par(mfrow=c(5,4))#,oma=c(0,6,0,0))
# exposure to outgroup
plot(r1$expOutgroup, r1$opClustering2 , xlab="", ylab="op. clustering", cex.lab=1.5, main="H=.6, beta op.distr.", pch=".")
abline(lm(r1$opClustering2~r1$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r2$expOutgroup, r2$opClustering2, xlab="", ylab="", cex.lab=1.5, main="H=.6, uniform op.distr.", pch=".")
abline(lm(r2$opClustering2~r2$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r3$expOutgroup, r3$opClustering2, xlab="", ylab="", cex.lab=1.5, main="H=.9, beta op.distr.", pch=".")
abline(lm(r3$opClustering2~r3$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r4$expOutgroup, r4$opClustering2, xlab="", ylab="", cex.lab=1.5, main="H=.9, uniform op.distr.", pch=".")
abline(lm(r4$opClustering2~r4$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
# NWAexpout
plot(r1$NWAexpOutgroup, r1$opClustering2, xlab="", ylab="op. clustering", cex.lab=1.5, main="", pch=".")
abline(lm(r1$opClustering2~r1$NWAexpOutgroup), col="blue")
mtext("exp.out. of minorities", side=1, line=2, cex=1)
plot(r2$NWAexpOutgroup, r2$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$opClustering2~r2$NWAexpOutgroup), col="blue")
mtext("exp.out. of minorities", side=1, line=2, cex=1)
plot(r3$NWAexpOutgroup, r3$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$opClustering2~r3$NWAexpOutgroup), col="blue")
mtext("exp.out. of minorities", side=1, line=2, cex=1)
plot(r4$NWAexpOutgroup, r4$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$opClustering2~r4$NWAexpOutgroup), col="blue")
mtext("exp.out. of minorities", side=1, line=2, cex=1)
# NATexpout
plot(r1$NATexpOutgroup, r1$opClustering2, xlab="", ylab="op. clustering", cex.lab=1.5, main="", pch=".")
abline(lm(r1$opClustering2~r1$NATexpOutgroup), col="blue")
mtext("exp.out. of majority", side=1, line=2, cex=1)
plot(r2$NATexpOutgroup, r2$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$opClustering2~r2$NATexpOutgroup), col="blue")
mtext("exp.out. of majority", side=1, line=2, cex=1)
plot(r3$NATexpOutgroup, r3$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$opClustering2~r3$NATexpOutgroup), col="blue")
mtext("exp.out. of majority", side=1, line=2, cex=1)
plot(r4$NATexpOutgroup, r4$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$opClustering2~r4$NATexpOutgroup), col="blue")
mtext("exp.out. of majority", side=1, line=2, cex=1)
# proportion non-western
plot(r1$p_nwa, r1$opClustering2, xlab="", ylab="op. clustering", cex.lab=1.5, main="", pch=".")
abline(lm(r1$opClustering2~r1$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r2$p_nwa, r2$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$opClustering2~r2$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r3$p_nwa, r3$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$opClustering2~r3$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r4$p_nwa, r4$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$opClustering2~r4$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 2)
plot(r1$grSegregation2, r1$opClustering2, xlab="", ylab="op. clustering", cex.lab=1.5, main="", pch=".")
abline(lm(r1$opClustering2~r1$grSegregation2), col="blue")
mtext("segregation", side=1, line=2, cex=1)
plot(r2$grSegregation2, r2$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$opClustering2~r2$grSegregation2), col="blue")
mtext("segregation", side=1, line=2, cex=1)
plot(r3$grSegregation2, r3$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$opClustering2~r3$grSegregation2), col="blue")
mtext("segregation", side=1, line=2, cex=1)
plot(r4$grSegregation2, r4$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$opClustering2~r4$grSegregation2), col="blue")
mtext("segregation", side=1, line=2, cex=1)
dev.off()


# 3) Segregation -> opinion alignment ---------
png(
  filename = "./outputGraphics/3.png",
  width = 1000,
  height = 1000,
  units = "px"
)
par(mfrow=c(5,4))#,oma=c(0,6,0,0))
plot(r1$expOutgroup, r1$relativeAlignment , xlab="", ylab="op. clustering", cex.lab=1.5, main="H=.6, beta op.distr.", pch=".")
abline(lm(r1$relativeAlignment~r1$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r2$expOutgroup, r2$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="H=.6, uniform op.distr.", pch=".")
abline(lm(r2$relativeAlignment~r2$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r3$expOutgroup, r3$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="H=.9, beta op.distr.", pch=".")
abline(lm(r3$relativeAlignment~r3$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r4$expOutgroup, r4$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="H=.9, uniform op.distr.", pch=".")
abline(lm(r4$relativeAlignment~r4$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
# NWAexpout
plot(r1$NWAexpOutgroup, r1$relativeAlignment, xlab="", ylab="op. alignment", cex.lab=1.5, main="", pch=".")
abline(lm(r1$relativeAlignment~r1$NWAexpOutgroup), col="blue")
mtext("exp.out. of minorities", side=1, line=2, cex=1)
plot(r2$NWAexpOutgroup, r2$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$relativeAlignment~r2$NWAexpOutgroup), col="blue")
mtext("exp.out. of minorities", side=1, line=2, cex=1)
plot(r3$NWAexpOutgroup, r3$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$relativeAlignment~r3$NWAexpOutgroup), col="blue")
mtext("exp.out. of minorities", side=1, line=2, cex=1)
plot(r4$NWAexpOutgroup, r4$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$relativeAlignment~r4$NWAexpOutgroup), col="blue")
mtext("exp.out. of minorities", side=1, line=2, cex=1)
# NATexpout
plot(r1$NATexpOutgroup, r1$relativeAlignment, xlab="", ylab="op. alignment", cex.lab=1.5, main="", pch=".")
abline(lm(r1$relativeAlignment~r1$NATexpOutgroup), col="blue")
mtext("exp.out. of majority", side=1, line=2, cex=1)
plot(r2$NATexpOutgroup, r2$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$relativeAlignment~r2$NATexpOutgroup), col="blue")
mtext("exp.out. of majority", side=1, line=2, cex=1)
plot(r3$NATexpOutgroup, r3$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$relativeAlignment~r3$NATexpOutgroup), col="blue")
mtext("exp.out. of majority", side=1, line=2, cex=1)
plot(r4$NATexpOutgroup, r4$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$relativeAlignment~r4$NATexpOutgroup), col="blue")
mtext("exp.out. of majority", side=1, line=2, cex=1)
# proportion non-western
plot(r1$p_nwa, r1$relativeAlignment, xlab="", ylab="op. alignment", cex.lab=1.5, main="", pch=".")
abline(lm(r1$relativeAlignment~r1$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r2$p_nwa, r2$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$relativeAlignment~r2$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r3$p_nwa, r3$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$relativeAlignment~r3$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r4$p_nwa, r4$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$relativeAlignment~r4$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 2)
plot(r1$grSegregation2, r1$relativeAlignment, xlab="", ylab="op. alignment", cex.lab=1.5, main="", pch=".")
abline(lm(r1$relativeAlignment~r1$grSegregation2), col="blue")
mtext("segregation", side=1, line=2, cex=1)
plot(r2$grSegregation2, r2$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$relativeAlignment~r2$grSegregation2), col="blue")
mtext("segregation", side=1, line=2, cex=1)
plot(r3$grSegregation2, r3$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$relativeAlignment~r3$grSegregation2), col="blue")
mtext("segregation", side=1, line=2, cex=1)
plot(r4$grSegregation2, r4$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$relativeAlignment~r4$grSegregation2), col="blue")
mtext("segregation", side=1, line=2, cex=1)
dev.off()















r <- subset(GsimResults, GsimResults$initialOpinionDistribution=="beta")
r1 <- subset(r, r$H == 0.6 & r$distanceDecay == 1)
r2 <- subset(r, r$H == 0.6 & r$distanceDecay == 2)
r3 <- subset(r, r$H == 0.6 & r$distanceDecay == 3)
r4 <- subset(r, r$H == 0.9 & r$distanceDecay == 1)
r5 <- subset(r, r$H == 0.9 & r$distanceDecay == 2)
r6 <- subset(r, r$H == 0.9 & r$distanceDecay == 3)

# 0) Descriptives   ---------
png(
  filename = "./outputGraphics/macroDescriptives1.png",
  width = 1000,
  height = 800,
  units = "px"
)
par(mfrow=c(4,6))#,oma=c(0,6,0,0))
# abs difference between groups' mean opinion
hist(abs(r1$meanOpinionG1 - r1$meanOpinionG2), xlab="", ylab="Frequency", cex.lab=1.5, main="H=.6, d-decay function 1")
mtext("Between-group\nop. difference", side=1, line=5, cex=1)
hist(abs(r2$meanOpinionG1 - r2$meanOpinionG2), xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 2", cex.main=1.5)
mtext("Between-group\nop. difference", side=1, line=5, cex=1)
hist(abs(r3$meanOpinionG1 - r3$meanOpinionG2), xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 3", cex.main=1.5)
mtext("Between-group\nop. difference", side=1, line=5, cex=1)
hist(abs(r4$meanOpinionG1 - r4$meanOpinionG2), xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 1", cex.main=1.5)
mtext("Between-group\nop. difference", side=1, line=5, cex=1)
hist(abs(r5$meanOpinionG1 - r5$meanOpinionG2), xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 2", cex.main=1.5)
mtext("Between-group\nop. difference", side=1, line=5, cex=1)
hist(abs(r6$meanOpinionG1 - r6$meanOpinionG2), xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 3", cex.main=1.5)
mtext("Between-group\nop. difference", side=1, line=5, cex=1)
# op polarization
hist(r1$polarizationIndex, xlab="", ylab="Frequency", cex.lab=1.5, main="")
mtext("op. polarization", side=1, line=3, cex=1)
hist(r2$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. polarization", side=1, line=3, cex=1)
hist(r3$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. polarization", side=1, line=3, cex=1)
hist(r4$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. polarization", side=1, line=3, cex=1)
hist(r5$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. polarization", side=1, line=3, cex=1)
hist(r6$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. polarization", side=1, line=3, cex=1)
# opinion clustering
hist(r1$opClustering2, xlab="", ylab="Frequency", cex.lab=1.5, main="")
mtext("op. clustering 2", side=1, line=3, cex=1)
hist(r2$opClustering2, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. clustering 2", side=1, line=3, cex=1)
hist(r3$opClustering2, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. clustering 2", side=1, line=3, cex=1)
hist(r4$opClustering2, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. clustering 2", side=1, line=3, cex=1)
hist(r5$opClustering2, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. clustering 2", side=1, line=3, cex=1)
hist(r6$opClustering2, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. clustering 2", side=1, line=3, cex=1)
# alignment
hist(r1$opAlignment2, xlab="", ylab="Frequency", cex.lab=1.5, main="")
mtext("op. alignment 2", side=1, line=3, cex=1)
hist(r2$opAlignment2, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. alignment 2", side=1, line=3, cex=1)
hist(r3$opAlignment2, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. alignment 2", side=1, line=3, cex=1)
hist(r4$opAlignment2, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. alignment 2", side=1, line=3, cex=1)
hist(r5$opAlignment2, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. alignment 2", side=1, line=3, cex=1)
hist(r6$opAlignment2, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. alignment 2", side=1, line=3, cex=1)
dev.off()


# 1) Segregation -> polarization ---------
png(
  filename = "./outputGraphics/macroResults 1a - segr-pol.png",
  width = 1000,
  height = 1000,
  units = "px"
)
par(mfrow=c(5,6))#,oma=c(0,6,0,0))
# exposure to outgroup
plot(r1$expOutgroup, r1$polarizationIndex, xlab="", ylab="op. polarization", cex.lab=1.5, main="H=.6, d-decay function 1", pch=".")
abline(lm(r1$polarizationIndex~r1$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r2$expOutgroup, r2$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 2", pch=".")
abline(lm(r2$polarizationIndex~r2$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r3$expOutgroup, r3$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 3", pch=".")
abline(lm(r3$polarizationIndex~r3$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r4$expOutgroup, r4$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 1", pch=".")
abline(lm(r4$polarizationIndex~r4$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r5$expOutgroup, r5$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 2", pch=".")
abline(lm(r5$polarizationIndex~r5$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r6$expOutgroup, r6$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 3", pch=".")
abline(lm(r6$polarizationIndex~r6$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
# proportion non-western
plot(r1$p_nwa, r1$polarizationIndex, xlab="", ylab="op. polarization", cex.lab=1.5, main="", pch=".")
abline(lm(r1$polarizationIndex~r1$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r2$p_nwa, r2$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$polarizationIndex~r2$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r3$p_nwa, r3$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$polarizationIndex~r3$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r4$p_nwa, r4$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$polarizationIndex~r4$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r5$p_nwa, r5$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$polarizationIndex~r5$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r6$p_nwa, r6$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$polarizationIndex~r6$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 1)
plot(r1$grSegregation1, r1$polarizationIndex, xlab="", ylab="op. polarization", cex.lab=1.5, main="", pch=".")
abline(lm(r1$polarizationIndex~r1$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r2$grSegregation1, r2$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$polarizationIndex~r2$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r3$grSegregation1, r3$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$polarizationIndex~r3$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r4$grSegregation1, r4$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$polarizationIndex~r4$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r5$grSegregation1, r5$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$polarizationIndex~r5$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r6$grSegregation1, r6$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$polarizationIndex~r6$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 2)
plot(r1$grSegregation2, r1$polarizationIndex, xlab="", ylab="op. polarization", cex.lab=1.5, main="", pch=".")
abline(lm(r1$polarizationIndex~r1$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r2$grSegregation2, r2$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$polarizationIndex~r2$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r3$grSegregation2, r3$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$polarizationIndex~r3$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r4$grSegregation2, r4$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$polarizationIndex~r4$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r5$grSegregation2, r5$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$polarizationIndex~r5$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r6$grSegregation2, r6$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$polarizationIndex~r6$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 3)
plot(r1$grSegregation3, r1$polarizationIndex, xlab="", ylab="op. polarization", cex.lab=1.5, main="", pch=".")
abline(lm(r1$polarizationIndex~r1$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r2$grSegregation3, r2$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$polarizationIndex~r2$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r3$grSegregation3, r3$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$polarizationIndex~r3$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r4$grSegregation3, r4$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$polarizationIndex~r4$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r5$grSegregation3, r5$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$polarizationIndex~r5$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r6$grSegregation3, r6$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$polarizationIndex~r6$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
dev.off()

#
png(
  filename = "./outputGraphics/macroResults 1b - segr-absOp.png",
  width = 1000,
  height = 1000,
  units = "px"
)
par(mfrow=c(5,6))#,oma=c(0,6,0,0))
# exposure to outgroup
plot(r1$expOutgroup, r1$meanAbsOpinion, xlab="", ylab="mean abs. op.", cex.lab=1.5, main="H=.6, d-decay function 1", pch=".")
abline(lm(r1$meanAbsOpinion~r1$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r2$expOutgroup, r2$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 2", pch=".")
abline(lm(r2$meanAbsOpinion~r2$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r3$expOutgroup, r3$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 3", pch=".")
abline(lm(r3$meanAbsOpinion~r3$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r4$expOutgroup, r4$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 1", pch=".")
abline(lm(r4$meanAbsOpinion~r4$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r5$expOutgroup, r5$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 2", pch=".")
abline(lm(r5$meanAbsOpinion~r5$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r6$expOutgroup, r6$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 3", pch=".")
abline(lm(r6$meanAbsOpinion~r6$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
# proportion non-western
plot(r1$p_nwa, r1$meanAbsOpinion, xlab="", ylab="mean abs. op.", cex.lab=1.5, main="", pch=".")
abline(lm(r1$meanAbsOpinion~r1$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r2$p_nwa, r2$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$meanAbsOpinion~r2$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r3$p_nwa, r3$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$meanAbsOpinion~r3$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r4$p_nwa, r4$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$meanAbsOpinion~r4$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r5$p_nwa, r5$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$meanAbsOpinion~r5$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r6$p_nwa, r6$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$meanAbsOpinion~r6$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 1)
plot(r1$grSegregation1, r1$meanAbsOpinion, xlab="", ylab="mean abs. op.", cex.lab=1.5, main="", pch=".")
abline(lm(r1$meanAbsOpinion~r1$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r2$grSegregation1, r2$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$meanAbsOpinion~r2$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r3$grSegregation1, r3$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$meanAbsOpinion~r3$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r4$grSegregation1, r4$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$meanAbsOpinion~r4$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r5$grSegregation1, r5$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$meanAbsOpinion~r5$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r6$grSegregation1, r6$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$meanAbsOpinion~r6$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 2)
plot(r1$grSegregation2, r1$meanAbsOpinion, xlab="", ylab="mean abs. op.", cex.lab=1.5, main="", pch=".")
abline(lm(r1$meanAbsOpinion~r1$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r2$grSegregation2, r2$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$meanAbsOpinion~r2$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r3$grSegregation2, r3$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$meanAbsOpinion~r3$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r4$grSegregation2, r4$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$meanAbsOpinion~r4$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r5$grSegregation2, r5$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$meanAbsOpinion~r5$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r6$grSegregation2, r6$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$meanAbsOpinion~r6$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 3)
plot(r1$grSegregation3, r1$meanAbsOpinion, xlab="", ylab="mean abs. op.", cex.lab=1.5, main="", pch=".")
abline(lm(r1$meanAbsOpinion~r1$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r2$grSegregation3, r2$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$meanAbsOpinion~r2$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r3$grSegregation3, r3$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$meanAbsOpinion~r3$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r4$grSegregation3, r4$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$meanAbsOpinion~r4$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r5$grSegregation3, r5$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$meanAbsOpinion~r5$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r6$grSegregation3, r6$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$meanAbsOpinion~r6$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
dev.off()




#
png(
  filename = "./outputGraphics/macroResults 1c - segr-propExtremists.png",
  width = 1000,
  height = 1000,
  units = "px"
)
par(mfrow=c(5,6))#,oma=c(0,6,0,0))
# exposure to outgroup
plot(r1$expOutgroup, r1$propExtremists, xlab="", ylab="prop. extremists", cex.lab=1.5, main="H=.6, d-decay function 1", pch=".")
abline(lm(r1$propExtremists~r1$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r2$expOutgroup, r2$propExtremists, xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 2", pch=".")
abline(lm(r2$propExtremists~r2$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r3$expOutgroup, r3$propExtremists, xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 3", pch=".")
abline(lm(r3$propExtremists~r3$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r4$expOutgroup, r4$propExtremists, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 1", pch=".")
abline(lm(r4$propExtremists~r4$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r5$expOutgroup, r5$propExtremists, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 2", pch=".")
abline(lm(r5$propExtremists~r5$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r6$expOutgroup, r6$propExtremists, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 3", pch=".")
abline(lm(r6$propExtremists~r6$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
# proportion non-western
plot(r1$p_nwa, r1$propExtremists, xlab="", ylab="prop. extremists", cex.lab=1.5, main="", pch=".")
abline(lm(r1$propExtremists~r1$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r2$p_nwa, r2$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$propExtremists~r2$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r3$p_nwa, r3$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$propExtremists~r3$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r4$p_nwa, r4$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$propExtremists~r4$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r5$p_nwa, r5$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$propExtremists~r5$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r6$p_nwa, r6$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$propExtremists~r6$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 1)
plot(r1$grSegregation1, r1$propExtremists, xlab="", ylab="prop. extremists", cex.lab=1.5, main="", pch=".")
abline(lm(r1$propExtremists~r1$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r2$grSegregation1, r2$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$propExtremists~r2$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r3$grSegregation1, r3$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$propExtremists~r3$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r4$grSegregation1, r4$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$propExtremists~r4$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r5$grSegregation1, r5$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$propExtremists~r5$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r6$grSegregation1, r6$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$propExtremists~r6$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 2)
plot(r1$grSegregation2, r1$propExtremists, xlab="", ylab="prop. extremists", cex.lab=1.5, main="", pch=".")
abline(lm(r1$propExtremists~r1$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r2$grSegregation2, r2$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$propExtremists~r2$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r3$grSegregation2, r3$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$propExtremists~r3$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r4$grSegregation2, r4$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$propExtremists~r4$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r5$grSegregation2, r5$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$propExtremists~r5$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r6$grSegregation2, r6$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$propExtremists~r6$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 3)
plot(r1$grSegregation3, r1$propExtremists, xlab="", ylab="prop. extremists", cex.lab=1.5, main="", pch=".")
abline(lm(r1$propExtremists~r1$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r2$grSegregation3, r2$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$propExtremists~r2$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r3$grSegregation3, r3$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$propExtremists~r3$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r4$grSegregation3, r4$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$propExtremists~r4$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r5$grSegregation3, r5$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$propExtremists~r5$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r6$grSegregation3, r6$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$propExtremists~r6$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
dev.off()

#
png(
  filename = "./outputGraphics/macroResults 1d - segr-opVar.png",
  width = 1000,
  height = 1000,
  units = "px"
)
par(mfrow=c(5,6))#,oma=c(0,6,0,0))
# exposure to outgroup
plot(r1$expOutgroup, r1$varOpinionGlobal, xlab="", ylab="op. variance", cex.lab=1.5, main="H=.6, d-decay function 1", pch=".")
abline(lm(r1$varOpinionGlobal~r1$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r2$expOutgroup, r2$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 2", pch=".")
abline(lm(r2$varOpinionGlobal~r2$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r3$expOutgroup, r3$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 3", pch=".")
abline(lm(r3$varOpinionGlobal~r3$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r4$expOutgroup, r4$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 1", pch=".")
abline(lm(r4$varOpinionGlobal~r4$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r5$expOutgroup, r5$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 2", pch=".")
abline(lm(r5$varOpinionGlobal~r5$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r6$expOutgroup, r6$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 3", pch=".")
abline(lm(r6$varOpinionGlobal~r6$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
# proportion non-western
plot(r1$p_nwa, r1$varOpinionGlobal, xlab="", ylab="op. variance", cex.lab=1.5, main="", pch=".")
abline(lm(r1$varOpinionGlobal~r1$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r2$p_nwa, r2$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$varOpinionGlobal~r2$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r3$p_nwa, r3$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$varOpinionGlobal~r3$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r4$p_nwa, r4$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$varOpinionGlobal~r4$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r5$p_nwa, r5$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$varOpinionGlobal~r5$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r6$p_nwa, r6$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$varOpinionGlobal~r6$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 1)
plot(r1$grSegregation1, r1$varOpinionGlobal, xlab="", ylab="op. variance", cex.lab=1.5, main="", pch=".")
abline(lm(r1$varOpinionGlobal~r1$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r2$grSegregation1, r2$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$varOpinionGlobal~r2$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r3$grSegregation1, r3$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$varOpinionGlobal~r3$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r4$grSegregation1, r4$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$varOpinionGlobal~r4$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r5$grSegregation1, r5$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$varOpinionGlobal~r5$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r6$grSegregation1, r6$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$varOpinionGlobal~r6$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 2)
plot(r1$grSegregation2, r1$varOpinionGlobal, xlab="", ylab="op. variance", cex.lab=1.5, main="", pch=".")
abline(lm(r1$varOpinionGlobal~r1$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r2$grSegregation2, r2$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$varOpinionGlobal~r2$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r3$grSegregation2, r3$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$varOpinionGlobal~r3$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r4$grSegregation2, r4$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$varOpinionGlobal~r4$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r5$grSegregation2, r5$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$varOpinionGlobal~r5$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r6$grSegregation2, r6$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$varOpinionGlobal~r6$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 3)
plot(r1$grSegregation3, r1$varOpinionGlobal, xlab="", ylab="op. variance", cex.lab=1.5, main="", pch=".")
abline(lm(r1$varOpinionGlobal~r1$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r2$grSegregation3, r2$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$varOpinionGlobal~r2$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r3$grSegregation3, r3$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$varOpinionGlobal~r3$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r4$grSegregation3, r4$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$varOpinionGlobal~r4$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r5$grSegregation3, r5$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$varOpinionGlobal~r5$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r6$grSegregation3, r6$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$varOpinionGlobal~r6$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
dev.off()







# 2) Segregation -> opinion clustering ---------
png(
  filename = "./outputGraphics/macroResults 2 - segr-clu.png",
  width = 1000,
  height = 1000,
  units = "px"
)
par(mfrow=c(5,6))#,oma=c(0,6,0,0))
# exposure to outgroup
plot(r1$expOutgroup, r1$opClustering1 , xlab="", ylab="op. clustering 1", cex.lab=1.5, main="H=.6, d-decay function 1", pch=".")
abline(lm(r1$opClustering1~r1$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r2$expOutgroup, r2$opClustering1, xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 2", pch=".")
abline(lm(r2$opClustering1~r2$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r3$expOutgroup, r3$opClustering1, xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 3", pch=".")
abline(lm(r3$opClustering1~r3$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r4$expOutgroup, r4$opClustering1, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 1", pch=".")
abline(lm(r4$opClustering1~r4$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r5$expOutgroup, r5$opClustering1, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 2", pch=".")
abline(lm(r5$opClustering1~r5$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r6$expOutgroup, r6$opClustering1, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 3", pch=".")
abline(lm(r6$opClustering1~r6$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
# proportion non-western
plot(r1$p_nwa, r1$opClustering1, xlab="", ylab="op. clustering 1", cex.lab=1.5, main="", pch=".")
abline(lm(r1$opClustering1~r1$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r2$p_nwa, r2$opClustering1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$opClustering1~r2$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r3$p_nwa, r3$opClustering1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$opClustering1~r3$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r4$p_nwa, r4$opClustering1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$opClustering1~r4$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r5$p_nwa, r5$opClustering1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$opClustering1~r5$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r6$p_nwa, r6$opClustering1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$opClustering1~r6$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 1)
plot(r1$grSegregation1, r1$opClustering1, xlab="", ylab="op. clustering 1", cex.lab=1.5, main="", pch=".")
abline(lm(r1$opClustering1~r1$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r2$grSegregation1, r2$opClustering1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$opClustering1~r2$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r3$grSegregation1, r3$opClustering1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$opClustering1~r3$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r4$grSegregation1, r4$opClustering1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$opClustering1~r4$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r5$grSegregation1, r5$opClustering1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$opClustering1~r5$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r6$grSegregation1, r6$opClustering1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$opClustering1~r6$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 2)
plot(r1$grSegregation2, r1$opClustering2, xlab="", ylab="op. clustering 2", cex.lab=1.5, main="", pch=".")
abline(lm(r1$opClustering2~r1$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r2$grSegregation2, r2$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$opClustering2~r2$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r3$grSegregation2, r3$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$opClustering2~r3$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r4$grSegregation2, r4$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$opClustering2~r4$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r5$grSegregation2, r5$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$opClustering2~r5$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r6$grSegregation2, r6$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$opClustering2~r6$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 3)
plot(r1$grSegregation3, r1$opClustering3, xlab="", ylab="op. clustering 3", cex.lab=1.5, main="", pch=".")
abline(lm(r1$opClustering3~r1$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r2$grSegregation3, r2$opClustering3, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$opClustering3~r2$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r3$grSegregation3, r3$opClustering3, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$opClustering3~r3$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r4$grSegregation3, r4$opClustering3, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$opClustering3~r4$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r5$grSegregation3, r5$opClustering3, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$opClustering3~r5$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r6$grSegregation3, r6$opClustering3, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$opClustering3~r6$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
dev.off()


# 3) Segregation -> opinion alignment ---------
png(
  filename = "./outputGraphics/macroResults 3 - segr-ali.png",
  width = 1000,
  height = 1000,
  units = "px"
)
par(mfrow=c(5,6))#,oma=c(0,6,0,0))
# exposure to outgroup
plot(r1$expOutgroup, r1$opAlignment1, xlab="", ylab="op. alignment 1", cex.lab=1.5, main="H=.6, d-decay function 1", pch=".")
abline(lm(r1$opAlignment1~r1$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r2$expOutgroup, r2$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 2", pch=".")
abline(lm(r2$opAlignment1~r2$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r3$expOutgroup, r3$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 3", pch=".")
abline(lm(r3$opAlignment1~r3$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r4$expOutgroup, r4$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 1", pch=".")
abline(lm(r4$opAlignment1~r4$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r5$expOutgroup, r5$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 2", pch=".")
abline(lm(r5$opAlignment1~r5$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r6$expOutgroup, r6$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 3", pch=".")
abline(lm(r6$opAlignment1~r6$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
# proportion non-western
plot(r1$p_nwa, r1$opAlignment1, xlab="", ylab="op. alignment 1", cex.lab=1.5, main="", pch=".")
abline(lm(r1$opAlignment1~r1$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r2$p_nwa, r2$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$opAlignment1~r2$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r3$p_nwa, r3$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$opAlignment1~r3$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r4$p_nwa, r4$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$opAlignment1~r4$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r5$p_nwa, r5$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$opAlignment1~r5$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r6$p_nwa, r6$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$opAlignment1~r6$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 1)
plot(r1$grSegregation1, r1$opAlignment1, xlab="", ylab="op. alignment 1", cex.lab=1.5, main="", pch=".")
abline(lm(r1$opAlignment1~r1$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r2$grSegregation1, r2$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$opAlignment1~r2$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r3$grSegregation1, r3$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$opAlignment1~r3$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r4$grSegregation1, r4$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$opAlignment1~r4$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r5$grSegregation1, r5$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$opAlignment1~r5$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r6$grSegregation1, r6$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$opAlignment1~r6$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 2)
plot(r1$grSegregation2, r1$opAlignment2, xlab="", ylab="op. alignment 2", cex.lab=1.5, main="", pch=".")
abline(lm(r1$opAlignment2~r1$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r2$grSegregation2, r2$opAlignment2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$opAlignment2~r2$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r3$grSegregation2, r3$opAlignment2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$opAlignment2~r3$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r4$grSegregation2, r4$opAlignment2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$opAlignment2~r4$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r5$grSegregation2, r5$opAlignment2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$opAlignment2~r5$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r6$grSegregation2, r6$opAlignment2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$opAlignment2~r6$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 3)
plot(r1$grSegregation3, r1$opAlignment3, xlab="", ylab="op. alignment 3", cex.lab=1.5, main="", pch=".")
abline(lm(r1$opAlignment3~r1$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r2$grSegregation3, r2$opAlignment3, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$opAlignment3~r2$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r3$grSegregation3, r3$opAlignment3, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$opAlignment3~r3$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r4$grSegregation3, r4$opAlignment3, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$opAlignment3~r4$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r5$grSegregation3, r5$opAlignment3, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$opAlignment3~r5$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r6$grSegregation3, r6$opAlignment3, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$opAlignment3~r6$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
dev.off()








hist(abs(r$meanOpinionG1 - r$meanOpinionG2))
hist(r$polarizationIndex)

hist(r$opClustering1)
hist(r$opClustering2)
hist(r$opClustering3)

hist(r$opAlignment1)
hist(r$opAlignment2)
hist(r$opAlignment3)


# 1) Segregation -> opinion polarization
#
# Global level: 
plot(r$expOutgroup, r$polarizationIndex, xlab="Exposure to outgroup", ylab="Polarization index")
plot(r$NATexpOutgroup, r$polarizationIndex)
plot(r$NWAexpOutgroup, r$polarizationIndex)

plot(r$p_nwa, r$polarizationIndex)
plot(r$grClustering1, r$polarizationIndex)
plot(r$grClustering2, r$polarizationIndex)
plot(r$grClustering3, r$polarizationIndex)
plot(r$grSegregation1, r$polarizationIndex)
plot(r$grSegregation2, r$polarizationIndex)
plot(r$grSegregation3, r$polarizationIndex)

# Local level: 
#plot(rw$exposureOutgroup, mean(abs(rw$opinion)))
hist(
  rw$opinion,
  breaks = 20,
  main = NULL,
  xlim = c(-1,1),
  xlab = "Opinion distribution",
  yaxt = "n",
  ylab = "Frequency"
)

plotHeatMap( # default function. Saves to file.
  data=rw,
  x=rw$exposureOutgroup,
  y=rw$timeFirstPol,
  xlab="exposure to outgroup",
  ylab="timeFirstPol",
  fileName = "test"
)

heatmap( # Custom function. Prints on screen.
  x = rw$exposureOutgroup,
  y = rw$timeFirstPol,
  bins=10,
  type = "count",#"prop",
  xlab = "exposure to outgroup",
  ylab = "# interactions before becoming extremist"
)















plotHeatMap(
  data=rw,
  x=rw$exposureOutgroup,
  y=rw$opinion,
  xlab="exposure to outgroup",
  ylab="opinion",
  fileName = "micro_expOutgroup_opinion2"
)
plotHeatMap(
  data=rw,
  x=rw$exposureOutgroup,
  y=rw$timeFirstPol,
  xlab="exposure to outgroup",
  ylab="# interactions before becoming extremist",
  fileName = "micro_groupClustering_timeFirstPol"
)


  
  
  
ggplot(rw, aes(x=rw$exposureOutgroup, y=abs(rw$opinion))) +
  #geom_point(shape=20, alpha= 0.01)+
  geom_smooth(method='lm')



a <- rw$exposureOutgroup[1:1000]#runif(1000, min=0, max=1)
b <- rw$opinion[1:1000]#rnorm(1000, mean=0, sd=1)
#b <- (a*0.5)+b
t <- as.data.frame(cbind(a,b))

ggplot(t, aes(x=a, y=abs(b))) + geom_point(shape=20,alpha=0.1)+
  geom_smooth(method='lm')





# 2) Segregation -> opinion clustering
#
# Global level: 
plot(r$expOutgroup, r$opClustering1, xlab="exposure to outgroup", ylab="spatial clustering of opinions")
plot(r$NATexpOutgroup, r$opClustering1)
plot(r$NWAexpOutgroup, r$opClustering1)


plot(r$p_nwa, r$opAlignment1)
plot(r$grClustering1, r$opClustering1)
plot(r$grClustering2, r$opClustering2)
plot(r$grClustering3, r$opClustering3)
plot(r$grSegregation1, r$opClustering1)
plot(r$grSegregation2, r$opClustering2)
plot(r$grSegregation3, r$opClustering3)


# Local level: 
plotHeatMap(
  data=rw,
  x=rw$exposureOutgroup,
  y=rw$opClustering,
  bins=40,
  xlab="exposure to outgroup",
  ylab="opinion clustering",
  fileName = "micro_expOutgroup_opClustering"
)


# 3) Segregation -> alignment
#
# Global level: 
plot(r$expOutgroup, r$opAlignment, xlab="exposure to outgroup", ylab="group-opinion alignment")
plot(r$NATexpOutgroup, r$opAlignment)
plot(r$NWAexpOutgroup, r$opAlignment)


plot(r$p_nwa, r$opAlignment)
plot(r$grClustering1, r$opAlignment1)
plot(r$grClustering2, r$opAlignment2)
plot(r$grClustering3, r$opAlignment3)
plot(r$grSegregation1, r$opAlignment1)

# Local level: 
plotHeatMap(
  data=rw,
  x=rw$exposureOutgroup,
  y=rw$opAlignment,
  bins=40,
  xlab="exposure to outgroup",
  ylab="group-opinion alignment",
  fileName = "micro_expOutgroup_opAlignment"
)




test <- subset(r, r$H == .3)
summary(lm(r$opClustering ~ r$WK_CODE))
summary(lm(r$opClustering ~ r$globalI_biv_count + r$p_nwa))

plot(r$WK_CODE, r$opClustering)



hist(r$opClustering)
hist(abs(rw$opAlignment))

hist(r$polarizationIndex)



################################################################################
table(names(r) %in% names(rw))
names(r)

# Adding independent variables to the complete agentlist rw
temp  <- r[,-c()]
r <- base::merge(
  x = rw,
  y = temp,
  by = "seed"
)
rm(temp)

















plot(r$expOutgroup, r$polarizationIndex)
summary(lm(r$polarizationIndex ~ r$expOutgroup))

plot(r$expOutgroup, abs(r$meanOpinionG2))
summary(lm(r$polarizationIndex ~ abs(r$meanOpinionG2)))

plot(r$expOutgroup, abs(r$meanOpinionG1 - r$meanOpinionG2))
summary(lm(r$polarizationIndex ~ abs(r$meanOpinionG1 - r$meanOpinionG2)))



run(
  seed = 12345,
  indexParameters = 99999,
  wijk = 1,
  initialOpinionDistribution = "uniform",
  H = 0.3,
  distanceDecay = 2,
  timeMax = 10,
  printOpinionHistogram = TRUE,
  exportOutput = FALSE
)
hist(test[[2]]$opinion)



# ==============================================================================
# Old implementation (do not run anything that follows)

# Time to save everything into a single file, that we shall name
# "completeDataset". This may well take a long while (and a lot of memory)
parameterSpace <- GparameterSpace
parameterSpace$printOpinionHistogram <- parameterSpace$exportOutput <- NULL
simResults <- GsimResults
simW <- GsimW
save(
  parameterSpace, simResults, simW,
  file = "./simOutput/completeDataset.RDATA"
  #file = "./simOutput/completeBaselineDataset.RDATA"
)

# Loading results, battery by battery. We combine all the results into
# global objects (GsimResults and GsimW)
GparameterSpace <- data.frame()
GsimResults <- data.frame()
GsimW <- list()
for (i in 1:length(files)){
  print(paste("Loading file", i, "of", length(files)))
  file <- files[i]
  load(paste0("./simOutput/peregrine/", files[i]))
  if (i == 1){
    GparameterSpace <- parameterSpace
  } else {
    if (any(GparameterSpace != parameterSpace)){
      stop(print("Can't load results batteries with inconsistent parameterSpace."))
    }
  }
  #simResults <- simResults[49:60,] ###################
  #simW <- simW[49:60]              ###################
  GsimResults <- as.data.frame(rbind(GsimResults,simResults))
  #GsimW <- c(GsimW, simW)# this is just too massive : it's 5GB when compressed..
}

# Time to save everything into a single file, that we shall name
# "completeDataset". This may well take a long while (and a lot of memory)
parameterSpace <- GparameterSpace
parameterSpace$printOpinionHistogram <- parameterSpace$exportOutput <- NULL
simResults <- GsimResults
simW <- GsimW
save(
  parameterSpace, simResults, simW,
  file = "./simOutput/completeDataset.RDATA"
  #file = "./simOutput/completeBaselineDataset.RDATA"
)
rm(files, GparameterSpace, GsimResults, GsimW)


#
# 
# Here we prepare the datasets for the data analyses.

# We import the settings and the output of the battery of simulation runs.
rm (list = ls( )) 
load("./cityData/geodata_Rotterdam.RData")
rm(worldList)
#load("./simOutput/sims_1.RData") #############################################
load("./simOutput/completeBaselineDataset.RDATA")

# We join the dataframes with parameter settings, simulation results, and 
# district descriptives.
parameterSpace$indexParameters <- c(1:nrow(parameterSpace))
r <- base::merge(
  x = simResults,
  y = parameterSpace,
  by = "indexParameters"
)
r$order <- c(1:nrow(r)) # Because joining dataframes can shuffle the ordering of
# the simulation runs, we create an index variable that
# we can later use to sort the dataframe in the correct
# order
citySummary$WK_NR <- c(1:nrow(citySummary))
r <- base::merge(
  x = r,
  y = citySummary,
  by.x = "wijk",
  by.y = "WK_NR"
)
r <- r[order(r$order),] # We sort the dataframe, and rid of the ordering index.
#r <- r[,-1]


# Adding measures of opinion clustering and alignment.
#
# We retrieve the opinion vector from the simulated district, we use it to
# calculate uni- and bivariate, global and local I.
# Then we merge the information to the results dataframes, r and rw:
#
#  > r: each row is a summary of a simulation run - contains global variables.
#  > rw: each row is an agent. Contains local variables. Variable rw$simRun 
#          nests agents by the simulation run they were in.
r$opClustering <- r$opAlignment <- NA
#rw <- as.data.frame(NA)
for (i in 1:nrow(r)){
  w <- simW[[i]]
  w$simRun <- i
  proxmat <- NA
  if (r$distanceDecay[i] == 1) {proxmat <- proximityList1[[r$wijk[i]]]}
  if (r$distanceDecay[i] == 2) {proxmat <- proximityList2[[r$wijk[i]]]}
  if (r$distanceDecay[i] == 3) {proxmat <- proximityList4[[r$wijk[i]]]}
  dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[r$wijk[i]])
  dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
  ops <- c(NA, length = length(dat))
  for (l in 1:length(dat)){
    cell <- subset(w, w$location == dat$OBJECTID[l])
    ops[l] <- mean(cell$opinion) ################  <  missing local measures here
  }
  dat$opinion <- ops
  print(paste0(
    "Run ", i, " of ", nrow(r), ". Calculating Moran's I on various attributes."
  ))
  
  # Global opinion clustering
  r$opClustering[i] <- MoranI(
    x = dat$opinion,
    proxmat = proxmat
  )
  
  # Global opinion alignment
  r$opAlignment[[i]] <- MoranI(
    x = dat$pnwal2014,
    y = dat$opinion,
    proxmat = proxmat
  )
  
  if(FALSE){ #########################
    # Local opinion clustering
    Ic <- MoranI(
      x=dat$opinion,
      proxmat=proxmat,
      type = "local"
    )
    
    # Local opinion alignment
    Ia <- MoranI(
      x=dat$pnwal2014,
      y=dat$opinion,
      proxmat=proxmat,
      type = "local"
    )
    ################################### presi da qui e messi in calib: Is, Iec
    
    for (a in 1:nrow(w)){
      w$opClustering[a] <- Ic[w$index[a]]
      w$opAlignment[a] <- Ia[w$index[a]]
      #w$localI_biv_count[a] <- Is[w$index[a]]
      #w$localI_count[a] <- Iec[w$index[a]]
      #w$p_nwa[a] <- dat$pnwal2014[w$index[a]]
    }
    
    #if(i<6){print(paste(r$opClustering[i], r$opAlignment[i], Ic[w$index[i]],Ia[w$index[i]]))}#########
    rm(Ic, Ia)
    ifelse(
      i != 1,
      rw <<- as.data.frame(rbind(rw, w)),
      rw <<- w
    )
  } ###################
}

# Saving the results datasets.
save(
  r,
  rw,
  file = "./simOutput/analyses.RDATA"
)


# Built-in function to plot heatmaps
#
plotHeatMap <- function(
  data,
  x,
  y,
  bins=20,
  xlab,
  ylab,
  fileName,
  legendTitle = "density"
){
  png(paste0("outputGraphics/", fileName,  ".png"), 
      width = 800, height = 600,
      units = "px", pointsize = 10, 
      res = NA)
  print(ggplot(
    data,
    stat = "density",
    aes(
      x=x,
      y=y
    )
  ) + #xlim(c(0.01, 1)) +
    scale_x_continuous(limits = c(NA,NA), expand = c(0, 0)) +
    scale_y_continuous(limits = c(NA,NA), expand = c(0, 0)) +
    geom_bin2d(bins = bins) +#, color ="white")+
    labs(
      fill = legendTitle,
      x = xlab,
      y = ylab
    ) +
    scale_fill_gradient(low =  "white", high = "black") +
    theme(
      axis.line=element_blank(),
      axis.ticks.length = unit(.25, "cm"),
      #axis.ticks = element_line(),
      axis.text = element_text(size=20), ####
      axis.title.y = element_text(size=25),
      axis.title.x = element_text(size = 25),
      title = element_text(size=25),
      legend.text = element_text(size=15),
      panel.background=element_blank(), plot.background=element_blank(),
      panel.border = element_blank()
    ))
  dev.off()
}



load("./testProximityMatrices.RDATA")

distance <- function1 <- function2 <- function3 <- c()
index <- 1
for (x in c(100, 200, 300, 500, 1000, 2000, 5000)){
  distance[index] <- x
  function1[index] <- sum(distances1 <= x) /length(distances1)
  function2[index] <- sum(distances2 <= x)/length(distances2)
  function3[index] <- sum(distances3 <= x)/length(distances3)
  
  index <- index + 1
}
test <- as.data.frame(cbind(distance,function1,function2,function3))
test
rm(test)


par(mfrow=c(1,3))
hist(
  distances1,
  col="red",
  xlim = c(0,2000),
  breaks=c(0:1000)*50,
  main="Distance decay function 1",
  yaxt = "n",
  ylab=NULL,
  xlab="distance between i and j (in meters)"
)
mtext(paste0("Max distance: ", round(max(distances1)), "m"), side=3)
hist(
  distances2,
  col="blue",
  xlim = c(0,2000),
  breaks=c(0:1000)*50,
  main="Distance decay function 2",
  yaxt = "n",
  ylab=NULL,
  xlab="distance between i and j (in meters)"
)
mtext(paste0("Max distance: ", round(max(distances2)), "m"), side=3)
hist(
  distances3,
  col="green",
  xlim = c(0,2000),
  breaks=c(0:1000)*50,
  main="Distance decay function 3",
  yaxt = "n",
  ylab=NULL,
  xlab="distance between i and j (in meters)"
)
mtext(paste0("Max distance: ", round(max(distances3)), "m"), side=3)



# measuring baseline alignment at t=0
d <- ali1 <- ali2 <- ali3 <- aliCon <- aliExtr <- c()
for (r in 0:9){
  print(paste0("replication #", r + 1 ," of 10"))
  for(i in 1:12){
    d[12*r + i] <- i
    w <- worldList[[i]]
    w$opinion <- rbeta(nrow(w), 3, 3, ncp = 0)
    w$opinion <- w$opinion * 2 - 1
    dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[i])
    dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
    ops <- c(NA, length = length(dat))
    for (l in 1:length(dat)){ # for every cell
      cell <- subset(w, w$location == dat$OBJECTID[l])
      ops[l] <- mean(cell$opinion)
    }
    if(FALSE){ ######################
    dat$opinion <- ops
    ali1[12*r + i] <- mean(abs(MoranI( # global alignment
      x = dat$pnwal2014,
      y = dat$opinion,
      proxmat = proximityList1[[i]],
      type = "local"
    )))
    ali2[12*r + i] <- mean(abs(MoranI( # global alignment
      x = dat$pnwal2014,
      y = dat$opinion,
      proxmat = proximityList2[[i]],
      type = "local"
    )))
    ali3[12*r + i] <- mean(abs(MoranI( # global alignment
      x = dat$pnwal2014,
      y = dat$opinion,
      proxmat = proximityList3[[i]],
      type = "local"
    )))
    }##########################
    
    dat$opinion <- 0
    
    aliCon[12*r + i] <- mean(abs(MoranI( # global alignment
      x = dat$pnwal2014,
      y = dat$opinion,
      proxmat = proximityList1[[i]],
      type = "local"
    )))
    dat$opinion <- 1
    aliExtr[12*r + i] <- mean(abs(MoranI( # global alignment
      x = dat$pnwal2014,
      y = dat$opinion,
      proxmat = proximityList1[[i]],
      type = "local"
    )))
  }
}





# Counting interactin events between random agents.
rm (list = ls( )) 
#dist0 <- dist100 <- dist500 <- dist1000 <- c()
# Loading resources
source("script NI&PA.R")
source("util.R")
require(spdep)
require(geosphere)
load("./cityData/geodata_Rotterdam.RData")
run (
  timeMax = 0,
  wijk = 1,
  exportOutput = FALSE
)
dat <- cbs100_rot[cbs100_rot$WK_CODE==districtsList[1],
                  c("nauto2014", "nnwal2014")]
dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
ncells <- length(dat)
distmat <- matrix(NA, nrow=ncells, ncol=ncells)
for (i in 1:ncells) {
  distmat[i,] <- distVincentyEllipsoid(coordinates(dat)[i,], coordinates(dat)) 
}

proxmat <- proximityList2[[1]]
egos <- alters <- distances <- c()
nCells <- length(table(agents$index))
popDensity <- c()
for (i in 1:nCells) {
  popDensity[i] <- nrow(subset(agents, agents$index == i))
}
probmat <- matrix(NA, nrow=nCells, ncol=nCells)
for (i in 1:nCells){
  probmat[,i] <- proxmat[,i] * popDensity[i]
}

for (ego in sample(nrow(agents), size=1)) {
  j <- sample(nrow(agents), size=1000)
  print(paste("counting encounters between", ego, "and 1000 potential alters"))
  for (iterations in 1:500000){
    egoCell <- agents$index[ego]
    targetCell <- sample(c(1:nCells), 1, prob=probmat[egoCell,])
    distance <- distmat[egoCell, targetCell]
    #distances1 <<- append(distances1, distmat[egoCell, targetCell]
    repeat{
      alter <- sample(c(1:popDensity[targetCell]),1)
      alter <- as.numeric(rownames(agents[which(agents$index==targetCell),][alter,]))
      if (ego != alter) {break}
    }
    if (alter %in% j){
      #print("found one")
      egos <- append(egos, ego)
      alters <- append(alters, alter)
      distances <- append(distances, distance)
    }
  }
  #print(paste("i=",ego,"from",egoCell,"j=",alter,"from",targetCell,"distance", distance))
}
cont3 <- as.data.frame(cbind(egos,alters,distances))

dist <- tally <- c()
for (i in unique(cont3$alters)){
  x <- subset(cont3, cont3$alters==i)
  dist <- append(dist, x$distances[1])
  tally <- append(tally, nrow(x))
  print(paste(i, "dist:", x$distances[1], "tally", nrow(x)))
}
test <- as.data.frame(cbind(unique(cont3$alters), dist, tally))
#print(test)
print(head(test[order(test$dist),]))
plot(dist, tally, xlim=c(0,1500))



#

rm (list = ls( )) 

source("script NI&PA.R")
test<-run (
  timeMax = 10,
  wijk = 9,
  exportOutput = TRUE
)






#rm (list = ls( )) 
#load("./simOutput/completeDataset.RDATA")



#==============================================================================
rm (list = ls( )) 


source("script NI&PA.R")
source("util.R")
run(timeMax = 0, wijk = 1)
i <- 1

dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[i])
dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
dat$observed1 <- dat$nnwal2014 - mean(dat$nnwal2014)
dat$grClustering1 <- MoranI(
  x=dat$nnwal2014,
  proxmat=proximityList1[[i]],
  type = "local"
)
hist(dat$grClustering1)
plot(dat$observed1, dat$grClustering1)

for (l in 1:nrow(dat)){
  ifelse(
    dat$j[l] < mean(dat$j),
    dat$test[l] <- 0,
    dat$test[l] <- 1
  )
}
dat$observed2 <- dat$test - mean(dat$test)
dat$testClu <- MoranI(
  x=dat$test,
  proxmat=proximityList1[[i]],
  type = "local"
)
hist(dat$testClu)
plot(dat$observed2,dat$testClu)


#============================================================================
# making sure MoranI() works as intended
run(timeMax = 0, wijk = 1)
i <- 1
dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[i])
dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
dat$gr <- dat$nnwal2014 - mean(dat$nnwal2014)
ops <- c(NA, length = length(dat))
for (l in 1:length(dat)){ # for every cell
  cell <- subset(agents, agents$location == dat$OBJECTID[l])
  ops[l] <- mean(cell$opinion)
}
dat$opinion <- ops

dat$segr <- MoranI(
  x=dat$gr,
  y=dat$opinion,
  proxmat=proximityList1[[i]],
  type = "local"
)
plot(dat$opinion, dat$segr)



for(a in 1:nrow(agents)){
  ifelse(
    agents$group[a] == 1,
    agents$opinion[a] <- 0.7,
    agents$opinion[a] <- -0.7
  )
}
ops <- c(NA, length = length(dat))
for (l in 1:length(dat)){ # for every cell
  cell <- subset(agents, agents$location == dat$OBJECTID[l])
  ops[l] <- mean(cell$opinion)
}
dat$opinion <- ops

dat$segr <- MoranI(
  x=dat$gr,
  y=dat$opinion,
  proxmat=proximityList1[[i]],
  type = "local"
)
plot(dat$opinion, dat$segr)




rr <- subset(
  r,
  r$initialOpinionDistribution == "groupBias" & 
    r$H == 0.6 &
    r$distanceDecay == 2
)
rri <- subset(
  ri,
  ri$initialOpinionDistribution == "groupBias" & 
    ri$H == 0.6 &
    ri$distanceDecay == 2
)
png(
  filename = "./outputGraphics/inspectionAlignmentMeasures_macro_uniform.png",
  width = 600,
  height = 1000,
  units = "px"
)
par(mfrow=c(6,4))
for (wijk in 1:12){
  a <- rr[rr$wijk==wijk,]
  for (i in 1:100){
    a$meanRawScores[i] <- mean(rri[rri$wijk==a$wijk[i]&rri$seed==a$seed[i],]$opAlignment2)
    a$meanAbsScores[i] <- mean(abs(rri[rri$wijk==a$wijk[i]&rri$seed==a$seed[i],]$opAlignment2))
  }
  print(plot(a$intuitiveAlignment, a$meanRawScores, main=districtsNames[wijk],
       xlab="global alignment", ylab="average *raw* alignment scores"))
  print(plot(a$intuitiveAlignment, a$meanAbsScores, #main=districtsNames[wijk],
       xlab="global alignment", ylab="average *absolute* alignment scores"))
}
dev.off()



# Add district-level global-alignment measure to agent-level data
temp <- base::merge( 
  x = rri,
  y = rr[,c(6,76)],
  by = "seed"
)

png(
  filename = "./outputGraphics/inspectionAlignmentMeasures_new.png",
  width = 600,
  height = 500,
  units = "px"
)
par(mfrow=c(3,4))
for (wijk in 1:12){
  a <- temp[temp$wijk==wijk,]
  x <- y <- c()
  for (i in 1:100){
    b <- a[a$seed == unique(a$seed)[i],]
    x[i] <- mean(b$intuitiveAlignment)
    y[i] <- abs(mean(b$opAlignment2))
  }
  print(plot(
    x,
    y,
    main=districtsNames[wijk],
    xlab="global alignment", ylab="*raw* alignment scores"
  ))
}
dev.off()



png(
  filename = "./outputGraphics/inspectionAlignmentMeasures_micro.png",
  width = 600,
  height = 1000,
  units = "px"
)
par(mfrow=c(6,4))
for (wijk in 1:12){
  a <- temp[temp$wijk==wijk,]
  # extracting only one run:
  a <- a[a$seed == a$seed[1],]
  #for (i in 1:100){
  #  a$meanRawScores[i] <- mean(rri[rri$wijk==a$wijk[i]&rri$seed==a$seed[i],]$opAlignment2)
  #  a$meanAbsScores[i] <- mean(abs(rri[rri$wijk==a$wijk[i]&rri$seed==a$seed[i],]$opAlignment2))
  #}
  print(plot(a$intuitiveAlignment, a$opAlignment2, main=districtsNames[wijk],
             xlab="global alignment", ylab="*raw* alignment scores"))
  print(plot(a$intuitiveAlignment, abs(a$opAlignment2), #main=districtsNames[wijk],
             xlab="global alignment", ylab="*absolute* alignment scores"))
}
dev.off()





wijk <- 1
agents <- worldList[[wijk]]
G1 <- which(agents$group == 1, arr.ind = TRUE)
G2 <- which(agents$group != 1, arr.ind = TRUE)
o1 <- rbeta(length(G1), 3, 3.5, ncp = 0)
o2 <- rbeta(length(G2), 3.5, 3, ncp = 0)
for (i in 1:nrow(agents)){
  if (agents$group[i] == 1) {
    agents$opinion[i] <- o1[1]
    o1 <- o1[-1]
  } else {
    agents$opinion[i] <- o2[1]
    o2 <- o2[-1]
  }
}
agents$opinion <- rbeta(nrow(agents), 3,3,ncp=0)
agents$opinion <- agents$opinion * 2 - 1
#agents[agents$group==-1,]$opinion <- 1
#agents[agents$group==1,]$opinion <- -1
dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[wijk])
dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
ops <- c(NA, length = length(dat))
for (l in 1:length(dat)){ # for every cell
  cell <- subset(agents, agents$location == dat$OBJECTID[l])
  ops[l] <- mean(cell$opinion)
}
dat$opinion <- ops
opAlignment2 <- MoranI( # opinion-group alignment
  x = dat$pnwal2014,
  y = dat$opinion,
  proxmat = proximityList2[[wijk]],
  type = "local"
)
agents <- base::merge(
  x=agents,
  y=as.data.frame(cbind(dat$OBJECTID,opAlignment2)),
  by.x="location",
  by.y = "V1"
)
lab = "no group bias"
par(mfrow=c(1,2))
hist(agents$opAlignment2, main=lab, xlab="raw local alignment")
hist(abs(agents$opAlignment2), main=lab, xlab="abs local alignment")


wijk <- 9
temp <- base::merge( 
  x = rri,
  y = rr[,c(6,76)],
  by = "seed"
)
par(mfrow=c(1,2))
a <- temp[temp$wijk==wijk,]
#a <- a[a$seed == a$seed[1],]

print(plot(a$intuitiveAlignment, a$opAlignment2, main=districtsNames[wijk],
           xlab="global alignment", ylab="*raw* alignment scores",pch=1,cex=0.1))
print(plot(a$intuitiveAlignment, abs(a$opAlignment2), #main=districtsNames[wijk],
           xlab="global alignment", ylab="*absolute* alignment scores",pch=1,cex=0.1))

summary(lm(a$intuitiveAlignment~abs(a$opAlignment2)))













x <- c()
for (i in 1:12){
  print(paste("District",i,"of 12"))
  a <- worldList[[i]]
  a$group <- -1
  a[sample(1:nrow(a), size=citySummary$n_nwa[i]),]$group <- 1
  
  dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[i])
  dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
  
  maxE <- calcExposure (
    agents = a,
    proxmat = proximityList2[[i]],
    dat = dat
  )[,3]
  
  x <- cbind(x, maxE)
  citySummary$maxExposure[i] <- mean(maxE)
}
print(mean(x))

length(x) == sum(citySummary$n_pop)






condition <- "beta"

rr <- subset(
  r,
  r$initialOpinionDistribution == condition & 
    r$H == 0.6 &
    r$distanceDecay == 2
)
rri <- subset(
  ri,
  ri$initialOpinionDistribution == condition  & 
    ri$H == 0.6 &
    ri$distanceDecay == 2
)

temp <- base::merge( 
  x = rri,
  y = rr[,c(6, 9, 33, 76)],
  by = "seed"
)

png(
  filename = paste0("./outputGraphics/zzz_",condition,".png"),
  width = 580,
  height = 500,
  units = "px"
)
par(mfrow=c(3,4))
for (wijk in 1:12){
  a <- temp[temp$wijk==wijk,]
  x <- y <- c()
  for (i in 1:100){
    b <- a[a$seed == unique(a$seed)[i],]
    x[i] <- mean(b$intuitiveAlignment)
    #y[i] <- abs(mean(b$opAlignment2))
    #y[i] <- mean(b$opNetAli2)
    y[i] <- mean(b$opAlignment2)
  }
  print(plot(
      x,
      y,
      #main=districtsNames[wijk],
      xlab="global alignment", ylab="bivariate Moran I)", #ylab="abs(mean(raw local alignment))",
      xlim=c(0,2),
      ylim=c(0,1)
    ),
    title(districtsNames[wijk], line=0.5)
  )
  title(
    paste0("initial opinion distribution:", condition),
    outer=TRUE, adj=0, cex.main=2, line=-2
  )
}
dev.off()



######################################################################

# Strong local alignment without global alignment
rr <- subset(
  r,
  r$initialOpinionDistribution == "groupBias" &
    r$wijk == 1 &
    r$H == 0.6 &
    r$distanceDecay == 2
)
rri <- subset(
  ri,
  ri$initialOpinionDistribution == "groupBias"  & 
    ri$wijk == 1 &
    ri$H == 0.6 &
    ri$distanceDecay == 2
)

hist(rri$opAlignment2)
hist(rr$intuitiveAlignment)
hist(rr$opNetAli2)

# find run with highest abs(mean(local alingment))
s <- NA
x <- 0
for (i in 1:100){
  temp <- rri[rri$seed==unique(rr$seed)[i],]
  xx <- abs(mean(temp$opAlignment2))
  if (xx > x) {
    x <- xx
    s <- unique(rr$seed)[i]
  }
}
#abs(mean(rri[rri$seed==s,]$opAlignment2))



mean(rri[rri$seed==s & rri$group!=1,]$opinion)
rr[rr$seed==s,]$intuitiveAlignment


library(compiler)
enableJIT(1)
test <- run (
  timeMax = 20,
  seed  =  s,
  initialOpinionDistribution = "groupBias",
  wijk = 1,
  H = 0.6,
  distanceDecay = 2,
  
  exportOutput = FALSE,
  exportTimeSeries = TRUE
)
enableJIT(0)


meanOpG1 <- meanOpG2 <- meanOpGlob <- sdOpG1 <- sdOpG2 <- sdOpGlob <- ali <- pol <- c()
for (i in 1:length(test[[1]])){
  temp <- test[[1]][[i]]
  tempi <- test[[2]][[i]]
  meanOpG1[i] <- temp$meanOpinionG1
  meanOpG2[i] <- temp$meanOpinionG2
  meanOpGlob[i] <- temp$meanOpinionGlobal
  sdOpG1[i] <- sd(tempi[tempi$group==1,]$opinion)
  sdOpG2[i] <- sd(tempi[tempi$group!=1,]$opinion)
  sdOpGlob[i] <- sd(tempi$opinion)
  ali[i] <- abs(mean(tempi$opAlignment2))
  pol [i] <- temp$polarizationIndex
}


temp <- as.matrix(cbind(meanOpG1,meanOpG2,meanOpGlob,sdOpG1,sdOpG2,sdOpGlob,ali,pol))
par(mar=c(5.1, 2.1, 4.1, 8.8), xpd=TRUE)
matplot(
  temp, type = c("b"),pch=1,col = 1:nrow(temp),
  main="Stadscentrum - groupBias",
  xlab="simulation steps", ylab=""
)
legend(
  "topright", inset=c(-0.41,0),#"bottomleft",
  pch=1,
  col=1:nrow(temp),
  legend = c(
    "mean group -1", "mean group 1", "mean global",
    "sd group -1", "sd group 1", "sd global", "abs(mean(raw_ali))", "pol"
  )
)



# Expectation 2b) 
# Districts with higher levels of mean outgroup exposure exhibit higher average
# local alignment.
print(districtViolins(
  depVar="opNetAli2",
  indepVar="expOutgr2",
  depVarLabel="average local alignment (s=100)",
  indepVarLabel="average outgroup exposure (s=100)"))


opAlignment2

print(districtViolins(
  depVar="opAlignment2",
  indepVar="expOutgr2",
  depVarLabel="bivariate I (global score, s=100)",
  indepVarLabel="average outgroup exposure (s=100)"))










# ==============================================================================
rm(list=ls())

# Here is our Moran I function, as of our last meeting.
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


# Let's simulate some individual agents, and save them in a dataframe "di"
set.seed(1234)
n <- 50
di <- data.frame(
  
  # I give them a random location out of the ten available. Here I assume the
  # world is one-dimensional: that is, locations/vierkanten are points on a
  # line.
  position = sample(1:10, n, replace=TRUE),
  
  # And a random opinion.
  opinion = runif (n, -1, 1)
)

# Then I give agents a group (either +1 or -1). I correlate the group with the
# location: the closer to location 1, the higher the probability of belonging to
# group -1. In other words, we are modeling some degree of group segregation.
for(i in 1:n){
  di$group[i] <- sample(
    c(-1,1),
    size = 1,
    prob = c(di$position[i] / 10, 1 - (di$position[i] / 10))
  )
}

# Now we calculate the aggregate statistics; for my convenience I do this in a
# function. This way I can quickly update the aggregate stats after I modify
# some agents.
calcAlignment <- function(di){
  
  # The aggregate, "vierkanten" stats are stored in a dataframe "d".
  # d$dens will contain the population density of each of the 10 locations:
  d <- data.frame(dens=as.vector(table(di$position) / n))
  
  # Then, for each location, I calculate the average opinion and the proportion
  # of agents who belong to group -1 (the nwal).
  for(l in 1:10){
    d$opinion[l] = mean(di[di$position == l,]$opinion)
    d$nwal[l] = sum(di$group[di$position==l]==-1) / nrow(di[di$position==l,])
  }
  
  # Now I calculate the distance & proximity matrix.
  distmat <- matrix(NA, nrow=10, ncol=10)
  for(r in 1:10) {for (c in 1:10) {distmat[r,c] <- abs(r - c) }}
  proxmat <- exp(-distmat)
  diag(proxmat) <- 0.5
  w <- proxmat <- proxmat/rowSums(proxmat)
  
  # There's all we need to calculate alignment:
  print(paste(
    "alignment =",
    abs(moranI(
      x = d$opinion,
      y = d$nwal,
      dens = d$dens,
      proxmat = w,
      N = n
    )$globalI)
  ))
}

# So let's see how much alignment there is. Opinions were given at random, so
# we should find little to no alignment:
calcAlignment(di)
# We obtained 0.149... - I guess it's a low figure.

# Now let's try to find the max alignment. So far we thought that max alignment
# is when all agents of one group have one opinion, and the agents of the other
# grup have a different opinion.
# So let's try just that:
for (i in 1:n){
  ifelse(
    di$group[i] == 1,
    di$opinion[i] <- -1,
    di$opinion[i] <- 1
  )
}

# Now I can easily re-calculate the aggregated statistics for the locations, and
# print alignment again.
calcAlignment(di)
# We got alignment = 0.607
# It's higher, as it should be. So far, so good.

# Now, if 0.607 is the max alignment, there are no changes to agents' opinion
# that can result in higher alignemnt scores.
#
# However, after a bit of fiddling, I found a case where we get an even higher
# alignment score. All we have to do is to assign a different opinion to agent
# 35:
di$opinion[35] <- 0.5
calcAlignment(di)
# Now alignment = 0.615, ever so slightly higher than our "max" alignment.





# _______________________________________________
rm(list = ls())
require(spdep)

# fake data
#set.seed(86235)
agent_id <- 1:1000
opinions <- runif(1000, -1, 1)
groups <- opinions + rnorm(n = 1000, mean = 0, sd = 1)
cor(opinions, groups)

# y define coordinates for each agent
x <- matrix(NA, nrow = 1000, ncol = 2)

#for (r in 1:nrow(x)){
#  x[r,1]<- sample (1:30, replace=TRUE, size =1)
#}
x[,2] <- 1:1000
x[,1] <- 1000:1
distance <- as.matrix(dist(x, diag = TRUE, upper = TRUE))

proximity <- exp(-3 * (distance)/max(distance))
diag(proximity) <- 0

proximity2 <- proximity
# equaldistances for everone
proximity2[proximity2 != 0] <- 1

coords <- coordinates(x)

proximity <- proximity/rowSums(proximity)

weightslist <- list()
for (i in 1:nrow(proximity)) {
  weightslist[[i]] <- proximity[i, ][-i]
}

nbs <- knn2nb(knearneigh(coords, k = 4))

nblist <- nb2listw(neighbours = nbs, weightslist)


# local
moranI <- function(x, y = NULL, proxmat, dens = NULL, N = length(x)) {
  # Adapted from Anselin (1995, eq. 7, 10, 11)
  # https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1538-4632.1995.tb00338.x
  
  # dens: the proportion of individuals in each cell over the district population if individual level
  # data dens is.null and N is simply length of input if we have aggregate data then N should be total
  # population size (or actually just a large number)
  
  # the N is actually arbitrary, thus simply make sufficiently large
  if (is.null(y)) {
    y <- x
  }
  if (is.null(dens)) {
    dens <- rep(1/N, times = N)
  }
  
  # correct scaling of opinions for densities
  v1dens_ind <- rep(x, times = (dens * N))
  v1dens <- (x - mean(v1dens_ind))/sd(v1dens_ind)
  v2dens_ind <- rep(y, times = (dens * N))
  v2dens <- (y - mean(v2dens_ind))/sd(v2dens_ind)
  
  # (density) weighted proximity matrix
  w <- proxmat
  wdens <- t(dens * t(w))
  wdens <- wdens/rowSums(wdens)
  
  # density and proximity weighted locals
  localI <- (v1dens * wdens %*% v2dens)  #formula 7
  
  # correct the normalization constants
  m2 <- sum(v1dens^2 * dens)
  S0 <- N  #we know the weight matrix for the individual level should add up to N
  ydens <- S0 * m2
  globalI <- sum(localI * dens * N)/ydens  # formula 10/11
  
  return(list(globalI = globalI, localI = as.numeric(localI), opinion = v1dens_ind, v1dens = v1dens, 
              proximity = wdens))
}



# local
moranIsig <- function(x, y = NULL, proxmat, dens = NULL, N = length(x), sig = TRUE, nsim = 2000) {
  # Adapted from Anselin (1995, eq. 7, 10, 11)
  # https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1538-4632.1995.tb00338.x
  
  # dens: the proportion of individuals in each cell over the district population if individual level
  # data dens is.null and N is simply length of input if we have aggregate data then N should be total
  # population size (or actually just a large number)
  
  # the N is actually arbitrary, thus simply make sufficiently large
  if (is.null(y)) {
    y <- x
  }
  if (is.null(dens)) {
    dens <- rep(1/N, times = N)
  }
  
  # correct scaling of opinions for densities
  v1dens_ind <- rep(x, times = (dens * N))
  v1dens <- (x - mean(v1dens_ind))/sd(v1dens_ind)
  v2dens_ind <- rep(y, times = (dens * N))
  v2dens <- (y - mean(v2dens_ind))/sd(v2dens_ind)
  
  # (density) weighted proximity matrix
  w <- proxmat
  wdens <- t(dens * t(w))
  wdens <- wdens/rowSums(wdens)
  
  # density and proximity weighted locals
  localI <- (v1dens * wdens %*% v2dens)  #formula 7
  
  # correct the normalization constants
  m2 <- sum(v1dens * v1dens * dens)
  S0 <- N  #we know the weight matrix for the individual level should add up to N
  ydens <- S0 * m2
  globalI <- sum(localI * dens * N)/ydens  # formula 10/11
  
  if (sig == TRUE) {
    # start the permutation tests.  how many individuals in each cell?
    ncells <- dens * N
    # assign a groupvariable for individuals to which cell they belong
    cell_id <- rep(1:length(dens), times = ncells)
    MI_sim <- MI_sim2 <- rep(NA, nsim)
    for (i in 1:nsim) {
      # the logic is, to randomly distribute the observed individuals across space
      permutation <- sample(1:N)
      v2dens_ind_new <- v2dens_ind[permutation]
      v1dens_ind_new <- v1dens_ind[permutation]
      # note that we are NOT affecting the density for each cell, the mean and sd of x and y these
      # individuals with there characteristics have to be aggregated again to the cell levels
      groups_cell <- split(v2dens_ind_new, f = cell_id)
      # calculate the mean score per cell
      v2dens_new <- sapply(groups_cell, FUN = mean)
      # make sure to scale again
      v2dens_new <- (v2dens_new - mean(v2dens_ind_new))/sd(v2dens_ind_new)
      groups_cell <- split(v1dens_ind_new, f = cell_id)
      # calculate the mean score per cell
      v1dens_new <- sapply(groups_cell, FUN = mean)
      # make sure to scale again
      v1dens_new <- (v1dens_new - mean(v1dens_ind_new))/sd(v1dens_ind_new)
      
      # calculate simulated MI
      localIsim <- (v1dens_new * wdens %*% v2dens_new)
      # correct the normalization constants
      m2 <- sum(v1dens_new * v1dens_new * dens)
      S0 <- N  #we know the weight matrix for the individual level should add up to N
      ydens <- S0 * m2
      MI_sim[i] <- sum(localIsim * dens * N)/ydens
      MI_sim2[i] <- sum(localIsim)
    }
    # calculate significance
    MI_p2 <- ifelse(sum(localI) > 0, sum(MI_sim2 > sum(localI))/nsim, sum(MI_sim2 < sum(localI))/nsim)
    MI_p <- ifelse(globalI > 0, sum(MI_sim > globalI)/nsim, sum(MI_sim < globalI)/nsim)
  }
  
  
  return(list(globalI = globalI, localI = as.numeric(localI), opinion = v1dens_ind, group = v2dens_ind, 
              group2 = v2dens_ind_new, v1dens = v1dens, proximity = wdens, MI_sim = MI_sim, MI_sim2 = MI_sim2, 
              MI_p = MI_p, MI_p2 = MI_p2))
}

# uni
require(spdep)
print("moranI")
moranI(opinions, proxmat = proximity)
print("moran from spdep")
moran(scale(opinions), n = 100, S0 = sum(proximity), nblist)
print("moran.mc from spdep")
moran.mc(scale(opinions), nblist, nsim = 100)
print("moranIsig")
moranIsig(opinions, proxmat = proximity, nsim = 200)













# fake data cell_level
set.seed(86235)
agent_id <- 1:5
opinions <- runif(5, -1, 1)
groups <- runif(5, 0, 1)  # a percentage

# define coordinates for the cells
x <- matrix(NA, nrow = 5, ncol = 2)
x[1, ] <- c(1, 1)
x[2, ] <- c(1, 2)
x[3, ] <- c(2, 1)
x[4, ] <- c(2, 2)
x[5, ] <- c(3, 2)
distance <- as.matrix(dist(x, diag = TRUE, upper = TRUE))

popsize <- c(10, 20, 30, 40, 50)
popN <- sum(popsize)
densities <- popsize/popN

proximity <- exp(-3 * (distance)/max(distance))
diag(proximity) <- 0

# disaggregate define coordinates for each individual
x_coord <- rep(c(1, 1, 2, 2, 3), times = popsize)
y_coord <- rep(c(1, 2, 1, 2, 2), times = popsize)
x_i <- as.matrix(cbind(x_coord, y_coord))

opinions_i <- rep(opinions, times = popsize)
groups_i_prob <- rep(groups, times = popsize)  #the percentage
groups_i <- NA
for (i in 1:length(groups_i_prob)) {
  groups_i[i] <- sample(c(-1, 1), size = 1, prob = c(groups_i_prob[i], 1 - groups_i_prob[i]))
}

distance_i <- as.matrix(dist(x_i, diag = TRUE, upper = TRUE))

proximity_i <- exp(-3 * (distance_i)/max(distance_i))
proximity_i[proximity_i == 1] <- 0


coords <- coordinates(x_i)

proximity_i <- proximity_i/rowSums(proximity_i)

weightslist <- list()

for (i in 1:nrow(proximity_i)) {
  weightslist[[i]] <- proximity_i[i, ][-i]
}

nbs <- knn2nb(knearneigh(coords, k = popN - 1))

nblist <- nb2listw(neighbours = nbs, weightslist)


# aggregate data, disaggregated, no constant densities
dis <- moranI(opinions_i, proxmat = proximity_i)
dis2 <- moranIsig(opinions_i, proxmat = proximity_i, nsim = 200)
dis3 <- moran.mc(scale(opinions_i), nblist, nsim = 200)


# now on aggregate data
ag <- moranI(opinions, dens = densities, N = popN, proxmat = proximity)
ag2 <- moranIsig(opinions, dens = densities, N = popN, proxmat = proximity, nsim = 200)

dis$globalI



################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


for(c in 1:length(unique(ri$H))){
  temp <- subset(
    ri,
    ri$distanceDecay==2 &
      ri$initialOpinionDistribution=="groupBias" &
      ri$H==unique(ri$H)[c]
  )
  print(violinPlotOld(
    depVar = log10(temp$nIntFirstExtr),
    indepVar = temp$expOutgr2,
    panels = temp$group,
    depVarLabel = "Time of first extremization (log_10)",
    indepVarLabel = "Outgroup exposure (s=100)"
  ) + facet_wrap(cols=vars(temp$group), labeller = as_labeller(c(
    "1"="Group 1: non-western residents",
    "-1" = "Group -1: western residents"
  ))) + ggtitle(paste0("H=",unique(ri$H)[c])))
}
for(c in 1:length(unique(ri$H))){
  temp <- subset(
    ri,
    ri$distanceDecay==2 &
      ri$initialOpinionDistribution=="groupBias" &
      ri$H==unique(ri$H)[c]
  )
  print(violinPlotOld(
    depVar = log10(temp$nIntFirstExtr),
    indepVar = temp$expOutgr2,
    panels = temp$group,
    depVarLabel = "Time of first extremization (log_10)",
    indepVarLabel = "Outgroup exposure (s=100)"
  ) + facet_wrap(cols=vars(temp$group), labeller = as_labeller(c(
    "1"="Group 1: non-western residents",
    "-1" = "Group -1: western residents"
  ))) + ggtitle(paste0("H=",unique(ri$H)[c])))
}














test <- subset(
  ri,
  ( 
    ri$distanceDecay == 2 &#2 &
      ri$initialOpinionDistribution == "beta" &
      ri$H == 0.6
  )
)
par(mfrow=c(4,3))
for (d in 1:12){
  #for (d in c(1,2,11)){
  #temp <- subset(rri, rri$wijk == d)
  temp <- subset(test, test$wijk == d)
  p <-violinPlotOld(
    data = temp,
    depVar = temp$opAlignment2,
    #depVar = abs(temp$opAlignment2),
    indepVar = temp$exposureOutgroup,
    panels = temp$group,
    #bins=5,
    #depVarLabel = "Alignment (agent level, s=100)",
    depVarLabel = "Alignment (raw scores agent level, s=100)",
    indepVarLabel = "Outgroup exposure (s=100)"
  ) + facet_wrap(temp$group, labeller = as_labeller(c(
    "1"="Group 1: non-western",
    "-1" = "Group -1: western"
  ))) + ggtitle(districtsNames[d])
  print(p)
}

par(mfrow=c(2,1))
hist(subset(rri, rri$group==1)$exposureOutgroup, xlim=c(0,1))
hist(subset(rri, rri$group==-1)$exposureOutgroup, xlim=c(0,1))

dia <- subset(di, di$group == -1)
dib <- subset(di, di$group == 1)
par(mfrow=c(2,1))
ggplot(
  dia,
  aes(
    x=cut_interval(x=dia$expOutgr2, length=0.1),
    y=abs(dia$opAlignment2))
) + 
  geom_boxplot() + 
  xlab ("Exposure to outgroup (s=100)") +
  ylab ("Alignment (s=100)") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )
ggplot(
  dib,
  aes(
    x=cut_interval(x=dib$expOutgr2, length=0.1),
    y=abs(dib$opAlignment2))
) + 
  geom_boxplot() + 
  xlab ("Exposure to outgroup (s=100)") +
  ylab ("Alignment (s=100)") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )

violinPlotOld(
  data = di,
  depVar = abs(di$opAlignment2),
  indepVar = di$expOutgr2,
  panels = di$group,
  bins=7,
  depVarLabel = "Alignment (agent level, s=100)",
  indepVarLabel = "Outgroup exposure (s=100)"
) + facet_wrap(di$group, labeller = as_labeller(c(
  "1"="Group 1: non-western",
  "-1" = "Group -1: western"
)))





violinPlotOld(
  depVar = log10(rri$timeFirstPol),
  indepVar = rri$exposureOutgroup,
  panels = rri$wijk,
  depVarLabel = "Time of first extremization (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)",
  bins=5
) + facet_wrap(rri$wijk, labeller = as_labeller(c(
  "1" = "Stadscentrum",
  "2" = "Delfshaven",
  "3" = "Overschie",
  "4" = "Noord",
  "5" = "Hillegersberg-Schiebroek",
  "6" = "Kralingen-Crooswijk",
  "7" = "Feijenoord",
  "8" = "IJsselmonde",
  "9" = "Pernis",
  "10" = "Prins Alexander",
  "11" = "Charlois",
  "12" = "Hoogvliet"
)))



#head(  cut(di$expOutgr2, breaks=(c(0:10)/10))   )
#table(  as.numeric(cut(di$expOutgr2, breaks=(c(0:10)/10)))   )

# Exploring single runs from dataset ab
load(file="./ab.RData")

a[[1]][[200]]$polarizationIndex
b[[1]][[200]]$polarizationIndex
c[[1]][[50]]$polarizationIndex
d[[1]]$polarizationIndex
e[[1]]$polarizationIndex

abs(a[[1]][[200]]$meanOpinionG1 - a[[1]][[200]]$meanOpinionG2)
abs(b[[1]][[200]]$meanOpinionG1 - b[[1]][[200]]$meanOpinionG2)
abs(c[[1]][[50]]$meanOpinionG1 - c[[1]][[50]]$meanOpinionG2)
abs(d[[1]]$meanOpinionG1 - d[[1]]$meanOpinionG2)
abs(e[[1]]$meanOpinionG1 - e[[1]]$meanOpinionG2)
# polarization and strong global alignment in all three runs

ai <- a[[2]][[200]]
bi <- b[[2]][[200]]
ci <- c[[2]][[50]]
di <- d[[2]]
ei <- e[[2]]
plot(di$timeFirstExtr, di$nIntFirstExtr)

table(ai$nIntFirstExtr < ai$timeFirstExtr)
table(bi$nIntFirstExtr < bi$timeFirstExtr)
table(ci$nIntFirstExtr < ci$timeFirstExtr)
table(di$nIntFirstExtr < di$timeFirstExtr)
table(ei$nIntFirstExtr < ei$timeFirstExtr)

# time first extremization
violinPlot(
  data = ai,
  depVar = log10(ai$timeFirstExtr),
  indepVar = ai$expOutgr2,
  depVarLabel = "Time of first extremization (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)"#,
)
violinPlot(
  data = bi,
  depVar = log10(bi$timeFirstExtr),
  indepVar = bi$expOutgr2,
  depVarLabel = "Time of first extremization (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)"#,
)
violinPlot(
  data = ci,
  depVar = log10(ci$timeFirstExtr),
  indepVar = ci$expOutgr2,
  depVarLabel = "Time of first extremization (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)"#,
)
violinPlot(
  data = di,
  depVar = log10(di$timeFirstExtr),
  indepVar = di$expOutgr2,
  depVarLabel = "Time of first extremization (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)"#,
)
violinPlot(
  data = ei,
  depVar = log10(ei$timeFirstExtr),
  indepVar = ei$expOutgr2,
  depVarLabel = "Time of first extremization (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)"#,
)

# N interatctions before first extremization
violinPlot(
  data = ai,
  depVar = log10(ai$nIntFirstExtr),
  indepVar = ai$expOutgr2,
  depVarLabel = "N int. before first extr. (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)",
  fill = "purple"
)
violinPlot(
  data = bi,
  depVar = log10(bi$nIntFirstExtr),
  indepVar = bi$expOutgr2,
  depVarLabel = "N int. before first extr. (log_10))",
  indepVarLabel = "Outgroup exposure (s=100)",
  fill = "purple"
)
violinPlot(
  data = ci,
  depVar = log10(ci$nIntFirstExtr),
  indepVar = ci$expOutgr2,
  depVarLabel = "N int. before first extr. (log_10))",
  indepVarLabel = "Outgroup exposure (s=100)",
  fill = "purple"
)
violinPlot(
  data = di,
  depVar = log10(di$nIntFirstExtr),
  indepVar = di$expOutgr2,
  depVarLabel = "N int. before first extr. (log_10))",
  indepVarLabel = "Outgroup exposure (s=100)",
  fill = "purple"
)

#by group
violinPlot(
  data = ai,
  depVar = log10(ai$timeFirstExtr),
  indepVar = ai$expOutgr2,
  panels = ai$group,
  depVarLabel = "Time of first extremization (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)"
) + facet_wrap(ai$group, labeller = as_labeller(c(
  "1"="Group 1: non-western residents",
  "-1" = "Group -1: western residents"
)))
violinPlot(
  data = bi,
  depVar = log10(bi$timeFirstExtr),
  indepVar = bi$expOutgr2,
  panels = bi$group,
  depVarLabel = "Time of first extremization (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)"
) + facet_wrap(bi$group, labeller = as_labeller(c(
  "1"="Group 1: non-western residents",
  "-1" = "Group -1: western residents"
)))
violinPlot(
  data = ci,
  depVar = log10(ci$timeFirstExtr),
  indepVar = ci$expOutgr2,
  panels = ci$group,
  depVarLabel = "Time of first extremization (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)"
) + facet_wrap(ci$group, labeller = as_labeller(c(
  "1"="Group 1: non-western residents",
  "-1" = "Group -1: western residents"
)))
violinPlot(
  data = di,
  depVar = log10(di$timeFirstExtr),
  indepVar = di$expOutgr2,
  panels = di$group,
  depVarLabel = "Time of first extremization (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)"
) + facet_wrap(di$group, labeller = as_labeller(c(
  "1"="Group 1: non-western residents",
  "-1" = "Group -1: western residents"
)))
violinPlot(
  data = ai,
  depVar = log10(ai$nIntFirstExtr),
  indepVar = ai$expOutgr2,
  panels = ai$group,
  depVarLabel = "N int. before first extr. (log_10))",
  indepVarLabel = "Outgroup exposure (s=100)",
  fill = "purple"
) + facet_wrap(ai$group, labeller = as_labeller(c(
  "1"="Group 1: non-western residents",
  "-1" = "Group -1: western residents"
)))
violinPlot(
  data = bi,
  depVar = log10(bi$nIntFirstExtr),
  indepVar = bi$expOutgr2,
  panels = bi$group,
  depVarLabel = "N int. before first extr. (log_10))",
  indepVarLabel = "Outgroup exposure (s=100)",
  fill = "purple"
) + facet_wrap(bi$group, labeller = as_labeller(c(
  "1"="Group 1: non-western residents",
  "-1" = "Group -1: western residents"
)))
violinPlot(
  data = ci,
  depVar = log10(ci$nIntFirstExtr),
  indepVar = ci$expOutgr2,
  panels = ci$group,
  depVarLabel = "N int. before first extr. (log_10))",
  indepVarLabel = "Outgroup exposure (s=100)",
  fill = "purple"
) + facet_wrap(ci$group, labeller = as_labeller(c(
  "1"="Group 1: non-western residents",
  "-1" = "Group -1: western residents"
)))
violinPlot(
  data = di,
  depVar = log10(di$nIntFirstExtr),
  indepVar = di$expOutgr2,
  panels = di$group,
  depVarLabel = "N int. before first extr. (log_10))",
  indepVarLabel = "Outgroup exposure (s=100)",
  fill = "purple"
) + facet_wrap(di$group, labeller = as_labeller(c(
  "1"="Group 1: non-western residents",
  "-1" = "Group -1: western residents"
)))





#

require(plyr)
require(ggplot2)
#rri$plot <- mapvalues(rri$group, c(-1,1), c("group -1", "b"))
#panelLabel <- labeller(c(
#    "1" = "group 1",
#    "-1" = "group -1"
#  )
#)
violinPlot(
  depVar = log10(rri$timeFirstPol),
  indepVar = rri$exposureOutgroup,
  depVarLabel = "Time of first extremization (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)"#,
)

violinPlot(
  depVar = log10(rri$timeFirstPol),
  indepVar = rri$exposureOutgroup,
  panels = rri$group,
  depVarLabel = "Time of first extremization (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)"
) + facet_wrap(rri$group, labeller = as_labeller(c(
  "1"="Group 1: non-western residents",
  "-1" = "Group -1: western residents"
))) 


violinPlot(
  depVar = log10(rri$timeFirstPol),
  indepVar = rri$exposureOutgroup,
  panels = rri$wijk,
  depVarLabel = "Time of first extremization (log_10)",
  indepVarLabel = "Outgroup exposure (s=100)",
  bins=5
) + facet_wrap(rri$wijk, labeller = as_labeller(c(
  "1" = "Stadscentrum",
  "2" = "Delfshaven",
  "3" = "Overschie",
  "4" = "Noord",
  "5" = "Hillegersberg-Schiebroek",
  "6" = "Kralingen-Crooswijk",
  "7" = "Feijenoord",
  "8" = "IJsselmonde",
  "9" = "Pernis",
  "10" = "Prins Alexander",
  "11" = "Charlois",
  "12" = "Hoogvliet"
)))


load(file= "oldRunsResults.RData")

print(districtViolins(
  depVar="opNetAli2",
  indepVar="expOutgr2",
  depVarLabel="Average local alignment (s=100)",
  indepVarLabel="outgroup exposure (s=100)"))

print(districtViolins(
  depVar="opAlignment2",
  indepVar="expOutgr2",
  depVarLabel="Average local alignment (s=100)",
  indepVarLabel="outgroup exposure (s=100)"))




# old measures of alignment

rr$oldAlignment <- rr$varOpinionGlobal - mean(c(rr$varOpinionG1, rr$varOpinionG2))
hist(rr$oldAlignment)
districtViolins(
  depVar="oldAlignment",
  indepVar="expOutgr2",
  depVarLabel="Alignment (difference between\nglobal- and group-variance",
  indepVarLabel="Outgroup exposure (s=100)")

rr$intuitiveAlignment <- abs(rr$varOpinionG1 - rr$varOpinionG2)
hist(rr$intuitiveAlignment)
districtViolins(
  depVar="intuitiveAlignment",
  indepVar="expOutgr2",
  depVarLabel="Alignment (abs. difference between\nthe average attitude of the two groups",
  indepVarLabel="Outgroup exposure (s=100)")


# We find and save separately the runs that were ran with the 'baseline'
# parameter configuration.
rr <- subset(
  r,
  ( 
    r$distanceDecay == 2 &
      r$initialOpinionDistribution == "groupBias" &
      r$H == 0.6
  )
)
for(i in 1:nrow(rr))(rr$wkName[i]<-districtsNames[rr$wijk[i]])


rri <- subset(
  ri,
  ( 
    ri$distanceDecay == 2 &#2 &
      ri$initialOpinionDistribution == "groupBias" &
      ri$H == 0.6
  )
)




# 1) segregation -> polarization

# Preliminary question: did we get some polarization at all?
# Here we use histograms to plot the distribution of the polarization index
# at the end of the simulation run. We look into the baseline configuration (rr)
# as well as all configurations together (r)
for(i in 1:nrow(r)){if(r$polarizationIndex[i]>1){r$polarizationIndex[i]<-1}}
par(mar=c(4,1,1,0.1))
hist(
  r$polarizationIndex, breaks=10,
  main="All parameter configurations",
  xlab="Polarization index", yaxt="n", ylab="", yaxt='n',
  col="black", border="white"
)
abline(v=0.3, col= "red")
title(ylab="Relative frequency", line=0)

for(i in 1:nrow(rr)){if(rr$polarizationIndex[i]>1){rr$polarizationIndex[i]<-1}}
par(mar=c(4,1,1,0.1))
hist(
  rr$polarizationIndex, breaks=10,
  main="Baseline parameter configuration",
  xlab="Polarization index", yaxt="n", ylab="", yaxt='n',
  col="black", border="white"
)
abline(v=0.3, col= "red")
title(ylab="Relative frequency", line=0)

for(i in 1:nrow(rr)){if(rr$polarizationIndex[i]>1){rr$polarizationIndex[i]<-1}}
par(mar=c(4,1,1,0.1))
hist(
  rri$opinion, breaks=10,
  main="",
  xlab="Attitude", yaxt="n", ylab="", yaxt='n',
  col="black", border="black"
)
title(ylab="Relative frequency", line=0)

par(mar=c(4,1,1,0.1))
hist(
  rri$timeFirstPol, breaks=10,
  main="",
  xlab="Time of first extremization", yaxt="n", ylab="", yaxt='n',
  col="black", border="black"
)
title(ylab="Relative frequency", line=0)

par(mar=c(4,1,1,0.1))
hist(
  rri$expOutgr2, breaks=10,
  main="",
  xlab="Exposure to outgroup (s=100)", yaxt="n", ylab="", yaxt='n',
  col="black", border="white"
)
title(ylab="Relative frequency", line=0)



# Plotting micro-level
#
# First we need to fetch the resources and polygon data.
load("./cityData/importedCBS.RData")

cit <- subset(nedb,nedb$WK_CODE %in% districtsList)
cit@data$WK_CODE <- droplevels(cit@data$WK_CODE)
cit <- spTransform(cit,CRS("+proj=longlat +ellps=WGS84"))

a <- "gray80"
b <- "gray60"
c <- "gray40"
d <- "gray20"
color <- c(a,b,c,d,b,c,b,a,c,d,d,a)
rm(a,b,c,d)

cit@data$COLOR <- NA
for (i in 1:length(cit@data$COLOR)){
  c <- citySummary
  c$color <- color
  l <- subset(c, c$WK_CODE == cit@data$WK_CODE[i])
  cit@data$COLOR[i] <- l$color
  #cit@polygons$COLOR[i] <- l$color
}
par(mar=c(0,0,0,0))
plot(cit, col=cit$COLOR, border=cit$COLOR)

ifelse(
  district == 0,
  dat <- subset(cbs100_rot, cbs100_rot$WK_CODE %in% districtsList),
  dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[district])
)
dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]


dat$tfp <- dat$absop <- #pol <- 
  dat$clus1 <- dat$clus2 <- dat$clus3 <-
  dat$li1 <- dat$ali2 <- dat$ali3 <-
  dat$exp1 <- dat$exp2 <- dat$exp3 <- NA
for (l in 1:length(dat)){ # for every cell
  cell <- subset(w, w$location == dat$OBJECTID[l])
  cell$timeFirstPol[is.na(cell$timeFirstPol)] <- 201
  dat$tfp[l] <- mean(cell$timeFirstPol, na.rm=TRUE)
  dat$absop[l] <- mean(abs(cell$opinion), na.rm=TRUE)
  #pol[l] <- mean(cell$pol, na.rm=TRUE)
  dat$clus1[l] <- mean(cell$opClustering1, na.rm=TRUE)
  dat$clus2[l] <- mean(cell$opClustering2, na.rm=TRUE)
  dat$clus3[l] <- mean(cell$opClustering3, na.rm=TRUE)
  dat$ali1[l] <- mean(abs(cell$opAlignment1), na.rm=TRUE)
  dat$ali2[l] <- mean(abs(cell$opAlignment2), na.rm=TRUE)
  dat$ali3[l] <- mean(abs(cell$opAlignment3), na.rm=TRUE)
  dat$exp2[l] <- mean(cell$exposureOutgroup, na.rm=TRUE)
}

cPalette <- colorRampPalette(c("blue", "red")) # low, high
cPaletteTrio <-colorRampPalette(c("blue", "white", "red"))
cPaletteReverse <- colorRampPalette(c("red", "blue")) # low, high
x <- log(dat$tfp +1, 10) 
dat$colo <- cPaletteReverse(20)[as.numeric(cut(x, breaks=20))]
par(mar=c(0,0,1,0))
plot(cit, col=cit$COLOR, border=cit$COLOR, main="Time of first polarization")
points(dat$l, dat$j, pch=16, col = dat$colo, cex=0.5)
par(mar=c(5,5,5,5))






# Time of first polarization
ggplot(
  rri,
  aes(
    x=cut_interval(x=rri$expOutgr2, length=0.1),
    y=log10(rri$nIntFirstExtr))
) + 
  geom_boxplot() + 
  xlab ("Exposure to outgroup (s=100)") +
  ylab ("Time of first extremization (log_10)") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )

require(dplyr)
require(tidyr)
ggplot(
  rri,
  aes(
    x=cut_interval(x=rri$expOutgr2, length=0.2),
    y=log10(rri$nIntFirstExtr))
) + 
  geom_boxplot() + 
  facet_wrap(rri$wijk) +
  xlab ("Exposure to outgroup (s=100)") +
  ylab ("Time of first extremization (log_10)") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    axis.line = element_line(colour = "black")
  )


# alignment

ggplot(
  rri,
  aes(
    x=cut_interval(x=rri$exposureOutgroup, length=0.1),
    y=abs(rri$opAlignment2))
) + 
  geom_boxplot() + 
  xlab ("Exposure to outgroup (s=100)") +
  ylab ("Alignment (s=100)") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )



# Absolute opinion
dat$colo <- cPalette(20)[as.numeric(cut(dat$absop, breaks=20))]
par(mar=c(0,0,1,0))
plot(cit, col=cit$COLOR, border=cit$COLOR, main="Absolute opinion")
points(dat$l, dat$j, pch=16, col = dat$colo, cex=0.5)
par(mar=c(5,5,5,5))

# Opinion clustering
x<-dat$clus2
for(i in 1:length(x)){
  if(is.na(x[i])){x[i] <- 0}
  if(x[i]>4){x[i]<-4}
}
dat$colo <- cPaletteTrio(20)[as.numeric(cut(x, breaks=seq(-4,4,8/20)))]
par(mar=c(0,0,1,0))
plot(cit, col=cit$COLOR, border=cit$COLOR, main="Opinion clustering (s=100)")
points(dat$l, dat$j, pch=16, col = dat$colo, cex=0.5)
par(mar=c(5,5,5,5))

heatmap(
  x = rri$exposureOutgroup,#rbeta(1000,3,9),
  y = x,#rbeta(1000,9,3),
  type = "proportion",
  bins = 11,
  xlab = "Exposure to outgroup (s=100)",
  ylab = "Opinion clustering (s=100)",
  legend = FALSE
)

# Opinion alignment
x<-dat$ali2
for(i in 1:length(x)){
  if(is.na(x[i])){x[i] <- 0}
  if(x[i]>4){x[i]<-4}
}
dat$colo <- cPalette(20)[as.numeric(cut(x, breaks=20))]
par(mar=c(0,0,1,0))
plot(cit, col=cit$COLOR, border=cit$COLOR, main="Group-attitude alignment (s=100)")
points(dat$l, dat$j, pch=16, col = dat$colo, cex=0.5)
par(mar=c(5,5,5,5))

heatmap(
  x = rri$exposureOutgroup,#rbeta(1000,3,9),
  y = x,#rbeta(1000,9,3),
  type = "count",
  bins = 11,
  xlab = "Exposure to outgroup (s=100)",
  ylab = "Group-attitude alignment (s=100)",
  legend = FALSE
)


# outgroup exposure
x<-dat$exp2
for(i in 1:length(x)){if(is.na(x[i])){x[i] <- 0}}
dat$colo <- cPalette(20)[as.numeric(cut(x, breaks=20))]
par(mar=c(0,0,1,0))
plot(cit, col=cit$COLOR, border=cit$COLOR, main="Outgroup exposure (s=100)")
points(dat$l, dat$j, pch=16, col = dat$colo, cex=0.5)
par(mar=c(5,5,5,5))

# proportion nwa
x<-dat$pnwal2014
#for(i in 1:length(x)){if(is.na(x[i])){x[i] <- 0}}
dat$colo <- cPalette(20)[as.numeric(cut(x, breaks=20))]
par(mar=c(0,0,1,0))
plot(cit, col=cit$COLOR, border=cit$COLOR, main="Proportion non-western residents")
points(dat$l, dat$j, pch=16, col = dat$colo, cex=0.5)
par(mar=c(5,5,5,5))





# clustering - alignment
x<-dat$clus2
for(i in 1:length(x)){
  if(is.na(x[i])){x[i] <- 0}
  if(x[i]>4){x[i]<-4}
}
y<-dat$ali2
for(i in 1:length(y)){
  if(is.na(y[i])){y[i] <- 0}
  if(y[i]>4){y[i]<-4}
}

heatmap(
  x = x,#rbeta(1000,3,9),
  y = y,#rbeta(1000,9,3),
  type = "prop",
  bins = 11,
  xlab = "Attitude clustering (s=100)",
  ylab = "Group-attitude alignment (s=100)",
  legend = FALSE
)



heatmap(
  x = dat$pnwal2014,#rbeta(1000,3,9),
  y = y,#rbeta(1000,9,3),
  type = "prop",
  bins = 11,
  xlab = "Attitude clustering (s=100)",
  ylab = "Group-attitude alignment (s=100)",
  legend = FALSE
)



print(districtViolins(
  depVar="opAlignment1",
  indepVar="expOutgr1",
  depVarLabel="Group-attitude alignment (s=100)",
  indepVarLabel="NWA outgroup exposure (s=100)"))




#plot(rri$exposureOutgroup, log(rri$timeFirstPol + 1,10))
y <- scale(log(rri$timeFirstPol + 1,10))
y <- y + min(y)


for(i in 1:nrow(rri)){if(is.na(rri$timeFirstPol[i])){rri$timeFirstPol[i]<-201}}

heatmap(
  x = rri$exposureOutgroup,#rbeta(1000,3,9),
  y = log(rri$timeFirstPol + 1,10),#rbeta(1000,9,3),
  type = "count",
  bins = 11,
  xlab = "Exposure to outgroup (s=100)",
  ylab = "Time of first polarization (log_10)",
  legend = FALSE
)

for(i in 1:nrow(rri)){if(is.na(rri$grSegregation2[i])){rri$grSegregation2[i]<-0}}
heatmap(
  x = rri$grSegregation2,#rbeta(1000,3,9),
  y = log(rri$timeFirstPol + 1,10),#rbeta(1000,9,3),
  type = "proportion",
  bins = 10,
  xlab = "Group segregation (s=100)",
  ylab = "Time of first polarization (log_10)",
  legend = FALSE
)


x <- log(dat$tfp +1, 10) 
dat$colo <- cPalette(20)[as.numeric(cut(dat$tfp, breaks=20))]
par(mar=c(0,0,1,0))
plot(
  cit, col=cit$COLOR, border=cit$COLOR,
  main="Time of first polarization"
)
points(dat$l, dat$j, pch=20, col = dat$colo, cex=1)


par(mar=c(5,5,5,5))


dat$colo <- cPalette(20)[as.numeric(cut(dat$tfp, breaks=20))]
plot(dat$l, dat$j, cex=2, pch=15, col=dat$colodistrict, main = 
       paste("time of first polarization \nin", placeName))
points(dat$l, dat$j, pch=20, col = dat$colo, cex=1)




dat$tfp <- tfp
dat$absop <- absop
dat$clus <- clus
dat$clus[dat$clus > 1.1] <- 1.1 #####
dat$clus[dat$clus < 0] <- 0############
#dat$pol <-pol
#dat$cclus <- cclus
#dat$aali <- aali
dat$ali <- ali
dat$ali[dat$ali > 1.1] <-1 #############
dat$exp <- exp


#dat$colodistrict <- cPaletteD(12)[as.numeric(cut(as.numeric(dat$WK_CODE), breaks=12))]







# Time first polarization
dat$colo <- cPalette(20)[as.numeric(cut(dat$tfp, breaks=20))]
plot(dat$l, dat$j, cex=2, pch=15, col=dat$colodistrict, main = 
       paste("time of first polarization \nin", placeName))
points(dat$l, dat$j, pch=20, col = dat$colo, cex=1)






# Plotting micro-level
#
# First we need to fetch the resources and polygon data.
load("./cityData/importedCBS.RData")
require(maptools)
require(rgeos)
require(tidyverse)
require(broom)



























district <- 0 # if 0, whole Rotterdam
placeName <- NA
ifelse(district==0, placeName<-"Rotterdam", placeName<-districtsNames[district])
cPalette <- colorRampPalette(c("blue", "red")) # low, high
cPaletteD <- colorRampPalette(c("beige","orange", "yellow","cyan", "green", "azure3", "black")) # low, high
#wijk <- rr[example,]$wijk
#load(paste0("./simOutput/peregrine/", rr[example,]$fileName))
#w <- worldList[[rr$indexParameters[example]]]
#w <- subset(rri, rri$wijk==example)
w <- rri
#w <- s
#safe <- w

polarizedSeeds <- subset(r, r$polarizationIndex > 0.7)$seed
w <- w[w$seed %in% polarizedSeeds,]


#w$cclus <- abs(w$opinion) * w$opClustering2
#w$aali <- abs(w$opinion) * w$opAlignment2
# plotting all districts together
#
ifelse(
  district == 0,
  dat <- subset(cbs100_rot, cbs100_rot$WK_CODE %in% districtsList),
  dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[district])
)
dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]


tfp <- absop <- pol <- clus <- ali <- exp <- c()
for (l in 1:length(dat)){ # for every cell
  cell <- subset(w, w$location == dat$OBJECTID[l])
  cell$timeFirstPol[is.na(cell$timeFirstPol)] <- 201
  tfp[l] <- mean(cell$timeFirstPol, na.rm=TRUE)
  absop[l] <- mean(abs(cell$opinion), na.rm=TRUE)
  #pol[l] <- mean(cell$pol, na.rm=TRUE)
  clus[l] <- mean(cell$opClustering2, na.rm=TRUE)
  ali[l] <- mean(abs(cell$opAlignment2), na.rm=TRUE)
  exp[l] <- mean(cell$exposureOutgroup, na.rm=TRUE)
}

dat$tfp <- tfp
dat$absop <- absop
dat$clus <- clus
dat$clus[dat$clus > 1.1] <- 1.1 #####
dat$clus[dat$clus < 0] <- 0############
#dat$pol <-pol
#dat$cclus <- cclus
#dat$aali <- aali
dat$ali <- ali
dat$ali[dat$ali > 1.1] <-1 #############
dat$exp <- exp


dat$colodistrict <- cPaletteD(12)[as.numeric(cut(as.numeric(dat$WK_CODE), breaks=12))]







# Time first polarization
dat$colo <- cPalette(20)[as.numeric(cut(dat$tfp, breaks=20))]
plot(dat$l, dat$j, cex=2, pch=15, col=dat$colodistrict, main = 
       paste("time of first polarization \nin", placeName))
points(dat$l, dat$j, pch=20, col = dat$colo, cex=1)


#plot(dat$l, dat$j, pch=20, cex=0.5)
#paste("time of first polarization \n", districtsNames[example]))


# Abs opinion
dat$colo <- cPalette(10)[as.numeric(cut(dat$absop, breaks=10))]
plot(dat$l, dat$j, pch=20, col = dat$colo, cex=1, main = 
       paste("abs. opinion \nin", placeName))
#paste("absolute opinion\n", districtsNames[example]))


# op clustering
dat$colo <- cPalette(10)[as.numeric(cut(dat$clus, breaks=10))]
plot(dat$l, dat$j, pch=15, col = dat$colodistrict, cex=2, main = 
       paste("op. clustering (s=100)\nin", placeName))

points(dat$l, dat$j, cex=1, pch=20, col=dat$colo)


# op alignment
dat$colo <- cPalette(20)[as.numeric(cut(dat$ali, breaks=20))]
plot(dat$l, dat$j, pch=15, col = dat$colodistrict, cex=2, main = 
       paste("op. alignment (s=100)\nin", placeName))
points(dat$l, dat$j, pch=20, col = dat$colo, cex=1)


# Abs opinion * clustering
dat$colo <- cPalette(10)[as.numeric(cut(dat$absop * dat$clus, breaks=10))]
plot(dat$l, dat$j, pch=20, col = dat$colo, cex=0.5, main = 
       paste("abs. op * clustering (s=100)\nin", placeName))
#paste("alignment (s=100)\n", districtsNames[example]))


# Abs opinion * clustering * alignment
dat$colo <- cPalette(10)[as.numeric(cut(dat$absop * dat$clus * dat$clus, breaks=10))]
plot(dat$l, dat$j, pch=20, col = dat$colo, cex=0.5, main = 
       paste("abs. op * clust. * align. (s=100)\nin", placeName))
#paste("alignment (s=100)\n", districtsNames[example]))


# Outgroup exposure
dat$colo <- cPalette(10)[as.numeric(cut(dat$exp, breaks=10))]
plot(dat$l, dat$j, pch=20, col = dat$colo, cex=0.5, main = 
       paste("exposure outgroup\n", placeName))


# proportion non-western
dat$colo <- cPalette(10)[as.numeric(cut(dat$pnwal2014, breaks=10))]
plot(dat$l, dat$j, pch=20, col = dat$colo, cex=0.5, main = 
       paste("proportion non-western\nin", placeName))



















# differences between districts
#
# 1) polarization
#boxplot(rr$polarizationIndex ~ rr$wkName)
ggplot(rr, aes(factor(rr$wkName), rr$polarizationIndex)) +
  ylab("Polarization index") +
  geom_violin(scale = "width",draw_quantiles = 0.5, fill="#ff8800") +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )

# 2) Clustering
ggplot(rr, aes(factor(rr$wkName), rr$opNetClus2)) +
  ylab("Attitude clustering (s=100)") +
  geom_violin(scale = "width",draw_quantiles = 0.5, fill="#ff8800") +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )

# 3) alignment
ggplot(rr, aes(factor(rr$wkName), rr$opNetAli2)) +
  ylab("Alignment (s=100)") +
  geom_violin(scale = "width",draw_quantiles = 0.5, fill="#ff8800") +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )





rr <- subset(
  r,
  ( 
    r$distanceDecay == 2 &
      r$initialOpinionDistribution == "groupBias" &
      r$H == 0.6
  )
)
for(i in 1:nrow(rr))(rr$wkName[i]<-districtsNames[rr$wijk[i]])


print(districtViolins(
  depVar="polarizationIndex",
  indepVar="expOutgr1",
  depVarLabel="Polarization index",
  indepVarLabel="exposure to outgroup (s=10)"))

print(districtViolins(
  depVar="polarizationIndex",
  indepVar="expOutgr2",
  depVarLabel="Polarization index",
  indepVarLabel="exposure to outgroup (s=100)"))

print(districtViolins(
  depVar="polarizationIndex",
  indepVar="expOutgr3",
  depVarLabel="Polarization index",
  indepVarLabel="exposure to outgroup (s=1000)"))

print(districtViolins(
  depVar="polarizationIndex",
  indepVar="p_nwa",
  depVarLabel="Polarization index",
  indepVarLabel="proportion of non-western residents"))



print(districtViolins(
  depVar="opNetClus2",
  indepVar="expOutgr2",
  depVarLabel="Opinion clustering (s=100)",
  indepVarLabel="exposure to outgroup (s=100)"))

print(districtViolins(
  depVar="opNetAli2",
  indepVar="expOutgr2",
  depVarLabel="Average local alignment (s=100)",
  indepVarLabel="outgroup exposure (s=100)"))



# plotting micro stats










test<- subset(rr, rr$polarizationIndex > 0.7)

ggplot(test, aes(factor(test$wkName), test$opNetClus2)) +
  ylab("Attitude clustering (s=100)") +
  geom_violin(scale = "width",draw_quantiles = 0.5, fill="#ff8800") +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )





# Regression tables =======================================================

rr <- subset(r,r$distanceDecay==2)

hist(r$polarizationIndex - mean(r$polarizationIndex))
hist(rr$polarizationIndex - mean(rr$polarizationIndex))

hist(r$opNetClus2 - mean(r$opNetClus2))
hist(rr$opNetClus2 - mean(rr$opNetClus2))

hist(r$opNetAli2 - mean(r$opNetAli2))
hist(rr$opNetAli2 - mean(rr$opNetAli2))


# 1) segr -> polarization

summary(lm(formula = rr$polarizationIndex ~ rr$expOutgr2 + 
             as.factor(rr$initialOpinionDistribution) + rr$H ))

summary(lm(formula = rr$polarizationIndex ~ rr$p_nwa +
             rr$grSegregation2 + rr$initialOpinionDistribution + rr$H ))

summary(lm(formula = rr$polarizationIndex ~ rr$p_nwa + 
             rr$grSegregation2 + rr$expOutgr2 +  rr$initialOpinionDistribution + rr$H ))

summary(lm(formula = rr$polarizationIndex ~ as.factor(rr$wijk)))


# 2) segr -> op clustering

summary(lm(formula = rr$opNetClus2 ~ rr$expOutgr2 + 
             as.factor(rr$initialOpinionDistribution) + rr$H ))

summary(lm(formula = rr$opNetClus2 ~ rr$p_nwa +
             rr$grSegregation2 + rr$initialOpinionDistribution + rr$H ))

summary(lm(formula = rr$opNetClus2 ~ rr$p_nwa + 
             rr$grSegregation2 + rr$expOutgr2 +  rr$initialOpinionDistribution + rr$H ))

summary(lm(formula = rr$opNetClus2 ~ as.factor(rr$wijk)))

# 3) segr -> alignment

summary(lm(formula = rr$opNetAli2 ~ rr$expOutgr2 + 
             as.factor(rr$initialOpinionDistribution) + rr$H ))

summary(lm(formula = rr$opNetAli2 ~ rr$p_nwa +
             rr$grSegregation2 + rr$initialOpinionDistribution + rr$H ))

summary(lm(formula = rr$opNetAli2 ~ rr$p_nwa + 
             rr$grSegregation2 + rr$expOutgr2 +  rr$initialOpinionDistribution + rr$H ))

summary(lm(formula = rr$opNetAli2 ~ as.factor(rr$wijk)))

fit<-lm(formula = rr$polarizationIndex ~ rr$p_nwa + 
          rr$grSegregation2 + rr$expOutgr2 +  rr$initialOpinionDistribution + rr$H )
par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(fit)  # Plot the model information




#rr <- subset(GsimResults, GsimResults$initialOpinionDistribu tion!="groupBias")
rr <- subset(rr, rr$distanceDecay == 2)


hist(subset(rr,rr$initialOpinionDistribution=="beta")$polarizationIndex)
hist(subset(rr,rr$initialOpinionDistribution=="uniform")$polarizationIndex)

# 1)
summary(lm(rr$meanAbsOpinion~rr$p_nwa))
summary(lm(rr$meanAbsOpinion~rr$p_nwa + rr$expOutgroup))
summary(lm(
  rr$meanAbsOpinion ~ rr$p_nwa + rr$expOutgroup +
    rr$initialOpinionDistribution +
    rr$H +
    rr$wijk
))
summary(lm(
  rr$polarizationIndex ~ rr$p_nwa + rr$expOutgroup +
    rr$initialOpinionDistribution +
    rr$H +
    rr$wijk
))

# 2)
summary(lm(
  rr$opClusteringA2 ~ rr$p_nwa + rr$expOutgroup +
    rr$initialOpinionDistribution +
    rr$H +
    rr$wijk
))

# 3)
# measuring baseline alignment at t=0
source("./script NI&PA.R")
load("./cityData/geodata_Rotterdam.RData")
d <- ali2 <- c()
for (r in 0:9){
  print(paste0("replication #", r + 1 ," of 10"))
  for(i in 1:12){
    d[12*r + i] <- i
    w <- worldList[[i]]
    w$opinion <- rbeta(nrow(w), 3, 3, ncp = 0)
    w$opinion <- w$opinion * 2 - 1
    dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[i])
    dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
    ops <- c(NA, length = length(dat))
    for (l in 1:length(dat)){ # for every cell
      cell <- subset(w, w$location == dat$OBJECTID[l])
      ops[l] <- mean(cell$opinion)
    }
    dat$opinion <- ops
    ali2[12*r + i] <- mean(abs(MoranI( # global alignment
      x = dat$pnwal2014,
      y = dat$opinion,
      proxmat = proximityList2[[i]],
      type = "local"
    )))
  }
}
test <- as.data.frame(cbind(d,ali2))
for(i in 1:nrow(rr)){
  rr$relativeAlignment[i] <- rr$opAlignment2[i] - mean(subset(test, test$d==rr$wijk[i])$ali2)
}
rr$absoluteAlignment <- rr$opAlignment2
summary(lm(
  rr$relativeAlignment ~ rr$p_nwa + rr$expOutgroup +
    rr$initialOpinionDistribution +
    rr$H +
    rr$wijk
))
summary(lm(
  rr$absoluteAlignment ~ rr$p_nwa + rr$expOutgroup +
    rr$initialOpinionDistribution +
    rr$H +
    rr$wijk
))

r1 <- subset(rr, rr$H == 0.6 & rr$initialOpinionDistribution == "beta")
r2 <- subset(rr, rr$H == 0.6 & rr$initialOpinionDistribution != "beta")
r3 <- subset(rr, rr$H == 0.9 & rr$initialOpinionDistribution == "beta")
r4 <- subset(rr, rr$H == 0.9 & rr$initialOpinionDistribution != "beta")

#0)
png(
  filename = "./outputGraphics/0.png",
  width = 1000,
  height = 800,
  units = "px"
)
par(mfrow=c(4,4))#,oma=c(0,6,0,0))
hist(r1$polarizationIndex, xlab="", ylab="Frequency", cex.lab=1.5, main="H=.6, beta op.distr.",cex.main=1.5, xlim=c(0,1))
mtext("op. polarization", side=1, line=5, cex=1)
hist(r2$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="H=.6, uniform op.distr.", cex.main=1.5, xlim=c(0,1))
mtext("op. polarization", side=1, line=5, cex=1)
hist(r3$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="H=.9, beta op.distr.", cex.main=1.5, xlim=c(0,1))
mtext("op. polarization", side=1, line=5, cex=1)
hist(r4$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="H=.9, uniform op.distr.", cex.main=1.5, xlim=c(0,1))
mtext("op. polarization", side=1, line=5, cex=1)
#op clustering
hist(r1$opClustering2, xlab="", ylab="Frequency", cex.lab=1.5, main="", xlim=c(-0.2,0.2))
mtext("op. clustering 2", side=1, line=3, cex=1)
hist(r2$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", xlim=c(-0.2,0.2))
mtext("op. clustering 2", side=1, line=3, cex=1)
hist(r3$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", xlim=c(-0.2,0.2))
mtext("op. clustering 2", side=1, line=3, cex=1)
hist(r4$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", xlim=c(-0.2,0.2))
#op alignment absolute
hist(r1$absoluteAlignment, xlab="", ylab="Frequency", cex.lab=1.5, main="", xlim=c(-0.3,0.3))
mtext("op. alignment 2", side=1, line=3, cex=1)
hist(r2$absoluteAlignment, xlab="", ylab="", cex.lab=1.5, main="", xlim=c(-0.3,0.3))
mtext("op. alignment 2", side=1, line=3, cex=1)
hist(r3$absoluteAlignment, xlab="", ylab="", cex.lab=1.5, main="", xlim=c(-0.3,0.3))
mtext("op. alignment 2", side=1, line=3, cex=1)
hist(r4$absoluteAlignment, xlab="", ylab="", cex.lab=1.5, main="", xlim=c(-0.3,0.3))
mtext("op. alignment 2", side=1, line=3, cex=1)
#op alignment relative
hist(r1$relativeAlignment, xlab="", ylab="Frequency", cex.lab=1.5, main="", xlim=c(-0.25,0.25))
mtext("relative op. alignment 2", side=1, line=3, cex=1)
hist(r2$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", xlim=c(-0.25,0.25))
mtext("relative op. alignment 2", side=1, line=3, cex=1)
hist(r3$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", xlim=c(-0.25,0.25))
mtext("relative op. alignment 2", side=1, line=3, cex=1)
hist(r4$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", xlim=c(-0.25,0.25))
mtext("relative op. alignment 2", side=1, line=3, cex=1)
dev.off()

#1a)
png(
  filename = "./outputGraphics/1a.png",
  width = 1000,
  height = 1000,
  units = "px"
)
par(mfrow=c(5,4))#,oma=c(0,6,0,0))
# exposure to outgroup
plot(r1$expOutgroup, r1$polarizationIndex, xlab="", ylab="op. polarization", cex.lab=1.5, main="H=.6, beta op.distr.", pch=".")
abline(lm(r1$polarizationIndex~r1$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r2$expOutgroup, r2$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="H=.6, uniform op.distr.", pch=".")
abline(lm(r2$polarizationIndex~r2$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r3$expOutgroup, r3$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="H=.9, beta op.distr.", pch=".")
abline(lm(r3$polarizationIndex~r3$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r4$expOutgroup, r4$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="H=.9, uniform op.distr.", pch=".")
abline(lm(r4$polarizationIndex~r4$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
# proportion non-western
plot(r1$p_nwa, r1$polarizationIndex, xlab="", ylab="op. polarization", cex.lab=1.5, main="", pch=".")
abline(lm(r1$polarizationIndex~r1$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r2$p_nwa, r2$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$polarizationIndex~r2$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r3$p_nwa, r3$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$polarizationIndex~r3$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r4$p_nwa, r4$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$polarizationIndex~r4$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 1)
plot(r1$grSegregation1, r1$polarizationIndex, xlab="", ylab="op. polarization", cex.lab=1.5, main="", pch=".")
abline(lm(r1$polarizationIndex~r1$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r2$grSegregation1, r2$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$polarizationIndex~r2$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r3$grSegregation1, r3$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$polarizationIndex~r3$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r4$grSegregation1, r4$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$polarizationIndex~r4$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 2)
plot(r1$grSegregation2, r1$polarizationIndex, xlab="", ylab="op. polarization", cex.lab=1.5, main="", pch=".")
abline(lm(r1$polarizationIndex~r1$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r2$grSegregation2, r2$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$polarizationIndex~r2$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r3$grSegregation2, r3$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$polarizationIndex~r3$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r4$grSegregation2, r4$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$polarizationIndex~r4$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 3)
plot(r1$grSegregation3, r1$polarizationIndex, xlab="", ylab="op. polarization", cex.lab=1.5, main="", pch=".")
abline(lm(r1$polarizationIndex~r1$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r2$grSegregation3, r2$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$polarizationIndex~r2$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r3$grSegregation3, r3$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$polarizationIndex~r3$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r4$grSegregation3, r4$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$polarizationIndex~r4$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
dev.off()
#1b)
png(
  filename = "./outputGraphics/1b.png",
  width = 1000,
  height = 1000,
  units = "px"
)
par(mfrow=c(5,4))#,oma=c(0,6,0,0))
# exposure to outgroup
plot(r1$expOutgroup, r1$meanAbsOpinion, xlab="", ylab="mean abs. op.", cex.lab=1.5, main="H=.6, beta op.distr.", pch=".")
abline(lm(r1$meanAbsOpinion~r1$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r2$expOutgroup, r2$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="H=.6, uniform op.distr.", pch=".")
abline(lm(r2$meanAbsOpinion~r2$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r3$expOutgroup, r3$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="H=.9, beta op.distr.", pch=".")
abline(lm(r3$meanAbsOpinion~r3$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r4$expOutgroup, r4$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="H=.9, uniform op.distr.", pch=".")
abline(lm(r4$meanAbsOpinion~r4$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
# proportion non-western
plot(r1$p_nwa, r1$meanAbsOpinion, xlab="", ylab="mean abs. op.", cex.lab=1.5, main="", pch=".")
abline(lm(r1$meanAbsOpinion~r1$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r2$p_nwa, r2$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$meanAbsOpinion~r2$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r3$p_nwa, r3$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$meanAbsOpinion~r3$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r4$p_nwa, r4$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$meanAbsOpinion~r4$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 1)
plot(r1$grSegregation1, r1$meanAbsOpinion, xlab="", ylab="mean abs. op.", cex.lab=1.5, main="", pch=".")
abline(lm(r1$meanAbsOpinion~r1$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r2$grSegregation1, r2$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$meanAbsOpinion~r2$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r3$grSegregation1, r3$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$meanAbsOpinion~r3$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r4$grSegregation1, r4$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$meanAbsOpinion~r4$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 2)
plot(r1$grSegregation2, r1$meanAbsOpinion, xlab="", ylab="mean abs. op.", cex.lab=1.5, main="", pch=".")
abline(lm(r1$meanAbsOpinion~r1$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r2$grSegregation2, r2$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$meanAbsOpinion~r2$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r3$grSegregation2, r3$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$meanAbsOpinion~r3$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r4$grSegregation2, r4$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$meanAbsOpinion~r4$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 3)
plot(r1$grSegregation3, r1$meanAbsOpinion, xlab="", ylab="mean abs. op.", cex.lab=1.5, main="", pch=".")
abline(lm(r1$meanAbsOpinion~r1$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r2$grSegregation3, r2$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$meanAbsOpinion~r2$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r3$grSegregation3, r3$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$meanAbsOpinion~r3$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r4$grSegregation3, r4$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$meanAbsOpinion~r4$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
dev.off()

#1c)
png(
  filename = "./outputGraphics/1c.png",
  width = 1000,
  height = 1000,
  units = "px"
)
par(mfrow=c(5,4))#,oma=c(0,6,0,0))
# exposure to outgroup
plot(r1$expOutgroup, r1$meanOpinionGlobal, xlab="", ylab="mean op.", cex.lab=1.5, main="H=.6, beta op.distr.", pch=".")
abline(lm(r1$meanOpinionGlobal~r1$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r2$expOutgroup, r2$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="H=.6, uniform op.distr.", pch=".")
abline(lm(r2$meanOpinionGlobal~r2$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r3$expOutgroup, r3$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="H=.9, beta op.distr.", pch=".")
abline(lm(r3$meanOpinionGlobal~r3$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r4$expOutgroup, r4$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="H=.9, uniform op.distr.", pch=".")
abline(lm(r4$meanOpinionGlobal~r4$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
# proportion non-western
plot(r1$p_nwa, r1$meanOpinionGlobal, xlab="", ylab="mean op.", cex.lab=1.5, main="", pch=".")
abline(lm(r1$meanOpinionGlobal~r1$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r2$p_nwa, r2$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$meanOpinionGlobal~r2$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r3$p_nwa, r3$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$meanOpinionGlobal~r3$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r4$p_nwa, r4$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$meanOpinionGlobal~r4$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 1)
plot(r1$grSegregation1, r1$meanOpinionGlobal, xlab="", ylab="mean op.", cex.lab=1.5, main="", pch=".")
abline(lm(r1$meanOpinionGlobal~r1$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r2$grSegregation1, r2$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$meanOpinionGlobal~r2$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r3$grSegregation1, r3$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$meanOpinionGlobal~r3$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r4$grSegregation1, r4$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$meanOpinionGlobal~r4$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 2)
plot(r1$grSegregation2, r1$meanOpinionGlobal, xlab="", ylab="mean op.", cex.lab=1.5, main="", pch=".")
abline(lm(r1$meanOpinionGlobal~r1$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r2$grSegregation2, r2$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$meanOpinionGlobal~r2$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r3$grSegregation2, r3$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$meanOpinionGlobal~r3$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r4$grSegregation2, r4$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$meanOpinionGlobal~r4$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 3)
plot(r1$grSegregation3, r1$meanOpinionGlobal, xlab="", ylab="mean op.", cex.lab=1.5, main="", pch=".")
abline(lm(r1$meanOpinionGlobal~r1$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r2$grSegregation3, r2$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$meanOpinionGlobal~r2$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r3$grSegregation3, r3$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$meanOpinionGlobal~r3$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r4$grSegregation3, r4$meanOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$meanOpinionGlobal~r4$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
dev.off()

#2)
png(
  filename = "./outputGraphics/2.png",
  width = 1000,
  height = 1000,
  units = "px"
)
par(mfrow=c(5,4))#,oma=c(0,6,0,0))
# exposure to outgroup
plot(r1$expOutgroup, r1$opClustering2 , xlab="", ylab="op. clustering", cex.lab=1.5, main="H=.6, beta op.distr.", pch=".")
abline(lm(r1$opClustering2~r1$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r2$expOutgroup, r2$opClustering2, xlab="", ylab="", cex.lab=1.5, main="H=.6, uniform op.distr.", pch=".")
abline(lm(r2$opClustering2~r2$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r3$expOutgroup, r3$opClustering2, xlab="", ylab="", cex.lab=1.5, main="H=.9, beta op.distr.", pch=".")
abline(lm(r3$opClustering2~r3$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r4$expOutgroup, r4$opClustering2, xlab="", ylab="", cex.lab=1.5, main="H=.9, uniform op.distr.", pch=".")
abline(lm(r4$opClustering2~r4$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
# NWAexpout
plot(r1$NWAexpOutgroup, r1$opClustering2, xlab="", ylab="op. clustering", cex.lab=1.5, main="", pch=".")
abline(lm(r1$opClustering2~r1$NWAexpOutgroup), col="blue")
mtext("exp.out. of minorities", side=1, line=2, cex=1)
plot(r2$NWAexpOutgroup, r2$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$opClustering2~r2$NWAexpOutgroup), col="blue")
mtext("exp.out. of minorities", side=1, line=2, cex=1)
plot(r3$NWAexpOutgroup, r3$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$opClustering2~r3$NWAexpOutgroup), col="blue")
mtext("exp.out. of minorities", side=1, line=2, cex=1)
plot(r4$NWAexpOutgroup, r4$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$opClustering2~r4$NWAexpOutgroup), col="blue")
mtext("exp.out. of minorities", side=1, line=2, cex=1)
# NATexpout
plot(r1$NATexpOutgroup, r1$opClustering2, xlab="", ylab="op. clustering", cex.lab=1.5, main="", pch=".")
abline(lm(r1$opClustering2~r1$NATexpOutgroup), col="blue")
mtext("exp.out. of majority", side=1, line=2, cex=1)
plot(r2$NATexpOutgroup, r2$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$opClustering2~r2$NATexpOutgroup), col="blue")
mtext("exp.out. of majority", side=1, line=2, cex=1)
plot(r3$NATexpOutgroup, r3$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$opClustering2~r3$NATexpOutgroup), col="blue")
mtext("exp.out. of majority", side=1, line=2, cex=1)
plot(r4$NATexpOutgroup, r4$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$opClustering2~r4$NATexpOutgroup), col="blue")
mtext("exp.out. of majority", side=1, line=2, cex=1)
# proportion non-western
plot(r1$p_nwa, r1$opClustering2, xlab="", ylab="op. clustering", cex.lab=1.5, main="", pch=".")
abline(lm(r1$opClustering2~r1$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r2$p_nwa, r2$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$opClustering2~r2$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r3$p_nwa, r3$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$opClustering2~r3$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r4$p_nwa, r4$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$opClustering2~r4$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 2)
plot(r1$grSegregation2, r1$opClustering2, xlab="", ylab="op. clustering", cex.lab=1.5, main="", pch=".")
abline(lm(r1$opClustering2~r1$grSegregation2), col="blue")
mtext("segregation", side=1, line=2, cex=1)
plot(r2$grSegregation2, r2$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$opClustering2~r2$grSegregation2), col="blue")
mtext("segregation", side=1, line=2, cex=1)
plot(r3$grSegregation2, r3$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$opClustering2~r3$grSegregation2), col="blue")
mtext("segregation", side=1, line=2, cex=1)
plot(r4$grSegregation2, r4$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$opClustering2~r4$grSegregation2), col="blue")
mtext("segregation", side=1, line=2, cex=1)
dev.off()


# 3) Segregation -> opinion alignment ---------
png(
  filename = "./outputGraphics/3.png",
  width = 1000,
  height = 1000,
  units = "px"
)
par(mfrow=c(5,4))#,oma=c(0,6,0,0))
plot(r1$expOutgroup, r1$relativeAlignment , xlab="", ylab="op. clustering", cex.lab=1.5, main="H=.6, beta op.distr.", pch=".")
abline(lm(r1$relativeAlignment~r1$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r2$expOutgroup, r2$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="H=.6, uniform op.distr.", pch=".")
abline(lm(r2$relativeAlignment~r2$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r3$expOutgroup, r3$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="H=.9, beta op.distr.", pch=".")
abline(lm(r3$relativeAlignment~r3$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r4$expOutgroup, r4$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="H=.9, uniform op.distr.", pch=".")
abline(lm(r4$relativeAlignment~r4$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
# NWAexpout
plot(r1$NWAexpOutgroup, r1$relativeAlignment, xlab="", ylab="op. alignment", cex.lab=1.5, main="", pch=".")
abline(lm(r1$relativeAlignment~r1$NWAexpOutgroup), col="blue")
mtext("exp.out. of minorities", side=1, line=2, cex=1)
plot(r2$NWAexpOutgroup, r2$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$relativeAlignment~r2$NWAexpOutgroup), col="blue")
mtext("exp.out. of minorities", side=1, line=2, cex=1)
plot(r3$NWAexpOutgroup, r3$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$relativeAlignment~r3$NWAexpOutgroup), col="blue")
mtext("exp.out. of minorities", side=1, line=2, cex=1)
plot(r4$NWAexpOutgroup, r4$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$relativeAlignment~r4$NWAexpOutgroup), col="blue")
mtext("exp.out. of minorities", side=1, line=2, cex=1)
# NATexpout
plot(r1$NATexpOutgroup, r1$relativeAlignment, xlab="", ylab="op. alignment", cex.lab=1.5, main="", pch=".")
abline(lm(r1$relativeAlignment~r1$NATexpOutgroup), col="blue")
mtext("exp.out. of majority", side=1, line=2, cex=1)
plot(r2$NATexpOutgroup, r2$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$relativeAlignment~r2$NATexpOutgroup), col="blue")
mtext("exp.out. of majority", side=1, line=2, cex=1)
plot(r3$NATexpOutgroup, r3$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$relativeAlignment~r3$NATexpOutgroup), col="blue")
mtext("exp.out. of majority", side=1, line=2, cex=1)
plot(r4$NATexpOutgroup, r4$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$relativeAlignment~r4$NATexpOutgroup), col="blue")
mtext("exp.out. of majority", side=1, line=2, cex=1)
# proportion non-western
plot(r1$p_nwa, r1$relativeAlignment, xlab="", ylab="op. alignment", cex.lab=1.5, main="", pch=".")
abline(lm(r1$relativeAlignment~r1$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r2$p_nwa, r2$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$relativeAlignment~r2$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r3$p_nwa, r3$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$relativeAlignment~r3$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r4$p_nwa, r4$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$relativeAlignment~r4$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 2)
plot(r1$grSegregation2, r1$relativeAlignment, xlab="", ylab="op. alignment", cex.lab=1.5, main="", pch=".")
abline(lm(r1$relativeAlignment~r1$grSegregation2), col="blue")
mtext("segregation", side=1, line=2, cex=1)
plot(r2$grSegregation2, r2$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$relativeAlignment~r2$grSegregation2), col="blue")
mtext("segregation", side=1, line=2, cex=1)
plot(r3$grSegregation2, r3$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$relativeAlignment~r3$grSegregation2), col="blue")
mtext("segregation", side=1, line=2, cex=1)
plot(r4$grSegregation2, r4$relativeAlignment, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$relativeAlignment~r4$grSegregation2), col="blue")
mtext("segregation", side=1, line=2, cex=1)
dev.off()















r <- subset(GsimResults, GsimResults$initialOpinionDistribution=="beta")
r1 <- subset(r, r$H == 0.6 & r$distanceDecay == 1)
r2 <- subset(r, r$H == 0.6 & r$distanceDecay == 2)
r3 <- subset(r, r$H == 0.6 & r$distanceDecay == 3)
r4 <- subset(r, r$H == 0.9 & r$distanceDecay == 1)
r5 <- subset(r, r$H == 0.9 & r$distanceDecay == 2)
r6 <- subset(r, r$H == 0.9 & r$distanceDecay == 3)

# 0) Descriptives   ---------
png(
  filename = "./outputGraphics/macroDescriptives1.png",
  width = 1000,
  height = 800,
  units = "px"
)
par(mfrow=c(4,6))#,oma=c(0,6,0,0))
# abs difference between groups' mean opinion
hist(abs(r1$meanOpinionG1 - r1$meanOpinionG2), xlab="", ylab="Frequency", cex.lab=1.5, main="H=.6, d-decay function 1")
mtext("Between-group\nop. difference", side=1, line=5, cex=1)
hist(abs(r2$meanOpinionG1 - r2$meanOpinionG2), xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 2", cex.main=1.5)
mtext("Between-group\nop. difference", side=1, line=5, cex=1)
hist(abs(r3$meanOpinionG1 - r3$meanOpinionG2), xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 3", cex.main=1.5)
mtext("Between-group\nop. difference", side=1, line=5, cex=1)
hist(abs(r4$meanOpinionG1 - r4$meanOpinionG2), xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 1", cex.main=1.5)
mtext("Between-group\nop. difference", side=1, line=5, cex=1)
hist(abs(r5$meanOpinionG1 - r5$meanOpinionG2), xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 2", cex.main=1.5)
mtext("Between-group\nop. difference", side=1, line=5, cex=1)
hist(abs(r6$meanOpinionG1 - r6$meanOpinionG2), xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 3", cex.main=1.5)
mtext("Between-group\nop. difference", side=1, line=5, cex=1)
# op polarization
hist(r1$polarizationIndex, xlab="", ylab="Frequency", cex.lab=1.5, main="")
mtext("op. polarization", side=1, line=3, cex=1)
hist(r2$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. polarization", side=1, line=3, cex=1)
hist(r3$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. polarization", side=1, line=3, cex=1)
hist(r4$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. polarization", side=1, line=3, cex=1)
hist(r5$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. polarization", side=1, line=3, cex=1)
hist(r6$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. polarization", side=1, line=3, cex=1)
# opinion clustering
hist(r1$opClustering2, xlab="", ylab="Frequency", cex.lab=1.5, main="")
mtext("op. clustering 2", side=1, line=3, cex=1)
hist(r2$opClustering2, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. clustering 2", side=1, line=3, cex=1)
hist(r3$opClustering2, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. clustering 2", side=1, line=3, cex=1)
hist(r4$opClustering2, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. clustering 2", side=1, line=3, cex=1)
hist(r5$opClustering2, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. clustering 2", side=1, line=3, cex=1)
hist(r6$opClustering2, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. clustering 2", side=1, line=3, cex=1)
# alignment
hist(r1$opAlignment2, xlab="", ylab="Frequency", cex.lab=1.5, main="")
mtext("op. alignment 2", side=1, line=3, cex=1)
hist(r2$opAlignment2, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. alignment 2", side=1, line=3, cex=1)
hist(r3$opAlignment2, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. alignment 2", side=1, line=3, cex=1)
hist(r4$opAlignment2, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. alignment 2", side=1, line=3, cex=1)
hist(r5$opAlignment2, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. alignment 2", side=1, line=3, cex=1)
hist(r6$opAlignment2, xlab="", ylab="", cex.lab=1.5, main="")
mtext("op. alignment 2", side=1, line=3, cex=1)
dev.off()


# 1) Segregation -> polarization ---------
png(
  filename = "./outputGraphics/macroResults 1a - segr-pol.png",
  width = 1000,
  height = 1000,
  units = "px"
)
par(mfrow=c(5,6))#,oma=c(0,6,0,0))
# exposure to outgroup
plot(r1$expOutgroup, r1$polarizationIndex, xlab="", ylab="op. polarization", cex.lab=1.5, main="H=.6, d-decay function 1", pch=".")
abline(lm(r1$polarizationIndex~r1$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r2$expOutgroup, r2$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 2", pch=".")
abline(lm(r2$polarizationIndex~r2$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r3$expOutgroup, r3$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 3", pch=".")
abline(lm(r3$polarizationIndex~r3$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r4$expOutgroup, r4$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 1", pch=".")
abline(lm(r4$polarizationIndex~r4$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r5$expOutgroup, r5$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 2", pch=".")
abline(lm(r5$polarizationIndex~r5$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r6$expOutgroup, r6$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 3", pch=".")
abline(lm(r6$polarizationIndex~r6$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
# proportion non-western
plot(r1$p_nwa, r1$polarizationIndex, xlab="", ylab="op. polarization", cex.lab=1.5, main="", pch=".")
abline(lm(r1$polarizationIndex~r1$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r2$p_nwa, r2$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$polarizationIndex~r2$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r3$p_nwa, r3$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$polarizationIndex~r3$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r4$p_nwa, r4$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$polarizationIndex~r4$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r5$p_nwa, r5$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$polarizationIndex~r5$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r6$p_nwa, r6$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$polarizationIndex~r6$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 1)
plot(r1$grSegregation1, r1$polarizationIndex, xlab="", ylab="op. polarization", cex.lab=1.5, main="", pch=".")
abline(lm(r1$polarizationIndex~r1$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r2$grSegregation1, r2$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$polarizationIndex~r2$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r3$grSegregation1, r3$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$polarizationIndex~r3$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r4$grSegregation1, r4$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$polarizationIndex~r4$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r5$grSegregation1, r5$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$polarizationIndex~r5$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r6$grSegregation1, r6$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$polarizationIndex~r6$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 2)
plot(r1$grSegregation2, r1$polarizationIndex, xlab="", ylab="op. polarization", cex.lab=1.5, main="", pch=".")
abline(lm(r1$polarizationIndex~r1$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r2$grSegregation2, r2$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$polarizationIndex~r2$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r3$grSegregation2, r3$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$polarizationIndex~r3$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r4$grSegregation2, r4$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$polarizationIndex~r4$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r5$grSegregation2, r5$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$polarizationIndex~r5$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r6$grSegregation2, r6$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$polarizationIndex~r6$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 3)
plot(r1$grSegregation3, r1$polarizationIndex, xlab="", ylab="op. polarization", cex.lab=1.5, main="", pch=".")
abline(lm(r1$polarizationIndex~r1$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r2$grSegregation3, r2$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$polarizationIndex~r2$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r3$grSegregation3, r3$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$polarizationIndex~r3$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r4$grSegregation3, r4$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$polarizationIndex~r4$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r5$grSegregation3, r5$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$polarizationIndex~r5$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r6$grSegregation3, r6$polarizationIndex, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$polarizationIndex~r6$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
dev.off()

#
png(
  filename = "./outputGraphics/macroResults 1b - segr-absOp.png",
  width = 1000,
  height = 1000,
  units = "px"
)
par(mfrow=c(5,6))#,oma=c(0,6,0,0))
# exposure to outgroup
plot(r1$expOutgroup, r1$meanAbsOpinion, xlab="", ylab="mean abs. op.", cex.lab=1.5, main="H=.6, d-decay function 1", pch=".")
abline(lm(r1$meanAbsOpinion~r1$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r2$expOutgroup, r2$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 2", pch=".")
abline(lm(r2$meanAbsOpinion~r2$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r3$expOutgroup, r3$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 3", pch=".")
abline(lm(r3$meanAbsOpinion~r3$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r4$expOutgroup, r4$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 1", pch=".")
abline(lm(r4$meanAbsOpinion~r4$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r5$expOutgroup, r5$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 2", pch=".")
abline(lm(r5$meanAbsOpinion~r5$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r6$expOutgroup, r6$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 3", pch=".")
abline(lm(r6$meanAbsOpinion~r6$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
# proportion non-western
plot(r1$p_nwa, r1$meanAbsOpinion, xlab="", ylab="mean abs. op.", cex.lab=1.5, main="", pch=".")
abline(lm(r1$meanAbsOpinion~r1$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r2$p_nwa, r2$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$meanAbsOpinion~r2$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r3$p_nwa, r3$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$meanAbsOpinion~r3$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r4$p_nwa, r4$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$meanAbsOpinion~r4$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r5$p_nwa, r5$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$meanAbsOpinion~r5$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r6$p_nwa, r6$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$meanAbsOpinion~r6$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 1)
plot(r1$grSegregation1, r1$meanAbsOpinion, xlab="", ylab="mean abs. op.", cex.lab=1.5, main="", pch=".")
abline(lm(r1$meanAbsOpinion~r1$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r2$grSegregation1, r2$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$meanAbsOpinion~r2$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r3$grSegregation1, r3$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$meanAbsOpinion~r3$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r4$grSegregation1, r4$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$meanAbsOpinion~r4$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r5$grSegregation1, r5$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$meanAbsOpinion~r5$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r6$grSegregation1, r6$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$meanAbsOpinion~r6$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 2)
plot(r1$grSegregation2, r1$meanAbsOpinion, xlab="", ylab="mean abs. op.", cex.lab=1.5, main="", pch=".")
abline(lm(r1$meanAbsOpinion~r1$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r2$grSegregation2, r2$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$meanAbsOpinion~r2$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r3$grSegregation2, r3$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$meanAbsOpinion~r3$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r4$grSegregation2, r4$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$meanAbsOpinion~r4$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r5$grSegregation2, r5$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$meanAbsOpinion~r5$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r6$grSegregation2, r6$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$meanAbsOpinion~r6$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 3)
plot(r1$grSegregation3, r1$meanAbsOpinion, xlab="", ylab="mean abs. op.", cex.lab=1.5, main="", pch=".")
abline(lm(r1$meanAbsOpinion~r1$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r2$grSegregation3, r2$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$meanAbsOpinion~r2$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r3$grSegregation3, r3$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$meanAbsOpinion~r3$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r4$grSegregation3, r4$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$meanAbsOpinion~r4$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r5$grSegregation3, r5$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$meanAbsOpinion~r5$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r6$grSegregation3, r6$meanAbsOpinion, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$meanAbsOpinion~r6$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
dev.off()




#
png(
  filename = "./outputGraphics/macroResults 1c - segr-propExtremists.png",
  width = 1000,
  height = 1000,
  units = "px"
)
par(mfrow=c(5,6))#,oma=c(0,6,0,0))
# exposure to outgroup
plot(r1$expOutgroup, r1$propExtremists, xlab="", ylab="prop. extremists", cex.lab=1.5, main="H=.6, d-decay function 1", pch=".")
abline(lm(r1$propExtremists~r1$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r2$expOutgroup, r2$propExtremists, xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 2", pch=".")
abline(lm(r2$propExtremists~r2$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r3$expOutgroup, r3$propExtremists, xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 3", pch=".")
abline(lm(r3$propExtremists~r3$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r4$expOutgroup, r4$propExtremists, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 1", pch=".")
abline(lm(r4$propExtremists~r4$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r5$expOutgroup, r5$propExtremists, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 2", pch=".")
abline(lm(r5$propExtremists~r5$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r6$expOutgroup, r6$propExtremists, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 3", pch=".")
abline(lm(r6$propExtremists~r6$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
# proportion non-western
plot(r1$p_nwa, r1$propExtremists, xlab="", ylab="prop. extremists", cex.lab=1.5, main="", pch=".")
abline(lm(r1$propExtremists~r1$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r2$p_nwa, r2$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$propExtremists~r2$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r3$p_nwa, r3$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$propExtremists~r3$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r4$p_nwa, r4$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$propExtremists~r4$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r5$p_nwa, r5$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$propExtremists~r5$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r6$p_nwa, r6$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$propExtremists~r6$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 1)
plot(r1$grSegregation1, r1$propExtremists, xlab="", ylab="prop. extremists", cex.lab=1.5, main="", pch=".")
abline(lm(r1$propExtremists~r1$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r2$grSegregation1, r2$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$propExtremists~r2$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r3$grSegregation1, r3$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$propExtremists~r3$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r4$grSegregation1, r4$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$propExtremists~r4$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r5$grSegregation1, r5$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$propExtremists~r5$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r6$grSegregation1, r6$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$propExtremists~r6$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 2)
plot(r1$grSegregation2, r1$propExtremists, xlab="", ylab="prop. extremists", cex.lab=1.5, main="", pch=".")
abline(lm(r1$propExtremists~r1$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r2$grSegregation2, r2$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$propExtremists~r2$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r3$grSegregation2, r3$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$propExtremists~r3$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r4$grSegregation2, r4$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$propExtremists~r4$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r5$grSegregation2, r5$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$propExtremists~r5$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r6$grSegregation2, r6$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$propExtremists~r6$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 3)
plot(r1$grSegregation3, r1$propExtremists, xlab="", ylab="prop. extremists", cex.lab=1.5, main="", pch=".")
abline(lm(r1$propExtremists~r1$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r2$grSegregation3, r2$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$propExtremists~r2$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r3$grSegregation3, r3$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$propExtremists~r3$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r4$grSegregation3, r4$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$propExtremists~r4$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r5$grSegregation3, r5$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$propExtremists~r5$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r6$grSegregation3, r6$propExtremists, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$propExtremists~r6$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
dev.off()

#
png(
  filename = "./outputGraphics/macroResults 1d - segr-opVar.png",
  width = 1000,
  height = 1000,
  units = "px"
)
par(mfrow=c(5,6))#,oma=c(0,6,0,0))
# exposure to outgroup
plot(r1$expOutgroup, r1$varOpinionGlobal, xlab="", ylab="op. variance", cex.lab=1.5, main="H=.6, d-decay function 1", pch=".")
abline(lm(r1$varOpinionGlobal~r1$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r2$expOutgroup, r2$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 2", pch=".")
abline(lm(r2$varOpinionGlobal~r2$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r3$expOutgroup, r3$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 3", pch=".")
abline(lm(r3$varOpinionGlobal~r3$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r4$expOutgroup, r4$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 1", pch=".")
abline(lm(r4$varOpinionGlobal~r4$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r5$expOutgroup, r5$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 2", pch=".")
abline(lm(r5$varOpinionGlobal~r5$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r6$expOutgroup, r6$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 3", pch=".")
abline(lm(r6$varOpinionGlobal~r6$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
# proportion non-western
plot(r1$p_nwa, r1$varOpinionGlobal, xlab="", ylab="op. variance", cex.lab=1.5, main="", pch=".")
abline(lm(r1$varOpinionGlobal~r1$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r2$p_nwa, r2$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$varOpinionGlobal~r2$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r3$p_nwa, r3$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$varOpinionGlobal~r3$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r4$p_nwa, r4$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$varOpinionGlobal~r4$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r5$p_nwa, r5$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$varOpinionGlobal~r5$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r6$p_nwa, r6$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$varOpinionGlobal~r6$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 1)
plot(r1$grSegregation1, r1$varOpinionGlobal, xlab="", ylab="op. variance", cex.lab=1.5, main="", pch=".")
abline(lm(r1$varOpinionGlobal~r1$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r2$grSegregation1, r2$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$varOpinionGlobal~r2$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r3$grSegregation1, r3$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$varOpinionGlobal~r3$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r4$grSegregation1, r4$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$varOpinionGlobal~r4$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r5$grSegregation1, r5$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$varOpinionGlobal~r5$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r6$grSegregation1, r6$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$varOpinionGlobal~r6$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 2)
plot(r1$grSegregation2, r1$varOpinionGlobal, xlab="", ylab="op. variance", cex.lab=1.5, main="", pch=".")
abline(lm(r1$varOpinionGlobal~r1$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r2$grSegregation2, r2$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$varOpinionGlobal~r2$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r3$grSegregation2, r3$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$varOpinionGlobal~r3$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r4$grSegregation2, r4$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$varOpinionGlobal~r4$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r5$grSegregation2, r5$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$varOpinionGlobal~r5$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r6$grSegregation2, r6$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$varOpinionGlobal~r6$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 3)
plot(r1$grSegregation3, r1$varOpinionGlobal, xlab="", ylab="op. variance", cex.lab=1.5, main="", pch=".")
abline(lm(r1$varOpinionGlobal~r1$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r2$grSegregation3, r2$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$varOpinionGlobal~r2$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r3$grSegregation3, r3$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$varOpinionGlobal~r3$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r4$grSegregation3, r4$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$varOpinionGlobal~r4$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r5$grSegregation3, r5$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$varOpinionGlobal~r5$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r6$grSegregation3, r6$varOpinionGlobal, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$varOpinionGlobal~r6$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
dev.off()







# 2) Segregation -> opinion clustering ---------
png(
  filename = "./outputGraphics/macroResults 2 - segr-clu.png",
  width = 1000,
  height = 1000,
  units = "px"
)
par(mfrow=c(5,6))#,oma=c(0,6,0,0))
# exposure to outgroup
plot(r1$expOutgroup, r1$opClustering1 , xlab="", ylab="op. clustering 1", cex.lab=1.5, main="H=.6, d-decay function 1", pch=".")
abline(lm(r1$opClustering1~r1$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r2$expOutgroup, r2$opClustering1, xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 2", pch=".")
abline(lm(r2$opClustering1~r2$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r3$expOutgroup, r3$opClustering1, xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 3", pch=".")
abline(lm(r3$opClustering1~r3$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r4$expOutgroup, r4$opClustering1, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 1", pch=".")
abline(lm(r4$opClustering1~r4$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r5$expOutgroup, r5$opClustering1, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 2", pch=".")
abline(lm(r5$opClustering1~r5$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r6$expOutgroup, r6$opClustering1, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 3", pch=".")
abline(lm(r6$opClustering1~r6$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
# proportion non-western
plot(r1$p_nwa, r1$opClustering1, xlab="", ylab="op. clustering 1", cex.lab=1.5, main="", pch=".")
abline(lm(r1$opClustering1~r1$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r2$p_nwa, r2$opClustering1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$opClustering1~r2$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r3$p_nwa, r3$opClustering1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$opClustering1~r3$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r4$p_nwa, r4$opClustering1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$opClustering1~r4$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r5$p_nwa, r5$opClustering1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$opClustering1~r5$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r6$p_nwa, r6$opClustering1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$opClustering1~r6$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 1)
plot(r1$grSegregation1, r1$opClustering1, xlab="", ylab="op. clustering 1", cex.lab=1.5, main="", pch=".")
abline(lm(r1$opClustering1~r1$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r2$grSegregation1, r2$opClustering1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$opClustering1~r2$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r3$grSegregation1, r3$opClustering1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$opClustering1~r3$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r4$grSegregation1, r4$opClustering1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$opClustering1~r4$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r5$grSegregation1, r5$opClustering1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$opClustering1~r5$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r6$grSegregation1, r6$opClustering1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$opClustering1~r6$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 2)
plot(r1$grSegregation2, r1$opClustering2, xlab="", ylab="op. clustering 2", cex.lab=1.5, main="", pch=".")
abline(lm(r1$opClustering2~r1$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r2$grSegregation2, r2$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$opClustering2~r2$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r3$grSegregation2, r3$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$opClustering2~r3$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r4$grSegregation2, r4$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$opClustering2~r4$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r5$grSegregation2, r5$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$opClustering2~r5$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r6$grSegregation2, r6$opClustering2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$opClustering2~r6$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 3)
plot(r1$grSegregation3, r1$opClustering3, xlab="", ylab="op. clustering 3", cex.lab=1.5, main="", pch=".")
abline(lm(r1$opClustering3~r1$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r2$grSegregation3, r2$opClustering3, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$opClustering3~r2$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r3$grSegregation3, r3$opClustering3, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$opClustering3~r3$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r4$grSegregation3, r4$opClustering3, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$opClustering3~r4$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r5$grSegregation3, r5$opClustering3, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$opClustering3~r5$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r6$grSegregation3, r6$opClustering3, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$opClustering3~r6$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
dev.off()


# 3) Segregation -> opinion alignment ---------
png(
  filename = "./outputGraphics/macroResults 3 - segr-ali.png",
  width = 1000,
  height = 1000,
  units = "px"
)
par(mfrow=c(5,6))#,oma=c(0,6,0,0))
# exposure to outgroup
plot(r1$expOutgroup, r1$opAlignment1, xlab="", ylab="op. alignment 1", cex.lab=1.5, main="H=.6, d-decay function 1", pch=".")
abline(lm(r1$opAlignment1~r1$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r2$expOutgroup, r2$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 2", pch=".")
abline(lm(r2$opAlignment1~r2$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r3$expOutgroup, r3$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="H=.6, d-decay function 3", pch=".")
abline(lm(r3$opAlignment1~r3$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r4$expOutgroup, r4$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 1", pch=".")
abline(lm(r4$opAlignment1~r4$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r5$expOutgroup, r5$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 2", pch=".")
abline(lm(r5$opAlignment1~r5$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
plot(r6$expOutgroup, r6$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="H=.9, d-decay function 3", pch=".")
abline(lm(r6$opAlignment1~r6$expOutgroup), col="blue")
mtext("exposure to outgroup", side=1, line=2, cex=1)
# proportion non-western
plot(r1$p_nwa, r1$opAlignment1, xlab="", ylab="op. alignment 1", cex.lab=1.5, main="", pch=".")
abline(lm(r1$opAlignment1~r1$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r2$p_nwa, r2$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$opAlignment1~r2$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r3$p_nwa, r3$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$opAlignment1~r3$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r4$p_nwa, r4$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$opAlignment1~r4$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r5$p_nwa, r5$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$opAlignment1~r5$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
plot(r6$p_nwa, r6$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$opAlignment1~r6$p_nwa), col="blue")
mtext("prop. non-western", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 1)
plot(r1$grSegregation1, r1$opAlignment1, xlab="", ylab="op. alignment 1", cex.lab=1.5, main="", pch=".")
abline(lm(r1$opAlignment1~r1$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r2$grSegregation1, r2$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$opAlignment1~r2$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r3$grSegregation1, r3$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$opAlignment1~r3$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r4$grSegregation1, r4$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$opAlignment1~r4$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r5$grSegregation1, r5$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$opAlignment1~r5$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
plot(r6$grSegregation1, r6$opAlignment1, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$opAlignment1~r6$grSegregation1), col="blue")
mtext("segregation 1", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 2)
plot(r1$grSegregation2, r1$opAlignment2, xlab="", ylab="op. alignment 2", cex.lab=1.5, main="", pch=".")
abline(lm(r1$opAlignment2~r1$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r2$grSegregation2, r2$opAlignment2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$opAlignment2~r2$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r3$grSegregation2, r3$opAlignment2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$opAlignment2~r3$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r4$grSegregation2, r4$opAlignment2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$opAlignment2~r4$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r5$grSegregation2, r5$opAlignment2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$opAlignment2~r5$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
plot(r6$grSegregation2, r6$opAlignment2, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$opAlignment2~r6$grSegregation2), col="blue")
mtext("segregation 2", side=1, line=2, cex=1)
# group segregation (bivar LISA with distance decay function 3)
plot(r1$grSegregation3, r1$opAlignment3, xlab="", ylab="op. alignment 3", cex.lab=1.5, main="", pch=".")
abline(lm(r1$opAlignment3~r1$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r2$grSegregation3, r2$opAlignment3, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r2$opAlignment3~r2$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r3$grSegregation3, r3$opAlignment3, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r3$opAlignment3~r3$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r4$grSegregation3, r4$opAlignment3, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r4$opAlignment3~r4$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r5$grSegregation3, r5$opAlignment3, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r5$opAlignment3~r5$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
plot(r6$grSegregation3, r6$opAlignment3, xlab="", ylab="", cex.lab=1.5, main="", pch=".")
abline(lm(r6$opAlignment3~r6$grSegregation3), col="blue")
mtext("segregation 3", side=1, line=2, cex=1)
dev.off()








hist(abs(r$meanOpinionG1 - r$meanOpinionG2))
hist(r$polarizationIndex)

hist(r$opClustering1)
hist(r$opClustering2)
hist(r$opClustering3)

hist(r$opAlignment1)
hist(r$opAlignment2)
hist(r$opAlignment3)


# 1) Segregation -> opinion polarization
#
# Global level: 
plot(r$expOutgroup, r$polarizationIndex, xlab="Exposure to outgroup", ylab="Polarization index")
plot(r$NATexpOutgroup, r$polarizationIndex)
plot(r$NWAexpOutgroup, r$polarizationIndex)

plot(r$p_nwa, r$polarizationIndex)
plot(r$grClustering1, r$polarizationIndex)
plot(r$grClustering2, r$polarizationIndex)
plot(r$grClustering3, r$polarizationIndex)
plot(r$grSegregation1, r$polarizationIndex)
plot(r$grSegregation2, r$polarizationIndex)
plot(r$grSegregation3, r$polarizationIndex)

# Local level: 
#plot(rw$exposureOutgroup, mean(abs(rw$opinion)))
hist(
  rw$opinion,
  breaks = 20,
  main = NULL,
  xlim = c(-1,1),
  xlab = "Opinion distribution",
  yaxt = "n",
  ylab = "Frequency"
)

plotHeatMap( # default function. Saves to file.
  data=rw,
  x=rw$exposureOutgroup,
  y=rw$timeFirstPol,
  xlab="exposure to outgroup",
  ylab="timeFirstPol",
  fileName = "test"
)

heatmap( # Custom function. Prints on screen.
  x = rw$exposureOutgroup,
  y = rw$timeFirstPol,
  bins=10,
  type = "count",#"prop",
  xlab = "exposure to outgroup",
  ylab = "# interactions before becoming extremist"
)















plotHeatMap(
  data=rw,
  x=rw$exposureOutgroup,
  y=rw$opinion,
  xlab="exposure to outgroup",
  ylab="opinion",
  fileName = "micro_expOutgroup_opinion2"
)
plotHeatMap(
  data=rw,
  x=rw$exposureOutgroup,
  y=rw$timeFirstPol,
  xlab="exposure to outgroup",
  ylab="# interactions before becoming extremist",
  fileName = "micro_groupClustering_timeFirstPol"
)





ggplot(rw, aes(x=rw$exposureOutgroup, y=abs(rw$opinion))) +
  #geom_point(shape=20, alpha= 0.01)+
  geom_smooth(method='lm')



a <- rw$exposureOutgroup[1:1000]#runif(1000, min=0, max=1)
b <- rw$opinion[1:1000]#rnorm(1000, mean=0, sd=1)
#b <- (a*0.5)+b
t <- as.data.frame(cbind(a,b))

ggplot(t, aes(x=a, y=abs(b))) + geom_point(shape=20,alpha=0.1)+
  geom_smooth(method='lm')





# 2) Segregation -> opinion clustering
#
# Global level: 
plot(r$expOutgroup, r$opClustering1, xlab="exposure to outgroup", ylab="spatial clustering of opinions")
plot(r$NATexpOutgroup, r$opClustering1)
plot(r$NWAexpOutgroup, r$opClustering1)


plot(r$p_nwa, r$opAlignment1)
plot(r$grClustering1, r$opClustering1)
plot(r$grClustering2, r$opClustering2)
plot(r$grClustering3, r$opClustering3)
plot(r$grSegregation1, r$opClustering1)
plot(r$grSegregation2, r$opClustering2)
plot(r$grSegregation3, r$opClustering3)


# Local level: 
plotHeatMap(
  data=rw,
  x=rw$exposureOutgroup,
  y=rw$opClustering,
  bins=40,
  xlab="exposure to outgroup",
  ylab="opinion clustering",
  fileName = "micro_expOutgroup_opClustering"
)


# 3) Segregation -> alignment
#
# Global level: 
plot(r$expOutgroup, r$opAlignment, xlab="exposure to outgroup", ylab="group-opinion alignment")
plot(r$NATexpOutgroup, r$opAlignment)
plot(r$NWAexpOutgroup, r$opAlignment)


plot(r$p_nwa, r$opAlignment)
plot(r$grClustering1, r$opAlignment1)
plot(r$grClustering2, r$opAlignment2)
plot(r$grClustering3, r$opAlignment3)
plot(r$grSegregation1, r$opAlignment1)

# Local level: 
plotHeatMap(
  data=rw,
  x=rw$exposureOutgroup,
  y=rw$opAlignment,
  bins=40,
  xlab="exposure to outgroup",
  ylab="group-opinion alignment",
  fileName = "micro_expOutgroup_opAlignment"
)




test <- subset(r, r$H == .3)
summary(lm(r$opClustering ~ r$WK_CODE))
summary(lm(r$opClustering ~ r$globalI_biv_count + r$p_nwa))

plot(r$WK_CODE, r$opClustering)



hist(r$opClustering)
hist(abs(rw$opAlignment))

hist(r$polarizationIndex)



################################################################################
table(names(r) %in% names(rw))
names(r)

# Adding independent variables to the complete agentlist rw
temp  <- r[,-c()]
r <- base::merge(
  x = rw,
  y = temp,
  by = "seed"
)
rm(temp)

















plot(r$expOutgroup, r$polarizationIndex)
summary(lm(r$polarizationIndex ~ r$expOutgroup))

plot(r$expOutgroup, abs(r$meanOpinionG2))
summary(lm(r$polarizationIndex ~ abs(r$meanOpinionG2)))

plot(r$expOutgroup, abs(r$meanOpinionG1 - r$meanOpinionG2))
summary(lm(r$polarizationIndex ~ abs(r$meanOpinionG1 - r$meanOpinionG2)))



run(
  seed = 12345,
  indexParameters = 99999,
  wijk = 1,
  initialOpinionDistribution = "uniform",
  H = 0.3,
  distanceDecay = 2,
  timeMax = 10,
  printOpinionHistogram = TRUE,
  exportOutput = FALSE
)
hist(test[[2]]$opinion)



# ==============================================================================
# Old implementation (do not run anything that follows)

# Time to save everything into a single file, that we shall name
# "completeDataset". This may well take a long while (and a lot of memory)
parameterSpace <- GparameterSpace
parameterSpace$printOpinionHistogram <- parameterSpace$exportOutput <- NULL
simResults <- GsimResults
simW <- GsimW
save(
  parameterSpace, simResults, simW,
  file = "./simOutput/completeDataset.RDATA"
  #file = "./simOutput/completeBaselineDataset.RDATA"
)

# Loading results, battery by battery. We combine all the results into
# global objects (GsimResults and GsimW)
GparameterSpace <- data.frame()
GsimResults <- data.frame()
GsimW <- list()
for (i in 1:length(files)){
  print(paste("Loading file", i, "of", length(files)))
  file <- files[i]
  load(paste0("./simOutput/peregrine/", files[i]))
  if (i == 1){
    GparameterSpace <- parameterSpace
  } else {
    if (any(GparameterSpace != parameterSpace)){
      stop(print("Can't load results batteries with inconsistent parameterSpace."))
    }
  }
  #simResults <- simResults[49:60,] ###################
  #simW <- simW[49:60]              ###################
  GsimResults <- as.data.frame(rbind(GsimResults,simResults))
  #GsimW <- c(GsimW, simW)# this is just too massive : it's 5GB when compressed..
}

# Time to save everything into a single file, that we shall name
# "completeDataset". This may well take a long while (and a lot of memory)
parameterSpace <- GparameterSpace
parameterSpace$printOpinionHistogram <- parameterSpace$exportOutput <- NULL
simResults <- GsimResults
simW <- GsimW
save(
  parameterSpace, simResults, simW,
  file = "./simOutput/completeDataset.RDATA"
  #file = "./simOutput/completeBaselineDataset.RDATA"
)
rm(files, GparameterSpace, GsimResults, GsimW)


#
# 
# Here we prepare the datasets for the data analyses.

# We import the settings and the output of the battery of simulation runs.
rm (list = ls( )) 
load("./cityData/geodata_Rotterdam.RData")
rm(worldList)
#load("./simOutput/sims_1.RData") #############################################
load("./simOutput/completeBaselineDataset.RDATA")

# We join the dataframes with parameter settings, simulation results, and 
# district descriptives.
parameterSpace$indexParameters <- c(1:nrow(parameterSpace))
r <- base::merge(
  x = simResults,
  y = parameterSpace,
  by = "indexParameters"
)
r$order <- c(1:nrow(r)) # Because joining dataframes can shuffle the ordering of
# the simulation runs, we create an index variable that
# we can later use to sort the dataframe in the correct
# order
citySummary$WK_NR <- c(1:nrow(citySummary))
r <- base::merge(
  x = r,
  y = citySummary,
  by.x = "wijk",
  by.y = "WK_NR"
)
r <- r[order(r$order),] # We sort the dataframe, and rid of the ordering index.
#r <- r[,-1]


# Adding measures of opinion clustering and alignment.
#
# We retrieve the opinion vector from the simulated district, we use it to
# calculate uni- and bivariate, global and local I.
# Then we merge the information to the results dataframes, r and rw:
#
#  > r: each row is a summary of a simulation run - contains global variables.
#  > rw: each row is an agent. Contains local variables. Variable rw$simRun 
#          nests agents by the simulation run they were in.
r$opClustering <- r$opAlignment <- NA
#rw <- as.data.frame(NA)
for (i in 1:nrow(r)){
  w <- simW[[i]]
  w$simRun <- i
  proxmat <- NA
  if (r$distanceDecay[i] == 1) {proxmat <- proximityList1[[r$wijk[i]]]}
  if (r$distanceDecay[i] == 2) {proxmat <- proximityList2[[r$wijk[i]]]}
  if (r$distanceDecay[i] == 3) {proxmat <- proximityList4[[r$wijk[i]]]}
  dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[r$wijk[i]])
  dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
  ops <- c(NA, length = length(dat))
  for (l in 1:length(dat)){
    cell <- subset(w, w$location == dat$OBJECTID[l])
    ops[l] <- mean(cell$opinion) ################  <  missing local measures here
  }
  dat$opinion <- ops
  print(paste0(
    "Run ", i, " of ", nrow(r), ". Calculating Moran's I on various attributes."
  ))
  
  # Global opinion clustering
  r$opClustering[i] <- MoranI(
    x = dat$opinion,
    proxmat = proxmat
  )
  
  # Global opinion alignment
  r$opAlignment[[i]] <- MoranI(
    x = dat$pnwal2014,
    y = dat$opinion,
    proxmat = proxmat
  )
  
  if(FALSE){ #########################
    # Local opinion clustering
    Ic <- MoranI(
      x=dat$opinion,
      proxmat=proxmat,
      type = "local"
    )
    
    # Local opinion alignment
    Ia <- MoranI(
      x=dat$pnwal2014,
      y=dat$opinion,
      proxmat=proxmat,
      type = "local"
    )
    ################################### presi da qui e messi in calib: Is, Iec
    
    for (a in 1:nrow(w)){
      w$opClustering[a] <- Ic[w$index[a]]
      w$opAlignment[a] <- Ia[w$index[a]]
      #w$localI_biv_count[a] <- Is[w$index[a]]
      #w$localI_count[a] <- Iec[w$index[a]]
      #w$p_nwa[a] <- dat$pnwal2014[w$index[a]]
    }
    
    #if(i<6){print(paste(r$opClustering[i], r$opAlignment[i], Ic[w$index[i]],Ia[w$index[i]]))}#########
    rm(Ic, Ia)
    ifelse(
      i != 1,
      rw <<- as.data.frame(rbind(rw, w)),
      rw <<- w
    )
  } ###################
}

# Saving the results datasets.
save(
  r,
  rw,
  file = "./simOutput/analyses.RDATA"
)


# Built-in function to plot heatmaps
#
plotHeatMap <- function(
  data,
  x,
  y,
  bins=20,
  xlab,
  ylab,
  fileName,
  legendTitle = "density"
){
  png(paste0("outputGraphics/", fileName,  ".png"), 
      width = 800, height = 600,
      units = "px", pointsize = 10, 
      res = NA)
  print(ggplot(
    data,
    stat = "density",
    aes(
      x=x,
      y=y
    )
  ) + #xlim(c(0.01, 1)) +
    scale_x_continuous(limits = c(NA,NA), expand = c(0, 0)) +
    scale_y_continuous(limits = c(NA,NA), expand = c(0, 0)) +
    geom_bin2d(bins = bins) +#, color ="white")+
    labs(
      fill = legendTitle,
      x = xlab,
      y = ylab
    ) +
    scale_fill_gradient(low =  "white", high = "black") +
    theme(
      axis.line=element_blank(),
      axis.ticks.length = unit(.25, "cm"),
      #axis.ticks = element_line(),
      axis.text = element_text(size=20), ####
      axis.title.y = element_text(size=25),
      axis.title.x = element_text(size = 25),
      title = element_text(size=25),
      legend.text = element_text(size=15),
      panel.background=element_blank(), plot.background=element_blank(),
      panel.border = element_blank()
    ))
  dev.off()
}



load("./testProximityMatrices.RDATA")

distance <- function1 <- function2 <- function3 <- c()
index <- 1
for (x in c(100, 200, 300, 500, 1000, 2000, 5000)){
  distance[index] <- x
  function1[index] <- sum(distances1 <= x) /length(distances1)
  function2[index] <- sum(distances2 <= x)/length(distances2)
  function3[index] <- sum(distances3 <= x)/length(distances3)
  
  index <- index + 1
}
test <- as.data.frame(cbind(distance,function1,function2,function3))
test
rm(test)


par(mfrow=c(1,3))
hist(
  distances1,
  col="red",
  xlim = c(0,2000),
  breaks=c(0:1000)*50,
  main="Distance decay function 1",
  yaxt = "n",
  ylab=NULL,
  xlab="distance between i and j (in meters)"
)
mtext(paste0("Max distance: ", round(max(distances1)), "m"), side=3)
hist(
  distances2,
  col="blue",
  xlim = c(0,2000),
  breaks=c(0:1000)*50,
  main="Distance decay function 2",
  yaxt = "n",
  ylab=NULL,
  xlab="distance between i and j (in meters)"
)
mtext(paste0("Max distance: ", round(max(distances2)), "m"), side=3)
hist(
  distances3,
  col="green",
  xlim = c(0,2000),
  breaks=c(0:1000)*50,
  main="Distance decay function 3",
  yaxt = "n",
  ylab=NULL,
  xlab="distance between i and j (in meters)"
)
mtext(paste0("Max distance: ", round(max(distances3)), "m"), side=3)



# measuring baseline alignment at t=0
d <- ali1 <- ali2 <- ali3 <- aliCon <- aliExtr <- c()
for (r in 0:9){
  print(paste0("replication #", r + 1 ," of 10"))
  for(i in 1:12){
    d[12*r + i] <- i
    w <- worldList[[i]]
    w$opinion <- rbeta(nrow(w), 3, 3, ncp = 0)
    w$opinion <- w$opinion * 2 - 1
    dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[i])
    dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
    ops <- c(NA, length = length(dat))
    for (l in 1:length(dat)){ # for every cell
      cell <- subset(w, w$location == dat$OBJECTID[l])
      ops[l] <- mean(cell$opinion)
    }
    if(FALSE){ ######################
      dat$opinion <- ops
      ali1[12*r + i] <- mean(abs(MoranI( # global alignment
        x = dat$pnwal2014,
        y = dat$opinion,
        proxmat = proximityList1[[i]],
        type = "local"
      )))
      ali2[12*r + i] <- mean(abs(MoranI( # global alignment
        x = dat$pnwal2014,
        y = dat$opinion,
        proxmat = proximityList2[[i]],
        type = "local"
      )))
      ali3[12*r + i] <- mean(abs(MoranI( # global alignment
        x = dat$pnwal2014,
        y = dat$opinion,
        proxmat = proximityList3[[i]],
        type = "local"
      )))
    }##########################
    
    dat$opinion <- 0
    
    aliCon[12*r + i] <- mean(abs(MoranI( # global alignment
      x = dat$pnwal2014,
      y = dat$opinion,
      proxmat = proximityList1[[i]],
      type = "local"
    )))
    dat$opinion <- 1
    aliExtr[12*r + i] <- mean(abs(MoranI( # global alignment
      x = dat$pnwal2014,
      y = dat$opinion,
      proxmat = proximityList1[[i]],
      type = "local"
    )))
  }
}





# Counting interactin events between random agents.
rm (list = ls( )) 
#dist0 <- dist100 <- dist500 <- dist1000 <- c()
# Loading resources
source("script NI&PA.R")
source("util.R")
require(spdep)
require(geosphere)
load("./cityData/geodata_Rotterdam.RData")
run (
  timeMax = 0,
  wijk = 1,
  exportOutput = FALSE
)
dat <- cbs100_rot[cbs100_rot$WK_CODE==districtsList[1],
                  c("nauto2014", "nnwal2014")]
dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
ncells <- length(dat)
distmat <- matrix(NA, nrow=ncells, ncol=ncells)
for (i in 1:ncells) {
  distmat[i,] <- distVincentyEllipsoid(coordinates(dat)[i,], coordinates(dat)) 
}

proxmat <- proximityList2[[1]]
egos <- alters <- distances <- c()
nCells <- length(table(agents$index))
popDensity <- c()
for (i in 1:nCells) {
  popDensity[i] <- nrow(subset(agents, agents$index == i))
}
probmat <- matrix(NA, nrow=nCells, ncol=nCells)
for (i in 1:nCells){
  probmat[,i] <- proxmat[,i] * popDensity[i]
}

for (ego in sample(nrow(agents), size=1)) {
  j <- sample(nrow(agents), size=1000)
  print(paste("counting encounters between", ego, "and 1000 potential alters"))
  for (iterations in 1:500000){
    egoCell <- agents$index[ego]
    targetCell <- sample(c(1:nCells), 1, prob=probmat[egoCell,])
    distance <- distmat[egoCell, targetCell]
    #distances1 <<- append(distances1, distmat[egoCell, targetCell]
    repeat{
      alter <- sample(c(1:popDensity[targetCell]),1)
      alter <- as.numeric(rownames(agents[which(agents$index==targetCell),][alter,]))
      if (ego != alter) {break}
    }
    if (alter %in% j){
      #print("found one")
      egos <- append(egos, ego)
      alters <- append(alters, alter)
      distances <- append(distances, distance)
    }
  }
  #print(paste("i=",ego,"from",egoCell,"j=",alter,"from",targetCell,"distance", distance))
}
cont3 <- as.data.frame(cbind(egos,alters,distances))

dist <- tally <- c()
for (i in unique(cont3$alters)){
  x <- subset(cont3, cont3$alters==i)
  dist <- append(dist, x$distances[1])
  tally <- append(tally, nrow(x))
  print(paste(i, "dist:", x$distances[1], "tally", nrow(x)))
}
test <- as.data.frame(cbind(unique(cont3$alters), dist, tally))
#print(test)
print(head(test[order(test$dist),]))
plot(dist, tally, xlim=c(0,1500))



#

rm (list = ls( )) 

source("script NI&PA.R")
test<-run (
  timeMax = 10,
  wijk = 9,
  exportOutput = TRUE
)






#rm (list = ls( )) 
#load("./simOutput/completeDataset.RDATA")



#==============================================================================
rm (list = ls( )) 


source("script NI&PA.R")
source("util.R")
run(timeMax = 0, wijk = 1)
i <- 1

dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[i])
dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
dat$observed1 <- dat$nnwal2014 - mean(dat$nnwal2014)
dat$grClustering1 <- MoranI(
  x=dat$nnwal2014,
  proxmat=proximityList1[[i]],
  type = "local"
)
hist(dat$grClustering1)
plot(dat$observed1, dat$grClustering1)

for (l in 1:nrow(dat)){
  ifelse(
    dat$j[l] < mean(dat$j),
    dat$test[l] <- 0,
    dat$test[l] <- 1
  )
}
dat$observed2 <- dat$test - mean(dat$test)
dat$testClu <- MoranI(
  x=dat$test,
  proxmat=proximityList1[[i]],
  type = "local"
)
hist(dat$testClu)
plot(dat$observed2,dat$testClu)


#============================================================================
# making sure MoranI() works as intended
run(timeMax = 0, wijk = 1)
i <- 1
dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[i])
dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
dat$gr <- dat$nnwal2014 - mean(dat$nnwal2014)
ops <- c(NA, length = length(dat))
for (l in 1:length(dat)){ # for every cell
  cell <- subset(agents, agents$location == dat$OBJECTID[l])
  ops[l] <- mean(cell$opinion)
}
dat$opinion <- ops

dat$segr <- MoranI(
  x=dat$gr,
  y=dat$opinion,
  proxmat=proximityList1[[i]],
  type = "local"
)
plot(dat$opinion, dat$segr)



for(a in 1:nrow(agents)){
  ifelse(
    agents$group[a] == 1,
    agents$opinion[a] <- 0.7,
    agents$opinion[a] <- -0.7
  )
}
ops <- c(NA, length = length(dat))
for (l in 1:length(dat)){ # for every cell
  cell <- subset(agents, agents$location == dat$OBJECTID[l])
  ops[l] <- mean(cell$opinion)
}
dat$opinion <- ops

dat$segr <- MoranI(
  x=dat$gr,
  y=dat$opinion,
  proxmat=proximityList1[[i]],
  type = "local"
)
plot(dat$opinion, dat$segr)




# Extra exploration: do opinion distributions of the two groups differ for 
# low or high values of outgroup exposure?
if(FALSE){ #######################
  f=rr[rr$polarizationIndex<0.5 &
         rr$polarizationIndex>0.3 & rr$wijk==2,][4,]$fileName
  s=rr[rr$polarizationIndex<0.5 &
         rr$polarizationIndex>0.3 & rr$wijk==2,][4,]$seed
  load(paste0("./simOutput/peregrine/", f))
  temp <- simResults[simResults$seed==s,]$indexParameters
  s<-simW[[temp]]
  
  par(mfrow=c(2,2))
  hist(
    subset(s, s$group==-1 & s$expOutgr2<0.4)$opinion,
    breaks=10,xlim=c(-1,1),xlab="attitude", yaxt="n",ylab="",
    yaxt='n',col="black", border="white",
    main="outgroup exposure in [0,0.4]\n\nwestern (majority)"
  )
  title(ylab="relative frequency", line=0)
  hist(
    subset(s, s$group==1 & s$expOutgr2<0.4)$opinion,
    breaks=10,xlim=c(-1,1),xlab="attitude", yaxt="n",ylab="",yaxt='n',
    col="black", border="white",
    main="\n\n\nnon-western (minority)"
  )
  title(ylab="", line=0)
  hist(
    subset(s, s$group==-1 & s$expOutgr2>0.4 & s$expOutgr2<=0.6)$opinion,
    breaks=10,xlim=c(-1,1),xlab="attitude", yaxt="n",ylab="",yaxt='n',
    col="black", border="white",
    main="outgroup exposure in (0.4,0.6]\n\nwestern (majority)"
  )
  title(ylab="relative frequency", line=0)
  hist(
    subset(s, s$group==1 & s$expOutgr2>0.4 & s$expOutgr2<=0.6)$opinion,
    breaks=10,xlim=c(-1,1),xlab="attitude", yaxt="n",ylab="",yaxt='n',
    col="black", border="white",
    main="\n\n\nnon-western (minority)"
  )
  title(ylab="", line=0)
}####################################



# __________________________
# Polarization index vs SD
polarizationIndex <- function(x){
  diff <- c()
  for (i in 1:length(x)){
    for (j in 1:length(x)) {
      if(i != j) {
        oD <- abs(x[i] - x[j])
        diff <- append (diff, oD)
      }
    }
  }
  return(var(diff))
}


x <- rbeta(n = 100, 0.3, 0.3) / 10
#x <- rnorm(n = 100, mean = 0, sd = 0.01)
hist(x, xlim = c(-1,1))
print(paste0("sd = ", sd(x), "; pol = ", polarizationIndex(x)))





#_______________________________________________________________________________
# Plotting extreme local alignment scores.
rm(list=ls())
source("util.R")
load("./cityData/geodata_Rotterdam.RData")
load("./simOutput/completeDataset.RDATA")
r$opAlignment1 <- abs(r$opAlignment1)
r$opAlignment2 <- abs(r$opAlignment2)
r$opAlignment3 <- abs(r$opAlignment3)
r$iniAli1 <- abs(r$iniAli1)
r$iniAli2 <- abs(r$iniAli2)
r$iniAli3 <- abs(r$iniAli3)

ri$opAlignment1 <- abs(ri$opAlignment1)
ri$opAlignment2 <- abs(ri$opAlignment2)
ri$opAlignment3 <- abs(ri$opAlignment3)
ri$iniOpAlignment1 <- abs(ri$iniOpAlignment1)
ri$iniOpAlignment2 <- abs(ri$iniOpAlignment2)
ri$iniOpAlignment3 <- abs(ri$iniOpAlignment3)

# selecting the runs from the baseline parameter configuration.
rr <- subset(
  r,
  r$initialOpinionDistribution == "groupBias" & 
    r$H == 0.6 &
    r$distanceDecay == 2
)
rri <- subset(
  ri,
  ri$initialOpinionDistribution == "groupBias" & 
    ri$H == 0.6 &
    ri$distanceDecay == 2
)

runs <- c()
for (w in 1:12) {runs[w] <- unique(rri[rri$wijk == w,]$seed)[1]}
rri2 <- rri[rri$seed %in% runs,]

seedz <- unique(rri2$seed)
rr2 <- r[r$seed %in% seedz,]

param <- expand.grid(
  wijk = 1:12,
  distr = unique(rri2$initialOpinionDistribution)
)

maxAli <- list(
  param = param,
  maxAli2 = list()
)


for (p in 1:nrow(param)){
  print(paste0("parameter config. ", p, " of ", nrow(param),". ", Sys.time()))
  wijk <- param$wijk[p]
  distr <- as.character(param$distr[p])
  agents <- worldList[[wijk]]
  
  scores <-  rep(NA, times = citySummary$n_squ[wijk])
  
  G1 <- which(agents$group == 1, arr.ind = TRUE) ###############################
  G2 <- which(agents$group != 1, arr.ind = TRUE)
  
  agents$opinion[G1] <- 1
  agents$opinion[G2] <- -1
  
  dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[wijk])
  dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
  dat$dens <- dat$inw2014 / citySummary$n_pop[wijk]
  ops <- c(NA, length = length(dat))
  for (l in 1:length(dat)){ # for every cell
    cell <- subset(agents, agents$location == dat$OBJECTID[l])
    ops[l] <- mean(cell$opinion)
  }
  dat$opinion <- ops
  
  scores <- moranI( # opinion-group alignment
    x = dat$pnwal2014, y = dat$opinion, proxmat =proximityList2[[wijk]],
    dens = dat$dens, N = citySummary$n_pop[wijk]
  )$localI

  
  maxAli$maxAli2[[p]] <- scores
}

#save(iniAli, file = "./simOutput/iniAli.RDATA")
rri2$opAlignment2 <- NA

for (p in 1:nrow(maxAli$param)){
  print(paste0("Processing ", p, " of ", nrow(maxAli$param), ". ", Sys.time()))
  
  for(index in 1:length(maxAli$maxAli2[[p]])) {
    a <- which(
      rri2$wijk == maxAli$param$wijk[p] & 
        rri2$initialOpinionDistribution == as.character(maxAli$param$distr[p]) &
        rri2$index == index
    )
    rri2$iniOpAlignment2[a] <- maxAli$maxAli2[[p]][index]
  }
  
}


rri2 <- rri[rri$wijk %in% c(1, 4, 7),]

districtLabels <- c(
  "1" = "Stadscentrum",
  "2" = "Delfshaven",
  "3" = "Overschie",
  "4" = "Noord",
  "5" = "Hillegersberg-Schiebroek",
  "6" = "Kralingen-Crooswijk",
  "7" = "Feijenoord",
  "8" = "IJsselmonde",
  "9" = "Pernis",
  "10" = "Prins Alexander",
  "11" = "Charlois",
  "12" = "Hoogvliet"
)

ggplot(
  rri2,
  aes(
    cut(expOutgr2, breaks=(0:5 / 5)),
    abs(opAlignment2)
  )) +
  #ggtitle("3rd run per district") +
  ggtitle("opinion west = -1\nopinion non-west = 1") +
  ylab("local alignment (s=100)") +
  xlab("outgroup exposure (s=100)") +
  geom_violin( # Alignment at t=0 (white)
    aes(y = abs(iniOpAlignment2)),
    fill = "white", color = "#ababab", scale = "width", draw_quantiles = 0.5) +
  geom_violin( # Alignment at t = 200 (gray)
    fill = "gray",
    scale = "width",
    draw_quantiles = 0.5
  ) +
  facet_grid(
    rri2$wijk ~ (rri2$group - 1),
    labeller = as_labeller(
      c(districtLabels, c("0"="non-western", "-2" = "western")))) +
  theme(
    plot.margin=unit(c(0,0,0,0),"pt"),
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(color=NA,fill="gray97"),
    axis.line = element_line(colour = "black"),
    strip.text.y = element_text(angle = 0, hjust=0),
    strip.background.y = element_rect(colour=NA, fill=NA)
  )



#_______________________________________________________________________________
# Checking baseline results for Pernis (why did it not reach the 
# extreme polarization & extreme alignment like the other districts?)

load("./simOutput/completeDataset.RDATA") 
library(ggplot2)

# Selecting the baseline runs from Pernis:
rr <- subset(
  r,
  r$initialOpinionDistribution == "groupBias" & 
  #r$wijk == 9 &
    r$H == 0.6 &
    r$distanceDecay == 2
)
rri <- subset(
  ri,
  ri$initialOpinionDistribution == "groupBias" & 
  #r$wijk == 9 &
    ri$H == 0.6 &
    ri$distanceDecay == 2
)


rri2 <- rri[rri$wijk %in% c(9, 2, 1),]

ggplot(
  rri2,
  aes(
    cut(expOutgr2, breaks=(0:5 / 5)),
    opAlignment2
  )) +
  #ggtitle("3rd run per district") +
  #ggtitle("no group bias") +
  ylab("local alignment (s=100)") +
  xlab("outgroup exposure (s=100)") +
  geom_violin( # Alignment at t = 200 (orange)
    color = "darkorange", fill = "darkorange", alpha = 0.4,#fill="gray",
    scale = "width",
    draw_quantiles = 0.5,
    position = position_nudge(x = 0.1)
  ) +
  facet_grid(
    rri2$wijk ~ (rri2$group - 1),
    labeller = as_labeller(
      c(districtLabels, c("0"="non-western", "-2" = "western")))) +
  theme(
    plot.margin=unit(c(0,0,0,0),"pt"),
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    #panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(color=NA,fill="gray97"),
    axis.line = element_line(colour = "black"),
    strip.text.y = element_text(angle = 0, hjust=0),
    strip.background.y = element_rect(colour=NA, fill=NA)
  )



rri2 <- rri[rri$wijk == 9,]
table(round(rri2$opinion), rri2$group)
table(rri2$opinion > 0.999 | rri2$opinion < 0.999)


#actual:
sd(c(rep(1, times = 46251), rep(-1, times = 3749)))

#theoretical max:
sd(c(rep(1, times = 46000), rep(-1, times = 4000)))

# Checking actual runs with positive alignment score (in countertendency)
seednr = 15
seedid <- unique(rri2[rri2$opAlignment2 > 0,]$seed)[seednr]
temp <- r[r$seed == seedid,]
load(paste0("./simOutput/peregrine/",temp$fileName))
w <- simW[[temp$indexParameters]]

table(w$opinion > 0.999 | w$opinion < 0.999)
table(round(w$opinion), w$group)

plot(w[w$group == 1,]$expOutgr2, w[w$group == 1,]$opAlignment2)
abline(h = 0, col = "red")


#plot(rri2$expOutgr2, rri2$opAlignment2)
x <- rri2[rri2$group == 1 & rri2$seed == seedid,]
plot(x$expOutgr2, x$opAlignment2)
abline(h = 0, col = "red")

x <- expand.grid(
  x_coor = unique(rri2$x_coor),
  y_coor = unique(rri2$y_coor)
)
x$g <- x$o <- x$a <- NA
for(i in 1:nrow(x)){
  xi <- w[w$x_coor == x$x_coor[i] & w$y_coor == x$y_coor[i],]
  if(nrow(xi) > 0){
    x$g[i] <- mean(xi$group)
    x$o[i] <- mean(xi$opinion)
    x$a[i] <- sum(xi$opAlignment2 > 0)
  }
}
x <- x[!is.na(x$g),]
x$a[x$a == 0] <- NA

ggplot(data = x, aes(x = x_coor, y_coor)) +
  geom_point(aes(fill = g), shape = 22, size = 18) +
  geom_point(aes(color = a), shape = 16, size = 5) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA)



aa <- c()
for (wijk in 1:12) {
  print(paste(
    citySummary$district[wijk],
    round(
      sum(worldList[[wijk]]$expOutgr2 > 0.7) / nrow(worldList[[wijk]]) * 100,
      digits = 2
    )
  ))
  aa <- c(aa, worldList[[wijk]]$expOutgr2)
}

length(aa)
round(sum(aa > 0.7) / length(aa) * 100, digits = 2)




###################33
#addTiles(leaflet())
spanH = 0.00071
spanV = 0.00045
#spanH = 0.00071
library("leaflet")
library("mapview")


m <- leaflet() %>%
  addTiles() %>%  # Add default OpenStreetMap map tiles
  #addMarkers(lng=4.4777325, lat=51.9244201, popup="Rotterdam") %>%
  addProviderTiles(providers$Stamen.TonerBackground) %>%
  #fitBounds(4.375, 51.885, 4.4, 51.894)  %>% # nice boundaries
  setView(
    lng = mean(x$y_coor),
    lat = mean(x$x_coor),
    zoom = 16
  ) %>%
  addRectangles(
    lng1 = x$y_coor - spanH, lng2 = x$y_coor + spanH,
    lat1 = x$x_coor - spanV, lat2 = x$x_coor + spanV,
    color = "transparent", fillColor = "darkorange", opacity = x$g#1#x$g#"transparent"
  )
  #addCircles(
  #  lng = x$y_coor,
  #  lat = x$x_coor,
  #  fillColor = c('red'))
m

mapshot(m, file = "./outputGraphics/map.png")

colorRampPalette(c("white", "darkOrange"), x$g)

################




#get_map(
#  location=c(mean(x$x_coor), mean(x$y_coor)),
#        source="stamen", maptype="watercolor", crop=FALSE)
#ggmap(myMap)
#ggmap(myMap)


# Starting opinion of minority agents who misalign

source("script NI&PA.r")

run(
  timeMax = 0,
  resetWorld = TRUE,
  seed  =  158749486,
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

a <- agents$opinion

run(
  timeMax = 200,
  resetWorld = TRUE,
  seed  =  158749486,
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

b <- agents$opinion
hist(b)


simdata <- data.frame(
  iniop = a,
  finop = b,
  group = factor(
    agents$group, levels = c(-1, 1),
    labels = c("western\nmajority", "non-western\nminority")),
  time = agents$nIntFirstExtr,
  location = agents$location
)
for(i in 1:nrow(simdata)){
  if (simdata$finop[i] <= -0.999) {simdata$finop[i] <- -1}
  if (simdata$finop[i] >= -0.999) {simdata$finop[i] <- 1}
}

library(ggplot2)

ggplot(data = simdata, aes(x = iniop, y = as.factor(finop))) +
  geom_boxplot() + 
  facet_wrap(~ group) +
  xlab("initial opinion\nt=0") + ylab("final opinion\nt=200")

ggplot(data = simdata, aes(x = time, y = as.factor(finop))) +
  geom_boxplot() + 
  facet_wrap(~ group) +
  xlab("# interactions before first extremization") +
  ylab("final opinion\nt=200")

table(simdata$group, simdata$finop)



loc <- data.frame(
  loc = unique(simdata[simdata$group == "non-western\nminority",]$location)
)
for (l in 1:nrow(loc)){
  x <- simdata[
    simdata$group == "non-western\nminority" & simdata$location == loc$loc[l],]
  loc$nNonWest[l] <- nrow(x)
  loc$negativeop[l] <- sum(x$finop == -1)
  loc$positiveop[l] <- sum(x$finop == 1)
}
names(loc) <- c(
  "sq.unit ID",
  "# non-west.",
  "# non-misaligned",
  "# misaligned"
)
write.csv(loc[loc$positiveop >= 1,], file = "./other/misaligned.csv")
print(loc[loc$"# misaligned" >= 1,])

#load("./cityData/geodata_Rotterdam.RData")
#wijk <- 9
#dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[wijk])
#dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
#dat$dens <- dat$inw2014 / citySummary$n_pop[wijk]
#ops <- c(NA, length = length(dat))
#for (l in 1:length(dat)){ # for every cell
#  cell <- subset(agents, agents$location == dat$OBJECTID[l])
#  ops[l] <- mean(cell$opinion)
#}
#dat$opinion <- ops
#
#opAlignment2 <- moranI( # opinion-group alignment
#  x = dat$pnwal2014,
#  y = dat$opinion,
#  proxmat = proximityList2[[wijk]],
#  dens = dat$dens,
#  N = citySummary$n_pop[wijk]
#)
#
#as <- base::merge(
#  x = agents[,1:2],
#  y = as.data.frame(cbind(
#    dat$OBJECTID,
#    opAlignment2 = opAlignment2$localI,
#  )),
#  by.x="location",
#  by.y = "V1"
#)
#W$opAlignment2 <- as$opAlignment2




ggplot(rr, aes(factor(expOutgr2), abs(intuitiveAlignment))) +
  geom_segment(aes(x = 0.5, xend = 12.5, y = 2, yend = 2), color = "black") +
  geom_violin( # alignment at t=0
    aes(y = iniIntuitiveAlignment),
    fill = "white", color = "#ababab",
    scale = "width",
    draw_quantiles = 0.5
  ) +
  geom_violin( # alignment at t=200
    color = "darkorange", fill = "darkorange", alpha = 0.4,#fill="gray",
    scale = "width",
    draw_quantiles = 0.5
  ) +
  ggtitle("districts ordered by\naverage outgroup exposure (s=100)") +
  ylab("between-groups attitude difference") +
  scale_x_discrete(labels = labs) +
  scale_y_continuous(
    #trans = scaleSquash(cut_from, cut_to, 50), ### squashing the y axis
    limits = c(0,2.1),
    expand = c(0,0)#,
    #breaks = breaks[breaks <= cut_from | breaks >= cut_to]
  ) +
  #geom_rect( # Highlighting the squashed portion of the chart
  #  aes(xmin = 0, xmax = 13, ymin = cut_from, ymax = cut_to),
  #  fill = "#fafafa"
  #) +
  #geom_segment( # highlighting the squashed portion on the Y axis
  #  aes(x = 0, xend = 0, y = 0, yend = cut_from), size = 0.5, color = "black") +
  #geom_segment( # highlighting the squashed portion on the Y axis
  #  aes(x = 0, xend = 0, y = cut_to, yend = 2), size = 0.5, color = "black") +
  #annotate("text", x = 0, y = cut_from, label = "\\", angle = "300", size = 4) +
  #annotate("text", x = 0, y = cut_to, label = "\\", angle = "300", size = 4) +
  #coord_cartesian(clip = "off") +
  theme(
    plot.margin=unit(c(0,0,0,40),"pt"),
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_blank()
  )