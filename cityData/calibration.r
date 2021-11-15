 # Calibration script
#
# Although this script can be sourced, I advise against it: it would take very
# long to run. The script is divided in sections which progressively build on
# one another. At the end of each section, the script saves the output that
# serves as input for the following sections: for this reason, each section
# can be run independently just as long as the output file is accessible
# ("./cityData/geodata_Rotterdam.RData").
#
#
# These are the sections:
#
# 1) Data import
# 2) Measures of segregation
# 3) World generation 
# 4) Miscellanea



# 1) Data import________________________________________________________________

# Clear working environment
rm (list = ls( )) 

# Necessary packages
require(rgdal)
require(dplyr)
require(sp)
require(spdep)
require(seg)
require(ggmap)
require(geosphere)

# Set script address as working directory.
# I noticed that this function doesn't work in some environments (possibly depending
# on the operating system). If it doesn't work, enter the working directory manually
# by replacing the command with: setwd("<current script directory>")
#setwd(dirname(parent.frame(2)$ofile))
setwd("D:/Dropbox/Progetto/w/calibration study/Scripts")
#setwd("C:/Users/Thomas/Dropbox/Progetto/w/calibration study/Scripts")

# Districts in Rotterdam
# (excluding industrial neighborhoods, empty land and municipality exclaves):
#
# Stadscentrum		WK059901
# Delfshaven    	WK059903
# Overschie	      WK059904
# Noord		        WK059905
# Hillegersberg-Schiebroek	WK059906
# Kralingen-Crooswijk	WK059908
# Feijenoord	    WK059910
# IJsselmonde	    WK059912
# Pernis				  WK059913
# Prins Alexander	WK059914
# Charlois	      WK059915
# Hoogvliet	      WK059916
#
# Dropping everything but the above-listed districts
districtsList <- c(
  "WK059901",
  "WK059903",
  "WK059904",
  "WK059905",
  "WK059906",
  "WK059908",
  "WK059910",
  "WK059912",
  "WK059913",
  "WK059914",
  "WK059915",
  "WK059916"
)
districtsNames <- c(
  "Stadscentrum",
  "Delfshaven",
  "Overschie",
  "Noord",
  "Hillegersberg-Schiebroek",
  "Kralingen-Crooswijk",
  "Feijenoord",
  "IJsselmonde",
  "Pernis",
  "Prins Alexander",
  "Charlois",
  "Hoogvliet"
)


districtsList <- c(
  "WK059901", # Rotterdam
  "WK059903",
  "WK059904",
  "WK059905",
  "WK059906",
  "WK059908",
  "WK059910",
  "WK059912",
  "WK059913",
  "WK059914",
  "WK059915",
  "WK059916",
  
  "WK036300", # Amsterdam
  "WK036301",
  "WK036302",
  "WK036303",
  "WK036304",
  "WK036305",
  "WK036306",
  "WK036307",
  
  "WK051804",  # Den Haag
  "WK051807",
  "WK051818",
  "WK051825",
  "WK051828",
  "WK051833",
  "WK051836",
  "WK051838",
  "WK051840",
  "WK051842",
  "WK051844"
)
districtsNames <- c(
  "Stadscentrum", # Rotterdam
  "Delfshaven",
  "Overschie",
  "Noord",
  "Hillegersberg-Schiebroek",
  "Kralingen-Crooswijk",
  "Feijenoord",
  "IJsselmonde",
  "Pernis",
  "Prins Alexander",
  "Charlois",
  "Hoogvliet",
  
  "Amsterdam-Centrum",   # Amsterdam
  "Amsterdam-Westpoort",
  "Amsterdam-West Westerpark",
  "Amsterdam-West Oud-West",
  "Amsterdam-Oost Zeeburg",
  "Amsterdam-West Bos en Lommer",
  "Amsterdam-West De Baarsjes",
  "Amsterdam-Noord",
  
  "Benoordenhout",  # Den Haag
  "Scheveningen",
  "Waldeck",
  "Mariahoeve en Marlot",
  "Centrum",
  "Bouwlust en Vrederust",
  "Moerwijk",
  "Laakkwartier en Spoorwijk ",
  "Wateringse Veld",
  "Ypenburg",
  "Leidschenveen"
)





# Function to convert between coordinates systems
rd2wgs84 <- function(X, Y)
{
  #http://www.dekoepel.nl/pdf/Transformatieformules.pdf
  #basispunten definieren
  X0 <- 155000.00 
  Y0 <- 463000.00 
  j0 <- 52.15517440
  l0 <- 5.38720621
  #coefficienten definieren
  K01<- 3235.65389 
  K20<- -32.58297
  K02<- -0.24750 
  K21<- -0.84978 
  K03<- -0.06550 
  K22<- -0.01709 
  K10<- -0.00738 
  K40<- 0.00530 
  K23<- -0.00039 
  K41<- 0.00033 
  K11<- -0.00012 
  
  L10<- 5260.52916
  L11<- 105.94684
  L12<- 2.45656
  L30<- -0.81885
  L13<- 0.05594
  L31<- -0.05607
  L01<- 0.01199
  L32<- -0.00256
  L14<- 0.00128
  L02<- 0.00022
  L20<- -0.00022
  L50<- 0.00026
  
  dX <- (X - X0)*10^-5
  dY <- (Y - Y0)*10^-5 
  {
    j <- j0 + 
      (
        K01*dX^0*dY^1 +
          K02*dX^0*dY^2 +
          K03*dX^0*dY^3 +
          K10*dX^1*dY^0 +
          K20*dX^2*dY^0 +
          K21*dX^2*dY^1 +
          K22*dX^2*dY^2 +
          K23*dX^1*dY^3 +
          K40*dX^2*dY^0 +
          K41*dX^2*dY^1 
      )/3600
  }
  
  {
    l <- l0 + 
      (
        L10*dX^1*dY^0 +
          L11*dX^1*dY^1 +
          L12*dX^1*dY^2 +
          L30*dX^3*dY^0 +
          L13*dX^1*dY^3 +
          L31*dX^3*dY^1 +
          L01*dX^0*dY^1 +
          L32*dX^3*dY^2 +
          L14*dX^1*dY^4 +
          L02*dX^0*dY^2 +
          L20*dX^2*dY^0 +
          L50*dX^5*dY^0 
      )/3600
  }
  wgs84<-cbind(j,l)
  return(wgs84)
}

# Load squares map. This takes really long (approx 10').
print("Loading CBS squares data. This will take approximately 10 minutes.")
cbs <- readOGR(dsn = "./2014-cbs-vierkant-100m", layer = "CBSvierkant100m201410")
cbs100 <- cbs@data

# Load map of districts & neighborhoods (takes about 2').
print("Loading CBS data on districts and neighborhoods. This will take approximately 2 minutes.")
nedb <- readOGR(dsn = "./shape 2014 versie 30/uitvoer_shape", layer = "buurt_2014")
#rot <- subset(nedb,nedb$GM_NAAM=="Rotterdam")

# Districts in Rotterdam
# (excluding industrial neighborhoods, empty land and municipality exclaves):
#
# Stadscentrum		WK059901
# Delfshaven    	WK059903
# Overschie	      WK059904
# Noord		        WK059905
# Hillegersberg-Schiebroek	WK059906
# Kralingen-Crooswijk	WK059908
# Feijenoord	    WK059910
# IJsselmonde	    WK059912
# Pernis				  WK059913
# Prins Alexander	WK059914
# Charlois	      WK059915
# Hoogvliet	      WK059916
#
# Dropping everything but the above-listed districts
rot <- subset(nedb,nedb$WK_CODE %in% districtsList)
#rot$WK_CODE
#table(cbs100$INW2014, useNA="always")
#ams<-subset(nedb, nedb$GM_NAAM == "Amsterdam")
#hag<-subset(nedb, nedb$GM_NAAM == "'s-Gravenhage")
#subset(table(ams$WK_CODE), table(ams$WK_CODE) > 0)
#subset(table(hag$WK_CODE), table(hag$WK_CODE) > 0)



# Only selecting inhabited gridcells
# (we include -99998 for now)
cbs <- subset(cbs,cbs$INW2014!=0)
cbs100 <- subset(cbs100,cbs100$INW2014!=0)

# Extract coordinate from the CBS format
cbs100$E <- as.numeric(substr(as.character(cbs100$C28992R100), 2, 5))
cbs100$N <- as.numeric(substr(as.character(cbs100$C28992R100), 7, length(as.character(cbs100$C28992R100))))

# Making centroids by hand (coordinate refers to south west corner, so I simply add 50m to x and 50m to y)
cbs100$east<-(cbs100$E*100)+50
cbs100$north<-(cbs100$N*100)+50

# Transform the cbs coordinate system to the nedb/alm coordinate system
cbs100<- cbind(cbs100, rd2wgs84(cbs100$east, cbs100$north))

# Make an identifier for each gridcell
cbs100$ncel<-1
cbs100$celid <- 1:length(cbs100$ncel)
head(cbs100)

# Make sure same coordinate system
proj4string(rot)
rot <- spTransform(rot,CRS("+proj=longlat +ellps=WGS84")) 
proj4string(rot)

# Select griddcells in Rotterdam
coordinates(cbs100) = c("l","j")
llCRS<-CRS("+proj=longlat +ellps=WGS84")
# Making spatialpoints object
loc_g<-SpatialPoints(cbs100,proj4string=llCRS)
proj4string(loc_g)<-proj4string(rot)
# Match info from Rotterdam polygon to our grids
locg2<-over(loc_g, rot)
#names(cbs100)
cbs100$BU_CODE <- locg2[,1]
cbs100$GM_NAAM <- locg2[,5]
cbs100$WK_CODE <- locg2[,3]
#names(cbs100)
#table(cbs100$GM_NAAM)

cbs100_rot <- subset(cbs100, cbs100$WK_CODE %in% districtsList)
cbs100_rot<- cbind(cbs100_rot, rd2wgs84(cbs100_rot$east, cbs100_rot$north))
#str(cbs100_rot)
#head(cbs100_rot)
#plot(cbs100_rot)
#plot(subset(cbs100_rot, cbs100_rot$WK_CODE =="WK059916"))
#cbs_rot <- subset(cbs, cbs$OBJECTID %in% cbs100_rot$OBJECTID)
#str(cbs_rot)

# From here on, we recode the variables on the demographic composition of the squares.

pautrecode_var <- c(0.00000, 97.37508, 82.06138, 69.92389, 50.40706 ,23.88861)
pwalrecode_var <- c(0.000000, 51.931746, 31.865649, 16.231917, 10.027004 , 5.405292)
pnwalrecode_var <- c(0.000000, 75.998914, 53.030019, 30.792600, 15.197071,  4.109998)
#make numeric first
cbs100_rot$p_auto2014 <-cbs100_rot$P_AUTO2014
cbs100_rot$p_wal2014 <-cbs100_rot$P_WAL2014
#cbs100_alm$p_nwal2014 <-cbs100_alm$P_NWAL2014
cbs100_rot$p_nwal2014 <-as.numeric(cbs100_rot$P_NWAL2014)
cbs100_rot$p_auto2014 <- as.numeric(cbs100_rot$p_auto2014)
cbs100_rot$p_wal2014 <- as.numeric(cbs100_rot$p_wal2014)
cbs100_rot$p_nwal2014 <- as.numeric(cbs100_rot$p_nwal2014)

#table(cbs100_alm$p_auto2014)
#table(cbs100_alm$P_AUTO2014, useNA="always")

cbs100_rot$inw2014 <- cbs100_rot$INW2014
cbs100_rot$inw2014[cbs100_rot$inw2014==-99998] <- 0

#recategorize auto
cbs100_rot$pauto2014[cbs100_rot$p_auto2014==5] <- 0
cbs100_rot$pauto2014[cbs100_rot$p_auto2014==4] <- pautrecode_var[2]
cbs100_rot$pauto2014[cbs100_rot$p_auto2014==3] <- pautrecode_var[3]
cbs100_rot$pauto2014[cbs100_rot$p_auto2014==2] <- pautrecode_var[4]
cbs100_rot$pauto2014[cbs100_rot$p_auto2014==1] <- pautrecode_var[5]
cbs100_rot$pauto2014[cbs100_rot$p_auto2014==7] <- pautrecode_var[6]
cbs100_rot$nauto2014 <- cbs100_rot$inw2014 * (cbs100_rot$pauto2014/100)
cbs100_rot$nauto2014[cbs100_rot$inw2014==0] <- 0
cbs100_rot$nauto2014b <- cbs100_rot$nauto2014

#table(cbs100_alm$nauto2014b, useNA="always")
#decide what to do with category 'geheim', should be below 5 and above 0
cbs100_rot$nauto2014b[is.na(cbs100_rot$nauto2014b)] <- 2
cbs100_rot$nauto2014mis <- ifelse(is.na(cbs100_rot$nauto2014), 1, 0)
cbs100_rot$nauto2014[is.na(cbs100_rot$nauto2014)] <- 0
#sum(cbs100_alm$nauto2014, na.rm=T)

#recategorize non-western
cbs100_rot$pnwal2014[cbs100_rot$p_nwal2014==5] <- 0
cbs100_rot$pnwal2014[cbs100_rot$p_nwal2014==4] <- pnwalrecode_var[2]
cbs100_rot$pnwal2014[cbs100_rot$p_nwal2014==3] <- pnwalrecode_var[3]
cbs100_rot$pnwal2014[cbs100_rot$p_nwal2014==2] <- pnwalrecode_var[4]
cbs100_rot$pnwal2014[cbs100_rot$p_nwal2014==1] <- pnwalrecode_var[5]
cbs100_rot$pnwal2014[cbs100_rot$p_nwal2014==7] <- pnwalrecode_var[6]
cbs100_rot$nnwal2014 <- cbs100_rot$inw2014 * (cbs100_rot$pnwal2014/100)
cbs100_rot$nnwal2014[cbs100_rot$inw2014==0] <- 0
cbs100_rot$nnwal2014b <- cbs100_rot$nnwal2014

cbs100_rot$nnwal2014b[is.na(cbs100_rot$nnwal2014b)] <- 2
cbs100_rot$nnwal2014mis <- ifelse(is.na(cbs100_rot$nnwal2014), 1, 0)
cbs100_rot$nnwal2014[is.na(cbs100_rot$nnwal2014)] <- 0

#I simply use nauto2014
cbs100_rot$nauto2014_trunc <- trunc(cbs100_rot$nauto2014)
cbs100_rot$nallo2014_trunc <- cbs100_rot$inw2014 - cbs100_rot$nauto2014_trunc
cbs100_rot$nnwal2014_trunc <- trunc(cbs100_rot$nnwal2014)
cbs100_rot$autoENwest2014_trunc <- cbs100_rot$inw2014 - cbs100_rot$nnwal2014_trunc

# Lastly, we summarize in a dataframe some district-level summary information.
#
# How many wijks in Rotterdam
cit <- subset(nedb,nedb$GM_NAAM == "Rotterdam")
cit <- subset(nedb,nedb$GM_NAAM == "Rotterdam" | nedb$GM_NAAM == "Amsterdam" | nedb$GM_NAAM == "'s-Gravenhage")
length(unique(cit$WK_CODE))

# Finding the number of 100*100m square units per district
n_squ <- c()
for (i in 1:length(districtsList)){
  cit <- subset(cbs100_rot,cbs100_rot$WK_CODE == districtsList[i])
  print (
    paste(
      "Wijk:",
      districtsList[i],
      "Square units:",
      length(unique(cit$OBJECTID))
    )
  )
  n_squ <- append(n_squ, length(unique(cit$OBJECTID)))
}
print(paste("Average", mean(n_squ)))

# Finding the population size of each district
n_pop <- n_nwa <- c()
for (i in 1:length(districtsList)){
  cit <- subset(nedb,nedb$WK_CODE == districtsList[i])
  pop <- sum(as.numeric(as.character(cit$AANT_INW)))
  # number of non-western immigrants
  cit$nwa <- c()
  for (n in 1:length(cit$P_N_W_AL)){
    if(as.numeric(as.character(cit$P_N_W_AL[n])) < 0) {cit$P_N_W_AL[n]<-0}
    cit$nwa[n] <- round(
      as.numeric(as.character(cit$AANT_INW[n])) *
        as.numeric(as.character(cit$P_N_W_AL[n])) / 100 )
  }
  nwa <- round(sum(as.numeric(as.character(cit$nwa))))
  print (
    paste(
      "Wijk:",
      districtsList[i],
      "Population:",
      pop,
      "Non-western population:",
      nwa
    )
  )
  n_pop <- append(n_pop, pop)
  n_nwa <- append(n_nwa, nwa)
}
print(paste("Population total:", sum(n_pop), "Average:", mean(n_pop)))

#Creating a dataframe with the summary information of each district
citySummary <- as.data.frame(cbind(
  districtsNames,
  districtsList,
  n_squ,
  n_pop,
  n_nwa
))
citySummary$n_pop <- as.numeric(as.character(citySummary$n_pop))
citySummary$n_squ <- as.numeric(as.character(citySummary$n_squ))
citySummary$n_nwa <- as.numeric(as.character(citySummary$n_nwa))
citySummary$p_nwa <- citySummary$n_nwa / citySummary$n_pop
names(citySummary) <- c(
  "district",
  "WK_CODE",
  "n_squ",
  "n_pop",
  "n_nwa"
)
citySummary <- citySummary[,c(1:5)]


# Interaction neighborhood
# NOTE: this is for backward compatibility - we now use a different definition
# of interaction neighborhood (see script sections 2 and 3)
#
# For each square in Rotterdam, we compute the distance
# (and associated interaction weight) with all other squares
# in the same district.
llCRS <- CRS("+proj=longlat +ellps=WGS84")
proj4string(cbs100_rot) <- llCRS
coords <- coordinates(cbs100_rot)
dlist <- 0
for (wijk in districtsList){
  focalDistrict <- subset(cbs100_rot, cbs100_rot$WK_CODE == wijk)
  
        # Note to self: make sure distance is in 100m
  #nblist <- dnearneigh(focalDistrict, 0, 1, row.names = focalDistrict$OBJECTID)
  nblist <- dnearneigh(focalDistrict, 0, 2)
  #col.nb.0.all <- dnearneigh(cbs100_alm, 0, 1.5)
  
  #include.self(nblist)
  dlist <- nbdists(nblist, focalDistrict)
}
dlist_decay <- lapply(dlist, function(x) 1/x)


# Saving script output:
save(
  cbs,
  nedb,
  file="./cityData/importedCBS.RData"
)

save(
  districtsList,
  districtsNames,
  rd2wgs84,
  cbs100_rot,
  nblist,
  dlist,
  dlist_decay,
  citySummary,
  file="geodata_Rotterdam.RData"
)
#save(cbs100_rot, districtsList, file="cityData/geodata_Rotterdam.RData")




# 2) Measures of segregation____________________________________________________
# (Global or at the level of the square cell)

# Clear working environment
rm (list = ls( )) 

# Necessary packages
require(rgdal)
require(dplyr)
require(sp)
require(spdep)
require(seg)
require(raster)
require(ggmap)
require(geosphere)

# Set script address as working directory.
#setwd(dirname(parent.frame(2)$ofile))
setwd("D:/Dropbox/Progetto/w/calibration study/Scripts")
#setwd("C:/Users/Thomas/Dropbox/Progetto/w/calibration study/Scripts")
#setwd("C:\\Users\\u838156\\surfdrive\\Shared\\ABM calibration paper\\scripts\\")

load("./cityData/geodata_Rotterdam.RData")



# White's index
#
# We calculate this proximity index for each district, using the square cells
# as base unit.
WhiteIndex <- c()
for(wijk in 1:length(districtsList)){
  print(paste0("[",wijk,"/", length(districtsList),"]"))
  print(paste0(
    "Calculating White's Index in: ",
    districtsNames[wijk],
    ". It may take a while."
  ))
  
  # We start by subsetting the square cells of the given district, making sure
  # that we discard the empty ones:
  dat <- cbs100_rot[cbs100_rot$WK_CODE==districtsList[wijk],
                    c("nauto2014", "nnwal2014")]
  dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
  
  # Next, we calculate the distance matrix, in meters.
  # We start by creating an empty matrix. Then we fill it in, row by row.
  ncells <- length(dat)
  distmat <- matrix(NA, nrow=ncells, ncol=ncells)
  for (i in 1:ncells) {
    distmat[i,] <- distVincentyEllipsoid(coordinates(dat)[i,], coordinates(dat)) 
  }
  #distmat <- distmat/10 #Here we adjust the unit of measurement, if needed.
  
  # Now we can calculate White's index:
  WhiteI <- seg::isp(dat,data=round(dat@data), nb=distmat)
  print(paste("White's Index =", WhiteI))
  WhiteIndex <- append(WhiteIndex, WhiteI)
}

# We add White's Index to the district-level summary statistics.
citySummary <- cbind(citySummary, WhiteIndex)



# Moran's I
#
# Calculated at the global or local level.

# Preparation of a fake world dataframe to test the index on.
if (FALSE){
  cityData <- subset(cbs100_rot, cbs100_rot$WK_CODE == districtsList[1])
  cityData <- cityData[cityData$nauto2014!=0 & cityData$nnwal2014!=0,]
  agents <- worldList[[1]]
  G1 <<- which(agents$group == 1, arr.ind = TRUE)
  G2 <<- which(agents$group != 1, arr.ind = TRUE)
  #world$opinion <<- rbeta(populationSize, 3, 3, ncp = 0)
  o1 <- rbeta(length(G1), 3, 3.5, ncp = 0)
  o2 <- rbeta(length(G2), 3.5, 3, ncp = 0)
  for (i in 1:nrow(world)){
    if (agents$group[i] == 1) {
      agents$opinion[i] <- o1[1]
      o1 <- o1[-1]
    } else {
      agents$opinion[i] <- o2[1]
      o2 <- o2[-1]
    }
  }
  agents$opinion <- agents$opinion * 2 - 1
  world <- agents
  proxmat <- proximityList[[1]]
  ops <- c(NA, length = length(cityData))
  for (i in 1:length(cityData)){
    cell <- subset(world, world$index == i)
    ops[i] <- mean(cell$opinion)
  }
  cityData$opinion <- ops
}




MoranI <- function(x, y=NULL, proxmat, type="global") {
  x <- scale(x)
  N <- nrow(proxmat)
  xm <- mean(x)
  if(is.null(y)){
    diag(proxmat) <- 0
    proxmat <- proxmat / rowSums(proxmat)
    y <- x
    ym <- xm
  } else {
    y <- scale(y)
    ym <- mean(y)
  }
  if (N != length(x) | length(x)!=length(y)) {
    stop("Data and/or proximity matrix have different lengths.")
  }
  
  if (type == "global"){
    
    # Global Moran's I
    W <- sum(proxmat)
    denom <- sum((x - xm) ^ 2)
    nom <- 0
    for (i in 1:N) {
      for (j in 1:N) {nom <- nom + proxmat[i,j]*(x[i] - xm)*(y[j] - ym)}
    }
    return((N*nom)/(W*denom))
  } else {
    
    # Local Moran's I - see https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1538-4632.1995.tb00338.x
    m2 <- sum((x - xm)^2) / N
    I <- c(NA, length = N)
    for (i in 1:N) {
      Zi <- x[i] - xm
      nom <- 0
      for (j in 1:N) {nom <- nom + (proxmat[i,j]*(y[j] - ym))}
      I[i] <- (Zi / m2) * nom
    }
    return(I)
  }
}

# With the following I quickly checked whether the spdep function "localmoran()"
# yields the same result as my function. Spoiler: it does.
#temp <- proxmat
#diag(temp) <- 0
#
# Whe test for local Moran's I
#test<-spdep::localmoran(x=cityData$pnwal2014, listw=mat2listw(temp))
#test[,1] == MoranI(x=cityData$pnwal2014, proxmat=proxmat, type = "local")
#
# And we test for the global index
#spdep::moran.test(x=cityData$pnwal2014, listw=mat2listw(temp))$estimate[1]
#MoranI(x=cityData$pnwal2014, proxmat=proxmat, type = "global")




# Moran scatterplot
MoranScatterplot <- function(x, lag_y, I, xlab="", ylab=""){
  x <- scale(x)
  plot(
    x,
    lag_y,
    xlab = xlab,
    ylab = ylab,
    main = paste("Moran's I =", round(I, digits = 5)),
    pch = 1#20
  )
  abline(h=0, lty = 2)
  abline(v=0, lty = 2)
  abline(0,I)
}

# We test out the univariate Moran I
MoranScatterplot(
  cityData$pnwal2014,
  MoranI(x=cityData$pnwal2014, proxmat=proxmat, type = "local"),
  MoranI(x=cityData$pnwal2014, proxmat=proxmat, type = "global"),
  xlab="Proportion of non-western immigrants",
  ylab="Lag: proportion of non-western immigrants"
)

# .. which is identical to what we obtain from the spdep functions:
temp <- proxmat
diag(temp) <-0
MoranScatterplot(
  cityData$pnwal2014,
  spdep::localmoran(x=cityData$pnwal2014, listw=mat2listw(temp))[,1],
  spdep::moran.test(x=cityData$pnwal2014, listw=mat2listw(temp))$estimate[1],
  xlab="Proportion of non-western immigrants",
  ylab="Lag: proportion of non-western immigrants"
)
rm(temp)


# And then we try out the bivariate Moran's I
MoranScatterplot(
  cityData$pnwal2014,
  MoranI(x=cityData$pnwal2014, y=cityData$opinion, proxmat, type = "local"),
  MoranI(x=cityData$pnwal2014, y=cityData$opinion, proxmat, type = "global"),
  xlab="Non-western immigrants",
  ylab="Lag: opinion"
)




if(FALSE){ # inspecting a different possible implementation
# ***************************** Temporary **************************************
MoranI_Barbosa <- function(x, y = NULL, W){
  if(is.null(y)){
    y = x
    diag(W) <- 0
  }
  xp <- scale(x)[,1]#(x - mean(x, na.rm=T))/sd(x, na.rm=T)
  yp <- scale(y)[,1]#(y - mean(y, na.rm=T))/sd(y, na.rm=T)
  W[which(is.na(W))] <- 0
  W <- W / rowSums(W)#############################################
  n <- nrow(W)
  
  global <- (xp%*%W%*%yp)/(n - 1)
  #global <- (xp%*%W%*%yp)/(n )
  local  <- (xp*W%*%yp)
  
  list(global = global, local  = as.numeric(local))
}

simula_moran <- function(x, y = NULL, W, nsims = 1000){
  
  if(is.null(y)) y = x
  
  n   = nrow(W)
  IDs = 1:n
  
  xp <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  W[which(is.na(W))] <- 0
  
  global_sims = NULL
  local_sims  = matrix(NA, nrow = n, ncol=nsims)
  
  ID_sample = sample(IDs, size = n*nsims, replace = T)
  
  y_s = y[ID_sample]
  y_s = matrix(y_s, nrow = n, ncol = nsims)
  y_s <- (y_s - apply(y_s, 1, mean))/apply(y_s, 1, sd)
  
  global_sims  <- as.numeric( (xp%*%W%*%y_s)/(n - 1) )
  local_sims  <- (xp*W%*%y_s)
  
  list(global_sims = global_sims,
       local_sims  = local_sims)
}
simula_moran(x=cityData$pnwal2014, W=proxmat)


MoranI(x=cityData$pnwal2014, y=cityData$opinion, proxmat=proxmat, type = "global")
MoranI_Barbosa(x=cityData$pnwal2014, y=cityData$opinion, W=proxmat)$global

MoranI(x=cityData$pnwal2014, proxmat=proxmat, type = "local")
MoranI_Barbosa(x=cityData$pnwal2014, W=proxmat)$local


MoranI(x=cityData$pnwal2014, y=cityData$opinion, proxmat=proxmat, type = "local")
MoranI_Barbosa(x=cityData$pnwal2014, y=cityData$opinion, W=proxmat)$local



test <- simW[[1]]
dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[1])
dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
proxmat <- proximityList2[[1]]
a <- MoranI(x=dat$pnwal2014, y=r$opinion, proxmat=proxmat, type = "local")
b <- MoranI_Barbosa(x=dat$pnwal2014, y=r$opinion, W=proxmat)$local

hist(a)
hist(b)

table(abs(a)>1)
table(abs(b)>1)


temp <- proxmat
diag(temp) <- 0
spdep::moran.test(x=cityData$pnwal2014, listw=mat2listw(temp))$estimate[1]
MoranI(x=cityData$pnwal2014, proxmat=proxmat, type = "global")
MoranI_Barbosa(x=cityData$pnwal2014, W=proxmat)$global
rm(temp)


temp<-as.data.frame(temp)



# ****************************** ******* **************************************
}


# Saving script output:
save(
  districtsList,
  districtsNames,
  rd2wgs84,
  MoranI,
  MoranScatterplot,
  cbs100_rot,
  nblist,
  dlist,
  dlist_decay,
  citySummary,
  file="./cityData/geodata_Rotterdam.RData"
)





# 3) World generation __________________________________________________________
#
# Here we create a simulated world of agents calibrated on the demographics of
# Rotterdam.

# Clear working environment
rm (list = ls( )) 

# Necessary packages
require(spdep)
require(geosphere)

# Set script address as working directory.
#setwd(dirname(parent.frame(2)$ofile))
setwd("D:/Dropbox/Progetto/w/calibration study/Scripts")
#setwd("C:/Users/Thomas/Dropbox/Progetto/w/calibration study/Scripts")
#setwd("C:\\Users\\u838156\\surfdrive\\Shared\\ABM calibration paper\\scripts\\")

load("./cityData/importedCBS.RData")
load("./cityData/geodata_Rotterdam.RData")

# Function to plot matrix for quick visual inspection:
printMat <- function(x){
  return(
    image(x, col=heat.colors(10000), zlim=c(min(x, na.rm=T),max(x, na.rm=T)))
  )
}

# Let's try out distance decay functions:
plotDistanceDecay <- function(){
  plot(exp(- (c(1:1000) / 100)),
       type = "l",
       col = "blue",
       xlab = "distance in meters",
       ylab = "probability of interaction")
  lines(exp(- (c(1:1000) / 1000)), col = "green")
  lines(exp(- (c(1:1000) / 10)), col = "red")
  lines(exp(- c(1:1000)), col = "black")
  legend(600,
         1,
         legend=c(
           "p=exp(-distance/1000)",
           "p=exp(-distance/100)",
           "p=exp(-distance/10)",
           "p=exp(-distance)"
         ),
         col=c("green", "blue", "red", "black"),
         lty=1
  )
}
plotDistanceDecay()
rm(plotDistanceDecay)
# I would say that distance/100 looks more plausible, so I proceed with that one.



createDistrict <- function(district){
  print(paste("District:", districtsNames[district]))
  # We select a district and a population size. We filter out square cells with
  # no inhabitants.
  
  wijk <- districtsList[district]
  #wijk <- districtsList[14] #*******************************************************************************
  scalePopulation <- 1#0.1
  cityData <- subset(cbs100_rot, cbs100_rot$WK_CODE == wijk)
  cityData <- cityData[cityData$nauto2014!=0 & cityData$nnwal2014!=0,]
  
  # We create the agents as the rows of a dataframe "world". To do so, we loop 
  # through all cells that belong to the district we selected.
  print("Creating world dataframe")
  world <- data.frame()
  for (i in 1:length(cityData@data$inw2014)){
    
    # First, we calculate the total number of agents we need to reproduce the
    # cell.
    cell <- cityData@data[i,]
    cell$inw2014 <- round(cell$inw2014 * scalePopulation)
    cell$nnwal2014_trunc <- round(cell$nnwal2014_trunc * scalePopulation)
    
    # Then we create agents' features using vectors as long as the number of
    # cell inhabitants.
    # We start with the group, assigning value +1 to non western immigrants, and
    # -1 to the rest.
    group <- rep(
      c(-1, 1),
      times = c((cell$inw2014-cell$nnwal2014_trunc), cell$nnwal2014_trunc))
    
    # Then we create vectors with the cell information. This includes: the cell
    # ID (location) and the coordinates (j, l, east, north).  Technically we
    # would only need the vector "location" to retrieve the rest of the
    # information in the cell dataset (cbs100_rot). However, for convenience we
    # merge it already here.
    location <- rep(cell$OBJECTID, times = cell$inw2014)
    index <- rep(i, times = cell$inw2014) # The index is just the number of the 
                                          # row in the cell dataset. We save it
                                          # in this vector so that later we can
                                          # use it to access the distance matrix
                                          # more quickly.
    x_coor <- rep(cell$j, times = cell$inw2014)
    y_coor <- rep(cell$l, times = cell$inw2014)
    east <- rep(cell$east, times = cell$inw2014)
    north <- rep(cell$north, times = cell$inw2014)
    
    # Finally, we put this information in the "world" dataset.
    world <- rbind(world, (data.frame(cbind(
      location,
      index,
      x_coor,
      y_coor,
      east,
      north,
      group
    ))))
  }
  populationSize <- nrow(world)
  
  # Next, we give agents an opinion.
  #This is randomly generated, so we need a seed.
  seed  <-  sample.int(999999999, size = 1)
  #seed  <-  123456789
  set.seed (seed)
  world$opinion <- rbeta(populationSize, 3, 3, ncp = 0)
  world$opinion <- world$opinion * 2 - 1
  #world$opinion=runif(populationSize, -1, 1)
  
  #world$opinion[world$index==1] <- -1
  #world$opinion[world$index<10] <- -1
  #world$opinion[world$index>80] <- 1
  #hist(world$opinion)
  
  # The world is now complete. What we need next is a proximity matrix between 
  # cells, that we will use to calculate the probabilities of interactions 
  # between agents.
  #
  # We first prepare the ingredients that we need to calculate a distance matrix
  print("Calculating proximity matrix")
  dat <- cityData[cityData, c("nauto2014", "nnwal2014")]
  ncells <- length(dat)
  distmat <- matrix(NA, nrow=ncells, ncol=ncells)
  proxmat1 <- matrix(NA, nrow=ncells, ncol=ncells)
  proxmat2 <- matrix(NA, nrow=ncells, ncol=ncells)
  proxmat3 <- matrix(NA, nrow=ncells, ncol=ncells)
  
  # Next, we fill in the matrix row by row with the distances (in meters)
  for (i in 1:ncells) {
    distmat[i,] <- distVincentyEllipsoid(coordinates(dat)[i,], coordinates(dat)) 
  }
  #proxmat <- proxmat/10 #Here we can adjust the unit of measurement, if needed.
  #printMat(proxmat)
  
  
  # We apply the distance decay function, so that the matrix expresses
  # proximity instead of distance.
  # NOTE: For agents who live in the same cell, we want to avoid to assume their
  # distance is 0, because their proximity would be maximal. Instead, we
  # assume that their proximity is the average distance between all points
  # in a square sized 100*100 meters. That would be about 52.14m. Therefore,
  # we assume that the distance between a cell and itself is about 52.14m.
  proxmat1 <- exp(-distmat/10)
  proxmat2 <- exp(-distmat/100)
  proxmat3 <- exp(-distmat/1000)
  diag(proxmat1) <- exp(-52.140543316/10)   #######
  diag(proxmat2) <- exp(-52.140543316/100)  #######
  diag(proxmat3) <- exp(-52.140543316/1000) #######
  #printMat(proxmat)
  
  # Normalize rows:
  #sum <- rowSums(proxmat)
  for (i in 1:length(distmat[1,])){
    #proxmat1[i,] <- proxmat1[i,] / sum[i]
    #proxmat2[i,] <- proxmat2[i,] / sum[i]
    #proxmat3[i,] <- proxmat3[i,] / sum[i]
    proxmat1[i,] <- proxmat1[i,] / sum(proxmat1[i,])
    proxmat2[i,] <- proxmat2[i,] / sum(proxmat2[i,])
    proxmat3[i,] <- proxmat3[i,] / sum(proxmat3[i,])
  }
  
  
  # NOTE: this section of commented code is switched off because made redundant.
  # Now exposure is calculated more fully at a later stage. Did not remove from
  # here for reference
  # (OLD implementation) Exposure to ingroup vs outgroup agents
  #
  # At this point we have all we need to calculate the local measures of
  # segregation. We start by calculating agent's exposure to other ingroup or
  # ougroup agents.
  # We do that by calculating two matrices, where we store the exposure to
  # natives and western immigrants (eAllo), and the exposure to non-western
  # immigrants (eNWI).
  #proxmat <- proxmat2 #**********************************************************
  #l <- nrow(proxmat)
  #eAllo <- matrix(NA, nrow=l, ncol=l)
  #eNWI <- matrix(NA, nrow=l, ncol=l)
  #eAllo_ingr <- matrix(NA, nrow=l, ncol=l)
  #eNWI_ingr <- matrix(NA, nrow=l, ncol=l)
  #for (i in 1:l){
  #  for (j in 1:l){
  #    
  #    # The exposure of a cell i to a group (Allo v NWI) in cell j is simply the
  #    # product of the relative proximity between the two cells, and the number
  #    # of members of the target group in cell j.
  #    eAllo[i,j] <- proxmat[i,j] * (cityData$inw2014[j] - cityData$nnwal2014_trunc[j])
  #    eNWI[i,j] <- proxmat[i,j] * cityData$nnwal2014_trunc[j]
  #    
  #    if (i==j){
  #      eAllo_ingr[i,j] <- proxmat[i,j] * (cityData$inw2014[j] - cityData$nnwal2014_trunc[j] - 1)
  #      eNWI_ingr[i,j] <- proxmat[i,j] * (cityData$nnwal2014_trunc[j] - 1)
  #    } else {
  #      eAllo_ingr[i,j] <- eAllo[i,j]
  #      eNWI_ingr[i,j] <- eNWI[i,j] 
  #    }
  #  }
  #}
  #
  # We normalize
  #vAllo <- rowSums(eAllo)
  #vNWI <- rowSums(eNWI)
  #vAllo_ingr <- rowSums(eAllo_ingr)
  #vNWI_ingr <- rowSums(eNWI_ingr)
  #sums_egoAllo <- vAllo_ingr + vNWI
  #sums_egoNWI <- vAllo + vNWI_ingr
  #eAllo_egoAllo <- eAllo_ingr / matrix(sums_egoAllo, nrow=l, ncol=l)
  #eAllo_egoNWI <- eAllo / matrix(sums_egoNWI, nrow=l, ncol=l)
  #eNWI_egoAllo <- eNWI / matrix(sums_egoAllo, nrow=l, ncol=l)
  #eNWI_egoNWI <- eNWI_ingr / matrix(sums_egoNWI, nrow=l, ncol=l)
  #
  # And finally we extract the vectorized indices of exposure to the two groups.
  #vAllo_egoAllo <- rowSums(eAllo_egoAllo)
  #vAllo_egoNWI <- rowSums(eAllo_egoNWI)
  #vNWI_egoAllo <- rowSums(eNWI_egoAllo)
  #vNWI_egoNWI <- rowSums(eNWI_egoNWI)
  #
  # The last thing we need to do is to assign the correct value for the exposure
  # (to ingroup and outgroup) to every agent.
  # For agents from group 1 (non-western), the ingroup exposure is of course the
  # exposure to other non-western (vNWI). For group -1, it's the other way
  # around.
  #for (i in 1:populationSize){
  #  if (world$group[i] == 1){
  #    world$exposureIngroup[i] <- vNWI_egoNWI[world$index[i]]
  #    world$exposureOutgroup[i] <- vAllo_egoNWI[world$index[i]]
  #  } else {
  #    world$exposureIngroup[i] <- vAllo_egoAllo[world$index[i]]
  #    world$exposureOutgroup[i] <- vNWI_egoAllo[world$index[i]]
  #  }
  #}
  #rm(
  #  l,
  #  eAllo, eAllo_ingr, eAllo_egoAllo, eAllo_egoNWI,
  #  eNWI, eNWI_ingr, eNWI_egoAllo, eNWI_egoNWI,
  #  vAllo, vAllo_ingr, vAllo_egoAllo, vAllo_egoNWI,
  #  vNWI, vNWI_ingr, vNWI_egoAllo, vNWI_egoNWI,
  #sums_egoAllo, sums_egoNWI
  #)
  
  print("Calculating local measures of segregation")
  
  # Local Moran's I.
  #
  # We calculate the index at the level of the square cell. Depending on "type",
  # this function returns the global value (type = "global"), or a vector of
  # values of I measured at each spatial unit (type = "local")
  
  # Local group clustering (local I_prop)
  I <- MoranI(
    x=cityData$pnwal2014,
    proxmat=proxmat,
    type = "local"
  )
  
  # Local group segregation (local biv_I_count)
  Is <- MoranI(
    x=cityData$inw2014 - cityData$nnwal2014,
    y=cityData$nnwal2014,
    proxmat=proxmat,
    type = "local"
  )
  
  # Local group clustering (local I_count)
  Iec <- MoranI(
    x=cityData$nnwal2014,
    proxmat=proxmat,
    type = "local"
  )
  
  
  for (i in 1:populationSize){
    world$localI_prop[i] <- I[world$index[i]]
    world$localI_biv_count[i] <- Is[world$index[i]]
    world$localI_count[i] <- Iec[world$index[i]]
    world$p_nwa[i] <- cityData$pnwal2014[world$index[i]]
  }
  rm(I,Is,Iec)
  
  
  agents <<- world
  proxmat1 <<- proxmat1
  proxmat2 <<- proxmat2
  proxmat3 <<- proxmat3
}

#createDistrict(9)



# With this code we generate a world dataframe and a proximity matrix for each
# district. We store this information in two lists, worldList and proximityList.
# NOTE: the execution of this code takes quite long (about 15 minutes).
worldList <- list()
proximityList1 <- proximityList2 <- proximityList3 <- list()

for (i in 1:length(districtsList)){
#for (i in c(1,3)){
  #world <- data.frame()
  #proxmat <- matrix()
  print(paste0("[", i, "/", length(districtsList), "]"))
  createDistrict(i)
  worldList[[i]] <- agents
  proximityList1[[i]] <- proxmat1
  proximityList2[[i]] <- proxmat2
  proximityList3[[i]] <- proxmat3
  print("District creation complete.")
}


# Now that we have generated the districts, we can calulate Moran's I  for every
# district and add it to the district descriptives:
globalI_propNWA <- c()
globalI_countNWA <- c()
#globalMoranI_nat <- c()
globalI_countNAT <-c()
globalI_biv_prop <- c()
globalI_biv_count <- c()
globalMoranI_opinion <- c()

for(wijk in 1:length(districtsList)){
  print(paste0("[",wijk,"/", length(districtsList),"]"))
  print(paste0(
    "Calculating global segregation measures in: ",
    districtsNames[wijk]
  ))
  
  # We subset the square cells of the given district, making sure
  # that we discard the empty ones:
  dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[wijk])
  dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
  agents <- worldList[[wijk]]
  proxmat <- proximityList2[[wijk]] #****************************************************
  
  # And we assign a group-biased opinion, using the same procedure that we use
  # in the ABM, when initialOpinionDistribution = "groupBias"
  G1 <<- which(agents$group == 1, arr.ind = TRUE)
  G2 <<- which(agents$group != 1, arr.ind = TRUE)
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
  agents$opinion <- agents$opinion * 2 - 1
  
  # We calculate the average opinion in each square unit
  ops <- c(NA, length = length(dat))
  for (i in 1:length(dat)){
    cell <- subset(agents, agents$index == i)
    ops[i] <- mean(cell$opinion)
  }
  dat$opinion <- ops
  
  # Finally, we calculate the Morans'I on both the group variable and the
  # opinion variable.
  Inwa <- MoranI(x=dat$pnwal2014, proxmat=proxmat, type = "global")
  IcountNWA <- MoranI(x=dat$nnwal2014, proxmat=proxmat, type = "global")
  #Inat <- MoranI(x=100 - dat$pnwal2014, proxmat=proxmat, type = "global")
  IcountNAT <- MoranI(x=dat$inw2014 - dat$nnwal2014, proxmat=proxmat, type = "global")
  IgroupbivP <- MoranI(x=100 - dat$pnwal2014, y=dat$pnwal2014, proxmat=proxmat, type = "global")
  IgroupbivC <- MoranI(x=dat$inw2014 - dat$nnwal2014, y=dat$nnwal2014, proxmat=proxmat, type = "global")
  
  o <- MoranI(x=dat$opinion, proxmat=proxmat, type = "global")
  globalI_propNWA <- append(globalI_propNWA, Inwa)
  globalI_countNWA <- append(globalI_countNWA, IcountNWA)
  #globalMoranI_nat <- append(globalMoranI_nat, Inat)
  globalI_countNAT <- append(globalI_countNAT, IcountNAT)
  globalI_biv_prop <- append(globalI_biv_prop, IgroupbivP)
  globalI_biv_count <- append(globalI_biv_count, IgroupbivC)
  globalMoranI_opinion <- append(globalMoranI_opinion, o)
}


# We conclude by adding the precalculated segregation measures to the district
# descriptives.
citySummary$p_nwa <- citySummary$n_nwa / citySummary$n_pop
citySummary$globalI_propNWA <- globalI_propNWA
citySummary$globalI_biv_prop <- globalI_biv_prop
citySummary$globalI_countNWA <- globalI_countNWA
citySummary$globalI_countNAT <- globalI_countNAT
citySummary$globalI_biv_count <- globalI_biv_count
 


# Here we measure exposure:
# Calculating exposure index for every district, for the two groups together and
# separately, for the three distance decay functions.
library(compiler)
enableJIT(1)
calcExposure <- function(agents, dat, proxmat){
  l <- nrow(proxmat)
  eAllo <- matrix(NA, nrow=l, ncol=l)
  eNWI <- matrix(NA, nrow=l, ncol=l)
  eAllo_ingr <- matrix(NA, nrow=l, ncol=l)
  eNWI_ingr <- matrix(NA, nrow=l, ncol=l)
  for (i in 1:l){
    for (j in 1:l){
      
      # The exposure of a cell i to a group (Allo v NWI) in cell j is simply the
      # product of the relative proximity between the two cells, and the number
      # of members of the target group in cell j.
      eAllo[i,j] <- proxmat[i,j] * (dat$inw2014[j] - dat$nnwal2014_trunc[j])
      eNWI[i,j] <- proxmat[i,j] * dat$nnwal2014_trunc[j]
      
      if (i==j){
        eAllo_ingr[i,j] <-
          proxmat[i,j] * (dat$inw2014[j] - dat$nnwal2014_trunc[j] - 1)
        eNWI_ingr[i,j] <-
          proxmat[i,j] * (dat$nnwal2014_trunc[j] - 1)
      } else {
        eAllo_ingr[i,j] <- eAllo[i,j]
        eNWI_ingr[i,j] <- eNWI[i,j] 
      }
    }
  }
  
  # We normalize
  vAllo <- rowSums(eAllo)
  vNWI <- rowSums(eNWI)
  vAllo_ingr <- rowSums(eAllo_ingr)
  vNWI_ingr <- rowSums(eNWI_ingr)
  sums_egoAllo <- vAllo_ingr + vNWI
  sums_egoNWI <- vAllo + vNWI_ingr
  eAllo_egoAllo <- eAllo_ingr / matrix(sums_egoAllo, nrow=l, ncol=l)
  eAllo_egoNWI <- eAllo / matrix(sums_egoNWI, nrow=l, ncol=l)
  eNWI_egoAllo <- eNWI / matrix(sums_egoAllo, nrow=l, ncol=l)
  eNWI_egoNWI <- eNWI_ingr / matrix(sums_egoNWI, nrow=l, ncol=l)
  
  # And finally we extract the vectorized indices of exposure to the two groups.
  vAllo_egoAllo <- rowSums(eAllo_egoAllo)
  vAllo_egoNWI <- rowSums(eAllo_egoNWI)
  vNWI_egoAllo <- rowSums(eNWI_egoAllo)
  vNWI_egoNWI <- rowSums(eNWI_egoNWI)
  
  # The last thing we need to do is to assign the correct value for the exposure
  # (to ingroup and outgroup) to every agent.
  # For agents from group 1 (non-western), the ingroup exposure is of course the
  # exposure to other non-western (vNWI). For group -1, it's the other way
  # around.
  for (i in 1:nrow(agents)){
    if (agents$group[i] == 1){
      agents$exposureIngroup[i] <- vNWI_egoNWI[agents$index[i]]
      agents$exposureOutgroup[i] <- vAllo_egoNWI[agents$index[i]]
    } else {
      agents$exposureIngroup[i] <- vAllo_egoAllo[agents$index[i]]
      agents$exposureOutgroup[i] <- vNWI_egoAllo[agents$index[i]]
    }
  }
  return(as.data.frame(cbind(
    agents$group,
    agents$exposureIngroup,
    agents$exposureOutgroup
  )))
}
#dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[1])
#dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
#test <- calcExposure(worldList[[1]], dat, proximityList1[[1]])

# initializing variables
citySummary$expIngr1 <- citySummary$expIngr2 <- citySummary$expIngr3 <- NA
citySummary$expOutgr1 <- citySummary$expOutgr2 <- citySummary$expOutgr3 <- NA
citySummary$NWAexpIngr1 <- citySummary$NWAexpIngr2 <- citySummary$NWAexpIngr3 <- NA
citySummary$NWAexpOutgr1 <- citySummary$NWAexpOutgr2 <-citySummary$NWAexpOutgr3 <- NA
citySummary$NATexpIngr1 <- citySummary$NATexpIngr2 <- citySummary$NATexpIngr3 <- NA
citySummary$NATexpOutgr1 <- citySummary$NATexpOutgr2 <- citySummary$NATexpOutgr3 <- NA


for (wijk in 1:length(districtsList)){
  #for(wijk in 9){
  print(paste("Calculating exposure in district", districtsNames[wijk]))
  dat <- subset(cbs100_rot, cbs100_rot$WK_CODE==districtsList[wijk])
  dat <- dat[dat$nauto2014!=0 & dat$nnwal2014!=0,]
  agents <- worldList[[wijk]]
  
  # with distance-decay 1
  proxmat <- proximityList1[[wijk]]
  exposure <- calcExposure(agents, dat, proxmat)
  worldList[[wijk]]$expIngr1 <- exposure[,2]
  worldList[[wijk]]$expOutgr1 <- exposure[,3]
  citySummary$expIngr1[wijk] <- mean(exposure[,2])
  citySummary$expOutgr1[wijk] <- mean(exposure[,3])
  citySummary$NWAexpIngr1[wijk] <- mean(subset(exposure,exposure[,1]==1)[,2])
  citySummary$NWAexpOutgr1[wijk] <- mean(subset(exposure,exposure[,1]==1)[,3])
  citySummary$NATexpIngr1[wijk] <- mean(subset(exposure,exposure[,1]!=1)[,2])
  citySummary$NATexpOutgr1[wijk] <- mean(subset(exposure,exposure[,1]!=1)[,3])
  
  # with distance-decay 2  
  proxmat <- proximityList2[[wijk]]
  exposure <- calcExposure(agents, dat, proxmat)
  worldList[[wijk]]$expIngr2 <- exposure[,2]
  worldList[[wijk]]$expOutgr2 <- exposure[,3]
  citySummary$expIngr2[wijk] <- mean(exposure[,2])
  citySummary$expOutgr2[wijk] <- mean(exposure[,3])
  citySummary$NWAexpIngr2[wijk] <- mean(subset(exposure,exposure[,1]==1)[,2])
  citySummary$NWAexpOutgr2[wijk] <- mean(subset(exposure,exposure[,1]==1)[,3])
  citySummary$NATexpIngr2[wijk] <- mean(subset(exposure,exposure[,1]!=1)[,2])
  citySummary$NATexpOutgr2[wijk] <- mean(subset(exposure,exposure[,1]!=1)[,3])
  
  # with distance-decay 3
  proxmat <- proximityList3[[wijk]]
  exposure <- calcExposure(agents, dat, proxmat)
  worldList[[wijk]]$expIngr3 <- exposure[,2]
  worldList[[wijk]]$expOutgr3 <- exposure[,3]
  citySummary$expIngr3[wijk] <- mean(exposure[,2])
  citySummary$expOutgr3[wijk] <- mean(exposure[,3])
  citySummary$NWAexpIngr3[wijk] <- mean(subset(exposure,exposure[,1]==1)[,2])
  citySummary$NWAexpOutgr3[wijk] <- mean(subset(exposure,exposure[,1]==1)[,3])
  citySummary$NATexpIngr3[wijk] <- mean(subset(exposure,exposure[,1]!=1)[,2])
  citySummary$NATexpOutgr3[wijk] <- mean(subset(exposure,exposure[,1]!=1)[,3])
}
enableJIT(0)


# Cleaning citySummary from redundant variables
#for(i in 1:length(districtsList)){
#  worldList[[i]]$exposureIngroup <- NULL
#  worldList[[i]]$exposureOutgroup <- NULL
#}
#citySummary$NWAexpIngroup <- citySummary$NWAexpOutgroup <-
#  citySummary$NATexpIngroup <- citySummary$NATexpOutgroup <-
#  citySummary$expIngroup <- citySummary$expOutgroup <- NULL

# Calculating each district's proportion of agents with very high exposure to
# the outgroup
citySummary$propHighlyExposed1 <-
  citySummary$propHighlyExposed2 <-
  citySummary$propHighlyExposed3 <- NA
for (w in 1:nrow(citySummary)){
  #threshold = 0.5
  a <- worldList[[w]]
  pop <- nrow(a)
  threshold1 <- mean(a$expOutgr1) + sd(a$expOutgr1)
  threshold2 <- mean(a$expOutgr2) + sd(a$expOutgr2)
  threshold3 <- mean(a$expOutgr3) + sd(a$expOutgr3)
  he1 <- length(a$expOutgr1[a$expOutgr1 > threshold1])
  he2 <- length(a$expOutgr2[a$expOutgr2 > threshold2])
  he3 <- length(a$expOutgr3[a$expOutgr3 > threshold3])
  citySummary$propHighlyExposed1[w] <- he1 / pop
  citySummary$propHighlyExposed2[w] <- he2 / pop
  citySummary$propHighlyExposed3[w] <- he3 / pop
}
par(mfrow=c(1,3))
plot(citySummary$expOutgr1, citySummary$propHighlyExposed1)
plot(citySummary$expOutgr2, citySummary$propHighlyExposed2)
plot(citySummary$expOutgr3, citySummary$propHighlyExposed3)

#Calculating variables for table 1
citySummary$temp <- round(citySummary$propHighlyExposed2 * 100, digits=2)
number_NWA <- c(NA)
exposures <- c()
for (w in 1:nrow(citySummary)){
  a <- worldList[[w]]
  number_NWA[w] <- nrow(a[a$group==1,])
  exposures <- c(exposures, a$expOutgr2)
}
# % non western in Rotterdam
round(sum(number_NWA) / sum(citySummary$n_pop) * 100, digits = 2)
# Mean exposure (s=100) and % agents with extremely high exposure (Rotterdam)
mean(exposures)
round(sum(exposures > mean(exposures + sd(exposures))) / length(exposures) * 100, digits = 2)

  
  
# Saving script output:
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

# Inspecting correlation between segregation measures
#
# Global level: here we compare districts
require(ggplot2)
ggplot(
  citySummary,
  aes(
    x = globalMoranI_group,
    y = WhiteIndex,
    color = expIngroup
  )
) +
  geom_point( size = 4) +
  scale_color_gradient(
    low = "blue",
    high = "red"
  ) +
  labs(
    color = "Exposure to ingroup",
    x = "Global Moran I",
    y = "White Index"
  )

plot(
  citySummary$WhiteIndex, 
  citySummary$expIngroup,
  main = "Districts: segregation measures",
  xlab = "White Index",
  ylab = "Average exposure to ingroup"
)

plot(
  citySummary$globalMoranI_group, 
  citySummary$expIngroup,
  main = "Districts: segregation measures",
  xlab = "Global Moran I: group clustering",
  ylab = "Average exposure to ingroup"
)

plot(
  citySummary$globalMoranI_group, 
  citySummary$expOutgroup,
  main = "Districts: segregation measures",
  xlab = "Global Moran I: group clustering",
  ylab = "Average exposure to outgroup"
)

plot(
  citySummary$globalMoranI_opinion, 
  citySummary$expIngroup,
  main = "Districts: segregation measures",
  xlab = "Global Moran I: opinion clustering",
  ylab = "Average exposure to ingroup"
)


citySummary <- citySummary1  # Steep decay
citySummary <- citySummary2  # Soft distance decay


plot(citySummary$globalI_propNWA,citySummary$globalI_countNWA)#,xlim=c(-0.05,0.6),ylim=c(-0.05,0.6))
plot(citySummary$globalI_propNWA,citySummary$globalI_countNAT)#,xlim=c(-0.05,0.6),ylim=c(-0.05,0.6))
plot(citySummary$globalI_countNWA,citySummary$globalI_countNAT)#,xlim=c(-0.05,0.6),ylim=c(-0.05,0.6))

plot(citySummary$globalI_biv_count,citySummary$globalI_countNAT)
plot(citySummary$globalI_biv_prop,citySummary$globalI_biv_count)
plot(citySummary$expOutgroup,citySummary$globalI_biv_prop)
plot(citySummary$expIngroup,citySummary$globalI_biv_prop)
plot(citySummary$expOutgroup,citySummary$globalI_biv_count)
plot(citySummary$expIngroup,citySummary$globalI_biv_count)

plot(citySummary$globalI_countNWA,citySummary$expOutgroup)

plot(citySummary$NWAexpIngroup,citySummary$globalI_biv_count)  ####
plot(citySummary$NWAexpOutgroup,citySummary$globalI_biv_count) ####

plot(citySummary$NATexpIngroup,citySummary$globalI_biv_count)
plot(citySummary$NATexpOutgroup,citySummary$globalI_biv_count)

plot(citySummary$NATexpOutgroup, citySummary$WhiteIndex)
 

plot(citySummary$p_nwa,citySummary$globalI_biv_count)
plot(citySummary$p_nwa,citySummary$NATexpOutgroup)
plot(citySummary$p_nwa,citySummary$NWAexpIngroup)
plot(citySummary$p_nwa,citySummary$expOutgroup)

plot(citySummary$p_nwa,globalI_countNWA)


summary(lm(
  citySummary$NWAexpIngroup ~ citySummary$globalI_biv_count
  + citySummary$p_nwa
  + citySummary$globalI_countNWA
))


if(FALSE){
  save(
    citySummary1,
    citySummary2,
    file="./cityData/desc_distance_decay.RData"
  )
  load("./cityData/desc_distance_decay.RData")
}



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
# Local level: here we compare agents in a give district.
# The quickest way to create the agent dataset is to have the ABM create it.
# So we load the ABM:
source("script NI&PA.R")

# And then we have it initialize the model based on a given wijk (in ths case,
# 1= Stadscentrum). The ABM outputs the agents dataframe.
run (wijk = 1, timeMax = 0)

# We calculate the average opinion in the cells
cityData <- subset(cbs100_rot, cbs100_rot$WK_CODE == districtsList[1])
cityData <- cityData[cityData$nauto2014!=0 & cityData$nnwal2014!=0,]
ops <- c(NA, length = length(cityData))
for (i in 1:length(cityData)){
  cell <- subset(agents, agents$index == i)
  ops[i] <- mean(cell$opinion)
}
cityData$opinion <- ops
#cityData$opinion[1:10] <- 1
#cityData$opinion 100:110 <- -1

# Based on the cell average opinion, we calculate the local (cell-level) Moran I.
# We then assign the value for local I to each agent in the agentset.
I <- MoranI(x=cityData$opinion, proxmat=proxmat, type = "local")
for (i in 1:populationSize){agents$localMoranI_opinion[i] <- I[agents$index[i]]}
rm(I)

# Now we are ready to plot exposure and local Moran's I at the agent level.
plot(
  agents$localMoranI, 
  agents$exposureIngroup,
  col = ifelse(agents$group == 1, "orange", "black"),
  main = "Agents: segregation measures",
  xlab = "Local Moran I: group clustering",
  ylab = "Exposure to ingroup"
)






o1 <- rbeta(length(G1), 3, 10, ncp = 0)
o2 <- rbeta(length(G2), 10, 3, ncp = 0)
for (i in 1:populationSize){
  if (agents$group[i] == 1) {
    agents$opinion[i] <- o1[1]
    o1 <- o1[-1]
  } else {
    agents$opinion[i] <- o2[1]
    o2 <- o2[-1]
  }
}
agents$opinion <- agents$opinion * 2 - 1
I <- MoranI(x=dat$opinion, proxmat=proxmat, type = "local")
for (i in 1:populationSize){agents$localMoranI_opinion[i] <- I[agents$index[i]]}
rm(I)

plot(
  agents$localMoranI_opinion, 
  agents$exposureIngroup,
  col = ifelse(agents$group == 1, "orange", "black"),
  main = "Agents: segregation measures",
  xlab = "Local Moran I: opinion clustering",
  ylab = "Exposure to ingroup"
)


table(subset(agents, agents$index==233)$group)
table(subset(agents, agents$index==212)$group)
table(subset(agents, agents$index==238)$group)
table(subset(agents, agents$index==234)$group)
table(subset(agents, agents$index==204)$group)
table(subset(agents, agents$index==232)$group)
table(subset(agents, agents$index==240)$group)




# 4) Miscellanea________________________________________________________________

# Clear working environment
rm (list = ls( )) 

# Necessary packages
require(rgdal)
require(dplyr)
require(sp)
require(spdep)
require(seg)
require(raster)
require(ggmap)
require(geosphere)

# Set script address as working directory.
#setwd(dirname(parent.frame(2)$ofile))
setwd("D:/Dropbox/Progetto/w/calibration study/Scripts")
#setwd("C:/Users/Thomas/Dropbox/Progetto/w/calibration study/Scripts")
#setwd("C:\\Users\\u838156\\surfdrive\\Shared\\ABM calibration paper\\scripts\\")

load("./cityData/geodata_Rotterdam.RData")



# Moran's I (old implementation)
#
# Once again, we define the function to calculate the index. This time we won't
# take neighborhoods as base unit, but individual agents.
# Like in the previous implementation, the function has an argument "nsize", which
# expresses the threshold (in km) used to define distance-based neighborhood.
# "w" is the dataframe with spatial features that will be passed by the calling command.
MoranI <- function(w, nsize=.1) {

  # We format the coordinates of each agent as spatial points.
  dataworld <- as.data.frame(cbind(w$l, w$j, w$group))
  names(dataworld) <- c("l", "j", "x")
  coords <- SpatialPoints(
    cbind(dataworld$l, dataworld$j),
    proj4string = CRS("+proj=longlat +ellps=WGS84")
  )
  
  # Then, we construct the proximity-neigborhoods.
  col.nb.0.all <- dnearneigh(coords, 0, nsize) #, longlat=TRUE
  
  # We measure the distance between neighbors, and apply a distance decay function.
  # We use the same distance-decay function that we adopt in the ABM.
  distanceDecayWeights <- nbdists(col.nb.0.all, coords) #, longlat=TRUE
  distanceDecayWeights <- lapply(distanceDecayWeights, function(x) 1/x)
  nb <- nb2listw(
    col.nb.0.all,
    glist=distanceDecayWeights,
    style="W",
    zero.policy=TRUE
  )
  print(paste0("Distance: ", nsize, "km"))
  
  # And now we have all we need to write the part of the function
  # that actually calculates Moran's I.
  sp.corg <- moran.test (
    dataworld$x,
    nb,
    zero.policy=TRUE
  )
  return(sp.corg)
}

# We initialize the variables we need for the next step
correlograms <- c()
distances <- c()
wijks <- c()

# And we call the function MoranI for every district and for
# a number of distance values (200m to 1.2km, with steps of 100m).
for(w in districtsList){
  print(paste("Calculating correlogram for district code: ", w))
  focalDistrict <<- subset(world, world$wijk == w)
  #for (m in 2:12){
  for (m in 2:12) {####################################################################
    nsize <- m / 10
    correlogram <- MoranI(focalDistrict, nsize)
    correlograms <- rbind(correlograms, correlogram$estimate)
    distances <- c(distances, nsize)
    wijks <- c(wijks, w)
  }
}

# We put the estimates and parameters into a dataframe, transforming to the right class
# where needed
correlograms <- cbind(correlograms, distances, wijks)
correlograms <- as.data.frame(correlograms)
names(correlograms) <- c("I", "expected", "variance", "distance", "district")
correlograms$I <- as.numeric(levels(correlograms$I))[correlograms$I]
correlograms$expected <- as.numeric(levels(correlograms$expected))[correlograms$expected]
correlograms$variance <- as.numeric(levels(correlograms$variance))[correlograms$variance]
correlograms$distance <- as.numeric(levels(correlograms$distance))[correlograms$distance]
correlograms$sdev <- sqrt(correlograms$variance)
#head(correlograms)

# We plot the correlograms from the dataframe.
# The plot won't be printed on screen, but exported to png.
png(
  filename = "./outputGraphics/MoranI_individualLevel.png",
  width = 800,
  height = 800,
  units = "px"
)
axisRange <- range(c(-0.1, 0.5))
par(mfrow=c(4,3), mar=c(4,4,3,0), oma=c(4,4,0,0))
for(wijk in districtsList){
  cpd <- subset(correlograms, correlograms$district == wijk)
  plot(cpd$distance, cpd$I,
       ylim=c(0, 0.2),
       #ylim=axisRange,
       xlim=c(0.2, 1.2),
       pch=19,
       xlab="",
       ylab="",
       #xlab="Distance (km)",
       #ylab="Moran's I",
       main=wijk
  )
}
mtext("Distance (km)",side=1,line=0,outer=TRUE,cex=1.3)
mtext("Moran's I",side=2,line=0,outer=TRUE,cex=1.3,las=0)
dev.off()



# (3)
# Dissimilarity Index
# 
# We start by rasterizing the simulated city ("world").

# First, we need to define the grids that we will overlay.
# We locate their south-west corner, and determine how big the cells
# should be, and therefore how many we need to cover the city.
minX <- min(world$x_coor)
minY <- min(world$y_coor)
maxX <- max(world$x_coor)
maxY <- max(world$y_coor)
spanX <- maxX - minX
spanY <- maxY - minY
nCells200m_X <- spanX %/% 200
if (spanX %% 200 > 0) {nCells200m_X <- nCells200m_X + 1}
nCells200m_Y <- spanY %/% 200
if (spanY %% 200 > 0) {nCells200m_Y <- nCells200m_Y + 1}
nCells500m_X <- spanX %/% 500
if (spanX %% 500 > 0) {nCells500m_X <- nCells500m_X + 1}
nCells500m_Y <- spanY %/% 500
if (spanY %% 500 > 0) {nCells500m_Y <- nCells500m_Y + 1}

# Now ce calculate the raster that we will later fill in with data.
# We create 2 rasters for two cell sizes: 200*200 and 500*500 meters.
# Note: the documentation on the function "raster()" does no clarify how it deals
# with space spans that are not exact multiples of the cell size. I did a quick 
# test, and it seems that the function trims off some areas, therefore excluding
# some agents. To prevent this, we pass to the function a space span that is the
# minimum exact multiple of the cell size, such that all agents' locations are
# accounted for. Rounding errors can potentially break this fix - because of this,
# after having created the rasters we will test whether they have the correct size.
rastMaxX200m <- nCells200m_X * 200 + minX
rastMaxY200m <- nCells200m_Y * 200 + minY
rastMaxX500m <- nCells500m_X * 500 + minX
rastMaxY500m <- nCells500m_Y * 500 + minY
raster200m <- raster(
  #ncols=nCell200m_X,
  #nrows=nCells200m_Y,
  xmn=minX,
  xmx=rastMaxX200m,
  ymn=minY,
  ymx=rastMaxY200m,
  resolution=200
)
raster500m <- raster(
  #ncols=nCell500m_X,
  #nrows=nCells500m_Y,
  xmn=minX,
  xmx=rastMaxX500m,
  ymn=minY,
  ymx=rastMaxY500m,
  resolution=500
)

# Here we make sure that the rasters have the correct size.
if (
  raster200m@ncols != nCells200m_X |
  raster200m@nrows != nCells200m_Y |
  raster500m@ncols != nCells500m_X |
  raster500m@nrows != nCells500m_Y )
{warning("Error with raster generation.")}


# Next, we overlay raster and simulated city : in other words, for each square in the
# newly created rasters we add data on the group composition of their agents.
agentsCoords <- cbind(world$x_coor, world$y_coor)
grid200m <- rasterize(
  agentsCoords,
  raster200m,
  world$group,
  fun=mean
)
grid500m <- rasterize(
  agentsCoords,
  raster500m,
  world$group,
  fun=mean
)
#plot(grid200m)
#plot(grid500m)

# We convert our raster data to the class "SpatialPolygonsDataFrame".
grid200m <- rasterToPolygons(grid200m)
grid200m@data$ID <- 1:nrow(grid200m)
grid500m <- rasterToPolygons(grid500m)
grid500m@data$ID <- 1:nrow(grid500m)

# Lastly, we assign agents to their cell in the corresponding grid (100m, 200m, 500m).
world$cell100m <- world$location
agents <- SpatialPointsDataFrame(agentsCoords, world)
overlay <- over(agents, grid200m)
world$cell200m <- overlay$ID
overlay <- over(agents, grid500m)
world$cell500m <- overlay$ID


# We have rasterized the simulated city.
# Now we can proceed to calculate the index of dissimilarity across raster cells.
# The library "seg" comes with a function, "dissim()", that's supposed to do just
# that. However, I found it easier to implement it from the scratch than to make
# theirs work. So here we go:
dissimilarityIndex <- function(grid, inputWorld) {
  
  # We calculate the denominators of the dissimilarity index:
  sumWestern <- sum(inputWorld$western)
  sumNonWestern <- sum (inputWorld$nonWestern)
  
  # We force R to interpret "grid" as the name of the variable we pass to the 
  # function as argument. For example, if we request dissimilarityIndex(cell500m, world)
  # the script will interpret "inputWorld$grid" as "world$cell500m".
  inputWorld$grid <- eval(substitute(grid), inputWorld)
  
  # Next, we mesure the composition of each individual locality.
  diss <- c()
  for (i in unique(inputWorld$grid)){
    w <- subset(inputWorld, inputWorld$grid == i)
    localWestern <- sum (w$western)
    localNonWestern <- sum (w$nonWestern)
    diss <- append(
      diss,
      abs((localWestern / sumWestern) - (localNonWestern / sumNonWestern))
    )
  }
  
  # And finally we have all the ingredients to calculate the index:
  return (sum(diss) * 0.5)
}
#dissimilarityIndex(cell500m, world)

# Now we proceed ad usual: we initialize the variables we need, and we calculate
# the index for each raster within every district. This isn't the most elegant way
# to loop trhough the degrees of freedom, but it works.
dissimilarities100 <- c()
dissimilarities200 <- c()
dissimilarities500 <- c()
wijks <- c()
localities <- c()
for(w in districtsList){
  print(paste("Calculating Dissimilaity Index for district code: ", w))
  focalDistrict <<- subset(world, world$wijk == w)
  dissimilarities100 <- c(dissimilarities100, dissimilarityIndex(cell100m, focalDistrict))
  dissimilarities200 <- c(dissimilarities200, dissimilarityIndex(cell200m, focalDistrict))
  dissimilarities500 <- c(dissimilarities500, dissimilarityIndex(cell500m, focalDistrict))
  wijks <- c(wijks, w)
}

# We also calculate the dissimilarity index at the level of the whole city. Note
# that here we use the rasters AND the districts map. This means that we calculate the
# dissimilarity index taking as spatial units individual districts as well as raster cells.
print("Calculating Dissimilaity Index on the whole city")
dissimilarities100 <- c(dissimilarities100, dissimilarityIndex(cell100m, world))
dissimilarities200 <- c(dissimilarities200, dissimilarityIndex(cell200m, world))
dissimilarities500 <- c(dissimilarities500, dissimilarityIndex(cell500m, world))
dissimilarityWholeCity <- dissimilarityIndex(wijk, world)
wijks <- c(wijks, "Whole city")

# We stitch together all the data generated into a single dataframe:
acrossDistricts <- c()
acrossDistricts[length(wijks)] <- dissimilarityWholeCity
dissimilarityIndices <- as.data.frame( cbind(
  wijks,
  dissimilarities100,
  dissimilarities200,
  dissimilarities500,
  acrossDistricts
))
dissimilarityIndices


# Finally, we proceed to plot the dissimilarity index. 
# Note: the plot will be printed on file, and not on screen.
xAxisLabels <- c("Raster 100m", "Raster 200m", "Raster 500m", "Across districts")
png(
  filename = "./outputGraphics/DissimilarityIndex_local.png",
  width = 1200,
  height = 2000,
  units = "px",
  res=150
)
#axisRange <- range(c(-0.1, 0.5))
par(mfrow=c(5,3), mar=c(7,4,3,0), oma=c(4,4,0,0))
for(i in dissimilarityIndices$wijks){
  cpd <- subset(dissimilarityIndices, dissimilarityIndices$wijks == i)
  cpd$dissimilarities100 <- as.numeric(levels(cpd$dissimilarities100))[cpd$dissimilarities100]
  cpd$dissimilarities200 <- as.numeric(levels(cpd$dissimilarities200))[cpd$dissimilarities200]
  cpd$dissimilarities500 <- as.numeric(levels(cpd$dissimilarities500))[cpd$dissimilarities500]
  cpd$acrossDistricts <- as.numeric(levels(cpd$acrossDistricts))[cpd$acrossDistricts]
  x <- c(1,2,3,4)
  y <- c(
    cpd$dissimilarities100,
    cpd$dissimilarities200,
    cpd$dissimilarities500,
    cpd$acrossDistricts
  )
  cpd$dissimilarities100
  #print(length(dissimilarities100))
  plot(
    x,
    y,
    xlim=c(0,5),
    ylim=c(0,1),
    pch=19,
    xlab="",
    ylab="",
    xaxt ='n',
    main=i
  )
  axis(
    1,
    at=1:4,
    labels=FALSE
  )
  text(
    seq(0.5, 3.5, by=1),
    par("usr")[3] - 0.3,
    labels = xAxisLabels,
    srt = 65,
    pos = 1,
    xpd = TRUE
  )
}
#mtext("Distance (km)",side=1,line=0,outer=TRUE,cex=1.3)
mtext("Dissimilarity Index",side=2,line=0,outer=TRUE,cex=1.3,las=0)
dev.off()





# Plot squares on map
#
# Here we do two things. We:
#   (1) plot the squares' centroids on the basemap of Rotterdam.
#   (2) plot agents (the ones that we used to calculate the micro-level measuers
#       of segregation) on the basemap of Rotterdam.


# (1)
# We first download the basemap, focused on the mean latitute and longitude of Rotterdam
exportMapPlot <- FALSE # If TRUE, saves the plot in a png file in "./output graphics/"
                       # Else, it just prints the plots on screen.
allSquareMapPoints <- cbind.data.frame(cbs100_rot$l, cbs100_rot$j)
meanLon <- mean(allSquareMapPoints[,1])
meanLat <- mean(allSquareMapPoints[,2])
areaToPlot <- get_map(
  location = c(lon = meanLon, lat = meanLat),
  #source = "google",
  #maptype = "satellite",
  #maptype = "terrain-background",
  source = "stamen",
  #maptype = "toner",
  maptype = "toner-background",
  #maptype = "terrain",
  color = "bw", 
  zoom = 11, # default for city level = 12
  scale = 2
)

# We then define a function to plot square centroids on our basemap
plotDistrict <- function(){
  ggmap(areaToPlot) +
    geom_point(
      data = allSquareMapPoints,
      aes(x = allSquareMapPoints[,1], y = allSquareMapPoints[,2], alpha = 0.3),
      color ="darkorange",
      size = 1, #0.05,
      shape = 16 #20
    ) +
    geom_point(
      data = squaresMapPoints,
      aes(x = squaresMapPoints[,1], y = squaresMapPoints[,2], alpha = 1),
      color ="orangered",
      size = 1, #0.05,
     shape = 16 #20
    ) +
    guides(fill=FALSE, alpha=FALSE, size=FALSE) +
    labs(x = "Longitude", y = "Latitude") +
    ggtitle(paste("Rotterdam - district code:", wijk))
}

# Then, we loop throught all districts: for each of them, we plot the
# square centroids that belong to the given district. This loop produces 
# as many plots as there are districts.
for(wijk in districtsList){
  focalDistrict <- subset(cbs100_rot, cbs100_rot$WK_CODE == wijk)
  squaresMapPoints <- cbind.data.frame(focalDistrict$l, focalDistrict$j)
  
  # If export is toggled on, this saves the plot in png format. Prints on screen otherwise.
  if (exportMapPlot == TRUE){
    png(
      filename = paste0("./output graphics/Rotterdam_wijk_", wijk, ".png"),
      width = 500,
      height = 500,
      units = "px"
    )
  } 
  print(plotDistrict())
  if (exportMapPlot == TRUE){dev.off()}
}
#squaresMapPoints <- cbind.data.frame(cbs100_rot$l, cbs100_rot$j)
#plotDistrict()



# (2)
# Here we plot agents. We do the same as before:
# We first download the basemap:
plotAgents <- function(){
  areaToPlot <- get_map(
    location = c(lon = mean(world$l), lat = mean(world$j)),
    #source = "google",
    #maptype = "satellite",
    #maptype = "terrain-background",
    source = "stamen",
    #maptype = "toner",
    maptype = "toner-background",
    #maptype = "terrain",
    color = "bw", 
    zoom = 12, # default for city level = 12
    scale = 2
  )
  
  # Then, we plot the basemap along with the agents.
  ggmap(areaToPlot) +
    geom_point(
      data = world,
      
      # (note how we color the two demographic groups differently)
      aes(x = world$l, y = world$j, alpha = 0.1, color = as.factor(world$group)),
      size = 0.05, #0.05,
      shape = 16 #20
    ) +  
    #scale_colour_manual(world$group=c("-1" = "red", "1" = "blue")) +
    guides(fill=FALSE, alpha=FALSE, size=FALSE)
}
plotAgents()



##########################   *** working ***   #################################
#Preparing the labels for the plot by adding the district names to the dataset cit.

#load("./cityData/geodata_Rotterdam.RData")
cit <- subset(nedb,nedb$WK_CODE %in% districtsList)
cit@data$WK_CODE <- droplevels(cit@data$WK_CODE)
cit <- spTransform(cit,CRS("+proj=longlat +ellps=WGS84"))

cit@data$WK_NAAM <- NA
cit@data$WK_NAAM <- as.factor(cit@data$WK_NAAM)
levels(cit@data$WK_NAAM) <- levels(citySummary$district)

for (i in 1:length(cit@data$WK_CODE)){
  l <- subset(citySummary, citySummary$WK_CODE == cit@data$WK_CODE[i])
  l <- l$district
  cit@data$WK_NAAM[i] <- l
}
#cit@data$WK_NAAM


# Joining neighborhood (buurten) polygon into a district polygon
require(maptools)
require(rgeos)
cit <- unionSpatialPolygons(SpP = cit, IDs = cit$WK_NAAM)
labCoords <- as.data.frame(gCentroid(cit, byid = TRUE))

# Here I reposition some of the labels by hand,
# so that they do not overlap in the final plot
labCoords[2,2] <- labCoords[2,2] - 0.005 # Delfshaven
labCoords[4,2] <- labCoords[4,2] + 0.005 # Hillegersberg-Schiebroek
labCoords[7,2] <- labCoords[7,2] - 0.002 # Kralingen-Crooswijk


# Downloading basemap
basemap <- get_map(
  location = "Rotterdam",
  #source = "stamen",
  source = "google",
  #maptype = "toner-background",
  maptype = "satellite",
  #color = "bw", 
  zoom = 11,
  scale = 2
)

# Plotting districts
ggmap(basemap, extent = "device") +
  geom_polygon(
    data = cit,
    aes(x=long, y=lat, group=group),
    color="dark orange",
    size = 1.2,
    fill="orange",
    alpha=0.1
  ) +
  geom_label(
    data = labCoords,
    aes(x=x, y=y, label=row.names(labCoords)),
    color = "medium purple 4",
    alpha = 0.7,
    fontface = "bold",
    lineheight = 0
    #label.size = 0
  )


## Plot opinion distributions
#
#
# Here we plot the two beta distrubitons that we use as initial
# opinion distribution.
x <- (c(1:1000) / 1000)
# These are our density functions: 
beta <- dbeta(x, 3, 3, ncp = 0)
unif <- dbeta(x, 1, 1, ncp = 0)
groupBias1 <- dbeta(x, 3, 3.5, ncp = 0)
groupBias2 <- dbeta(x, 3.5, 3, ncp = 0) 

par(mfrow=c(1,3), mar = c(2,1,1,1))
plot((x * 2)-1,unif,
     xlim=c(-1,1),
     yaxt='n', xaxt='n', #cex.axis=2,
     bty="n",
     type="l",
     xlab="",
     ylab="",
     col="black"#"purple"#, main="α=β=3"
)
axis(side=1, at=c(-1,0,1), lwd=3)
title(expression(paste(alpha,"=",beta,"=1")), col.main="black")#"purple")
#
plot((x * 2)-1,beta,
     xlim=c(-1,1),
     yaxt='n', xaxt='n', #cex.axis=2,
     bty="n",
     type="l",
     xlab="",
     ylab="",
     col="black"#"purple"#, main="α=β=3"
)
axis(side=1, at=c(-1,0,1), lwd=3)
title(expression(paste(alpha,"=",beta,"=3")), col.main="black")#"purple")
#
plot((x * 2)-1,groupBias1,
     xlim=c(-1,1),
     yaxt='n', xaxt='n', #cex.axis=2,
     bty="n",
     type="l",
     xlab="",
     ylab="",
     lty=3,
     col="black"#"blue"#, main="α=β=3"
)
points((x * 2)-1,groupBias2, type="l", lty=5, col="black")#"red")
axis(side=1, at=c(-1,0,1), lwd=3)
title(expression(paste(alpha,"=3, ",beta,"=3.5")), col.main="black", adj=0)
title(expression(paste(alpha,"=3.5, ",beta,"=3")), col.main="black", adj=1)

rm(unif, beta, groupBias1, groupBias2)


# Here we do a similar thing, by plotting histograms of the sampled distribution
hst <- function (data, yt, t){
  ggplot(
    data=data,
    aes(x=data[,1]),
    stat(count())
  ) +
    #geom_histogram(binwidth=0.1, colour="black", alpha=0) +
    geom_density(color="white", fill="black", position = "stack") +
    ylim(0,0.77) +
    labs(
      x = "Attitude scale",
      y = yt,
      title = t
    ) +
    theme(axis.line.y=element_blank(),axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.ticks.x=element_line(size=1),
          axis.ticks.length=unit(0.3, "cm"),
          axis.line.x=element_line(size=1.5),
          axis.text.x=element_text(size=20),
          #axis.title.x=element_blank(),
          #axis.title =  element_text(angle = 0, hjust = 0),
          axis.title.y = element_text(size=25),
          axis.title.x = element_text(size = 25),
          title = element_text(size=25),
          #axis.title.x.top =  element_text(hjust = 0),
          #axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(), plot.background=element_blank(),
          panel.border = element_blank()
    )
}
#hst(as.data.frame(rbeta(10000, 3, 3, ncp = 0) * 2 - 1), t="α=β=3")
#hst(as.data.frame(rbeta(10000, 1, 1, ncp = 0) * 2 - 1), t="α=β=1")

png("outputGraphics/initial opinion distribution.png", 
    width = 600, height = 300,
    units = "px", pointsize = 10, 
    res = NA)
grid.arrange(
  hst(as.data.frame(rbeta(1000000, 2, 2, ncp = 0) * 2 - 1),yt="Frequency", t="α=β=3"),
  hst(as.data.frame(rbeta(1000000, 1, 1, ncp = 0) * 2 - 1), yt="", t="α=β=1"),
  nrow = 1
)
dev.off()


## Plot stylized examples of outcome measures
#
#
# We start by creating the stylized examples of districts:
require(ggplot2)
set.seed(1234567890)
crtFakeDistrict <- function(n,relativeGroupSize=0.5, bins=8){
  x <- sample.int(bins, size = n, replace = TRUE)
  y <- sample.int(bins, size = n, replace = TRUE)
  d <- as.data.frame(cbind(x,y))
  d <- unique(d)
  if (nrow(d) != n) {
    while(nrow(d) < n) {
      l <- n - nrow(d)
      x <- sample.int(bins, size = l, replace = TRUE)
      y <- sample.int(bins, size = l, replace = TRUE)
      d <- rbind(d, as.data.frame(cbind(x,y)))
      d <- unique(d)
    }
  }
  d$g <- rbinom(n, 1, relativeGroupSize)
  return(d)
}



# This is the function to plot the agents in the fake district
plt  <- function(data, marg=margin(0,0,0,0, "pt")){
  ggplot(
    data,
    aes(
      x = x,
      y = y,
      fill = o
    )
  ) +
    geom_point(size = 15, shape = 21) +
    scale_fill_gradient(
      low = "black",
      high = "white",
      limits = c(-1, 1)
    ) +
    xlim(0.7, 8.3) + ylim(0.7,8.3) +
    #scale_x_continuous(position = "top") +
    labs(
      #x = x_lab1,
      fill = "Attitude scale"
    ) +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          plot.margin = marg,
          #axis.title.x=element_blank(),
          #axis.title =  element_text(angle = 0, hjust = 0),
          axis.title = element_text(size = 25),
          axis.title.x.top =  element_text(hjust = 0),
          #axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(), plot.background=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    )
}

hst <- function(d, y=NULL, marg=margin(0,0,0,0, "pt")) {
  if (is.null(y)){
    ggplot(data=d, aes(d$o)) +
      geom_histogram(
        breaks=c(((0:7)/7*2)-1),
        col="white",
        fill="black"
      ) +
      scale_x_continuous(
        breaks=c(-1,0,1),
        limits=c(-1,1)
      ) +
      scale_y_continuous(breaks=NULL) +
      theme(
        plot.margin = marg,
        axis.line.x = element_line(size = 1, linetype = "solid"),
        axis.text.x= element_text(size=20),
        axis.title = element_text(size = 25),
        legend.position="none",
        panel.background=element_blank(), plot.background=element_blank()
      )
  } else {
    ggplot(data=d, aes(d$o, fill=as.factor(d$g))) +
      geom_histogram(
        breaks=c(((0:7)/7*2)-1),
        col="white",
        position = "stack"
        #fill="black"
      ) +
      scale_fill_manual(values=c("#ff8800", "#4c00ff")) +
      scale_x_continuous(
        breaks=c(-1,0,1),
        limits=c(-1,1)
      ) +
      scale_y_continuous(breaks=NULL) +
      theme(
        plot.margin = marg,
        axis.line.x = element_line(size = 1, linetype = "solid"),
        axis.text.x= element_text(size=20),
        axis.title = element_text(size = 25),
        legend.position="none",
        panel.background=element_blank(), plot.background=element_blank()
      )
  }
}
#hst(d=d9, y=TRUE)
#hst(d)
#hst(d, marg=margin(0,20,0,0,"pt"))

# Here's a function to artificially increase clustering.
findNewSpot <- function (d, i, q){
  old_x <- d$x[i]
  old_y <- d$y[i]
  #print(paste(old_x, old_y, t))
  it <- 1
  repeat{
    d$x[i] <- sample.int(4, size = 1, replace = TRUE) + q
    d$x[i] <- sample.int(4, size = 1, replace = TRUE) + q
    #print(paste("tweaked - it", it))
    if(nrow(unique(d)) == 35){break}
    if(it == 5){
      d$x[i] <- old_x
      d$y[i] <- old_y
      #print("no change")
      break
    } else {it <- it + 1}
  }
}

# Now we can create the 4 fake districts based on the dimensions we need to show
# (clustering and polarization)
d <- crtFakeDistrict(35)
d$o <- rbeta(35, 0.02, 0.02, ncp = 0) * 2 - 1 # Y pol, N clustering
plt(d)

d2 <- crtFakeDistrict(35)
for(i in 1: 35){                            # Y pol, Y clustering
  t <- (d2$x[i] + d2$y[i]) / 16
  d2$o[i] <- rbeta(1, 0.1 * t, 0.1 * (1 - t), ncp=0) * 2 - 1
  if (d2$o[i] > 0 & t < 0.7){findNewSpot(d2, i, 4)}
  if (d2$o[i] < 0 & t > 0.3){findNewSpot(d2, i, 0)}
}
plt(d2)

d3 <- crtFakeDistrict(35)
d3$o <- rbeta(35, 5, 5, ncp = 0) * 2 - 1    # N pol, N clustering
plt(d3)

d4 <- crtFakeDistrict(35)
for(i in 1: 35){                            # N pol, Y clustering
  t <- (d4$x[i] + d4$y[i]) / 16
  #ifelse(
  #  t<0.5,
  #  d4$o[i] <- rbeta(1, 1.5, 1, ncp=0)- 1,
  #  d4$o[i] <- rbeta(1, 1, 1.5, ncp=0)
  #)
  d4$o[i] <- rbeta(1, 10 * t, 10 * (1 - t), ncp=0) * 2 - 1
  if (d4$o[i] > 0 & t < 0.7){findNewSpot(d4, i, 4)}
  if (d4$o[i] < 0 & t > 0.3){findNewSpot(d4, i, 0)}
}
plt(d4)







# Polarization
#
#
require(gtable)
require(grid)
require(gridExtra)
png("outputGraphics/polarized district.png", 
    width = 850, height = 800,
    units = "px", pointsize = 5, 
    res = NA)
grid.arrange(
  hst(d3, marg=margin(0,35,0,0,"pt"))+ labs(x="", y="Attitude frequency"),
  hst(d, marg=margin(0,0,0,35,"pt"))+ labs(x="", y=""),
  plt(d3, marg=margin(0,35,10,0,"pt"))+ labs(x="Weak polarization", y=""),
  plt(d, marg=margin(0,0,10,35,"pt"))+ labs(x="Strong polarization", y=""),
  nrow=2#,
  #plt(d2)+ labs(x="", y=""),
  #plt(d4)+ labs(x="Strong polarization", y=""),
  #top="Clustering",
  #left="Polarization"
)
dev.off()



# Polarization v clustering
#
# 
require(gtable)
require(grid)
require(gridExtra)
png("./outputGraphics/polarization v clustering.png", 
    width = 800, height = 800,
    units = "px", pointsize = 5, 
    res = NA)
grid.arrange(
  plt(d)+ labs(x="", y="Strong polarization"),
  plt(d2)+ labs(x="", y=""),
  plt(d3)+ labs(x="Weak clustering", y="Weak polarization"),
  plt(d4)+ labs(x="Strong clustering", y=""),
  nrow=2#,
  #top="Clustering",
  #left="Polarization"
)
dev.off()



# Clustering: local v global:
#
#
# Here we do the same as before: we start by generating two fake district that
# differ on the clustering dimension:

d5 <- crtFakeDistrict(35) # global clustering
d5$o <- rbeta(35, 0.1, 0.1, ncp = 0) * 2 - 1
for(i in 1:35){
  if(d5$x[i] < 5){d5$o[i] <- rbeta(1, 0.5, 5, ncp = 0) -1}
  if(d5$x[i] > 5){d5$o[i] <- rbeta(1, 5, 0.5, ncp = 0)}
}
plt(d5)

d6 <- crtFakeDistrict(35) # local clustering
d6$o <- rbeta(35, 0.01, 0.01, ncp = 0) * 2 - 1
for(i in 1:35){
  if(d6$y[i] < 5){
    d6$o[i] <- rbeta(1, 5, 0.2, ncp = 0)
    if(d6$x[i] < 5) {
      d6$o[i] <- rbeta(1, 0.2, 5, ncp = 0) -1
    }
  }
  if(d6$y[i] > 5){
    d6$o[i] <- rbeta(1, 5, 0.2, ncp = 0)
    if(d6$x[i] > 5){
      d6$o[i] <- rbeta(1, 0.2, 5, ncp = 0) -1
    }
  }
}
plt(d6)


# And then we wrap them in an output plot.
png("outputGraphics/local v global clustering.png", 
    width = 1000, height = 400,
    units = "px", pointsize = 5, 
    res = NA)
grid.arrange(
  plt(d, marg=margin(0,5,10,0,"pt")) +
    labs(x="Polarization\nwithout clustering", y=NULL),
  plt(d6, marg=margin(0,5,10,5,"pt")) + 
    labs(x="Polarization and\nlocal clustering", y=NULL),
  plt(d5, marg=margin(0,0,10,5,"pt")) + 
    labs(x="Polarization,\nlocal and global clustering", y=NULL),
  nrow=1#,
  #top="Clustering",
  #left="Polarization"
)
dev.off()


# Alignment:
#
#
# We define a function to plot the fake district which is also able to plot the
# group membership of agents.
plt2  <- function(data, marg=margin(0,0,0,0,"pt")){
  ggplot(
    data,
    aes(
      x = x,
      y = y,
      fill = o,
      col = g
    )
  ) +
    geom_point(size = 14, shape = 21,  stroke = 3) +
    scale_fill_gradient(
      low = "black",
      high = "white",
      limits = c(-1, 1)
    ) +
    scale_color_gradient(
      low="#ff8800",
      high="#4c00ff"#"black"
    ) +
    xlim(min(data$x)-0.3, max(data$x)+0.3) + ylim(min(data$y)-0.3, max(data$y)+0.3) +
    #scale_x_continuous(position = "top") +
    labs(
      #x = x_lab1,
      fill = "Attitude scale"
    ) +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          plot.margin = marg,
          axis.title = element_text(size = 25),
          axis.title.x.top =  element_text(hjust = 0),
          legend.position="none",
          panel.background=element_blank(), plot.background=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    )
}
plt2(d9)



# With this dummy district we show that high local levels of exposure do not 
# necessarily lead to high global averages.
require(ggplot2)
require(gridExtra)
test <- crtFakeDistrict(16, relativeGroupSize=0, bins=4)
test$o <- 1
test$g<-0#[sample(1:nrow(test), 2)] <- 1 
proxmat <- exp(-distmat/10)
for (i in 1:length(distmat[1,])){
  proxmat[i,] <- proxmat[i,] / sum(proxmat[i,])
}

plt2(test)

test$index <- c(1:nrow(test))
test$inw2014 <- 1
test$nnwal2014_trunc <- test$group <- test$g
distmat <- as.matrix(dist(test, diag=TRUE, upper=TRUE))
test$exposureOutgroup <- calcExposure(test, test, proxmat)$V3

print(paste(
  "Test 1 - Population average outgroup exposure:",
  round(mean(test$exposureOutgroup), digits=2), ";",
  "Average outgroup exposure of blues:",
  round(mean(subset(test, test$g==1)$exposureOutgroup), digits=2)
))

test2<-test
test2$g <- 0
test2[c(1,11),]$g <- 1
#test2$g <- 0
plt2(test2)
test2$index <- c(1:nrow(test2))
test2$inw2014 <- 1
test2$nnwal2014_trunc <- test2$group <- test2$g
distmat <- as.matrix(dist(test2, diag=TRUE, upper=TRUE))
test2$exposureOutgroup <- calcExposure(test2, test2, proxmat)$V3

print(paste(
  "Test 2 - Population average outgroup exposure:",
  round(mean(test2$exposureOutgroup), digits=2), ";",
  "Average outgroup exposure of blues:",
  round(mean(subset(test2, test2$g==1)$exposureOutgroup), digits=2)
))

# Exposure
png("outputGraphics/exposure.png", 
    width = 800, height = 400,
    units = "px", pointsize = 5, 
    res = NA)
grid.arrange(
  plt2(test, marg=margin(0,30,10,0,"pt")) + 
    labs(x="",y=""),
    #labs(x="Lower,\nBlues have higher outgr. exp.", y=""),
  plt2(test2, marg=margin(0,30,10,0,"pt")) + 
    labs(x="",y=""),
    #labs(x="Higher average outgr. exp.,\nBlues have lower outgr. exp.", y=""),
  nrow=1
)
dev.off()















# Next we generate the needed fake districts, in this case 4.

d7 <- crtFakeDistrict(35) # Y clustering, N alignment
d7$o <- rbeta(35, 0.01, 0.01, ncp = 0) * 2 - 1
for(i in 1:35){
  if(d7$y[i] < 5){
    d7$o[i] <- rbeta(1, 5, 0.2, ncp = 0)
    if(d7$x[i] < 5) {
      d7$o[i] <- rbeta(1, 0.2, 5, ncp = 0) -1
    }
  }
  if(d7$y[i] > 5){
    d7$o[i] <- rbeta(1, 5, 0.2, ncp = 0)
    if(d7$x[i] > 5){
      d7$o[i] <- rbeta(1, 0.2, 5, ncp = 0) -1
    }
  }
}
plt2(d7)

d8 <- crtFakeDistrict(35) # N clustering, N alignment
d8$o <- rbeta(35, 0.01, 0.01, ncp = 0) * 2 - 1
plt2(d8)

d9 <- crtFakeDistrict(35) # Y clustering, Y alignment
d9$o <- rbeta(35, 0.01, 0.01, ncp = 0) * 2 - 1
for(i in 1:35){
  if(d9$y[i] < 5){
    d9$o[i] <- rbeta(1, 5, 0.2, ncp = 0)
    if(d9$x[i] < 5) {
      d9$o[i] <- rbeta(1, 0.2, 5, ncp = 0) -1
    }
  }
  if(d9$y[i] > 5){
    d9$o[i] <- rbeta(1, 5, 0.2, ncp = 0)
    if(d9$x[i] > 5){
      d9$o[i] <- rbeta(1, 0.2, 5, ncp = 0) -1
    }
  }
  p <- abs(d9$o[i])
  if (d9$o[i] > 0){
    if (d9$g[i] <- rbinom(1, size=1, prob=p)){d9$g[i] <- 1} else {d9$g[i] = -1}
  } else {
    if (d9$g[i] <- rbinom(1, size=1, prob=p)){d9$g[i] <- -1} else {d9$g[i] = 1}
  }
}
plt2(d9)

d0 <- crtFakeDistrict(35) # N clustering, Y alignment
d0$o <- rbeta(35, 0.01, 0.01, ncp = 0) * 2 - 1
for(i in 1:35){
  if (d0$o[i] > 0){
    if (d0$g[i] <- rbinom(1, size=1, prob=p)){d0$g[i] <- 1} else {d0$g[i] = -1}
  } else {
    if (d0$g[i] <- rbinom(1, size=1, prob=p)){d0$g[i] <- -1} else {d0$g[i] = 1}
  }
}
plt2(d0)


d10 <- crtFakeDistrict(50) # N alignment
#d10$o <- rbeta(35, 0.15, 0.15, ncp = 0) * 2 - 1
d10$g <- -1
for(i in 1:50){
  if(d10$y[i] > 4 & d10$x[i] > 4){
    d10$g[i] <- rbinom(1, size=1, prob=0.9)
  } else {
    if(d10$y[i] < 5 & d10$x[i] < 5){
      d10$g[i] <- rbinom(1, size=1, prob=0.9)
    } else {
      d10$g[i] <- rbinom(1, size=1, prob=0.1)
    }
  }
  d10$o[i] <- rbeta(1, 0.15, 0.15, ncp = 0) * 2 - 1
}
plt2(d10)
hst(d10, y=TRUE)


d11 <- crtFakeDistrict(50) # Y local alignment
#d10$o <- rbeta(35, 0.15, 0.15, ncp = 0) * 2 - 1
d11$g <- -1
for(i in 1:50){
  if(d11$y[i] > 4 & d11$x[i] > 4){
    d11$g[i] <- rbinom(1, size=1, prob=0.9)
  } else {
    if(d11$y[i] < 5 & d11$x[i] < 5){
      d11$g[i] <- rbinom(1, size=1, prob=0.9)
    } else {
      d11$g[i] <- rbinom(1, size=1, prob=0.1)
    }
  }
  ifelse(
    d11$y[i] < 5,
    d11$o[i] <- rbeta(1, 1, 0.2, ncp = 0) * 2 - 1,
    d11$o[i] <- rbeta(1, 0.2, 1, ncp = 0) * 2 - 1
  )
}
plt2(d11)
hst(d11, y=TRUE)


d12 <- crtFakeDistrict(50) # Y global and local alignment
#d10$o <- rbeta(35, 0.15, 0.15, ncp = 0) * 2 - 1
d12$g <- -1
for(i in 1:50){
  if(d12$y[i] > 4 & d12$x[i] > 4){
    d12$g[i] <- rbinom(1, size=1, prob=0.9)
  } else {
    if(d12$y[i] < 5 & d12$x[i] < 5){
      d12$g[i] <- rbinom(1, size=1, prob=0.9)
    } else {
      d12$g[i] <- rbinom(1, size=1, prob=0.1)
    }
  }
  ifelse(
    d12$g[i] == 1,
    d12$o[i] <- rbeta(1, 1, 0.2, ncp = 0) * 2 - 1,
    d12$o[i] <- rbeta(1, 0.2, 1, ncp = 0) * 2 - 1
  )
}
plt2(d12)
hst(d12, y=TRUE)



# Alignment
png("outputGraphics/Alignment.png", 
    width = 1100, height = 720,
    units = "px", pointsize = 5, 
    res = NA)
grid.arrange(
  hst(d=d10, y=TRUE, marg=margin(0,30,0,0,"pt")) + labs(x="", y="Attitude frequency"),
  hst(d=d11, y=TRUE, marg=margin(0,30,0,00,"pt")) + labs(x="", y=""),
  hst(d=d12, y=TRUE, marg=margin(0,30,0,00,"pt")) + labs(x="", y=""),
  plt2(d10, marg=margin(0,30,10,0,"pt")) + 
    labs(x="Polarization\nwithout alignment", y=""),
  plt2(d11, marg=margin(0,30,10,0,"pt")) + 
    labs(x="Polarization and\nlocal alignment", y=""),
  plt2(d12, marg=margin(0,30,10,0,"pt")) + 
    labs(x="Polarization,\nlocal and global alignment", y=""),
  nrow=2
)
dev.off()


# Clustering v alignment
png("outputGraphics/clustering v alignment.png", 
    width = 800, height = 800,
    units = "px", pointsize = 5, 
    res = NA)
grid.arrange(
  plt2(d7) + labs(x="", y="Strong (local) attitude clustering"),
  plt2(d9) + labs(x="", y=""),
  plt2(d8) + labs(x="Weak alignment", y="Weak attitude clustering"),
  plt2(d0) + labs(x="Strong alignment", y=""),
  nrow=2
)
dev.off()




# ************ This is too slow to be doable *****************************
# At this point we have all we need to calculate the local measures of
# segregation. We start by calculating agent's exposure to other ingroup or
# ougroup agents.
# We do that by calculing two matrices, one with the proximity between ingroup
# agents, and one between outgroup agents.
prob_xx <- matrix(NA, nrow=populationSize , ncol=populationSize)
prob_xy <- matrix(NA, nrow=populationSize , ncol=populationSize)

# To fill in the top triangle of the matrices, we loop through them:
for (i in 1:populationSize){
  for (j in i:populationSize){
    
    # For agents who live in the same cell, we want to avoid to assume their
    # distance is 0, because their proximity would be maximal. Instead, we
    # assume that their proximity is the average distance between all points
    # in a square sized 100*100 meters. That would be about 52.14m
    if (world$index[i]==world$index[j]) {
      prob_xx[i,j] <- exp(-52.140543316/100) * (world$group[i] == world$group[j])
      prob_xy[i,j] <- exp(-52.140543316/100) * (world$group[i] != world$group[j])
    } else {
      
      # To find the proximity between agents from different cells, we resort to
      # the matrix proxmat, which contains the proximities between all cells.
      prob_xx[i,j] <- proxmat[world$index[i],world$index[j]] *
        (world$group[i] == world$group[j])
      prob_xy[i,j] <- proxmat[world$index[i],world$index[j]] *
        (world$group[i] != world$group[j])
    }
  }
}

# We complete the matrix by erasing the diagonal and mirroring the triangle.
diag(prob_xx) <- diag(prob_xy) <- NA
prob_xx[lower.tri(prob_xx)] <- (t(prob_xx))[lower.tri(prob_xx)]
prob_xy[lower.tri(prob_xy)] <- (t(prob_xy))[lower.tri(prob_xy)]
#printMat(prob_xx)

# Next, we find the normalizing constant and we normalize the matrices.
sums <- rowSums(prob_xx, na.rm=T) + rowSums(prob_xy, na.rm=T)
prob_xx <- prob_xx / matrix(sums, nrow=length(sums), ncol=length(sums)) 
prob_xy <- prob_xy / matrix(sums, nrow=length(sums), ncol=length(sums)) 

# Lastly, we add the agent-level ingroup- and outgroup-exposure to the world
# dataframe.
exposureIngroup <- rowSums(prob_xx, na.rm=T)
exposureOutgroup <- rowSums(prob_xy, na.rm=T)
world <- cbind(world, exposureIngroup, exposureOutgroup)
# ************ This is too slow to be doable *****************************



# Global Moran I (old method.) - The function we currently use is in section 3.
#
# First we define the function, which takes as input the argument "nsize": this
# is the threshold (in kilometers) used to defined the distance-based neighborhood.
# "w" is the district information that will be passed by the calling command.
MoranI <- function(w, nsize=.1) {
  
  # We format the district squares (centroids) as spatial points
  dataworld <- as.data.frame(cbind(w$l, w$j, w$nnwal2014))
  names(dataworld) <- c("l", "j", "x")
  coords <- SpatialPoints(
    cbind(dataworld$l, dataworld$j),
    proj4string = CRS("+proj=longlat +ellps=WGS84")
  )
  
  # Then, we find the neighbors of each square...
  col.nb.0.all <- dnearneigh(coords, 0, nsize) #, longlat=TRUE
  print(paste0("Distance: ", nsize, "km"))
  
  # ... and calculate the correlogram.
  # This function isn't optimal, as for some reason "order = 1" 
  # is not a valid argument, and every time we execute this command
  # we get the estimates for an unnecessary second lag (order=2).
  # But it's not a big deal: we will just ignore the estimates for
  # the second lag.
  sp.corg <- sp.correlogram(
    col.nb.0.all, dataworld$x,
    order=2, # *bug*: won't accept the default argument: order=1
    method="I",
    randomisation=FALSE,
    zero.policy=TRUE
  )
  return(sp.corg)
}

# We initialize the variables we need for the next step
correlograms <- c()
distances <- c()
wijks <- c()

# And we call the function MoranI for every district and for
# a number of distance values (200m to 1.2km, with steps of 100m).
for(wijk in districtsList){
  print(paste("Calculating correlogram for district code: ", wijk))
  focalDistrict <- subset(cbs100_rot, cbs100_rot$WK_CODE == wijk)
  for (m in 2:12){
    nsize <- m / 10
    correlogram <- MoranI(focalDistrict, nsize)
    correlograms <- rbind(correlograms, correlogram$res[1,])
    distances <- c(distances, nsize)
    wijks <- c(wijks, wijk)
  }
}

# We put the estimates and parameters into a dataframe, transforming to the right class
# where needed
correlograms <- cbind(correlograms, distances, wijks)
correlograms <- as.data.frame(correlograms)
names(correlograms) <- c("I", "expected", "variance", "distance", "district")
correlograms$I <- as.numeric(levels(correlograms$I))[correlograms$I]
correlograms$expected <- as.numeric(levels(correlograms$expected))[correlograms$expected]
correlograms$variance <- as.numeric(levels(correlograms$variance))[correlograms$variance]
correlograms$distance <- as.numeric(levels(correlograms$distance))[correlograms$distance]
correlograms$sdev <- sqrt(correlograms$variance)
#head(correlograms)

# We plot the correlograms from the dataframe.
# The plot won't be printed on screen, but exported to png.
png(
  filename = "./output graphics/MoranI.png",
  width = 800,
  height = 800,
  units = "px"
)
axisRange <- range(c(-0.1, 0.5))
par(mfrow=c(4,3))
for(wijk in districtsList){
  cpd <- subset(correlograms, correlograms$district == wijk)
  plot(cpd$distance, cpd$I,
       ylim=axisRange,
       xlim=c(0.2, 1.3),
       pch=19,
       xlab="Distance",
       ylab="Moran's I",
       main=wijk
  )
  arrows(
    x0 = cpd$distance,
    y0 = cpd$I-cpd$sdev,
    x1 = cpd$distance,
    y1 = cpd$I+cpd$sdev,
    code = 3,
    length=0.04,
    angle=90
  )
}
dev.off()




