rm (list = ls( ))
library(gridExtra)
library(ggplot2)
library(reshape2)
library(scales)


# Figure 1 _____________________________________________________________________
# Example districts with higher peaks of outgroup exposure vs higher average
# outgroup exposure.
#
# We want to create a bunch of points to position on a mock-up district.
d <- expand.grid(
  x = 1:4,
  y = 1:3,
  col = 0
)
d1 <- d2 <- d
d1$panel <- "district A"
d2$panel <- "district B"

# Marking two corner points in d1
d1$col[which((d1$x == 1 & d1$y == 1) | (d1$x == 4 & d1$y == 3))] <- 1

# Marking the middle two points in d2
d2$col[which(d2$x %in% 2:3 & d2$y == 2)] <- 1

d <- rbind(d1,d2)
d$col <- as.factor(d$col)


png(
  filename = "./outputGraphics/figure 1 - outgroup exposure example.png",
  width = 1000, height = 500, res = 300, units = "px"
)
ggplot (d, aes(x = x, y = y, fill = col)) +
  geom_point(shape = 21, size = 9, color = "#383838", stroke = 1.1) +
  facet_grid(cols = vars(panel)) +
  scale_x_continuous(expand = c(0.15, 0.15)) +
  scale_y_continuous(expand = c(0.15, 0.15)) +
  scale_fill_manual(values = c("0" = "#f2f2f2", "1" = "#383838")) +
  theme(
    panel.background =  element_rect(color = NA, fill = "#dedede"),
    panel.border = element_blank(),
    panel.spacing = unit(1, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "none",
    strip.background = element_blank()
  )
dev.off()
rm(d, d1, d2)


# Figure 2 _____________________________________________________________________
# Initial opinion distributions.
#
d1 <- data.frame(
  o = 0:1000 / 1000,
  y = dbeta(0:1000/1000, shape1 = 1, shape2 = 1),
  panel = 1,
  group = 1
)
d2 <- data.frame(
  o = 0:1000 / 1000,
  y = dbeta(0:1000/1000, shape1 = 3, shape2 = 3),
  panel = 2,
  group = 1
)
d3 <- data.frame(
  o = 0:1000 / 1000,
  y = dbeta(0:1000/1000, shape1 = 3, shape2 = 3.5),
  panel = 3,
  group = 2
)
d4 <- data.frame(
  o = 0:1000 / 1000,
  y = dbeta(0:1000/1000, shape1 = 3.5, shape2 = 3),
  panel = 3,
  group = 3
)
d <- rbind(d1, d2, d3, d4)
d$group <- as.factor(d$group)

panelLabels <- c(
  "uniform, no group bias",
  "bell-shaped, no group bias",
  "bell-shaped, group bias"
)
d$panel <- factor(
  d$panel,
  levels = 1:3,
  labels = panelLabels
)
d$o <- d$o * 2 - 1

labs <- data.frame(
  label = c(
    "\u03B1=\u03B2=1",
    "\u03B1=\u03B2=3",
    "\u03B1=3, \u03B2=3.5      \u03B1=3.5, \u03B2=3"),
  panel = factor(
    1:3,
    levels = 1:3,
    labels = panelLabels
  ),
  x = 0,
  y = c(1.2,2.1,2.1)
)

png(
  filename = "./outputGraphics/figure 2 - initial opinion distributions.png",
  width = 1800, height = 700, res = 300, units = "px"
)
ggplot(data = d, aes(x = o, y = y)) +
  geom_line(aes(linetype = group), size = 1) + #, linetype = group
  xlab("attitude") + #ylab ("probability") +
  scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1, 1), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,2.25)) +
  geom_text(data = labs, aes(x = x, y = y, label = label)) +
  facet_wrap(~panel) +
  theme(
    legend.position = "none",
    panel.background =  element_blank(),
    panel.border = element_blank(),
    panel.spacing = unit(1.2, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "#808080"),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank() 
  )
dev.off()
rm(d, d1, d2, d3, d4, panelLabels, labs)



# Figure 3 _____________________________________________________________________
# Plotting distance decay functions.
#
ddf <- function(distances = 1:50000/10, s) {return(exp( - distances / s))}

d <- data.frame(
  dist = 1:50000/10,
  steep = ddf(s = 10),
  medium = ddf(s = 100),
  mild = ddf(s = 1000)
)
d <- reshape2::melt(d, id.vars = c("dist"))

labs <- data.frame(
  label = c("steep (s=10)", "medium (s=100)", "mild (s=1000)"),
  x = c(10, 100, 1000),
  y = c(0.5, 0.5, 0.5)
)

png(
  filename = "./outputGraphics/figure 3 - distance decay.png",
  width = 1000, height = 850, res = 300, units = "px"
)
ggplot(data = d) +
  geom_line(aes(x = dist, y = value, linetype = variable), size = 1) +
  xlab("distance (meters)") + ylab ("probability of interaction") +
  scale_x_continuous(
    trans='log10', breaks = c(1, 10, 100, 1000, 5000),
    limits = c(1, 5000), expand = c(0,0)
  ) +
  #geom_label(data = labs, aes(x = x, y = y, label = label), fill = NA, label.size	= 0) +
  geom_text(data = labs, aes(x = x, y = y, label = label), angle = 293) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,1), expand = c(0,0)) +
  theme(
    legend.position = "none",
    panel.background =  element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "#808080"),
    axis.text.x=element_text(hjust=0.8)
  )
dev.off()


#_______________________________________________________________________________
################################################################################
#################### Plotting simulation results ###############################
################################################################################

# Defining district labels that will be used throughout:
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

# A function to plot violin plots:
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

# Loading simulation results file:
#load("./cityData/geodata_Rotterdam.RData")
load("./simOutput/completeDataset.RDATA") ##
#load("./simOutput/iniAli.RDATA")
# This has loaded two dataframes into memory:
#   >"r": each row is a simulation run; the columns are parameters and 
#    district-level outcome metrics;
#   >"ri": each row is an agent; the columns are parameters and agent-level
#    metrics.
# The key of these two dataset is the variable "seed", the unique ID of each
# simulation run.


# recoding variables.
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





# Figure 4 _____________________________________________________________________
# Expectation 1a)
# Agents who are more exposed to outgroup agents develop extreme attitudes after
# fewer interactions.
#
# We filter out agents who never developed extreme attitudes:
rri2 <- rri[!is.na(rri$timeFirstExtr),]
png(
  filename = "./outputGraphics/figure 4 - expectation_1a.png",
  width = 1800, height = 1400, res = 300, units = "px"
)
ggplot(
  rri2,
  aes(
    cut(expOutgr2, breaks=(c(0:5)/5)),
    log10(nIntFirstExtr)
  )) +
  ylab("number of interactions to\nfirst extremization (log_10)") +
  xlab("outgroup exposure (s=100)") +
  geom_violin(
    fill = "gray",
    draw_quantiles = c(0.5)#,
    #bw = "bcv"
  ) +
  facet_wrap(rri2$wijk, labeller = as_labeller(districtLabels)) +
  theme(
    plot.margin=unit(c(0,0,0,0),"pt"),
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )
dev.off()
rm(rri2)

# Figure 5 _____________________________________________________________________
# Expectation 1b)
# Districts with higher levels of mean outgroup exposure exhibit a higher degree
# of polarization at any given point of time (i.e. measured as the average
# number of interaction events per agent).
#
#png(
#  filename = "./outputGraphics/figure 5 - expectation_1b.png",
#  width = 1300, height = 1200, res = 300, units = "px"
#)
#print(districtViolins(
#  depVar="SDopinions",
#  depVar2 = "iniSDopinions",
#  indepVar="expOutgr2",
#  depVarLabel="attitude polarization",
#  indepVarLabel="average outgroup exposure (s=100)"))
#dev.off()

labs_prep <- citySummary[order(citySummary$expOutgr2),]$district
labs <- c()
for (i in 1:length(labs_prep)){
  labs[i] <- paste0(
    labs_prep[i],
    " (",
    round(
      subset(citySummary, citySummary$district==labs_prep[i])$expOutgr2,
      digits = 3
    ),
    ")"
  )
}

png(
  filename = "./outputGraphics/figure 5 - expectation_1b.png",
  width = 1300, height = 1200, res = 300, units = "px"
)

ggplot(rr, aes(factor(expOutgr2), SDopinions)) +
  geom_col( # Max polarization
    data = citySummary,
    #aes(x = as.factor(expOutgr2), y = maxAlignment2),
    aes(x = as.factor(expOutgr2), y = maxSDopinions),
    fill = "#ededed"
  ) +
  geom_violin( # Initial polarization
    aes(y = iniSDopinions),
    fill = "white", color = "#ababab",
    scale = "width", width = 0.7,
    draw_quantiles = 0.5
  ) +
  geom_violin( # Final polarization
    fill="gray",
    scale = "width", width = 0.7,
    draw_quantiles = 0.5
  ) + 
  scale_x_discrete(labels=labs) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
  ggtitle("districts ordered by\naverage outgroup exposure (s=100)") +
  ylab("attitude polarization") +
  theme(
    plot.margin=unit(c(0,0,0,40),"pt"),
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )

dev.off()
  



# Figure # _____________________________________________________________________
# Showing non-relevant runs from the baseline:
#png(
#  filename = "./outputGraphics/figure X - polarized runs.png",
#  width = 1200, height = 1000, res = 300, units = "px"
#)
#hist(
#  rr$SDopinions, breaks=10,
#  main="baseline parameter configuration",
#  xlab="opinion polarization", yaxt="n", ylab="", yaxt='n',
#  col="black", border="white"
#)
#abline(v=0.3, col= "red")
#title(ylab="relative frequency", line=0)
#dev.off()
#
# And here are all runs (not only the baseline)
#png(
#  filename = "./outputGraphics/filtering polarized runs _ all runs.png",
#  width = 1200, height = 1000, res = 300, units = "px"
#)
#hist(
#  r$polarizationIndex, breaks=10,
#  main="all parameter configurations",
#  xlab="polarization index", yaxt="n", ylab="", yaxt='n',
#  col="black", border="white"
#)
#abline(v=0.3, col= "red")
#title(ylab="relative frequency", line=0)
#dev.off()
#
#
#
# Removing non-relevant (i.e. non-polarized) runs:
#rr <- subset(rr, rr$polarizationIndex >= 0.3)
#rri <- subset(rri, rri$seed %in% rr$seed)




# Figure 6 _____________________________________________________________________
# Expectation 2a)
# Agents who are more exposed to their outgroup exhibit higher scores of local
# alignment.
#
# First we only select few representative districts to plot in Figure 6:
rri2 <- rri[rri$wijk %in% c(1, 4, 7),]

# If we want to only plot one run per district:
#runs <- c()
#for (w in 1:12) {runs[w] <- unique(rri[rri$wijk == w,]$seed)[3]}
#rri2 <- rri2[rri2$seed %in% runs,]

# Picking an alternative initial configuration:
#rri2 <- subset(
#  ri,
#  ri$initialOpinionDistribution == "beta" &
#    ri$H == 0.6 &ri$distanceDecay == 2 & ri$wijk %in% c(1, 4, 7)
#)

png(
  filename = "./outputGraphics/figure 6 - expectation_2a.png",
  width = 1100, height = 1200, res = 300, units = "px"
)
ggplot(
  rri2,
  aes(
    cut(expOutgr2, breaks=(0:5 / 5)),
    abs(opAlignment2)
  )) +
  #ggtitle("3rd run per district") +
  #ggtitle("no group bias") +
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
dev.off()
rm(rri2)

# Then we save a plot for the appendix (Figure A13) that contains all districts:
png(
  filename = "./outputGraphics/figure A1 - expectation_2a.png",
  width = 1100, height = 3000, res = 300, units = "px"
)
ggplot(
  rri,
  aes(
    cut(expOutgr2, breaks=(0:5 / 5)),
    abs(opAlignment2)#, fill=as.factor(rri$group)
  )) +
  ylab("local alignment (s=100)") +
  xlab("outgroup exposure (s=100)") +
  geom_violin( # Alignment at t=0 (white)
    aes(y = abs(rri$iniOpAlignment2)),
    fill = "white", color = "#ababab", scale = "width", draw_quantiles = 0.5) +
  geom_violin( # Alignment at t = 200 (gray)
    fill = "gray",
    scale = "width",
    draw_quantiles = 0.5
  ) +
  facet_grid(
    rri$wijk ~ (rri$group - 1),
    labeller = as_labeller(
      c(districtLabels, c("0"="non-western", "-2" = "western")))) +
  theme(
    plot.margin=unit(c(0,0,0,00),"pt"),
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(color=NA,fill="gray97"),
    axis.line = element_line(colour = "black"),
    strip.text.y = element_text(angle = 0, hjust=0),
    strip.background.y = element_rect(colour=NA, fill=NA)
  )
dev.off()



# Figure 7 _____________________________________________________________________
# Expectation 2b) 
# Districts with higher levels of mean outgroup exposure exhibit higher average
# local alignment.

#png(
#  filename = "./outputGraphics/figure 7 - expectation_2b.png",
#  width = 1300, height = 1200, res = 300, units = "px"
#)
#print(districtViolins(
#  data = rr,
#  depVar = "opAlignment2",
#  indepVar = "expOutgr2",
#  depVar2 = "iniAli2",
#  depVarLabel = "average local alignment (s=100)",
#  indepVarLabel = "average outgroup exposure (s=100)"
#))
#dev.off()

labs_prep <- citySummary[order(citySummary$expOutgr2),]$district
labs <- c()
for (i in 1:length(labs_prep)){
  labs[i] <- paste0(
    labs_prep[i],
    " (",
    round(
      subset(citySummary, citySummary$district==labs_prep[i])$expOutgr2,
      digits = 3
    ),
    ")"
  )
}

png(
  filename = "./outputGraphics/figure 7 - expectation_2b.png",
  width = 1300, height = 1200, res = 300, units = "px"
)

ggplot(rr, aes(factor(expOutgr2), opAlignment2)) +
  geom_col( # Max alignment
    data = citySummary,
    aes(x = as.factor(expOutgr2), y = maxAlignment2),
    fill = "#ededed"#, width = 1.2
  ) +
  geom_violin( # Initial alignment
    aes(y = iniAli2),
    fill = "white", color = "#ababab",
    scale = "width", width = 0.7,
    draw_quantiles = 0.5
  ) +
  geom_violin( # Final alignment
    fill="gray",
    scale = "width", width = 0.7,
    draw_quantiles = 0.5
  ) +
  ylab("average local alignment (s=100)") +
  ggtitle("districts ordered by\naverage outgroup exposure (s=100)") +
  scale_x_discrete(labels=labs) +
  scale_y_continuous(expand = c(0,0)) +
  theme(
    plot.margin=unit(c(0,0,0,40),"pt"),
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )

dev.off()



# Figure 8 _____________________________________________________________________
# Expectation 2c) 
# There is an inverted U-shaped effect of average outgroup exposure
# in the district and the difference between the average attitude of
# the two groups (global level alignment).
#png(
#  filename = "./outputGraphics/figure 8 - expectation_2c.png",
#  width = 1300, height = 1200, res = 300, units = "px"
#)
#districtViolins(
#  depVar = "intuitiveAlignment",
#  depVar2 = "iniIntuitiveAlignment",
#  indepVar = "expOutgr2",
#  depVarLabel = "attitude difference between groups",
#  indepVarLabel = "average outgroup exposure (s=100)")
#dev.off()

labs_prep <- citySummary[order(citySummary[,"expOutgr2"]),]$district
labs <- c()
for (i in 1:length(labs_prep)){
  labs[i] <- paste0(
    labs_prep[i],
    " (",
    round(
      subset(citySummary, citySummary$district==labs_prep[i])[,"expOutgr2"],
      digits = 3
    ),
    ")"
  )
}

# Points will be concentrated between very low values of alignment (t=0) and
# very high (at t=200), and nothing in-between. To better plot these data, we
# can "squash" the empty portion on the y axis.
# adapted from https://stackoverflow.com/a/35514861
scaleSquash <- function(from, to, factor) {
  tr <- function(x) {
    if (any(is.na(x))) return(x)
    isq <- x > from & x < to
    ito <- x >= to
    x[isq] <- from + (x[isq] - from) / factor
    x[ito] <- from + (to - from) / factor + (x[ito] - to)
    return(x)
  }
  
  inv <- function(x) {
    if (any(is.na(x))) return(x)
    isq <- x > from & x < from + (to - from) / factor
    ito <- x >= from + (to - from) / factor
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from) / factor))
    return(x)
  }
  
  return(trans_new("scaleSquash", tr, inv))
}

cut_from = 0.26
cut_to = 1.66
breaks = 0:20/10

png(
  filename = "./outputGraphics/figure 8 - expectation_2c.png",
  width = 1300, height = 1200, res = 300, units = "px"
)
ggplot(rr, aes(factor(expOutgr2), abs(intuitiveAlignment))) +
  geom_violin( # alignment at t=0
  aes(y = iniIntuitiveAlignment),
  fill = "white", color = "#ababab",
  scale = "width",
  draw_quantiles = 0.5
  ) +
  geom_violin( # alignment at t=200
    fill="gray",
    scale = "width",
    draw_quantiles = 0.5
  ) +
  ggtitle("districts ordered by\naverage outgroup exposure (s=100)") +
  ylab("between-groups attitude difference") +
  scale_x_discrete(labels = labs) +
  scale_y_continuous(
    trans = scaleSquash(cut_from, cut_to, 50), ### squashing the y axis
    limits = c(0,2.1),
    expand = c(0,0),
    breaks = breaks[breaks <= cut_from | breaks >= cut_to]
  ) +
  geom_rect( # Highlighting the squashed portion of the chart
    aes(xmin = 0, xmax = 13, ymin = cut_from, ymax = cut_to),
    fill = "#fafafa"
  ) +
  geom_segment( # highlighting the squashed portion on the Y axis
    aes(x = 0, xend = 0, y = 0, yend = cut_from), size = 0.5, color = "black") +
  geom_segment( # highlighting the squashed portion on the Y axis
    aes(x = 0, xend = 0, y = cut_to, yend = 2), size = 0.5, color = "black") +
  #geom_segment( # highlighting the squashed portion on the Y axis
  #  aes(x = 0, xend = 0, y = cut_from, yend = cut_to),
  #  size = 1, color = "#bdbdbd") +
  #geom_text(aes(x = 0, y = cut_from, label = "test", angle = 45)) +
  #geom_label(aes(x = 0, y = cut_from, label = "test")) +
  annotate("text", x = 0, y = cut_from, label = "\\", angle = "300", size = 4) +
  annotate("text", x = 0, y = cut_to, label = "\\", angle = "300", size = 4) +
  coord_cartesian(clip = "off") +
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
dev.off()


# Robustness ===================================================================
#
#

# Robustness to H
rr <- subset(
  r,
  r$initialOpinionDistribution == "groupBias" & 
    r$distanceDecay == 2
)
rr$opAlignment1 <- abs(rr$opAlignment1)
rr$opAlignment2 <- abs(rr$opAlignment2)
rr$opAlignment3 <- abs(rr$opAlignment3)
rr$iniAli1 <- abs(rr$iniAli1)
rr$iniAli2 <- abs(rr$iniAli2)
rr$iniAli3 <- abs(rr$iniAli3)
rri <- subset(
  ri,
  ri$initialOpinionDistribution == "groupBias" & 
    ri$distanceDecay == 2
)
# 1a)
png(
  filename = "./outputGraphics/robustness_H_1a.png",
  width = 2000, height = 1000, res = 300, units = "px"
)
violinPlotOld(
  depVar = log10(rri$nIntFirstExtr),
  indepVar = rri$expOutgr2,
  depVarLabel = "time of first extremization (log_10)",
  indepVarLabel = "outgroup exposure (s=100)"
) + facet_grid(
  cols=vars(rri$H), labeller=as_labeller(c("0.6"="H=0.6","0.9"="H=0.9")))
dev.off()

# 1b)
png(
  filename = "./outputGraphics/robustness_H_1b.png",
  width = 2000, height = 1000, res = 300, units = "px"
)
print(districtViolins(
  depVar="polarizationIndex",
  indepVar="expOutgr2",
  depVarLabel="polarization index",
  indepVarLabel="average outgroup exposure (s=100)")) + 
  facet_grid(
    cols=vars(rr$H), labeller=as_labeller(c("0.6"="H=0.6","0.9"="H=0.9")))
dev.off()

# 2a)
#rr <- subset(rr, rr$polarizationIndex >= 0.3)
#rri <- subset(rri, rri$seed %in% rr$seed)
#png(
#  filename = "./outputGraphics/robustness_H_2a.png",
#  width = 2000, height = 2000, res = 300, units = "px"
#)
#violinPlotOld(
#  depVar = abs(rri$opAlignment2),
#  indepVar = rri$expOutgr2,
#  depVarLabel = "local alignment (s=100)",
#  indepVarLabel = "outgroup exposure (s=100)"
#) + facet_grid(
#  rri$H ~ rri$group,
#  labeller = as_labeller(
#    c("0.6"="H=0.6","0.9"="H=0.9","-1"="western","1"="non-western")))
#dev.off()
#
# 2b)
#print(districtViolins(
#  depVar="opNetAli2",
#  indepVar="expOutgr2",
#  depVarLabel="average local alignment (s=100)",
#  indepVarLabel="average outgroup exposure (s=100)")) + 
#  facet_grid(rr$H, labeller=as_labeller(c("0.6"="H=0.6","0.9"="H=0.9")))
#
# 2c)
#districtViolins(
#  depVar="intuitiveAlignment",
#  indepVar="expOutgr2",
#  depVarLabel="global alignment (abs. difference between
#  the average attitude of the two groups)",
#  indepVarLabel="average outgroup exposure (s=100)") + 
#  facet_grid(rr$H, labeller=as_labeller(c("0.6"="H=0.6","0.9"="H=0.9")))



# Robustness to initialOpinionDistribution
rr <- subset(
  r,
  r$H == "0.6" & 
    r$distanceDecay == 2
)
rr$opAlignment1 <- abs(rr$opAlignment1)
rr$opAlignment2 <- abs(rr$opAlignment2)
rr$opAlignment3 <- abs(rr$opAlignment3)
rr$iniAli1 <- abs(rr$iniAli1)
rr$iniAli2 <- abs(rr$iniAli2)
rr$iniAli3 <- abs(rr$iniAli3)
rri <- subset(
  ri,
  ri$H == 0.6 & 
    ri$distanceDecay == 2
)

# 1a)
png(
  filename = "./outputGraphics/robustness_opDistr_1a.png",
  width = 2000, height = 1000, res = 300, units = "px"
)
violinPlotOld(
  depVar = log10(rri$nIntFirstExtr),
  indepVar = rri$expOutgr2,
  depVarLabel = "time of first extremization (log_10)",
  indepVarLabel = "outgroup exposure (s=100)"
) + facet_grid(
  cols=vars(rri$initialOpinionDistribution),
  labeller=as_labeller(c(
    "uniform"="uniform",
    "groupBias"="beta\nwith group bias",
    "beta"="beta\nno group bias")))
dev.off()


# 1b)
png(
  filename = "./outputGraphics/robustness_opDistr_1b.png",
  width = 2200, height = 1100, res = 300, units = "px"
)
print(districtViolins(
  depVar="polarizationIndex",
  indepVar="expOutgr2",
  depVarLabel="polarization index",
  indepVarLabel="average outgroup exposure (s=100)")) + 
  facet_grid(
    cols = vars(rr$initialOpinionDistribution), 
    labeller=as_labeller(c(
      "uniform"="uniform",
      "groupBias"="beta\nwith group bias",
      "beta"="beta\nno group bias")))
dev.off()

# 2a)
rr <- subset(rr, rr$polarizationIndex >= 0.3)
rri <- subset(rri, rri$seed %in% rr$seed)
png(
  filename = "./outputGraphics/robustness_opDistr_2a.png",
  width = 1800, height = 1400, res = 300, units = "px"
)
violinPlotOld(
  depVar = abs(rri$opAlignment2),
  indepVar = rri$expOutgr2,
  depVarLabel = "local alignment (s=100)",
  indepVarLabel = "outgroup exposure (s=100)"
) + facet_grid(
  rri$initialOpinionDistribution ~rri$group,
  labeller=as_labeller(c(
    "-1"="western","1"="non-western",
    "uniform"="uniform",
    "groupBias"="beta\nwith group bias",
    "beta"="beta\nno group bias")))
dev.off()

# 2b)
png(
  filename = "./outputGraphics/robustness_opDistr_2b.png",
  width = 2200, height = 1100, res = 300, units = "px"
)
print(districtViolins(
  depVar="opAlignment2",
  indepVar="expOutgr2",
  depVar2="iniAli2",
  depVarLabel="global alignment (s=100)",
  indepVarLabel="average outgroup exposure (s=100)")) + 
  facet_grid(
    cols = vars(rr$initialOpinionDistribution), 
    labeller=as_labeller(c(
      "uniform"="uniform",
      "groupBias"="beta\nwith group bias",
      "beta"="beta\nno group bias")))
dev.off()

# 2c)
png(
  filename = "./outputGraphics/robustness_opDistr_2c.png",
  width = 2200, height = 1100, res = 300, units = "px"
)
districtViolins(
  depVar="intuitiveAlignment",
  indepVar="expOutgr2",
  depVarLabel="attitude difference between groups",
  indepVarLabel="average outgroup exposure (s=100)") + 
  facet_grid(
    cols = vars(rr$initialOpinionDistribution), 
    labeller=as_labeller(c(
      "uniform"="uniform",
      "groupBias"="beta\nwith group bias",
      "beta"="beta\nno group bias")))
dev.off()





# Robustness to distance decay function (in the interaction)
rr <- subset(
  r,
  r$initialOpinionDistribution == "groupBias" & 
    r$H == 0.6
)
rr$opAlignment1 <- abs(rr$opAlignment1)
rr$opAlignment2 <- abs(rr$opAlignment2)
rr$opAlignment3 <- abs(rr$opAlignment3)
rr$iniAli1 <- abs(rr$iniAli1)
rr$iniAli2 <- abs(rr$iniAli2)
rr$iniAli3 <- abs(rr$iniAli3)
rri <- subset(
  ri,
  ri$initialOpinionDistribution == "groupBias" & 
    ri$H == 0.6
)

# 1a)
png(
  filename = "./outputGraphics/robustness_DDFint_1a.png",
  width = 2000, height = 1000, res = 300, units = "px"
)
violinPlotOld(
  depVar = log10(rri$nIntFirstExtr),
  indepVar = rri$expOutgr2,
  depVarLabel = "time of first extremization (log_10)",
  indepVarLabel = "outgroup exposure (s=100)"
) + facet_grid(
  cols=vars(rri$distanceDecay),
  labeller=as_labeller(c(
    "1"="s=10",
    "2"="s=100",
    "3"="s=1000")))
dev.off()

# 1b)
png(
  filename = "./outputGraphics/robustness_DDFint_1b.png",
  width = 2200, height = 1100, res = 300, units = "px"
)
print(districtViolins(
  depVar="polarizationIndex",
  indepVar="expOutgr2",
  depVarLabel="polarization index",
  indepVarLabel="average outgroup exposure (s=100)")) + 
  facet_grid(
    cols=vars(rr$distanceDecay), 
    labeller=as_labeller(c(
      "1"="s=10",
      "2"="s=100",
      "3"="s=1000")))
dev.off()

# 2a)
rr <- subset(rr, rr$polarizationIndex >= 0.3)
rri <- subset(rri, rri$seed %in% rr$seed)
png(
  filename = "./outputGraphics/robustness_DDFint_2a.png",
  width = 1800, height = 1400, res = 300, units = "px"
)
violinPlotOld(
  depVar = abs(rri$opAlignment2),
  indepVar = rri$expOutgr2,
  depVarLabel = "local alignment (s=100)",
  indepVarLabel = "outgroup exposure (s=100)"
) + facet_grid(
  rri$distanceDecay ~(rri$group - 1),
  labeller=as_labeller(c(
    "-2"="western","0"="non-western",
    "1"="s=10",
    "2"="s=100",
    "3"="s=1000")))
dev.off()

# 2b)
png(
  filename = "./outputGraphics/robustness_DDFint_2b.png",
  width = 2200, height = 1100, res = 300, units = "px"
)
print(districtViolins(
  depVar="opAlignment2",
  indepVar="expOutgr2",
  depVar2="iniAli2",
  depVarLabel="global alignment (s=100)",
  indepVarLabel="average outgroup exposure (s=100)"))+ 
  facet_grid(
    cols=vars(rr$distanceDecay), 
    labeller=as_labeller(c(
      "1"="s=10",
      "2"="s=100",
      "3"="s=1000")))
dev.off()

# 2c)
png(
  filename = "./outputGraphics/robustness_DDFint_2c.png",
  width = 2200, height = 1100, res = 300, units = "px"
)
districtViolins(
  depVar="intuitiveAlignment",
  indepVar="expOutgr2",
  depVarLabel="attitude difference between groups",
  indepVarLabel="average outgroup exposure (s=100)") + 
  facet_grid(
    cols=vars(rr$distanceDecay), 
    labeller=as_labeller(c(
      "1"="s=10",
      "2"="s=100",
      "3"="s=1000")))
dev.off()






# Robustness to distance decay function (in the outcome measures)
rr <- subset(
  r,
  r$initialOpinionDistribution == "groupBias" & 
    r$H == 0.6 &
    r$distanceDecay == 2
)
rr$opAlignment1 <- abs(rr$opAlignment1)
rr$opAlignment2 <- abs(rr$opAlignment2)
rr$opAlignment3 <- abs(rr$opAlignment3)
rr$iniAli1 <- abs(rr$iniAli1)
rr$iniAli2 <- abs(rr$iniAli2)
rr$iniAli3 <- abs(rr$iniAli3)
rri <- subset(
  ri,
  ri$initialOpinionDistribution == "groupBias" & 
    ri$H == 0.6 &
    ri$distanceDecay == 2
)

# 1a
violinPlotOld(
  depVar = log10(rri$nIntFirstExtr),
  indepVar = rri$expOutgr1,
  depVarLabel = "time of first extremization (log_10)",
  indepVarLabel = "outgroup exposure (s=10)"
)
violinPlotOld(
  depVar = log10(rri$nIntFirstExtr),
  indepVar = rri$expOutgr3,
  depVarLabel = "time of first extremization (log_10)",
  indepVarLabel = "outgroup exposure (s=1000)"
)

#1b)
print(districtViolins(
  depVar="polarizationIndex",
  indepVar="expOutgr1",
  depVarLabel="polarization index",
  indepVarLabel="average outgroup exposure (s=10)"))
print(districtViolins(
  depVar="polarizationIndex",
  indepVar="expOutgr3",
  depVarLabel="polarization index",
  indepVarLabel="average outgroup exposure (s=1000)"))

#2a)
rr <- subset(rr, rr$polarizationIndex >= 0.3)
rri <- subset(rri, rri$seed %in% rr$seed)

violinPlotOld(
  depVar = abs(rri$opAlignment3),
  indepVar = rri$expOutgr1,
  depVarLabel = "local alignment (s=10)",
  indepVarLabel = "outgroup exposure (s=10)"
)
violinPlotOld(
  depVar = abs(rri$opAlignment3),
  indepVar = rri$expOutgr3,
  depVarLabel = "local alignment (s=1000)",
  indepVarLabel = "outgroup exposure (s=1000)"
)

#2b)
print(districtViolins(
  depVar="opAlignment1",
  indepVar="expOutgr1",
  depVar2="iniAli2",
  depVarLabel="global alignment (s=10)",
  indepVarLabel="average outgroup exposure (s=10)"))
print(districtViolins(
  depVar="opAlignment3",
  indepVar="expOutgr3",
  depVar2="iniAli2",
  depVarLabel="global alignment (s=1000)",
  indepVarLabel="average outgroup exposure (s=1000)"))



# 2c)
districtViolins(
  depVar="intuitiveAlignment",
  indepVar="expOutgr1",
  depVarLabel="attitude difference between groups",
  indepVarLabel="average outgroup exposure (s=10)")
districtViolins(
  depVar="intuitiveAlignment",
  indepVar="expOutgr3",
  depVarLabel="attitude difference between groups",
  indepVarLabel="average outgroup exposure (s=1000)")

















