########################################################################################################
#
# This is a very quick and (very) dirty piece of code to compute probabilities of alignments with stars 
# in a tumulus. Inspired by: https://hal.science/hal-03223357v3
#
# The general strategy is to load in memory some information about a tumulus (namely, the azimuths of 
# stones as seen from the center of the tumulus) and randomly shuffle these positions to check if the 
# alignments observed could happen by chance.
#
# The two things that the program needs are: 
# 
#   (1) A file describing the position of stones in the tumulus (tumulus_azimuths.csv)
#   (2) A file with the list of azimuths to consider (typically, azimuth for the rise/set of notable stars): objects.csv
#
#

# Change the next command (setwd) to point to the folder on your computer where you have downloaded the two csv files 
# that come with this code.
setwd("C:/Users/jf/Documents/perso/2024/Tumulus")

RESOLUTION <- 0.1 # Steps in resolution for angles calculation.

# Let's start by defining two functions that will be needed for these simulations.

################################################################################
# randomizePositions
#
# This function randomizes the position of stones on the circle.
# It takes as parameters: tz and minSpace.
# 
#   tz: a data frame containing the list of stones in the tumulus, with at least the following columns:
#       - az1 & az2: azimuths of the two ends of the stone, as seen from the center
#       - angSize: the angular size of the stone, as seen from the center (can be obtained easily by doing the difference of az2 and az1)
#
#   minSpace: the minimum angular separation between two stones. Default value is 3.0 degrees.
randomizePositions <- function(tz, minSpace = 3.0){
  
  # I start by making a copy of the tumulus data.frame
  tzr <- tz
  
  # I'm going to generate a list of all positions available in increment of 0.1 degree.
  # After positioning each stone, I will remove the positions occupied by this stone to the list of possible positions for the next random draw (the next stone)
  positionsAvailable <- seq(from=0, to=360, by = RESOLUTION)
  
  # For each stone, I will randomly draw a new position and remove this position from the list 
  # of available positions (we don't stack stones on top of each other... In fact, there is a parameter 
  # 'minSpace' to require a minimum distance (in azimuth) between two stones.)
  for(i in (1:nrow(tzr))){
    initialAz1 <- tz$az1[i]
    angZise <- tz$angSize[i]
    newAz1 <- sample(positionsAvailable, 1)
    newAz2 <- (newAz1 + angZise)%%360
    tzr$az1[i] <- newAz1
    tzr$az2[i] <- newAz2
    
    startRemovePosition <- (newAz1 - minSpace) %% 360
    endRemovePosition <- (newAz2 + minSpace) %% 360
    widthToRemove <- (endRemovePosition - startRemovePosition) %% 360
    
    vr <- c() # vr will contain the list of positions occupied by this stone (and therefore to remove for the 
    # next random draw).
    # If the positions span the absolute North (zero azimuth)
    if( endRemovePosition < startRemovePosition ){
      v1 <- seq(from = startRemovePosition, to = 360, by = RESOLUTION)
      v2 <- seq(from = 0, to = endRemovePosition, by = RESOLUTION)
      vr <- c(v1, v2)
    } else {
      vr <- seq(from = startRemovePosition, to = endRemovePosition, by = RESOLUTION)
    }
    positionsAvailable <- positionsAvailable[which( is.element(positionsAvailable, vr)==F)]
  }
  tzr
}


#####################################################################################################
# getListOfPositionsCovered
#
# This function return the list of positions that would be considered as being part of an alignment
#
# It takes two parameters: t and precision
#
#   t: a data frame with the list of stones' positions (typically, returned by randomizePositions)
#
#   precision: similar to the 'résolution' parameter in https://hal.science/hal-03223357v3
#   It consider the star to be aligned if it misses the stone by no more than 'precision'
#   Default value: 1.0 degrees.
getListOfPositionsCovered <- function(t, precision = 1.0){
  vpos <- c()
  for(i in (1:nrow(t))){
    az1 <- (t$az1[i] - precision) %% 360
    az2 <- (t$az2[i] + precision) %% 360
    if( az1 > az2 ){ # If the positions covered by this stone span the azimuth zero
      v1 <- seq(from = az1, to = 360.0, by = RESOLUTION)
      v2 <- seq(from = 0.0, to = az2, by = RESOLUTION)
      vpos <- c(vpos, v1, v2)
    } else {
      v <- seq(from = az1, to = az2, by = RESOLUTION)
      vpos <- c(vpos, v)
    }
  }
  vpos
}

# We start by loading in memory the information about the tumulus and the list of azimuths of interest.
tz <- read.csv("tumulus_azimuths.csv")
tobj <- read.csv("objects.csv")

NB_RANDS <- 10000 # This is the number of random draws to be done. Going much above 10,000 will take a while, with 
# very little benefit because we likely end up re-sampling the same sets of positions.

# Preparing a data frame (tRes) to store the resuls for each random draw
# nbRiseOrSet: number of objects for which there is an alignment for either rise or set
# nbRise: number of objects for which there is an alignment for the rise
# nbRiseAndSet: number of objects for which there is an alignment for both rise and set
# All these are initialised to zero
tRes <- data.frame( round = seq(from=1, to = NB_RANDS, by = 1) )
tRes$nbRiseOrSet <- 0
tRes$nbRise <- 0
tRes$nbRiseAndSet <- 0

# This is the main loop, which will perform a new random draw at each iteration and store the results.
# It should take about 1 minute to run for 10,000 random draws (unless you have a super fast/slow computer)
for( i in (1:NB_RANDS)){
  tzr <- randomizePositions(tz, minSpace = 3.0)
  positionsCovered <- round( getListOfPositionsCovered(t = tzr, precision = 1.0) , 1 )
  tobj$riseMatch <- is.element(tobj$rise, positionsCovered)
  tobj$setMatch <- is.element(tobj$set, positionsCovered)
  
  tRes$nbRise[i] <- length(which(tobj$riseMatch==T))
  tRes$nbRiseOrSet[i] <- length(which(tobj$riseMatch==T | tobj$setMatch==T))
  tRes$nbRiseAndSet[i] <- length(which(tobj$riseMatch==T & tobj$setMatch==T))
}

cat("Probability to observe at least one object with both rise and set in an alignment: ", 
    length(which(tRes$nbRiseAndSet>0)) / nrow(tRes), 
    "\n"
)

cat("Mean number of objects with an alignment for their rise: ", 
    mean(tRes$nbRise), " (out of: ", nrow(tobj), ")\n")


hist(tRes$nbRise)
hist(tRes$nbRiseAndSet)

