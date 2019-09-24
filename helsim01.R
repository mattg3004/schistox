#helsim01 - R version of the IBM in Medley 1989 thesis, and Anderson&Medley 1985
#Graham Medley, July 2013
#
#Aspects missing:
# host population dynamics (aging, birth, death)
# age-related exposure function
# schistsome free-living populations (currently only two free-living: infective and uninfective)
# multiple immature stages in host (currently only zero or one)
# more esoteric chemotherapy (currently only mass random)

#clear the decks
rm(list=ls(all.names=TRUE))
closeAllConnections()

#define the functions
source ( "helsimFNS01.R" )

#read parameters
params <- readParams ( "params.txt")

#set some defaults
maxStep = 1/52 #time steps never exceed this for deterministic update of freeliving popns

#set up the output and recording
outFile <- file ( paste(params$runName,"Out.txt",sep="" ), "w" )
outInt <- 1 / params$outFreq
outTimes <- seq ( 0, params$maxTime, outInt )
resWorms <- array(0,dim=c(params$numReps,length(outTimes),params$N))
resFL <- array(0,dim=c(params$numReps,length(outTimes),2))

for ( repNum in 1:params$numReps ) {
  cat ( "repetition",repNum )
  #setup the simulation
  simData <- setupSD ( params )
  time <- 0
  FLlast <- 0

  #record first time point
  resWorms[repNum,1,] <- simData$worms[,1]
  resFL[repNum,1,] <- simData$freeLiving
  outRes ( outFile, repNum, 0, simData )
  outCount <- 2
  nextOutTime <- outInt
  chemoTimes <- params$chemos[,2]
  nextChemoTime <- min(chemoTimes)
  nextChemo <- which.min(chemoTimes)
  nextStep <- min ( nextOutTime, time+maxStep, nextChemoTime )

  while ( time<params$maxTime ) {
    rates <- calcRates ( params, simData )
    tstep <- rexp(1,sum(rates))
    if ((time+tstep)<nextStep) {
      time <- time+tstep
      simData <- doEvent ( rates, params, simData )
    } else {
      simData <- doFreeLive ( params, simData, nextStep-FLlast )
      FLlast <- nextStep
      if (time+tstep>nextChemoTime) { #chemotherapy
        cat("c")
        simData <- doChemo ( nextChemo, params, simData )
        chemoTimes[nextChemo] <- params$maxTime+1000
        nextChemoTime <- min(chemoTimes)
        nextChemo <- which.min(chemoTimes)
      }
      if ((time+tstep)>nextOutTime){ #record
        resWorms[repNum,outCount,] <- simData$worms[,1]
        resFL[repNum,outCount,] <- simData$freeLiving
        outRes ( outFile, repNum, nextOutTime, simData )
        outCount <- outCount+1
        nextOutTime <- nextOutTime+outInt
        cat ( ".")
      }
      time <- nextStep
      nextStep <- min ( nextOutTime, time+maxStep, nextChemoTime )

    }
  } #while time
  cat ( "\n" )
} #for repetitions

#plot means
means = apply ( resWorms, c(1,2), mean )
plot(c(0,tail(outTimes,1)),c(0,max(means)),"n",main="Mean Adult Burdens",xlab="Time")
for ( ii in 1:params$numReps ) lines ( outTimes, means[ii,] )

#plot variances
vars = apply ( resWorms, c(1,2), var )
plot ( c(0,tail(outTimes,1)),c(0,max(vars)),"n",main="Variance Adult Burdens",xlab="Time")
for ( ii in 1:params$numReps ) lines ( outTimes, vars[ii,])

#plot individual burdens
numToPlot = min ( 10, params$N ) #number of individuals to plot
plot ( c(0,tail(outTimes,1)), c(0,max(resWorms[params$numReps,,1:numToPlot])),"n",
       main="Indiv Adult Burdens", xlab="Time", ylab="Worms", sub="last rep")
for ( ii in 1:numToPlot ) lines ( outTimes, resWorms[params$numReps,,ii] )

#plot distribution of worms
hist(resWorms[params$numReps,length(outTimes),],main="Distribution",xlab="# adult worms",
     sub="last rep, last time")

#plot correlation between susceptibility and burden if relevant
if ( params$hostHet[1]>0 ) plot (simData$si,resWorms[params$numReps,length(outTimes),],
                              main="susceptibility",xlab="s_i",ylab="#worms",sub="last rep, last time")

#plot free living stages
plot(c(0,tail(outTimes,1)),c(0,max(resFL[params$numReps,,])),"n",ylog=TRUE,
     main="Freeliving populations",xlab="Time")
lines(outTimes,resFL[params$numReps,,1])
lines(outTimes,resFL[params$numReps,,2],col="red")

close ( outFile )
