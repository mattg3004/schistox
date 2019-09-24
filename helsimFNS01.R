#reads parameters from file
readParams <- function (fileName) {
  #define the connection
  inFile <- file ( fileName, "r")
  #read the simulation size variables
  runName <- strsplit ( readLines ( inFile, n=1 ), split=" ")[[1]][1]
  numReps <- as.integer(strsplit ( readLines ( inFile, n=1), split="-" )[[1]][1])
  maxTime <- as.integer(strsplit ( readLines ( inFile, n=1), split="-" )[[1]][1])
  outFreq <- as.integer(strsplit ( readLines ( inFile, n=1), split="-" )[[1]][1])
  N <- as.integer(strsplit ( readLines ( inFile, n=1), split="-" )[[1]][1])
  #read the simulation control variables
  readList <- function () { #reads a set of numbers from one line
    ip <- strsplit ( readLines ( inFile, n=1), split="-" )[[1]][1]
    ip <- as.double ( strsplit ( ip, split=" ")[[1]])
  }
  hostHet <- readList()
  pickUp <- readList()
  devStages <- readList()
  numWormPops <- devStages[1]+1
  mu <- as.double(strsplit ( readLines ( inFile, n=1), split="-" )[[1]][1])
  freeStages <- readList()
  eggProd <- readList()
  contact <- readList()
  numChemos <- as.integer(strsplit ( readLines ( inFile, n=1), split="-" )[[1]][1])
  if (numChemos==0) {
    chemos <- NULL
  } else {
    chemos <- array(0,dim=c(numChemos,4))
    for (ii in 1:numChemos) {
      chemos[ii,] <- readList()
    }
  } 
  
  #output something for the user's entertainment
  cat ( "Running", runName, "with", numReps, "repetitions", "over", maxTime, "years\n")
  if (pickUp[1]==0){
    u <- 1
  } else {
    if (pickUp[1]==1) {
      u <- -log(pickUp[2])*pickUp[2]/(1-pickUp[2])
      cat ( "   hetergeneous pickup with mean", 1/u, "\n")
    }
  }

  
  close ( inFile )
  pars <- list ( runName=runName, numReps=numReps, maxTime=maxTime, outFreq=outFreq, N=N,
                 hostHet=hostHet, pickUp=pickUp, contact=contact, u=u,
                 numWormPops=numWormPops, devStages=devStages, freeStages=freeStages,
                 eggProd=eggProd, mu=mu,
                 numChemos=numChemos, chemos=chemos )
  return (pars)
}

#sets up the simulation to initial conditions based on analytical equilibria
setupSD <-  function ( pars ) {
  
  #host contact heterogeneity
  if(pars$hostHet[1]==0) { 
    si <- rep(1,pars$N)
  }
  if(pars$hostHet[1]==1){
    si <- rgamma(pars$N,pars$hostHet[2],rate=pars$hostHet[2])
  }
  
  #pickup function
  calcLoCDF <- function ( n, p ){ #return cdf of logistic distribution
    q = 1 - p
    dist <- rep ( 0, n )
    dist[ 1 ] <- -q / log( p )
    for ( ii in 2:(n-1) ) {
      dist[ ii ] <- dist[ii-1]*q*(ii-1)/ii
    }
    dist <- cumsum ( dist )
    dist[n]<- 1
    return ( dist )
  }
  if ( pars$pickUp[1]==0 ) {
    pickUpCDF <- NULL #dummy
  }
  if ( pars$pickUp[1]==1 ) {
    pickUpCDF <- calcLoCDF ( 500, pars$pickUp[2] ) #logarithmic distribution
  }
  
  #worm populations - [N,1] are the adults
  worms <- array(10,dim=c(pars$N,pars$numWormPops))
  
  #freeliving populations: uninfective and infective
  freeLiving <- c(1000,1000)
  
  SD <- list ( ID=1:pars$N, si=si, worms=worms, freeLiving=freeLiving,
               pickUpCDF=pickUpCDF )
  
  return(SD)
}

#output the simulation
outRes <- function ( f, r, t, SD ) {
  cat ( file=f, r, t, SD$worms[,1], "\n") #adults only
  return()
}

#calculate the event rates
calcRates <- function ( pars, SD ) {
  
  #worm death
  r1 <- pars$mu*sum(SD$worms)
  
  #worm contact
  if ( pars$contact[1]==0 ) r2 <- sum(SD$si)*pars$contact[2]*pars$u
  if ( pars$contact[1]==1 ) r2 <- 
    sum(SD$si)*pars$contact[2]*SD$freeLiving[2]*pars$u/pars$N
  
  #worm maturation - only works with one
  if ( pars$devStage[1]==0) {
    r3 <- 0
  } else {
    r3 <- pars$devStage[2]*sum(SD$worms[,2]) }
  
  rates <- c ( r1, r2, r3 )
  
  return ( rates )
} #calcRates ****

#enact an event
doEvent <- function ( rates, pars, SD ) {
  
  #determine which event
  event <- which((runif(1)*sum(rates))<cumsum(rates))[1]
  
  if (event==1) { #worm death
    indiv <- which ( runif(1)*sum(SD$worms)<cumsum(rowSums(SD$worms)))[1]
    popn <- which ( runif(1)*sum(SD$worms[indiv,])<cumsum(SD$worms[indiv,]))[1]
    SD$worms[indiv,popn] <- SD$worms[indiv,popn]-1
    return ( SD )
  }
  if (event==2) { #worm contact and pickup
    indiv <- which ( (runif(1)*sum(SD$si))<cumsum(SD$si))[1]
    if (pars$pickUp[1]==0) {
      SD$worms[indiv,pars$numWormPops] <- SD$worms[indiv,pars$numWormPops]+pars$pickUp[2]
      SD$freeLiving[2] <- SD$freeLiving[2]-pars$pickUp[2]
      return ( SD )
    }
    if (pars$pickUp[1]==1) {
      pickup <- which (runif(1)<SD$pickUpCDF)[1]
      SD$worms[indiv,pars$numWormPops] <- SD$worms[indiv,pars$numWormPops]+pickup
      SD$freeLiving[2] <- SD$freeLiving[2]-pickup
      return ( SD )
    }
  }
  if (event==3) { #worm maturation
    indiv <- which ( (runif(1)*sum(SD$worms[,2]))<cumsum(SD$worms[,2]))[1]
    SD$worms[indiv,1] <- SD$worms[indiv,1]+1
    SD$worms[indiv,2] <- SD$worms[indiv,2]-1
    return ( SD )
  }

} #doEvent *****

#update the freeliving populations deterministically
doFreeLive <- function ( pars, SD, ts ) {
  eggs <- ts*sum(pars$eggProd[1]*SD$worms[,1]*exp(-SD$worms[,1]*pars$eggProd[2]))
  leaving1 <- SD$freeLiving[1]*sum(pars$freeStages[1:2])*ts
  SD$freeLiving[1] <- SD$freeLiving[1] - leaving1 + eggs
  mature1 <- leaving1 * pars$freeStages[2]/sum(pars$freeStages[1:2])
  SD$freeLiving[2] <- SD$freeLiving[2]*(1-pars$freeStages[3]*ts) + mature1
  return ( SD )
} #freeLive ****

#chemotherapy
doChemo <- function ( cs, pars, SD ) {
  if ( pars$chemos[cs,1] == 1) { #mass random
    toTreat <- runif(pars$N) <= pars$chemos[cs,3]
    SD$worms[toTreat,1] <- rbinom(size=SD$worms[toTreat,1],n=sum(toTreat),p=1-pars$chemos[cs,4])
  }
  return ( SD )
  
} #doChemo ****
