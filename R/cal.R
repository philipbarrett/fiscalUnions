#############################################################################
# cal.R
#
# Script to calibrate the fiscal unions model
# Philip Barrett, Chicago
# Created: 04jul2016
#############################################################################
rm(list=ls())
library(readr)
library(vars)
library(mvtnorm)
library(markovchain)
library(expm)
library(reshape)
library(mFilter)

Rcpp::sourceCpp('code/2016/fiscalUnions/R/rcpp/dists.cpp')
source('code/2016/fiscalUnions/R/AR1disc.R')

## 0. Script variables
gdp.file <- '~/Dropbox/data/2016/fiscalUnions/gdpDefANA.csv'
    # OECD annual GDP data
alp.file <- '~/Dropbox/data/2016/fiscalUnions/labProd.csv'
    # Labor productivity
tax.file <- '~/Dropbox/data/2016/fiscalUnions/govtNomTaxANA.csv'
    # This data source (OECD ANAs) gives consistent measures of government 
    # revenue, expenditure, and deficits.  The only downside is that it has no 
    # German data before 1995.  For this, I need to stitch on the canned 
    # revenue/GDP from the OECD taxation database.
    # QQ: Can I use W German data prior to 1995??
cons.file <- '~/Dropbox/data/2016/fiscalUnions/govtConsExpd.csv'
    # Government consumption expenditure data
int.file <- '~/Dropbox/data/2016/fiscalUnions/interestRates.csv'
    # Interest rates
debt.file <- '~/Dropbox/data/2016/fiscalUnions/govtDebt.csv'
    # Debt/GDP ratios
pop.file <- '~/Dropbox/data/2016/fiscalUnions/population.csv'
    # Population data
cts <- c("Germany", "France")
alp.code <- c( 'T_GDPHRS', 'C')       # GDP per hour worked, current prices
expd.code <- c( 'GTE', 'GS13')        # Total expenditure, general government
nom.gdp.code <- c( 'B1_GA', 'C' )     # GDP, current prices
gdp.def.code <- c( 'B1_GA', 'DOB' )   # GDP, deflator (could use expenditure deflator here)

y.min <- 1980     # 1960
nn <- 5000        # Simulation periods
n.Z <- 15         # Number of discretized single-country points
n.Z.joint <- 30   # For the joint country process
n.sim <- 2e6      # Simulation length for cross-check
real.int.inf <- 'CPI'       # Can also use 'gdpdef'
save.file <- '~/Dropbox/data/2016/fiscalUnions/taxCoeffs.rdata'


## 1. Controls
x.bar <- 2/3
    # Average hours per person
g.expd <- FALSE
    # Set to TRUE to use total govt expenditure, otherwise uses government
    # consumption
alt.g <- FALSE
    # Alternative trend path for g

## 2. Read in the data
names(cts) <- cts

### 2.1 GDP per hour worked ###
df.nom.alp <- read_csv( alp.file )
names( df.nom.alp )[8] <- 'Year'
df.nom.alp <- subset( df.nom.alp, SUBJECT==alp.code[1] & MEASURE==alp.code[2] &
                        Country %in% cts & Year >= y.min )
df.nom.alp <- df.nom.alp[, c('Country', 'Subject', 'Measure', 'Year', 'Value')]
df.nom.alp$Year <- as.numeric(df.nom.alp$Year)
df.nom.alp$Value <- as.numeric(df.nom.alp$Value)
    # Drop needless columns
nom.alp <- reshape( df.nom.alp, timevar='Country', 
                    idvar=names(df.nom.alp)[2:(ncol(df.nom.alp)-1)], direction='wide' )
    # The nominal GDP dataframe
rownames(nom.alp) <- NULL
names(nom.alp)[ ncol(nom.alp)-0:1] <- sapply( names(nom.alp)[ ncol(nom.alp)-0:1], function(x) substr(x,7,nchar(x)) )
    # Adjust the names

### 2.3 Nominal GDP ###
df.nom.gdp <- read_csv( gdp.file )
df.nom.gdp <- subset( df.nom.gdp, TRANSACT==nom.gdp.code[1] & MEASURE==nom.gdp.code[2] &
                        Country %in% cts & Year >= y.min )
df.nom.gdp <- df.nom.gdp[, c('Country', 'Transaction', 'Measure', 'Year', 'Unit', 'PowerCode', 'Value')]
df.nom.gdp$Year <- as.numeric(df.nom.gdp$Year)
df.nom.gdp$Value <- as.numeric(df.nom.gdp$Value)
    # Drop needless columns
nom.gdp <- reshape( df.nom.gdp, timevar='Country', 
                    idvar=names(df.nom.gdp)[2:(ncol(df.nom.gdp)-1)], direction='wide' )
    # The nominal GDP dataframe
rownames(nom.gdp) <- NULL
names(nom.gdp)[ ncol(nom.gdp)-0:1] <- sapply( names(nom.gdp)[ ncol(nom.gdp)-0:1], function(x) substr(x,7,nchar(x)) )
    # Adjust the names

### 2.4 GDP deflator ###
df.gdp.def <- read_csv( gdp.file )
df.gdp.def <- subset( df.gdp.def, TRANSACT==gdp.def.code[1] & MEASURE==gdp.def.code[2] &
                        Country %in% cts & Year >= y.min )
df.gdp.def <- df.gdp.def[, c('Country', 'Transaction', 'Measure', 'Year', 'Unit', 'PowerCode', 'Value')]
df.gdp.def$Year <- as.numeric(df.gdp.def$Year)
df.gdp.def$Value <- as.numeric(df.gdp.def$Value)
    # Drop needless columns
gdp.def <- reshape( df.gdp.def, timevar='Country', 
                    idvar=names(df.nom.gdp)[2:(ncol(df.nom.gdp)-1)], direction='wide' )
    # The GDP deflator dataframe
rownames(gdp.def) <- NULL
names(gdp.def)[ ncol(gdp.def)-0:1] <- sapply( names(gdp.def)[ ncol(gdp.def)-0:1], function(x) substr(x,7,nchar(x)) )
    # Adjust the names

### 2.5 Nominal expenditures ###
df.expd <- read_csv( tax.file )
df.expd <- subset( df.expd, TRANSACT==expd.code[1] & SECTOR==expd.code[2] &
                     Country %in% cts & Year >= y.min )
df.expd <- df.expd[, c('Country', 'Transaction', 'Sector', 'Measure', 'Year', 'Unit', 'PowerCode', 'Value')]
df.expd$Year <- as.numeric(df.expd$Year)
df.expd$Value <- as.numeric(df.expd$Value)
    # Drop needless columns
expd <- reshape( df.expd, timevar='Country', idvar=names(df.expd)[2:(ncol(df.expd)-1)], direction='wide' )
    # The nominal expenditure dataframe
rownames(expd) <- NULL
names(expd)[ ncol(expd)-0:1] <- sapply( names(expd)[ ncol(expd)-0:1], function(x) substr(x,7,nchar(x)) )
    # Adjust the names

### 2.6 Govt consumption expenditures ###
df.gcons <- read.table( cons.file, sep="\t", header = TRUE )
    # For some reason read_csv gives the wrong data back!!!
df.gcons <- df.gcons[ -( nrow(df.gcons) - 0:4 ), ]
df.gcons <- df.gcons[, -c(1,2,4) ]
names( df.gcons )[1] <- 'Country'
names( df.gcons )[-1] <- sapply( names( df.gcons )[-1], 
                                 function(x) substr(x,2,5) )
df.gcons <- subset( df.gcons, Country %in% cts )
df.gcons[,-1] <- apply(df.gcons[,-1],2,as.numeric)
    # Clean the data
df.gcons.melt <- melt(df.gcons, id='Country' ) 
names( df.gcons.melt )[2:3] <- c('Year', 'GovtCons')
for( i in 2:3 ) df.gcons.melt[,i] <- as.numeric( as.character( df.gcons.melt[,i ] ) )
    # Melt down the data
gcons <- reshape( subset( df.gcons.melt, Year >= y.min & Country %in% cts ), timevar='Country', 
                  idvar='Year', direction='wide' )
    # Reshape as required
rownames(gcons) <- NULL
names(gcons)[ ncol(gcons)-0:1] <- sapply( names(gcons)[ ncol(gcons)-0:1], function(x) substr(x,10,nchar(x)) )
gcons[,cts] <- gcons[,cts] / 1e6
    # Adjust the names

### 2.7 Population ###
df.pop <- read_csv( pop.file )
df.pop <- df.pop[ !is.na(df.pop[['Country Code']]), ]
df.pop <- df.pop[, -(2:4) ]
names( df.pop ) <- sapply( names( df.pop ), 
                           function(x) strsplit(x, ' ')[[1]][1] )
df.pop <- subset( df.pop, Country %in% cts )
    # Clean the data
pop.melt <- melt(df.pop, id='Country' ) 
names( pop.melt )[2:3] <- c('Year', 'Population')
for( i in 2:3 ) pop.melt[,i] <- as.numeric( as.character( pop.melt[,i ] ) )
    # Melt down the data
pop <- reshape( subset( pop.melt, Year >= y.min ), timevar='Country', 
                idvar='Year', direction='wide' )
    # Reshape as required
rownames(pop) <- NULL
names(pop)[ ncol(pop)-0:1] <- sapply( names(pop)[ ncol(pop)-0:1], function(x) substr(x,12,nchar(x)) )
    # Adjust the names

## 3. Create single-country VARs

### 3.1 Create and filter real ave Labor Prod'y
log.real.alp <- data.frame( Year=nom.alp$Year, log( nom.alp[,cts]/ gdp.def[,cts] * 100 ) )
    # Log real ave lab prod
real.alp.hp <- apply( log.real.alp[,cts], 2, hpfilter,
                      freq=6.25, type="lambda", drift=FALSE )
real.alp.trend <- data.frame( Year=nom.alp$Year, sapply( real.alp.hp, function(x) x$trend ) )
real.alp.cycle <- data.frame( Year=nom.alp$Year, sapply( real.alp.hp, function(x) x$cycle ) )
    # Trend and cycle ALP
par(mfrow=c(2,1))
with( log.real.alp, plot( Year, France, ylim=range(log.real.alp[,cts]), 
                  col='blue', lwd=2, type='l', ylab='Log labor productivity' ) )
with( log.real.alp, lines( Year, Germany, ylim=range(real.alp.trend[,cts]), 
                   col='red', lwd=2 ) )
with( real.alp.trend, lines( Year, Germany, col='red', lwd=2, lty=2 ) )
with( real.alp.trend, lines( Year, France, col='blue', lwd=2, lty=2 ) )
    # Quick plot
with( real.alp.cycle, plot( Year, France, ylim=range(real.alp.cycle[,cts]), 
                  col='blue', lwd=2, type='l', ylab='Log labor productivity (cycle)' ) )
with( real.alp.cycle, lines( Year, Germany, ylim=range(real.alp.cycle[,cts]), 
                           col='red', lwd=2 ) )
abline(h=0, lwd=.5)
par(mfrow=c(1,1))

### 3.2 Create and filter government expenditure
df.g <- if( g.expd ) expd else gcons
g.shr <- data.frame( Year = df.g$Year, df.g[,cts] / nom.gdp[,cts] )
    # Share of expenditure in 
with( g.shr, plot( Year, France, type='l', lwd=2, col='blue', 
       ylim=range(g.shr[,cts], na.rm = TRUE ), ylab='Govt share of GDP' ) )
with( g.shr, lines( Year, Germany, lwd=2, col='red' ) )
    # Quick plot
ave.g.shr <- apply(g.shr[,cts], 2, mean, na.rm=T )
    # Average share of G in GDP
A.bar <- 1 / ( 1 - x.bar )
g.bar <- ave.g.shr
    # So average output is unity and g is relative to GDP
real.g <- data.frame( Year = df.g$Year, df.g[,cts] / gdp.def[,cts] * 100 )
    # NB: SHOULD UPDATE TO ALLOW FOR DIFFERENT DEFLATOR FOR EXPENDITURE
real.g <- real.g[ apply(real.g,1, function(x) all(!is.na(x))), ]
log.real.g <- data.frame( Year=real.g$Year, log(real.g[,cts]) )
log.real.g.hp <- apply( log.real.g[,cts], 2, hpfilter,
                   freq=6.25, type="lambda", drift=FALSE )
real.g.trend <- data.frame( Year=real.g$Year, sapply( log.real.g.hp, function(x) x$trend ) )
real.g.cycle <- data.frame( Year=real.g$Year, sapply( log.real.g.hp, function(x) x$cycle ) )
    # Trend and cycle log real GDP
if( alt.g ){
  real.g.trend <- data.frame( Year=real.alp.trend$Year,
                  t( ( apply( log.real.g[,cts], 2, mean ) / 
                    apply( log.real.alp[,cts], 2, mean ) ) * t( real.alp.trend[,cts] ) ) )
  real.g.cycle <- data.frame( Year=log.real.g$Year,
                    log.real.g[,cts] - real.g.trend[1:nrow(log.real.g), cts] )
}
    # Alternative specification for g
par(mfrow=c(2,1))
with( log.real.g, plot( Year, France, ylim=range(log.real.g[,cts]), 
                          col='blue', lwd=2, type='l', ylab='Log govt spending' ) )
with( log.real.g, lines( Year, Germany, ylim=range(real.g.trend[,cts]), 
                           col='red', lwd=2 ) )
with( real.g.trend, lines( Year, Germany, col='red', lwd=2, lty=2 ) )
with( real.g.trend, lines( Year, France, col='blue', lwd=2, lty=2 ) )
with( real.g.cycle, plot( Year, France, ylim=range(real.g.cycle[,cts]), 
                            col='blue', lwd=2, type='l', ylab='Log govt spending (cycle)' ) )
with( real.g.cycle, lines( Year, Germany, ylim=range(real.g.cycle[,cts]), 
                             col='red', lwd=2 ) )
abline(h=0, lwd=.5)
par(mfrow=c(1,1))
    # Quick plot

n.pds <- min( sapply( list( real.g.cycle, real.alp.cycle ), nrow ) )
    # Number of periods

## 4. Create the VARs ##
var.A.g <- lapply( cts, function(cty) 
              VAR( cbind( alp=real.alp.cycle[1:n.pds,cty], 
                          g=real.g.cycle[1:n.pds,cty] ), type='none' ) )
    # The VAR (gives non-zero constants because one period is dropped when lagging)
cor.A.g <- sapply( cts, function(cty) cor( real.alp.cycle[1:n.pds,cty], 
                                            real.g.cycle[1:n.pds,cty] ) )
    # The correlation
cov.A.g <- lapply( cts, function(cty) var( cbind( real.alp.cycle[1:n.pds,cty], 
                                           real.g.cycle[1:n.pds,cty] ) ) )
    # The covariance
phi <- lapply( cts, function(cty) t( sapply( c('alp','g'),
                     function(vble) var.A.g[[cty]]$varresult[[vble]]$coefficients ) ) )
    # The matrix of auto-regressive coefficients
sig <- lapply( cts, function(cty) summary(var.A.g[[cty]])$covres )
    # The Covariance of the residuals
lam <- lapply( cts, function(cty) lyapunovEq( phi[[cty]], sig[[cty]] ) )
    # The covariance of the stationary distribution
X <- lapply( lam, function(x) rmvnorm( nn, sigma=x ) )
eps <- lapply( sig, function(x) rmvnorm( nn, sigma=x ) )
    # Random variables for integration
X.prime.e <- lapply( cts, function(cty) X[[cty]] %*% t(phi[[cty]]) )
    # Expected value of X prime
Z <- lapply( X, function(x) kmeans( x, n.Z )$centers )
    # The Z points
T.p <- lapply( cts, function(cty) trans_prob( phi[[cty]], X[[cty]], 
                                              eps[[cty]], Z[[cty]] ) )
    # The transition probabilities
T.vals <- lapply(cts, 
                 function(cty) t( c( A.bar, g.bar[cty] ) * t(exp(Z[[cty]])) ) )
    # Return to the level of the shocks (not logs)
ergodic <- lapply( T.p, function(x) (x %^% 100)[1,] )
    # The ergodic distribution
plot( T.vals[[1]], cex = ergodic[[1]] * 40, main="Ergodic distribution",
      xlim=range(sapply(T.vals, function(x) x[,1] ) ), col='red', lwd=2,
      ylim=range(sapply(T.vals, function(x) x[,2] ) ), xlab='A', ylab='g' ) 
points( T.vals[[2]], cex = ergodic[[2]] * 40, col='blue', lwd=2 ) 
    # Plot the ergodic distribution
for(cty in cts) heatmap(T.p[[cty]], Colv ="Rowv", main=cty )
    # Heatmap of transition probabilities

## 5. Simulate the discretized process and compare to the VAR ##
mc.indiv <- lapply( cts, function(cty) 
  new( "markovchain", byrow=TRUE, transitionMatrix=T.p[[cty]], name=cty ) )
    # The markov chain object
set.seed(1234)
sim.indiv <- lapply( cts, function(cty)
         T.vals[[cty]][as.numeric(rmarkovchain( n.sim, mc.indiv[[cty]], 
                                                useRCpp=TRUE ) ), ] )
for( cty in cts ) rownames(sim.indiv[[cty]]) <- NULL
for( cty in cts ) colnames(sim.indiv[[cty]]) <- c('A', 'g') 
    # Simulated series

### 5.1 Time series plot ###
n.pds <- nrow( real.alp.cycle )
with( real.alp.cycle, plot( 1:n.pds, France, ylim=range(real.alp.cycle[,cts]), 
                            col='blue', lwd=2, type='l', ylab='Log labor productivity (cycle)' ) )
with( real.alp.cycle, lines( 1:n.pds, Germany, ylim=range(real.alp.cycle[,cts]), 
                             col='red', lwd=2 ) )
for( i in 1:2) 
  lines( 1:n.pds, log( sim.indiv[[cts[i]]][1:n.pds,'A'] / A.bar ), lwd=2, lty=2, 
         col=c('red','blue')[i] )
abline( h=0, lwd=.5 )
legend( 'bottomright', cts, lwd=2, col=c('red','blue'), bty='n' )
    # Simulation
n.pds <- nrow(real.g.cycle)
with( real.g.cycle, plot( 1:n.pds, France, ylim=range(real.g.cycle[,cts]), 
                            col='blue', lwd=2, type='l', ylab='Log labor productivity (cycle)' ) )
with( real.g.cycle, lines( 1:n.pds, Germany, ylim=range(real.g.cycle[,cts]), 
                             col='red', lwd=2 ) )
for( i in 1:2) 
  lines( 1:n.pds, log( sim.indiv[[cts[i]]][1:n.pds,'g'] / g.bar[i] ), lwd=2, lty=2, 
         col=c('red','blue')[i] )
abline( h=0, lwd=.5 )
legend( 'bottomright', cts, lwd=2, col=c('red','blue'), bty='n' )

cov.mc <- lapply( sim.indiv, function(x) var(log(x)) )
    # simulated covariance

var.mc <- lapply( sim.indiv, function(x) VAR( log( x ), type='none' ) )
    # Create the VAR on the simulated data
par(mfrow = c(1, 2))
for( cty in cts ){
  barplot( rbind( cov.A.g[[cty]][c(1,2,4)], 
                  cov.mc[[cty]][c(1,2,4)] ), 
           beside=TRUE, col=c('blue', 'red'), main=paste0(cty, ' covariances'),
           legend.text=c( 'Data', 'Markov Chain \nSimulation' ),
           args.legend =list(x='topleft', bty='n') )
}
for( cty in cts ){
  barplot( rbind( c( sapply( var.A.g[[cty]]$varresult, function(x) x$coefficients ) ),
                  c( sapply( var.mc[[cty]]$varresult, function(x) x$coefficients ) ) ), 
           beside=TRUE, col=c('blue', 'red'), main=paste0(cty, ' VAR coefficients'),
           legend.text=c('Data','Markov Chain \nSimulation'), 
           args.legend =list(x='bottomleft', bty='n') ) 
}
par(mfrow = c(1, 1))


## 6. Create the joint VAR ##
real.alp.g.cycle <- cbind( Year = real.g.cycle$Year,
                           alp=real.alp.cycle[1:nrow(real.g.cycle),cts], 
                           g=real.g.cycle[, cts ])
var.joint <- VAR( real.alp.g.cycle[,-1], type='none' )
    # The VAR (gives non-zero constants because one period is dropped when lagging)
cov.joint <- var( real.alp.g.cycle[,-1] )
    # The covariance
phi.joint <- t(sapply( var.joint$varresult, function(x) x$coefficients ))
    # The matrix of auto-regressive coefficients
sig.joint <- summary(var.joint)$covres
    # The covariance of the residuals
lam.joint <- lyapunovEq( phi.joint, sig.joint )
    # The covariance of the stationary distribution
X <- rmvnorm( nn, sigma=lam.joint )
eps <- rmvnorm( nn, sigma=sig.joint )
    # Random variables for integration
X.prime.e <- X %*% t(phi.joint)
    # Expected value of X prime
Z <- kmeans( X, n.Z.joint )$centers
    # The Z points
T.p.joint <- trans_prob( phi.joint, X, eps, Z )
    # The transition probabilities
T.vals.joint <- exp(Z) * ( rep(1,nrow(Z)) %*% t(c( rep( A.bar, length(cts)), g.bar ) ) )
    # Return to the level of the shocks (not logs).  NO SCALING HERE!!!
colnames( T.vals.joint ) <- paste0( rep( c("A.","g."), each=length(cts)), cts )
    # Label the columns
ergodic.joint <- ( T.p.joint %^% 100)[1,]
    # The ergodic distribution
plot( T.vals.joint[,c(1,3)], cex = ergodic.joint * 80, main="Ergodic distribution",
      xlim=range( T.vals.joint[,1:2] ), col='red', lwd=2,
      ylim=range( T.vals.joint[,3:4] ), xlab='A', ylab='g' ) 
points( T.vals.joint[,c(2,4)], cex = ergodic.joint * 80, col='blue', lwd=2 ) 
    # Plot the ergodic distribution
heatmap(T.p.joint, Colv ="Rowv", main='Joint transition probabilities' )
    # Heatmap of transition probabilities

## 7. Simulate the discretized process and compare to the VAR ##
mc.joint <- new( "markovchain", byrow=TRUE, transitionMatrix=T.p.joint, name="Joint" )
    # The markov chain object
set.seed(1234)
sim.joint <- T.vals.joint[as.numeric(rmarkovchain( n.sim, mc.joint, useRCpp=TRUE ) ), ]
rownames(sim.joint) <- NULL
colnames(sim.joint) <- colnames( T.vals.joint ) 
    # Simulated taxes

### 7.1 Plot the covariance ###
cov.mc.joint <- list( Germany=cov(log(sim.joint[,c(1,3)])), 
                      France=cov(log(sim.joint[,c(2,4)])) )
par(mfrow = c(1, 2))
for( cty in cts ){
  barplot( rbind( cov.A.g[[cty]][c(1,2,4)], 
                  cov.mc.joint[[cty]][c(1,2,4)] ), 
           beside=TRUE, col=c('blue', 'red'), main=paste0(cty, ' covariances'),
           legend.text=c( 'Data', 'Markov Chain \nSimulation (joint)' ),
           args.legend =list(x='topleft', bty='n') )
}
par(mfrow=c(1,1))


### 8. Scaling ###
rho <- mean( pop[,cts[1]] / ( apply( pop[,cts], 1, sum ) ) )
    # Population shares
gdp.shr <- data.frame( nom.gdp$Year, nom.gdp[,cts] / apply( nom.gdp[,cts],1, sum ) )
    # Share of total GDP
ave.gdp.shr <- apply( gdp.shr[,cts], 2, mean )
    # Average shares
A.bar.joint <- 1/ ( 1 - x.bar ) * ave.gdp.shr / c( rho, 1-rho )
g.bar.joint <- ave.gdp.shr * ave.g.shr
g.bar.joint.chk <- g.bar * A.bar.joint / A.bar * c(rho, 1-rho )
if(!all( g.bar.joint == g.bar.joint.chk ) ) warning( 'g.bar.joint check fails' )
    # Scaling
T.vals.joint <- exp(Z) * ( rep(1,nrow(Z)) %*% t(c( A.bar.joint, g.bar.joint ) ) )

## TO ADD: SCALING AND SAVING ###