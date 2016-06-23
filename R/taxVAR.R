#############################################################################
# taxVAR.R
#
# Script to calibrate the fiscal unions model
# Philip Barrett, Chicago
# Created: 10jun2016
#############################################################################
rm(list=ls())
library(readr)
library(vars)
library(mvtnorm)
library(markovchain)
library(expm)
library(reshape)
# library(tidyr)
# library(plyr)
# library(lubridate)

Rcpp::sourceCpp('code/2016/fiscalUnions/R/rcpp/dists.cpp')
source('code/2016/fiscalUnions/R/AR1disc.R')

## 0. Script variables
tax.file <- '~/Dropbox/data/2016/fiscalUnions/govtNomTaxANA.csv'
    # This data source (OECD ANAs) gives consistent measures of government
    # revenue, expenditure, and deficits.  The only downside is that it has no
    # German data before 1995.  For this, I need to stitch on the canned
    # revenue/GDP from the OECD taxation database
gdp.file <- '~/Dropbox/data/2016/fiscalUnions/gdpDefANA.csv'
    # From OECD annual GDP data
check.file <- '~/Dropbox/data/2016/fiscalUnions/taxRevenue.csv'
    # From the OECD taxation database
pop.file <- '~/Dropbox/data/2016/fiscalUnions/population.csv'
    # Population data
cons.file <- '~/Dropbox/data/2016/fiscalUnions/govtConsExpd.csv'
    # Government consumption expenditure data (required for pre-1995 Germany)
int.file <- '~/Dropbox/data/2016/fiscalUnions/interestRates.csv'
    # Interest rates
debt.file <- '~/Dropbox/data/2016/fiscalUnions/govtDebt.csv'
    # Debt/GDP ratios
cts <- c("Germany", "France")
tax.code <- c( 'GTR', 'GS13')         # Total revenue, general government
expd.code <- c( 'GTE', 'GS13')         # Total expenditure, general government
nom.gdp.code <- c( 'B1_GA', 'C' )     # GDP, current prices
gdp.def.code <- c( 'B1_GA', 'DOB' )   # GDP, deflator (could use expenditure deflator here)
tax.gdp.code <- 'NES'                 # Total tax revenue
y.min <- 1980     # 1960
nn <- 2000        # Simulation periods
n.Z <- 20         # Number of discretized periods
n.ar1 <- 10        # Number of AR(1) points
n.sim <- 2e6      # Simulation length for cross-check
real.int.inf <- 'CPI'       # Can also use 'gdpdef'
save.file <- '~/Dropbox/data/2016/fiscalUnions/taxCoeffs.rdata'

## 1. Read the data
names(cts) <- cts

### 1.1 Nominal taxes ###
df.tax <- read_csv( tax.file )
df.tax <- subset( df.tax, TRANSACT==tax.code[1] & SECTOR==tax.code[2] &
                    Country %in% cts & Year >= y.min )
df.tax <- df.tax[, c('Country', 'Transaction', 'Sector', 'Measure', 'Year', 'Unit', 'PowerCode', 'Value')]
df.tax$Year <- as.numeric(df.tax$Year)
df.tax$Value <- as.numeric(df.tax$Value)
    # Drop needless columns
tax <- reshape( df.tax, timevar='Country', idvar=names(df.tax)[2:(ncol(df.tax)-1)], direction='wide' )
    # The nominal tax dataframe
rownames(tax) <- NULL
names(tax)[ ncol(tax)-0:1] <- sapply( names(tax)[ ncol(tax)-0:1], function(x) substr(x,7,nchar(x)) )
    # Adjust the names

### 1.2 Nominal GDP ###
df.nom.gdp <- read_csv( gdp.file )
df.nom.gdp <- subset( df.nom.gdp, TRANSACT==nom.gdp.code[1] & MEASURE==nom.gdp.code[2] &
                      Country %in% cts & Year >= y.min )
df.nom.gdp <- df.nom.gdp[, c('Country', 'Transaction', 'Measure', 'Year', 'Unit', 'PowerCode', 'Value')]
df.nom.gdp$Year <- as.numeric(df.nom.gdp$Year)
df.nom.gdp$Value <- as.numeric(df.nom.gdp$Value)
    # Drop needless columns
nom.gdp <- reshape( df.nom.gdp, timevar='Country', 
                    idvar=names(df.nom.gdp)[2:(ncol(df.nom.gdp)-1)], direction='wide' )
    # The nominal tax dataframe
rownames(nom.gdp) <- NULL
names(nom.gdp)[ ncol(nom.gdp)-0:1] <- sapply( names(nom.gdp)[ ncol(nom.gdp)-0:1], function(x) substr(x,7,nchar(x)) )
    # Adjust the names

### 1.3 GDP deflator ###
df.gdp.def <- read_csv( gdp.file )
df.gdp.def <- subset( df.gdp.def, TRANSACT==gdp.def.code[1] & MEASURE==gdp.def.code[2] &
                        Country %in% cts & Year >= y.min )
df.gdp.def <- df.gdp.def[, c('Country', 'Transaction', 'Measure', 'Year', 'Unit', 'PowerCode', 'Value')]
df.gdp.def$Year <- as.numeric(df.gdp.def$Year)
df.gdp.def$Value <- as.numeric(df.gdp.def$Value)
    # Drop needless columns
gdp.def <- reshape( df.gdp.def, timevar='Country', 
                    idvar=names(df.nom.gdp)[2:(ncol(df.nom.gdp)-1)], direction='wide' )
    # The nominal tax dataframe
rownames(gdp.def) <- NULL
names(gdp.def)[ ncol(gdp.def)-0:1] <- sapply( names(gdp.def)[ ncol(gdp.def)-0:1], function(x) substr(x,7,nchar(x)) )
    # Adjust the names

### 1.4 Tax/GDP ratios ###
df.tax.gdp <- read_csv( check.file )
df.tax.gdp <- subset( df.tax.gdp, GOV==tax.gdp.code & Country %in% cts & 
                        Tax == 'Total tax revenue' & Year >= y.min )
df.tax.gdp <- df.tax.gdp[, c('Government', 'Tax', 'Variable', 'Country', 'Year', 'Unit', 'Value')]
df.tax.gdp[ df.tax.gdp$Value > 100, ]$Value <- df.tax.gdp[ df.tax.gdp$Value > 100, ]$Value / 1000
    # Clean the data (some values mis-recorded)
df.tax.gdp$Year <- as.numeric(df.tax.gdp$Year)
df.tax.gdp$Value <- as.numeric(df.tax.gdp$Value)
tax.gdp <- reshape( df.tax.gdp, timevar='Country', 
                    idvar=c('Year', 'Government', 'Tax', 'Variable', 'Unit'), 
                    direction='wide' )
    # The taxes dataframe
names(tax.gdp)[ ncol(tax.gdp)-0:1] <- sapply( names(tax.gdp)[ ncol(tax.gdp)-0:1], function(x) substr(x,7,nchar(x)) )
    # Adjust the names

### 1.5 Population ###
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

### 1.6 Nominal expenditures ###
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

### 1.7 Govt consumption expenditures ###
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

### 1.8 Real interest rates ###
df.int <- read_csv( int.file )
names(df.int)[-1] <- c( 'Interest.France', 'Interest.Germany', 
                        'CPI.France', 'gdpdef.France',
                        'CPI.Germany', 'gdpdef.Germany' )
df.int.a <- aggregate( df.int[,-1], list(format(df.int$DATE, '%Y' )), mean )
names(df.int.a)[1] <- 'Year'
    # Annual averages
real.int <- data.frame( Year=df.int.a$Year, 
                        sapply(cts, function(cty) 
                          df.int.a[[paste0('Interest.', cty)]] - 
                            df.int.a[[paste0( real.int.inf, '.', cty)]] ) )
for( i in 1:ncol(real.int) ) real.int[,i] <- as.numeric( as.character( real.int[,i] ) )
real.int <- subset( real.int, Year >= y.min )
    # Realized real rates 

### 1.9 Debt levels ###
debt.gdp <- read_csv( debt.file )
names( debt.gdp ) <- c( 'Year', 'France', 'Germany' )
debt.gdp[ debt.gdp==0 ] <- NA
debt.gdp$Year <- format( debt.gdp$Year, '%Y' )
for( i in 1:ncol(debt.gdp) ) debt.gdp[,i] <- as.numeric(debt.gdp[,i])
debt.gdp <- subset( debt.gdp, Year >= y.min )
    # The debt/GDP ratios

## 2. Compare the canned and home-made tax/GDP ratios (as a check) ##
tax.gdp.hand <- data.frame( Year = tax$Year, tax[,cts] / nom.gdp[,cts] * 100 )
    # Create by hand
plot( tax.gdp$Year, tax.gdp[, cts[1]] / tax.gdp.hand[ 1:nrow(tax.gdp),cts[1]], type='l', col='blue', lwd=2, 
      xlab='Year', ylab='Ratio of measures', ylim=c(.75,.9) )
lines( tax.gdp$Year, tax.gdp[, cts[2]] / tax.gdp.hand[ 1:nrow(tax.gdp),cts[2]], col='red', lwd=2 )
    # Tax plot
plot( expd$Year, gcons[, cts[1]] / expd[ , cts[1]], type='l', col='blue', lwd=2, 
      xlab='Year', ylab='Gov exp/cons ratio', ylim=c(0,1) )
lines( expd$Year, gcons[, cts[2]] / expd[ ,cts[2]], col='red', lwd=2 )
    # Expenditure/consumption plot. Looks bad for France, but fine for Germany. Why?

### 2.1 If Germany in cts, stitch on the OECD tax database data
if( 'Germany' %in% cts ){
  scale.factor.tax <- mean( tax.gdp.hand$Germany[1:nrow(tax.gdp)] / tax.gdp$Germany, 
                        na.rm = TRUE ) / 100
  replace.tax <- is.na( tax$Germany )
  scale.factor.cons <- mean( gcons$Germany / expd$Germany, na.rm = TRUE )
  replace.tax <- is.na( tax$Germany )
  replace.cons <- is.na( expd$Germany )
      # The rows to replace
  tax$Germany[replace.tax] <- scale.factor.tax * tax.gdp$Germany[replace.tax] * nom.gdp$Germany[replace.tax]
  expd$Germany[replace.cons] <- gcons$Germany[replace.cons] / scale.factor.cons
}

## 2.2 replot ##
plot( tax.gdp$Year, tax.gdp[, cts[1]], type='l', ylim = c(30,55), col='blue', lwd=2, 
      xlab='Year', ylab='Tax/GDP' )
lines( tax.gdp$Year, tax.gdp[ ,cts[2]], type='l', col='red', lwd=2 )
for(i in 1:2){
  lines( tax$Year, tax[ ,cts[i]] / nom.gdp[ ,cts[i]] * 100, 
         type='l', col=c('blue','red')[i], lwd=2, lty=2 )
}
legend( 'topleft', paste( rep( cts, 2), rep( c( 'OECD', 'Constructed'), each=2 ) ), 
        lwd=2, col=c('blue','red'), lty=c(1,1,2,2), bty='n' )
    # The two series


## 3. Create the real tax series and measure the trend growth rates ##
real.tax <- data.frame( Year=tax$Year, tax[,cts] / ( gdp.def[,cts] / 100 ) )
    # In 2014 Euros
lm.tax.ratio <- lm( log( Germany / France ) ~ Year, data=real.tax )
l.lm.tax.cts <- lapply( cts, function(cty) lm( paste( 'log(', cty, ' ) ~ Year' ), 
                                                      data=real.tax ) )
    # Measuring the relative growth rates of the two real tax series
gam <- exp( mean( sapply( l.lm.tax.cts, function(x) x$coefficients['Year'] ) ) ) - 1
    # The mean annual growth rate of taxes.  Require the exp and subtraction to
    # convert from log levels to annual growth rates (v. small diff though).
tax.dt <- data.frame( Year=tax$Year, log( real.tax[,cts] ) -
                      sapply( cts, function(cty) predict( l.lm.tax.cts[[cty]] ) ) )
    # Real detrended taxes


## 4. The demeaned version ##
var.tax <- VAR( tax.dt[,cts], type='none' )
    # The VAR (gives non-zero constants because one period is dropped when lagging)
phi <- t(sapply( 1:length(cts), function(i) var.tax$varresult[[cts[i]]]$coefficients ))
    # The matrix of auto-regressive coefficients
sig <- summary(var.tax)$covres
    # The Covariance of the residuals
lam <- lyapunovEq( phi, sig )
    # The covariance of the stationary distribution
X <- rmvnorm( nn, sigma=lam )
eps <- rmvnorm( nn, sigma=sig )
    # Random variables for integration
X.prime.e <- X %*% t(phi)
    # Expected value of X prime
Z <- kmeans( X, n.Z )$centers
    # The Z points
T.p <- trans_prob( phi, X, eps, Z )
    # The transition probabilities
T.vals <- exp(Z)
    # Return to the true growth rate version (not logs)
ergodic <- (T.p %^% 100)[1,]
    # The ergodic distribution
plot( T.vals, cex = ergodic * 40, main="Ergodic distribution" )
    # Plot the ergodic distribution
heatmap(T.p, Colv ="Rowv" )
    # Heatmap of transition probabilities


## 5. Simulate the discretized process and compare to the VAR ##
mc.tax <- new( "markovchain", byrow=TRUE, transitionMatrix=T.p, name="Tax" )
    # The markov chain object
set.seed(1234)
sim.tax <- T.vals[as.numeric(rmarkovchain( n.sim, mc.tax, useRCpp=TRUE ) ), ]
rownames(sim.tax) <- NULL
colnames(sim.tax) <- cts
    # Simulated taxes

### 5.1 Time series plot ###
n.pds <- max( tax.dt$Year ) - y.min + 1
par(mfrow = c(1, 2))
  plot( tax.dt$Year, tax.dt[, cts[1]], type='l', col='blue', lwd=2, main='Data',
        xlab='Year', ylab='Detrended real taxes' )
  lines( tax.dt$Year,  tax.dt[, cts[2]], col='red', lwd=2 )
  abline( h=0, lwd=.5 )
  legend( 'bottomright', cts, lwd=2, col=c('blue','red'), bty='n' )
      # Data
  plot( 1:n.pds, log( sim.tax[ 1:n.pds, 1 ] ), type='l', col='blue', ylim=range( tax.dt[, cts]),
        xlab='Sim pd', ylab='Detrended real taxes', lwd=2, main='Simulation' )
  lines( 1:n.pds, log( sim.tax[ 1:n.pds, 2 ] ), type='l', col='red', lwd=2 )
  abline( h=0, lwd=.5 )
  legend( 'bottomright', cts, lwd=2, col=c('blue','red'), bty='n' )
      # Simulation
par(mfrow = c(1, 1))

var.mc <- VAR( log( sim.tax ), type='none' )
    # Create the VAR on the simulated data
par(mfrow = c(1, 2))
for( cty in cts ){
  barplot( rbind( var.mc$varresult[[cty]]$coefficients, 
                var.tax$varresult[[cty]]$coefficients ), 
         beside=TRUE, col=c('blue', 'red'), main=paste0(cty, ' VAR coefficients'),
         legend.text=c('Markov Chain \nSimulation', 'Data'), 
         args.legend =list(x=if(cty=='Germany') 'topright' else 'topleft', bty='n') )
}
par(mfrow = c(1, 1))


## 6. Individual countries ##
l.ar1 <- lapply( cts, function(cty) ar(tax.dt[[cty]], order.max=1 ) )
    # The list of AR(1) autoregressions
l.indiv <- lapply( cts, function(cty) ar1.disc( n.ar1, l.ar1[[cty]]$ar, 
                            sqrt(l.ar1$Germany$var.pred ), ext=FALSE) )
    # The individual countries' discretized processes


## 7. Population growth and per capita real government spending ##
l.lm.pop.cts <- lapply( cts, function(cty) lm( paste( 'log(', cty, ' ) ~ Year' ), 
                                               data=pop ) )
    # Mesauring the average growth rates of the countries' populations
nn <- exp( mean( sapply( l.lm.pop.cts, function(x) x$coefficients['Year'] ) ) ) - 1
    # The mean annual growth rate of population  Require the exp and subtraction to
    # convert from log levels to annual growth rates (v. small diff though).
real.expd <- data.frame( Year=expd$Year, expd[,cts] / ( gdp.def[,cts] / 100 ) )
    # Real government expenditure
real.expd.pc <- data.frame( Year=expd$Year, real.expd[,cts] / pop[,cts] * 1e6 )
    # Real per capita expenditure.  Useful for a quick plot
real.expd.t <- data.frame( Year=expd$Year, real.expd[,cts] / 
                             exp( sapply( l.lm.tax.cts, predict ) ) )
    # Expenditure relative to trend taxes
# sd.expd <- apply( real.expd.t[,-1], 2, sd, na.rm=T )
sd.expd.t <- apply( subset(real.expd.t, Year != 1995 )[,-1], 2, sd, na.rm=T )
    # The standard deviation of expenditure.  Drop 1995 as it looks wrong.
    # Use this for calibrating sigma.  (Relative to trend taxes)
mu.expd.t <- apply( subset(real.expd.t, Year != 1995 )[,-1], 2, mean, na.rm=T )
    # Mean expenditure.  Use for calibrating beta (equivalent to avg defecit) 
    # Wait, really?  Doesn't g > T imply that debt is negative in the long run? 
    # Or that gamma > r, right?
with( real.expd.t, plot( Year, France, lwd=2, col='red', typ='l', ylim=c(.9,1.2 ),
                         xlab='Year', ylab='Expenditure to tax trend ratio' ) )
with( real.expd.t, lines( Year, Germany, lwd=2, col='blue', typ='l' ) )
abline( h=mu.epxd, col=c('red', 'blue'), lty=2)
    # Plot the data to have a quick look at the two variances

## 8. Interest rates ##
rr <- mean(as.matrix(real.int[,cts])) / 100
    # The average real interest rate
r.diff <- real.int[,cts[1]] - real.int[,cts[2]]
rr.diff.ave <- c( abs = mean(abs(r.diff)), mean=mean(r.diff), 
                  median=median(r.diff) )
    # The average difference in the real interest rate

## 9. Debt levels ##
debt.t <- data.frame( Year=debt.gdp$Year,
                      ( debt.gdp[,cts] * nom.gdp[,cts] / gdp.def[,cts] ) / 
                        exp( sapply( l.lm.tax.cts, predict ) ) )
    # Debt/Tax trend ratio
# Try to hit these on average (kind of dumb)

## 10. Save the results ##
T.vals.temp <- matrix( 0, nrow(T.vals), ncol(T.vals) )
    # For some reason, reading T.vals directly in julia gives an error
for( i in 1:length(T.vals) ) T.vals.temp[i] <- T.vals[i]
    # Fill the values by hand
save( l.indiv, sd.expd, mu.expd, gam, nn, rr, T.vals.temp, T.p, n.Z, n.ar1, file = save.file )  # l.indiv
