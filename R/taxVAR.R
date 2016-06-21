library(readr)
library(vars)
library(mvtnorm)
library(markovchain)
library(expm)
# library(tidyr)
# library(dplyr)
# library(lubridate)

Rcpp::sourceCpp('code/2016/fiscalUnions/R/rcpp/dists.cpp')

## 0. Script variables
file <- '~/Dropbox/data/2016/fiscalUnions/taxRevenue.csv'
cts <- c("Germany", "France")
gov <- 'NES'      # Total tax revenue
y.min <- 1980     # 1960
nn <- 1000        # Simulation periods
n.Z <- 20         # Number of discretized periods
n.sim <- 1e8   # Simulation length for cross-check

## 1. Read the data
df <- read_csv( file )
df <- subset( df, GOV==gov & Country %in% cts & Tax == 'Total tax revenue' )
df <- df[, c('Government', 'Tax', 'Variable', 'Country', 'Year', 'Unit', 'Value')]
df[ df$Value > 100, ]$Value <- df[ df$Value > 100, ]$Value / 1000
    # Clean the data (some values mis-recorded)
df$Year <- as.numeric(df$Year)
df$Value <- as.numeric(df$Value)
tax <- reshape( df, timevar='Country', idvar=c('Year', 'Government', 'Tax', 'Variable', 'Unit'), 
                  direction='wide' )
    # The taxes dataframe
names(tax)[ ncol(tax)-0:1] <- sapply( names(tax)[ ncol(tax)-0:1], function(x) substr(x,7,nchar(x)) )
    # Adjust the names

## 2. A quick plot ##
plot( tax$Year[ tax$Year >= y.min ], tax[tax$Year >= y.min ,cts[1]], type='l', ylim = c(30,45), col='blue', lwd=2, 
      xlim=c(y.min,max(tax$Year)), xlab='Year', ylab='Tax/GDP' )
lines( tax$Year[tax$Year >= y.min ], tax[tax$Year >= y.min ,cts[2]], type='l', col='red', lwd=2 )
legend( 'bottomright', cts, lwd=2, col=c('blue','red'), bty='n' )
var.tax <- VAR( subset(tax, Year > y.min)[,cts] )

## 3. The demeaned version ##
tax.dm <- tax
tax.means <- apply( subset(tax, Year > y.min)[,cts], 2, mean )
tax.dm[,cts] <- tax[,cts] - rep( tax.means, each=nrow(tax) )
var.tax.dm <- VAR( subset(tax.dm, Year > y.min)[,cts],type='none'  )
    # Need to remove the constant because regress on *lagged* sample so mean of
    # regressors not quite zero
phi <- t(sapply( 1:length(cts), function(i) var.tax.dm$varresult[[cts[i]]]$coefficients ))
    # The matrix of auto-regressive coefficients
sig <- summary(var.tax.dm)$covres
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

## 4. Return to the version with means ##
T.vals <- Z + rep( tax.means, each=n.Z )
    # Re-mean
ergodic <- (T.p %^% 100)[1,]
    # The ergodic distribution
plot( T.vals, cex = ergodic * 10, main="Ergodic distribution" )
    # Plot the ergodic distribution
heatmap(T.p, Colv ="Rowv" )
    # Heatmap of transition probabilities

## 5. Simulate the discretized process and compare to the VAR ##
mc.tax <- new( "markovchain", byrow=TRUE, transitionMatrix=T.p, name="Tax" )
    # The markov chain object
sim.tax <- T.vals[as.numeric(rmarkovchain( n.sim, mc.tax ) ), ]
rownames(sim.tax) <- NULL
colnames(sim.tax) <- cts
    # Simulated taxes
plot( 1:40, sim.tax[ 1:40, 1 ], type='l', col='blue', ylim = c(30,45),
      xlab='Sim pd', ylab='Tax/GDP', lwd=2 )
lines( 1:40, sim.tax[ 1:40, 2 ], type='l', col='red', lwd=2 )
legend( 'bottomright', cts, lwd=2, col=c('blue','red'), bty='n' )
    # Quick plot
var.mc <- VAR( sim.tax )
    # Create the VAR on the simulated data
par(mfrow = c(1, 2))
for( cty in cts ){
  barplot( rbind( var.mc$varresult[[cty]]$coefficients, 
                var.tax$varresult[[cty]]$coefficients ), 
         beside=TRUE, col=c('blue', 'red'), main=paste0(cty, ' VAR coefficients'),
         legend.text=c('Markov Chain \nSimulation', 'Data'), 
         args.legend =list(x='topleft',bty='n'), log='y' )
}
par(mfrow = c(1, 1))

### Q: Can I use markovchainFit on the nearest-neighbour version of a simulation
### from the VAR to generate an answer?



