library(readr)
library(vars)
library(mvtnorm)
# library(tidyr)
# library(dplyr)
# library(lubridate)

# Rcpp::sourceCpp('code/2016/fiscalUnions/R/rcpp/dists.cpp')

## 0. Script variables
file <- '~/Dropbox/data/2016/fiscalUnions/taxRevenue.csv'
cts <- c("Germany", "France")
gov <- 'NES'  # Total tax revenue
y.min <- 1980 # 1960
nn <- 1000   # Simulation periods
n.Z <- 10     # Number of discretized periods

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
plot( tax$Year, tax[,cts[1]], type='l', ylim = c(30,45), col='blue', lwd=2 )
lines( tax$Year, tax[,cts[2]], type='l', col='red', lwd=2 )
var.tax <- VAR( subset(tax, Year > y.min)[,cts] )

## 3. The demeaned version ##
tax.dm <- tax
tax.means <- rep( apply( subset(tax, Year > y.min)[,cts], 2, mean ), each=nrow(tax) )
tax.dm[,cts] <- tax[,cts] - tax.means
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
T.vals <- Z + tax.means
    # Re-mean

# X prime itself needs to be an array of depth number of countries

