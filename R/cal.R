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

Rcpp::sourceCpp('code/2016/fiscalUnions/R/rcpp/dists.cpp')
source('code/2016/fiscalUnions/R/AR1disc.R')

## 0. Script variables
gdp.file <- '~/Dropbox/data/2016/fiscalUnions/gdpDefANA.csv'
    # OECD annual GDP data
hours.pw.file <- '~/Dropbox/data/2016/fiscalUnions/hrsPerWorker.csv'
    # Annual hours per worker
emp.file <- '~/Dropbox/data/2016/fiscalUnions/totalEmployment.csv'
    # Total employment
exp.file <- '~/Dropbox/data/2016/fiscalUnions/govtNomTaxANA.csv'
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
cts <- c("Germany", "France")
expd.code <- c( 'GTE', 'GS13')         # Total expenditure, general government
nom.gdp.code <- c( 'B1_GA', 'C' )     # GDP, current prices
gdp.def.code <- c( 'B1_GA', 'DOB' )   # GDP, deflator (could use expenditure deflator here)

y.min <- 1980     # 1960
nn <- 2000        # Simulation periods
n.Z <- 30         # Number of discretized tax levels
n.ar1 <- 10        # Number of AR(1) points
n.sim <- 2e6      # Simulation length for cross-check
real.int.inf <- 'CPI'       # Can also use 'gdpdef'
save.file <- '~/Dropbox/data/2016/fiscalUnions/taxCoeffs.rdata'


## 1. Controls
x.bar <- 2/3
    # Average hours per person
exp <- TRUE
    # Set to TRUE to use total govt expenditure, otherwise uses government
    # consumption


## 2. Read in the data
names(cts) <- cts

### 2.1 Hours per worker ###
cts.hrs <- if( 'Germany' %in% cts ) c(cts, 'West Germany'='West Germany' ) else cts
df.hrs <- read_csv( hours.pw.file )
names( df.hrs )[8] <- 'Year'
df.hrs$Year <- as.numeric(as.character(df.hrs$Year))
df.hrs <- subset( df.hrs, Country %in% cts.hrs & Year >= y.min )
df.hrs <- df.hrs[, c('Country', 'Employment status', 'Year', 'Unit', 'Value')]
df.hrs$Value <- as.numeric(df.hrs$Value)
    # Drop needless columns
hrs <- reshape( df.hrs, timevar='Country', 
                    idvar=names(df.hrs)[2:(ncol(df.hrs)-1)], direction='wide' )
    # The nominal GDP dataframe
rownames(hrs) <- NULL
names(hrs)[ 4:ncol(hrs)] <- sapply( names(hrs)[ 4:ncol(hrs)], function(x) substr(x,7,nchar(x)) )
    # Adjust the names
if( 'Germany' %in% cts ){
  hrs.ratio <- mean( hrs$Germany / hrs$`West Germany`, na.rm=TRUE )
  hrs$Germany[is.na(hrs$Germany)] <- hrs.ratio * hrs$`West Germany`[is.na(hrs$Germany)]
      # Paste in the W German hours
  with( hrs, plot( Year, France, type='l', lwd=2, col='red', ylab='Hours/worker',
                   ylim=range( hrs[,4:6], na.rm = T ) ) )
  with( hrs, lines( Year, Germany, type='l', lwd=2, col='blue' ) )
  with( hrs, lines( Year, `West Germany`, type='l', lwd=2, col='blue', lty=2 ) )
      # Quick plot
}


### 2.2 Nominal GDP ###
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

### 2.3 GDP deflator ###
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