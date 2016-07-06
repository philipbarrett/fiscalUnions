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
cts <- c("Germany", "France")
alp.code <- c( 'T_GDPHRS', 'C')       # GDP per hour worked, current prices
expd.code <- c( 'GTE', 'GS13')        # Total expenditure, general government
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

## 3. Create single-country VARs

## 3.1 Create and filter real ave Labor Prod'y
log.real.alp <- data.frame( Year=nom.alp$Year, log( nom.alp[,cts]/ gdp.def[,cts] * 100 ) )
    # Log real ave lab prod
real.alp.hp <- apply( log.real.alp[,cts], 2, hpfilter,
                      freq=6.25, type="lambda", drift=FALSE )
real.alp.trend <- data.frame( Year=nom.alp$Year, sapply( real.alp.hp, function(x) x$trend ) )
real.alp.cycle <- data.frame( Year=nom.alp$Year, sapply( real.alp.hp, function(x) x$cycle ) )
    # Trend and cycle ALP




