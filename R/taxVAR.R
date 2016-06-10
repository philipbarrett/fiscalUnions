library(readr)
# library(tidyr)
# library(dplyr)
library(vars)
# library(lubridate)

## 0. Script variables
file <- '~/Dropbox/data/2016/fiscalUnions/taxRevenue.csv'
cts <- c("Germany", "France")
gov <- 'FED'

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

plot( tax$Year, tax[,cts[1]], type='l', ylim = c(10,20), col='blue', lwd=2 )
lines( tax$Year, tax[,cts[2]], type='l', col='red', lwd=2 )

VAR( tax[,cts] )
