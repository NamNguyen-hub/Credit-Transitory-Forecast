# Decomposition Codes
## 1. Load library ----
rm(list=ls())
library(mFilter)
library(xts)
library(tsbox) # Convert time series types

library(neverhpfilter)  #Hamilton filter

library(smooth) # MA forecast
library(Mcomp)

## BN decomposition filter
# library(devtools)
# devtools::install_github("KevinKotze/tsm")
library(tsm)

## 
library(forecast)
library(vars)
library(tseries)
# library(dynlm) AIC / BIC model selections packages
# library(zoo)
# library(dyn)

## 2. Load data
## Set working directory
library('rstudioapi')
setwd(dirname(getActiveDocumentContext()$path))

getwd()

# Window of sample data
startdate='1988-12-31'
enddate ='2020-01-01'
startdate_ts = c(1989,1)
startdate_hamilton = c(1986,2) #hamilton filter (need extra 11 periods:  4 lags + 8 periods ahead forecast)
country = 'US'

# Importing file
filepath1 = ('../HPCredit/Data Collection/1.Latest/Paper2')
filepath2 = sprintf('/MergedData_%s.txt',country)
filepath = paste(filepath1, filepath2, sep='')

df <- read.table(filepath, header=TRUE, sep=',')
df = subset(df, date > as.Date(startdate)) # Limit series data to after 1990
df = subset(df, date < as.Date(enddate)) # Limit series data to before 2020

#Cycles var name list
varlist = c("ID", "date", "Credit_log", "HPIndex_log")
df = df[varlist]
credit <- xts(df[,-c(1,2,4)], order.by=as.Date(df[,"date"], "%Y-%m-%d"))
# hpi <- xts(df[,-c(1,2,3)], order.by=as.Date(df[,"date"], "%Y-%m-%d"))
credit0 <-ts_ts(credit)

# Prep data for hamilton filter (need extra 11 periods:  4 lags + 8 periods ahead forecast)
df <- read.table(filepath, header=TRUE, sep=',')
df = subset(df, date > as.Date(startdate_hamilton)) 
varlist = c("ID", "date", "Credit_log", "HPIndex_log")
df = df[varlist]
credit_hamilton <- xts(df[,-c(1,2,4)], order.by=as.Date(df[,"date"], "%Y-%m-%d"))
credit_hamilton1 <- ts_ts(credit_hamilton)

dy = diff(credit0)
dy.true=dy


# BIS credit gap
filepath1 = ('../HPCredit/Data Collection/1.Latest/Paper2')
filepath2 = sprintf('/Credit_BISgap_%s.txt',country)
filepath = paste(filepath1, filepath2, sep='')

df <- read.table(filepath, header=TRUE, sep=',')
df = subset(df, date > as.Date(startdate)) # Limit series data to after 1990
df = subset(df, date < as.Date(enddate)) # Limit series data to before 2020
varlist = c("date", "obs_value")
df = df[varlist]
credit_BIS <- xts(df[,-c(1)], order.by=as.Date(df[,"date"], "%Y-%m-%d"))
credit_BIS <- ts_ts(credit_BIS)


## VAR Models
## dy.true ~ dy.L1 + cy.L1
#### dy = Credit to household first differenced
#### cy = Credit to household cyclical component decomposed using different methods

  n.end=20 #Initial sample estimation
  t=length(dy) #Full Sample size
  n=t-n.end-3 #Forecast sample
  
  
invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))
  
library(dlm)

filterHP <- function(series,lambda=1600){
  # Function description: HP filter using DLM package.
  if(!"ts" %in% class(series)) stop("series is not a \"ts\" object.")
  
  # Set state priors
  level <- series[1]
  slope <- mean(diff(series),na.rm=TRUE)
  
  # Set up HP filter in a DLM model
  model <- function(param){
    trend <- dlmModPoly(dV = 1,
                        dW = c(0,1/lambda),
                        m0 = c(level,slope),
                        C0 = 2 * diag(2))
    
    # AR(2) model
    cycle <- dlmModARMA(ar     = ARtransPars(c(0,0)),
                        sigma2 = 1e-07)      
    
    hp_filter_dlm <- trend  + cycle
    
    return(hp_filter_dlm)
  }
  
  # MLE Estimation
  MLE       <- dlmMLE(series,c(0.5,0.4),model)
  # Estimated parameters
  EstParams <- model(MLE$par)
  # # Smoothed series
  # Smooth_Estimates <- dlmSmooth(series,EstParams)
  # Filtered series
  Smooth_Estimates = dlmFilter(series,EstParams)
  
  # Trend and Cycle
  trend <- Smooth_Estimates$m[,1]
  cycle <- series - trend
  
  # Plot the data ---
  par(mfrow = c(2,1),
      oma = c(1,3,0,0) + 0.1,
      mar = c(1.5,1,1,1) + 0.1)
  plot(series,las=1,col="black")
  lines(trend,col="blue")
  legend("topleft",legend=c("Observed","Trend"),border=FALSE,bty="n",col=c("black","blue"),lwd=1)
  title(main="HP Filter - Trend")
  par(mar = c(1,1,1.5,1) + 0.1)
  plot(cycle,las=1,col="red")
  title(main="HP Filter - Cycle")
  abline(h=0)
  par(mfrow = c(1,1),
      oma = c(0,0,0,0),
      mar = c(5.1,4.1,4.1,2.1))  
  
  # Return the data
  data <- ts.union(series,trend,cycle)
  return(data)
}

# Load some data
data(USecon)

# Compare two implementations
filterHP(USecon[,"M1"],lambda=1600)

filterHP(credit0, lambda=1600)

library(mFilter)
par(mfrow=c(1,1))

out1=filterHP(credit0, lambda=1600)[,"cycle"]

plot(hpfilter(credit0,freq=1600)$cycle, ylim=c(-5,5), col=1,ylab="")
lines(out1, col=2)
legend("topleft",legend=c("HPfilter 2 sided", "HPfilter 1 sided"),
       col=1:2,lty=rep(1,2),ncol=2)

# The data used is annual, so arguably lambda = 1600 should really be 6.25.