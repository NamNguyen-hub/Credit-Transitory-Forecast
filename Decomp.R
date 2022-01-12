# Decomposition Codes
## 1. Load data ----
rm(list=ls())
library(mFilter)
library(xts)
library(tsbox)
library(neverhpfilter)
## BN decomposition filter
# library(devtools)
# devtools::install_github("KevinKotze/tsm")
library(tsm)

## Set working directory
library('rstudioapi')
setwd(dirname(getActiveDocumentContext()$path))
getwd()

# Window of sample data
startdate='1988-12-31'
enddate ='2020-01-01'
startdate_ts = c(1989,1)

country = 'US'

# Importing file
filepath1 = ('../HPCredit/Data Collection/1.Latest/Paper2')
filepath2 = sprintf('/MergedData_%s.txt',country)
filepath = paste(filepath1, filepath2, sep='')

df <- read.table(filepath, header=TRUE, sep=',')
df = subset(df, date > as.Date(startdate)) # Limit series data to after 1990
df = subset(df, date < as.Date(enddate)) # Limit series data to after 1990

#Cycles var name list
varlist = c("ID", "date", "Credit_log", "HPIndex_log")
df = df[varlist]
credit <- xts(df[,-c(1,2,4)], order.by=as.Date(df[,"date"], "%Y-%m-%d"))
hpi <- xts(df[,-c(1,2,3)], order.by=as.Date(df[,"date"], "%Y-%m-%d"))
credit <-ts_ts(credit)
hpi <-ts_ts(hpi)

## 5 Mfilter functions
credit.hp <- mFilter(credit,filter="HP", type = "lambda", freq = 1600)  # Hodrick-Prescott filter
credit.hp3k <- mFilter(credit,filter="HP", type = "lambda", freq = 3000)  # Hodrick-Prescott filter
credit.hp400k <- mFilter(credit,filter="HP", type = "lambda", freq = 400000)  # Hodrick-Prescott filter
# credit.bk <- mFilter(credit,filter="BK")  # Baxter-King filter, I do not use due to tail sampling
# credit.cf <- cffilter(credit, pl=1.5, pu=39.8, type="fixed")  # Christiano-Fitzgerald filter
# credit.cf <- cffilter(credit, pl=6, pu=32, type="asymmetric")
credit.hamilton <- yth_filter(credit, h = 8, p = 4)


credit.bw <- bwfilter(credit, drift=FALSE)  # Butterworth filter
# credit.tr <- mFilter(credit,filter="TR")  # Trigonometric regression filter

credit.linear <- tslm(credit ~ trend) # Linear trend decomp
c.linear <- credit - fitted(credit.linear)
credit.quad <- tslm(credit ~ trend + I(trend^2)) # Quadratic trend decomp
c.quad <- credit - fitted(credit.quad)
bn.decomp <- bnd(credit, nlag = 3) # Beveridge-Nelson decomposition
t.bn <- ts(bn.decomp[, 1], start = startdate_ts, frequency = 4) 
c.bn <- ts(bn.decomp[, 2], start = startdate_ts, frequency = 4) 


#Spectral decomposition
#Wavelet decomp
# library(tsm)
# library(waveslim)
# inf.d4 <- modwt(credit, "d4", n.levels = 3)
# names(inf.d4) <- c("w1", "w2", "w3", "v3")
# inf.d4 <- phase.shift(inf.d4, "d4")
# par(mfrow = c(5, 1))
# plot.ts(credit, axes = FALSE, ylab = "actual", main = "")
# for (i in 1:4) {
#   plot.ts(inf.d4[[i]], axes = FALSE, ylab = names(inf.d4)[i])
# }
# axis(side = 1, at = c(seq(1, (length(credit)), by = 48)), 
#      labels = c("1990", "2000", "2010", "2020"))
# t.wavelet <- ts(inf.d4$v3, start = c(1989, 2), frequency = 4)
# c.wavelet <- credit - t.wavelet
# plot.ts(credit, ylab = "credit")
# lines(t.wavelet, col = "red")
#Dirichlet decomp

credit.cf.cycle = na.omit(credit.cf$cycle)
time(credit.cf.cycle)

## filter cycles to cut off at startdate_uc
## filter enddate 2016-Q4 as CF and 
startdate_uc=c(1989,2)

c.hp<- credit.hp$cycle
c.hp3k <- credit.hp3k$cycle
c.hp400k <-credit.hp400k$cycle
c.cf<- credit.cf$cycle
c.bw<- credit.bw$cycle
c.tr<- credit.tr$cycle

## BK and CF filters both use symmetric sample of 12 periods or 3 years
c.hp = window(c.hp, start=startdate_uc)
c.hp3k = window(c.hp3k, start=startdate_uc)
c.hp400k <- window(c.hp400k, start=startdate_uc)
c.cf = window(c.cf, start=startdate_uc)
c.bw = window(c.bw, start=startdate_uc)
c.tr = window(c.tr, start=startdate_uc)
c.linear = window(c.linear, start=startdate_uc)
c.quad = window(c.quad, start=startdate_uc)
c.bn = window(c.bn, start=startdate_uc)


## Import UC cycle components
df3 <- read.table(sprintf("../HPCredit/Regression/VAR_2_crosscycle/Output/OutputData/uc_yc_%s.txt",country), header=FALSE, sep=",")
## BK and CF filters both use symmetric sample of 12 periods
df3=df3[1:(nrow(df3)),1:2] 
c.uc = ts(df3[,2], start=1989, frequency=4)
t.uc = ts(df3[,1], start=1989, frequency=4)
par(mfrow=c(2,1),mar=c(3,3,2,1),cex=.8)
c.uc = window(c.uc, start=startdate_uc)

plot(credit,main="Credit Series & Estimated Trend", col=1, ylab="")
lines(credit.hp$trend,col=2)
lines(credit.hp3k$trend,col=3)
lines(credit.hp400k$trend,col=4)
lines(credit.bw$trend,col=5)
lines(credit.tr$trend,col=6)
lines(fitted(credit.linear), col=7)
lines(fitted(credit.quad), col=8)
lines(t.bn, col =9)
lines(t.uc, col=10)

legend("topleft",legend=c("series", "HP", "HP3k", "HP400k","BW","TR","linear","quad","BN","UC"),
       col=1:10,lty=rep(1,10),ncol=2)

plot(c.hp,main="Estimated Cyclical Component",
     ylim=c(-20,20),col=2,ylab="")
lines(c.hp3k,col=3)
lines(c.hp400k,col=4)
lines(c.bw,col=5)
lines(c.tr,col=6)
lines(c.linear, col=7)
lines(c.quad, col=8)
lines(c.bn, col=9)
lines(c.uc,col=10)