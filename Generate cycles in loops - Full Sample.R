# Decomposition Codes
## 1. Load library ----
rm(list=ls())

country = 'JP'

library(mFilter)
library(xts)
library(tsbox) # Convert time series types

library(neverhpfilter)  #Hamilton filter

library(smooth) # MA forecast
library(Mcomp)
library(dplyr)
## BN decomposition filter
# library(devtools)
# devtools::install_github("KevinKotze/tsm")
library(tsm)

## 
library(forecast)
library(vars)
library(tseries)


## 2. Load data
## Set working directory
library('rstudioapi')
setwd(dirname(getActiveDocumentContext()$path))
getwd()

## One sided HP filter function
source("HPfilters/OneSidedHPfilterfunc.R")


# Window of sample data
startdate_var_input='1988-10-01'
startdate_var='1989-01-01'
enddate ='2021-04-01'
startdate_ts = c(1988,4)
startdate_diff = c(1989,1) #sample date is pushed forward because of diff operator
startdate_hamilton_xts = '1983-01-01'
startdate_hamilton = c(1983,1) #hamilton filter (need extra 23 periods:  4 lags + 20 periods ahead forecast)

# Importing file
filepath1 = ('../HPCredit/Data Collection/1.Latest/Paper2')
filepath2 = sprintf('/credit_%s.txt',country)
filepath = paste(filepath1, filepath2, sep='')

df <- read.table(filepath, header=TRUE, sep=',')
# Prep data for hamilton filter (need extra 23 periods:  4 lags + 20 periods ahead forecast)
df1 <- df %>%
  filter(grepl("^(A)", cg_dtype))
df1 = subset(df1, date >= as.Date(startdate_hamilton_xts)) 
df1 = subset(df1, date <= as.Date(enddate)) # Limit series data to before 2020
varlist = c("date", "obs_value")
df1 = df1[varlist]
credit_hamilton <- xts(df1[,-c(1)], order.by=as.Date(df1[,"date"], "%Y-%m-%d"))

#Log transforming
# credit_hamilton = 100*log(credit_hamilton)
credit_hamilton1 <- ts_ts(credit_hamilton)


# Prep data for dy (diff(credit))
credit0=window(credit_hamilton1, start=startdate_var_input)
dy = diff(credit_hamilton1)
dy = window(dy, start=startdate_diff)
dy.true=dy

## HP filter series

# df1 <- df %>%
#   filter(grepl("^(A)", cg_dtype))
# varlist = c("date", "obs_value")
# df1 = df1[varlist]
# df1 = subset(df1, date >= as.Date(startdate_hamilton_xts))
# df1 = subset(df1, date <= as.Date(enddate))
# credit <- xts(df1[,-c(1)], order.by=as.Date(df1[,"date"], "%Y-%m-%d"))
credit1 <- credit_hamilton1
c.hp=filterHP(credit1)[,"cycle"]
c.hp3k=filterHP(credit1, lambda=3000)[,"cycle"]
c.hp400k=filterHP(credit1, lambda=400000)[,"cycle"]

x<-list(c.hp,c.hp3k,c.hp400k)
x<-lapply(x, window, start=startdate_diff)
c.hp1=x[1][[1]]
c.hp3k1=x[2][[1]]
c.hp400k1=x[3][[1]]
## VAR Models
## dy.true ~ dy.L1 + cy.L1
#### dy = Credit to household first differenced
#### cy = Credit to household cyclical component decomposed using different methods

t=length(dy) #Full Sample size

# to account for first difference of dependent variable diff(credit)
## UC components
## Import UC cycle components
df3 <- read.table(sprintf("../HPCredit/Regression/VAR_2/Output/OutputData/uc_yc_%s.txt",country), header=FALSE, sep=",")
## BK and CF filters both use symmetric sample of 12 periods
df3=df3[,1:2] 
c.uc = ts(df3[,2], start=startdate_hamilton, frequency=4)
c.uc1 = window(c.uc, start=startdate_diff)
i=1
#   

c.hamilton2=matrix(0,t,1) #store of 1 sided cycle decomp
c.linear2=matrix(0,t,1)
c.quad2=matrix(0,t,1)
c.bn2=matrix(0,t,1)


for(i in 1:t){
  
  
  credit_hamilton = credit_hamilton1[(1):(i-1+24)]
  credit_hamilton = ts(credit_hamilton, start = startdate_hamilton, frequency =4)
  credit_xts = ts_xts(credit_hamilton)
  
  #### Cycle creation
  # credit.hp <- mFilter(credit,filter="HP", type = "lambda", freq = 1600)  # Hodrick-Prescott filter
  # credit.hp3k <- mFilter(credit,filter="HP", type = "lambda", freq = 3000)  # Hodrick-Prescott filter
  # credit.hp400k <- mFilter(credit,filter="HP", type = "lambda", freq = 400000)  # Hodrick-Prescott filter
  
  
  c.hamilton4 <- yth_filter(credit_xts, h = 20, p = 4)$value.cycle # Hamilton filter
  c.hamilton2[(1+i-1)]=c.hamilton4[i-1+24]
  
  # credit.bw <- bwfilter(credit1, drift=FALSE)  # Butterworth filter
  credit.linear <- tslm(credit_hamilton ~ trend) # Linear trend decomp
  c.linear <- credit_hamilton - fitted(credit.linear)
  c.linear2[1+i-1]=c.linear[i-1+24]
  
  credit.quad <- tslm(credit_hamilton ~ trend + I(trend^2)) # Quadratic trend decomp
  c.quad <- credit_hamilton - fitted(credit.quad)
  c.quad2[1+i-1]=c.quad[i-1+24]
  
  bn.decomp <- bnd(credit_hamilton, nlag = 2) # Beveridge-Nelson decomposition
  c.bn <- ts(bn.decomp[, 2], start = startdate_hamilton, frequency = 4) 
  c.bn2[1+i-1]=c.bn[i-1+24]
  
} 

x=list(c.hp1,c.hp3k1,c.hp400k1,c.uc1,c.hamilton2, c.linear2 , c.quad2, c.bn2)
# x<-lapply(x, function(x) x=x[1:(n.end+i-1)])
x<-lapply(x, ts, start=startdate_diff, frequency=4)

c.hp=x[1][[1]]
c.hp3k=x[2][[1]]
c.hp400k=x[3][[1]]
c.uc=x[4][[1]]
c.hamilton = x[5][[1]]
c.linear = x[6][[1]]
c.quad = x[7][[1]]
c.bn = x[8][[1]]


dy=window(dy, start=startdate_diff)
dy.true=dy

t=length(dy)
## Save data into csv file

## Load dy
startdate_ts = startdate_diff # assign differenced values to 1 period after sample startdate_var, 

var_1=ts(cbind(dy,c.hp), start=startdate_ts, frequency =4)
var_2=ts(cbind(dy,c.hp3k), start=startdate_ts, frequency =4)
var_3=ts(cbind(dy,c.hp400k), start=startdate_ts, frequency =4)
var_4=ts(cbind(dy,c.hamilton), start=startdate_ts, frequency =4)
# var_5=ts(cbind(dy,c.bw), start=startdate_ts, frequency =4)
var_5=ts(cbind(dy,c.linear), start=startdate_ts, frequency =4)
var_6=ts(cbind(dy,c.quad), start=startdate_ts, frequency =4)
var_7=ts(cbind(dy,c.bn), start=startdate_ts, frequency =4)
var_8=ts(cbind(dy,c.uc), start=startdate_ts, frequency =4)


  i=1

n.end=20
n=t-n.end-3


# set matrix for storage
pred_var1=matrix(0,n,4)
pred_var2=matrix(0,n,4)
pred_var3=matrix(0,n,4)
pred_var4=matrix(0,n,4)
pred_var5=matrix(0,n,4)
pred_var6=matrix(0,n,4)
pred_var7=matrix(0,n,4)
pred_var8=matrix(0,n,4)


pred_ar=matrix(0,n,4)
pred_comb=matrix(0,n,4)

w_bg1=matrix(0,n,8)#need separate weights for each horizon in Bates-Granger method
w_bg2=matrix(0,n,8)
w_bg3=matrix(0,n,8)
w_bg4=matrix(0,n,8)
DMW=matrix(0,n,1)
c.combined=matrix(NA,t,1)
c.df = matrix(NA,t,9)
c.weight = matrix(NA,t,8)
  
for(i in 1:n){
  

#### Run regression
  x_var1 = var_1[1:(n.end+i-1),]
  x_var2 = var_2[1:(n.end+i-1),]
  x_var3 = var_3[1:(n.end+i-1),]
  x_var4 = var_4[1:(n.end+i-1),]
  x_var5 = var_5[1:(n.end+i-1),]
  x_var6 = var_6[1:(n.end+i-1),]
  x_var7 = var_7[1:(n.end+i-1),]
  x_var8 = var_8[1:(n.end+i-1),]

  x_ar=dy[1:(n.end+i-1)]
    ## set up for AR(1) regression
  # 
  # info.crit1=VARselect(x_var21,lag.max=4,type="const")
  # info.crit2=VARselect(var_22,lag.max=4,type="const")
  # info.crit3=VARselect(var_23,lag.max=4,type="const")
  # info.crit4=VARselect(var_24,lag.max=4,type="const")
  # 
  # info.crit31=VARselect(var_31,lag.max=4,type="const")
  # info.crit32=VARselect(var_32,lag.max=4,type="const") 
  # info.crit33=VARselect(var_33,lag.max=4,type="const")
  # 
  # n_1=info.crit1$selection[3]
  # n_2=info.crit2$selection[3]  
  # n_3=info.crit3$selection[3]  
  # n_4=info.crit4$selection[3]    
  # n_31=info.crit31$selection[3]
  # n_32=info.crit32$selection[3]
  # n_33=info.crit33$selection[3]
  
  model.var1=VAR(x_var1,type="const",ic="SC")
  for_var1=predict(model.var1,n.ahead=4,se.fit=FALSE)
  pred_var1[i,1:4]=for_var1$fcst$dy[1:4] 
  
  model.var2=VAR(x_var2,type="const",ic="SC")
  for_var2=predict(model.var2,n.ahead=4,se.fit=FALSE)
  pred_var2[i,1:4]=for_var2$fcst$dy[1:4] 
  
  model.var3=VAR(x_var3,type="const",ic="SC")
  for_var3=predict(model.var3,n.ahead=4,se.fit=FALSE)
  pred_var3[i,1:4]=for_var3$fcst$dy[1:4] 
  
  model.var4=VAR(x_var4,type="const",ic="SC")
  for_var4=predict(model.var4,n.ahead=4,se.fit=FALSE)
  pred_var4[i,1:4]=for_var4$fcst$dy[1:4] 
  
  model.var5=VAR(x_var5,type="const",ic="SC")
  for_var5=predict(model.var5,n.ahead=4,se.fit=FALSE)
  pred_var5[i,1:4]=for_var5$fcst$dy[1:4]
  
  model.var6=VAR(x_var6,type="const",ic="SC")
  for_var6=predict(model.var6,n.ahead=4,se.fit=FALSE)
  pred_var6[i,1:4]=for_var6$fcst$dy[1:4] 
  
  model.var7=VAR(x_var7,type="const",ic="SC")
  for_var7=predict(model.var7,n.ahead=4,se.fit=FALSE)
  pred_var7[i,1:4]=for_var7$fcst$dy[1:4] 
  
  model.var8=VAR(x_var8,type="const",ic="SC")
  for_var8=predict(model.var8,n.ahead=4,se.fit=FALSE)
  pred_var8[i,1:4]=for_var8$fcst$dy[1:4] 
  
  # model.var9=VAR(x_var9,type="const",ic="SC")
  # for_var9=predict(model.var9,n.ahead=4,se.fit=FALSE)
  # pred_var9[i,1:4]=for_var9$fcst$dy[1:4] 
  
  model.ar=arima(x_ar,order=c(1,0,0),method="ML")
  pred_ar[i,1:4]=predict(model.ar,n.ahead=4,se.fit=FALSE)[1:4]
  
  ###Bates-Granger
  # sam_size=n.end+i-1
  
  mse1_1 = mean((dy.true[(n.end+1):(t-n+i-3)] - pred_var1[1:i,1])^2)
  mse2_1 = mean((dy.true[(n.end+2):(t-n+i-2)] - pred_var1[1:i,2])^2)
  mse3_1 = mean((dy.true[(n.end+3):(t-n+i-1)] - pred_var1[1:i,3])^2)
  mse4_1 = mean((dy.true[(n.end+4):(t-n+i)] - pred_var1[1:i,4])^2)
  
  mse1_2 = mean((dy.true[(n.end+1):(t-n+i-3)] - pred_var2[1:i,1])^2)
  mse2_2 = mean((dy.true[(n.end+2):(t-n+i-2)] - pred_var2[1:i,2])^2)
  mse3_2 = mean((dy.true[(n.end+3):(t-n+i-1)] - pred_var2[1:i,3])^2)
  mse4_2 = mean((dy.true[(n.end+4):(t-n+i)] - pred_var2[1:i,4])^2)
  
  mse1_3 = mean((dy.true[(n.end+1):(t-n+i-3)] - pred_var3[1:i,1])^2)
  mse2_3 = mean((dy.true[(n.end+2):(t-n+i-2)] - pred_var3[1:i,2])^2)
  mse3_3 = mean((dy.true[(n.end+3):(t-n+i-1)] - pred_var3[1:i,3])^2)
  mse4_3 = mean((dy.true[(n.end+4):(t-n+i)] - pred_var3[1:i,4])^2)
  
  mse1_4 = mean((dy.true[(n.end+1):(t-n+i-3)] - pred_var4[1:i,1])^2)
  mse2_4 = mean((dy.true[(n.end+2):(t-n+i-2)] - pred_var4[1:i,2])^2)
  mse3_4 = mean((dy.true[(n.end+3):(t-n+i-1)] - pred_var4[1:i,3])^2)
  mse4_4 = mean((dy.true[(n.end+4):(t-n+i)] - pred_var4[1:i,4])^2)
  
  mse1_5 = mean((dy.true[(n.end+1):(t-n+i-3)] - pred_var5[1:i,1])^2)
  mse2_5 = mean((dy.true[(n.end+2):(t-n+i-2)] - pred_var5[1:i,2])^2)
  mse3_5 = mean((dy.true[(n.end+3):(t-n+i-1)] - pred_var5[1:i,3])^2)
  mse4_5 = mean((dy.true[(n.end+4):(t-n+i)] - pred_var5[1:i,4])^2)
  
  mse1_6 = mean((dy.true[(n.end+1):(t-n+i-3)] - pred_var6[1:i,1])^2)
  mse2_6 = mean((dy.true[(n.end+2):(t-n+i-2)] - pred_var6[1:i,2])^2)
  mse3_6 = mean((dy.true[(n.end+3):(t-n+i-1)] - pred_var6[1:i,3])^2)
  mse4_6 = mean((dy.true[(n.end+4):(t-n+i)] - pred_var6[1:i,4])^2)
  
  mse1_7 = mean((dy.true[(n.end+1):(t-n+i-3)] - pred_var7[1:i,1])^2)
  mse2_7 = mean((dy.true[(n.end+2):(t-n+i-2)] - pred_var7[1:i,2])^2)
  mse3_7 = mean((dy.true[(n.end+3):(t-n+i-1)] - pred_var7[1:i,3])^2)
  mse4_7 = mean((dy.true[(n.end+4):(t-n+i)] - pred_var7[1:i,4])^2)
  
  mse1_8 = mean((dy.true[(n.end+1):(t-n+i-3)] - pred_var8[1:i,1])^2)
  mse2_8 = mean((dy.true[(n.end+2):(t-n+i-2)] - pred_var8[1:i,2])^2)
  mse3_8 = mean((dy.true[(n.end+3):(t-n+i-1)] - pred_var8[1:i,3])^2)
  mse4_8 = mean((dy.true[(n.end+4):(t-n+i)] - pred_var8[1:i,4])^2)
  
  # mse1_9 = mean((dy.true[(n.end+1):(t-n+i-3)] - pred_var9[1:i,1])^2)
  # mse2_9 = mean((dy.true[(n.end+2):(t-n+i-2)] - pred_var9[1:i,2])^2)
  # mse3_9 = mean((dy.true[(n.end+3):(t-n+i-1)] - pred_var9[1:i,3])^2)
  # mse4_9 = mean((dy.true[(n.end+4):(t-n+i)] - pred_var9[1:i,4])^2)
  

  
  wbgstr_11=1/mse1_1
  wbgstr_12=1/mse1_2
  wbgstr_13=1/mse1_3
  wbgstr_14=1/mse1_4
  wbgstr_15=1/mse1_5
  wbgstr_16=1/mse1_6
  wbgstr_17=1/mse1_7
  wbgstr_18=1/mse1_8
  # wbgstr_19=1/mse1_9


  
  wbgstr_21=1/mse2_1
  wbgstr_22=1/mse2_2
  wbgstr_23=1/mse2_3
  wbgstr_24=1/mse2_4
  wbgstr_25=1/mse2_5
  wbgstr_26=1/mse2_6
  wbgstr_27=1/mse2_7
  wbgstr_28=1/mse2_8
  # wbgstr_29=1/mse2_9

  
  
  wbgstr_31=1/mse3_1
  wbgstr_32=1/mse3_2
  wbgstr_33=1/mse3_3
  wbgstr_34=1/mse3_4
  wbgstr_35=1/mse3_5
  wbgstr_36=1/mse3_6
  wbgstr_37=1/mse3_7
  wbgstr_38=1/mse3_8
  # wbgstr_39=1/mse3_9


  
  wbgstr_41=1/mse4_1
  wbgstr_42=1/mse4_2
  wbgstr_43=1/mse4_3
  wbgstr_44=1/mse4_4
  wbgstr_45=1/mse4_5
  wbgstr_46=1/mse4_6
  wbgstr_47=1/mse4_7
  wbgstr_48=1/mse4_8
  # wbgstr_49=1/mse4_9


  
  wbgstr1_sum=(wbgstr_11+wbgstr_12+wbgstr_13+wbgstr_14+wbgstr_15+wbgstr_16+wbgstr_17+wbgstr_18)
  wbgstr2_sum=(wbgstr_21+wbgstr_22+wbgstr_23+wbgstr_24+wbgstr_25+wbgstr_26+wbgstr_27+wbgstr_28)
  wbgstr3_sum=(wbgstr_31+wbgstr_32+wbgstr_33+wbgstr_34+wbgstr_35+wbgstr_36+wbgstr_37+wbgstr_38)
  wbgstr4_sum=(wbgstr_41+wbgstr_42+wbgstr_43+wbgstr_44+wbgstr_45+wbgstr_46+wbgstr_47+wbgstr_48)
  
  
  wbg_11=wbgstr_11/wbgstr1_sum
  wbg_12=wbgstr_12/wbgstr1_sum
  wbg_13=wbgstr_13/wbgstr1_sum
  wbg_14=wbgstr_14/wbgstr1_sum  
  wbg_15=wbgstr_15/wbgstr1_sum
  wbg_16=wbgstr_16/wbgstr1_sum  
  wbg_17=wbgstr_17/wbgstr1_sum  
  wbg_18=wbgstr_18/wbgstr1_sum  
  # wbg_19=wbgstr_19/wbgstr1_sum  
  
  
  wbg_21=wbgstr_21/wbgstr2_sum
  wbg_22=wbgstr_22/wbgstr2_sum
  wbg_23=wbgstr_23/wbgstr2_sum
  wbg_24=wbgstr_24/wbgstr2_sum 
  wbg_25=wbgstr_25/wbgstr2_sum
  wbg_26=wbgstr_26/wbgstr2_sum  
  wbg_27=wbgstr_27/wbgstr2_sum  
  wbg_28=wbgstr_28/wbgstr2_sum  
  # wbg_29=wbgstr_29/wbgstr2_sum  
  
  wbg_31=wbgstr_31/wbgstr3_sum
  wbg_32=wbgstr_32/wbgstr3_sum
  wbg_33=wbgstr_33/wbgstr3_sum
  wbg_34=wbgstr_34/wbgstr3_sum 
  wbg_35=wbgstr_35/wbgstr3_sum
  wbg_36=wbgstr_36/wbgstr3_sum  
  wbg_37=wbgstr_37/wbgstr3_sum  
  wbg_38=wbgstr_38/wbgstr3_sum  
  # wbg_39=wbgstr_39/wbgstr3_sum  
  
  wbg_41=wbgstr_41/wbgstr4_sum
  wbg_42=wbgstr_42/wbgstr4_sum
  wbg_43=wbgstr_43/wbgstr4_sum
  wbg_44=wbgstr_44/wbgstr4_sum 
  wbg_45=wbgstr_45/wbgstr4_sum
  wbg_46=wbgstr_46/wbgstr4_sum  
  wbg_47=wbgstr_47/wbgstr4_sum  
  wbg_48=wbgstr_48/wbgstr4_sum  
  # wbg_49=wbgstr_49/wbgstr4_sum  
  
  w_bg1[i,1:8]=c(wbg_11,wbg_12,wbg_13,wbg_14,wbg_15,wbg_16,wbg_17,wbg_18)
  w_bg2[i,1:8]=c(wbg_21,wbg_22,wbg_23,wbg_24,wbg_25,wbg_26,wbg_27,wbg_28)
  w_bg3[i,1:8]=c(wbg_31,wbg_32,wbg_33,wbg_34,wbg_35,wbg_36,wbg_37,wbg_38)
  w_bg4[i,1:8]=c(wbg_41,wbg_42,wbg_43,wbg_44,wbg_45,wbg_46,wbg_47,wbg_48) 

  ## Calculate combined cycle

  c.combined[i+n.end] = (c.hp[i+n.end]*w_bg1[i,1]+c.hp3k[i+n.end]*w_bg1[i,2]+c.hp400k[i+n.end]*w_bg1[i,3]
                         +c.hamilton[i+n.end]*w_bg1[i,4]+c.linear[i+n.end]*w_bg1[i,5]+c.quad[i+n.end]*w_bg1[i,6]
                   +c.bn[i+n.end]*w_bg1[i,7]+c.uc[i+n.end]*w_bg1[i,8])
 
  

  
  
  
   ## Storage of cyclical components
  c.df[i+n.end-1,] = c(c.hp[i+n.end],c.hp3k[i+n.end],c.hp400k[i+n.end],
                       c.hamilton[i+n.end],c.linear[i+n.end],c.quad[i+n.end],
                       c.bn[i+n.end],c.uc[i+n.end], c.combined[i+n.end])
 
  
   c.weight[i+n.end-1,] = c(mean(w_bg1[i,1],w_bg2[i,1],w_bg3[i,1],w_bg4[i,1]),
                            mean(w_bg1[i,2],w_bg2[i,2],w_bg3[i,2],w_bg4[i,2]),
                            mean(w_bg1[i,3],w_bg2[i,3],w_bg3[i,3],w_bg4[i,3]),
                            mean(w_bg1[i,4],w_bg2[i,4],w_bg3[i,4],w_bg4[i,4]),
                            mean(w_bg1[i,5],w_bg2[i,5],w_bg3[i,5],w_bg4[i,5]),
                            mean(w_bg1[i,6],w_bg2[i,6],w_bg3[i,6],w_bg4[i,6]),
                            mean(w_bg1[i,7],w_bg2[i,7],w_bg3[i,7],w_bg4[i,7]),
                            mean(w_bg1[i,8],w_bg2[i,8],w_bg3[i,8],w_bg4[i,8]))
  
  
   c.combined[i+n.end] =  (c.hp[i+n.end]*c.weight[i+n.end-1,1]
                           +c.hp3k[i+n.end]*c.weight[i+n.end-1,2]
                           +c.hp400k[i+n.end]*c.weight[i+n.end-1,3]
                           +c.hamilton[i+n.end]*c.weight[i+n.end-1,4]
                           +c.linear[i+n.end]*c.weight[i+n.end-1,5]
                           +c.quad[i+n.end]*c.weight[i+n.end-1,6]
                           +c.bn[i+n.end]*c.weight[i+n.end-1,7]
                           +c.uc[i+n.end]*c.weight[i+n.end-1,8])
   
   
  # c.hp1[i+n.end-1] = c.hp[i+n.end-1]
  # c.hp3k1[i+n.end-1] = c.hp3k[i+n.end-1]
  # c.hp400k1[i+n.end-1] = c.hp400k[i+n.end-1]
  # c.hamilton1[i+n.end-1] = c.hamilton[i+n.end-1]
  # c.bw1[i+n.end-1] = c.bw[i+n.end-1]
  # c.linear1[i+n.end-1] = c.linear[i+n.end-1]
  # c.quad1[i+n.end-1] = c.quad[i+n.end-1]
  # c.bn1[i+n.end-1] = c.bn[i+n.end-1]
  # c.uc1[i+n.end-1] = c.uc[i+n.end-1]
  
  print(c.combined)
  print(i)
# ## Calcalte DMW time series values
# ### Squared error of Model A: AR(1) forecast  
# e_sqr_ar1=(dy.true[(i+1)]-pred_ar[i,1])^2#1-step ahead forecast error
# ### Squared error of Model B: combination avarage forecast
# pred_c1=(pred_var1[i,1]+pred_var2[i,1]+pred_var2[i,1]+pred_var4[i,1]+pred_var5[i,1]+pred_var6[i,1])/6
# e_sqr_c1=(dy.true[(i+1)]-pred_c1)^2
# ### DMW time series:
# DMW[i]=e_sqr_c1-e_sqr_ar1
  
}
### Extract MSE from lm function ----
  
## Create combined cyclical components
### Save all MSE and calculate normalized weights
### Create combined cyclical component
#### Read notes and problem set from 735: combined forecast
pred_bg1=(w_bg1[,1]*pred_var1[,1]+w_bg1[,2]*pred_var2[,1]
          +w_bg1[,3]*pred_var3[,1]+w_bg1[,4]*pred_var4[,1]
          +w_bg1[,5]*pred_var5[,1]+w_bg1[,6]*pred_var6[,1]
          +w_bg1[,7]*pred_var7[,1]+w_bg1[,8]*pred_var8[,1])

pred_bg2=(w_bg2[,1]*pred_var1[,2]+w_bg2[,2]*pred_var2[,2]
          +w_bg2[,3]*pred_var3[,2]+w_bg2[,4]*pred_var4[,2]
          +w_bg2[,5]*pred_var5[,2]+w_bg2[,6]*pred_var6[,2]
          +w_bg2[,7]*pred_var7[,2]+w_bg2[,8]*pred_var8[,2])

pred_bg3=(w_bg3[,1]*pred_var1[,3]+w_bg3[,2]*pred_var2[,3]
          +w_bg3[,3]*pred_var3[,3]+w_bg3[,4]*pred_var4[,3]
          +w_bg3[,5]*pred_var5[,3]+w_bg3[,6]*pred_var6[,3]
          +w_bg3[,7]*pred_var7[,3]+w_bg3[,8]*pred_var8[,3])

pred_bg4=(w_bg4[,1]*pred_var1[,4]+w_bg4[,2]*pred_var2[,4]
          +w_bg4[,3]*pred_var3[,4]+w_bg4[,4]*pred_var4[,4]
          +w_bg4[,5]*pred_var5[,4]+w_bg4[,6]*pred_var6[,4]
          +w_bg4[,7]*pred_var7[,4]+w_bg4[,8]*pred_var8[,4])
pred_bg=cbind(pred_bg1,pred_bg2,pred_bg3,pred_bg4)

pred_c=(pred_var1+pred_var2+pred_var3+pred_var4+pred_var5+pred_var6+pred_var7+pred_var8)/8


## Test forecast retrospect and prospect
### Final
### Quasi-real time (QTR)




# prediction errors
# compute prediction errors

e1_varcomb=dy.true[(n.end+1):(t-3)]-pred_c[1:n,1]#1-step ahead forecast error
e2_varcomb=dy.true[(n.end+2):(t-2)]-pred_c[1:n,2] #2-step ahead forecast error
e3_varcomb=dy.true[(n.end+3):(t-1)]-pred_c[1:n,3]#1-step ahead forecast error
e4_varcomb=dy.true[(n.end+4):t]-pred_c[1:n,4] #2-step ahead forecast error
e14_varcomb=(e1_varcomb+e2_varcomb+e3_varcomb+e4_varcomb)/4

rmse1_varcomb=sqrt(mean(e1_varcomb^2))
rmse2_varcomb=sqrt(mean(e2_varcomb^2))
rmse3_varcomb=sqrt(mean(e3_varcomb^2))
rmse4_varcomb=sqrt(mean(e4_varcomb^2))
rmse14_varcomb=sqrt(mean(e14_varcomb^2))

e1_ar=dy.true[(n.end+1):(t-3)]-pred_ar[1:n,1]#1-step ahead forecast error
e2_ar=dy.true[(n.end+2):(t-2)]-pred_ar[1:n,2] #2-step ahead forecast error
e3_ar=dy.true[(n.end+3):(t-1)]-pred_ar[1:n,3]#1-step ahead forecast error
e4_ar=dy.true[(n.end+4):t]-pred_ar[1:n,4] #2-step ahead forecast error
e14_ar=(e1_ar+e2_ar+e3_ar+e4_ar)/4

rmse1_ar=sqrt(mean(e1_ar^2))
rmse2_ar=sqrt(mean(e2_ar^2))
rmse3_ar=sqrt(mean(e3_ar^2))
rmse4_ar=sqrt(mean(e4_ar^2))
rmse14_ar=sqrt(mean(e14_ar^2))

e1_varbg=dy.true[(n.end+1):(t-3)]-pred_bg[1:n,1]#1-step ahead forecast error
e2_varbg=dy.true[(n.end+2):(t-2)]-pred_bg[1:n,2] #2-step ahead forecast error
e3_varbg=dy.true[(n.end+3):(t-1)]-pred_bg[1:n,3]#3-step ahead forecast error
e4_varbg=dy.true[(n.end+4):t]-pred_bg[1:n,4] #4-step ahead forecast error
e14_varbg=(e1_varbg+e2_varbg+e3_varbg+e4_varbg)/4

rmse1_varbg=sqrt(mean(e1_varbg^2))
rmse2_varbg=sqrt(mean(e2_varbg^2))
rmse3_varbg=sqrt(mean(e3_varbg^2))
rmse4_varbg=sqrt(mean(e4_varbg^2))
rmse14_varbg=sqrt(mean(e14_varbg^2))


## HP filter rmse
e1_var1=dy.true[(n.end+1):(t-3)]-pred_var1[1:n,1]#1-step ahead forecast error
e2_var1=dy.true[(n.end+2):(t-2)]-pred_var1[1:n,2] #2-step ahead forecast error
e3_var1=dy.true[(n.end+3):(t-1)]-pred_var1[1:n,3]#3-step ahead forecast error
e4_var1=dy.true[(n.end+4):t]-pred_var1[1:n,4] #4-step ahead forecast error
e14_var1=(e1_var1+e2_var1+e3_var1+e4_var1)/4

rmse1_var1=sqrt(mean(e1_var1^2))
rmse2_var1=sqrt(mean(e2_var1^2))
rmse3_var1=sqrt(mean(e3_var1^2))
rmse4_var1=sqrt(mean(e4_var1^2))
rmse14_var1=sqrt(mean(e14_var1^2))

##HP filter lambda = 3k rmse
e1_var2=dy.true[(n.end+1):(t-3)]-pred_var2[1:n,1]#1-step ahead forecast error
e2_var2=dy.true[(n.end+2):(t-2)]-pred_var2[1:n,2] #2-step ahead forecast error
e3_var2=dy.true[(n.end+3):(t-1)]-pred_var2[1:n,3]#3-step ahead forecast error
e4_var2=dy.true[(n.end+4):t]-pred_var2[1:n,4] #4-step ahead forecast error
e14_var2=(e1_var2+e2_var2+e3_var2+e4_var2)/4

rmse1_var2=sqrt(mean(e1_var2^2))
rmse2_var2=sqrt(mean(e2_var2^2))
rmse3_var2=sqrt(mean(e3_var2^2))
rmse4_var2=sqrt(mean(e4_var2^2))
rmse14_var2=sqrt(mean(e14_var2^2))

##HP filter lambda = 400k rmse
e1_var3=dy.true[(n.end+1):(t-3)]-pred_var3[1:n,1]#1-step ahead forecast error
e2_var3=dy.true[(n.end+2):(t-2)]-pred_var3[1:n,2] #2-step ahead forecast error
e3_var3=dy.true[(n.end+3):(t-1)]-pred_var3[1:n,3]#3-step ahead forecast error
e4_var3=dy.true[(n.end+4):t]-pred_var3[1:n,4] #4-step ahead forecast error
e14_var3=(e1_var3+e2_var3+e3_var3+e4_var3)/4

rmse1_var3=sqrt(mean(e1_var3^2))
rmse2_var3=sqrt(mean(e2_var3^2))
rmse3_var3=sqrt(mean(e3_var3^2))
rmse4_var3=sqrt(mean(e4_var3^2))
rmse14_var3=sqrt(mean(e14_var3^2))


##Hamilton filter
e1_var4=dy.true[(n.end+1):(t-3)]-pred_var4[1:n,1]#1-step ahead forecast error
e2_var4=dy.true[(n.end+2):(t-2)]-pred_var4[1:n,2] #2-step ahead forecast error
e3_var4=dy.true[(n.end+3):(t-1)]-pred_var4[1:n,3]#3-step ahead forecast error
e4_var4=dy.true[(n.end+4):t]-pred_var4[1:n,4] #4-step ahead forecast error
e14_var4=(e1_var4+e2_var4+e3_var4+e4_var4)/4

rmse1_var4=sqrt(mean(e1_var4^2))
rmse2_var4=sqrt(mean(e2_var4^2))
rmse3_var4=sqrt(mean(e3_var4^2))
rmse4_var4=sqrt(mean(e4_var4^2))
rmse14_var4=sqrt(mean(e14_var4^2))

##linear filter
e1_var5=dy.true[(n.end+1):(t-3)]-pred_var5[1:n,1]#1-step ahead forecast error
e2_var5=dy.true[(n.end+2):(t-2)]-pred_var5[1:n,2] #2-step ahead forecast error
e3_var5=dy.true[(n.end+3):(t-1)]-pred_var5[1:n,3]#3-step ahead forecast error
e4_var5=dy.true[(n.end+4):t]-pred_var5[1:n,4] #4-step ahead forecast error
e14_var5=(e1_var5+e2_var5+e3_var5+e4_var5)/4

rmse1_var5=sqrt(mean(e1_var5^2))
rmse2_var5=sqrt(mean(e2_var5^2))
rmse3_var5=sqrt(mean(e3_var5^2))
rmse4_var5=sqrt(mean(e4_var5^2))
rmse14_var5=sqrt(mean(e14_var5^2))

##quad filter
e1_var6=dy.true[(n.end+1):(t-3)]-pred_var6[1:n,1]#1-step ahead forecast error
e2_var6=dy.true[(n.end+2):(t-2)]-pred_var6[1:n,2] #2-step ahead forecast error
e3_var6=dy.true[(n.end+3):(t-1)]-pred_var6[1:n,3]#3-step ahead forecast error
e4_var6=dy.true[(n.end+4):t]-pred_var6[1:n,4] #4-step ahead forecast error
e14_var6=(e1_var6+e2_var6+e3_var6+e4_var6)/4

rmse1_var6=sqrt(mean(e1_var6^2))
rmse2_var6=sqrt(mean(e2_var6^2))
rmse3_var6=sqrt(mean(e3_var6^2))
rmse4_var6=sqrt(mean(e4_var6^2))
rmse14_var6=sqrt(mean(e14_var6^2))
##BN filter
e1_var7=dy.true[(n.end+1):(t-3)]-pred_var7[1:n,1]#1-step ahead forecast error
e2_var7=dy.true[(n.end+2):(t-2)]-pred_var7[1:n,2] #2-step ahead forecast error
e3_var7=dy.true[(n.end+3):(t-1)]-pred_var7[1:n,3]#3-step ahead forecast error
e4_var7=dy.true[(n.end+4):t]-pred_var7[1:n,4] #4-step ahead forecast error
e14_var7=(e1_var7+e2_var7+e3_var7+e4_var7)/4

rmse1_var7=sqrt(mean(e1_var7^2))
rmse2_var7=sqrt(mean(e2_var7^2))
rmse3_var7=sqrt(mean(e3_var7^2))
rmse4_var7=sqrt(mean(e4_var7^2))
rmse14_var7=sqrt(mean(e14_var7^2))
##UC filter
e1_var8=dy.true[(n.end+1):(t-3)]-pred_var8[1:n,1]#1-step ahead forecast error
e2_var8=dy.true[(n.end+2):(t-2)]-pred_var8[1:n,2] #2-step ahead forecast error
e3_var8=dy.true[(n.end+3):(t-1)]-pred_var8[1:n,3]#3-step ahead forecast error
e4_var8=dy.true[(n.end+4):t]-pred_var8[1:n,4] #4-step ahead forecast error
e14_var8=(e1_var8+e2_var8+e3_var8+e4_var8)/4

rmse1_var8=sqrt(mean(e1_var8^2))
rmse2_var8=sqrt(mean(e2_var8^2))
rmse3_var8=sqrt(mean(e3_var8^2))
rmse4_var8=sqrt(mean(e4_var8^2))
rmse14_var8=sqrt(mean(e14_var8^2))

# summary(DMW)
# rolling_samp=20
# tstat=matrix(0,(n-rolling_samp),1)
# for(i in 1:(n-rolling_samp)){
# DMW1=DMW[i:(i+rolling_samp),1]
# DMWols = lm(DMW1~1)
# summary(DMWols)
# tstat[i]=coef(summary(DMWols))[, "t value"]
# }
# plot(tstat)

# 
# plot(c.hp,main="Estimated Cyclical Component",
#      ylim=c(-12,15),col=1,ylab="")
# lines(c.hp3k,col=2)
# lines(c.hp400k,col=3)
# lines(c.hamilton,col=4)
# lines(c.bw,col=5)
# lines(c.linear, col=6)
# lines(c.quad, col=7)
# lines(c.bn, col=8)
# lines(c.uc,col=2, lty =2)
# lines(c.combined, col=3, lty=6)
# lines(diff(credit), col=1, lty=5)
# legend("topleft",legend=c("HP", "HP3k", "HP400k", "Hamilton","BW","linear","quad","BN","UC", "combined", "diff(credit)"),
#        col=c(1:8,2:3,1), lty=c(rep(1,8),2,3,5),ncol=2)

# Baseline forecast is AR(1)

country

# filepath = sprintf('RMSE_ratio_%s.txt',country)
# logfile = file(filepath)
# sink(logfile, append = TRUE, type="output")
# 
# 
# 'Bates-Granger weight combination forecast'
# 'rmse1_varbg/rmse1_ar; rmse2_varbg/rmse2_ar;rmse3_varbg/rmse3_ar; rmse4_varbg/rmse4_ar;rmse14_varbg/rmse14_ar'
# rmse1_varbg/rmse1_ar; rmse2_varbg/rmse2_ar;rmse3_varbg/rmse3_ar; rmse4_varbg/rmse4_ar;rmse14_varbg/rmse14_ar
# 
# 
# 'Average Combination RMSE ratio'
# 'rmse1_varcomb/rmse1_ar; rmse2_varcomb/rmse2_ar;rmse3_varcomb/rmse3_ar; rmse4_varcomb/rmse4_ar; rmse14_varcomb/rmse14_ar'
# rmse1_varcomb/rmse1_ar; rmse2_varcomb/rmse2_ar;rmse3_varcomb/rmse3_ar; rmse4_varcomb/rmse4_ar; rmse14_varcomb/rmse14_ar
# 
# 'BN filter'
# 'rmse1_var7/rmse1_ar; rmse2_var7/rmse2_ar;rmse3_var7/rmse3_ar; rmse4_var7/rmse4_ar; rmse14_var7/rmse14_ar'
# rmse1_var7/rmse1_ar; rmse2_var7/rmse2_ar;rmse3_var7/rmse3_ar; rmse4_var7/rmse4_ar; rmse14_var7/rmse14_ar
# 
# 'HP filter with lambda = 1600'
# 'rmse1_var1/rmse1_ar; rmse2_var1/rmse2_ar;rmse3_var1/rmse3_ar; rmse4_var1/rmse4_ar; rmse14_var1/rmse14_ar'
# rmse1_var1/rmse1_ar; rmse2_var1/rmse2_ar;rmse3_var1/rmse3_ar; rmse4_var1/rmse4_ar; rmse14_var1/rmse14_ar
# 
# 'HP filter with lambda = 3,000'
# 'rmse1_var2/rmse1_ar; rmse2_var2/rmse2_ar;rmse3_var2/rmse3_ar; rmse4_var2/rmse4_ar; rmse14_var2/rmse14_ar'
# rmse1_var2/rmse1_ar; rmse2_var2/rmse2_ar;rmse3_var2/rmse3_ar; rmse4_var2/rmse4_ar; rmse14_var2/rmse14_ar
# 
# 
# 'HP filter with lambda = 400,000'
# 'rmse1_var2/rmse1_ar; rmse2_var2/rmse2_ar;rmse3_var2/rmse3_ar; rmse4_var2/rmse4_ar; rmse14_var2/rmse14_ar'
# rmse1_var3/rmse1_ar; rmse2_var3/rmse2_ar;rmse3_var3/rmse3_ar; rmse4_var3/rmse4_ar; rmse14_var3/rmse14_ar
# 
# 'Hamilton filter'
# 'rmse1_var4/rmse1_ar; rmse2_var4/rmse2_ar;rmse3_var4/rmse3_ar; rmse4_var4/rmse4_ar; rmse14_var4/rmse14_ar'
# rmse1_var4/rmse1_ar; rmse2_var4/rmse2_ar;rmse3_var4/rmse3_ar; rmse4_var4/rmse4_ar; rmse14_var4/rmse14_ar
# 
# 'linear filter'
# 'rmse1_var5/rmse1_ar; rmse2_var5/rmse2_ar;rmse3_var5/rmse3_ar; rmse4_var5/rmse4_ar; rmse14_var5/rmse14_ar'
# rmse1_var5/rmse1_ar; rmse2_var5/rmse2_ar;rmse3_var5/rmse3_ar; rmse4_var5/rmse4_ar; rmse14_var5/rmse14_ar
# 
# 'quad filter'
# 'rmse1_var6/rmse1_ar; rmse2_var6/rmse2_ar;rmse3_var6/rmse3_ar; rmse4_var6/rmse4_ar; rmse14_var6/rmse14_ar'
# rmse1_var6/rmse1_ar; rmse2_var6/rmse2_ar;rmse3_var6/rmse3_ar; rmse4_var6/rmse4_ar; rmse14_var6/rmse14_ar
# 
# 'UC filter'
# 'rmse1_var8/rmse1_ar; rmse2_var8/rmse2_ar;rmse3_var8/rmse3_ar; rmse4_var8/rmse4_ar; rmse14_var8/rmse14_ar'
# rmse1_var8/rmse1_ar; rmse2_var8/rmse2_ar;rmse3_var8/rmse3_ar; rmse4_var8/rmse4_ar; rmse14_var8/rmse14_ar
# 
# 
# closeAllConnections() # Close connection to log file
# # 

name1 <- c("HP", "HP3k", "HP400k", "Hamilton","linear","quad","BN","UC","combined cycle")
c.df=as.data.frame(c.df)
c.weight=as.data.frame(c.weight)
names(c.df)<-name1
name1 <- name1[-length(name1)]
names(c.weight)<-name1

# Graph 1sided cycle data for methods

datef = seq(as.Date("1989-04-01"), by="quarter", length.out=nrow(c.df))
c.df$date = datef
c.df$date <- as.Date(c.df$date)

c.hamilton1s = ts(c.df$Hamilton, start=startdate_diff, freq=4)
# c.bw1s = ts(c.df$BW, start=startdate_diff, freq=4)
c.linear1s = ts(c.df$linear, start=startdate_diff, freq=4)
c.quad1s = ts(c.df$quad, start=startdate_diff, freq=4)
c.bn1s = ts(c.df$BN, start=startdate_diff, freq=4)

min(c.linear1s)

pdfpath = sprintf('Combined_cycle_%s.pdf', country)

pdf(file = pdfpath,   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5.5) # The height of the plot in inches


par(mfrow=c(1,1),mar=c(3,3,2,1),cex=.8)

c.combined = ts(c.combined, start=1989, frequency=4)

maint = sprintf("%s Household Credit Differenced Series & Combined Cycles", country)
plot(dy.true,main=maint, col=1, ylab="", ylim=c(-12,15), lty=5)
lines(c.combined,col=3, lty=6)
lines(c.hp, col=2)
lines(c.hp400k,col=4)
legend("topleft",legend=c("diff(credit)", "combined forecast cycle", "HP cycle component total credit","HP cycle with lambda = 400k (BIS-credit-gap)"),
       col=c(1,3,2,4),lty=c(5,2,1,1),ncol=2)


plot(c.hp,main="Estimated Cyclical Component",
     ylim=c(-15,15),col=1,ylab="")
lines(c.hp3k,col=2)
lines(c.hp400k,col=3)
lines(c.hamilton,col=4)
lines(c.linear, col=5)
lines(c.quad, col=6)
lines(c.bn, col=7)
lines(c.uc,col=2, lty =2)
lines(c.combined, col=3, lty=6)
lines(dy.true, col=1, lty=5)
legend("topleft",legend=c("HP", "HP3k", "HP400k", "Hamilton","linear","quad","BN","UC", "combined", "diff(credit)"),
       col=c(1:7,2:3,1), lty=c(rep(1,7),2,3,5),ncol=2)


c.hpw = ts(c.weight$HP, start=startdate_diff, freq=4)
c.hp3kw=ts(c.weight$HP3k, start=startdate_diff, freq=4)
c.hp400kw=ts(c.weight$HP400k, start=startdate_diff, freq=4)
c.hamiltonw = ts(c.weight$Hamilton, start=startdate_diff, freq=4)
c.linearw = ts(c.weight$linear, start=startdate_diff, freq=4)
c.quadw = ts(c.weight$quad, start=startdate_diff, freq=4)
c.bnw = ts(c.weight$BN, start=startdate_diff, freq=4)
c.ucw = ts(c.weight$UC, start=startdate_diff, freq=4)
dev.off()



pdfpath = sprintf('1-4steps_average_weight_%s.pdf', country)

pdf(file = pdfpath,   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5.5) # The height of the plot in inches

c.weight = na.omit(c.weight)
maint = sprintf('%s Estimated Bates Granger weights of Cyclical Component', country)

plot(c.hpw,main=maint,
          sub="1-4 steps average",
     ylim=c(0.05,0.2),col=1,ylab="")
lines(c.hp3kw,col=2)
lines(c.hp400kw,col=3)
lines(c.hamiltonw,col=4)
lines(c.linearw, col=5)
lines(c.quadw, col=6)
lines(c.bnw, col=7)
lines(c.ucw,col=2, lty =2)
legend("topleft",legend=c("HP", "HP3k", "HP400k", "Hamilton","linear","quad","BN","UC"),
       col=c(1:7,2), lty=c(rep(1,7),2),ncol=2)

dev.off()

filepath = sprintf('GeneratedCycles_%s.csv',country)
c.df$HP<-c.hp
c.df$HP3k<-c.hp3k
c.df$HP400k<-c.hp400k
c.df$Hamilton<-c.hamilton
c.df$linear<-c.linear
c.df$quad<-c.quad
c.df$BN<-c.bn
c.df$UC<-c.uc
write.table(c.df, filepath, sep=',' , row.names = FALSE)

write.table(c.weight, "GeneratedWeights.csv", sep=',' )

#write table append

df_rmse1=as.data.frame(matrix(0,10,7))
names(df_rmse1)=c("Country","filter","1step","2steps","3steps","4steps","avg1-4steps")
df_rmse1[,1]=country

df_rmse1[1,2]="HP"
df_rmse1[1,3:7]=c(rmse1_var1/rmse1_ar, rmse2_var1/rmse2_ar,rmse3_var1/rmse3_ar, rmse4_var1/rmse4_ar, rmse14_var1/rmse14_ar)
df_rmse1[2,2]="HP3k"
df_rmse1[2,3:7]=c(rmse1_var2/rmse1_ar, rmse2_var2/rmse2_ar,rmse3_var2/rmse3_ar, rmse4_var2/rmse4_ar, rmse14_var2/rmse14_ar)
df_rmse1[3,2]="HP400k"
df_rmse1[3,3:7]=c(rmse1_var3/rmse1_ar, rmse2_var3/rmse2_ar,rmse3_var3/rmse3_ar, rmse4_var3/rmse4_ar, rmse14_var3/rmse14_ar)
df_rmse1[4,2]="Hamilton"
df_rmse1[4,3:7]=c(rmse1_var4/rmse1_ar, rmse2_var4/rmse2_ar,rmse3_var4/rmse3_ar, rmse4_var4/rmse4_ar, rmse14_var4/rmse14_ar)
df_rmse1[5,2]="Linear"
df_rmse1[5,3:7]=c(rmse1_var5/rmse1_ar, rmse2_var5/rmse2_ar,rmse3_var5/rmse3_ar, rmse4_var5/rmse4_ar, rmse14_var5/rmse14_ar)
df_rmse1[6,2]="Quad"
df_rmse1[6,3:7]=c(rmse1_var6/rmse1_ar, rmse2_var6/rmse2_ar,rmse3_var6/rmse3_ar, rmse4_var6/rmse4_ar, rmse14_var6/rmse14_ar)
df_rmse1[7,2]="BN"
df_rmse1[7,3:7]=c(rmse1_var7/rmse1_ar, rmse2_var7/rmse2_ar,rmse3_var7/rmse3_ar, rmse4_var7/rmse4_ar, rmse14_var7/rmse14_ar)
df_rmse1[8,2]="UC"
df_rmse1[8,3:7]=c(rmse1_var8/rmse1_ar, rmse2_var8/rmse2_ar,rmse3_var8/rmse3_ar, rmse4_var8/rmse4_ar, rmse14_var8/rmse14_ar)
df_rmse1[9,2]="avgcomb"
df_rmse1[9,3:7]=c(rmse1_varcomb/rmse1_ar, rmse2_varcomb/rmse2_ar,rmse3_varcomb/rmse3_ar, rmse4_varcomb/rmse4_ar, rmse14_varcomb/rmse14_ar)
df_rmse1[10,2]="Bates-Granger"
df_rmse1[10,3:7]=c(rmse1_varbg/rmse1_ar, rmse2_varbg/rmse2_ar,rmse3_varbg/rmse3_ar, rmse4_varbg/rmse4_ar,rmse14_varbg/rmse14_ar)

#write.table append
if (country=="US"){
write.table(df_rmse1, "RMSE_fullsample.csv", sep=',', row.names=FALSE)
}else {
  write.table(df_rmse1, "RMSE_fullsample.csv", sep=',', append=TRUE, row.names=FALSE, col.names = FALSE)
}
