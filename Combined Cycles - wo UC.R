## Calculate MSE
### Load variable

# c.hp<- credit.hp$cycle
# c.bk<- credit.bk$cycle
# c.cf<- c.credit.cf$cycle
# c.bw<- credit.bw$cycle
# c.tr<- credit.tr$cycle
# c.uc<- credit.uc

dy = diff(credit)
startdate = 1989

dy = window(dy, start=startdate, end=2016.75)
#dy.true=diff(credit)
dy.true=dy

library(forecast)
library(vars)
library(tseries)
library(dynlm)
library(zoo)
library(dyn)
#### y = HPI
#### c = Credit to household cyclical component
###### There are 6 cycles components
##### delta_y = alpha + beta_i c_{t-1} + delta_y_{t-1} + v_t
# m1 <- lm(delta_y ~ c_i + lag.delta_y)
## x6

var_1=ts(cbind(dy,c.hp), start=startdate, frequency =4)
var_2=ts(cbind(dy,c.bk), start=startdate, frequency =4)
var_3=ts(cbind(dy,c.cf), start=startdate, frequency =4)
var_4=ts(cbind(dy,c.bw), start=startdate, frequency =4)
var_5=ts(cbind(dy,c.tr), start=startdate, frequency =4)

n.end=20 #Initial sample estimation
t=length(dy) #Full Sample size
n=t-n.end-3 #Forecast sample

# set matrix for storage
pred_var1=matrix(0,n,4)
pred_var2=matrix(0,n,4)
pred_var3=matrix(0,n,4)
pred_var4=matrix(0,n,4)
pred_var5=matrix(0,n,4)

pred_ar=matrix(0,n,4)
pred_comb=matrix(0,n,4)

w_bg1=matrix(0,n,7)#need separate weights for each horizon in Bates-Granger method
w_bg2=matrix(0,n,7)
w_bg3=matrix(0,n,7)
w_bg4=matrix(0,n,7)

for(i in 1:n){
  
  #### Run regression
  x_var1 = var_1[1:(n.end+i-1),]
  x_var2 = var_2[1:(n.end+i-1),]
  x_var3 = var_3[1:(n.end+i-1),]
  x_var4 = var_4[1:(n.end+i-1),]
  x_var5 = var_5[1:(n.end+i-1),]

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
  
  
  model.ar=arima(x_ar,order=c(1,0,0),method="ML")
  pred_ar[i,1:4]=predict(model.ar,n.ahead=4,se.fit=FALSE)[1:4]
  
  ###Bates-Granger
  sam_size=n.end+i-1
  
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
  
  wbgstr_11=1/mse1_1
  wbgstr_12=1/mse1_2
  wbgstr_13=1/mse1_3
  wbgstr_14=1/mse1_4
  wbgstr_15=1/mse1_5
  
  
  wbgstr_21=1/mse2_1
  wbgstr_22=1/mse2_2
  wbgstr_23=1/mse2_3
  wbgstr_24=1/mse2_4
  wbgstr_25=1/mse2_5
  
  
  wbgstr_31=1/mse3_1
  wbgstr_32=1/mse3_2
  wbgstr_33=1/mse3_3
  wbgstr_34=1/mse3_4
  wbgstr_35=1/mse3_5
  
  
  wbgstr_41=1/mse4_1
  wbgstr_42=1/mse4_2
  wbgstr_43=1/mse4_3
  wbgstr_44=1/mse4_4
  wbgstr_45=1/mse4_5
  
  
  wbgstr1_sum=(wbgstr_11+wbgstr_12+wbgstr_13+wbgstr_14+wbgstr_15)
  wbgstr2_sum=(wbgstr_21+wbgstr_22+wbgstr_23+wbgstr_24+wbgstr_25)
  wbgstr3_sum=(wbgstr_31+wbgstr_32+wbgstr_33+wbgstr_34+wbgstr_35)
  wbgstr4_sum=(wbgstr_41+wbgstr_42+wbgstr_43+wbgstr_44+wbgstr_45)
  
  
  wbg_11=wbgstr_11/wbgstr1_sum
  wbg_12=wbgstr_12/wbgstr1_sum
  wbg_13=wbgstr_13/wbgstr1_sum
  wbg_14=wbgstr_14/wbgstr1_sum
  wbg_15=wbgstr_15/wbgstr1_sum  

  
  wbg_21=wbgstr_21/wbgstr2_sum
  wbg_22=wbgstr_22/wbgstr2_sum
  wbg_23=wbgstr_23/wbgstr2_sum
  wbg_24=wbgstr_24/wbgstr2_sum
  wbg_25=wbgstr_25/wbgstr2_sum  

  
  wbg_31=wbgstr_31/wbgstr3_sum
  wbg_32=wbgstr_32/wbgstr3_sum
  wbg_33=wbgstr_33/wbgstr3_sum
  wbg_34=wbgstr_34/wbgstr3_sum
  wbg_35=wbgstr_35/wbgstr3_sum  

  
  wbg_41=wbgstr_41/wbgstr4_sum
  wbg_42=wbgstr_42/wbgstr4_sum
  wbg_43=wbgstr_43/wbgstr4_sum
  wbg_44=wbgstr_44/wbgstr4_sum
  wbg_45=wbgstr_45/wbgstr4_sum  

  
  w_bg1[i,1:5]=c(wbg_11,wbg_12,wbg_13,wbg_14,wbg_15)
  w_bg2[i,1:5]=c(wbg_21,wbg_22,wbg_23,wbg_24,wbg_25)
  w_bg3[i,1:5]=c(wbg_31,wbg_32,wbg_33,wbg_34,wbg_35)
  w_bg4[i,1:5]=c(wbg_41,wbg_42,wbg_43,wbg_44,wbg_45) 
  
## Calcalte DMW time series values
### Squared error of Model A: AR(1) forecast  
e_sqr_ar1=(dy.true[(i+1)]-pred_ar[i,1])^2#1-step ahead forecast error
### Squared error of Model B: combination avarage forecast
pred_c1=(pred_var1[i,1]+pred_var2[i,1]+pred_var3[i,1]+pred_var4[i,1]+pred_var5[i,1])/5
e_sqr_c1=(dy.true[(i+1)]-pred_c1)^2
### DMW time series:
DMW[i]=e_sqr_c1-e_sqr_ar1
}

### Extract MSE from lm function

## Create combined cyclical components
### Save all MSE and calculate normalized weights
### Create combined cyclical component
#### Read notes and problem set from 735: combined forecast
pred_bg1=(w_bg1[,1]*pred_var1[,1]+w_bg1[,2]*pred_var2[,1]
          +w_bg1[,3]*pred_var3[,1]+w_bg1[,4]*pred_var4[,1]+w_bg1[,5]*pred_var5[,1])
pred_bg2=(w_bg2[,1]*pred_var1[,2]+w_bg2[,2]*pred_var2[,2]
          +w_bg2[,3]*pred_var3[,2]+w_bg2[,4]*pred_var4[,2]+w_bg2[,5]*pred_var5[,2])

pred_bg3=(w_bg3[,1]*pred_var1[,3]+w_bg3[,2]*pred_var2[,3]
          +w_bg3[,3]*pred_var3[,3]+w_bg3[,4]*pred_var4[,3]+w_bg3[,5]*pred_var5[,3])

pred_bg4=(w_bg4[,1]*pred_var1[,4]+w_bg4[,2]*pred_var2[,4]
          +w_bg4[,3]*pred_var3[,4]+w_bg4[,4]*pred_var4[,4]+w_bg4[,5]*pred_var5[,4])
pred_bg=cbind(pred_bg1,pred_bg2,pred_bg3,pred_bg4)

pred_c=(pred_var1+pred_var2+pred_var3+pred_var4+pred_var5)/5


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

# 
# summary(DMW)
# rolling_samp=20
# tstat=matrix(0,(n-rolling_samp),1)
# for(i in 1:(n-rolling_samp)){
#   DMW1=DMW[i:(i+rolling_samp),1]
#   # DMW1 = as.data.frame(DMW1)
#   DMWols = lm(DMW1~1)
#   summary(DMWols)
#   tstat[i]=coef(summary(DMWols))[, "t value"]
# }
# plot(tstat)

rmse1_varcomb/rmse1_ar; rmse2_varcomb/rmse2_ar;rmse3_varcomb/rmse3_ar; rmse4_varcomb/rmse4_ar; rmse14_varcomb/rmse14_ar
rmse1_varbg/rmse1_ar; rmse2_varbg/rmse2_ar;rmse3_varbg/rmse3_ar; rmse4_varbg/rmse4_ar;rmse14_varbg/rmse14_ar
