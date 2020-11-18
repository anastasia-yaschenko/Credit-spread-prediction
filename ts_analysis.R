library(xts)
library(tseries)
library(forecast)
library(urca)
library(rlist)
library(strucchange)
library(vars)
library(aTSA)
library(tsDyn)
library(ggplot2)
library(xtable)
Sys.setlocale(locale = "C")

wd<-"C:/Users/Dell/Desktop/Time-series project"
setwd(wd)

#Load data
data<-read.csv2("BAMLH0A0HYM2.csv", sep = ",", dec = ".", na.strings = ".")
data<-data[!is.na(data$BAMLH0A0HYM2), ] #delete NA
row.names(data)<-data$DATE
colnames(data)<-c('DATE', 'SPREAD')
spread_for_plot<-xts(data$SPREAD, as.Date(data$DATE))
plot(spread_for_plot, main = 'Credit spread')

sp500<-read.csv2("SP500.csv", sep = ",", dec = ".", na.strings = ".")
sp500<-sp500[!is.na(sp500$SP500), ]
row.names(sp500)<-sp500$DATE
sp_for_plot<-xts(sp500$SP500, as.Date(sp500$DATE))
plot(sp_for_plot, main = 'S&P500')

#data used for testing cointegration
total_df<-merge(data, sp500, by = 0)


#Make test set from last 11 days
hold_out<-data[seq(dim(data)[1]-10, dim(data)[1]), ]
#Train data
spread<-data[seq(1, dim(data)[1]-10),]$SPREAD
spread<-xts(spread, as.Date(data[seq(1, dim(data)[1]-10),]$DATE))

sp<-sp500$SP500

#DATA EXPLORATION
#We can take log of time series to reduce possible heteroscedasticity
log_spread<-log(spread)
log_sp<-log(sp)
log_sp<-xts(log_sp, as.Date(sp500$DATE))

total_df$SPREAD<-log(total_df$SPREAD)
total_df$SP500<-log(total_df$SP500)
hold_out$SPREAD<-log(hold_out$SPREAD)

plot(spread[seq(0,50)], type="l", main = 'Spread')#before transformation
plot(log_spread[seq(0,50)], type="l", main = 'Log(spread)')#after transformation

#Stationarity for log(spread)
#ACF and PACF plots
acf(log_spread, lag.max = 100) #judging by ACF, nonstationary - autocorrelation does not quickly decay to 0
pacf(log_spread, lag.max = 20)

#ADF test
#Choose by BIC since sample is large. Specification without trend but with drift
summary(ur.df(log_spread, type = "drift", lags = 30, selectlags = "BIC"))
#do not reject null hypothesis of non-stationarity->try to difference time series
#KPSS test
summary(ur.kpss(log_spread))#same conclusion

diff_spread<-diff(log_spread)[-1]*100 #multiply by 100 to avoid calculation problems

#ADF test
summary(ur.df(diff_spread, type = "drift", lags = 30, selectlags = "BIC"))
#Now null hypothesis of non-stationarity is rejected
#KPSS test
summary(ur.kpss(diff_spread))#stationary

#Stationarity for S&P500
#ACF and PACF plots
acf(log_sp, lag.max = 100) #judging by ACF, nonstationary - autocorrelation does not quickly decay to 0
pacf(log_sp, lag.max = 20)

#ADF test
summary(ur.df(log_sp, type = "drift", lags = 30, selectlags = "BIC"))
#do not reject null hypothesis of non-stationarity-> difference time series
#KPSS test
summary(ur.kpss(log_sp))#non-stationary

diff_sp<-diff(log_sp)[-1]*100

#ADF test
summary(ur.df(diff_sp, type = "drift", lags = 30, selectlags = "BIC"))
#Now null hypothesis of non-stationarity is rejected
#KPSS test
summary(ur.kpss(diff_sp))#stationary
#So, log(S&P)~I(1) and log(spread)~I(1) 



#MODEL FIT - try to fit ARMA model for differenced log(spread)
acf(diff_spread, main='ACF')#candidate lags = 1, 7, 10- MA part
pacf(diff_spread, main = 'PACF')#candidate lags = 2, 7, 10- AR part
#Estimate candidate models
BICs<-list()
AICs<-list()
AR<-list()
MA<-list()
for (i in c(1, 7, 10)){ 
  for (j in c(2, 7, 10)){
    AR<-list.append(AR, j)
    MA<-list.append(MA, i)
    model<-arima(diff_spread, c(j, 0, i), optim.control = list(maxit=500))
    BICs<-list.append(BICs, c(BIC(model)))
    AICs<-list.append(AICs, c(AIC(model)))
  }
}
BICs<-c(unlist(BICs))
AICs<-c(unlist(AICs))
AR<-c(unlist(AR))
MA<-c(unlist(MA))

models_df<-data.frame(AR, MA, BICs, AICs)
models_df
models_df[which.min(models_df[,"BICs"]),c("AR","MA")] #ARMA(2,1) by BIC
models_df[which.min(models_df[,"AICs"]),c("AR","MA")] #ARMA(7,10) by AIC
#Naturally, BIC chooses smaller number of lags

#Fit auto ARIMA to check if any other model is better
auto_model<-auto.arima(as.numeric(diff_spread), d=0, ic="bic", seasonal = FALSE)
summary(auto_model) #ARMA(2, 0)
BIC(auto_model)
AIC(auto_model)
#This model is better than previous one based on both AIC and BIC. 

#RESIDUAL DIAGNOSTICS - check that residuals are white noise
resid<-auto_model$residuals
acf(resid) #no significant autocorrelation
pacf(resid) #some significant partial autocorrelation
Box.test(resid, lag = 20, type = c('Ljung-Box'), fitdf = 3) #residuals are serially autocorrelated
#Residuals are not white noise, so let us try other models

mod21<-arima(as.numeric(diff_spread), c(2, 0, 1))
resid<-mod21$residuals
acf(resid) #no significant autocorrelation
pacf(resid) #some significant partial autocorrelation
Box.test(resid, lag = 20, type = c('Ljung-Box'), fitdf = 4)#still not white noise

mod710<-arima(as.numeric(diff_spread), c(7, 0, 10), optim.control = list(maxit=500))
resid<-mod710$residuals
acf(resid) #no significant autocorrelation
pacf(resid) #no significant partial autocorrelation
Box.test(resid, lag = 20, type = c('Ljung-Box'), fitdf = 18)#no autocorrelation
Box.test(resid, lag = 19, type = c('Ljung-Box'), fitdf = 18)#no autocorrelation
#This model produces white noise as residuals, so we choose it

# graphs of fitted values
fit <- xts(fitted(mod710), index(diff_spread))
par(mfrow = c(1,1))
plot(xts::last(diff_spread, 100), main = 'ARMA(7,10) Fitted values')
lines(xts::last(fit,100), type = "l", col = "red")

function_errors<-function(err){
  mean_err<-data.frame(c(mean(err^2), mean(abs(err))))
  row.names(mean_err)<-c("MSE", "MAE")
  colnames(mean_err)<-c("value")
  mean_err
}
err_arma<-diff_spread - fit
function_errors(err_arma)#in-sample fit

#FORECASTS
#Forecast for last 20 days
hold_out_spread<-xts(diff(hold_out$SPREAD)*100, as.Date(hold_out$DATE[-1]))

pred<-forecast(mod710, 10)[,2]
pred<-xts(pred, as.Date(hold_out$DATE[-1]))
par(mfrow = c(1,1))
plot(hold_out_spread, main = 'ARMA(7,10) prediction for 10 days')
lines(pred, type = "l", col = "red")
lines(pred*10, type = "l", col = "lightblue")
function_errors(hold_out_spread - pred)#out-of-sample fit
        

#STRUCTURAL BREAKS
all_spread<-diff(log(data$SPREAD))*100
spread_df<-data.frame(all_spread)           
row.names(spread_df)<-data$DATE[-1]

#Structural reak in levels
F_stats<- Fstats(all_spread ~ 1, from=0.15, to=0.85, data = spread_df)
sctest(F_stats)#no break

#Structural break in variance
spread_df$e2<-NA
spread_df$e2[11:nrow(spread_df)]<-(mod710$residuals)^2
F_stats2<- Fstats(e2 ~ 1, from=0.15, to=0.85, data = spread_df)
sctest(F_stats2)#there is a break
dat<-row.names(spread_df)[breakpoints(F_stats2)$breakpoints]
dat
par(mfrow = c(1,1))
plot(spread, main = 'Credit spread')
lines(spread[dat]*100, type='h',col = "red")

#fit new models accounting for break in variance
part1_diff_spread<-diff_spread['/2019-12-02']
part2_diff_spread<-diff_spread['2019-12-02/']

m710_part1<-arima(as.numeric(part1_diff_spread), c(7, 0, 10), optim.control = list(maxit=500))
m710_part2<-arima(as.numeric(part2_diff_spread), c(7, 0, 10), optim.control = list(maxit=500))
m710_part1
m710_part2
#check residuals are white noise
Box.test(m710_part1$residuals, lag = 20, type = c('Ljung-Box'), fitdf = 18)
Box.test(m710_part2$residuals, lag = 20, type = c('Ljung-Box'), fitdf = 18)
#they are white noise for both models

pred2<-forecast(m710_part2, 10)[, 2]
pred2<-xts(pred2, as.Date(hold_out$DATE[-1]))
par(mfrow = c(1,1))
plot(hold_out_spread, main = 'ARMA(7,10) prediction for 10 days (accounting for break in volatility)', ylim=c(min(hold_out_spread), 15))
lines(pred2, type = "l", col = "red") #a little bit better, but still not a good fit
lines(pred2*10, type = "l", col = "lightblue") #a little bit better, but still not a good fit
function_errors(as.numeric(part1_diff_spread) - as.numeric(fitted(m710_part1)))#in-sample fit
function_errors(as.numeric(part2_diff_spread) - as.numeric(fitted(m710_part2)))#in-sample fit
function_errors(as.numeric(hold_out_spread) - as.numeric(pred2))#out-of-sample fit


#VAR model
y1<-diff(total_df$SP500)[-1]*100
y1_cv<-y1[seq((length(y1) - 9), length(y1))]
y1<-y1[seq(1, length(y1) - 9)]

y2<-diff(total_df$SPREAD)[-1]*100
y2_cv<-y2[seq((length(y2) - 9), length(y2))]
y2<-y2[seq(1, length(y2) - 9)]

VARselect(cbind(y1,y2)) #Best model by AIC has 10 lags, Schwarz-2 lags 
var_bic <- VAR(cbind(y1,y2),2)
summary(var_bic)
var_aic <- VAR(cbind(y1,y2),10)
summary(var_aic)

#In-sample fit
var_bic_fit<-fitted(var_bic)[, 2]
var_aic_fit<-fitted(var_aic)[, 2]

err_var_bic<-var_bic_fit - y2[c(-1,-2)]
err_var_aic<-var_aic_fit - y2[seq(-10, -1)]

function_errors(err_var_bic)
function_errors(err_var_aic)


plot_var_predictions<-function(mod){
  var_pred<-predict(mod, 10)$fcst$y2[, 1]
  print(function_errors(y2_cv - var_pred))
  var_pred<-xts(var_pred, as.Date(total_df$DATE.y[seq((nrow(total_df) - 9), nrow(total_df))]))
  y2_cv<-xts(y2_cv, as.Date(total_df$DATE.y[seq((nrow(total_df) - 9), nrow(total_df))]))
  par(mfrow = c(1,1))
  plot(y2_cv, main = 'VAR prediction for 10 days', type = 'l',  ylim=c(min(y2_cv), 10))
  lines(var_pred, type = "l", col = "red") 
  lines(var_pred*10, type = "l", col = "lightblue") 
}
plot_var_predictions(var_bic)#a little bit better, but still not a good fit
plot_var_predictions(var_aic)#more volatile, but still not a good fit
#VAR(2) seems to produce better forecasts, so proceed with it

#Granger causality
causality(var_bic, "y1")$Granger #S&P500 Granger causes spreads
causality(var_bic, "y2")$Granger #spread Granger causes S&P500 at 5% significance level

#cointegration
y3<-total_df$SP500
y4<-total_df$SPREAD
summary(ca.jo(cbind(y3, y4)))#do not reject null hypothesis of no cointegration
coint.test(y3, y4)#do not reject null hypothesis of no cointegration

#Let us try 13-week Treasury bill rate
tres<-read.csv2("^IRX.csv", sep = ",", dec = ".", na.strings = "null")
row.names(tres)<-tres$Date
tres<-subset(tres, select=c("Adj.Close", "Date"))
tres_for_plot<-xts(tres$Adj.Close, as.Date(tres$Date))
plot(tres_for_plot, main = '13 week Treasury bill')
tres$log_close<-log(tres$Adj.Close)
tres<-tres[!is.na(tres$log_close), ]


total_df_tres<-merge(data, tres, by = 0)
y5<-total_df_tres$log_close
y6<-log(total_df_tres$SPREAD)

#check stationarity
#ADF test
summary(ur.df(y5, type = "drift", lags = 30, selectlags = "BIC"))
#do not reject null hypothesis of non-stationarity
#KPSS test
summary(ur.kpss(y5))#non - stationary
y5_diff<-diff(y5)
summary(ur.df(y5_diff, type = "drift", lags = 30, selectlags = "BIC"))
#do not reject null hypothesis of non-stationarity
#KPSS test
summary(ur.kpss(y5_diff))
#So, Treasury rate~I(1)

summary(ca.jo(cbind(y5, y6)))#reject null hypothesis of no cointegration at 5% significance level
coint.test(y5, y6)#do not reject null hypothesis of no cointegration for specification without trend
#We get contradictory results, but let us assume that Johansen test is right and there is cointegration
df_plot<-data.frame(log_close = log(total_df_tres$SPREAD), DATE = total_df_tres$DATE)#, as.Date(total_df_tres$Date))
ggplot(total_df_tres,aes(DATE,log_close,  group = 1))+geom_line(aes(color="13-week Treasury bill"), color='red')+
  geom_line(data=df_plot,aes(color="Credit spread"), color='black')+
  labs(title = "13-week Treasury bill and credit spread")



#Estimate VECM model
y5_cv<-y5[seq((length(y5) - 9), length(y5))]
y5<-y5[seq(1, length(y5) - 9)]
y6_cv<-y6[seq((length(y6) - 9), length(y6))]
y6<-y6[seq(1, length(y6) - 9)]

vecm <- VECM(cbind(y6, y5), 2, estim = "ML")
summary(vecm)

vecm_pred<-predict(vecm, n.ahead = 10)[, 2]
vecm_pred<-xts(vecm_pred, as.Date(total_df_tres$Date[seq((nrow(total_df_tres) - 9), nrow(total_df_tres))]))
y6_cv<-xts(y6_cv, as.Date(total_df_tres$Date[seq((nrow(total_df_tres) - 9), nrow(total_df_tres))]))
par(mfrow = c(1,1))
plot(y6_cv, main = 'VECM prediction for 10 days', type = 'l', ylim=c(-2.5, max(y6_cv)+1))
lines(vecm_pred, type = "l", col = "red") 
vecm_fit<-fitted(vecm)[, 1]
function_errors(vecm_fit - y6[c(-1,-2,-3)])#in-sample fit
function_errors(y6_cv - vecm_pred)#out-of-sample fit
