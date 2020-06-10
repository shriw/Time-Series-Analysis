library(Quandl)
library(dplyr)
library(xts)
library(ggplot2)
library(tseries)
library(forecast)
library(CADFtest)
library(fGarch)
library(lubridate)
# Univariate Data #

data <- read.csv("D:/Sem 3/TSA/sharePrices.csv")

#### EDA ####
head(data)
tail(data)
summary(data)
names(data) <- c("Date","SP")
data <- arrange(data,Date)
frequency(data)
data$Date <- as.Date(as.yearmon(data$Date, format="%Y-%m"))

head(data)

ts <- ts(data$SP,start = c(2001,4) ,frequency = 12)
autoplot(ts, xlab="Year", ylab = "Indian Share Price Index, 2015=100")
tsLog <- log(ts)
plot(tsLog)


#tsLog[which(!is.finite(tsLog))] <- NA
#tsLog <- tsLog[complete.cases(tsLog)]
#plot(tsLog) # trend, No Seasonality, non stationary, Variance is not exactly constant. Could be heteroskedastic

ts.components <- decompose(tsLog)
plot(ts.components)

# check for stationarity
adf.test(ts,alternative = "stationary")
# Failed to reject Null hypothesis. Data is non stationary.

max.lag<-round(sqrt(length(ts)))
CADFtest(ts, type= "trend", criterion= "BIC", max.lag.y=max.lag)
#pvalue is greater than 0.5 therefore fail to reject Ho and therefore stochastic trend


ts.diff <- diff(ts)
head(tsLog.diff)
autoplot(ts.diff)

adf.test(ts.diff,alternative = "stationary")
CADFtest(ts.diff, type= "trend", criterion= "BIC", max.lag.y=max.lag)
# Data is stationary

ggmonthplot(ts.diff)
ggsubseriesplot(ts.diff)
ggseasonplot(ts)


#Seasonality
ts.sdiff <- diff(diff(ts, lag=12))
autoplot(ts.sdiff)
ggmonthplot(ts.sdiff)
#ggseasonplot(ts.sdiff)
CADFtest(ts.sdiff, type= "trend", criterion= "BIC", max.lag.y=max.lag)

ggAcf(ts.sdiff)
ggPacf(ts.sdiff)

# Test for white noise
max.lag<-round(sqrt(length(tsLog.diff)))
Box.test(ts.sdiff, lag = max.lag, type = "Ljung-Box") #p-value is less than 0.05 therefore not white noise

#### Models #####

tsTest <- tail(ts, 10)
tsTrain <- head(ts, 214)

ts.arima <- auto.arima(tsTrain)

ts.arima1 <- Arima(tsTrain,c(2,1,1),c(0,1,0) )
ts.plot(ts.arima1$residuals) 
ggAcf(ts.arima1$residuals) # significant auto correlation at 12
Box.test(ts.arima1$residuals ,lag = max.lag, type="Ljung-Box") # p-val<0.05 model is not valid
ts.arima1 %>% forecast(h=10) %>% autoplot() +autolayer(tsTest)
tsfcst <- ts.arima1 %>% forecast(h=10)
accuracy(tsfcst, tsTest)


ts.arima2 <- Arima(tsTrain, c(2,1,0), c(0,1,0))
#ts.plot(ts.arima2$residuals) 
ggAcf(ts.arima2$residuals) # few autocorrelations, not significant
Box.test(ts.arima2$residuals ,lag = max.lag, type="Ljung-Box") # p-val<0.05 model is not valid
ts.arima2 %>% forecast(h=10) %>% autoplot() +autolayer(tsTest)
tsfcst <- ts.arima2 %>% forecast(h=10)
accuracy(tsfcst, tsTest)

ts.arima3 <- Arima(tsTrain, c(2,1,1),c(0,1,1))
ts.plot(ts.arima3$residuals) 
ggAcf(ts.arima3$residuals) # no significant auto correlation
Box.test(ts.arima3$residuals ,lag = max.lag, type="Ljung-Box") # p-val>0.05 model is valid
ts.arima3 %>% forecast(h=10) %>% autoplot() +autolayer(tsTest)
tsfcst <- ts.arima3 %>% forecast(h=10)
accuracy(tsfcst, tsTest)


ts.arima4 <- Arima(tsTrain, c(2,1,1), c(1,1,1))
ts.plot(ts.arima4$residuals) 
ggAcf(ts.arima4$residuals) # no significant auto correlation
Box.test(ts.arima3$residuals ,lag = max.lag, type="Ljung-Box") # p-val>0.05 model is valid
ts.arima4 %>% forecast(h=10) %>% autoplot() +autolayer(tsTest)
tsfcst <- ts.arima4 %>% forecast(h=10)
accuracy(tsfcst, tsTest)

ts.arima5 <- Arima(tsTrain, c(1,1,2), c(0,1,1))
ts.plot(ts.arima5$residuals) 
ggAcf(ts.arima5$residuals) # no significant auto correlation
Box.test(ts.arima5$residuals ,lag = max.lag, type="Ljung-Box") # p-val>0.05 model is valid
ts.arima5 %>% forecast(h=10) %>% autoplot() +autolayer(tsTest)
tsfcst <- ts.arima5 %>% forecast(h=10)
accuracy(tsfcst, tsTest)

ts.arima6 <- Arima(tsTrain, c(2,1,1))
#ts.plot(ts.arima5$residuals) 
ggAcf(ts.arima6$residuals) # no significant auto correlation
Box.test(ts.arima6$residuals ,lag = max.lag, type="Ljung-Box") # p-val<0.05 model is not valid
ts.arima6 %>% forecast(h=10) %>% autoplot() +autolayer(tsTest) 
tsfcst <- ts.arima6 %>% forecast(h=10)
accuracy(tsfcst,tsTest)


dm.test(residuals(ts.arima3), residuals(ts.arima4), h=10)  #
dm.test(residuals(ts.arima3), residuals(ts.arima5), h=10)
dm.test(residuals(ts.arima3), residuals(ts.arima6), h=10)
dm.test(residuals(ts.arima4), residuals(ts.arima5), h=10)
dm.test(residuals(ts.arima4), residuals(ts.arima6), h=10)
dm.test(residuals(ts.arima5), residuals(ts.arima6), h=10)


#### Heteroscedasticity ####
ggAcf(ts.arima3$residuals)
ggAcf(ts.arima6$residuals^2)


library(tseries)
tsgarch <- garch(tsTrain,trace = F)
tsgarch.res <- na.remove(tsgarch$res[-1])
acf(tsgarch.res)
acf(tsgarch.res^2)




fitGarch <- garchFit(~garch(1,1), data = tsTrain)
summary(fitGarch) # all except bet1 significant,  
plot(fitGarch) # Not normal, Acf shows autocorrelations model invalid

fitGarch <- garchFit(~arma(2,1)+garch(1,1), data = tsTrain) #, cond.dist = "QMLE")
summary(fitGarch)
plot(fitGarch) # Standardized and squared standardized residuals have no autocorrelation hence model is valid
# few periods of volatility 

fitGarch <- garchFit(~arma(2,1)+garch(1,1), data = tsTrain, cond.dist = "QMLE")
summary(fitGarch)
 plot(fitGarch)

fgarchfcst <- predict(fitGarch, n.ahead=10)
meanfgarchfcst <- fgarchfcst$meanForecast
accuracy(meanfgarchfcst,tsTest)
