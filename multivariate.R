library(dplyr)
library(zoo)
library(forecast)
library(vars)
library(broom)
library(xts)
library(ggplot2)
library(tseries)
library(CADFtest)
library(fGarch)
library(lubridate)

##### EDA ####
IrInf <- read.csv("D:/Sem 3/TSA/interestRates.csv")
names(IrInf) <- c("Time", "IR", "IF")
head(IrInf)
tail(IrInf)
IrInf$Time <- as.Date(as.yearmon(IrInf$Time,"%Y-%m"))
IrTs <- ts(IrInf$IR, start=c(2011, 12), frequency = 12)
infTs <- ts(IrInf$IF, start=c(2011,12),frequency= 12)
autoplot(cbind(IrTs,infTs), y= "Inflation and Interest Rates")
max.lag<-round(sqrt(length(IrTs)))

#### Stationarity ####
CADFtest(IrTs,type = "trend", criterion = "BIC",max.lag.y =max.lag )
# pvalue > 0.05 therefore non stationary like we saw in plot
CADFtest(infTs,type = "trend", criterion = "BIC",max.lag.y =max.lag )
# pvalue > 0.05 therefore non stationary. Need to take differences.

infTs.diff <- diff(infTs)
irTs.diff <- diff(IrTs)
CADFtest(infTs.diff,type = "trend", criterion = "BIC",max.lag.y = max.lag )
CADFtest(irTs.diff,type = "trend", criterion = "BIC",max.lag.y = max.lag )
# Differencing of first order results into stationarity.

autoplot(cbind(irTs.diff,infTs.diff), y="Differenced Time Series")

inf.sdiff <- diff(diff(infTs, lag = 12)) # Not stationary
ir.sdiff <- diff(diff(IrTs, lag = 12))
autoplot(cbind(inf.sdiff,ir.sdiff), y="Differenced Time Series")

ggmonthplot(inf.sdiff)
ggmonthplot(ir.sdiff)
# pvalue < 0.05 therefore stationary I(1)

#### LR on Diff and non Diff TS ####
fitLm <- lm(IrTs ~ infTs)
summary(fitLm) #Coefficients are significant
ggAcf(fitLm$residuals) # Significant Autocorrelations
Box.test(fitLm$residuals, lag = max.lag, type = "Ljung-Box") # p-value < 5%
# model not valid


fitLmD <- lm(irTs.diff ~ infTs.diff)
summary(fitLmD)
ggAcf(fitLmD$residuals)
Box.test(fitLmD$residuals, lag = max.lag, type = "Ljung-Box")
# No significant AutoCorr and model valid because p val > 5%


#### DLM ####
lag <- 4
n <- length(irTs.diff)
infTs.diff.0 <- infTs.diff[(lag+1):n]
infTs.diff.1 <- infTs.diff.0[lag:(n-1)]
infTs.diff.2 <- infTs.diff[(lag-1):(n-2)]
infTs.diff.3 <- infTs.diff[(lag-2):(n-3)]
infTs.diff.4 <- infTs.diff[(lag-3):(n-4)]

irTs.diff.0 <- irTs.diff[(lag+1):n]
irTs.diff.1 <- irTs.diff[lag:(n-1)]
irTs.diff.2 <- irTs.diff[(lag-1):(n-2)]
irTs.diff.3 <- irTs.diff[(lag-2):(n-3)]
irTs.diff.4 <- irTs.diff[(lag-3):(n-4)]



fit_dlm1 <- lm(infTs.diff.0 ~ irTs.diff.0+irTs.diff.1)
acf(fit_dlm1$residuals) # sig autocorrelations present
Box.test(fit_dlm1$residuals, lag = max.lag, type = "Ljung-Box") # pval < 5% model not valid

fit_dlm1 <- lm(irTs.diff.0 ~ infTs.diff.0+infTs.diff.1)
acf(fit_dlm1$residuals) # sig autocorrelations present
Box.test(fit_dlm1$residuals, lag = max.lag, type = "Ljung-Box") # pval > 5% model valid


fit_dlm2 <- lm(infTs.diff.0 ~ irTs.diff.0+irTs.diff.1 + irTs.diff.2)
acf(fit_dlm2$residuals) # sig autocorrelations present
Box.test(fit_dlm2$residuals, lag = max.lag, type = "Ljung-Box") # pval < 5% model not valid

fit_dlm2 <- lm(irTs.diff.0 ~ infTs.diff.0+infTs.diff.1 +infTs.diff.2)
acf(fit_dlm2$residuals) # sig autocorrelations present
Box.test(fit_dlm2$residuals, lag = max.lag, type = "Ljung-Box") # pval > 5% model not valid


fit_dlm3 <- lm(irTs.diff.0 ~ infTs.diff.0+infTs.diff.1 + infTs.diff.2 + infTs.diff.3)
acf(fit_dlm3$residuals) # sig autocorrelations present
Box.test(fit_dlm3$residuals, lag = max.lag, type = "Ljung-Box") # pval < 5% model not valid

fit_dlm3 <- lm(infTs.diff.0 ~ irTs.diff.0+irTs.diff.1 + irTs.diff.2 + irTs.diff.3)
acf(fit_dlm3$residuals) # sig autocorrelations present
Box.test(fit_dlm3$residuals, lag = max.lag, type = "Ljung-Box") # pval < 5% model not valid

fit_dlm4 <- lm(infTs.diff.0 ~ irTs.diff.0+irTs.diff.1 + irTs.diff.2 + irTs.diff.3 + irTs.diff.4)
acf(fit_dlm4$residuals) # sig autocorrelations present
Box.test(fit_dlm4$residuals, lag = max.lag, type = "Ljung-Box") # pval < 5% model not valid

fit_dlm4 <- lm(irTs.diff.0 ~ infTs.diff.0+infTs.diff.1 + infTs.diff.2 + infTs.diff.3 + infTs.diff.4)
ggAcf(fit_dlm4$residuals) # sig autocorrelations present
Box.test(fit_dlm4$residuals, lag = max.lag, type = "Ljung-Box") # pval > 5% model valid

#### ADLM ####

lag <- 4
n <- length(infTs.diff)
infTs.diff.0 <- infTs.diff[(lag+1):n]
infTs.diff.1 <- infTs.diff.0[lag:(n-1)]
infTs.diff.2 <- infTs.diff[(lag-1):(n-2)]
infTs.diff.3 <- infTs.diff[(lag-2):(n-3)]
infTs.diff.4 <- infTs.diff[(lag-3):(n-4)]
infTs.diff.5 <- infTs.diff[(lag-4):(n-5)]


irTs.diff.0 <- irTs.diff[(lag+1):n]
irTs.diff.1 <- irTs.diff[lag:(n-1)]
irTs.diff.2 <- irTs.diff[(lag-1):(n-2)]
irTs.diff.3 <- irTs.diff[(lag-2):(n-3)]
irTs.diff.4 <- irTs.diff[(lag-3):(n-4)]
irTs.diff.5 <- irTs.diff[(lag-4):(n-5)]
irTs.diff.6 <- irTs.diff[(lag-5):(n-6)]

fit_adlm4 <- lm(infTs.diff.0 ~ infTs.diff.1+infTs.diff.2+infTs.diff.3+infTs.diff.4+
                 irTs.diff.1+irTs.diff.2+irTs.diff.3+irTs.diff.4)
summary(fit_adlm4)
# 14% variability as per R2, Fstat not sign, regressors not jointly signif, infTs.diff.1 is signif

plot.ts(fit_adlm6$residuals)
ggAcf(fit_adlm4$residuals)
Box.test(fit_adlm4$residuals, lag = max.lag, type = "Ljung-Box")

#ADLM4 is valid. ADLM3 not valid. ADLM2 Valid Inf~IR
#ADLM2 not valid. ADLM3 not valid. ADLM4 not valid  IR ~ Inf



#### Granger Causality ####


fit_adlm_nox <- lm(infTs.diff.0 ~ infTs.diff.1+infTs.diff.2+infTs.diff.3+infTs.diff.4)
anova(fit_adlm4,fit_adlm_nox)

fit_adlm_nox <- lm(infTs.diff.0 ~ infTs.diff.1+infTs.diff.2)
anova(fit_adlm3,fit_adlm_nox)

# P Val > 5% Ho not rejected and therefore no explanotory power therefore no GC



#### Engle Granger ####

fit_ci <- lm(infTs ~ IrTs)
res_fit_ci <- fit_ci$residuals
CADFtest(res_fit_ci,type="drift",criterion="BIC",max.lag.y=max.lag)
# -1.94 > -3.41 , tf no cointegration. 
# tf no ecm. 

#### VAR ####

library(vars)
library(gridExtra)
dlogdata<-data.frame(infTs.diff,irTs.diff)
names(dlogdata)<-c("Inf","IR")

fit_var1<-VAR(dlogdata,type="const",p=1)
summary(fit_var1) # R2 =4% pval > 5% tf, regressors not jointly significant
#R2 of IR 12 % pvalue < 5% regressors are jointly significant. 
var1_residuals<-resid(fit_var1)

plt1 <- ggAcf(var1_residuals[,1])
plt2 <- ggAcf(var1_residuals[,2]) # White noise, tf , valid.
plt3 <- ggCcf(var1_residuals[,1],var1_residuals[,2])
grid.arrange(plt1, plt2, plt3 , ncol=2)


VARselect(dlogdata,lag.max=10,type="const")
fit_varautom<-VAR(dlogdata,type="const",p=1)
summary(fit_varautom)
varautom_residuals<-resid(fit_varautom)
par(mfrow=c(2,2))
acf(varautom_residuals[,1])
acf(varautom_residuals[,2])
ccf(varautom_residuals[,1],varautom_residuals[,2])

irf_var<-irf(fit_var1,ortho=F,boot=T)
autoplot(irf_var) # Significant Neg response at t+1 in both cases. 
library(sparsevar)
plot(irf_var)
impulse <- irf_var
lags <- c(1:11)

irf1<-data.frame(impulse$irf$IR[,1],impulse$Lower$IR[,1],
                 impulse$Upper$IR[,1],lags)
irf2<-data.frame(impulse$irf$IR[,2],impulse$Lower$IR[,2],
                 impulse$Upper$IR[,2])
number_ticks <- function(n) {function(limits) pretty(limits, n)}


PIB_PIB <- ggplot(data = irf1,aes(lags,impulse.irf.IR...1.)) +
  geom_line(aes(y = impulse.Upper.IR...1.), colour = 'red') +
  geom_line(aes(y = impulse.Lower.IR...1.), colour = 'red')+
  geom_line(aes(y = impulse.irf.IR...1.))+
  geom_ribbon(aes(x=lags, ymax=impulse.Upper.IR...1., ymin=impulse.Lower.IR...1.), fill="red", alpha=.1) +
  xlab("") + ylab("Inf") + ggtitle("Impulse Response from IR") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),                    
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(2,10,2,10), "mm"))+
  scale_x_continuous(breaks=number_ticks(10)) +
  geom_line(colour = 'black')



PIB_CON <- ggplot(data = irf2,aes(lags,impulse.irf.IR...2.)) +
  geom_line(aes(y = impulse.Upper.IR...2.), colour = 'red') +
  geom_line(aes(y = impulse.Lower.IR...2.), colour = 'red')+
  geom_line(aes(y = impulse.irf.IR...2.))+
  geom_ribbon(aes(x=lags, ymax=impulse.Upper.IR...2., ymin=impulse.Lower.IR...2.), fill="red", alpha=.1) +
  xlab("") + ylab("IR") + ggtitle("") +
  theme(axis.title.x=element_blank(),
        #           axis.text.x=element_blank(),                    
        #           axis.ticks.x=element_blank(),
        plot.margin = unit(c(-10,10,4,10), "mm"))+
  scale_x_continuous(breaks=number_ticks(10)) +
  geom_line(colour = 'black')

grid.arrange(PIB_PIB, PIB_CON, nrow=2)

#### Johansens Test ####
data <- data.frame(infTs,IrTs)
library(urca)
trace_test<-ca.jo(data,type="trace",K=1,ecdet="const",spec="transitory")
summary(trace_test)
maxeigen_test<-ca.jo(logdata,type="eigen",K=4,ecdet="const",spec="transitory")
summary(maxeigen_test)
