#Financial Econometrics Final Exam

#Chosen US stock : PerkinElmer, Inc.
#Listed on the NYSE (PKI), PerkinElmer Inc is a company specialized in
#conception of scientific measuring equipments.

'library(quantmod)
library(xts)
library(readr)
library(latex2exp)
library(summarytools)
library(qwraps2)
library(normtest)
library(nortest)
library(moments)
library(xtable)
library(sm)
library(astsa)
library(forecast)
library(portes)'

#STEP1 : collecting data

rm(list=ls()) #clearing environment

#quantmod
PKI<- getSymbols("PKI",from="2009-12-31", to="2019-12-31", auto.assign=FALSE) #extracting data from Yahoo

#zoo
PKI_core<-coredata(PKI) #converting into a matrix
PKI_index<-index(PKI)

#Some observations
first(PKI, "1 week")
last(PKI, "1 week")
periodicity(PKI) 

#Focusing on the adjusted price :
PKI.d<-PKI$PKI.Adjusted #here is our time series
rm(PKI) 
plot(PKI.d, type = 'l') #to observe the evolution of the adjusted price over the 10 last years

#Monthly : starting to register the last day of each month and then construct the dataset
last_day_of_month <- endpoints(PKI.d, on = "months")
PKI.m <-PKI.d[last_day_of_month]
plot(PKI.m, type = 'l') 

#STEP 2 : Let's now compute log-return of our both time series
ss_dates <- "20091231/20191231"
PKI.P.d <- PKI.d[ss_dates]
PKI.P.m <- PKI.m[ss_dates]

PKI.r.d <- diff(log(PKI.P.d))
PKI.r.m <- diff(log(PKI.P.m)) 

#Convert into data.frame
PKI.r.d.df <- cbind(index(PKI.r.d), data.frame(PKI.r.d))
colnames(PKI.r.d.df) <- c("date",'daily log-return')
PKI.r.d.df<-PKI.r.d.df[-1,]
PKI.r.m.df <- cbind(index(PKI.r.m), data.frame(PKI.r.m))
colnames(PKI.r.m.df) <- c("date",'monthly log-return')
PKI.r.m.df<-PKI.r.m.df[-1,]

#Now let's compute the summary table (using moments package as well)

summary(PKI.r.d.df)
summary(PKI.r.m.df)
sd(PKI.r.d.df[,2])
sd(PKI.r.m.df[,2])
skewness(PKI.r.d.df[,2])
skewness(PKI.r.m.df[,2])
kurtosis(PKI.r.d.df[,2])
kurtosis(PKI.r.m.df[,2])
kurtosis.norm.test(PKI.r.d.df[,2])
kurtosis.norm.test(PKI.r.m.df[,2])
quantile(PKI.r.d.df[,2], prob=.05)
quantile(PKI.r.m.df[,2], prob=.05)
quantile(PKI.r.d.df[,2], prob=.95)
quantile(PKI.r.m.df[,2], prob=.95)
jarque.test(PKI.r.d.df[,2])
jarque.test(PKI.r.m.df[,2])
lillie.test(PKI.r.d.df[,2])
lillie.test(PKI.r.m.df[,2])
D1=2*1.96*(0.01519349/sqrt(2263))
D2=2*1.96*(0.06223045/sqrt(108))
D1; D2


#Stylized fact 1 : prices are non-stationary
plot(PKI.d)
plot(PKI.m)
plot(log(PKI.d))
scatter.smooth(x = PKI.P.d[1:2515,],y = PKI.P.d[2:2516,], col=5, xlab = "log-P(t-1)", ylab = "log-P" )

I1<-as.data.frame(log(PKI.P.d))[1:2014,]
I2<-as.data.frame(log(PKI.P.d))[2015:2516,]
mean(I1); sd(I1); mean(I2); sd(I2)

acf(log(PKI.P.d))

#Stylized fact 2 & 3: returns are stationary, asymmetric

par(mfrow=c(1,2))

#log-returns over time
plot(x = PKI.r.d.df[,1], y = PKI.r.d.df[,2], type = 'l', 
     col='red', xlab='Date', ylab='Return', main=TeX('daily log-return PKI ($r_t$)'), yaxp=c(-.15,.15,10))
abline(h = 0, lty=2)
plot(x = PKI.r.m.df[,1], y = PKI.r.m.df[,2], type = 'l', 
     col='red', xlab='Date', ylab='Return', main=TeX('monthly log-return PKI ($r_t$)'), yaxp=c(-.15,.15,10))
abline(h = 0, lty=2)

#log-return histograms
hist(PKI.r.d.df[,2], 100, xlab = "daily log-return", main="", freq = FALSE)
X<-seq(-.2,.2,.01)
lines(x = X, y = (1/(sqrt(2*pi)*0.01551177))*exp(-((X-0.0006442)^2)/(2*(0.01551177^2))), col='red')
hist(PKI.r.m.df[,2], 15, xlab = "monthly log-return",main="", freq=FALSE)
lines(x = X, y = (1/(sqrt(2*pi)*0.06322579))*exp(-((X-0.01350)^2)/(2*(0.06322579^2))), col='red')

#Stylized fact 4 : Heavy tails

#QQ-plots
qqnorm(y = PKI.r.d.df[,2], pch=1, frame=FALSE)
qqline(PKI.r.d.df[,2], col = "steelblue", lwd = 2)
qqnorm(y = PKI.r.m.df[,2], pch=1, frame=FALSE)
qqline(PKI.r.m.df[,2], col = "steelblue", lwd = 2)
#With t-distribution
plot(y = quantile(PKI.r.d.df[,2], prob = seq(0,1,1/length(PKI.r.d.df[,2]))), x = qt(p = seq(0,1,1/length(PKI.r.d.df[,2])),df = 3), xlab='t quantiles, df=3', ylab='daily returns distribution quantiles')
plot(y = quantile(PKI.r.d.df[,2], prob = seq(0,1,1/length(PKI.r.d.df[,2]))), x = qt(p = seq(0,1,1/length(PKI.r.d.df[,2])),df = 5), xlab='t quantiles, df=5', ylab='')
plot(y = quantile(PKI.r.d.df[,2], prob = seq(0,1,1/length(PKI.r.d.df[,2]))), x = qt(p = seq(0,1,1/length(PKI.r.d.df[,2])),df = 10), xlab='t quantiles, df=10', ylab='')

#Stylized fact 5, lilliefors test
sortedreturns<-sort(PKI.r.m.df[,2])
seq.ind <- seq(1,length(sortedreturns),1)
emp.cdf <-seq.ind/length(sortedreturns)
emp.cdf.2 <-(seq.ind-1)/length(sortedreturns)
theor.cdf <- pnorm(sortedreturns,mean(sortedreturns),sd(sortedreturns))

plot(y=emp.cdf,  x=sortedreturns, col="blue", pch=1, lwd=1, lty=1, cex=0.5, 
     xlab='PKI sorted annual log returns', ylab='cdf')
points(y=theor.cdf, x=sortedreturns, col="red", pch=0,lwd=1, lty=1, cex=0.5)
legend("topleft", legend=c("Empirical cdf" , "Normal cdf"), pch=c(1,0), col=c("blue", "red"))
KS.L.stat1 = max(abs(emp.cdf-theor.cdf))
KS.L.stat2 = max(abs(emp.cdf.2-theor.cdf))
KS.L.stat = max(c(KS.L.stat1,KS.L.stat2))
KS.L.stat1
KS.L.stat2
KS.L.stat
plot(y=abs(emp.cdf.2-theor.cdf),x=sortedreturns, type='l',
     xlab='PKI sorted annual log returns', ylim=c(0,0.15), col='blue')
abline(h=0.805/sqrt(length(sortedreturns)),lwd = 4, lty = 1, col='orange')
abline(h=0.886/sqrt(length(sortedreturns)),lwd = 4, lty = 1, col='red')
abline(h=1.031/sqrt(length(sortedreturns)),lwd = 4, lty = 1, col='darkred')
text(x=-0.5, y=0.805/sqrt(length(sortedreturns))-0.006, TeX('10% crit. value $KS_L$ = 0.0983'),  adj = c(0,0), cex= 0.7)
text(x=-0.5, y=0.886/sqrt(length(sortedreturns))+0.002, TeX('5% crit. value $KS_L$ = 0.1082'),  adj = c(0,0), cex= 0.7)
text(x=-0.5, y=1.031/sqrt(length(sortedreturns))+0.002, TeX('1% crit. value $KS_L$ = 0.1259'),  adj = c(0,0), cex= 0.7)
text(x=-0.5, y=0.032, TeX('$ G(\\tilde{x}_t)-\\Phi(\\tilde{x}_t, \\hat{\\mu}, \\hat{\\sigma}^2)$'),  adj = c(0,0))

#Test for aggregational gaussianity
mean(PKI.r.m.df[,2])/mean(PKI.r.d.df[,2])

#Stylized fact 6 : ACF tests

par(mfrow=c(1,2))
lag.max.acf = 40;  lim.y.axes = c(-0.08,0.08)
data2plot = PKI.r.d.df[,2]; # daily returns 
Acf(data2plot, main=TeX('Daily returns : $r_t$'), lag.max = lag.max.acf, xlab = "lag in days", ylim=lim.y.axes, xlim = c(1,40))
data2plot = PKI.r.m.df[,2]; # monthly returns 
Acf(data2plot, main=TeX('Monthly returns : $r_t^m$'), lag.max = lag.max.acf, xlab = "lag in months", ylim=lim.y.axes, xlim = c(1,40))
Acf(data2plot)
my.data <- PKI.r.d.df[,2]
lags.all <- seq(1,25,1);
acf(my.data, lag.max <- max(lags.all), plot = FALSE)
# Barteltt interval
1.96/sqrt(length(my.data))
# test only first lag
Box.test(my.data, lag = 1, type = c("Box-Pierce"), fitdf = 0)
Box.test(my.data, lag = 1, type = c("Ljung-Box"), fitdf = 0)
# test all first 5 lags
Box.test(my.data, lag = 5, type = c("Ljung-Box"), fitdf = 0)
Box.test(my.data, lag = 5, type = c("Box-Pierce"), fitdf = 0)

my.max.lag      <- 25
lags.all        <- seq(1,my.max.lag,1)
my.acf          <- acf(my.data, lag.max = my.max.lag, plot = FALSE)
my.acf.diameter <- qnorm(0.975)/sqrt(length(my.data))
my.acf.tstat.0  <- (my.acf$acf[-1] - 0)/sqrt(1/length(my.data))
my.LjungBox     <- LjungBox(my.data, lags=lags.all)
my.BoxPierce    <- BoxPierce(my.data, lags=lags.all)
crit.value.5.BP <- qchisq(0.95,lags.all) 

my.table <- cbind(my.BoxPierce[,1],
                  my.acf$acf[-1],
                  my.acf.diameter,
                  my.acf.tstat.0,
                  my.BoxPierce[,2],
                  my.BoxPierce[,4],
                  my.LjungBox[,2],
                  my.LjungBox[,4],
                  crit.value.5.BP)
my.table.df <-as.data.frame(my.table)
names(my.table.df)  <- c("lag","acf","acf diam.","acf test","Box-Pierce stat","BP pval","LB stat","LB pval","crit")
rownames(my.table.df) <-c()
options(scipen = 999)
a <- data.matrix(my.table.df)
round(a, digits = 3)

#Stylized fact 7 : Volatility clustering

plot(x = PKI.r.d.df[,1], y = (PKI.r.d.df[,2])^2, type = 'l', 
     col='red', xlab='Date', ylab='Return', main=TeX('daily SQUARED log-return $r_t$^2'), yaxp=c(-.15,.15,10))

rt.d.subsamp <- PKI.r.d[ss_dates];
rt.d.subsamp.df <- cbind(index(rt.d.subsamp), data.frame(rt.d.subsamp)); names(rt.d.subsamp.df)[1] <- "date";

wind.length <- 252; # 1 year of daily data = 252 working days

# mean
roll.mom <- rollapply(data = rt.d.subsamp,
                      width = wind.length,
                      function(y){c(as.numeric(mean(y)),sd(y),skewness(y),kurtosis(y),
                                    moment(y, order = 4, central = TRUE))},
                      align = "right", 
                      by.column=FALSE)
mean.plot    <- roll.mom[,1];
mean.plot.ub <- roll.mom[,1]+1.96*roll.mom[,2]/sqrt(wind.length);
mean.plot.lb <- roll.mom[,1]-1.96*roll.mom[,2]/sqrt(wind.length);
a <-cbind(mean.plot, mean.plot.lb, mean.plot.ub)
mean.plot.all <- cbind(index(a), data.frame(a)); names(mean.plot.all)[1] <- "date";


data2plot <- mean.plot.all[complete.cases(mean.plot.all),]
seq.dates.x <- data2plot[,1]
plot(x = seq.dates.x, y = data2plot[,2]*100, type = 'l', col="blue"   , lty = 1, lwd = 2,
     xlab=""  , ylab=TeX("mean (in percentage)"),    main=TeX('Rolling mean(on 252 days)  \\%'),   xaxt ="none",
     ylim=c(-0.5,0.5))
lines(x=seq.dates.x, y=data2plot[,3]*100, col="red", lwd=1, lty=1)
lines(x=seq.dates.x, y=data2plot[,4]*100, col="red", lwd=1, lty=1)
seq_sel  <- endpoints(data2plot$date, on = 'years'); 
date_seq <- data2plot$date[seq_sel]; 
date_lab <- format(date_seq,"%y")
axis(1, at = date_seq, label = date_lab, las = 1, cex.axis=1.0); abline(0,0, lty = 1) # add zero line
dev.off()

#St. Dev
sd.plot    <- roll.mom[,2];
mu4        <- roll.mom[,5];
sd.plot.ub <- roll.mom[,2]+1.96*(1/(2*sd.plot)*sqrt(mu4-sd.plot^4))/sqrt(wind.length);
sd.plot.lb <- roll.mom[,2]-1.96*(1/(2*sd.plot)*sqrt(mu4-sd.plot^4))/sqrt(wind.length);
a <-cbind(sd.plot,sd.plot.lb,sd.plot.ub)
sd.plot.all <- cbind(index(a), data.frame(a)); names(sd.plot.all)[1] <- "date";

data2plot <- sd.plot.all[complete.cases(sd.plot.all),]
#postscript(file="../figures/SP500_stdev_rolling_1981_2018.eps",width=12,height=6,horizontal = FALSE, onefile=FALSE) 
plot(x = seq.dates.x, y = data2plot[,2]*100, type = 'l', col="blue"   , lty = 1, lwd = 2,
     xlab=""  , ylab=TeX("st.dev (in percentage)"),    main=TeX('Rolling st. dev.(on 252 days)  \\%'),   xaxt ="none",
     ylim=c(0,4))     # do not diaply x and y-axes labels/ticks) 
lines(x=seq.dates.x, y=data2plot[,3]*100, col="red", lwd=1, lty=1)
lines(x=seq.dates.x, y=data2plot[,4]*100, col="red", lwd=1, lty=1)
axis(1, at = date_seq, label = date_lab, las = 1, cex.axis=1.0); abline(0,0, lty = 1) # add zero line
dev.off()


#ACF r2 and |r|
par(mfrow=c(1,2))
lag.max.acf = 80;  lim.y.axes = c(-0.20,0.20)
data2plot = abs(PKI.r.d.df[,2]) ; 
Acf(data2plot, main=TeX('Daily absolute returns : $|r_t|$'), lag.max = lag.max.acf, xlab = "lag in days", ylim=lim.y.axes)
data2plot2 = abs(PKI.r.m.df[,2]) ;
Acf(data2plot2, main=TeX('Monthly absolute returns : $|r_t|$'), lag.max = lag.max.acf, xlab = "lag in months", ylim=lim.y.axes)

data2plot3 = (PKI.r.d.df[,2])^2 ; 
Acf(data2plot3, main=TeX('Daily Squared returns : $r_t$^2'), lag.max = lag.max.acf, xlab = "lag in days", ylim=lim.y.axes)
data2plot4 = (PKI.r.m.df[,2])^2 ;
Acf(data2plot4, main=TeX('Monthly Squared returns : $r_t$^2'), lag.max = lag.max.acf, xlab = "lag in months", ylim=lim.y.axes)


#STYLIZED FACT 8

ret   = (PKI.r.d.df[,2])   ; # daily returns 
ret2  = (PKI.r.d.df[,2])^2 ; # daily SQUARED returns 

ret  <- as.numeric(rt.d.subsamp.df); ret2 <- ret^2
ccf(ret, ret2, lag.max = 10, type = "correlation",  plot = TRUE, 
    main=TeX('Cross-coorrelation between daily $r_{t+j}$ and $r_t^2$ = corr($r_{t+j}$, $r_{t}^2$)'), 
    xlab = TeX('lag $j$ in days'), ylab=TeX('Cross-correlation'))

#Comparing with VIX index
VIX      <- getSymbols("^VIX",from="1951-12-31", to="2018-12-31", auto.assign=FALSE) 
VIXcsv   <-  as.xts(read.zoo("Data_VIX.csv", header = TRUE, format = "%Y-%m-%d"))
VIX.d.all  <- VIX$VIX.Adjusted ; names(VIX.d.all)  <- "VIX.d"

a <- merge(PKI.P.d,VIX.d.all) # merge two datasets (log prices and VIX ) with different legths --> NA are kept
b <- diff(a)               # compute changes in pt and VIX conmpared to previus period

a.df <- cbind(index(a), data.frame(a)); names(a.df)[1] <- "date";
b.df <- cbind(index(b), data.frame(b)); names(b.df)[1] <- "date";

a.comp <- a.df[complete.cases(a.df),] # remove all rows in which there is at least one NA -> generate balanced (complete) matrix
b.comp <- b.df[complete.cases(b.df),]

head(a.comp); tail(a.comp)
head(b.comp); tail(b.comp)

par(mfrow=c(1,1))
data2plot <- a.comp;
plot(x = data2plot[,1], y = data2plot[,2], type ="l", ylab=TeX("PKI (blue)"),    main=TeX('PKI (left, blue) and VIX (right, red) indexes'),   
     xlab="", xaxt ="none",   col = "blue", lwd =2)
par(new = TRUE)
plot(x = data2plot[,1], y = data2plot[,3], type = "l", xaxt = "n", yaxt = "n",
     ylab = "", xlab = "", col = "red", lty = 1)
axis(side = 4)
mtext("VIX", side = 4, line = 3)
seq_sel <- endpoints(data2plot$date, on = 'years'); date_seq = data2plot$date[seq_sel]; date_lab = format(date_seq,"%b-%y")
axis(1, at = date_seq, label = date_lab, las = 1, cex.axis=1.0)

par(mfrow=c(1,1))
data2plot <- b.comp;
plot(x = data2plot[,2], y = data2plot[,3], xlab=TeX("PKI log return"),    main=TeX('($VIX_t-VIX_{t-1}$) vs. $r_t$'),   
     ylab="VIX Change", col = "blue", lwd =1)
abline(lm(data2plot[,3]~data2plot[,2]), col="red", lwd=2) # regression line (y~x) 

#END