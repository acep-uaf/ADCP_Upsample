#install.packages('Langevin')
#install.packages('plotrix')
#install.packages('pracma')
#install.packages('tseries')
library(tseries)
library(ggplot2)
library(dplyr)
library("pracma")
library("Langevin")
library("plotrix")

#velData <- read.csv(file = 'vel_4.csv')
velData <- read.csv(file = '230713ad2cp_swvel4.csv')

head(velData)


# nome data is in 10 second time steps for 1 year. 

fs <- 1
#x <- slice(velData, 1:45000)
x <- velData[,1] 
#x <- x[1:45000]  # remove load values that do not have many values and throw analysis off

# the Nome dataset needs to be made stationary: possibly by subtracting a moving average from the data
# this will not result in different values dependent on load. 

# try dividing a moving average 
x_window <- 10 # 
x_total <- 45000
x_movavg <- movavg(x, x_window,type = 's')
#x_movavg <- c(x_movavg[round(x_window/2):length(x)],x[(length(x)-x_window+round(x_window/2)+2):length(x)]) # center the moving average
x_movavg <- x_movavg[round(x_window/2):(x_total+round(x_window/2)-1)] # center the moving average
x <- x[1:x_total]
x_std = sd(x_movavg)
x_minus_movavg <- (x - x_movavg) / x_std

#x_mean <- mean(x)
#x_std <- sd(x)
#x_minus_movavg <- (x - x_mean) / x_std

# test
par(mfrow = c(1, 1))
plot(x[1:1000], type = 'l', ylab = "Velocity (m/s)")
lines(x_movavg[1:1000],col='red')
legText = c("Data", "Moving Average")
legend(x = "topright", legend = legText, lty = c(1, 1), col = c('black', 'red'))

plot(x[44000:45000],type = 'l', ylab = "Velocity (m/s)")
lines(x_movavg[44000:45000],col='red')
legText = c("Data", "Moving Average")
legend(x = "topright", legend = legText, lty = c(1, 1), col = c('black', 'red'))

plot(x_minus_movavg[1:1000],type = 'l', ylab = "", main = "Normalized Difference from Moving Average")
plot(x_minus_movavg[44000:45000],type = 'l', ylab = "", main = "Normalized Difference from Moving Average")

plot(x,type = 'l', ylab = "Velocity (m/s)")
lines(x_movavg,col='red')
legText = c("Data", "Moving Average")
legend(x = "topright", legend = legText, lty = c(1, 1), col = c('black', 'red'))

plot(x_minus_movavg,type = 'l', ylab = "", main = "Normalized Difference from Moving Average")

x_den <- density(x_minus_movavg)
plot(x_den, main = "Probablility Density of Normalized Difference from Moving Average")

# calculate langevin coefficients 
bins <- 20 # number of bins to divide the velocities into 
steps <- c(1)#:2) # vector of steps to calculate conditional moments for different tau values (change in time). I am not entirely sure what effect this has.
# since we want 10 second data, I believe we are only interested in a time step of 1. But, I ran time steps of 1 - 10. 
ests_detrended <- Langevin1D(x_minus_movavg, bins, steps, sf=fs)
summary(ests_detrended)
plot(ests_detrended)

attach(ests_detrended)
par(mfrow = c(1, 2))
plotCI(mean_bin, D1, uiw = eD1, xlab = "x [a.u.]",
       ylab = "", cex = 2, pch = 20)
ylab.text = expression(paste("Drift coefficient ", D^(1), "(x) [a.u.]"))
mtext(ylab.text, side = 2, line = 2.5)
plotCI(mean_bin, D2, uiw = eD2, xlab = "x [a.u.]",
       ylab = "", cex = 2, pch = 20)
ylab.text = expression(paste("Diffusion coefficient ", D^(2), "(x) [a.u.]"))
mtext(ylab.text, side = 2, line = 2.5)

#linearModD1 <- lm(D1 ~ mean_bin + I(mean_bin^2) + I(mean_bin^3), weights = 1/eD1) # linear with no offset
linearModD1 <- lm(D1 ~ mean_bin, weights = 1/eD1) # linear with no offset
summary(linearModD1)
linearModD2 <- lm(D2 ~ mean_bin + I(mean_bin^2), weights = 1/eD2 ) # quadratic with no linear
summary(linearModD2)

plot.new()
par(mfrow = c(1, 2))
plotCI(mean_bin, D1, uiw = eD1, xlab = "x [a.u.]",
       ylab = "", cex = 2, pch = 20)
ylab.text = expression(paste("Drift coefficient ", D^(1), "(x) [a.u.]"))
mtext(ylab.text, side = 2, line = 2.5)
lines(mean_bin[as.numeric(names(predict(linearModD1)))],predict(linearModD1,na.rm = FALSE),col = 'red')
plotCI(mean_bin, D2, uiw = eD2, xlab = "x [a.u.]",
       ylab = "", cex = 2, pch = 20)
ylab.text = expression(paste("Diffusion coefficient ", D^(2), "(x) [a.u.]"))
mtext(ylab.text, side = 2, line = 2.5)
lines(mean_bin[as.numeric(names(predict(linearModD2)))],predict(linearModD2,na.rm = FALSE),col = 'red')

# generate time series from langevin coefficients
par(mfrow = c(1, 1))
set.seed(11)
calcVals <- timeseries1D(length(x), startpoint = 0,
                        #d13 = coefficients(linearModD1)[4],
                        #d12 = coefficients(linearModD1)[3],
                         d13 = 0,
                         d12 = 0,
                         d11 = coefficients(linearModD1)[2],
                         d10 = coefficients(linearModD1)[1],
                         d22 = coefficients(linearModD2)[3], 
                         d21 = coefficients(linearModD2)[2],
                        #d20 = 0, sf=fs)
                         d20 = coefficients(linearModD2)[1], sf=fs)
plot(calcVals)
#calcVals <- timeseries1D(length(x), d11 = -0.3538739,
#                         d22 = 0.3511254704, d21 = 0,
#                         d20 = 0.0001702648, sf=fs)

t <- 1:length(x)
plot(x_minus_movavg[1:1000], type = 'l', ylab = "", main = "Normalized Difference from Moving Average")
lines(calcVals[1:1000], col='orange')
legText = c("Measured", "Generated")
legend(x = "topright", legend = legText, lty = c(1, 1), col = c('black', 'red'))

plot(x_minus_movavg[44000:45000], type = 'l', ylab = "", main = "Normalized Difference from Moving Average")
lines(calcVals[44000:45000], col='red')
plot(x_minus_movavg, type = 'l', ylab = "", main = "Normalized Difference from Moving Average")
lines(calcVals, col='red')

x_calc <- x_std * calcVals + x_movavg
#x_calc <- x_std * calcVals + x_mean

plot(x_calc,type = 'l')
plot(x[1:1000], type = 'l')
lines(x_calc[1:1000],col='red')
#lines(x[1:1000],col='blue')
# this looks pretty good. Real data will only have average values every 10 min, so it may be beneficial to rerun analysis using average value every 10 min with linear 
# interpolations. 

x_fft <- fft(x_minus_movavg)
calc_fft <- fft(calcVals)

plot(Mod(x_fft), type = 'l')
plot(Mod(calc_fft), type = 'l')

# test if time series is stationary 
adf.test(x_minus_movavg) # p-value < 0.05 indicates the TS is stationary
kpss.test(x_minus_movavg) # p-value > 0.05 indicates TS is stationary 





###### try to upsample to 1 second, from 10 second base data ##########

# try to calculate 1 second values with sf = 0.1
#calcVals_1s = timeseries1D(length(x), startpoint = 0,d13 = 0,d12 = 0, d11 = coefficients(linearModD1)[1],
#                           d10 = 0,
#                           d22 = coefficients(linearModD2)[2], d21 = 0,
#                           d20 = coefficients(linearModD2)[1], sf=10)

#numSeconds = 3000
#plot(linspace(1,numSeconds,numSeconds/10),x_minus_movavg[1:(numSeconds/10)],type = 'l')
#lines(linspace(1,numSeconds,numSeconds/10),calcVals[1:(numSeconds/10)],col='red')
#lines(linspace(1,numSeconds,numSeconds),calcVals_1s[1:numSeconds],col='blue')
