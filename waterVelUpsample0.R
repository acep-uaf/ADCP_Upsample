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

# the dataset needs to be made stationary: possibly by subtracting a moving average from the data
# this will not result in different values dependent on load. 

# try dividing a moving average 
x_window <- 10 # 
x_total <- 57600
end_window <- c((x_total-1000): x_total)
x_movavg <- movavg(x, x_window,type = 's')
#x_movavg <- c(x_movavg[round(x_window/2):length(x)],x[(length(x)-x_window+round(x_window/2)+2):length(x)]) # center the moving average
x_movavg <- x_movavg[round(x_window/2):(x_total+round(x_window/2)-1)] # center the moving average
x <- x[1:x_total]
x_std = sd(x_movavg)
x_minus_movavg <- (x - x_movavg) / x_std

# plot first 1000 velocity time steps with moving average
par(mfrow = c(1, 1))
plot(x[1:1000], type = 'l', ylab = "Velocity (m/s)")
lines(x_movavg[1:1000],col='red')
legText = c("Data", "Moving Average")
legend(x = "bottomright", legend = legText, lty = c(1, 1), col = c('black', 'red'))

#plot(x[end_window],type = 'l', ylab = "Velocity (m/s)")
#lines(x_movavg[end_window],col='red')
#legText = c("Data", "Moving Average")
#legend(x = "bottomright", legend = legText, lty = c(1, 1), col = c('black', 'red'))

#plot the whole signal with its moving average
#plot(x,type = 'l', ylab = "Velocity (m/s)")
#lines(x_movavg,col='red')
#legText = c("Data", "Moving Average")
#legend(x = "topright", legend = legText, lty = c(1, 1), col = c('black', 'red'))

#plot the stationary signal
#plot(x_minus_movavg[1:1000],type = 'l', ylab = "", main = "Normalized Difference from Moving Average")
#plot(x_minus_movavg[end_window],type = 'l', ylab = "", main = "Normalized Difference from Moving Average")
#plot(x_minus_movavg,type = 'l', ylab = "", main = "Normalized Difference from Moving Average")

#plot the probablity density function (PDF) of the stationary signal
#x_den <- density(x_minus_movavg)
#plot(x_den, main = "Probablility Density of Normalized Difference from Moving Average")

# calculate langevin coefficients 
bins <- 20 # number of bins to divide the velocities into 
steps <- c(1) # vector of steps to calculate conditional moments for different tau values (change in time). 
# Since we upsample by a factor of 10 in the moving average, I believe we are only interested in a time step of 1.  
ests_detrended <- Langevin1D(x_minus_movavg, bins, steps, sf=fs)
summary(ests_detrended)
plot(ests_detrended)

attach(ests_detrended)

# Determine linear and quadratic models of drift and diffusion components
linearModD1 <- lm(D1 ~ mean_bin, weights = 1/eD1) # linear with no offset
summary(linearModD1)
linearModD2 <- lm(D2 ~ mean_bin + I(mean_bin^2), weights = 1/eD2 ) # quadratic with no linear
summary(linearModD2)

#plot Drift/Diffusion components with their models
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
                         d20 = coefficients(linearModD2)[1], sf=fs)

plot(x_minus_movavg[1:1000], type = 'l', ylab = "", main = "Normalized Difference from Moving Average")
lines(calcVals[1:1000], col='orange')
legText = c("Measured", "Generated")
legend(x = "topright", legend = legText, lty = c(1, 1), col = c('black', 'orange'))

#plot(x_minus_movavg[end_window], type = 'l', ylab = "", main = "Normalized Difference from Moving Average")
#lines(calcVals[end_window], col='red')
#plot(x_minus_movavg, type = 'l', ylab = "", main = "Normalized Difference from Moving Average")
#lines(calcVals, col='red')

x_calc <- x_std * calcVals + x_movavg

# Plot velocities
plot(x[1:1000], type = 'l', ylab = "",)
lines(x_calc[1:1000],col='orange')
legend(x = "bottomright", legend = legText, lty = c(1, 1), col = c('black', 'orange'))
#plot(x[end_window], type = 'l', ylab = "",)
#lines(x_calc[end_window],col='orange')
#legend(x = "bottomright", legend = legText, lty = c(1, 1), col = c('black', 'orange'))
#plot(x, type = 'l', ylab = "",)
#lines(x_calc,col='orange')
#legend(x = "bottomright", legend = legText, lty = c(1, 1), col = c('black', 'orange'))
# this looks pretty good.  

#calculate fast fourier transform
x_fft <- fft(x)
calc_fft <- fft(x_calc)

par(mfrow = c(2, 1))
y.lim = c(0, 300)
plot(Mod(x_fft[1:(x_total/2)]), type = 'l', ylim = y.lim, ylab = "", main = "FFT of Measured Velocities")
plot(Mod(calc_fft[1:(x_total/2)]), type = 'l', ylim = y.lim, ylab = "", xlab = 'Index', main = "FFT of Generated Velocities")

# test if time series is stationary 
adf.test(x_minus_movavg) # p-value < 0.05 indicates the TS is stationary
kpss.test(x_minus_movavg) # p-value > 0.05 indicates TS is stationary 





