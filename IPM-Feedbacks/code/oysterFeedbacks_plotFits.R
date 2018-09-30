## Jacob L Moore (jlmoor@ucdavis.edu)
## 07/01/17

## Plots fits to growth and fecundity data, and age-specific size distributions.
## Creates Figures 1, 2, and 5.

## Requires the following files:
	# growthSurvData_WB.csv
	# WB_sizeDists.Rdata
	# oysterFeedbacks_demography.R
	# oysterFeedbacks_estimateParams.R

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Misc setup
# -------------------------------------------------------------------
# -------------------------------------------------------------------

rm (list = ls ())
graphics.off()
options(warn=-1)

library(bbmle)

setwd('/Users/JLM/Documents/UC_Davis/Schreiber_Lab/OysterProject/IPM_Feedbacks/writeup/code-data/')
data.folder = '/Users/JLM/Documents/UC_Davis/Schreiber_Lab/OysterProject/IPM_Feedbacks/writeup/code-data/data'

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Specify miscellaneous model parameters
# -------------------------------------------------------------------
# -------------------------------------------------------------------
siteName = 'WB'

## Set type of feedback to incorporate. Options include: 
# ‘none’ (for original linear model), or ‘positive’ (for model with positive feedbacks)
feedbacks = 'none'

## Set minimum/maximum size and maximum age
min.size = 0
max.size = 250
max.age = 10
num.calcs = 150
plot.length = 40 	# how many time steps to plot

# Set up discetization for size (this evaluates the integral using the midpoint rule)
n = max.size
xs = seq(min.size, max.size, length=n+1)
deltax = xs[2]-xs[1] # dx increment in sizes
xs.mid.temp = .5 * (xs[-1] + xs[-(n+1)]) # midpoints
xs.mid =  c(xs.mid.temp, xs.mid.temp[n]+deltax)	# add extra point on end for evicted class

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Load data
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# load growth and survival data
d = read.csv(paste(data.folder, '/growthSurvData_WB.csv', sep=''), header=TRUE, sep=',')
d$age = d$timeStep+1

# remove oysters that die for growth fitting
d2 = d[-which(is.na(d$sizeNext)),]

# remove entries where no change to avoid error about non-finite initial values
if(length(which(d2$sizeChange==0)!=0)) {
	d2 = d2[-which(d2$sizeChange==0),]
}


source('oysterFeedbacks_demography.R')
source('oysterFeedbacks_estimateParams.R')

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Plot parameters
# -------------------------------------------------------------------
# -------------------------------------------------------------------

font.name = 'serif'
overall.scale = 1.75
main.scale = 1.25
axis.scale = 1
label.scale = 1.35
line.width.axes = 1.5

line.width = 2.5
line.width.eq = 1
line.width.pts = 3
xx = seq(0,max.size, by=0.05)

age.cols = rainbow(max.age)
age.ltys = c(1:max.age)

legend.text.age9 = c('Age=1', 'Age=2', 'Age=3', 'Age=4', 'Age=5', 'Age=6', 'Age=7', 'Age=8', 'Age=9')
legend.text.age10 = c('Age=1', 'Age=2', 'Age=3', 'Age=4', 'Age=5', 'Age=6', 'Age=7', 'Age=8', 'Age=9', 'Age=10')

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Growth and survival relationships
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# ---------- fits to growth data
quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
plot(d$size[which(d$timeStep==0)], d$log.sizeChange[which(d$timeStep==0)], main='Growth Regression', xlab='Size at time t (mm)', ylab='log(change in size) (mm)', lwd=line.width.pts, col='gray47', axes=T, xlim=c(0, max.size), ylim=c(-1,5))
points(d$size[which(d$timeStep==1)], d$log.sizeChange[which(d$timeStep==1)],lwd=line.width.pts, col='gray70')
sizeNext.fit = params$growth.int + xx*params$growth.size
lines(xx, sizeNext.fit, lwd=line.width, lty=1, col='black')
legend('topright', legend=c('Year 1', 'Year 2'), col=c('gray47', 'gray70'), lty=c(NA, NA), pch=c(1,1), lwd=c(line.width.pts, line.width.pts), cex=.7, horiz=FALSE)
box(lwd=line.width.axes)
dev.copy(pdf, paste(siteName, 'growthRegression1.pdf', sep=''))
dev.off()

# ---------- convert growth fit to sizeNext against sizeCurrent; Figure 1A
quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
plot(d$size[which(d$timeStep==0)], d$sizeNext[which(d$timeStep==0)], main='Growth', xlim=c(0,250), ylim=c(0,250), xlab='Size at time t (mm)', ylab='Size at time t+1 (mm)', lwd=line.width.pts, col='gray47', axes=T)
points(d$size[which(d$timeStep==1)], d$sizeNext[which(d$timeStep==1)],lwd=line.width.pts, col='gray70')
sizeNext.fit = xx + exp(params$growth.int + xx*params$growth.size)
lines(xx, sizeNext.fit, lwd=line.width, lty=1, col='black')
abline(a=0, b=1, lwd=line.width.axes, lty=3)
legend('bottomright', legend=c('Year 1', 'Year 2'), col=c('gray47', 'gray70'), lty=c(NA, NA), pch=c(1,1), lwd=c(line.width.pts, line.width.pts), cex=.7, horiz=FALSE)
box(lwd=line.width.axes)
dev.copy(pdf, paste(siteName, 'growthRegression2.pdf', sep=''))
dev.off()

# ---------- fits to survival data; Figure 1B
quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
plot(d$size[which(d$timeStep==0)], jitter(d$surv[which(d$timeStep==0)], factor=0.1), xlab='Size (mm)', ylab='Probability', main='Survival', xlim=c(0,200), lwd=line.width.pts, col='gray47', axes=T)
points(d$size[which(d$timeStep==1)], jitter(d$surv[which(d$timeStep==1)], factor=0.1), lwd=line.width.pts, col='gray70')
y.surv = predict(surv.reg, data.frame(size=xx), type='response')
lines(xx, y.surv, col='black', lty=1, lwd=line.width)
legend('bottomright', legend=c('Year 1', 'Year 2'), col=c('gray47', 'gray70'), lty=c(NA, NA), pch=c(1,1), lwd=c(line.width.pts, line.width.pts), cex=.7, horiz=FALSE)
box(lwd=line.width.axes)
dev.copy(pdf, paste(siteName, 'survivalRegressions.pdf', sep=''))
dev.off()

# ---------- size distribution of each age' Figure 5
load(paste(data.folder, '/WB_sizeDists.Rdata', sep=''))
quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
matplot(N.mat.P, type='l', lwd=line.width, lty=age.ltys, col=age.cols, xlab='Size', ylab='Frequency')
legend('topright', legend=legend.text.age10, lwd=line.width, lty=age.ltys, col=age.cols, cex=0.5)
box(lwd=line.width.axes)
dev.copy(pdf, paste('WB_ageSpecificSizeDists.pdf', sep=''))
dev.off()


# # -------------------------------------------------------------------
# # -------------------------------------------------------------------
# # Fecundity relationships
# # -------------------------------------------------------------------
# # -------------------------------------------------------------------	

# ---------- Size-dependent number of eggs; Figure 2B
quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
log.num.eggs = params$fec.int + params$fec.slope*log(xx)
plot(log(data.virg.fec$size.fec[may.dates]), log(data.virg.fec$num.eggs[may.dates]), col='gray47',main='Number of Offspring', xlab='log size (mm)', ylab='log number of eggs', axes=T, ylim=c(0, 15), xlim=c(2,5))
lines(log(xx), log.num.eggs, lwd=line.width)
box(lwd=line.width.axes)
dev.copy(pdf, paste(siteName, 'fecundityRegression1.pdf', sep=''))
dev.off()

# ---------- Number of offspring, non log scale
quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
num.eggs = exp(params$fec.int)*xx^params$fec.slope
plot(data.virg.fec$size.fec[may.dates], data.virg.fec$num.eggs[may.dates], col='gray47', xlim=c(0, 250), axes=T, main='Fecundity', xlab='Size (mm)', ylab='Number of eggs', ylim=c(0, max(num.eggs)))
lines(xx, num.eggs, lwd=line.width, col='black')
box(lwd=line.width.axes)
dev.copy(pdf, paste(siteName, 'fecundityRegression2.pdf', sep=''))
dev.off()


# ---------- Size distribution of recruits; Figure 2C
quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
hist(firstYearSizes, main='Recruit Size', xlab='Size (mm)', ylab='Density', freq=FALSE, xlim=c(0,35))
offspring.size = dnorm(xx, mean=params$recruit.size.mean, sd=params$recruit.size.sd)
lines(xx, offspring.size, lwd=line.width, type='l')
box(lwd=line.width.axes)
dev.copy(pdf, paste(siteName, 'recruitDistribution.pdf', sep=''))
dev.off()


# ---------- Size-specific sex ratios; Figure 2A
quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
predict.virg = params$sr.int + params$sr.slope*xx
predict.virg[which(predict.virg>1)]=1
plot(x.sr, sr, main='Sex Ratio', xlab='Size (mm)', ylab='Proportion female', xlim=c(0,250), ylim=c(0,1), pch=20, lwd=line.width.pts, col='gray47')
lines(xx, predict.virg, lwd=line.width)
box(lwd=line.width.axes)
dev.copy(pdf, paste(siteName, 'sexRatio.pdf', sep=''))
dev.off()


