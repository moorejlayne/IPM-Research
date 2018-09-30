## Jacob L Moore (jlmoor@ucdavis.edu); 10/15/15

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Setup and formatting
# -------------------------------------------------------------------
# -------------------------------------------------------------------

rm (list = ls ())
graphics.off()
options(warn=-1)
setwd('/Users/JLM/Documents/UC_Davis/Schreiber_Lab/OysterProject/IPM/writeup/EcologicalApplicationsResubmission/code_data/')

max.size = 300
max.age = 15

# Specify model to plot ('size' or 'age.size.int')
statsMethod = 'age.size.int'

# Do you want to include the evicted class in the plots? (yes=1, no=0)
plot.evict = 0

# ---------- Plot formatting specifications

line.width.pts  = 2.5
line.width.col = 4
line.width.bk = 5
line.width.axes = 3
line.width.contour = 3
line.type.low = 1
line.type.mid = 1
line.type.high = 1
line.color.low = 'black'	
line.color.mid = 'gray50'	# dark gray
line.color.high = 'gray76' 	# light gray

font.name = 'serif'
overall.scale = 1.75
main.scale = 1.5
axis.scale = 1.25
label.scale = 1.5

cols = rainbow(max.age)
line.typs = seq(from=1, to=max.age, by=1)
xx = seq(0,max.size, by=0.05)
plot.power=.8
x.age = line.typs

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Load data for statistical fitting plots
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# load growth and survival data
d = read.csv('growthsurvivalData.csv', header=TRUE, sep=',')
d$age = d$timeStep+1
# remove oysters that die for growth fitting
d2 = d[-which(is.na(d$sizeNext)),]
# remove entries where no change to avoid error about non-finite initial values
if(length(which(d2$sizeChange==0)!=0)) {
	d2 = d2[-which(d2$sizeChange==0),]
}
d2$log.sizeChange = log(d2$sizeChange)
# sample from data to plot less
data.idx = sample(seq(1:length(d2$size)), 1000)
d.sub = d[data.idx,]
source('oysterIPM_estimateParams.R')


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Plots for size-only model
# -------------------------------------------------------------------
# -------------------------------------------------------------------
if(statsMethod=='size') {
	
	filename.data = paste('size_maxage1_maxsize', max.size, sep='')
	filename = filename.data
	load(paste(filename.data, '_lambdamid.RData', sep=''))
	load(paste(filename.data, '_lambdahigh.RData', sep=''))		

# ---------- Reproductive values
	quartz()
	par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
	plot(xs.mid[1:(n+plot.evict)], repVals.size.overall.scaled.mid[1:(n+plot.evict)], type='l', lty=line.type.mid, lwd=line.width.bk, main='Reproductive Values, Scaled \n Size Only', xlab='Size (mm)', ylab='Normalized Reproductive Value', ylim=c(0, max(repVals.size.overall.scaled.high, repVals.size.overall.scaled.mid)), col=line.color.mid)
	lines(xs.mid[1:(n+plot.evict)], repVals.size.overall.scaled.high[1:(n+plot.evict)], type='l', lwd=line.width.bk, lty=line.type.high, col=line.color.high)
	legend('bottomright', legend=c(expression(lambda==1.0), expression(lambda==1.5)), col=c(line.color.mid, line.color.high), lty=c(1,1), lwd=c(2,2), cex=0.7)
	box(lwd=line.width.axes)
	dev.copy(pdf, paste(filename, '_repVals_Scaled.pdf', sep=''))
	dev.off()

# ---------- Stable size
	quartz()
	par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
	plot(xs.mid[1:(n+plot.evict)], stable.size.overall.mid[1:(n+plot.evict)], type='l', lty=line.type.mid, lwd=line.width.bk, main='Stable Size Distribution \n Size Only', xlab='Size (mm)', ylab='Density', ylim=c(0, max(stable.size.overall.high[1:(n+plot.evict)], stable.size.overall.mid[1:(n+plot.evict)])), col=line.color.mid)
	lines(xs.mid[1:(n+plot.evict)], stable.size.overall.high[1:(n+plot.evict)], type='l', lwd=line.width.bk, lty=line.type.high, col=line.color.high)
	legend('topright', legend=c(expression(lambda==1.0), expression(lambda==1.5)), col=c( line.color.mid, line.color.high), lty=c( 1,1), lwd=c( 2,2), cex=0.7)
	box(lwd=line.width.axes)
	dev.copy(pdf, paste(filename, '_stableSize.pdf', sep=''))
	dev.off()
	
# ---------- Survival elasticity
	quartz()
	par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
	plot(xs.mid[1:(n+plot.evict)], elas.P.size.high[1:(n+plot.evict)], type='l', lwd=line.width.bk, main=paste('Size Only Survival Elasticity'), xlab='Size (mm)', ylab='Density', axes=T, ylim=c(0, max(elas.P.size.high[1:(n+plot.evict)], elas.P.size.mid[1:(n+plot.evict)], na.rm=T)), lty=line.type.high, col=line.color.high)	
	lines(xs.mid[1:(n+plot.evict)], elas.P.size.mid[1:(n+plot.evict)], lty=line.type.mid, lwd=line.width.bk, col=line.color.mid)
	legend('topright', legend=c(expression(lambda==1.0), expression(lambda==1.5)), col=c(line.color.mid, line.color.high), lty=c(1,1), lwd=c(2,2), cex=0.7)	
	box(lwd=line.width.axes)
	dev.copy(pdf, paste(filename, '_survElasSize.pdf', sep=''))
	dev.off()
	
# ---------- Fecundity elasticity
	quartz()
	par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
	plot(xs.mid[1:(n+plot.evict)], elas.F.size.high[1:(n+plot.evict)], type='l', lwd=line.width.bk, main=paste('Size Only Fecundity Elasticity'), xlab='Size (mm)', ylab='Density', axes=T, ylim=c(0, max(elas.F.size.high[1:(n+plot.evict)], elas.F.size.mid[1:(n+plot.evict)], na.rm=T)), lty=line.type.high, col=line.color.high)	
	lines(xs.mid[1:(n+plot.evict)], elas.F.size.mid[1:(n+plot.evict)], lty=line.type.mid, lwd=line.width.bk, col=line.color.mid)
	legend('topleft', legend=c(expression(lambda==1.0), expression(lambda==1.5)), col=c(line.color.mid, line.color.high), lty=c(1,1), lwd=c(2,2), cex=0.7)	
	box(lwd=line.width.axes)
	dev.copy(pdf, paste(filename, '_fecElasSize.pdf', sep=''))
	dev.off()

# ---------- Statistical Fitting	
# # plot legend separately
quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
plot(0,0, col='white', axes=FALSE, xlab='', ylab='')
legend.txt = as.character(seq(from=1, to=max.age, by=1))
legend.txt = c('Year 1', 'Year 2', 'Year 3')
legend('center', legend=legend.txt, col=c('gray0', 'gray47', 'gray70'), lty=c(NA, NA, NA), pch=c(1,1,1), lwd=c(line.width.pts, line.width.pts, line.width.pts), cex=.7, horiz=FALSE)
dev.copy(pdf, 'sizeLegend.pdf')
dev.off()

# fits to growth data
quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
plot(d.sub$size[which(d.sub$timeStep==0)], log(d.sub$sizeChange[which(d.sub$timeStep==0)]), main='Growth Regression', xlab='Size at time t (mm)', ylab='log(change in size) (mm)', lwd=line.width.pts, col='gray0', axes=T, xlim=c(6, 200), ylim=c(0,5))
points(d.sub$size[which(d.sub$timeStep==1)], log(d.sub$sizeChange[which(d.sub$timeStep==1)]),lwd=line.width.pts, col='gray47')
points(d.sub$size[which(d.sub$timeStep==2)], log(d.sub$sizeChange[which(d.sub$timeStep==2)]),lwd=line.width.pts, col='gray70')
for(currentAge in 1) {
	sizeNext.fit = params$growth.int + xx*params$growth.size + currentAge*params$growth.age + xx*currentAge*params$growth.age.size
	lines(xx, sizeNext.fit, lwd=line.width.col, lty=line.typs[currentAge], col=cols[currentAge])
}
box(lwd=line.width.axes)
dev.copy(pdf, paste(filename, 'growthRegression1.pdf', sep=''))
dev.off()

# convert growth fit to sizeNext against sizeCurrent
quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
plot(d.sub$size[which(d.sub$timeStep==0)], d.sub$sizeNext[which(d.sub$timeStep==0)], main='Growth', xlim=c(6,200), ylim=c(0,250), xlab='Size at time t (mm)', ylab='Size at time t+1 (mm)', lwd=line.width.pts, col='gray0', axes=T)
points(d.sub$size[which(d.sub$timeStep==1)], d.sub$sizeNext[which(d.sub$timeStep==1)],lwd=line.width.pts, col='gray47')
points(d.sub$size[which(d.sub$timeStep==2)], d.sub$sizeNext[which(d.sub$timeStep==2)],lwd=line.width.pts, col='gray70')
for(currentAge in 1) {
	sizeNext.fit = xx + exp(params$growth.int + xx*params$growth.size + currentAge*params$growth.age + xx*currentAge*params$growth.age.size)
	lines(xx, sizeNext.fit, lwd=line.width.col, lty=line.typs[currentAge], col=cols[currentAge])
}
abline(a=0, b=1, lwd=line.width.axes, lty=3)
box(lwd=line.width.axes)
dev.copy(pdf, paste(filename, 'growthRegression2.pdf', sep=''))
dev.off()

# fits to survival data
quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
plot(d.sub$size[which(d.sub$timeStep==0)], jitter(d.sub$surv[which(d.sub$timeStep==0)], factor=0.1), xlab='Size (mm)', ylab='Probability', main='Survival', xlim=c(6,200), lwd=line.width.pts, col='black', axes=T)
points(d.sub$size[which(d.sub$timeStep==1)], jitter(d.sub$surv[which(d.sub$timeStep==1)], factor=0.1), lwd=line.width.pts, col='gray47')
points(d.sub$size[which(d.sub$timeStep==2)], jitter(d.sub$surv[which(d.sub$timeStep==2)], factor=0.1), lwd=line.width.pts, col='gray70')
for(currentAge in 1) {
	y.surv = predict(surv.reg, data.frame(size=xx, age=currentAge), type='response')
	lines(xx, y.surv, col=cols[currentAge], lty=line.typs[currentAge], lwd=line.width.col)
}
box(lwd=line.width.axes)
dev.copy(pdf, paste(filename, 'survivalRegressions.pdf', sep=''))
dev.off()


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Plots for age- and size-structured model
# -------------------------------------------------------------------
# -------------------------------------------------------------------	
} else {
	filename.data = paste('agesize_maxage', max.age, '_maxsize', max.size, sep='')
	filename = filename.data
	load(paste(filename.data, '_lambdalow.RData', sep=''))
	load(paste(filename.data, '_lambdamid.RData', sep=''))
	load(paste(filename.data, '_lambdahigh.RData', sep=''))

# ---------- Reproductive values
	# Size specific
	quartz()
	par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
	plot(xs.mid[1:(n+plot.evict)], repVals.size.overall.scaled.low[1:(n+plot.evict)], type='l', lty=line.type.low, lwd=line.width.bk, main=paste('Size-specific Reproductive Values\n Age x Size', sep=''), xlab='Size (mm)', ylab='Normalized Reproductive Value', ylim=c(0, max(repVals.size.overall.scaled.high[1:(n+plot.evict)], repVals.size.overall.scaled.mid[1:(n+plot.evict)], na.rm=T)), col=line.color.low)
	lines(xs.mid[1:(n+plot.evict)], repVals.size.overall.scaled.high[1:(n+plot.evict)], type='l', lwd=line.width.bk, lty=line.type.mid, col=line.color.mid)
	lines(xs.mid[1:(n+plot.evict)], repVals.size.overall.scaled.mid[1:(n+plot.evict)], type='l', lwd=line.width.bk, lty=line.type.high, col=line.color.high)
	legend('topleft', legend=c(expression(lambda==0.5), expression(lambda==1.0), expression(lambda==1.5)), col=c(line.color.low, line.color.mid, line.color.high), lty=c(1,1,1), lwd=c(2,2,2), cex=0.7)
	box(lwd=line.width.axes)
	dev.copy(pdf, paste(filename, '_repValsSize_Scaled.pdf', sep=''))
	dev.off()
	
	# Age specific
	quartz()
	par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
	plot(x.age,repVals.age.scaled.low, main=paste('Age-specific Reproductive Values, \n Age x Size', sep=''), ylab='Normalized Reproductive Value', xlab='Age (years)', ylim=c(0, max(repVals.age.scaled.low, repVals.age.scaled.mid, repVals.age.scaled.high, na.rm=T)), type='b', lty=1, lwd=line.width.bk, col=line.color.low, pch=19)
	points(x.age,repVals.age.scaled.high, type='b', lty=line.type.high, lwd=line.width.bk, col=line.color.high, pch=19)
	points(x.age,repVals.age.scaled.mid, type='b', lty=line.type.mid, lwd=line.width.bk, col=line.color.mid, pch=19)
	legend('topright', legend=c(expression(lambda==0.5), expression(lambda==1.0), expression(lambda==1.5)), col=c(line.color.low, line.color.mid, line.color.high), pch=c(19,19,19), lty=c(1,1,1), lwd=c(2,2,2), cex=0.7)
	box(lwd=line.width.axes)
	dev.copy(pdf, paste(filename, '_repValsAge_Scaled.pdf', sep=''))
	dev.off()
	
	
# ---------- Stable size
	# Size specific
	quartz()
	par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
	plot(xs.mid[1:(n+plot.evict)], stable.size.overall.low[1:(n+plot.evict)], type='l', lty=line.type.low, lwd=line.width.bk, main=paste('Stable Size Distribution \n Age x Size', sep=''), xlab='Size (mm)', ylab='Density', ylim=c(0, max(stable.size.overall.high[1:(n+plot.evict)], stable.size.overall.low[1:(n+plot.evict)], na.rm=T)), col=line.color.low)
	lines(xs.mid[1:(n+plot.evict)], stable.size.overall.mid[1:(n+plot.evict)], type='l', lwd=line.width.bk, lty=line.type.mid, col=line.color.mid)
	lines(xs.mid[1:(n+plot.evict)], stable.size.overall.high[1:(n+plot.evict)], type='l', lwd=line.width.bk, lty=line.type.high, col=line.color.high)
	legend('topright', legend=c(expression(lambda==0.5), expression(lambda==1.0), expression(lambda==1.5)), col=c(line.color.low, line.color.mid, line.color.high), lty=c(1,1,1), lwd=c(2,2,2), cex=0.7)
	box(lwd=line.width.axes)
	dev.copy(pdf, paste(filename, '_stableSize.pdf', sep=''))
	dev.off()
	
	# Age specific
	quartz()
	par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
	plot(x.age,stable.age.low, main=paste('Stable Age Distribution \n Age x Size', sep=''), ylab='Density', xlab='Age (years)', ylim=c(0, max(stable.age.low, stable.age.mid, stable.age.high, na.rm=T)), type='b', lty=1, lwd=line.width.bk, col=line.color.low, pch=19)
	points(x.age,stable.age.high, type='b', lty=line.type.high, lwd=line.width.bk, col=line.color.high, pch=19)
	points(x.age,stable.age.mid, type='b', lty=line.type.mid, lwd=line.width.bk, col=line.color.mid, pch=19)
	legend('top', legend=c(expression(lambda==0.5), expression(lambda==1.0), expression(lambda==1.5)), col=c(line.color.low, line.color.mid, line.color.high), pch=c(19,19,19), lty=c(1,1,1), lwd=c(2,2,2), cex=0.7)
	box(lwd=line.width.axes)
	dev.copy(pdf, paste(filename, '_stableAge.pdf', sep=''))
	dev.off()


# ---------- Survival elasticity
	# Size specific
	quartz()
	par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
	plot(xs.mid[1:(n+plot.evict)], elas.P.size.high[1:(n+plot.evict)], type='l', lwd=line.width.bk, main=paste('Size-specific Survival Elasticity \n Age x Size', sep=''), xlab='Size (mm)', ylab='Density', axes=T, ylim=c(0, max(elas.F.size.low[1:(n+plot.evict)], elas.F.size.high[1:(n+plot.evict)], elas.F.size.mid[1:(n+plot.evict)], elas.P.size.low[1:(n+plot.evict)], elas.P.size.high[1:(n+plot.evict)], elas.P.size.mid[1:(n+plot.evict)], na.rm=T)), lty=line.type.high, col=line.color.high)	
	lines(xs.mid[1:(n+plot.evict)], elas.P.size.low[1:(n+plot.evict)], lty=line.type.low, lwd=line.width.bk, col=line.color.low)
	lines(xs.mid[1:(n+plot.evict)], elas.P.size.mid[1:(n+plot.evict)], lty=line.type.mid, lwd=line.width.bk, col=line.color.mid)
	legend('topleft', legend=c(expression(lambda==0.5), expression(lambda==1.0), expression(lambda==1.5)), col=c(line.color.low, line.color.mid, line.color.high), lty=c(1,1,1), lwd=c(2,2,2), cex=0.7)
	box(lwd=line.width.axes)
	dev.copy(pdf, paste(filename, '_survElasSize.pdf', sep=''))
	dev.off()
	
	# Age specific
	quartz()
	par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
	plot(x.age,elas.P.age.low, main=paste('Age-specific Survival Elasticity \n Age x Size', sep=''), xlab='Age (years)', ylab='Density', axes=T, ylim=c(0,max(elas.F.age.low, elas.F.age.high, elas.F.age.mid, elas.P.age.low, elas.P.age.high, elas.P.age.mid, na.rm=T)), lwd=line.width.bk, lty=line.type.low, col=line.color.low, type='b', pch=19)
	points(x.age, elas.P.age.mid, type='b', lty=line.type.mid, lwd=line.width.bk, col=line.color.mid, pch=19)
	points(x.age, elas.P.age.high, type='b', lty=line.type.high, lwd=line.width.bk, col=line.color.high, pch=19)
	legend('bottomleft', legend=c(expression(lambda==0.5), expression(lambda==1.0), expression(lambda==1.5)), col=c(line.color.low, line.color.mid, line.color.high), pch=c(19,19,19), lty=c(1,1,1), lwd=c(2,2,2), cex=0.7)
	box(lwd=line.width.axes)
	dev.copy(pdf, paste(filename, '_survElasAge.pdf', sep=''))
	dev.off()

# ---------- Fecundity elasticity
	# Size specific
	quartz()
	par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
	plot(xs.mid[1:(n+plot.evict)], elas.F.size.low[1:(n+plot.evict)], type='l', lwd=line.width.bk, main=paste('Size-specific Fecundity Elasticity \n Age x Size', sep=''), xlab='Size (mm)', ylab='Density', axes=T, ylim=c(0, max(elas.F.size.low[1:(n+plot.evict)], elas.F.size.high[1:(n+plot.evict)], elas.F.size.mid[1:(n+plot.evict)], elas.P.size.low[1:(n+plot.evict)], elas.P.size.high[1:(n+plot.evict)], elas.P.size.mid[1:(n+plot.evict)], na.rm=T)), lty=line.type.low, col=line.color.low)	
	lines(xs.mid[1:(n+plot.evict)], elas.F.size.mid[1:(n+plot.evict)], lty=line.type.mid, lwd=line.width.bk, col=line.color.mid)
	lines(xs.mid[1:(n+plot.evict)], elas.F.size.high[1:(n+plot.evict)], lty=line.type.high, lwd=line.width.bk, col=line.color.high)
	legend('topleft', legend=c(expression(lambda==0.5), expression(lambda==1.0), expression(lambda==1.5)), col=c(line.color.low, line.color.mid, line.color.high), lty=c(1,1,1), lwd=c(2,2,2), cex=0.7)
	box(lwd=line.width.axes)
	dev.copy(pdf, paste(filename, '_fecElasSize.pdf', sep=''))
	dev.off()

	# Age specific
	quartz()	
	par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
	plot(x.age,elas.F.age.low, main=paste('Age-specific Fecundity Elasticity \n Age x Size', sep=''), xlab='Age (years)', ylab='Density', axes=T, ylim=c(0,max(elas.F.age.low, elas.F.age.high, elas.F.age.mid, elas.P.age.low, elas.P.age.high, elas.P.age.mid, na.rm=T)), lwd=line.width.bk, lty=line.type.low, col=line.color.low, type='b', pch=19)
	points(x.age, elas.F.age.mid, type='b', lty=line.type.mid, lwd=line.width.bk, col=line.color.mid, pch=19)
	points(x.age, elas.F.age.high, type='b', lty=line.type.high, lwd=line.width.bk, col=line.color.high, pch=19)
	legend('topleft', legend=c(expression(lambda==0.5), expression(lambda==1.0), expression(lambda==1.5)), col=c(line.color.low, line.color.mid, line.color.high), pch=c(19,19,19), lty=c(1,1,1), lwd=c(2,2,2), cex=0.7)
	box(lwd=line.width.axes)	
	dev.copy(pdf, paste(filename, '_fecElasAge.pdf', sep=''))
	dev.off()	


# ---------- Statistical Fitting
# # plot legend separately
quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
plot(0,0, col='white', axes=FALSE, xlab='', ylab='')
legend.txt = as.character(seq(from=1, to=max.age, by=1))
legend.txt = c('Year 1', 'Year 2', 'Year 3', '', 'Age=1', 'Age=2', 'Age=3', 'Age=4', 'Age=5', 'Age=6', 'Age=7', 'Age=8', 'Age=9', 'Age=10')
legend('center', legend=legend.txt, col=c('gray0', 'gray47', 'gray70', 'white', cols), lty=c(NA, NA, NA, NA, line.typs), pch=c(1,1,1,rep(NA, max.age+1)), lwd=c(line.width.pts, line.width.pts, line.width.pts, 2, rep(line.width.col, max.age)), cex=.7, horiz=FALSE)
dev.copy(pdf, 'agesizeLegend.pdf')
dev.off()

# fits to growth data
quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
plot(d.sub$size[which(d.sub$timeStep==0)], log(d.sub$sizeChange[which(d.sub$timeStep==0)]), main='Growth Regression', xlab='Size at time t (mm)', ylab='log(change in size) (mm)', lwd=line.width.pts, col='gray0', axes=T, xlim=c(6, 200), ylim=c(0,5))
points(d.sub$size[which(d.sub$timeStep==1)], log(d.sub$sizeChange[which(d.sub$timeStep==1)]),lwd=line.width.pts, col='gray47')
points(d.sub$size[which(d.sub$timeStep==2)], log(d.sub$sizeChange[which(d.sub$timeStep==2)]),lwd=line.width.pts, col='gray70')
for(currentAge in 1:max.age) {
	sizeNext.fit = params$growth.int + xx*params$growth.size + currentAge*params$growth.age + xx*currentAge*params$growth.age.size
	lines(xx, sizeNext.fit, lwd=line.width.col, lty=line.typs[currentAge], col=cols[currentAge])
}
box(lwd=line.width.axes)
dev.copy(pdf, paste(filename, 'growthRegression1.pdf', sep=''))
dev.off()

# convert growth fit to sizeNext against sizeCurrent
quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
plot(d.sub$size[which(d.sub$timeStep==0)], d.sub$sizeNext[which(d.sub$timeStep==0)], main='Growth', xlim=c(6,200), ylim=c(0,250), xlab='Size at time t (mm)', ylab='Size at time t+1 (mm)', lwd=line.width.pts, col='gray0', axes=T)
points(d.sub$size[which(d.sub$timeStep==1)], d.sub$sizeNext[which(d.sub$timeStep==1)],lwd=line.width.pts, col='gray47')
points(d.sub$size[which(d.sub$timeStep==2)], d.sub$sizeNext[which(d.sub$timeStep==2)],lwd=line.width.pts, col='gray70')
for(currentAge in 1:max.age) {
	sizeNext.fit = xx + exp(params$growth.int + xx*params$growth.size + currentAge*params$growth.age + xx*currentAge*params$growth.age.size)
	lines(xx, sizeNext.fit, lwd=line.width.col, lty=line.typs[currentAge], col=cols[currentAge])
}
abline(a=0, b=1, lwd=line.width.axes, lty=3)
box(lwd=line.width.axes)
dev.copy(pdf, paste(filename, 'growthRegression2.pdf', sep=''))
dev.off()

# fits to survival data
quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
plot(d.sub$size[which(d.sub$timeStep==0)], jitter(d.sub$surv[which(d.sub$timeStep==0)], factor=0.1), xlab='Size (mm)', ylab='Probability', main='Survival', xlim=c(6,200), lwd=line.width.pts, col='black', axes=T)
points(d.sub$size[which(d.sub$timeStep==1)], jitter(d.sub$surv[which(d.sub$timeStep==1)], factor=0.1), lwd=line.width.pts, col='gray47')
points(d.sub$size[which(d.sub$timeStep==2)], jitter(d.sub$surv[which(d.sub$timeStep==2)], factor=0.1), lwd=line.width.pts, col='gray70')
for(currentAge in 1:max.age) {
	y.surv = predict(surv.reg, data.frame(size=xx, age=currentAge), type='response')
	lines(xx, y.surv, col=cols[currentAge], lty=line.typs[currentAge], lwd=line.width.col)
}
box(lwd=line.width.axes)
dev.copy(pdf, paste(filename, 'survivalRegressions.pdf', sep=''))
dev.off()

}

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Fecundity-related plots (will plot regardless of statsMethod)
# -------------------------------------------------------------------
# -------------------------------------------------------------------	

# ---------- Size-dependent number of eggs
quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
log.num.eggs = params$fec.int + params$fec.slope*xx
offspring.size = dnorm(xx, mean=params$recruit.size.mean, sd=params$recruit.size.sd)
plot(xx, log.num.eggs, main='Number of Offspring', xlab='Size (mm)', ylab='log number of eggs', lwd=line.width.bk, type='l', xlim=c(0, 300), axes=T)
points(data.gig.fec$size.fec, log(data.gig.fec$num.eggs))
box(lwd=line.width.axes)
dev.copy(pdf, 'fecundityRegression.pdf')
dev.off()

# ---------- Size distribution of recruits
quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
hist(firstYearSizes, main='Offspring Size', xlab='Size (mm)', ylab='Density', freq=FALSE, ylim=c(0,0.07), xlim=c(0,65))
lines(xx, offspring.size, main='Offspring Size', lwd=line.width.bk, type='l')
box(lwd=line.width.axes)
dev.copy(pdf, 'recruitDistribution.pdf')
dev.off()

# ---------- Size-specific sex ratios
quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
# From Buroker 1983 (Gigas from WA)
data = read.csv(paste('sexRatiosWA.csv', sep=''), header=TRUE, sep=',')
sizes = rep(data$shellLength*10, 4)
propFemale = c(data$propF1, data$propF2, data$propF3, data$propF4)
m.sr.gigas = lm(propFemale~sizes)
#(Intercept) 0.0310907
#sizes       0.0043945
size.seq = seq(from=0, to=230, by=0.5)
predict.gigas = m.sr.gigas$coefficients[1] + m.sr.gigas$coefficients[2]*size.seq
predict.gigas[which(predict.gigas>1)]=1
plot(sizes, propFemale, main='Sex Ratio', xlab='Size', ylab='Proportion female', xlim=c(0,230), ylim=c(0,1), pch=20)
lines(size.seq, predict.gigas, lwd=line.width.bk)
box(lwd=line.width.axes)
dev.copy(pdf, 'sexRatio.pdf')
dev.off()


