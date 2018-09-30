## Jacob L Moore (jlmoor@ucdavis.edu); 10/15/15
## Using growth kernels from the age- and size-structured model, simulates individual trajectories of growth. This generates output given in Appendix A. 

## INPUT: 
# growthsurvivalData.csv

## OUTPUT:
# Creates two plots showing individual size at each time point for 1000 simulations; one plot for normal growth fitting (allowing for positive and negative growth), and one plot for lognormal growth fitting (allowing for only positive growth). 

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Misc setup
# -------------------------------------------------------------------
# -------------------------------------------------------------------

rm (list = ls ())
graphics.off()
options(warn=-1)

library(bbmle)
setwd('/Users/JLM/Documents/UC_Davis/Schreiber_Lab/OysterProject/IPM/writeup/Ecological Applications Resubmission/code_data/')

optimMethod='Nelder-Mead'
growthMethod='age.size.int'

min.size = 0
max.size = 300
max.age = 15
num.reps = 5000
time.steps=max.age

filename = paste('agesize_maxage', max.age, '_maxsize', max.size, sep='')

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Load data
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# load full data
d = read.csv('growthsurvivalData.csv', header=TRUE, sep=',')
d$age = d$timeStep+1
# remove oysters that die for growth fitting
d2 = d[-which(is.na(d$sizeNext)),]
# remove entries where no change to avoid error about non-finite initial values
if(length(which(d2$sizeChange==0)!=0)) {
	d2 = d2[-which(d2$sizeChange==0),]
}
d2$log.sizeChange = log(d2$sizeChange)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Fit models
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# distribution of offspring sizes	
firstYearSizes = d$size[which(d$timeStep==d$timeStep[1])]	# pull out data from first time point
recruit.size.mean=mean(firstYearSizes, na.rm=TRUE)
recruit.size.sd=sd(firstYearSizes, na.rm=TRUE)

start.size=rnorm(num.reps, recruit.size.mean, recruit.size.sd) 

# growth
growth.fitting = function(growthMethod, fitMethod, d2) {
	optimMethod = 'Nelder-Mead'
	# allows for positive and negative growth
	if(fitMethod=='normal') {
		m = mle2(sizeChange~dnorm(mean=a+b*age+d*size+h*size*age, sd=c), start=list(a=mean(d2$sizeChange), b=0, d=0, h=0, c=sd(d2$sizeChange)), data=d2, method=optimMethod)
	# only allows for positive growth
	} else { # 'lognormal'
		m = mle2(log.sizeChange~dnorm(mean=a+b*age+d*size+h*size*age, sd=c), start=list(a=mean(d2$log.sizeChange), b=0, d=0, h=0, c=sd(d2$log.sizeChange)), data=d2, method=optimMethod)		
	}
}

# survival
surv.fitting = function(growthMethod, d) {
	if(growthMethod =='size') {
		m = glm(surv~size, data=d, family=binomial)
	} else { 	# 'age.size.int'
		m = glm(surv~size+age+size*age, data=d, family=binomial)
	}
	return(m)
}
surv.reg = surv.fitting(growthMethod, d)

# Fit normal model
params=data.frame(
  	surv.int=0,
  	surv.age=0,
  	surv.size=0,
  	surv.age.size=0,
   	growth.int=0,
  	growth.var=0,
  	growth.age=0,
  	growth.size=0,
  	growth.age.size=0
)


growth.reg = growth.fitting(growthMethod, 'normal', d2)
growth.reg.normal = growth.reg
if(growthMethod=='age.size.int') {
	params$growth.int = coef(growth.reg)[1]				# a
	params$growth.age = coef(growth.reg)[2]				# b
	params$growth.size = coef(growth.reg)[3]			# d
	params$growth.age.size = coef(growth.reg)[4]			# h
	params$growth.var = coef(growth.reg)[5]				# c
} else {		# 'size'
	params$growth.int = coef(growth.reg)[1]				# a
	params$growth.size = coef(growth.reg)[2]			# b	
	params$growth.var = coef(growth.reg)[3]				# c		
}
if(growthMethod =='age.size.int') {
	params$surv.int = coefficients(surv.reg)[1]
	params$surv.size = coefficients(surv.reg)[2]	
	params$surv.age = coefficients(surv.reg)[3]	
	params$surv.age.size = coefficients(surv.reg)[4]			
} else {
	params$surv.int = coefficients(surv.reg)[1]
	params$surv.size = coefficients(surv.reg)[2]
}
params.normal = params

# Fit lognormal model
params=data.frame(
  	surv.int=0,
  	surv.age=0,
  	surv.size=0,
  	surv.age.size=0,
   	growth.int=0,
  	growth.var=0,
  	growth.age=0,
  	growth.size=0,
  	growth.age.size=0
)

growth.reg = growth.fitting(growthMethod, 'lognormal', d2)
growth.reg.lognormal = growth.reg
if(growthMethod=='age.size.int') {
	params$growth.int = coef(growth.reg)[1]				# a
	params$growth.age = coef(growth.reg)[2]				# b
	params$growth.size = coef(growth.reg)[3]			# d
	params$growth.age.size = coef(growth.reg)[4]			# h
	params$growth.var = coef(growth.reg)[5]				# c
} else {		# 'size'
	params$growth.int = coef(growth.reg)[1]				# a
	params$growth.size = coef(growth.reg)[2]			# b	
	params$growth.var = coef(growth.reg)[3]				# c		
}
if(growthMethod =='age.size.int') {
	params$surv.int = coefficients(surv.reg)[1]
	params$surv.size = coefficients(surv.reg)[2]	
	params$surv.age = coefficients(surv.reg)[3]	
	params$surv.age.size = coefficients(surv.reg)[4]			
} else {
	params$surv.int = coefficients(surv.reg)[1]
	params$surv.size = coefficients(surv.reg)[2]
}
params.lognormal = params


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Simulate growth iterations
# -------------------------------------------------------------------
# -------------------------------------------------------------------

simulate.growth = function(start.size, time.steps, fitMethod, params) {
	indiv.size = rep(0,time.steps)
	indiv.size[1] = start.size
	for(tstep in 2:time.steps) {
		currentSize = indiv.size[tstep-1]
		currentAge = tstep-1
		# grow
		if(currentSize!=0) {
			if(currentSize>=max.size) {
				nextSize=max.size
			} else {
				if(fitMethod=='normal') {
					nextSize = currentSize + rnorm(1, mean=params$growth.int + currentAge*params$growth.age + currentSize*params$growth.size + currentAge*currentSize*params$growth.age.size, sd=params$growth.var) 
				} else { # fitMethod=='lognormal'
					nextSize = currentSize + rlnorm(1, mean=params$growth.int + currentAge*params$growth.age + currentSize*params$growth.size + currentAge*currentSize*params$growth.age.size, sd=params$growth.var) 
				}				
				if(nextSize>max.size) nextSize=max.size
			}
		} else { nextSize=0 }		
		indiv.size[tstep]=nextSize
	}
	return(indiv.size)
}

indiv.growth.normal = sapply(start.size, simulate.growth, time.steps=time.steps, fitMethod='normal', params=params.normal)
indiv.growth.lognormal = sapply(start.size, simulate.growth, time.steps=time.steps, fitMethod='lognormal', params=params.lognormal)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Simulate growth & survival iterations
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# survival function (probability of surviving)
	s.x = function(currentSize, currentAge, params) {
		u = exp(params$surv.int + currentAge*params$surv.age + currentSize*params$surv.size + currentSize*currentAge*params$surv.age.size)
		return(u/(1+u))
	}

simulate.growth.surv = function(start.size, time.steps, fitMethod, params) {
	indiv.size = rep(0,time.steps)
	indiv.size[1] = start.size
	for(tstep in 2:time.steps) {
		currentSize = indiv.size[tstep-1]
		currentAge = tstep-1
			
		# survive
		surv.prob = s.x(currentSize,currentAge,params)
		dice = runif(1, min=0, max=1)
		if(dice>surv.prob) {currentSize=0} else {currentSize=currentSize}
		
		# grow
		if(currentSize!=0) {
			if(currentSize>=max.size) {
				nextSize=max.size
			} else {
				if(fitMethod=='normal') {
					nextSize = currentSize + rnorm(1, mean=params$growth.int + currentAge*params$growth.age + currentSize*params$growth.size + currentAge*currentSize*params$growth.age.size, sd=params$growth.var) 
				} else { # fitMethod=='lognormal'
					nextSize = currentSize + rlnorm(1, mean=params$growth.int + currentAge*params$growth.age + currentSize*params$growth.size + currentAge*currentSize*params$growth.age.size, sd=params$growth.var) 
				}				
				if(nextSize>max.size) nextSize=max.size
			}
		} else { nextSize=0 }
				
		indiv.size[tstep]=nextSize
	}
	return(indiv.size)
}
indiv.growth.surv.normal = sapply(start.size, simulate.growth.surv, time.steps=time.steps, fitMethod='normal', params=params.normal)
indiv.growth.surv.lognormal = sapply(start.size, simulate.growth.surv, time.steps=time.steps, fitMethod='lognormal', params=params.lognormal)



# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Make plots, no mortality
# -------------------------------------------------------------------
# -------------------------------------------------------------------

font.name = 'serif'
overall.scale = 1.75
main.scale = 1.5
axis.scale = 1.25
label.scale = 1.5
line.width.axes = 3
line.width.average = 5
line.width.sims = 0.5

# normal growth
quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
tvect = seq(from=1, to=time.steps, by=1)
avg.growth.normal = rowMeans(indiv.growth.normal)
matplot(tvect,indiv.growth.normal, type='l', lwd=line.width.sims, xlab='Time (years)', ylab='Size (mm)', ylim=c(0,max.size), main='Normal Growth Fitting', col='darkgray')
lines(tvect,avg.growth.normal,col='black',lwd=line.width.average)
box(lwd=line.width.axes)
dev.copy(pdf, paste(filename, '_normalTrajectories.pdf', sep=''))
dev.off()

# lognormal growth
quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
tvect = seq(from=1, to=time.steps, by=1)
avg.growth.lognormal = rowMeans(indiv.growth.lognormal)
matplot(tvect,indiv.growth.lognormal, type='l', lwd=line.width.sims, xlab='Time (years)', ylab='Size (mm)', ylim=c(0,max.size), main='Lognormal Growth Fitting', col='darkgray')
lines(tvect,avg.growth.lognormal,col='black',lwd=line.width.average)
box(lwd=line.width.axes)
dev.copy(pdf, paste(filename, '_lognormalTrajectories.pdf', sep=''))
dev.off()

# # # create histogram of final time point of growth trajectories
# x.seq = seq(from=min(indiv.growth.lognormal[max.age,]), to=max(indiv.growth.lognormal[max.age,]), by=0.1)
# y.seq=dnorm(x.seq, mean=mean(indiv.growth.lognormal[max.age,]), sd = sd(indiv.growth.lognormal[max.age,]))
# hist(indiv.growth.lognormal[max.age,],freq=F, breaks=50,xlab='size (mm)', ylab='density', main=paste('reps=', num.reps, ' age=', max.age))
# lines(x.seq,y.seq)

# overlay stable size distribution
load("/Users/JLM/Documents/UC_Davis/Schreiber_Lab/OysterProject/IPM/writeup/EcologicalApplicationsResubmission/code_data/output/agesize_maxage15_maxsize300_lambdalow.RData")

quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
hist(indiv.growth.lognormal[max.age,],freq=F, breaks=60,xlab='size (mm)', ylab='density', main=paste('reps=', num.reps, ' age=', max.age), xlim = c(150,275))
lines(xs.mid, stable.size.overall.low, lwd=2)
#box(lwd=line.width.axes)
legend('topright', legend='stable size, lambda=0.5', lwd=2, col='black', cex=0.4)
dev.copy(pdf, paste(filename, '_stableSize_Trajectory_overlay.pdf', sep=''))
dev.off()


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Make plots, mortality
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# lognormal growth
quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
tvect = seq(from=1, to=time.steps, by=1)
avg.growth.surv.lognormal = rowMeans(indiv.growth.surv.lognormal)
# remove those that die
indiv.growth.surv.lognormal.nozero = indiv.growth.surv.lognormal
indiv.growth.surv.lognormal.nozero[which(indiv.growth.surv.lognormal.nozero==0)]=NA
avg.growth.surv.lognormal.nozero = rowMeans(indiv.growth.surv.lognormal.nozero, na.rm=TRUE)

#percent.surv = 1-(length(which(indiv.growth.surv.lognormal[15,]==0)) / num.reps)	# 0.375
percent.surv = 1-(length(which(indiv.growth.surv.lognormal[max.age,]==0)) / num.reps)	# 0.375

matplot(tvect,indiv.growth.surv.lognormal, type='l', lwd=line.width.sims, xlab='Time (years)', ylab='Size (mm)', ylim=c(0,max.size), main='Lognormal Growth Fitting \nw/Mortality', col='darkgray')
lines(tvect,avg.growth.surv.lognormal.nozero,col='black',lwd=line.width.average)
lines(tvect,avg.growth.surv.lognormal,col='black',lwd=2, lty=2)
box(lwd=line.width.axes)
dev.copy(pdf, paste(filename, '_lognormalTrajectories_mortality.pdf', sep=''))
dev.off()

# overlay stable size distribution
load("/Users/JLM/Documents/UC_Davis/Schreiber_Lab/OysterProject/IPM/writeup/EcologicalApplicationsResubmission/code_data/output/agesize_maxage15_maxsize300_lambdalow.RData")

quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
hist(indiv.growth.surv.lognormal.nozero[max.age,],freq=F, breaks=60,xlab='size (mm)', ylab='density', main=paste('reps=', num.reps, ' age=', max.age), xlim = c(160,275))
lines(xs.mid, stable.size.overall.low, lwd=2)
#box(lwd=line.width.axes)
legend('topright', legend='stable size, lambda=0.5', lwd=2, col='black', cex=0.4)
dev.copy(pdf, paste(filename, '_stableSize_Trajectory_overlay_mortality.pdf', sep=''))
dev.off()

	
quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
tvect = seq(from=1, to=time.steps, by=1)
matplot(tvect,indiv.growth.surv.lognormal, type='l', lwd=line.width.sims, xlab='Time (years)', ylab='Size (mm)', ylim=c(0,max.size), main=paste('Lognormal Growth Fitting \nw/Mortality (', 100*percent.surv, '% surv)', sep=''), col='darkgray')
lines(tvect,avg.growth.surv.lognormal.nozero,col='black',lwd=line.width.average)
lines(tvect,avg.growth.surv.lognormal,col='black',lwd=2, lty=2)
lines(tvect,avg.growth.lognormal,col='red',lwd=line.width.average, lty=2)
legend('topleft', legend=c('growth+surv (living)', 'growth', 'growth+surv (all)'), lty=c(1,2,2), lwd=c(line.width.average, line.width.average,2), col=c('black', 'red', 'black'), cex=0.6)
box(lwd=line.width.axes)
dev.copy(pdf, paste(filename, '_lognormalTrajectories_mortalityComp.pdf', sep=''))
dev.off()

