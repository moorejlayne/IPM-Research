## Jacob L Moore (jlmoor@ucdavis.edu)
## 07/01/17

## Plots phase plane and trajectories for three initial conditions (starting from either equilibrium size distribution or harvested size distribution). Includes separatrix for each starting size distribution, as well as figure of each starting size distribution. 
## Creates Figure 3.

## Requires the following files:
	# growthSurvData_WB.csv
	# WB_eqVals_Linear.Rdata
	# oysterFeedbacks_demography.R
	# oysterFeedbacks_estimateParams.R
	# oysterFeedbacks_runIPM.R
	# WB_separatrixOutput_final.csv

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Misc setup
# -------------------------------------------------------------------
# -------------------------------------------------------------------

rm (list = ls ())
graphics.off()
options(warn=-1)

## load libraries
library(TeachingDemos)

setwd('/Users/JLM/Documents/UC_Davis/Schreiber_Lab/OysterProject/IPM_Feedbacks/writeup/code-data/')
data.folder = '/Users/JLM/Documents/UC_Davis/Schreiber_Lab/OysterProject/IPM_Feedbacks/writeup/code-data/data'

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Specify miscellaneous model parameters
# -------------------------------------------------------------------
# -------------------------------------------------------------------

feedbacks = 'positive'
rhoVal = 'low'
deltaVal = 'mid'
alphaVal = 'estimate'

## Designate minimum/maximum size and maximum age
min.size = 0
max.size = 250
max.age = 10
num.calcs = 150
siteName = 'WB'
plot.length = 50

# set up discetization for size (this evaluates the integral using the midpoint rule)
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

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Load demographic functions, and fit parameters
# -------------------------------------------------------------------
# -------------------------------------------------------------------

source('oysterFeedbacks_demography.R')
source('oysterFeedbacks_estimateParams.R')

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Specify restoration parameters
# -------------------------------------------------------------------
# -------------------------------------------------------------------

shell.add.val = 0
oyster.add.val = 0

# i=1 represents initial conditions, so to add in 'first' time step, need to add at t=2
time.add.shell = 2
time.add.oyster = 2
age.to.add = 1

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Run model with multiple initial conditions
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# these are run in pairs (e.g. for a given simulation, H=initH[i] and N=initN[i])
initHs = c(50000, 30000, 55000)
initNs = c(5E6, 9E6, 1.2E7)

# set initial size distribution
dists = c('distEq', 'distHarvest')	
dist.mat = eq.Mat[,,1] / sum(eq.Mat[,,1])	
dist.mat.Eq = dist.mat / sum(dist.mat)
dist.mat.Harv = dist.mat
dist.mat.Harv[which(xs.mid.temp>=75)[1]:dim(dist.mat.Harv)[1],] = dist.mat.Harv[which(xs.mid.temp>=75)[1]:dim(dist.mat.Harv)[1],] * 0
dist.mat.Harv = dist.mat.Harv / sum(dist.mat.Harv)

# build array to store trajectories for each model run
H.sum.all = array(NA, dim=c(length(dists)*length(initHs), num.calcs))
N.sum.all = H.sum.all

count = 1

# run model for each initial condition
for(dist.type in dists) {
			
	if(dist.type=='distEq') {
		dist.mat = dist.mat.Eq
	} else {
		dist.mat = dist.mat.Harv
	} 							
	shell.addition = shell.add.val
	total.oyster.add = oyster.add.val
	oyster.addition = array(0, dim=c(max.size+1, max.age, length(time.add.oyster)))
	oyster.addition[,age.to.add,1]=total.oyster.add*(stable.dist.matrix.linear[,age.to.add] / sum(stable.dist.matrix.linear[,age.to.add]	))	
	for(ss in 1:length(initHs)) {
		initN = initNs[ss]
		init.vect = initN*dist.mat
		initH = initHs[ss]
		init.vect.dead = initH
		source('oysterFeedbacks_runIPM.R')
		H.sum.all[count,] = H.sum
		N.sum.all[count,] = N.sum	
		count = count+1
	}
}

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Plot formatting
# -------------------------------------------------------------------
# -------------------------------------------------------------------

font.name = 'serif'
overall.scale = 2
main.scale = 1.4
axis.scale = 1.2
label.scale = 1.5
line.width.axes = 1
legend.scale = 0.9
legend.scale.mod = 0.5

line.width.eq = 2
line.type.eq = 1
line.col.eq = 'black'

line.width.harv = 3
line.type.harv = 6
line.col.harv = 'gray50'

pch.eq = 1 #19
pch.start = 8
pch.end = 15

pch.points.cex = 1.3

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Plots
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# load separatrix data
d.additions = read.table(paste(data.folder, '/WB_separatrixOutput_final.csv', sep=''), header=TRUE, sep=',')
oystersAdded = cbind(d.additions$oystersAdded[which(d.additions$popDist=='distEq')], d.additions$oystersAdded[which(d.additions$popDist=='distHarvest')])
initialH = cbind(d.additions$initH[which(d.additions$popDist=='distEq')], d.additions$initH[which(d.additions$popDist=='distHarvest')])

H.sum.eq = t(H.sum.all[1:3,])
H.sum.harv = t(H.sum.all[4:6,])
N.sum.eq = t(N.sum.all[1:3,])
N.sum.harv = t(N.sum.all[4:6,])

y.range = c(2E6, 1.4E7)
x.range = c(min(initialH)+10, max(initialH))


filename = paste('WB_positive_rholow_alphaestimate_deltamid', sep='')

# Figure 3
pdf(paste(filename, '_phasePlots.pdf', sep=''), width=10, height=7)
par(bg='white', family=font.name, cex = overall.scale, cex.main=main.scale, cex.lab=label.scale, cex.axis=axis.scale, mar=c(4.5,4.5,4,3))
layout(rbind(c(1,1,2), c(1,1,3)))

# plot phase plane
matplot(initialH[,1], oystersAdded, type='l', lty=c(1,4), lwd=1, ylab='Total oysters', col='gray50', xlim = x.range, ylim=y.range, xlab= expression("Substrate ("*""*m^2*")"), main='Phase Plane')
points(H.eqs, N.eqs, pch=c(pch.eq), col=c('black'), lwd=3, cex= pch.points.cex+0.5)
matplot(H.sum.eq[1:plot.length,], N.sum.eq[1:plot.length,], type='l', lwd=line.width.eq, col=line.col.eq, add=TRUE, lty=line.type.eq)
matplot(H.sum.harv[1:plot.length,], N.sum.harv[1:plot.length,], type='l', lwd=line.width.harv, col=line.col.harv, add=TRUE, lty=line.type.harv)
points(initHs, initNs, pch=pch.start, cex= pch.points.cex, lwd=2)
legend('topleft', legend=c('Unstable equilibrium', 'Initial condition', 'Initial population: equilibrium distribution', '    Simulation Trajectory', '    Separatrix', 'Initial population: harvested distribution', '    Simulation Trajectory', '    Separatrix'), col=c('black', 'black', 'white', line.col.eq, line.col.eq, 'white', line.col.harv, line.col.harv), pch=c(pch.eq, pch.start, NA, NA, NA, NA, NA, NA), lty=c(NA, NA, NA, line.type.eq, line.type.eq, NA, line.type.harv, line.type.harv), lwd=c(1, 1, 1, line.width.eq, 1, 1, line.width.harv, 1 ), cex=legend.scale, inset=c(0,0))
legend('topright', 'I', bty='n', inset=c(0,0.15), cex= legend.scale+legend.scale.mod)
legend('topright', 'II', bty='n', inset=c(0,0.55), cex=legend.scale+legend.scale.mod)
legend('bottomright', 'III', bty='n', inset=c(0, 0.1), cex=legend.scale+legend.scale.mod )

# plot equilibrium distribution
plot(xs.mid, rowSums(dist.mat.Eq), type='l', lty=line.type.eq, lwd=line.width.eq, col=line.col.eq, main='Equilibrium Distribution', xlab='Size', ylab='Frequency')

# plot harvested distribution
plot(xs.mid, rowSums(dist.mat.Harv), type='l', lty=line.type.harv, lwd=line.width.harv, col=line.col.harv, main='Harvested Distribution', xlab='Size', ylab='Frequency')
dev.off()
	