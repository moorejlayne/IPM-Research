
## Jacob L Moore (jlmoor@ucdavis.edu)
## 07/01/17

## Runs a single iteration of the IPM model with either no feedbacks, or positive feedbacks. 

## Requires the following files:
	# growthSurvData_WB.csv
	# oysterFeedbacks_demography.R
	# oysterFeedbacks_estimateParams.R
	# oysterFeedbacks_runIPM.R
	# oysterFeedbacks_makePlots.R
	
## Creates the following files:
	# WB_eqVals_Linear.Rdata (if run w/feedbacks='none')
	# Plots of initial size distributions, populations trajectories, and phase planes

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Misc setup
# -------------------------------------------------------------------
# -------------------------------------------------------------------

rm (list = ls ())
graphics.off()
options(warn=-1)

## load libraries
library(bbmle)
library(TeachingDemos)

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
feedbacks = 'positive'

## Set whether to used estimated value of local retention, or a low value of local retention. This is only relevant when feedbacks=='positive'.
rhoVal = 'low'		# 'low' OR 'estimate'

## Set whether to used estimated value of alpha, or a low or high value. This is only relevant when feedbacks=='positive'.
alphaVal = 'estimate' # 'estimate' OR 'low' OR 'high'

## Set whether to used high, mid, or low value for delta. This is only relevant when feedbacks=='positive'.
deltaVal = 'mid' # 'low' OR 'mid' OR 'high'


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

# # load growth and survival data
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

# Distribution of existing oyster population. 
dist.type = 'distHarvest' # 'distEq' OR 'distHarvest' OR 'distScorched'

eqMod.N = -0.05		# what percent away from N equilibrium to start 
eqMod.H = -0.05		# what percent away from H equilibrium to start
	# the above only relevant if feedbacks=='positive'

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Build & run IPM, save data
# -------------------------------------------------------------------
# -------------------------------------------------------------------

if(feedbacks=='none') {
	init.vect = matrix(0, nrow=(n+1), ncol=max.age)
	init.vect[,1]=dnorm(xs.mid, mean=params$recruit.size.mean, sd=params$recruit.size.sd)
	init.vect.dead = 0
	shell.addition = 0
	oyster.addition = array(0, dim=c(max.size+1, max.age, length(time.add.oyster)))
	
	source('oysterFeedbacks_runIPM.R')
	filename = paste('WB_Linear', sep='')
		
} else {	
	dist.mat = eq.Mat[,,1] / sum(eq.Mat[,,1])		
	if(dist.type=='distEq') {
		dist.mat = dist.mat / sum(dist.mat)
	} else if(dist.type=='distHarvest') {
		dist.mat[which(xs.mid.temp>=75)[1]:dim(dist.mat)[1],] = dist.mat[which(xs.mid.temp>=75)[1]:dim(dist.mat)[1],] * 0
		dist.mat = dist.mat / sum(dist.mat)
	} else { # dist.type=='distScorched'
		dist.mat = dist.mat*0
	}
				
	shell.addition = shell.add.val
	total.oyster.add = oyster.add.val
	oyster.addition = array(0, dim=c(max.size+1, max.age, length(time.add.oyster)))
	oyster.addition[,age.to.add,1]=total.oyster.add*(stable.dist.matrix.linear[,age.to.add] / sum(stable.dist.matrix.linear[,age.to.add]))

	initN = N.eqs + N.eqs*eqMod.N
	initH = H.eqs + H.eqs*eqMod.H
	init.vect = initN*dist.mat 
	init.vect.dead = initH
	source('oysterFeedbacks_runIPM.R')
	
	filename = paste(siteName, '_', feedbacks, '_rho', rhoVal, '_alpha', alphaVal, '_delta', deltaVal, '_initDist', dist.type, '_eqModN', eqMod.N, '_eqModH', eqMod.H, sep='')	
		
}

source('oysterFeedbacks_makePlots.R')
