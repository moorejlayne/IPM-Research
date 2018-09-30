## Jacob L Moore (jlmoor@ucdavis.edu)
## 07/01/17

## Projects forward the survival+growth kernel to determine the 
## size distribution of oysters of a particular age.

## Requires the following files:
	# growthSurvData_WB.csv
	# WB_eqVals_Linear.Rdata
	# oysterFeedbacks_demography.R
	# oysterFeedbacks_estimateParams.R
## Creates the following files:
	# WB_sizeDists.Rdata


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

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Load demographic functions and fit parameters
# -------------------------------------------------------------------
# -------------------------------------------------------------------

source('oysterFeedbacks_demography.R')
source('oysterFeedbacks_estimateParams.R')

load(paste(data.folder, '/WB_eqVals_Linear.Rdata', sep=''))

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Iterate growth-surv kernel
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# set initial population
init.vect = matrix(0, nrow=(n+1), ncol=max.age)
init.vect[,1]=20000*dnorm(xs.mid, mean=params$recruit.size.mean, sd=params$recruit.size.sd)

# set up kernels for growth and survival
S = array(0, dim=c((n+1),max.age))
G = array(0, dim=c((n+1),(n+1),max.age))
P = array(0, dim=c((n+1),(n+1),max.age))

stable.dens = array(NA, dim=c((max.size+1), max.age))
for(a in 1:max.age) {
	# build growth kernel
	G.temp = deltax*outer(xs.mid.temp,xs.mid.temp,g.yx, currentAge=a, params=params)	
	prob.growth.evict = 1-colSums(G.temp)
	G.temp = rbind(G.temp, prob.growth.evict)
	colAdd = c(rep(0,times=n), 1)
	G.temp = cbind(G.temp, colAdd)	
	G[,,a] = G.temp
	
	# build survival and shell kernel
	S.temp = s.x(xs.mid.temp,currentAge=a,params=params, max.age=max.age)
	S.temp = c(S.temp, S.temp[n])
	S[,a] = S.temp

	# build survival*growth kernel
	P.temp = G.temp; for(i in 1:n) P.temp[,i]=G.temp[,i]*S[i,a]
	P[,,a]=P.temp		 	
}
K = P
N.mat = array(NA, dim=c((n+1),max.age))
N.mat[,1] = init.vect[1:(n+1)]
for(a in 2:max.age) {
	N.mat[,a] = K[,,(a-1)] %*% N.mat[,(a-1)]
}
	
# rescale so each age sums to 1
for(i in 1:max.age) {
	N.mat[,i] = N.mat[,i] / sum(N.mat[,i])
	stable.dens[,i] = stable.dist.matrix.linear[,i]/sum(stable.dist.matrix.linear[,i])
}
N.mat.P = N.mat

save(N.mat.P, stable.dens, file=paste(siteName, '_sizeDists.Rdata', sep=''))



