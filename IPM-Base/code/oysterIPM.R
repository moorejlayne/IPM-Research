
## Jacob L Moore (jlmoor@ucdavis.edu); 10/15/15
## Creates and simulates an IPM for C. gigas populations that are recruitment limited and decreasing (lambda~=0.5), stable (lambda~=1), and increasing with high recruitment ((lambda~=1.5).

## INPUT:
# growthsurvivalData.csv
# fecundityData.csv
# sexRatiosWA.csv

## OUTPUT:
# Creates 3 .RData files containing model output (one for lambda=0.5, one for lambda=1.0, and one for lambda=1.5). These files can be plotted using oysterIPM_makePlots.R. 

## SOURCED CODE:
# oysterIPM_demography.R
# oysterIPM_estimateParams.R
# oysterIPM_runIPM.R
# oysterIPM_saveData.R

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
setwd('/Users/JLM/Documents/UC_Davis/Schreiber_Lab/OysterProject/IPM/writeup/EcologicalApplicationsResubmission/code_data/')

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Specify model and define associated parameters
# -------------------------------------------------------------------
# -------------------------------------------------------------------

## Designate minimum/maximum size and maximum age
min.size = 0
max.size = 300
max.age = 15

## Designate type of model used for statistical fitting. Options include:
# ‘size’ (for size-only model) or ‘age.size.int’ (for age- and size-structured model)
statsMethod = 'age.size.int'

## Defines the probability of recruitment. This is the proportion of produced larvae that survive and successfully enter the oyster population. Current values of estab.prob.low, estab.prob.mid, and estab.prob.high are set such that lambda.low~=0.5, lambda.mid~=1, and lambda.high~=1.5, respectively. This assumes the max.size=300, and max.age = 15. If these values are changed, the values listed below will no longer correspond to the appropriate lambda values.  
if(statsMethod=='size') {
	estab.prob.low = 2.45E-18 # lambda=0.96123
	estab.prob.mid = 3.74E-14
	estab.prob.high = 6.68E-12
	max.age = 1
} else {	# statsMethod=='age.size.int'
	estab.prob.low = 2.44E-1500 #2.44E-15
	estab.prob.mid = 1.00E-11
	estab.prob.high = 3.97E-10
}


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Load data
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


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Load demographic functions and fit parameters
# -------------------------------------------------------------------
# -------------------------------------------------------------------

source('oysterIPM_demography.R')
source('oysterIPM_estimateParams.R')


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Build & run IPM, save data
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# set up discetization for size (this evaluates the integral using the midpoint rule)
n = max.size
xs = seq(min.size, max.size, length=n+1)
deltax = xs[2]-xs[1] # dx increment in sizes
xs.mid.temp = .5 * (xs[-1] + xs[-(n+1)]) # midpoints
xs.mid =  c(xs.mid.temp, xs.mid.temp[n]+deltax)	# add extra point on end for evicted class


# Run the IPM for each value of lambda
lamvals = 'low' #c('low','mid', 'high')
for(lam in lamvals) {
	if(lam=='low') {
		params$establishment.prob = estab.prob.low
	} else if(lam=='mid') {
		params$establishment.prob = estab.prob.mid
	} else {
		params$establishment.prob = estab.prob.high
	}
	source('oysterIPM_runIPM.R')
}

