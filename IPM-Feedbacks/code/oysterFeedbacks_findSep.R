## Jacob L Moore (jlmoor@ucdavis.edu)
## 07/01/17

## Runs the bisectional search to find the separatrix between persistence and extinction for different initial population size distributions.

## Requires the following files:
	# growthSurvData_WB.csv
	# WB_sizeDists.Rdata
	# WB_eqVals_Linear.Rdata
	# oysterFeedbacks_demography.R
	# oysterFeedbacks_estimateParams.R
	# oysterFeedbacks_runIPM.R
	
## Creates the following files:
	# WB_separatrixOutput_final_<timestamp>.csv  [Note, '_<timestamp>' must be removed for this file to work with oysterFeedbacks_plotSep.R]

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

setwd('/Users/JLM/Documents/UC_Davis/Schreiber_Lab/OysterProject/IPM_Feedbacks/writeup/code-data/')
data.folder = '/Users/JLM/Documents/UC_Davis/Schreiber_Lab/OysterProject/IPM_Feedbacks/writeup/code-data/data'

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Specify miscellaneous model parameters
# -------------------------------------------------------------------
# -------------------------------------------------------------------
siteName = 'WB'
feedbacks = 'positive'

## Set minimum/maximum size and maximum age
min.size = 0
max.size = 250
max.age = 10
num.calcs = 150

# Set up discetization for size (this evaluates the integral using the midpoint rule)
n = max.size
xs = seq(min.size, max.size, length=n+1)
deltax = xs[2]-xs[1] # dx increment in sizes
xs.mid.temp = .5 * (xs[-1] + xs[-(n+1)]) # midpoints
xs.mid =  c(xs.mid.temp, xs.mid.temp[n]+deltax)	# add extra point on end for evicted class

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Load data and source files
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

load(paste(data.folder, '/WB_sizeDists.Rdata', sep=''))
source('oysterFeedbacks_demography.R')

# i=1 represents initial conditions, so to add in 'first' time step, need to add at t=2
time.add.shell = c(2)
time.add.oyster = c(2)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Set up conditions through which to iterate
# -------------------------------------------------------------------
# -------------------------------------------------------------------

tt = substr(Sys.time(), 12,19)
tt = chartr(':', '.', tt)

rhoVal = 'low'				# 'low' OR 'estimate'
alphaVal = 'estimate'		# 'low' OR 'estimate' OR 'high'
deltaVal = 'mid'				# 'low' OR 'mid' OR 'high'
popDists = c('distHarvest', 'distEq')		# starting population size distribution; 'distEq' OR 'distHarvest'

source('oysterFeedbacks_estimateParams.R')

# Because the simulation takes a while to run, the following lines allow you to 
# specify the span of initial substrate values on which to run the search. 
# The values below are interesting for alpha='estimate', delta='mid', 
# and popDist='distHarvest'
if(rhoVal=='low') {
	initH.seq = seq(60000, 70000, by=1000)
} else {
	initH.seq = seq(10,610,by=25)
}

# specify filenames
filename.temp.oyster = paste('WB_separatrixOutput_temp_', tt, sep='')
filename.final.oyster = paste('WB_separatrixOutput_final_', tt, sep='')

idx.count = 1

for(popDist in popDists) {
	

	# set so no oysters to begin
	initN = 0
	dist.mat = eq.Mat[,,1] / sum(eq.Mat[,,1])
	init.vect = initN*dist.mat
	
	# set up initial population size distribution
	if(popDist=='distEq') {
		size.dist = dist.mat / sum(dist.mat)
	} else {
		size.dist = dist.mat
		size.dist[which(xs.mid.temp>=75)[1]:dim(size.dist)[1],] = size.dist[which(xs.mid.temp>=75)[1]:dim(size.dist)[1],] * 0
		size.dist = size.dist / sum(size.dist)
	}
	
	# for each initial value of substrate, find the minimum number of oysters (distributed across sizes according to popDist) such that the population will persist
	for (initH in initH.seq) {
		init.vect.dead = 0 #initH
		print(paste('sub=', initH, sep=''))
		
		# start search with oysters above the equilibrium
		total.oyster.add = N.eqs*1.75
		last.oyster.added = total.oyster.add
		shell.addition=initH
		last.shell.addition=0
	
		lastOver = total.oyster.add
		lastUnder = 0
		tolCond = lastOver-lastUnder
		negPop = TRUE
	
		count=1	
		
		# initialize temporary file for troubleshooting
		if(idx.count==1) {
			write.table(t(c(rhoVal, alphaVal, deltaVal, N.eqs, H.eqs, initH, round(last.oyster.added, digits=2), popDist, round(tolCond, digits=2), round(lastUnder, digits=2), round(lastOver, digits=2), negPop, count)), file=paste(filename.temp.oyster, '.csv', sep=''), append=FALSE, sep=',', col.names=c('rhoVal', 'alphaVal', 'deltaVal', 'Neqs', 'Heqs', 'initH', 'oystersAdded', 'popDist', 'tolCond',  'lastUnder', 'lastOver', 'negPop', 'count'), row.names=FALSE)
		} else {
			write.table(t(c(rhoVal, alphaVal, deltaVal, N.eqs, H.eqs, initH, round(last.oyster.added, digits=2), popDist, round(tolCond, digits=2), round(lastUnder, digits=2), round(lastOver, digits=2), negPop, count)), file=paste(filename.temp.oyster, '.csv', sep=''), append=TRUE, sep=',', col.names=FALSE, row.names=FALSE)
		}	
		
		# while not close enough to equilibrium, and population is going toward extinction
		while(tolCond > N.eqs*0.00001 | negPop) {
			print(paste('tolCond=', round(tolCond, digits=2), ', oysters added=', round(last.oyster.added, digits=2), ', negPop=', negPop, sep=''))
								
			# set number of oysters added							
			oyster.addition = array(0, dim=c(max.size+1, max.age, length(time.add.oyster)))
			oyster.addition[,,1] = total.oyster.add*(size.dist)
	
			# run model
			source('oysterFeedbacks_runIPM.R')
			last.oyster.added = total.oyster.add	
	
			# check if over threshold (i.e. population persisting) or not (i.e. population declining to zero)
			if(N.sum[num.calcs]>N.eqs && N.sum[num.calcs]>N.sum[num.calcs-10] && H.sum[num.calcs]>H.eqs && H.sum[num.calcs]>H.sum[num.calcs-10]) {
				# if overthreshold
				lastOver = total.oyster.add
				total.oyster.add = (total.oyster.add + lastUnder) / 2
				negPop=FALSE	
				count=0			
			} else {
				# if under threshold
				lastUnder = total.oyster.add
				negPop=TRUE					
				if(count==1) { # didn't get high enough on first try
					total.oyster.add = total.oyster.add*3
					lastOver = total.oyster.add
				} else {
					total.oyster.add = (total.oyster.add + lastOver) / 2 
					count=0
				}
	
			}
			tolCond = lastOver-lastUnder
			
			# update temporary .csv file
			write.table(t(c(rhoVal, alphaVal, deltaVal, N.eqs, H.eqs, initH, round(last.oyster.added, digits=2), popDist, round(tolCond, digits=2), round(lastUnder, digits=2), round(lastOver, digits=2), negPop, count)), file=paste(filename.temp.oyster, '.csv', sep=''), append=TRUE, sep=',', col.names=FALSE, row.names=FALSE)
		}	
	
		# update final .csv file
		if(idx.count==1) {
			write.table(t(c(rhoVal, alphaVal, deltaVal, N.eqs, H.eqs, initH, last.oyster.added, popDist, round(tolCond, digits=2), num.calcs)), file=paste(filename.final.oyster, '.csv', sep=''), append=FALSE, sep=',', col.names=c('rhoVal', 'alphaVal', 'deltaVal', 'Neqs', 'Heqs', 'initH', 'oystersAdded', 'popDist', 'tolCond', 'num.calcs'), row.names=FALSE)	
		} else {
			write.table(t(c(rhoVal, alphaVal, deltaVal, N.eqs, H.eqs, initH, last.oyster.added, popDist, round(tolCond, digits=2), num.calcs)), file=paste(filename.final.oyster, '.csv', sep=''), append=TRUE, sep=',', col.names=FALSE, row.names=FALSE)
		}
		idx.count=idx.count+1		
	}
}