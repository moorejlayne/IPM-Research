## Jacob L Moore (jlmoor@ucdavis.edu)
## 07/01/17

## Implements bisectional search to determine total amount of substrate required to cross the threshold surface.

## Requires the following files:
	# growthSurvData_WB.csv
	# WB_eqVals_Linear.Rdata
	# oysterFeedbacks_demography.R
	# oysterFeedbacks_estimateParams.R
	# oysterFeedbacks_runIPM.R

## Creates the following files:
	# WB_substrateOutput_final_<timestamp>.csv  [Note, '_<timestamp>' must be removed for this file to work with oysterFeedbacks_plotAdditions.R]


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

## Iterate through desired parameter combinations
rhoVals = c('low', 'estimate') 			# 'low' OR 'estimate'
alphaVals = c('estimate') 				# 'low' OR 'estimate' OR 'high'
deltaVals = c('mid', 'low', 'high')		# 'low' OR 'mid' OR 'high'    

# additional specifications
initN.mod.seq = c(0.9, 1.1)					# percent of N.eq at which to begin population
dist.start.seq = c('distHarvest', 'distEq')	# size-distribution of initial population

# specify temporary and final filenames to save
filename.temp.substrate = paste(siteName, '_substrateOutput_temp_', tt, sep='')
filename.final.substrate = paste(siteName, '_substrateOutput_final_', tt, sep='')

# simulate substrate addition for all parameter combinations
idx.count = 1
for(RR in 1:length(rhoVals)) {
	rhoVal=rhoVals[RR]
	for(AA in 1:length(alphaVals)) {
		alphaVal=alphaVals[AA] 
		for(DD in 1:length(deltaVals)) {
			deltaVal=deltaVals[DD]
			
			source('oysterFeedbacks_estimateParams.R')

			for(DS in 1:length(dist.start.seq)) {
				# set up initial population size distribution
				dist.type=dist.start.seq[DS]
				dist.mat = eq.Mat[,,1] / sum(eq.Mat[,,1])
				if(dist.type=='distEq') {
					dist.mat = dist.mat
					dist.mat = dist.mat / sum(dist.mat)
				} else { # (dist.type=='distHarvest') 
					dist.mat[which(xs.mid.temp>=75)[1]:dim(dist.mat)[1],] = dist.mat[which(xs.mid.temp>=75)[1]:dim(dist.mat)[1],] * 0
					dist.mat = dist.mat / sum(dist.mat)
				} 
				
				# iterate through each initial population size
				for(ii in 1:length(initN.mod.seq)) {
					# display current parameter combo
					print(paste('delta=', deltaVal, ', alpha=', alphaVal, ', rho=', rhoVal, ', startDist=', dist.type, ', initN=', initN.mod.seq[ii], sep=''))				
					# set up population
					initN = initN.mod.seq[ii]*N.eqs
					init.vect = initN*dist.mat 
					oyster.addition = array(0, dim=c(max.size+1, max.age, length(time.add.oyster)))
					initH = 0
					init.vect.dead = initH
					
					# initialize bisection search
					shell.addition = H.eqs*1.1
					last.shell.addition = shell.addition					
					lastOver = shell.addition
					lastUnder = 0
					tolCond = lastOver-lastUnder
					negPop = TRUE					
					count = 1					
					
					# initialize temporary file for troubleshooting
					if(idx.count==1) {
						write.table(t(c(rhoVal, alphaVal, deltaVal, N.eqs, H.eqs, initN, dist.type, round(last.shell.addition, digits=2), round(tolCond, digits=2), round(lastUnder, digits=2), round(lastOver, digits=2), negPop, count)), file=paste(filename.temp.substrate, '.csv', sep=''), append=FALSE, sep=',', col.names=c('rhoVal', 'alphaVal', 'deltaVal', 'Neqs', 'Heqs', 'initN', 'startDist', 'substrateAdded', 'tolCond',  'lastUnder', 'lastOver', 'negPop', 'count'), row.names=FALSE)
					} else {
						write.table(t(c(rhoVal, alphaVal, deltaVal, N.eqs, H.eqs, initN, dist.type, round(last.shell.addition, digits=2), round(tolCond, digits=2), round(lastUnder, digits=2), round(lastOver, digits=2), negPop, count)), file=paste(filename.temp.substrate, '.csv', sep=''), append=TRUE, sep=',', col.names=FALSE, row.names=FALSE)
					}				
						
					# while not close enough to equilibrium, and population is going toward extinction				
					while(tolCond>H.eqs*0.001 | negPop) {
						# while not close enough to equilibrium, and population is going toward extinction
						print(paste('tolCond=', round(tolCond, digits=2), ', substrate added=', round(last.shell.addition, digits=2), ', negPop=', negPop, sep=''))
						
						# run model
						source('oysterFeedbacks_runIPM.R')
						last.shell.addition = shell.addition
						
						# check if over threshold (i.e. population persisting) or not (i.e. population declining to zero)
						if(N.sum[num.calcs]>N.eqs && N.sum[num.calcs]>N.sum[num.calcs-10] && H.sum[num.calcs]>H.eqs && H.sum[num.calcs]>H.sum[num.calcs-10]) {
							# if overthreshold
							lastOver = shell.addition
							shell.addition = (shell.addition + lastUnder) / 2
							negPop=FALSE	
							count=0			
						} else {
							# if under threshold
							lastUnder = shell.addition
							negPop=TRUE					
							if(count==1) { # didn't get high enough on first try
								shell.addition = shell.addition*5
								lastOver = shell.addition
							} else {
								shell.addition = (shell.addition + lastOver) / 2 
								count=0
							}
						}
						tolCond = lastOver-lastUnder	
						
						# update temporary .csv file
						write.table(t(c(rhoVal, alphaVal, deltaVal, N.eqs, H.eqs, initN, dist.type, round(last.shell.addition, digits=2), round(tolCond, digits=2), round(lastUnder, digits=2), round(lastOver, digits=2), negPop, count)), file=paste(filename.temp.substrate, '.csv', sep=''), append=TRUE, sep=',', col.names=FALSE, row.names=FALSE)		
					}
					
					# update final .csv file
					if(idx.count==1) {
						write.table(t(c(rhoVal, alphaVal, deltaVal, N.eqs, H.eqs, initN, dist.type, round(last.shell.addition, digits=2), round(tolCond, digits=2), num.calcs)), file=paste(filename.final.substrate, '.csv', sep=''), append=FALSE, sep=',', col.names=c('rhoVal', 'alphaVal', 'deltaVal', 'Neqs', 'Heqs', 'initN', 'startDist', 'substrateAdded', 'tolCond', 'num.calcs'), row.names=FALSE)	
					} else {
						write.table(t(c(rhoVal, alphaVal, deltaVal, N.eqs, H.eqs, initN, dist.type, round(last.shell.addition, digits=2), round(tolCond, digits=2), num.calcs)), file=paste(filename.final.substrate, '.csv', sep=''), append=TRUE, sep=',', col.names=FALSE, row.names=FALSE)
					}
					idx.count=idx.count+1
				}				
			}
		}
	}
}		