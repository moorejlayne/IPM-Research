
## Jacob L Moore (jlmoor@ucdavis.edu)
## 07/01/17

## Implements bisectional search to determine total number of oysters of a particular age that are required to cross the threshold surface, beginning from a scorched earth. 

## Requires the following files:
	# growthSurvData_WB.csv
	# WB_eqVals_Linear.Rdata
	# oysterFeedbacks_demography.R
	# oysterFeedbacks_estimateParams.R
	# oysterFeedbacks_runIPM.R

## Creates the following files:
	# WB_oysterOutput_final_<timestamp>.csv  [Note, '_<timestamp>' must be removed for this file to work with oysterFeedbacks_plotAdditions.R]

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

# age cohorts to add, and size distribution of each age
age.to.add.seq = c(1:9)	
size.dist = N.mat.P

# initial substrate amount
initH = 0

# Specify temporary and final filenames to save
filename.temp.oyster = paste(siteName, '_oysterOutput_temp_', tt, sep='')		# for troubleshooting
filename.final.oyster = paste(siteName, '_oysterOutput_final_', tt, sep='')

# simulate oyster addition for all parameter combinations
idx.count = 1
for(RR in 1:length(rhoVals)) {
	rhoVal=rhoVals[RR]
	for(AA in 1:length(alphaVals)) {
		alphaVal=alphaVals[AA] 
		for(DD in 1:length(deltaVals)) {
			deltaVal=deltaVals[DD]
			
			source('oysterFeedbacks_estimateParams.R')
			
			# set so no oysters to begin with
			initN = 0	
			dist.mat = eq.Mat[,,1] / sum(eq.Mat[,,1])
			init.vect = initN*dist.mat

			init.vect.dead = initH
			
			# iterate through each age 	
			for(aa in 1:length(age.to.add.seq)) {
				age.to.add = age.to.add.seq[aa]
				
				# display current parameter combo 
				print(paste('delta=', deltaVal, ', alpha=', alphaVal, ', rho=', rhoVal, ', age=', aa, sep=''))
				
				# initialize bisection search 															
				if(aa==1) {		# if on first age, start above the equilibrium
					total.oyster.add = N.eqs*1.5	
					last.oyster.added=total.oyster.add
				} else {			
					# if already found required amount for another age, start there (this helps sim run faster)
					total.oyster.add = last.oyster.added
				}						
				shell.addition=0
				last.shell.addition=0	
				lastOver = total.oyster.add
				lastUnder = 0
				tolCond = lastOver-lastUnder
				negPop = TRUE
		
				count=1
						
				# initialize temporary file for troubleshooting
				if(idx.count==1) {
					write.table(t(c(rhoVal, alphaVal, deltaVal, N.eqs, H.eqs, initH, round(last.oyster.added, digits=2), aa, round(tolCond, digits=2), round(lastUnder, digits=2), round(lastOver, digits=2), negPop, count)), file=paste(filename.temp.oyster, '.csv', sep=''), append=FALSE, sep=',', col.names=c('rhoVal', 'alphaVal', 'deltaVal', 'Neqs', 'Heqs', 'initH', 'oystersAdded', 'a', 'tolCond',  'lastUnder', 'lastOver', 'negPop', 'count'), row.names=FALSE)
				} else {
					write.table(t(c(rhoVal, alphaVal, deltaVal, N.eqs, H.eqs, initH, round(last.oyster.added, digits=2), aa, round(tolCond, digits=2), round(lastUnder, digits=2), round(lastOver, digits=2), negPop, count)), file=paste(filename.temp.oyster, '.csv', sep=''), append=TRUE, sep=',', col.names=FALSE, row.names=FALSE)
				}
						
				# while not close enough to equilibrium, and population is going toward extinction
				while(tolCond > N.eqs*0.001 | negPop) {
					# display current parameter combo 
					print(paste('tolCond=', round(tolCond, digits=2), ', oysters added=', round(last.oyster.added, digits=2), ', negPop=', negPop, sep=''))
					
					# set number of oysters added							
					oyster.addition = array(0, dim=c(max.size+1, max.age, length(time.add.oyster)))
					oyster.addition[,age.to.add,1] = total.oyster.add*(size.dist[,age.to.add])
			
					# run model
					source('oysterFeedbacks_runIPM.R')
					last.oyster.added = total.oyster.add	
							
					# check if over threshold (i.e. population persisting) or not (i.e. population declining to zero)
					if(N.sum[num.calcs]>N.eqs && N.sum[num.calcs]>N.sum[num.calcs-10] && H.sum[num.calcs]>H.eqs && H.sum[num.calcs]>H.sum[num.calcs-10]) {
						# if over threshold
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
					write.table(t(c(rhoVal, alphaVal, deltaVal, N.eqs, H.eqs, initH, round(last.oyster.added, digits=2), aa, round(tolCond, digits=2), round(lastUnder, digits=2), round(lastOver, digits=2), negPop, count)), file=paste(filename.temp.oyster, '.csv', sep=''), append=TRUE, sep=',', col.names=FALSE, row.names=FALSE)
				}

				# update final .csv file
				if(idx.count==1) {
					write.table(t(c(rhoVal, alphaVal, deltaVal, N.eqs, H.eqs, initH, last.oyster.added, aa,  round(tolCond, digits=2), num.calcs)), file=paste(filename.final.oyster, '.csv', sep=''), append=FALSE, sep=',', col.names=c('rhoVal', 'alphaVal', 'deltaVal', 'Neqs', 'Heqs', 'initH', 'oystersAdded', 'a','tolCond', 'num.calcs'), row.names=FALSE)	
				} else {
					write.table(t(c(rhoVal, alphaVal, deltaVal, N.eqs, H.eqs, initH, last.oyster.added, aa, round(tolCond, digits=2), num.calcs)), file=paste(filename.final.oyster, '.csv', sep=''), append=TRUE, sep=',', col.names=FALSE, row.names=FALSE)
				}
				idx.count=idx.count+1
			}								
		}
	}
}