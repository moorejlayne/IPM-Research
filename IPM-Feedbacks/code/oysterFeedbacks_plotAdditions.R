
## Jacob L Moore (jlmoor@ucdavis.edu)
## 07/01/17

## Plots data from substrate and oyster addition simulations, with analytic approximations. 
## Creates Figures 6 and 7.

## Requires the following files:
	# WB_oysterOutput_final.csv
	# WB_substrateOutput_final.csv
	# WB_sensitivities_rhoVallow_alphaValestimate_deltaValmid.Rdata
	# WB_sensitivities_rhoValestimate_alphaValestimate_deltaValmid.Rdata

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

siteName = 'WB'
max.age=10
n = 250

rhoVals =  c('low', 'estimate')		# need both 'low' AND 'estimate'
deltaVal = 'mid' 					# 'low' OR 'mid' OR 'high' 
alphaVal = 'estimate'			# 'low' OR 'estimate' OR 'high'

startDistType = c('distHarvest', 'distEq')		# need both 'distHarvest' AND 'distEq

oyster.d = read.table(paste(data.folder, '/WB_oysterOutput_FINAL.csv', sep=''), header=TRUE, sep=',')
substrate.d = read.table(paste(data.folder, '/WB_substrateOutput_FINAL.csv', sep=''), header=TRUE, sep=',')

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Plot formatting
# -------------------------------------------------------------------
# -------------------------------------------------------------------

font.name = 'serif'
overall.scale = 1.5
main.scale = 1.05
axis.scale = 0.9
label.scale = 1.05
label.scale.matplot = 1

line.width = 2.5
line.width.eq = 1

col.two = c('gray47', 'gray70')
col.three = c('black', 'gray25', 'gray90' )
col.three.bar = rep(col.three, each=2)

age.cols = rainbow(9)
age.cols.bar = rep(age.cols, each=2)
age.ltys = c(1:9)

legend.text.age9 = c('Age=1', 'Age=2', 'Age=3', 'Age=4', 'Age=5', 'Age=6', 'Age=7', 'Age=8', 'Age=9')
legend.text.age10 = c('Age=1', 'Age=2', 'Age=3', 'Age=4', 'Age=5', 'Age=6', 'Age=7', 'Age=8', 'Age=9', 'Age=10')

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Plot oyster additions (Figure 6)
# -------------------------------------------------------------------
# -------------------------------------------------------------------

filename.oyster = paste('WB_deltaVal', deltaVal, '_alphaVal', alphaVal, sep='')			

# load analytic approximation for rho=low
filename.sens = paste('WB_sensitivities_rhoVal', 'low', '_alphaVal', alphaVal, '_deltaVal', deltaVal, sep='')
yy = load(paste(data.folder, '/', filename.sens, '.Rdata', sep=''))
oysters.req.low = oysters.req
rm(list=yy)

# load analytic approximation for rho=estimate
filename.sens = paste('WB_sensitivities_rhoVal', 'estimate', '_alphaVal', alphaVal, '_deltaVal', deltaVal, sep='')
yy = load(paste(data.folder, '/', filename.sens, '.Rdata', sep=''))
oysters.req.est = oysters.req
rm(list=yy)					

# subset oyster data so only using data for specific parameters 
oyster.sub2 = oyster.d[which(oyster.d$deltaVal==deltaVal),]
oyster.sub3 = oyster.sub2[which(oyster.sub2$alphaVal==alphaVal),]

# for plotting multiple y-axes
oyster.rhoLow = oyster.sub3[which(oyster.sub3$rhoVal=='low'),]
oyster.rhoLow.vals = oyster.rhoLow$oystersAdded
oyster.rhoLow.prop = oyster.rhoLow.vals / max(oyster.rhoLow.vals)
H.eqs.rhoLow = oyster.rhoLow$Heqs[1]
N.eqs.rhoLow = oyster.rhoLow$Neqs[1]
oyster.rhoEst = oyster.sub3[which(oyster.sub3$rhoVal=='estimate'),]
oyster.rhoEst.vals = oyster.rhoEst$oystersAdded
oyster.rhoEst.prop = oyster.rhoEst.vals / max(oyster.rhoEst.vals)
H.eqs.rhoEst = oyster.rhoEst$Heqs[1]
N.eqs.rhoEst = oyster.rhoEst$Neqs[1]
oysters.req.low.prop = oysters.req.low / max(oyster.rhoLow.vals)
oysters.req.est.prop = oysters.req.est / max(oyster.rhoEst.vals)
oyster.rhoLow.prop2 = oyster.rhoLow.vals / N.eqs.rhoLow
oyster.rhoEst.prop2 = oyster.rhoEst.vals / N.eqs.rhoEst
oysters.req.low.prop2 = oysters.req.low / N.eqs.rhoLow
oysters.req.est.prop2 = oysters.req.est / N.eqs.rhoEst
		
leftAxisLabs = pretty(seq(0, max(oyster.rhoLow.vals), length.out=10))
leftAxisAt = leftAxisLabs / max(oyster.rhoLow.vals)
rightAxisLabs = pretty(seq(0, max(oyster.rhoEst.vals), length.out=10))
rightAxisAt = rightAxisLabs / max(oyster.rhoEst.vals)		
yaxismod = 0.01
if(max(leftAxisAt, rightAxisAt)>1) yaxismod = max(leftAxisAt, rightAxisAt)-1 + 0.01*max(leftAxisAt, rightAxisAt)

# Figure 6A
quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale, mar=c(5,4,4,6)+0.1)
df.bar = barplot(rbind(oyster.rhoLow.prop, oyster.rhoEst.prop), ylim=c(0, max(oyster.rhoLow.prop, oyster.rhoEst.prop)+0.1), beside=TRUE, yaxt='n', main=paste('Oysters Required', sep=''), xlab='Age', col=age.cols.bar, names.arg=c(1:(max.age-1)), density=c(rep(c(100, 30), (max.age-1))))
points(df.bar, rbind(oysters.req.low.prop[1:9], oysters.req.est.prop[1:9]), pch=20, lwd=2, col=c('black', 'gray40'), cex=0.6)
legend('topright', legend=c(expression(rho['low']), expression(rho['est']), 'Analytic'), col=c('black', 'gray47', 'black'), density=c(100, 30, 0), pch=c(NA, NA, 20), border=c('black', 'black', 'white'), cex=0.75)
abline(h=0, lty=1, lwd=1, col='black')
abline(h=N.eqs.rhoLow/max(oyster.rhoLow.vals), lty=1, lwd=1, col='black')
abline(h=N.eqs.rhoEst/max(oyster.rhoEst.vals), lty=2, lwd=1, col='black')
axis(2, at = leftAxisAt, labels=leftAxisLabs)
mtext(expression(paste('Number of oysters, ', rho['low']), sep=''), side=2, line=3, cex = 1.5)
axis(4, at = rightAxisAt, labels=rightAxisLabs)
mtext(expression(paste('Number of oysters, ', rho['est']), sep=''), side=4, line=3, cex=1.5, srt=10)
box()
dev.copy(pdf, paste(filename.oyster, '_oystersRequired.pdf', sep=''))
dev.off()	

# Figure 6B
quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale, mar=c(5,4,4,6)+0.1)
df.bar = barplot(rbind(oysters.req.low.prop2[1:9]/oyster.rhoLow.prop2, oysters.req.est.prop2[1:9]/oyster.rhoEst.prop2), ylim=c(0, 2), beside=TRUE, xlab='Age', col=age.cols.bar, names.arg=c(1:(max.age-1)), density=c(rep(c(100, 30), (max.age-1))), main='Ratio of Estimates', ylab='Oysters required: analytic / numeric')
legend('topright', legend=c(expression(rho['low']), expression(rho['est'])), col=c('black', 'gray47'), density=c(100, 30, 0), pch=c(NA, NA), border=c('black', 'black'), cex=0.75)
abline(h=0, lty=1, lwd=1, col='black')
abline(h=1, lty=2, lwd=1, col='black')
box()
dev.copy(pdf, paste(filename.oyster, '_oystersRequired_compRATIO.pdf', sep=''))
dev.off()	
		

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Plot substrate additions (Figure 7)
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# subset oyster data so only using data for specific parameters 
substrate.sub2 = substrate.d[which(substrate.d$deltaVal==deltaVal),]
substrate.sub3 = substrate.sub2[which(substrate.sub2$alphaVal==alphaVal),]

# make one plot for each type of start distribution (Figure 7A -> distEq; Figure 7B -> distHarvest)
for(startDist in startDistType) {
	filename.substrate = paste('WB_deltaVal', deltaVal, '_alphaVal', alphaVal, '_', startDist, sep='')
	
	# subset oyster data so only using data for specific parameters 
	substrate.sub = substrate.sub3[which(substrate.sub3$startDist==startDist),]

	# for plotting multiple y-axes				
	substrate.low = substrate.sub[which(substrate.sub$rhoVal=='low'),]
	substrate.low.vals = substrate.low$substrateAdded
	substrate.low.prop = substrate.low.vals / max(substrate.low.vals)
	H.eqs.low = substrate.low$Heqs[1]
	H.eqs.low.prop = H.eqs.low / max(substrate.low.vals)
	N.eqs.low = substrate.low$Neqs[1]	
	substrate.est = substrate.sub[which(substrate.sub$rhoVal=='estimate'),]
	substrate.est.vals = substrate.est$substrateAdded
	substrate.est.prop = substrate.est.vals / max(substrate.est.vals)
	H.eqs.est = substrate.est$Heqs[1]
	H.eqs.est.prop = H.eqs.est / max(substrate.est.vals)
	N.eqs.est = substrate.est$Neqs[1]
	
	leftAxisLabs = pretty(seq(0, max(substrate.low.vals), length.out=10))
	leftAxisAt = leftAxisLabs / max(substrate.low.vals)
	rightAxisLabs = pretty(seq(0, max(substrate.est.vals), length.out=10))
	rightAxisAt = rightAxisLabs / max(substrate.est.vals)
	x.names = c(0.9, 1.1)
										
	# Figure 7
	quartz()
	par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale, mar=c(5,4,4,6)+0.2)
	df.bar = barplot(rbind(substrate.low.prop, substrate.est.prop), ylim=c(0, max(substrate.low.prop, substrate.est.prop)+0.1), beside=TRUE, yaxt='n', main=paste('Substrate Required \n', startDist), xlab= expression(paste('Initial population size (% of ', hat(N), ')', sep='')), col='gray30', names.arg=c('90%', '110%'), density=c(100,30))
	abline(h=0, lty=1, lwd=1, col='black')
	abline(h=H.eqs.low.prop, lty=1, lwd=2, col='gray30')
	abline(h=H.eqs.est.prop, lty=2, lwd=2, col='gray30')
	axis(2, at = leftAxisAt, labels=leftAxisLabs)
	mtext(expression(paste('Amount of substrate (', m^2, '), ', rho['low']), sep=''), side=2, line=3, cex = 1.5)
	axis(4, at = rightAxisAt, labels=rightAxisLabs)
mtext(expression(paste('Amount of substrate (', m^2, '), ', rho['est']), sep=''), side=4, line=3, cex=1.5, srt=10)
	legend('topright', legend=c(expression(rho['low']), expression(rho['est'])), density=c(100, 30), cex=0.75)
	box()
	dev.copy(pdf, paste(filename.substrate, '_substrateRequired.pdf', sep=''))
	dev.off()								
}				

