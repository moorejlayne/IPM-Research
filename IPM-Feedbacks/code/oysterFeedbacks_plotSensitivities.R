## Jacob L Moore (jlmoor@ucdavis.edu)
## 07/01/17

## Plots sensitivity of the threshold surface to delta, rho, and alpha when rho=low and rho=est.
## Creates Figure 4.

## Requires the following files:
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

n = 250
max.age=10

deltaVal = 'mid'			# 'low' OR 'mid' OR 'high' [given that sensitivities have been run]
alphaVal = 'estimate'	# 'low' OR 'estimate' OR 'high' [given that sensitivities have been run]

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Plot formatting
# -------------------------------------------------------------------
# -------------------------------------------------------------------

font.name = 'serif'
overall.scale = 1.75
main.scale = 1.05
axis.scale = .8
label.scale = 1.5
label.scale.matplot = 1

line.width = 2.5
line.width.eq = 1

age.cols = rainbow(9)
age.cols.bar = rep(age.cols, each=2)
age.ltys = c(1:9)

legend.text.age9 = c('Age=1', 'Age=2', 'Age=3', 'Age=4', 'Age=5', 'Age=6', 'Age=7', 'Age=8', 'Age=9')
legend.text.age10 = c('Age=1', 'Age=2', 'Age=3', 'Age=4', 'Age=5', 'Age=6', 'Age=7', 'Age=8', 'Age=9', 'Age=10')

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Make plot
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# load in data for rho=low
rhoVal = 'low'
filename = paste('/WB_sensitivities_rhoVal', rhoVal, '_alphaVal', alphaVal, '_deltaVal', deltaVal, sep='')
yy = load(paste(data.folder, filename, '.Rdata', sep=''))
delta.low = deltaE
alpha.low = alphaE
rho.low = rhoE

# load in data for rho=est
rhoVal = 'estimate'
filename = paste('/WB_sensitivities_rhoVal', rhoVal, '_alphaVal', alphaVal, '_deltaVal', deltaVal, sep='')
yy = load(paste(data.folder, filename, '.Rdata', sep=''))
delta.est = deltaE
alpha.est = alphaE
rho.est = rhoE	

# filename for saving			
filename2 = paste('WB_sensitivities', '_deltaVal', deltaVal, '_alphaVal', alphaVal, sep='')

# Plot elasticities; Figure 4
quartz()
par(bg='white', family=font.name, cex = overall.scale, cex.main = 0.6, cex.lab =label.scale, cex.axis=axis.scale)
elas.mat1 = c(alpha.low[1], rho.low[1], delta.low[1])
elas.mat2 = c(alpha.est[1], rho.est[1], delta.est[1])
elas.mat = rbind(elas.mat1, elas.mat2)
barplot(elas.mat, names.arg=c(expression(alpha), expression(rho), expression(delta)), ylab='Elasticity', beside=TRUE, ylim=c(-3.5, 1.25))
abline(h=0)
legend('bottomright', legend=c(expression(rho['low']), expression(rho['est'])), col=c('gray10', 'gray70'), pch=c(20,20), cex=0.6)
box()
dev.copy(pdf, paste(filename2, '_Elasticity.pdf', sep=''))
dev.off()