## Jacob L Moore (jlmoor@ucdavis.edu)
## 07/01/17

## Creates plots of initial size/age distributions, trajectories of oyster number and substrate levels, 
## and phase planes for a single run of the IPM. 

## Called by oysterFeedbacks.R

graphics.off()

# ---------- Plot formatting specifications

font.name = 'serif'
overall.scale = 1.65
main.scale = 1.05
axis.scale = 0.9
label.scale = 1.05
line.width.axes = 1
legend.scale = 0.6

line.width = 2.5
line.width.eq = 1
line.type.above = 1
line.type.below = 2

pch.eq = 21
pch.start = 17
pch.end = 15

legend.text.age9 = c('Age=1', 'Age=2', 'Age=3', 'Age=4', 'Age=5', 'Age=6', 'Age=7', 'Age=8', 'Age=9')
legend.text.age10 = c('Age=1', 'Age=2', 'Age=3', 'Age=4', 'Age=5', 'Age=6', 'Age=7', 'Age=8', 'Age=9', 'Age=10')
				
if(feedbacks=='none') {
	
	quartz()
	par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
	plot(xs.mid, stable.size.dist, type='l', lwd=line.width, col='black', xlab='Size', ylab='Frequency', main=paste('Stable Size Distribution', sep=''))
	box(lwd=line.width.axes)
	dev.copy(pdf, paste(filename, '_sizeDist.pdf', sep=''))
	dev.off()
	
	quartz()
	par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
	plot(1:max.age, stable.age.dist, type='b', lwd=line.width, col='black', xlab='Age', ylab='Frequency', main=paste('Stable Age Distribution', sep=''))
	box(lwd=line.width.axes)
	dev.copy(pdf, paste(filename, '_ageDist.pdf', sep=''))
	dev.off()
				
} else {

##------ Plot starting size and age distributions
	quartz()
	par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
	plot(xs.mid, rowSums(init.vect), type='l', lwd=line.width, col='black', xlab='Size', ylab='Frequency', main=paste('Starting Size Distribution', sep=''))
	box(lwd=line.width.axes)
	dev.copy(pdf, paste(filename, '_startingSizeDist.pdf', sep=''))
	dev.off()
	
	quartz()
	par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
	plot(1:max.age, colSums(init.vect), type='b', lwd=line.width, col='black', xlab='Size', ylab='Frequency', main=paste('Starting Age Distribution', sep=''))
	box(lwd=line.width.axes)
	dev.copy(pdf, paste(filename, '_startingAgeDist.pdf', sep=''))
	dev.off()

##------ Plot trajectories of oyster numbers and substrate amounts
	quartz()
	par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
	plot(1:plot.length, N.sum[1:plot.length], main=paste('Oyster Trajectory'), xlab='Time', ylab='Oysters', lwd=line.width, lty=line.type.above, type='line', ylim=c(min(min(N.sum[1:plot.length]), N.eqs), max(max(N.sum[1:plot.length]), N.eqs)))
	abline(h=N.eqs, lwd=line.width.eq, lty=3, col='gray50')
	dev.copy(pdf, paste(filename, '_oysters.pdf', sep=''))
	dev.off()
		
	quartz()
	par(bg='white', family=font.name, cex = overall.scale, cex.main = main.scale, cex.lab = label.scale, cex.axis=axis.scale)
	plot(1:plot.length, H.sum[1:plot.length], main=paste('Substrate Trajectory'), xlab='Time', ylab='Substrate', lwd=line.width, lty=line.type.above, type='line', ylim=c(min(min(H.sum[1:plot.length]), H.eqs), max(max(H.sum[1:plot.length]), H.eqs)))
	abline(h=H.eqs, lwd=line.width.eq, lty=3, col='gray50')
	dev.copy(pdf, paste(filename, '_substrate.pdf', sep=''))
	dev.off()

##------ Plot phase plane of oyster numbers and substrate amounts
	if(dist.type=='distHarvest') inset.loc = 'topleft'
	if(dist.type=='distEq') inset.loc = 'topleft'

	quartz()
	par(bg='white', family=font.name, cex = overall.scale, cex.main=main.scale, cex.lab=label.scale, cex.axis=axis.scale)
	plot(H.sum[1:plot.length], N.sum[1:plot.length], type='l', lwd=line.width, lty=line.type.above, ylab='Total Oysters', xlab= expression("Substrate ("*""*m^2*")"), xlim=c(min(min(H.sum[1:plot.length]), H.eqs), max(max(H.sum[1:plot.length]), H.eqs)), ylim=c(min(min(N.sum[1:plot.length]), N.eqs), max(max(N.sum[1:plot.length]), N.eqs)), main='Phase Plane')
	points(H.eqs, N.eqs, pch=c(pch.eq), col=c('black'))
	points(H.sum[plot.length], N.sum[plot.length], lwd=line.width, pch=c(pch.end), col='black')
	points(H.sum[1], N.sum[1], lwd=line.width, pch=c(pch.start), col='black')
	legend('bottomright', legend=c('Unstable Equilibrium', 'time=0', paste('time=', plot.length, sep='')), col=c('black', 'black', 'black'), pch=c(pch.eq, 24, 22), cex=legend.scale)
	subplot(plot(xs.mid, rowSums(dist.mat), type='l', lwd=1, col='black', main='Initial Distribution', xlab='Size', ylab='Frequency'), type='plt', inset.loc, pars=list(cex=0.95, cex.main=0.75, cex.lab=.75, cex.axis=0.5, tcl=-0.25, mgp=c(1.1,0.25,0), bty='n', mai=c(1,1,0.25,1)), inset=c(0.1,0.08))
	dev.copy(pdf, paste(filename, '_phasePlane.pdf', sep=''))
	dev.off()
		
}
