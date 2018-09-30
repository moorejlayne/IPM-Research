## Jacob L Moore (jlmoor@ucdavis.edu)
## 07/01/17

## Contains functions for creating growth, survival, and fecundity kernels. Called by various scripts.

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Functions for life history 
# -------------------------------------------------------------------
# -------------------------------------------------------------------
	
# survival function (probability of surviving)
s.x = function(currentSize, currentAge, params, max.age) {
	# this ensures that all oysters that die at max age will enter dead shell class
	if(currentAge==max.age) {
		u=0*currentSize
	} else {
		u = exp(params$surv.int + currentAge*params$surv.age + currentSize*params$surv.size + currentSize*currentAge*params$surv.age.size)
	}
	return(u/(1+u))
}


# growth function 
g.yx = function(futureSize, currentSize, currentAge, params) {	
	lnorm.mean=params$growth.int + currentAge*params$growth.age + currentSize*params$growth.size + currentAge*currentSize*params$growth.age.size
	lnorm.sd=params$growth.var
	upperSize = futureSize+(deltax/2)
	lowerSize = futureSize-(deltax/2)
	dist = plnorm(upperSize-currentSize, mean=lnorm.mean, sd=lnorm.sd)-plnorm(lowerSize-currentSize, mean=lnorm.mean, sd=lnorm.sd)
}	

# fecundity function
f.yx = function(futureSize, currentSize, currentAge, params, params.shell, feedbacks, currentPop, currentShell) {
	sex.ratio = params$sr.int + params$sr.slope*currentSize
	sex.ratio[which(sex.ratio>1)]=1
	mean.num.offspring = exp(params$fec.int)*(currentSize^params$fec.slope)*sex.ratio
		
	if(feedbacks=='none') {
		total.offspring = mean.num.offspring * params$establishment.prob
	} else {
		total.offspring = mean.num.offspring * ((params.shell$rho * currentShell)/(params.shell$alpha + currentShell))		
	} 		
	return(total.offspring * dnorm(futureSize, mean=params$recruit.size.mean, sd=params$recruit.size.sd))
}  
