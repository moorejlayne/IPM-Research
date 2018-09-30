## Jacob L Moore (jlmoor@ucdavis.edu); 10/15/15

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# General functions for life history 
# -------------------------------------------------------------------
# -------------------------------------------------------------------
	
# survival function (probability of surviving)
	s.x = function(currentSize, currentAge, params, max.age=100) {
		u = exp(params$surv.int + currentAge*params$surv.age + currentSize*params$surv.size + currentSize*currentAge*params$surv.age.size)
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
	f.yx = function(futureSize, currentSize, currentAge, siteName, params) {
		# fecundity = establishment prob * mean num offspring * prob distribution of offspring size y for adult of size x 	
		sex.ratio = 0.0310907 + 0.0043945*currentSize   # relationship obtained by fitting data from Buroker 1983 (see data in sexRatiosWA.csv)
	 	if(sex.ratio>1) sex.ratio=1
		mean.num.offspring = exp(params$fec.int)*exp(params$fec.slope*currentSize)*sex.ratio
		return(params$establishment.prob * mean.num.offspring * dnorm(futureSize, mean=params$recruit.size.mean, sd=params$recruit.size.sd))
	}  
