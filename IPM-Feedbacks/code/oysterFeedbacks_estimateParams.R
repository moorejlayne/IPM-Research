## Jacob L Moore (jlmoor@ucdavis.edu)
## 07/01/17

# Estimates parameters used in the model. Called by various scripts.

if(feedbacks!='none') {
	# load in values from linear model
	load(paste(data.folder, '/WB_eqVals_Linear.Rdata', sep=''))	

	# define remaining shell parameters
	if(deltaVal=='low') {
		params.shell$delta = 0.05
	} else if(deltaVal=='mid') {
		params.shell$delta = 0.2
	} else { # = 'high'
		params.shell$delta = 0.4
	}
	
	params.shell$p.hat = params$establishment.prob
		
	if(alphaVal=='low') params.shell$alpha = 100
	if(alphaVal=='high') params.shell$alpha = 1000000
	# else alphaVal=='estimate' and value retained from WB_eqVals_Linear.Rdata
		
	if(rhoVal=='low') {
		params.shell$rho = params$establishment.prob+0.5*params$establishment.prob
	} else {		# = 'estimate' from Puckett 2016
		params.shell$rho = 0.002617898 
	}
		
	# find equilibrium values	
	source('oysterFeedbacks_findEq.R')
	
} else {
	
	params=data.frame(
	  	surv.int=0,
	  	surv.age=0,
	  	surv.size=0,
	  	surv.age.size=0,
	   	growth.int=0,
	  	growth.var=0,
	  	growth.age=0,
	  	growth.size=0,
	  	growth.age.size=0,
	  	fec.slope = 0,
	  	fec.int = 0,
	  	sr.int = 0,
	  	sr.slope = 0,
	  	recruit.size.mean=0,
	  	recruit.size.sd=0,
	  	establishment.prob=1.8953272480974699859e-05	
	  					   
		)
	
	# delta, rho, and p.hat will be specified when run code with feedbacks=='positive'
	# alpha will be specified below
	params.shell = data.frame(
		scaling.a = 3.442786E-6,
		scaling.b = 1.56,
		delta = 0,
		alpha = 0,
		rho = 0,
		p.hat = 0
	)
	
	# -------------------------------------------------------------------
	# -------------------------------------------------------------------
	# Survival function
	# -------------------------------------------------------------------
	# -------------------------------------------------------------------
	
	# Surival data is given in growthsurvivalData.csv, and represents a binary survival (1) or mortality (0) value. See text for details of data. 	
	surv.fitting = function(d) {
		m = glm(surv~size, data=d, family=binomial)
		return(m)
	}
	surv.reg = surv.fitting(d)
	params$surv.int = coefficients(surv.reg)[1]
	params$surv.size = coefficients(surv.reg)[2]

	
	
	# -------------------------------------------------------------------
	# -------------------------------------------------------------------
	# Growth function
	# -------------------------------------------------------------------
	# -------------------------------------------------------------------
	
	# Growth data is given in growthsurvivalData.csv. See text for details of data. 	
	growth.fitting = function(d2) {
		optimMethod = 'Nelder-Mead'
		m = mle2(log.sizeChange~dnorm(mean=a+b*size, sd=c), start=list(a=mean(d2$log.sizeChange), b=0, c=sd(d2$log.sizeChange)), data=d2, method=optimMethod)
		return(m)
	}
	growth.reg = growth.fitting(d2)
	params$growth.int = coef(growth.reg)[1]		# a
	params$growth.size = coef(growth.reg)[2]		# b	
	params$growth.var = coef(growth.reg)[3]		# c		

	
	
	# -------------------------------------------------------------------
	# -------------------------------------------------------------------
	# Fecundity function
	# -------------------------------------------------------------------
	# -------------------------------------------------------------------
			
	data.virg.fec.raw = read.csv(paste(data.folder, '/fecundityData_WB.csv', sep=''), header=TRUE, sep=',') 
		
	clean.fecdata.virg = function(data_fec) {
		idx = seq(1:dim(data_fec)[1])
		date = as.Date(data_fec$Date[idx],"%m/%d/%y") 
		species = rep('CV', length(idx))
		location = rep('WB',length(idx))
		size.fec=data_fec$LVL..mm[idx]
		num.eggs = data_fec$Mean...Eggs[idx]				
		d.fec = data.frame(date, species, location, size.fec, num.eggs)
		return(d.fec)
	}		
	data.virg.fec = clean.fecdata.virg(data.virg.fec.raw)	

	# fit as a scaling relationship
	fec.reg = lm(log(num.eggs)~log(size.fec), data=data.virg.fec)
	params$fec.int=coefficients(fec.reg)[1]
	params$fec.slope=coefficients(fec.reg)[2]
	
	# recruit size distribution	
	data.sizes = read.csv(paste(data.folder, '/recruitSizes_WB.csv', sep=''), header=TRUE, sep=',')	
	# only take sizes from August recruits
	firstYearSizes = c(data.sizes$LVL[which(data.sizes$Cohort=='june2006')], data.sizes$LVL[which(data.sizes$Cohort=='june2007')])	
	params$recruit.size.mean=mean(firstYearSizes, na.rm=TRUE)
	params$recruit.size.sd=sd(firstYearSizes, na.rm=TRUE)
	
	# sex ratio calculations, from May 2007 and May 2008 (Mroch et al 2012)
	x.sr = c(15, 30, 45, 60, 75, 90, 105, 120, 135, 150)
	sr = c(0.1558, 0.3874, 0.682, 0.7667, 0.8974, 0.9579, 1, 0.98, 0.9375, 1)
	sr.reg =  lm(sr~x.sr)
	params$sr.int=coefficients(sr.reg)[1]
	params$sr.slope=coefficients(sr.reg)[2]	
	

	# -------------------------------------------------------------------
	# -------------------------------------------------------------------
	# Alpha calculation (data from Puckett et al. 2016; see text)
	# -------------------------------------------------------------------
	# -------------------------------------------------------------------
	
	popSizeData = read.csv(paste(data.folder, '/popSizeData_WB.csv', sep=''), header=TRUE, sep=',')	
	
	xs.bins = popSizeData$size
	N.total = 3583233	
	H.total = 4134.5
	adults = N.total*popSizeData$freq
	sex.ratio = params$sr.int + params$sr.slope*xs.bins
	sex.ratio[which(sex.ratio>1)]=1
	mean.num.offspring = exp(params$fec.int)*(xs.bins^params$fec.slope)*sex.ratio
	num.larvae = sum(adults*mean.num.offspring)
	recruits = 4266804
	rhoEst = 0.002617898
	params.shell$alpha = (num.larvae*rhoEst*H.total/recruits) - H.total

}