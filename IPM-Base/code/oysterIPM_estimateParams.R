## Jacob L Moore (jlmoor@ucdavis.edu); 10/15/15

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
  	recruit.size.mean=0,
  	recruit.size.sd=0,
  	establishment.prob=0
)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Survival function
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# Surival data is given in growthsurvivalData.csv, and represents a binary survival (1) or mortality (0) value. See text for details of data. 

surv.fitting = function(statsMethod, d) {
	if(statsMethod =='size') {
		m = glm(surv~size, data=d, family=binomial)
	} else { 	# 'age.size.int'
		m = glm(surv~size+age+size*age, data=d, family=binomial)
	}
	return(m)
}
surv.reg = surv.fitting(statsMethod, d)
if(statsMethod =='age.size.int') {
	params$surv.int = coefficients(surv.reg)[1]
	params$surv.size = coefficients(surv.reg)[2]	
	params$surv.age = coefficients(surv.reg)[3]	
	params$surv.age.size = coefficients(surv.reg)[4]			
} else {
	params$surv.int = coefficients(surv.reg)[1]
	params$surv.size = coefficients(surv.reg)[2]
}


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Growth function
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# Growth data is given in growthsurvivalData.csv. See text for details of data. 

growth.fitting = function(statsMethod, d2) {
	optimMethod = 'Nelder-Mead'
	if(statsMethod=='size') {
		m = mle2(log.sizeChange~dnorm(mean=a+b*size, sd=c), start=list(a=mean(d2$log.sizeChange), b=0, c=sd(d2$log.sizeChange)), data=d2, method=optimMethod)
	} else { 	# 'age.size.int'
		m = mle2(log.sizeChange~dnorm(mean=a+b*age+d*size+h*size*age, sd=c), start=list(a=mean(d2$log.sizeChange), b=0, d=0, h=0, c=sd(d2$log.sizeChange)), data=d2, method=optimMethod)
	}
}
growth.reg = growth.fitting(statsMethod, d2)
if(statsMethod=='age.size.int') {
	params$growth.int = coef(growth.reg)[1]				# a
	params$growth.age = coef(growth.reg)[2]				# b
	params$growth.size = coef(growth.reg)[3]			# d
	params$growth.age.size = coef(growth.reg)[4]		# h
	params$growth.var = coef(growth.reg)[5]				# c
} else {
	params$growth.int = coef(growth.reg)[1]				# a
	params$growth.size = coef(growth.reg)[2]			# b	
	params$growth.var = coef(growth.reg)[3]				# c		
}


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Fecundity function
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# Fecundity data is given in fecundityData.csv. Data on the relationship between dry tissue weight (DTW) and the number of eggs produced in millions (num.eggs) was obtained from Kang et al. 2003. To get the relationship between shell length in mm (size.fec) and the number of eggs, we used to relationship given in Ren et al. 2003 to convert DTW to shell length: size = 10 * ( (DTW/4.22E-4) ^ (1/3.74) ). 

data.gig.fec = read.table('fecundityData.csv', header=TRUE, sep=',')
data.gig.fec$num.eggs =	data.gig.fec$num.eggs*1E6 
fec.reg = lm(log(num.eggs)~size.fec, data=data.gig.fec)
params$fec.int=coefficients(fec.reg)[1]
params$fec.slope=coefficients(fec.reg)[2]	

# distribution of offspring sizes	
firstYearSizes = d$size[which(d$timeStep==d$timeStep[1])]	# pull out data from first time point
params$recruit.size.mean=mean(firstYearSizes, na.rm=TRUE)
params$recruit.size.sd=sd(firstYearSizes, na.rm=TRUE)
