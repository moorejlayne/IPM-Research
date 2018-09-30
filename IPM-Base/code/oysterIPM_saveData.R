## Jacob L Moore (jlmoor@ucdavis.edu); 10/15/15

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Save data
# -------------------------------------------------------------------
# -------------------------------------------------------------------

if(statsMethod=='size') {
	filename.data = paste('size_maxage', max.age, '_maxsize', max.size, '_lambda', lam, sep='')
		if(lam=='low') {
		lambda.low = lambda
		stable.size.overall.low = stable.dist
		repVals.size.overall.low = rep.vals
		repVals.size.overall.scaled.low = repVals.size.overall.scaled
		elas.low = elas
		elas.size.low = elas.size
		elas.F.low = elas.F
		elas.P.low = elas.P
		sens.low = sens
		elas.P.size.low = elas.P.size
		elas.F.size.low = elas.F.size
		params.low = params
		F.mat.low = F.mat
		K.mat.low = K.mat
		P.mat.low = P.mat
		save(elas.size.low,repVals.size.overall.scaled.low, elas.F.low, elas.P.low, params.low, F.mat.low, K.mat.low, P.mat.low, n, max.size, max.age, xs.mid, lambda.low, stable.size.overall.low, repVals.size.overall.low, sens.low, elas.low, elas.P.size.low, elas.F.size.low, file=paste(filename.data, '.RData', sep=''))			
	} else if(lam=='mid') {
		lambda.mid = lambda
		stable.size.overall.mid = stable.dist
		repVals.size.overall.mid = rep.vals
		repVals.size.overall.scaled.mid = repVals.size.overall.scaled
		elas.mid = elas
		elas.size.mid = elas.size
		elas.P.mid = elas.P
		elas.F.mid = elas.F
		sens.mid = sens
		elas.P.size.mid = elas.P.size
		elas.F.size.mid = elas.F.size
		params.mid = params
		F.mat.mid = F.mat
		K.mat.mid = K.mat
		P.mat.mid = P.mat					
		save(elas.size.mid,repVals.size.overall.scaled.mid, elas.P.mid, elas.F.mid, F.mat.mid, K.mat.mid, P.mat.mid, params.mid, n, max.size, max.age, xs.mid, lambda.mid, stable.size.overall.mid, repVals.size.overall.mid, sens.mid, elas.mid, elas.P.size.mid, elas.F.size.mid, file=paste(filename.data, '.RData', sep=''))					
	} else {
		lambda.high = lambda
		stable.size.overall.high = stable.dist
		repVals.size.overall.high = rep.vals
		repVals.size.overall.scaled.high = repVals.size.overall.scaled
		elas.high = elas
		elas.size.high = elas.size
		elas.P.high = elas.P
		elas.F.high = elas.F
		sens.high = sens
		elas.P.size.high = elas.P.size
		elas.F.size.high = elas.F.size
		params.high = params
		F.mat.high = F.mat
		K.mat.high = K.mat
		P.mat.high = P.mat
		save(elas.size.high,repVals.size.overall.scaled.high, elas.P.high, elas.F.high, params.high, F.mat.high, K.mat.high, P.mat.high, n, max.size, max.age, xs.mid, lambda.high, stable.size.overall.high, repVals.size.overall.high, sens.high, elas.high, elas.P.size.high, elas.F.size.high, file=paste(filename.data, '.RData', sep=''))				
	}

} else { 	# age.size.int
	filename.data = paste('agesize_maxage', max.age, '_maxsize', max.size, '_lambda', lam, sep='')
	if(lam=='low') {
		lambda.low = lambda
		stable.dist.low = stable.dist
		rep.vals.low = rep.vals
		sens.low = sens
		elas.low = elas
		elas.P.low = elas.P
		elas.F.low = elas.F
		stable.size.overall.low = stable.size.overall
		stable.age.low = stable.age
		repVals.age.low = repVals.age
		repVals.size.overall.low = repVals.size.overall
		repVals.size.overall.scaled.low = repVals.size.overall.scaled
		repVals.age.scaled.low = repVals.age.scaled
		elas.P.size.low = elas.P.size
		elas.P.age.low = elas.P.age
		elas.F.size.low = elas.F.size
		elas.F.age.low = elas.F.age
		params.low = params
		F.mat.low = F.mat
		K.mat.low = K.mat
		P.mat.low = P.mat
		save(elas.P.size.low, elas.P.age.low, elas.F.size.low, elas.F.age.low, elas.P.low, elas.F.low, params.low, F.mat.low, K.mat.low, P.mat.low, n, max.age, max.size, xs.mid, stable.age.low, lambda.low, stable.dist.low, rep.vals.low, sens.low, elas.low, stable.size.overall.low, stable.age.low, repVals.age.low, repVals.size.overall.low, repVals.size.overall.scaled.low, repVals.age.scaled.low, file=paste(filename.data, '.RData', sep=''))

	} else if(lam=='mid') {
		lambda.mid = lambda
		stable.dist.mid = stable.dist
		rep.vals.mid = rep.vals
		sens.mid = sens
		elas.mid = elas
		elas.P.mid = elas.P
		elas.F.mid = elas.F
		stable.size.overall.mid = stable.size.overall
		stable.age.mid = stable.age
		repVals.age.mid = repVals.age
		repVals.size.overall.mid = repVals.size.overall
		repVals.size.overall.scaled.mid = repVals.size.overall.scaled
		repVals.age.scaled.mid = repVals.age.scaled
		elas.P.size.mid = elas.P.size
		elas.P.age.mid = elas.P.age
		elas.F.size.mid = elas.F.size
		elas.F.age.mid = elas.F.age
		params.mid = params
		F.mat.mid = F.mat
		K.mat.mid = K.mat
		P.mat.mid = P.mat
		save(elas.P.size.mid, elas.P.age.mid, elas.F.size.mid, elas.F.age.mid, elas.P.mid, elas.F.mid, params.mid, F.mat.mid, P.mat.mid, K.mat.mid, n, max.age, max.size, xs.mid, stable.age.mid, lambda.mid, stable.dist.mid, rep.vals.mid, sens.mid, elas.mid, stable.size.overall.mid, stable.age.mid, repVals.age.mid, repVals.size.overall.mid, repVals.size.overall.scaled.mid, repVals.age.scaled.mid, file=paste(filename.data, '.RData', sep=''))

	} else {
		lambda.high = lambda
		stable.dist.high = stable.dist
		rep.vals.high = rep.vals
		sens.high = sens
		elas.high = elas
		elas.P.high = elas.P
		elas.F.high = elas.F
		stable.size.overall.high = stable.size.overall
		stable.age.high = stable.age
		repVals.age.high = repVals.age
		repVals.size.overall.high = repVals.size.overall
		repVals.size.overall.scaled.high = repVals.size.overall.scaled
		repVals.age.scaled.high = repVals.age.scaled
		elas.P.size.high = elas.P.size
		elas.P.age.high = elas.P.age
		elas.F.size.high = elas.F.size
		elas.F.age.high = elas.F.age
		params.high = params
		F.mat.high = F.mat
		K.mat.high = K.mat
		P.mat.high = P.mat
		save(elas.F.high, elas.P.high, params.high, F.mat.high, K.mat.high, P.mat.high, n, max.age, max.size, xs.mid, stable.age.high, lambda.high, stable.dist.high, rep.vals.high, sens.high, elas.high, stable.size.overall.high, stable.age.high, repVals.age.high, repVals.size.overall.high, repVals.size.overall.scaled.high, repVals.age.scaled.high, elas.P.size.high, elas.P.age.high, elas.F.size.high, elas.F.age.high, file=paste(filename.data, '.RData', sep=''))
	}
}

save.image()
unlink(".RData")	
