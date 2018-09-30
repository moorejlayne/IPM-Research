## Jacob L Moore (jlmoor@ucdavis.edu); 10/15/15

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Build kernels / IPM
# -------------------------------------------------------------------
# -------------------------------------------------------------------

if(statsMethod=='size') {
	# build growth kernel
	G.temp = deltax*outer(xs.mid.temp,xs.mid.temp,g.yx,params=params, currentAge=1)		# growth
	prob.growth.evict = 1-colSums(G.temp)	# set up discrete class for size > maxSize
	G = rbind(G.temp, prob.growth.evict)
	colAdd = c(rep(0,times=n), 1)
	G = cbind(G, colAdd)
	
	# build survival*growth kernel
	S = s.x(xs.mid.temp,params=params, currentAge=1)	# survival
	S = c(S, S[n])
	P= G; for(i in 1:(n+1)) P[,i]=G[,i]*S[i]	# growth * survival
	P.mat = P	# rename
	
	# build fecundity kernel
	F.temp = deltax*outer(xs.mid.temp,xs.mid.temp,f.yx, currentAge=1, siteName, params=params)
	F = rbind(F.temp, F.temp[n,])
	F = cbind(F, F.temp[,n])
	F.mat = F	# rename

	# build complete kernel that includes growth, survival, and fecundity
	K.mat = P + F%*%P	
} else { #age.size.int	
	# Set up individual matrices for survival, growth, and fecundity	
	S = array(0, dim=c((n+1),max.age))
	G = array(0, dim=c((n+1),(n+1),max.age))
	P = array(0, dim=c((n+1),(n+1),max.age))
	F = array(0, dim=c((n+1),(n+1),max.age))
	
	for(a in 1:max.age) {
		# build growth kernel
		G.temp = deltax*outer(xs.mid.temp,xs.mid.temp,g.yx, currentAge=a, params=params)	
		prob.growth.evict = 1-colSums(G.temp)
		G.temp = rbind(G.temp, prob.growth.evict)
		colAdd = c(rep(0,times=n), 1)
		G.temp = cbind(G.temp, colAdd)	
		G[,,a] = G.temp
		
		# build survival*growth kernel
		S.temp = s.x(xs.mid.temp,currentAge=a,params=params, max.age=max.age)
		S.temp = c(S.temp, S.temp[n])
		P.temp = G.temp; for(i in 1:n) P.temp[,i]=G.temp[,i]*S.temp[i]
		P[,,a]=P.temp
	
		# build fecundity kernel
		F.temp = deltax*outer(xs.mid.temp,xs.mid.temp,f.yx, currentAge=a, siteName, params=params)
		F.temp = rbind(F.temp, F.temp[n,])
		F.temp = cbind(F.temp, F.temp[,n])
		F[,,a]=F.temp
	}
	
	# build complete kernel that includes growth, survival, and fecundity 
	zMat = matrix(0, nrow=(n+1), ncol=(n+1))	# filler matrix of zeros

	# row of fecundity matrices
	F.mat = zMat	# first age doesn't reproduce
	for(r in 2:max.age) {
		F.mat = cbind(F.mat, F[,,r])	
	}
	zMat2 = matrix(0, nrow=(dim(F.mat)[2]-(n+1)), ncol=dim(F.mat)[2])
	F.mat = rbind(F.mat, zMat2)
	
	# build up rest of matrix that contains growth*survival
	P.cols = array(dim=c((max.age-1)*(n+1),(n+1), max.age))
	# first column
	num.zeros.post = max.age-2
	temp.col = P[,,1]
	for (r.loc in 1:num.zeros.post) temp.col=rbind(temp.col,zMat)
	P.cols[,,1]=temp.col
	# middle columns
	for(col.loc in 2:(max.age-1)) {
		num.zeros.pre = col.loc-1
		num.zeros.post = max.age-col.loc-1
		temp.col=zMat
		if(num.zeros.pre>1) {
			for(v in 2:num.zeros.pre) temp.col=rbind(temp.col,zMat)
		}
		temp.col = rbind(temp.col, P[,,col.loc])
		if(num.zeros.post>0) {
			for(w in 1:num.zeros.post) temp.col=rbind(temp.col,zMat)
		}
		P.cols[,,col.loc]=temp.col
	}
	# last column
	temp.col=zMat
	for (r.loc in 2:(max.age-1)) temp.col=rbind(temp.col,zMat)
	P.cols[,,max.age]=temp.col
	# join columns together
	P.mat = P.cols[,,1]
	for(r in 2:max.age) P.mat = cbind(P.mat,P.cols[,,r])
	zMat3 = matrix(0, nrow=(n+1), ncol=dim(P.mat)[2])
	P.mat = rbind(zMat3, P.mat)
	
	# final matrix
	K.mat = P.mat + F.mat%*%P.mat	
}

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Evaluate IPM
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# Stable distribution
init.vect2 = rep(1, dim(K.mat)[1])	
num.calcs = 1000
for(i in 1:num.calcs){
	temp.vect2 = K.mat %*% init.vect2
	init.vect2 = (temp.vect2 / max(temp.vect2))[,1]
	}	
v = init.vect2	
stable.dist = v / sum(v)

# Long-term population growth rate
lambda = sum(K.mat %*% v) / sum(v)

# Reproductive values		
init.vect1 = rep(1, dim(K.mat)[1])	
num.calcs = 1000
for(i in 1:num.calcs){
	temp.vect1 =  init.vect1 %*% K.mat
	init.vect1 = (temp.vect1 / max(temp.vect1))[1,]
	}
w = init.vect1	
rep.vals = w / sum(stable.dist*w)

# Elasticity
v.dot.w = sum(stable.dist*rep.vals)
sens = outer(rep.vals, stable.dist) / v.dot.w
elas = sens*K.mat / lambda
elas.F = sens %*% t(P.mat) * F.mat / lambda
elas.P = ( sens + t(F.mat) %*% sens) * P.mat / lambda

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Normalize values, and break up into age and size
# -------------------------------------------------------------------
# -------------------------------------------------------------------

if(statsMethod=='size') {
	stable.size.overall = stable.dist[1:(n+1)]
	repVals.size.overall = rep.vals[1:(n+1)]
	repVals.size.overall.scaled = repVals.size.overall / sum(repVals.size.overall)
	elas.size = colSums(elas)
	elas.F.size = colSums(elas.F)
	elas.P.size = colSums(elas.P)
} else {	# age.size.int
	# ---------- Stable distributions
	stable.age = matrix(nrow=1, ncol=max.age)
	stable.size.temp = matrix(nrow=max.age, ncol=(n+1))
	stable.age[1]=sum(stable.dist[1:(n+1)])
	stable.size.temp[1,]=stable.dist[1:(n+1)]		
	stable.temp = stable.dist[1:(n+1)]
	stable.scaled = stable.temp / sum(stable.temp)
	for(currentAge in 2:max.age) {
		start.idx = (n+1)*(currentAge-1)+1
		end.idx = (n+1)*currentAge
		stable.temp = stable.dist[start.idx:end.idx]
		stable.scaled = stable.temp / sum(stable.temp)
		stable.age[currentAge]=sum(stable.dist[start.idx:end.idx])
		stable.size.temp[currentAge,]=stable.dist[start.idx:end.idx]
	}
	stable.size.overall = colSums(stable.size.temp)
	
	# ---------- Reproductive values
	repVals.age = matrix(nrow=1, ncol=(max.age))
	repVals.size.temp = matrix(nrow=(max.age), ncol=(n+1))
	repVals.age[1]=sum(rep.vals[(1):(n+1)])
	repVals.size.temp[1,]=rep.vals[(1):(n+1)]	
	for(currentAge in 2:max.age) {
		start.idx = (n+1)*(currentAge-1)+1
		end.idx = (n+1)*currentAge
		repVals.age[currentAge]=sum(rep.vals[start.idx:end.idx])
		repVals.size.temp[currentAge,]=rep.vals[start.idx:end.idx]
	}
	repVals.size.overall = colSums(repVals.size.temp)
	repVals.size.overall.scaled = repVals.size.overall / sum(repVals.size.overall)
	repVals.age.scaled = repVals.age / sum(repVals.age)
	
	# ---------- Survival elasticity
	elas.P.colsum = colSums(elas.P)
	elas.P.size = matrix(0, nrow=1, ncol=(n+1))
	elas.P.age = matrix(0, nrow=1,ncol=max.age)
	elas.P.size[(1:(n+1))] = elas.P.colsum[1:(n+1)]
	elas.P.age[1] = sum(elas.P.colsum[1:(n+1)])
	for(currentAge in 2:max.age) {
		start.idx = (n+1)*(currentAge-1)+1
		end.idx = (n+1)*currentAge
		elas.P.size[1:(n+1)] = elas.P.size + elas.P.colsum[start.idx:end.idx]
		elas.P.age[currentAge] = sum(elas.P.colsum[start.idx:end.idx])
	}
	
	# ---------- Fecundity elasticity
	elas.F.colsum = colSums(elas.F)
	elas.F.size = matrix(0, nrow=1, ncol=(n+1))
	elas.F.age = matrix(0, nrow=1,ncol=max.age)
	elas.F.size[(1:(n+1))] = elas.F.colsum[1:(n+1)]
	elas.F.age[1] = sum(elas.F.colsum[1:(n+1)])
	for(currentAge in 2:max.age) {
		start.idx = (n+1)*(currentAge-1)+1
		end.idx = (n+1)*currentAge
		elas.F.size[1:(n+1)] = elas.F.size + elas.F.colsum[start.idx:end.idx]
		elas.F.age[currentAge] = sum(elas.F.colsum[start.idx:end.idx])
	}
}

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Save data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
source('oysterIPM_saveData.R')
