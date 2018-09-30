## Jacob L Moore (jlmoor@ucdavis.edu) 
## 07/01/17

# Builds and runs the integral-projection model. Called by various scripts. 

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Build kernels / IPM
# -------------------------------------------------------------------
# -------------------------------------------------------------------

	# Set up individual matrices for survival, growth, and fecundity	
	S = array(0, dim=c((n+1),max.age))
	G = array(0, dim=c((n+1),(n+1),max.age))
	P = array(0, dim=c((n+1),(n+1),max.age))
	F = array(0, dim=c((n+1),(n+1),max.age))
	H.new.base = S
	
	for(a in 1:max.age) {
		# build growth kernel
		G.temp = deltax*outer(xs.mid.temp,xs.mid.temp,g.yx, currentAge=a, params=params)	
		prob.growth.evict = (1-colSums(G.temp))	
		G.temp = rbind(G.temp, prob.growth.evict)
		colAdd = c(rep(0,times=n), 1)
		G.temp = cbind(G.temp, colAdd)	
		G[,,a] = G.temp
		
		# build survival and shell kernel
		S.temp = s.x(xs.mid.temp,currentAge=a,params=params, max.age=max.age)
		S.temp = c(S.temp, S.temp[n])
		S[,a] = S.temp
		
		# amount of new shell added to dead shell variable (before multiplied by population size)
		H.new.base[,a] = params.shell$scaling.a*(xs.mid^params.shell$scaling.b) * (G[,,a] %*% (1-S[,a]))
				
		# build survival*growth kernel
		P.temp = G.temp; for(i in 1:n) P.temp[,i]=G.temp[,i]*S[i,a]
		P[,,a]=P.temp
	
		# build fecundity kernel
		F.temp = deltax*outer(xs.mid.temp,xs.mid.temp,f.yx, currentAge=a, params=params, params.shell=params.shell, feedbacks, currentPop=NA, currentShell=NA)
		F.temp = rbind(F.temp, F.temp[n,])
		F.temp = cbind(F.temp, F.temp[,n])
		F[,,a]=F.temp		 	
	}

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Simulate Population (Evaluate IPM)
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# set up matrix to iterate through time steps
N.mat = array(NA, dim=c((n+1),max.age, num.calcs))
	
# vectors to keep track of total population size and shell at each time step
N.sum = rep(NA, num.calcs)
H.sum = N.sum

# add in intial size vectors
N.mat[,,1] = init.vect
N.sum[1] = sum(N.mat[,,1])
H.sum[1] = init.vect.dead

# simulate population 
for (i in 2:num.calcs) {
	#print(i)

## ---- update shell substrate
	# new shell from growth/survival of oysters from previous census
	H.new = array(0, dim=c((n+1),max.age))		
	H.new[,1] = H.new.base[,1] * N.mat[,1,(i-1)]	
	for(a in 2:max.age) {
		H.new[,a] = H.new.base[,a] * N.mat[,a,(i-1)]
	}
	
	# decay of shell from previous year
	H.sum[i] = H.sum[i-1]*(exp(-params.shell$delta)) + sum(H.new)
	
	# new shell added
	if(time.add.shell==i) {
		H.sum[i] = H.sum[i] + shell.addition[1]
	}	
	currentShell = H.sum[i]

## ---- simulate population		
	currentPop = N.sum[i-1]	
	N.new = matrix(0, nrow=(n+1), ncol=max.age)
	for(a in 2:max.age) {
		if(feedbacks!='none') {
			# build fecundity kernel
			F.temp = deltax*outer(xs.mid.temp,xs.mid.temp,f.yx, currentAge=a, params=params, params.shell=params.shell, feedbacks, currentPop=currentPop, currentShell=currentShell)
			F.temp = rbind(F.temp, F.temp[n,])
			F.temp = cbind(F.temp, F.temp[,n])
			F[,,a]=F.temp
		}	
				
		# iterate population
		N.new[,a] = ((F[,,a] %*% P[,,(a-1)]) %*% N.mat[,(a-1),(i-1)])		
		N.mat[,a,i] = P[,,(a-1)] %*% N.mat[,(a-1),(i-1)]
	}
	N.mat[,1,i] = rowSums(N.new)
	
	# Add additional oysters
	if(time.add.oyster==i) {
		N.mat[,,i]=N.mat[,,i]+oyster.addition[,,1]
	}
	N.sum[i] = sum(N.mat[,,i])

	# rescale linear model
	if(feedbacks=='none') {
		N.mat[,,i] = N.mat[,,i] / sum(N.mat[,,i])
	}		
}

stable.dist.matrix = N.mat[,,num.calcs]
stable.size.dist = rowSums(N.mat[,,num.calcs]) / sum(rowSums(N.mat[,,num.calcs]))
stable.age.dist = colSums(N.mat[,,num.calcs]) / sum(colSums(N.mat[,,num.calcs]))	

# save linear data
if(feedbacks=='none') {
	lambda.linear = N.sum[num.calcs]
	stable.dist.matrix.linear = stable.dist.matrix
	Z = matrix(data=NA, nrow=(n+1), ncol=max.age)
	for(a in 1:max.age) {
		Z[,a] = H.new.base[,a] * stable.dist.matrix.linear[,a]
	}  
	Z.hat = sum(Z)
	save(Z.hat, H.new.base, params, params.shell, stable.dist.matrix.linear, lambda.linear, file=paste('WB_eqVals_Linear.Rdata', sep=''))
}

