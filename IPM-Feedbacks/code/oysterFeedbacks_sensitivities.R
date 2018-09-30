
## Jacob L Moore (jlmoor@ucdavis.edu)
## 07/01/17

## Calculations sensitivity of the threshold surface to delta, rho, and alpha.

## Requires the following files:
	# growthSurvData_WB.csv
	# WB_eqVals_Linear.Rdata
	# WB_sizeDists.Rdata
	# oysterFeedbacks_demography.R
	# oysterFeedbacks_estimateParams.R

## Creates the following files:
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

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Specify miscellaneous model parameters
# -------------------------------------------------------------------
# -------------------------------------------------------------------
siteName = 'WB'
feedbacks = 'positive'

## Set minimum/maximum size and maximum age
min.size = 0
max.size = 250
max.age = 10
num.calcs = 150

# Set up discetization for size (this evaluates the integral using the midpoint rule)
n = max.size
xs = seq(min.size, max.size, length=n+1)
deltax = xs[2]-xs[1] # dx increment in sizes
xs.mid.temp = .5 * (xs[-1] + xs[-(n+1)]) # midpoints
xs.mid =  c(xs.mid.temp, xs.mid.temp[n]+deltax)	# add extra point on end for evicted class

rhoVals = c('low', 'estimate')
deltaVals = 'mid' #c('mid', 'low', 'high')
alphaVals = 'estimate' #c('estimate', 'low', 'high')

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Load data and source files
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# load growth and survival data
d = read.csv(paste(data.folder, '/growthSurvData_WB.csv', sep=''), header=TRUE, sep=',')
d$age = d$timeStep+1

# remove oysters that die for growth fitting
d2 = d[-which(is.na(d$sizeNext)),]

# remove entries where no change to avoid error about non-finite initial values
if(length(which(d2$sizeChange==0)!=0)) {
	d2 = d2[-which(d2$sizeChange==0),]
}

source('oysterFeedbacks_demography.R')

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Calculate sensitivites
# -------------------------------------------------------------------
# -------------------------------------------------------------------

for(rhoVal in rhoVals) {
	for(deltaVal in deltaVals) {
		for(alphaVal in alphaVals) {

			source('oysterFeedbacks_estimateParams.R')

			# Set up individual matrices for survival, growth, and fecundity	
			S = array(0, dim=c((n+1),max.age))
			G = array(0, dim=c((n+1),(n+1),max.age))
			P = array(0, dim=c((n+1),(n+1),max.age))
			F = array(0, dim=c((n+1),(n+1),max.age))
			H.new.base = S
			
			# build kernels, and fecundity kernel when currentShell=H.eqs
			for(a in 1:max.age) {
				# build growth kernel
				G.temp = deltax*outer(xs.mid.temp,xs.mid.temp,g.yx, currentAge=a, params=params)	
				prob.growth.evict = 1-colSums(G.temp)
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
				F.temp = deltax*outer(xs.mid.temp,xs.mid.temp,f.yx, currentAge=a, params=params, params.shell=params.shell, feedbacks, currentPop=NA, currentShell=H.eqs)
				F.temp = rbind(F.temp, F.temp[n,])
				F.temp = cbind(F.temp, F.temp[,n])
				F[,,a]=F.temp		 	
			}
	
	
		# ------------- Build and evaluate sensitivity terms
			# equilbrium population
				dist.mat = eq.Mat[,,1] / sum(eq.Mat[,,1])
				eqN.vect = as.vector(N.eqs*dist.mat)	# length=2510
				
			# assemble A matrix
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
						for(vv in 2:num.zeros.pre) temp.col=rbind(temp.col,zMat)
					}
					temp.col = rbind(temp.col, P[,,col.loc])
					if(num.zeros.post>0) {
						for(ww in 1:num.zeros.post) temp.col=rbind(temp.col,zMat)
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
				A = P.mat + F.mat%*%P.mat	
				A.I = diag(dim(A)[1])
				
			# derivatives of f
				df.dn = as.vector(H.new.base)
				df.dH = exp(-params.shell$delta)
				df.ddelta = -params.shell$delta*H.eqs*exp(-params.shell$delta)
				df.drho = 0
				df.dalpha = 0
				
			# derivatives of A	
				# find "W" term
				W.mat.temp = F.mat*(params.shell$alpha + H.eqs) / (params.shell$rho*H.eqs)
				W.mat = W.mat.temp %*% P.mat	
					
				dA.dH = (params.shell$rho * params.shell$alpha * W.mat) / (params.shell$alpha + H.eqs)^2 
				dA.dH.N = dA.dH %*% eqN.vect
				
				dA.drho = (H.eqs * W.mat) / (params.shell$alpha + H.eqs) 
				dA.drho.N = dA.drho %*% eqN.vect
				
				dA.dalpha = (-params.shell$rho * H.eqs * W.mat) / (params.shell$alpha + H.eqs)^2 
				dA.dalpha.N = dA.dalpha %*% eqN.vect
				
				dA.ddelta = matrix(data=0, nrow=length(eqN.vect), ncol=length(eqN.vect))
				dA.ddelta.N = dA.ddelta %*% eqN.vect
				
				dF.dx.temp1 = cbind(A-A.I, dA.dH.N)
				dF.dx.temp2 = c(df.dn, df.dH-1) 
				dF.dx = rbind(dF.dx.temp1, dF.dx.temp2)
				dF.dx.inv = solve(dF.dx)
				
			# sensitivity to alpha
				dF.dalpha = rbind(dA.dalpha.N, df.dalpha)		
				alphaS = dF.dx.inv %*% -dF.dalpha
				N.alphaS = matrix(alphaS, nrow=(n+1), ncol=max.age)
				
			# sensitivity to rho
				dF.drho = rbind(dA.drho.N, df.drho)	
				rhoS = dF.dx.inv %*% -dF.drho
				N.rhoS = matrix(rhoS, nrow=(n+1), ncol=max.age)
				
			# sensitivity to delta
				dF.ddelta = rbind(dA.ddelta.N, df.ddelta)	
				deltaS = dF.dx.inv %*% -dF.ddelta
				N.deltaS = matrix(deltaS, nrow=(n+1), ncol=max.age)
	 
		# ------------- Elasticity terms		
			deltaE = deltaS * params.shell$delta / c(eqN.vect, H.eqs)
			alphaE = alphaS * params.shell$alpha / c(eqN.vect, H.eqs)
			rhoE = rhoS * params.shell$rho / c(eqN.vect, H.eqs)
			


	# -------------------------------------------------------------------
	# -------------------------------------------------------------------
	# Determine number of oysters of each age to add
	# -------------------------------------------------------------------
	# -------------------------------------------------------------------

	# ------------- load size distributions of each age
		load('~/Documents/UC_Davis/Schreiber_Lab/OysterProject/IPM_Feedbacks/data/WB_sizeDists.Rdata')

	# ------------- Linearization around the equilibrium
		# Jacobian
		J.temp1 = cbind(A, dA.dH.N)
		J.temp2 = c(df.dn, df.dH) 
		J = rbind(J.temp1, J.temp2)
		J.trans = t(J)
			
		# eigenvectors
		J.eigenvR = Re(eigen(J)$vectors[,1])
		J.eigenvL = Re(eigen(J.trans)$vectors[,1])	
		J.eigenvL.scaled = J.eigenvL / sum(J.eigenvL)
		J.eigenvL.mat.scaled = matrix(J.eigenvL.scaled, nrow=(n+1), ncol=max.age)


	# ------------- Geometry
		w = abs(J.eigenvL)
		w.length = sqrt(sum(w*w))
		
		xhat = c(as.vector(stable.dist.matrix.linear*N.eqs), H.eqs)
		xhat.length = sqrt(sum(xhat*xhat))
		
		w.xhat = sum(w*xhat)
		
		angleD = acos(w.xhat / (w.length * xhat.length))
		angleB = pi/2 - angleD
		
		v.all = c(as.vector(N.mat.P), 0) 
		oysters.req = rep(NA, max.age)

		# iterate through each age
		for(aa in 1:max.age) {
			start.idx = 1 + (aa-1)*251
			end.idx = 251*aa
		
			v.age = v.all*0
			v.age[start.idx:end.idx] = v.all[start.idx:end.idx]
			v.age.length = sqrt(sum(v.age*v.age))
			
			v.w = sum(v.age*w) 
		
			angleE = acos(v.w / (v.age.length*w.length))
			angleC = pi/2 - angleE
			
			b = xhat.length * sin(angleB) / sin(angleC)
			oysters.req[aa] = b / v.age.length		
		}

		# save data so no need to re-run	
		filename = paste(siteName, '_sensitivities', '_rhoVal', rhoVal, '_alphaVal', alphaVal, '_deltaVal', deltaVal, sep='')
		save(J.eigenvR, J.eigenvL, J.eigenvL.scaled, J.eigenvL.mat.scaled, J, deltaE, alphaE, rhoE, deltaS, alphaS, rhoS, N.deltaS, N.alphaS, N.rhoS, oysters.req, file=paste(filename, '.Rdata', sep='') )
			
 		}
	}
}
	


