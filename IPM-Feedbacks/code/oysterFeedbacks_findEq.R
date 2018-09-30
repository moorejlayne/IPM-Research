## Jacob L Moore (jlmoor@ucdavis.edu)
## 07/01/17

# Solves for equilibrium of positive feedback model. Called by oysterFeedbacks_estimateParams.R

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Find roots when positive only feedbacks
# -------------------------------------------------------------------
# -------------------------------------------------------------------

H.eqs = (params.shell$p.hat*params.shell$alpha) / (params.shell$rho - params.shell$p.hat)
N.eqs = (1-exp(-params.shell$delta)) * H.eqs / Z.hat

# size by age matrix of equilibirum values
eq.Mat = array(NA, dim=c(max.size+1,max.age,1))
for(a in 1:max.age) {	
	eq.Mat[,a,1] = N.eqs*stable.dist.matrix.linear[,a]
}	
