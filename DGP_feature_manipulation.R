library(MASS)
set.seed(235098)
N = 1000


X <- function(afBeta_K, afTheta_NK, Costs) {
  afBeta_K <- as.vector(afBeta_K)
  if (is.null(Costs)) {
    return_val = afTheta_NK
  } else {
    Cinv = getCinv(Costs)
    
    afCostGammaInv_N = getCostGammaInv(Costs)
    
    tiled_CostsMultipliedByBeta = matrix( t(Cinv %*% afBeta_K), nrow=nrow(afTheta_NK), ncol=length(afBeta_K), byrow=TRUE) #recycle this by row]
    return_val = afTheta_NK + ((1/afCostGammaInv_N) * tiled_CostsMultipliedByBeta)
  }
  
  return_val
  
}

getCinv <- function(Costs) {
  invertExceptConstant( Costs[["C"]] )
}


getCostGammaInv <- function(Costs) {
  Costs[["afCostGammaInv"]]
}

invertExceptConstant <- function(afCosts_KK) {
  
  stopifnot( NROW(afCosts_KK) == NCOL(afCosts_KK) )
  K = NROW(afCosts_KK)
  if ((sum(is.finite(afCosts_KK[0,0:K]))==0) & (sum(is.finite(afCosts_KK[0:K,0])) ==0)) {
    afCosts_Kminus1Kminus1 = afCosts_KK[2:K, 2:K]
    afCosts_Kminus1Kminus1inv = solve(afCosts_Kminus1Kminus1)
    
    return_val = rbind( matrix(rep(0,K),nrow=1,ncol=K), cbind(matrix(rep(0,K-1),nrow=K-1,ncol=1), afCosts_Kminus1Kminus1inv))
  }  else {
    return_val = solve( afCosts_KK )
  }
  return_val
}

# STREAM A NEW INDIVIDUAL EACH FUNCTION CALL
lossFreshIndiv <- function(afBeta_K, bIncludeManipulation=TRUE, N=1) {
  
  K = 4
  
  # PAYOUT FUNCTION
  afBeta_TRUE_K = c(0.2, 3, 0.1, 0.1)
  
  # Amount of idiosyncratic error in desired payouts
  fSigma = 0.5
  
  # UNDERLYING TYPES
  afTheta_Covariance_Kminus1Kminus1 = matrix(c(1, 1, 0.1, 1, 2, 1, 0.1, 1, 1),nrow=3,ncol=3)
  if (N==1) {
    afTheta_NK = matrix( c( 1, mvrnorm(n=N, mu=rep(0, K-1), Sigma=afTheta_Covariance_Kminus1Kminus1)), nrow=N, ncol=K)
  } else {
    afTheta_NK = cbind( rep(1, N), mvrnorm(n=N, mu=rep(0, K-1), Sigma=afTheta_Covariance_Kminus1Kminus1)) 
  }
  
  # COST FUNCTION
  bHeterogeneousCostsCorrelatedWithObservable = TRUE
  
  # cross costs
  afCosts_KK    = matrix(c(Inf, Inf, Inf, Inf,
                           Inf, 1, 0.1, 0.2,
                           Inf, 0.1, 2, 0.8,
                           Inf, 0.2, 0.8, 4),nrow=4,ncol=4)
  
  # individual multiplier
  afCostGammaInv_N = rep(1,N)
  if (bHeterogeneousCostsCorrelatedWithObservable) {
    afCostGammaInv_N[ afTheta_NK[ 2]  > 0.2] = 0.1
  }
  Costs = list(C=afCosts_KK, afCostGammaInv=afCostGammaInv_N)
  
  afEpsilon_N_fresh <- rnorm(N,0,fSigma)
  
  Y <- afTheta_NK %*% afBeta_TRUE_K + afEpsilon_N_fresh
  
  if (bIncludeManipulation) {
    afDeviation <- Y - X(afBeta_K, afTheta_NK, Costs) %*% afBeta_K
  } else {
    afDeviation <- Y - afTheta_NK %*% afBeta_K
  }
  
  (t(afDeviation) %*% afDeviation)/nrow(afDeviation)
}

# auxiliary 
lossMimicPanel <- function(afBeta_K, N, nSeed=12407, bIncludeManipulation=TRUE) {
  set.seed( nSeed )
  mean( replicate(N, lossFreshIndiv(afBeta_K, bIncludeManipulation=bIncludeManipulation)) )
}


