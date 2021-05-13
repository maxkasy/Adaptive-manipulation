library(memoise) # package to memoze functions - for posterior calculations

# X1 (n x d) input matrix
# X2 (m x d) input matrix
# theta (d+2 x 1) vector of hyperparameters. Lengthscale for each dimension d
# and prior variance sigma_0^2 and noise variance sigma_n^2.
BuildCovarianceMatrix = function(X1, X2, theta){
  X1 <- as.matrix(X1) 
  X2 <- as.matrix(X2)
  n <- dim(X1)[1]
  d <- dim(X1)[2]
  m <- dim(X2)[1]
  K <- matrix(0, n, m)
  if(d==1){
    Lambda_inv <- theta[1:d]^(-2) 
  } else {
    Lambda_inv <- diag(theta[1:d]^(-2))  
  }
  for(i in 1:n){
    for(j in 1:m){
      x_diff <- matrix(X1[i,]-X2[j,],1,d)
      r2 <- x_diff %*% Lambda_inv %*% t(x_diff)
      K[i,j] <- theta[d+1]*exp(-0.5*r2)
    }
  }
  return(K)
}


# X (n x d) points sampled
# theta (d+2 x 1) vector of hyperparameters
posterior_prep = function(X, Y, theta){
  X <- as.matrix(X)
  C <- BuildCovarianceMatrix(X, X, theta)
  K = C + diag(theta[dim(X)[2]+2], nrow(X))
  cho = chol(K)
  Kinv <- chol2inv(cho)
  logdetK = 2*sum(log(diag(cho)))
  KinvY <- Kinv%*%Y
  list(C = C, K=K, Kinv=Kinv, KinvY=KinvY, logdetK = logdetK)
}

# memoize this function, to avoid re-computation
memo_posterior_prep = memoise(posterior_prep)

# x (q x d) points to be predicted 
# X (n x d) points sampled
# Y (n x 1) value at points sampled
# theta (d+2 x 1) vector of hyperparameters
PosteriorMean = function(x, X, Y, theta){
  p_prep = memo_posterior_prep(X, Y, theta)
  Cs <- BuildCovarianceMatrix(x, X, theta)
  mu <- Cs%*%p_prep$KinvY
  return(mu)
}

# X (n x d) points sampled
# Y (n x 1) value at points sampled
# theta (d+2 x 1) vector of hyperparameters
marginal_llh = function(X, Y, theta){
  p_prep = memo_posterior_prep(X, Y, theta)
  
  ll <- -0.5*(t(Y)%*%p_prep$KinvY + p_prep$logdetK)
  return(ll)
}




# Computes gradient of ll wrt hyperparameters 
marginal_llh_grad = function(X, Y, theta){
  d <- dim(X)[2]
  n <- dim(X)[1]
  p_prep = memo_posterior_prep(X, Y, theta)
  C= p_prep$C
  K = p_prep$K
  Kinv = p_prep$Kinv
  KinvY = p_prep$KinvY
  
  # dK_s/dtheta_i 
  dK_sdl <- rep(list(matrix(0, n, n)),d)
  for(k in 1:d){
    for(i in 1:n){
      for(j in 1:n){
        dK_sdl[[k]][i,j] <- C[i,j]*(X[i,k]-X[j,k])^2/(theta[k]^3)
      }
    } 
  }
  dK_sdsigma2_f <- K/theta[d+1]
  dK_sdsigma2_n <- diag(n)
  
  # KinvY KinvY'-Kinv 
  KinvY_Kinv <- KinvY%*%t(KinvY)-Kinv 
  
  # dll/dtheta_i 
  grad_ll <- matrix(0, d+2, 1)
  for(k in 1:d){
    grad_ll[k] <- 0.5*sum(diag(KinvY_Kinv%*%dK_sdl[[k]]))
  }
  grad_ll[d+1] <- 0.5*sum(diag(KinvY_Kinv%*%dK_sdsigma2_f))
  grad_ll[d+2] <- 0.5*sum(diag(KinvY_Kinv%*%dK_sdsigma2_n))
  return(grad_ll)
}



# X (n x d) points sampled
# Y (n x 1) value at points sampled
# M (1 x 1) # of restarts
maximize_marginal_LLH_update = function(X, Y, theta){
  d <- dim(X)[2]
  theta_min <- rep(1e-3, d+2)
  theta_max <- rbind(matrix(rep(1e2, d),d,1), 
                     matrix(var(Y),2,1))
  
  res <- optim(par = theta,
               fn = function(theta) -marginal_llh(X, Y, theta),
               gr = function(theta) -marginal_llh_grad(X, Y, theta), 
               lower = theta_min,
               upper = theta_max,
               method = "L-BFGS-B",
               control = list(maxit =10, trace=0)) # MAX: do only 10 update steps each time
  
  theta_opt <- res$par
  return(theta_opt)
} 

# X (n x d) points sampled
# Y (n x 1) value at points sampled
# theta (d+2 x 1) vector of hyperparameters
# xmin (1 x d) min value x is allowed to obtain
# xmax (1 x d) max value x is allowed to obtain
# nFeatures (1 x 1) # of features in analytic approximation
ThompsonSpectralSampling = function(X, Y, theta, xmin, xmax, nFeatures, M){
  n <- dim(X)[1]
  d <- dim(X)[2]
  X <- t(X)
  xmin <- t(xmin)
  xmax <- t(xmax)
  
  W <- matrix(rnorm(nFeatures*d), nFeatures, d)*rep(theta[1:d,]^(-1), nFeatures)
  b <- 2*pi*runif(nFeatures)
  noise <- rnorm(nFeatures)
  
  Phi_t <- sqrt(2*theta[d+1]/nFeatures)*cos(W%*%X+matrix(rep(b,n), nFeatures, n)) # Phi_t is t(Phi) i.e. transpose of Phi 
  A <- t(Phi_t)%*%Phi_t+diag(n)*theta[d+2]
  Ainv <- chol2inv(chol(A))
  z <- Phi_t%*%Y/theta[d+2]
  m <- z - Phi_t%*%(Ainv%*%(t(Phi_t)%*%z))
  eig <- eigen(A)
  R <- (sqrt(eig$values)*(sqrt(eig$values)+sqrt(theta[d+2])))^(-1)
  omega <- noise - (Phi_t%*%(eig$vectors%*%(R*(t(eig$vectors)%*%(t(Phi_t)%*%noise)))))+m
  
  obj = function(x){
    t(omega)%*%(sqrt(2*theta[d+1]/nFeatures)*cos(W%*%x+matrix(b, nFeatures, 1)))
  }
  gradobj = function(x){
    -t(omega)%*%(sqrt(2*theta[d+1]/nFeatures)*(matrix(rep(sin(W%*%x+b),d), nFeatures, d)*W))
  }
  
  
# MAX: again, optionally provide starting value rather than number of restarts  
  x_0 <- matrix(runif(d*M, xmin[1], xmax[1]),d,M)
  x_list <- matrix(0, M, d)
  y_list <- matrix(0, M, 1)
  for(i in 1:(M)){
    res <- optim(par = matrix(x_0[,i], d, 1),
                 fn = obj,
                 gr = gradobj,
                 lower = xmin,
                 upper = xmax,
                 method = "L-BFGS-B")
    x_list[i,] <- res$par
    y_list[i] <- res$value
  }
  x_next <- matrix(x_list[which.min(y_list),],1,d)
  y_best <- min(y_list)
  list(x_next = x_next,
       y_best = y_best)
}




# X (n x d) points sampled
# Y (n x 1) value at points sampled
# theta (d+2 x 1) vector of hyperparameters
# xmin (1 x d) lower bound for x
# xmax (1 x d) upper bound for x
# M (1 x 1) # of restarts
MinimizePosteriorMean = function(X, Y, theta, xmin, xmax, M){
  d <- dim(X)[2]
  obj_list <- matrix(0, M, 1)
  x_list <- matrix(0, M, d)
  x_0 <- matrix(runif(d*M, xmin, xmax), M, d)
  for(i in 1:M){
    res <- optim(par = x_0[i,], 
                 fn = function(x) PosteriorMean(t(x), X, Y, theta),
                 method = "L-BFGS-B",
                 upper = xmax,
                 lower = xmin, 
                 control = list(trace = 0))
    obj_list[i] <- res$value
    x_list[i,] <- res$par
  }
  x_best <- x_list[which.min(obj_list),]
  y_best <- min(obj_list)
  list(x_best = x_best,
       y_best = y_best)
} 

# FUN objective function 
# X (n x d) points sampled
# Y (n x 1) value at points sampled
# theta (d+2 x 1) vector of hyperparameters
# xmin (1 x d) min value x is allowed to obtain
# xmax (1 x d) max value x is allowed to obtain
# nFeatures # of features 
# nEval (1 x 1) # of evaluations 
BayesianOptimization = function(FUN, nEval, X, Y, theta, xmin, xmax, nFeatures, M){
  theta_opt <- theta
  for(i in 1:nEval){
    ss <- ThompsonSpectralSampling(X, Y, theta_opt, xmin, xmax, nFeatures, M)
    x_next <- ss$x_next
    y_next <- as.matrix(FUN(as.vector(x_next)))
    X <- rbind(X, x_next) 
    Y <- rbind(Y, y_next)
    
    if(i%%10==0){ # update hyperparameters every 10th iteration 
        theta_opt <- maximize_marginal_LLH_update(X, Y, theta_opt)
    }
    cat(paste("Iteration:", i, "\n"))
    cat(paste("Minimum:", min(Y), "\n"))
    cat(paste("Y:", y_next, "\n"))
    cat(paste("X", "current:", x_next,"\n\n"))
  }
  min_posterior_mean <- MinimizePosteriorMean(X, Y, theta_opt, xmin, xmax, M)
  X_opt <- min_posterior_mean$x_best
  muhat_opt <- min_posterior_mean$y_best
  list(X = X,
       Y = Y,
       X_opt = X_opt,
       theta_opt = theta_opt,
       muhat_opt = muhat_opt)
} 


