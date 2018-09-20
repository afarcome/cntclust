library(ContaminatedMixt) 
library(mixture)
library(mclust)
library(gdata)   # used for the upper triangular matrix
library(Matrix)  # forceSymmetric()
library(mnormt)  # pd.solve
library(mvnmle)  # Maximum likelihood for the multivariate normal distribution
library(mvtnorm)                                        
library(tclust)  # used for initialization
library(gtools)  # permutations

########################################
## Density of the Contaminated Normal ##
########################################

dmcn <- function(x, mu, Sigma, alpha = 0.01, eta = 1.01){
  
  alpha*dmvnorm(x=x,mean=mu,sigma=eta*Sigma)+(1-alpha)*dmvnorm(x=x,mean=mu,sigma=Sigma)
  
}

##################################################
## Random generation from a Contaminated Normal ##
##################################################

rmcn <- function(n, mu = rep(0,d), Sigma, alpha = 0.01, eta = 1.01){
  
  if(missing(Sigma))
    stop("Sigma is missing")
  if(alpha<0 | alpha>1)
    stop("alpha must be in (0,1)")
  if(eta<1)
    stop("eta must be greater than 1")
  
  d <- if(is.matrix(Sigma)) 
    ncol(Sigma)
  else 1
  
  X   <- array(0,c(n,d),dimnames=list(1:n,paste("X.",1:d,sep="")))
  bad <- rbinom(n=n,size=1,prob=alpha)
  for(i in 1:n){
    if(bad[i]==1)
      X[i,] <- rmnorm(n = 1, mean = mu, varcov = eta*Sigma)
    else
      X[i,] <- rmnorm(n = 1, mean = mu, varcov = Sigma)
  }
  
  return(X)  
  
}

###############################################################################################
## Transform a (d x d x k) array of diagonal matrices in a (d x k) matrix of diagonal values ##
###############################################################################################

Lambda.Flat <- function(Lambda){
  
  apply(Lambda,3,diag)
  
}

# Eq. (3.2) in Fritz, Garcia-Escudero & Mayo-Iscar (2011)

#########################################################################
## Basic Primal Heuristic method: feasible solution                    ##
## Used to guarantee that each cluster has, at least, d+1 observations ##
## Gallegos & Ritter (2010) CSDA, Appendix B, page 652                 ##
#########################################################################

BPH <- function(bestD,dens,z,group,no.trim,k,size.min){
  
  # bestD: Mahalanobis distances of the untrimmed observations
  # dens: (n x k)-matrix with the numerators of the elements of z 
  # z: matrix of posterior probabilities
  # group: classification vector
  # no.trim: # number of observations which are considered as not outlying
  # k: number of clusters
  # size.min: minimum size for a cluster
  
  size    <- colSums(z)
  bad.vec <- which(size<size.min)
  bad.num <- length(bad.vec)
  
  # Step 4.alpha
  
  cont1 <- sum(size[bad.vec])
  cont2 <- no.trim
  
  while(cont1 < bad.num*size.min){
    
    pos <- bestD[cont2] # position of the observation to be reassigned
    
    if(group[pos] %in% bad.vec){
      
      cont2 <- cont2 - 1
      
    } else{
      
      ass.pos <- which.max(dens[pos,bad.vec])
      group[pos] <- bad.vec[ass.pos] # new group of membership for that observation
      z <- unmap(group,1:k)
      
      # ricalcoliamo
      
      size <- colSums(z)
      bad.vec <- which(size < size.min)
      bad.num <- length(bad.vec)
      
      cont1 <- sum(size[bad.vec])
      cont2 <- cont2 - 1
      
    }
    
  }
  
  return(list(
    z=z, # matrix of posterior probabilities 
    group=group # classification vector
    )
    )
  
}

###################################
## Penalized estimation of alpha ##
###################################

penalized.alpha <- function(z,v,alpha.max=0.5,alpha.init=NULL,lambda=0,penalty="LASSO",method="BFGS",maxit=1000,reltol=1e-15,trace=0){
  
  # z: matrix of posterior probabilities
  # v: matrix of posterior probabilities to be good or bad
  # alpha.max: maximum proportion of mild outliers in each cluster
  # alpha.init: initial values for the proportions of mild outliers in each cluster
  # lambda: positive penalty parameter
  # penalty: "LASSO" or "RIDGE"
  # method: the optimization method used; see optim() for details
  # maxit: maximum number of iterations; see optim() for details
  # reltol: relative convergence tolerance; see optim() for details
  # trace: Non-negative integer. If positive, tracing information on the progress of the optimization is produced; see optim() for details
  
  if(!is.matrix(z))
    z <- matrix(z, nrow = length(z), ncol = 1, dimnames = list(names(z), deparse(substitute(z))))
  if(!is.matrix(v))
    v <- matrix(v, nrow = length(v), ncol = 1, dimnames = list(names(v), deparse(substitute(v))))
  
  n <- nrow(z)
  k <- ncol(z)
  
  if(is.null(alpha.init))
    alpha.init <- runif(k,min=0,max=alpha.max)
  
  # from (0,alpha.max) to (0,1)
  
  temp <- alpha.init/alpha.max
  
  # from (0,1) to (-Inf,Inf)
  
  initial.values <- log(temp/(1-temp))  
    
  f <- function(par,z,v,alpha.max,lambda,penalty,n,k){
    
    # from (-Inf,Inf) to (0,1) 
    
    par.new <- exp(par)/(1+exp(par))
    
    # from (0,1) to (0,alpha.max) 
    
    alpha <- alpha.max*par.new
    
    if(penalty=="LASSO")
      pen <- lambda*log(n)*sum(alpha)
    if(penalty=="RIDGE")
      pen <- lambda*log(n)*sum(alpha^2)
    
    loglik <- sum(z*(v*matrix(log(alpha),nrow=n,ncol=k,byrow=TRUE)+(1-v)*matrix(log(1-alpha),nrow=n,ncol=k,byrow=TRUE)))
    
    obj <- loglik - pen
    
    return(obj)
    
  }
  
  res   <- optim(par=initial.values, fn=f, z=z, v=v, alpha.max=alpha.max, lambda=lambda, penalty=penalty, n=n, k=k, method=method, control=list(fnscale=-1,maxit=maxit,reltol=reltol,trace=trace))
  obj   <- res$value
  alpha <- res$par
  
  # from (-Inf,Inf) to (0,1) 
  
  alpha.new <- exp(alpha)/(1+exp(alpha))
  
  # from (0,1) to (alpha.min,1) 
  
  alpha <- alpha.max*alpha.new
  
  return(
    list(
      penalized.loglik = obj, # penalized log-likelihood
      alpha = alpha # proportions of mild outliers
    )
  ) 
    
}
  
#####################################################################################
## Trimmed version of Mixtures of Gaussian and contaminated Gaussian distributions 
## N.B. The trimming proportion (alpha0) is fixed   ##
#####################################################################################

cntclust <- function( 
  X,                 # data matrix
  k,                 # number of mixture components
  alpha0 = 0,        # trimming proportion
  density = "CN",    # "N" for normal or "CN" for contaminated normal
  initialization = "tclust", # or manual; in that case, a partition (init.clus) must be provided
  init.clus = NULL,  # initial partition to be used when initialization = "manual"
  penalty = NULL,    # can be: NULL (no penalization), "LASSO", or "RIDGE"
  lambda = 0,        # nonnegative penalty parameter
  Cstepz = TRUE,     # if true, then the CEM algorithm is applied; if FALSE, then the EM algorithm is considered
  restr = "eigen",   # "deter" and "sigma" do not work at the moment 
  restr.fact = 100,  # eigenvalues constraint
  alpha.max = 0.5,   # maximum proportion of mild outliers in each group for the contaminated normal distribution
  size.min = 0,      # minimum size for a cluster
  method = "BFGS",   # used when penalty is not NULL
  tol = 10^-20,      # one of the 2 stopping rules (see max.iter for the other)
  max.iter = 1000,   # maximum number of iterations
  trace = 0,         # Non-negative integer. If positive, tracing information on the progress of the optimization is produced; see optim() for details
  plotll = FALSE     # plot of the log-likelihood versus iterations
)
{
    
    eps = 0.001,       # to avoid problems with alpha
    Cstepv = FALSE,    # will be implemented in the future 

    if(is.data.frame(X))
    X <- data.matrix(X)
  if(!is.matrix(X))
    X <- matrix(X, nrow = length(X), ncol = 1, dimnames = list(names(X), deparse(substitute(X))))
  if(!is.numeric(X))
    stop("numeric matrix/vector expected for X")
  if(any(is.na(X)))  
    stop("No NAs allowed")
  if(is.null(k)) 
    stop("The number of groups k is NULL")
  if(!is.null(init.clus) & alpha0 != length(which(init.clus == 0))/nrow(X)){
    temp   <- length(which(init.clus == 0))
    alpha0 <- temp/nrow(X)
    cat("\n")
    cat(paste("alpha0 has been set to ",alpha0,", i.e. to the proportion of zeros in init.clus",sep=""))
    cat("\n")
  }
    
  #stop("alpha0 must be equal to the proportion of 0's in init.clus")
  
  # --------------------- #
  # Definition of objects #
  # --------------------- #
  
  n       <- nrow(X) # number of units
  d       <- ncol(X) # number of variables
  no.trim <- floor(n*(1-alpha0)) # number of observations which are considered as not outlying
  
  #size.min <- d + 1 # minimum size for each cluster 
  
  if(k*size.min > no.trim){
    stop("There are too many clusters w.r.t. the number of 'good' observations")
  }
  
  # E-step quantities
  
  z          <- array(0,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep="")))
  v          <- array(1/n,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep="")))
  mah        <- array(0,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep="")))
  zvgood     <- array(0,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep="")))
  zvbad      <- array(0,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep="")))
  correction <- array(1,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep=""))) 
  
  # M-step Parameters 
  
  prior  <- numeric(k) 
  mu     <- array(0,c(d,k),dimnames=list(paste("X.",1:d,sep=""),paste("comp.",1:k,sep="")))
  Sigma  <- array(0,c(d,d,k),dimnames=list(paste("X.",1:d,sep=""),paste("X.",1:d,sep=""),paste("comp.",1:k,sep="")))
  Lambda <- array(0,c(d,d,k),dimnames=list(paste("X.",1:d,sep=""),paste("X.",1:d,sep=""),paste("comp.",1:k,sep="")))
  Gamma  <- array(0,c(d,d,k),dimnames=list(paste("X.",1:d,sep=""),paste("X.",1:d,sep=""),paste("comp.",1:k,sep="")))
  if(density=="N"){
    alpha <- rep(0,k)
    eta   <- rep(1,k)
  }
  if(density=="CN"){
    alpha <- rep(1/n,k)
    eta   <- rep((n+1)/n,k)
  }
  
  # Distribution 
  
  dens  <- array(0,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep=""))) # weights*fX
  densX <- array(0,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep=""))) # fX
  
  # ------------------------------- #
  # Initialization of the algorithm #
  # ------------------------------- #
  
  if(initialization == "manual"){
    
    temp   <- length(which(init.clus == 0))
    alpha0 <- temp/n
    sizes  <- table(init.clus[which(init.clus>0)])
    prior  <- sizes/sum(sizes)
    bestD  <- which(init.clus > 0)
    z[bestD,] <- unmap(init.clus[bestD],1:k)
    for(j in 1:k){
      mu[,j] <- colSums(z[bestD,j]*correction[bestD,j]/sum(z[bestD,j]*correction[bestD,j])*X[bestD,])
      temp   <- crossprod(sqrt(z[bestD,j]*correction[bestD,j]/sum(z[bestD,j]))*(X[bestD,]-matrix(rep(mu[,j],no.trim),no.trim,d,byrow=TRUE)))
      Sigma[,,j] <- as.matrix(forceSymmetric(x=temp))
    }
    temp <- list()
    for(j in 1:k){
      if(density=="N") 
        densX[bestD,j] <- dmvnorm(x=X[bestD,],mean=mu[,j],sigma=Sigma[,,j])
      if(density=="CN") 
        densX[bestD,j] <- dmcn(x=X[bestD,],mu=mu[,j],Sigma=Sigma[,,j],alpha=alpha[j],eta=eta[j])
      
      dens[bestD,j] <- prior[j]*densX[bestD,j]   
    }    
    
    # ------------------------------------- # 
    # Global - Observed-data log-likelihood # 
    # ------------------------------------- #
    
    temp$obj <- sum(log(rowSums(dens[bestD,]))) 
    
  }
  
  # if(initialization == "gpcm"){
  #   
  #   init <- gpcm(data=X, G=k, mnames="VVV", start=0)
  #   prior <- init$gpar$pi
  #   for(j in 1:k){
  #     mu[,j] <- init$gpar[[j]]$mu
  #     Sigma[,,j] <- init$gpar[[j]]$sigma
  #   }
  #   temp <- list()
  #   temp$obj <- init$BIC[,,1]
  #   
  #   print(temp$obj)
  #   
  # }    
    
  if(initialization == "tclust"){
    
    temp   <- tclust(x = X, k = k, alpha = alpha0, restr = restr, restr.fact = restr.fact)
    prior  <- temp$weights
    mu     <- temp$centers
    Sigma  <- temp$cov
    sizes  <- temp$size
    
  }
  
  # Penalized Log-Likelihood (obj) on alpha_1,...,alpha_k
  
  if(is.null(penalty)){
    obj <- temp$obj
  }  else{
    if(penalty == "LASSO"){
      obj <- temp$obj - lambda*log(n)*sum(alpha)
    }
    if(penalty == "RIDGE"){
      obj <- temp$obj - lambda*log(n)*sum(alpha^2)
    }
  }
  
  # EM algorithm --------------------------------------------------
  
  # Preliminary definition of convergence criterions
  
  check     <- 0
  iteration <- 1
  
  while(check<1){
    
    # ++++++ #
    # E-Step #
    # ++++++ #
    
    # ---------------------------------------------------------- #
    # Step 2.1: Trimming and cluster assignments (E and C-steps) #
    # ---------------------------------------------------------- #
    
    # We work on the whole set of n observations
    
    for(j in 1:k){
      if(density=="N") 
        densX[,j] <- dmvnorm(x=X, mean=mu[,j], sigma=Sigma[,,j], log=FALSE)
      # dmnorm(x=X,mean=mu[,j],varcov=Sigma[,,j])
      if(density=="CN") 
        densX[,j] <- dmcn(x=X,mu=mu[,j],Sigma=Sigma[,,j],alpha=alpha[j],eta=eta[j])
      
      dens[,j] <- prior[j]*densX[,j]    
    } 
    
    # Define the subset of the untrimmed observations
    # Garcìa-Escudero et al. (2008): The Annals of Statistics
    
    D       <- apply(dens,1,max)          # step 2.1: Discriminant functions
    temp    <- order(D,decreasing = TRUE) # from the worst to the best
    bestD   <- temp[1:no.trim]            # the best no.trim values
    trimmed <- temp[-(1:no.trim)]
    
    # posterior probabilities on the whole set of n observations
    
    z <- dens/matrix(rep(rowSums(dens),k),ncol=k)
    
    # ++++++ #
    # C-Step # Step 2.2 in Garcìa-Escudero et al. (2008): The Annals of Statistics
    # ++++++ #
    
    if(Cstepz){
      group <- map(z)
      z     <- unmap(group,1:k)
      
      if(min(colSums(z)) < size.min){
        
        print(paste("The minimum cluster size should be ",size.min,", instead the cluster sizes are:",sep=""))
        print(colSums(z))
        
        temp <- BPH(bestD=bestD,dens=dens,z=z,group=group,no.trim=no.trim,k=k,size.min=size.min)
        group <- temp$group
        z <- temp$z
        
        print("With the basic primal heuristic solution the cluster sizes become:")
        print(colSums(z))
        cat("\n")
        
      }
      
    }
    
    if(density == "CN"){
      for(j in 1:k){ 
        #v[,j] <- (alpha[j]*dmnorm(x=X,mean=mu[,j],varcov=eta[j]*Sigma[,,j]))/dmcn(x=X,mu=mu[,j],Sigma=Sigma[,,j],alpha=alpha[j],eta=eta[j])
        v[,j] <- (alpha[j]*dmvnorm(x=X,mean=mu[,j],sigma=eta[j]*Sigma[,,j]))/dmcn(x=X,mu=mu[,j],Sigma=Sigma[,,j],alpha=alpha[j],eta=eta[j])
      }
      v[is.nan(v)] = 1 # to avoid 0/0 in v
      
      if(Cstepv){
        v <- round(v)
      }
      
      zvbad  <- z*v 
      zvgood <- z*(1-v)
      correction <- (1-v)+v*matrix(rep(1/eta),n,k,byrow=TRUE)
    }
    
    # ------------------------------------ #
    # Step 2.2: Update parameters (M-step) #
    # ------------------------------------ #
    
    # ++++++ #
    # M-Step #
    # ++++++ #
    
    # Parameters
    
    prior <- colMeans(z[bestD,])
    sizes <- colSums(z[bestD,])
    for(j in 1:k){
      
      #pippo <- colMeans(X[which(z[bestD,j] == 1),])
      mu[,j] <- colSums(z[bestD,j]*correction[bestD,j]/sum(z[bestD,j]*correction[bestD,j])*X[bestD,])
      temp   <- crossprod(sqrt(z[bestD,j]*correction[bestD,j]/sum(z[bestD,j]))*(X[bestD,]-matrix(rep(mu[,j],no.trim),no.trim,d,byrow=TRUE)))
      Sigma[,,j] <- as.matrix(forceSymmetric(x=temp))
      
    }
      
    ## Impose constraints on the eigenvalues
    
    if(restr == "eigen"){
      
      for(j in 1:k){
      
        temp        <- eigen(Sigma[,,j])
        Lambda[,,j] <- diag(temp$values)
        Gamma[,,j]  <- temp$vectors
        
      }
      LambdaFlat <- Lambda.Flat(Lambda)
      LambdaFlatnew <- tclust:::.restr2_eigenv(autovalues=LambdaFlat, ni.ini=sizes, restr.fact=restr.fact, 1e-16)
      for(j in 1:k){
        Lambda[,,j] <- diag(LambdaFlatnew[,j])
        temp <- Gamma[,,j] %*% Lambda[,,j] %*% t(Gamma[,,j])
        Sigma[,,j] <- as.matrix(forceSymmetric(x=temp))
      }
      
    }
      
    # CN
    
    if(density == "CN"){
      
      for(j in 1:k){
        
        # direct method
        
        alpha[j] <- min(alpha.max,sum(zvbad[bestD,j])/sum(z[bestD,j]))
        alpha[which(alpha==alpha.max)] = (alpha.max-eps)
        
        # indirect method
        
        # alpha[j] <- optimize(g,c(alpha.min,1),maximum = TRUE,j=j,z=z.trim,v=v.trim)$maximum
        
        aj           <- sum(zvbad[bestD,j])
        mah[bestD,j] <- mahalanobis(x=X[bestD,], center=mu[,j], cov=Sigma[,,j], inverted=FALSE)
        bj           <- sum(zvbad[bestD,j]*mah[bestD,j])
        
        eta[j]       <- max(1.001,bj/(d*aj))
        
      }
      
      if(!is.null(penalty))
        alpha <- penalized.alpha(z=z,v=v,alpha.max=alpha.max,alpha.init=alpha,lambda=lambda,penalty=penalty,method=method,maxit=max.iter,reltol=tol,trace=trace)$alpha 
        
    }
    
    # --------- #
    # densities #
    # --------- #
    
    for(j in 1:k){
      if(density=="N") 
        densX[bestD,j] <- dmvnorm(x=X[bestD,],mean=mu[,j],sigma=Sigma[,,j])
        #densX[bestD,j] <- dmnorm(x=X[bestD,],mean=mu[,j],varcov=Sigma[,,j])
      if(density=="CN") 
        densX[bestD,j] <- dmcn(x=X[bestD,],mu=mu[,j],Sigma=Sigma[,,j],alpha=alpha[j],eta=eta[j])
      
      dens[bestD,j] <- prior[j]*densX[bestD,j]   
    }    
    
    # ------------------------------------- # 
    # Global - Observed-data log-likelihood # 
    # ------------------------------------- #
    
    llvalues <- sum(log(rowSums(dens[bestD,]))) # CONTROLLARE la giusta likelihood di convergenza
    #loglik <- c(loglik,llvalues)
    
    # Penalized Log-Likelihood
    
    if(is.null(penalty)){
      obj <- c(obj,llvalues)
    } else{
      if(penalty=="LASSO"){
        obj <- c(obj,(llvalues-lambda*log(n)*sum(alpha)))
      }
      if(penalty=="RIDGE"){
        obj <- c(obj,(llvalues-lambda*log(n)*sum(alpha^2)))
      }
    }
    
    # cat("*")
    iteration <- iteration + 1
    
    # if(iteration == max.iter | (loglik[iteration]-loglik[iteration-1])<tol)
    if(iteration == max.iter | (obj[iteration]-obj[iteration-1])<tol)
      check <- 1

  }
  
  # cat("\n")
  
  # ------------------------ #
  # Penalized Log-likelihood #
  # ------------------------ #
  
  #finalloglik <- loglik[iteration] 
  finalobj <- obj[iteration] 
  
  # ---------------------------------------------------------- #
  # plot to check monotonicity of the penalized log-likelihood #
  # ---------------------------------------------------------- #
  
  if(plotll){
    
    par(mai=c(0.84,0.8,0.012,0.004))
    par(las = 3)
    par(cex.axis=0.7)
    par(cex.lab=1.2)
    plot(0:(iteration-1),obj[1:iteration],type="l",axes = FALSE,xlab="iterations",ylab="Objective function",lwd=2)
    axis(1, at = 0:(iteration-1),label = 0:(iteration-1)) 
    axis(2)
    box(col = "black")
    
  }
  
  # The EM-algorithm is finished #
  
  # --------------------- #
  # Classification Matrix #
  # --------------------- #
  
  clas.num <- rep(0,n)
  clas.chr <- rep("trimmed",n)
  clas.num[bestD] <- apply(z[bestD,],1,which.max) 
  for(j in 1:k){
    temp <- bestD[which(clas.num[bestD]==j)]
    clas.chr[temp] <- ifelse(v[temp,j]>0.5,"bad","*") 
  }
  classification <- data.frame(group=clas.num, detection=clas.chr)
  
  # -------------------- #
  # Number of parameters #
  # -------------------- #
  
  nparX <- k*(d + d*(d+1)/2)
  if(density == "CN") 
    nparX <- nparX + 2*k
  npar <- (k-1) + nparX
  
  # -------------------- #
  # Information criteria #
  # -------------------- #
  
  AIC <- 2*finalobj - 2*npar            
  BIC <- 2*finalobj - npar*log(no.trim) 
  
  result <- list(
    X              = X,
    k              = k,            
    n              = n,
    no.trim        = no.trim,
    density        = density,
    classification = classification,
    npar           = npar,
    sizes          = round(sizes,5),
    prior          = round(prior,5),
    mu             = round(mu,5),
    Sigma          = round(Sigma,5),
    alpha          = round(alpha,5),
    eta            = round(eta,5),
    z              = round(z,5),
    iter.stop      = iteration,
    obj            = round(finalobj,5),
    AIC            = round(AIC,5),
    BIC            = round(BIC,5),
    call           = match.call()    
  )
  class(result) <- "cntclust"
  
  return(result)
  
}

##################################
# Random generator from Mixtures #
# of N, t, and CN distributions  #
##################################

rM <- function(N, # number of observations
               model="N", # mixture component density ("N", "t", or "CN")
               prior, # mixture weights
               mu, # (d x k)-matrix of mean vectors
               Sigma, # 3D-array of scalar matrices
               alpha, # vector of proportion of mild outliers (when model="CN")
               eta, # vector of inflation parameters (when model="CN")
               df, # vector of degrees of freedom (when model="t")
               seed=NULL # seed for random generation
               ){
  
  set.seed(seed)
  
  P <- nrow(mu)
  K <- length(prior)
  
  X      <- array(0,c(N,P),dimnames=list(1:N,paste("X.",1:P,sep="")))
  z      <- t(rmultinom(N, 1, prob = prior))
  groups <- map(z)
  for(i in 1:N){
    if(model=="N")
      X[i,] <- rmnorm(1,mean=mu[,groups[i]],varcov=Sigma[,,groups[i]])
    if(model=="t")
      X[i,] <- rmt(1,mean=mu[,groups[i]],S=Sigma[,,groups[i]],df=df[groups[i]])
    if(model=="CN")
      X[i,] <- rmcn(1, mu=mu[,groups[i]],Sigma=Sigma[,,groups[i]],alpha=alpha[groups[i]],eta=eta[groups[i]])
  }
  
  return(
    list(
      X = X, # data matrix 
      z = z, # (n x k)-matrix of posterior probabilities
      groups = groups # classification vector
    )
  )  
  
}
