Code for the paper "Robust model based clustering with mild and gross outliers" 

This project reports R code for implementing the methods and reproducing
data analysis and simulation studies in the paper "Robust model based 
clustering with mild and gross outliers"

The following libraries shall be installed in R before using code from
this project:

library(ContaminatedMixt) 

library(mixture)

library(mclust)

library(gdata)   

library(Matrix)  

library(mnormt)  

library(mvnmle)  

library(mvtnorm)                                        

library(tclust)  

library(gtools)  

The main function for implementation is called cntclust. 
Note that the code is commented.

The cntclust function takes in input: 

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

  trace = 0,         # Non-negative integer. If positive, tracing information on
 the progress of the optimization is produced; see optim() for details

  plotll = FALSE     # plot of the log-likelihood versus iterations

and gives in output a list with elements:  

    X, k, n, mu, Sigma, alpha, eta, AIC, BIC, call              

    no.trim: n*(1-\alpha0)     

    density: "N" or "CN" as inputed         

    classification: matrix indicating the group of each observation
    and whether it belongs to the inflated CN component or not 

    npar: number of parameters            

    sizes: cluster sizes           

    prior: mixture weights        

    z:  (n x k)-matrix of posterior probabilities         

    iter.stop: number of iterations at convergence 

    obj: likelihood or penalized likelihood at convergence             

