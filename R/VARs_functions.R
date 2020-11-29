#' @export
### Root2coef
### This is function which maps roots (in L) of the characteristic function of an AR process as inout
### to the AR coefficients
###   
### input
###   
### p:   nummber of lags
### r_p  (optional) p-vector of roots outside the unit circle 
###
### output 
### a_p  p-vector of lag-coefficients
### 
### the recursion follows from the equation:    (1-a[3,1]L^1-a[3,2]L^2-a[3,3]L^3)(L-r_4) 
###                                           = (1-a[4,1]L^1-a[4,2]L^2-a[4,3]L^3-a[4,4]L^4)r_4 
###
###
###
Roots2coef = function(p,r_p) {
   if (missing(r_p)) r_p <- 0.5/(runif(p)-0.5)  #random number outside unit circle
   if (min(abs(r_p)) < 1) {
   	return("r_p is within the unit circle")
   	#print("r_p is within the unit circle") 
   }
   a = matrix(0,p,p)
   a[1,1] = 1/r_p[1]
   if ( p>1 ) {
   for (i in 2:p) {
      for ( j in 1:i ) {
         if (j == 1)           a[i,j] = a[i-1,j] + 1/r_p[i] 
         if ((j > 1) & (j<i))  a[i,j] = a[i-1,j] - a[i-1,j-1]/r_p[i]
         if (j == i)           a[i,j] = -a[i-1,i-1]/r_p[i]
      }
   }
   }   
   
   #R = matrix(0,p,p)
   #for (i in 1: p)     {
   #   for (j in 1: p ) {
   #        R[i,j] = r_p[j]^i
   #   }
   #}
   #a[p,]%*%R[,]
   return(a[p,])
}



#' This function will generate iid multivariate normal random time series.
#'
#' @param T     : length of the generated time series
#' @param sigma : (n x n) covariance matrix of the normal series
#' @return      : T x n matrix of iid normal time series
#' @export
rnormSIGMA = function(T,sigma) {
    # generate random numbers from iid multivariate normal distribution with covariance matrix Sigma 
    n = dim(sigma)[1]
    U = rnorm(T*n)
    dim(U) = c(T,n)
    U = U%*%chol(sigma)
    return(U)
}



#'  This function generates random numbers from iid multivariate conditional normal distribution with covariance matrix Sigma, given i-th component has the value of v, this will be an (n-1) dimensional random number  
#'
#' @param T     : length of the generated time series
#' @param sigma : (n x n) covariance matrix of the normal series
#' @param I     : index of the conditioning component
#' @param v     : the value of the  conditioning component
#' @switch      : switch = 1, gives the conditional random series, while switch = 0, gives the expected values. 
#' @return      : T x (n-1) matrix of iid conditional normal time series
#' @export
rnormSIGMA_cond = function(T,sigma,I,v,switch) {
    # generate random numbers from iid multivariate conditional normal distribution with covariance matrix Sigma, given 
    # i-th component has the value of v, this will be an (n-1) dimensional random number  
      sigma_cond = as.matrix(sigma[-I,-I])-sigma[-I,I]%*%solve(sigma[I,I])%*%sigma[I,-I]
      mu_cond = sigma[-I,I]%*%solve(sigma[I,I])%*%v
      U = rnormSIGMA(T,sigma_cond)*switch+as.vector(c(1:T)/c(1:T))%*%t(mu_cond)
      return(U)
}


#'  This function selects randomly N elements out of a set of T elements 
#'
#' @param N     : the number to be selected
#' @param T     : the total number of elements 
#' @export
NoutofT = function(N,T) {
  unique(round(runif(3*N)*(T-1)))[1:N]+1
}


#' This function generate impulse response functions for VAR,CIVAR,MRVAR MRCIVAR, also for GVAR, CIGVAR, MRGVAR and MRCIGVAR.
#' For the later four classes of models it also provides the functionalities to calculate the global, regional and country-specific shocks.
#' It also calculates global and regional responses.    
#' @export
irf_B_sigma = function (B, sigma, nstep, comb, irf = c("gen", "chol", "chol1","gen1","genN1", "comb1","smat"),G=NA,smat=NA) 
{
    neq <- dim(B)[1]
    nvar <- dim(B)[2]
    lags <- dim(B)[3]
    n = dim(sigma)[1]

    if (irf == "smat") {
        smat = (smat)
    }

    if (irf == "chol") 	   {        smat = chol(sigma)    					    }
    if (irf == "PTdecomp") {        smat = chol(G%*%sigma%*%t(G))          		    }
    if (irf == "gen")      {        smat = t(sigma %*% diag(1/(sqrt(diag(sigma)))))     }
    if (irf == "chol1")    {        smat = chol(sigma) %*% diag(1/diag(chol(sigma)))    }
    if (irf == "gen1")     {        smat = t(sigma %*% diag(1/(diag(sigma))))           }
    if (irf == "genN1")    {        smat = t(sigma %*% diag(-1/(diag(sigma))))          }
    if ( irf == "concerts1") { smat  = t((sigma%*%diag(1/(diag(sigma))))%*%comb)        };
    if ( irf == "concerts0") { smat  = t((sigma%*%diag(1/(sqrt(diag(sigma)))))%*%comb)  };

    if ( irf == "concertc") { 
	   c     = as.numeric(!(comb[,1]==0));
         smat  = matrix(0,n,n); 
         for (i in 1: n) { smat[i,] = sigma[i,]%*%diag(c)%*%INVI(sigma,c,i)%*%comb; }
         smat = t(smat);
    };

    if (irf == "comb") {
        DD = diag(t(comb) %*% sigma %*% comb)
        for (i in 1:length(DD)) {
            if (DD[i] > 0) 
                DD[i] = sqrt(1/DD[i])
        }
        DD = diag(DD)
        smat = t(sigma %*% comb %*% DD)
    }
    if (irf == "comb1") {
        DD = diag(t(comb) %*% sigma %*% comb)
        for (i in 1:length(DD)) {
            if (DD[i] > 0) 
                DD[i] = (1/DD[i])
        }
        DD = diag(DD)
        smat = t(sigma %*% comb %*% DD)
    }
    if (dim(smat)[2] != dim(B)[2]) 
        stop("B and smat conflict on # of variables")
    response <- array(0, dim = c(neq, nvar, nstep))
    response[, , 1] <- t(smat)
    for (it in 2:nstep) {
        for (ilag in 1:min(lags, it - 1)) response[, , it] <- response[, 
            , it] + B[, , ilag] %*% response[, , it - ilag]
    }
    dimnames(response) <- list(dimnames(B)[[2]], dimnames(smat)[[1]], 
        as.character(0:(nstep - 1)))
    return(response)
}

