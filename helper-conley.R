## -----------------------------------------------------------------------------
## helper-conley.R
## Kyle Butts, CU Boulder Economics 
## 
## This creates R functions for Conley (1999) Spatial-HAC standard errors. 
## For panel or cross-sectional data. For cross-sectional just don't include time and id vectors 
## 
## -----------------------------------------------------------------------------

library(Rcpp)
library(RcppArmadillo)

# setwd("")
Rcpp::sourceCpp("helper-conley.cpp")


#' Conley Spatial-HAC Variance-Covariance Matrix
#' 
#' Formula from: Conley (1999) and Robust Inference with Multi-way Clustering
#' Original code from: https://github.com/darinchristensen/conley-se
#' 
#' @param X - nxk matrix of covariates (don't include lat and long if not in regression; if using fixed effects, they should not be included and X should be centered. See the `demeaned` option in `fixest::feols` if you want centered Xs)
#' @param e - vector of residuals
#' @param coords - nx2 matrix of lat and long (order doesn't matter)
#' @param dist_cutoff - Cutoff distance in km for spatial correlation
#' @param lag_cutoff - Cutoff time in units of your time variable for serial correlation
#' @param cores - multicore support with parallel::mclapply
#' @param verbose - if TRUE, prints out a lot of details
#' @param id - nx1 vector of id variable
#' @param time - nx1 vector of time variable (can be in any units)
#'
#' @note could be speed up if you put the lapply within C++ code. I calculate the N^2 distances for each time period
conley_ses <- function(X, e, coords, dist_cutoff, lag_cutoff, cores = 1, verbose = FALSE, id = NULL, time = NULL) {
	
	# Load multicore
	if(cores > 1) library(parallel)
	
	n <- nrow(X)
	k <- ncol(X)
	
	panel <- !(is.null(id) | is.null(time))
	if(!panel) print("No Panel dimension provided. Using just Spatial adjustment")
	
	# Spatial Correlation ------------------------------------------------------
	
	if(panel) {
		# Unique time indices
		time_unique <- unique(time)
		time_n <- length(time_unique)
		
		# For each time t, get XeeXh for that time 
		if(cores == 1) {
			XeeXhs <- lapply(time_unique, function(t) {
				sub_X <- X[time == t, ]
				sub_e <- e[time == t]
				sub_coords <- coords[time == t, ]
				n1 <- nrow(sub_X)
				
				XeeX_spatial(sub_coords, dist_cutoff, sub_X, sub_e, n1, k)
			})
		} else {
			XeeXhs <- mclapply(time_unique, function(t) {
				sub_X <- X[time == t, ]
				sub_e <- e[time == t]
				sub_coords <- coords[time == t, ]
				n1 <- nrow(sub_X)
				
				XeeX_spatial(sub_coords, dist_cutoff, sub_X, sub_e, n1, k)
			}, mc.cores = cores)
		}
		
		# Reduce across time
		XeeX_spatial <- Reduce("+",  XeeXhs)
	} else {
		XeeX_spatial <- XeeX_spatial(coords, dist_cutoff, X, e, n, k)
	}
	
	
	
	
	# Serial Correlation -------------------------------------------------------
	
	if(panel) XeeX_serial <- XeeX_serial(time, id, lag_cutoff, X, e, k)
	
	
	# Cov-Var Sandwich ---------------------------------------------------------
	
	# Bread of Sandwich
	invXX <- solve(t(X) %*% X)
	
	# Sandwiches w/ and w/o autocorrelation
	V_spatial <- invXX %*% (XeeX_spatial) %*% t(invXX)
	if(panel) V_spatial_HAC  <- invXX %*% (XeeX_spatial + XeeX_serial) %*% t(invXX)
	
	
	# Return Variance-Covariance Matrix
	if(panel) {
		return(list(
			"Spatial" = V_spatial,
			"Spatial_HAC" = V_spatial_HAC
		))
	} else {
		return(list(
			"Spatial" = V_spatial
		))
	}
} 
