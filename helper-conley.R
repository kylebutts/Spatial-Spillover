## -----------------------------------------------------------------------------
## helper-conley.R
## Kyle Butts, CU Boulder Economics 
## 
## This creates R functions for Conley (1999) Spatial-HAC standard errors. 
## For panel or cross-sectional data. For cross-sectional just don't include time and id vectors 
## -----------------------------------------------------------------------------

library(Rcpp)
library(RcppArmadillo)

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
#' @param lag_cutoff - (optional if panel) Cutoff time in units of your time variable for serial correlation
#' @param id - (optional if panel) nx1 vector of id variable
#' @param time - (optional if panel) nx1 vector of time variable (can be in any units)
#'
#' @note could be speed up if you put the lapply within C++ code. I calculate the N^2 distances for each time period
conley_ses <- function(X, e, coords, dist_cutoff, lag_cutoff = 0, id = NULL, time = NULL) {
	
	X <- as.matrix(X)
	coords <- as.matrix(coords)
	e <- as.matrix(e)
	
	n <- nrow(X)
	k <- ncol(X)
	
	panel <- !(is.null(id) | is.null(time))
	if(!panel) print("No Panel dimension provided. Using just Spatial adjustment")
	if(panel) {
		time = as.matrix(time)
		id = as.matrix(id)
	}
	
	# Spatial Correlation ------------------------------------------------------
	
	if(panel) {
		# Unique time indices
		time_unique <- unique(time)
		
		# For each time t, get XeeX_spatial for that time 
		XeeXhs <- lapply(time_unique, function(t) {
			sub_X <- as.matrix(X[time == t, ])
			sub_e <- as.matrix(e[time == t])
			sub_coords <- as.matrix(coords[time == t, ])
			n1 <- nrow(sub_X)
			
			XeeX_spatial(coords=sub_coords, cutoff=dist_cutoff, X=sub_X, e=sub_e, k=k)
		})
		
		# Reduce across time
		XeeXh_spatial <- Reduce("+",  XeeXhs)
	} else {
		XeeXh_spatial <- XeeX_spatial(coords=coords, cutoff=dist_cutoff, X=X, e=e, k=k)
	}
	
	
	
	
	# Serial Correlation -------------------------------------------------------
	
	if(panel) XeeXh_serial <- XeeX_serial(time=time, id=id, cutoff=lag_cutoff, X=X, e=e, k=k)
	
	
	# Heteroskedastic Robust ---------------------------------------------------
	
	XeeXh_robust <- XeeX_robust(X, e, k)
	
	
	
	# Cov-Var Sandwich ---------------------------------------------------------
	# Bread of Sandwich
	invXX <- solve(t(X) %*% X)
	
	# Stata , robust SEs
	V_robust <- n/(n-k) * (invXX %*% XeeXh_robust %*% invXX)
	
	# Spatial SEs
	V_spatial <- invXX %*% (XeeXh_spatial) %*% invXX
	
	# Spatial HAC SEs
	if(panel) V_spatial_HAC  <- invXX %*% (XeeXh_spatial + XeeXh_serial - XeeXh_robust) %*% invXX
	
	
	# Return Variance-Covariance Matrix
	if(panel) {
		return(list(
			"Robust" = V_robust,
			"Spatial" = V_spatial,
			"Spatial_HAC" = V_spatial_HAC
		))
	} else {
		return(list(
			"Robust" = V_robust,
			"Spatial" = V_spatial
		))
	}
} 
