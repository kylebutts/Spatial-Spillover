
# Winsorize helper function ----------------------------------------------------

# https://gist.github.com/maifeng/3e5a1f0dd65bf0fd9c2b8f4ac8f5fab3
winsorize_x = function(x, cut = 0.01){
	cut_point_top <- quantile(x, 1 - cut, na.rm = T)
	cut_point_bottom <- quantile(x, cut, na.rm = T)
	i = which(x >= cut_point_top) 
	x[i] = cut_point_top
	j = which(x <= cut_point_bottom) 
	x[j] = cut_point_bottom
	return(x)
}


# Format Regression tex cell ---------------------------------------------------
# Regression Star helper function
reg_format <- function(pt, se) {
	pt_str <- paste0("$", sprintf("%0.4f", pt))
	
	t_stat <- abs(pt/se)
	
	if(t_stat > 2.576) {
		pt_str <- paste0(pt_str, "^{***}$")
	}else if(t_stat > 1.960) {
		pt_str <- paste0(pt_str, "^{**}$")
	}else if(t_stat > 1.645) {
		pt_str <- paste0(pt_str, "^{*}$")
	}else {
		pt_str <- paste0(pt_str, "$")
	}
	
	return(pt_str)
}



# Oaxaca-Blinder Regression ----------------------------------------------------

#' @description This function returns the Oaxaca-Blinder Estimand for the Average Treatment Effect on the Treated, following Patrick Kline - Oaxaca-Blinder as a Reweighting Estimator (2011)
#'
#' @param data dataframe to estimate with
#' @param formula Formula for linear regression specification. Either formula or string object. Do not include the treatment variable in this
#' @param treat string for treatment variable, must be 1 = treat, 0 = control
#' @param cluster optional - string for cluster variable. Variable "indicating which observations belong to the same cluster". Passed along to `clubSandwich::vcovCR()`
oaxaca_estimate <- function(data, formula, treat, cluster = NA) {
	
	## Regression using Control units ------------------------------------------
	
	# Prepare formula
	if(is.character(formula)) formula <- as.formula(formula)
	Y_var <- all.vars(formula[[2]])
	X_vars <- all.vars(formula)[-1]

	
	# Subset by control and treated
	data_control <- data[data[[treat]] == 0, ]
	data_treat <- data[data[[treat]] == 1, ]
	control_reg <- lm(formula, data = data_control) 
	
	## Mean of TE --------------------------------------------------------------
	
	# Counterfactual Mean
	pred <- predict(control_reg, data)
	
	# Y_i(1) - \hat{Y}_i(0)
	diff <- data[[Y_var]] - pred
	
	# E[Y_i(1) - \hat{Y}_i(0) | D_i = 1]
	te <- mean(diff[data[[treat]] == 1], na.rm = TRUE)
	
	
	## Variance of TE ----------------------------------------------------------
	
	# Stata Robust Var-Cov matrix
	if(is.na(cluster)) {
		V0 <- sandwich::vcovHC(control_reg, type = "HC1")
	} else {
		V0 <- clubSandwich::vcovCR(control_reg, cluster = data_control[[cluster]], type = "CR1S")
	}
	
	# Var(te) comes from out of sample prediction variance formula
	X <- model.matrix(formula, data = data)
	D <- data[[treat]]
	V1 <- var(data_treat[[Y_var]]) / sum(D)
	Vdiff <- V1 + (t(D) %*% X %*% V0 %*% t(X) %*% D)/(sum(D)^2)
	
	se <- as.vector(sqrt(Vdiff))
	
	return(list("te" = te, "se" = se))	
}


