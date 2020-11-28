// https://github.com/darinchristensen/conley-se

#include <iostream>
#include <math.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;


// Helper functions ------------------------------------------------------------

// Degrees to Radians ----------------------------------------------------------
double to_radians_cpp(double degrees){
	return(degrees * 3.141593 / 180);
}

// Haversine Distance ----------------------------------------------------------
double haversine_cpp(double lat1, double long1,
                     double lat2, double long2,
                     std::string unit="km"){
	int radius = 6378;
	double delta_phi = to_radians_cpp(lat2 - lat1);
	double delta_lambda = to_radians_cpp(long2 - long1);
	double phi1 = to_radians_cpp(lat1);
	double phi2 = to_radians_cpp(lat2);
	double term1 = pow(sin(delta_phi / 2), 2);
	double term2 = cos(phi1) * cos(phi2) * pow(sin(delta_lambda/2), 2);
	double the_terms = term1 + term2;
	double delta_sigma = 2 * atan2(sqrt(the_terms), sqrt(1-the_terms));
	double distance = radius * delta_sigma;
	
	/* if it is anything *but* km it is miles */
	if(unit != "km"){
		return(distance*0.621371);
	}
	
	return(distance);
}

// Distance Formula used by Solomon Hsiang -------------------------------------
double sh_cpp(double lat1, double long1,
              double lat2, double long2) {
	double distance = pow(pow(111 * (lat1 - lat2), 2) + pow(cos(lat1 * 3.1415927 / 180) * 111 * (long1 - long2), 2), .5);
	return(distance);
}



// Cross-sectional spatial clustering ------------------------------------------

//' Cross-sectional XeeX matrix 
//' Lg because it doesn't store distance matrices
//' 
//' @param M - nx2 matrix of lat-long points 
//' @param cutoff - cutoff distance in kilometers
//' @param X - nxk model matrix from regression
//' @param e - nx1 vector or residuals
//' @param n1 - 
//' @param k - number of covariates
//' @param kernel - bartlett or any string "" for 0/1 Kernel
//' @param dist_fn - Haversine or any string "" for other distance function. Use Haversine
//  [[Rcpp::export]]
arma::mat XeeX_spatial(arma::mat M, double cutoff,
                    arma::mat X, arma::vec e, int n1, int k,
                    std::string kernel="bartlett",
                    std::string dist_fn="Haversine"){
	
	long long int nrow = M.n_rows;
	
	// Initialize XeeXh matrix
	arma::mat XeeXh(k, k, fill::zeros);
	arma::mat temp;
	
	for( long long int i = 0; i < nrow; i++ ){
		for( long long int j = 0; j < nrow; j++ ){
			double d;
			
			// Distance Function
			if(dist_fn == "Haversine") {
				d = haversine_cpp(M(i,0), M(i,1), M(j,0), M(j,1));
			} else {
				d = sh_cpp(M(i,0), M(i,1), M(j,0), M(j,1));
			}
			
			if(d <= cutoff) {
				// Kernel 
				int K = (kernel == "bartlett") ? (1 - d / cutoff) : 1;
				
			
				// Iteration (i,j): XeeXh += X_i' X_j e_i e_j K(d_ij)
				// Formula from Robust Inference with Multi-way Clustering
				XeeXh += (e[i] * e[j] * K) * (X.row(i).t() * X.row(j));
			}
		}
	}
	
	// Return kxk matrix X'ee'X
	return XeeXh;
}



//' Create Variance-Covariance matrix of serial correlations
//' 
//' @param time - nx1 vector of time periods
//' @param id - nx1 vector of IDs  
//' @param cutoff - cutoff distance in unit of time 
//' @param X - nxk model matrix from regression
//' @param e - nx1 vector or residuals
//' @param k - number of covariates
// [[Rcpp::export]]
arma::mat XeeX_serial(arma::vec time, arma::vec id, double cutoff,
                   arma::mat X, arma::vec e, int k){
	
	long long int nrow = time.n_elem;
	
	arma::mat XeeXh(k, k, fill::zeros);
	
	
	for( long long int i = 0; i < nrow; i++ ){
		for( long long int j = 0; j < nrow; j++ ) {
			double t_diff = abs(time[i] - time[j]);
			bool same_id = id[i] == id[j];
			
			if(t_diff <= cutoff && same_id) {
				// Iteration (i,j): XeeXh += X_i' X_j e_i e_j K(i,j)
				// K(i,j) = 1(t_i - t_j <= cutoff) * 1(id_i == id_j)
				// Formula from Robust Inference with Multi-way Clustering
				XeeXh += (e[i] * e[j]) * (X.row(i).t() * X.row(j));
			}
		}
	}
	
	return XeeXh;
}





/*** R
cli::cli_alert_success("Loaded function {.code XeeX_spatial} and {.code XeeX_serial}")
*/
