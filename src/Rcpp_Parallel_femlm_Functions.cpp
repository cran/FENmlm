#include <cmath>
#include <Rcpp.h>
#include <stdio.h>
#include <omp.h>
#include <Rmath.h>
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]

// This file contains femlm functions that will be parallelized with the omp library

// [[Rcpp::export]]
void set_omp(int n_cpu){
	// function used to set the number of threads used by omp
	// it avoids using it again later on
	omp_set_dynamic(0);
	omp_set_num_threads(n_cpu);
}

double cpppar_abs(double x){
	//simple function to compute the absolute value
	if(x > 0){
		return(x);
	} else {
		return(-x);
	}
}

// [[Rcpp::export]]
NumericVector cpppar_DichotomyNR(int K, int family, double theta, double epsDicho, NumericVector lhs, NumericVector mu, NumericVector borne_inf, NumericVector borne_sup, IntegerVector obsCluster, IntegerVector tableCluster){
	// ARGUMENTS:
	// N: number of observations
	// K: number of clusters classes
	// family:
	// 	- 2: negative binomial
	// 	- 3: logit
	// theta: used only for the NB, but necessary
	// epsDicho: the stopping criterion
	// lhs: the dependent variable
	// mu: the vector of the value of mu
	// borne_inf, init, borne_sup => their length is equal to the number of classes
	// obsCluster: the integer vector that is equal to order(dum[[g]])
	// tableCluster: it gives for each classes the number of element

	// INFO:
	// this algorithm applies dichotomy to find out the cluster coefficients
	// this is for the logit and for the negative binomial
	// after 10 iterations: full dichotomy

	int itermax = 100, iterFullDicho = 10;
	int k; // classical index
	NumericVector res(K);

	IntegerVector start_vector(K), end_vector(K);
	for(k=0 ; k<K ; k++){

		if(k == 0){
			start_vector[k] = 0;
			end_vector[k] = tableCluster[k];
		} else {
			start_vector[k] = end_vector[k-1];
			end_vector[k] = end_vector[k-1] + tableCluster[k];
		}
	}

	#pragma omp parallel
	{
		double x1, x0;
		bool ok;
		int k, i, iter;
		double value, derivee = 0;
		double obs, exp_mu;
		double lower_bound, upper_bound;
		int start, end;

		#pragma omp for private(x1, x0, ok, i, iter, value, derivee, obs, exp_mu, lower_bound, upper_bound, start, end)
		for(k=0 ; k<K ; k++){
			// we loop over each cluster

			// initialisation:
			x1 = 0;
			iter = 0;
			ok = true;

			// the indices of the observations
			start = start_vector[k];
			end = end_vector[k];

			// the bounds
			lower_bound = borne_inf(k);
			upper_bound = borne_sup(k);

			// Update of the value if it goes out of the boundaries
			// because we dont know ex ante if 0 is within the bounds
			if(x1 >= upper_bound || x1 <= lower_bound){
				x1 = (lower_bound + upper_bound)/2;
			}

			while(ok){
				iter++;

				// 1st step: initialisation des bornes
				value = 0;

				for(i=start ; i<end ; i++){
					obs = obsCluster(i); // obscluster must be in cpp index
					if(family == 2){
						// negbin
						value += lhs(obs) - (theta + lhs(obs)) / (1 + theta*exp(-x1 - mu(obs)));
					} else {
						// logit
						value += lhs(obs) - 1 / (1 + exp(-x1 - mu(obs)));
					}

				}

				// update of the bounds.
				if(value > 0){
					lower_bound = x1;
				} else {
					upper_bound = x1;
				}

				// 2nd step: NR iteration or Dichotomy
				x0 = x1;
				if(value == 0){
					ok = false;
				} else if(iter <= iterFullDicho){
					derivee = 0;

					for(i=start ; i<end ; i++){
						obs = obsCluster(i);
						exp_mu = exp(x1 + mu(obs));
						if(family == 2){
							// negbin
							derivee += - theta * (theta + lhs(obs)) / ( (theta/exp_mu + 1) * (theta + exp_mu) );
						} else {
							// logit
							derivee += - 1 / ( (1/exp_mu + 1) * (1 + exp_mu) );
						}
					}

					x1 = x0 - value / derivee;

					// 3rd step: dichotomy (if necessary)
					// Update of the value if it goes out of the boundaries
					if(x1 >= upper_bound || x1 <= lower_bound){
						x1 = (lower_bound + upper_bound)/2;
					}
				} else {
					x1 = (lower_bound + upper_bound)/2;
				}

				// the stopping criteria
				if(iter == itermax){
					ok = false;
				}

				if( ((x0-x1) < 0 && x1-x0 < epsDicho) || ((x0-x1) > 0 && x0-x1 < epsDicho) ){
					ok = false;
				}

			}

			// res[k] = iter;
			// res[k] = value;
			res[k] = x1;
		}
	}

	return(res);
}



// [[Rcpp::export]]
NumericVector cpppar_tapply_vsum(int K, NumericVector x, IntegerVector obsCluster, IntegerVector tableCluster){
	// K: nber of classes
	// x: a vector
	// dum: the N vector of clusters

	int k;

	NumericVector res(K);

	IntegerVector start_vector(K), end_vector(K);
	for(k=0 ; k<K ; k++){

		if(k == 0){
			start_vector[k] = 0;
			end_vector[k] = tableCluster[k];
		} else {
			start_vector[k] = end_vector[k-1];
			end_vector[k] = end_vector[k-1] + tableCluster[k];
		}
	}

	#pragma omp parallel
	{
		int k, i;
		#pragma omp for private(i)
		for(k=0 ; k<K ; k++){
			for(i=start_vector[k] ; i<end_vector[k] ; i++){
				res[k] += x[obsCluster(i)];
			}
		}
	}

	return(res);
}

// [[Rcpp::export]]
NumericMatrix cpppar_PartialDerivative(int Q, int N, int V, double epsDeriv, NumericVector ll_d2,	NumericMatrix F, NumericVector init, IntegerMatrix obsCluster, IntegerVector tableCluster, IntegerVector nbCluster){
	// takes in:
	// init: the initialisation of the sum of derivatives vector
	// ll_d2: the second derivative
	// F: the jacobian matrix
	// V: the number of variables
	// tableCluster: unlist(tableCluster_all)

	int iterMax = 1000;

	int q, k;
	int K;
	int index;
	int sum_cases=0;
	IntegerVector start_cluster(Q), end_cluster(Q);
	// start/end_cluster: the beginning and end of the whole clusters

	for(q=0 ; q<Q ; q++){
		sum_cases += nbCluster(q);
		if(q == 0){
			start_cluster(q) = 0;
			end_cluster(q) = nbCluster(q);
		} else {
			start_cluster(q) = start_cluster(q-1) + nbCluster(q-1);
			end_cluster(q) = end_cluster(q-1) + nbCluster(q);
		}
	}

	// start/end_cluster_coef: the start and end, observation wise, of cluster coefficients in the obsCluster matrix
	IntegerVector start_cluster_coef(sum_cases), end_cluster_coef(sum_cases);

	for(q=0 ; q<Q ; q++){
		K = end_cluster[q] - start_cluster[q];
		for(k=0 ; k<K ; k++){

			index = start_cluster[q] + k;
			if(k == 0){
				start_cluster_coef[index] = 0;
				end_cluster_coef[index] = tableCluster[index];
			} else {
				start_cluster_coef[index] = end_cluster_coef[index - 1];
				end_cluster_coef[index] = end_cluster_coef[index - 1] + tableCluster[index];
			}
		}
	}

	// Rprintf("start_cluster_coef: ok\n");

	NumericMatrix clusterDeriv(sum_cases, V); // the 'temporary' derivatives
	NumericVector sum_lld2(sum_cases);
	// the result to return
	NumericMatrix S(N, V);

	// Creation of the sum_lld2
	for(q=0 ; q<Q ; q++){
		K = end_cluster[q] - start_cluster[q];
		#pragma omp parallel
		{
			double value_thread;
			int obs, k, index, i;

			#pragma omp for private(index, i, value_thread, obs)
			for(k=0 ; k<K ; k++){
				index = start_cluster[q] + k;
				value_thread = 0;
				for(i=start_cluster_coef[index] ; i<end_cluster_coef[index] ; i++){
					obs = obsCluster(i, q);
					value_thread += ll_d2[obs];
				}
				sum_lld2[index] = value_thread;
			}
		}
	}

	#pragma omp parallel
	{
		// initialisation of S
		#pragma omp for collapse(2)
		for(int v=0 ; v<V ; v++){
			for(int i=0 ; i<N ; i++){
				S(i,v) = init(i,v);
			}
		}
	}


	#pragma omp parallel
	{

		// variables
		int q, i, k, obs, index, K;
		int iter;
		double value, max_diff;
		bool ok;

		#pragma omp for private(q, i, k, obs, index, K, iter, value, max_diff, ok)
		for(int v=0 ; v<V ; v++){
			// we loop over each variable

			ok = true;
			iter = 0;
			while( ok & (iter<iterMax) ){
				iter++;
				ok = false;

				// Rprintf("\n iter=%i -- ", iter);

				for(q=0 ; q<Q ; q++){

					// Rprintf(" q=%i ", q);

					// sum for each cluster
					K = end_cluster[q] - start_cluster[q];

					// computing deriv
					for(k=0 ; k<K ; k++){
						index = start_cluster[q] + k;
						value = 0;
						for(i=start_cluster_coef[index] ; i<end_cluster_coef[index] ; i++){
							obs = obsCluster(i, q);
							value += (F(obs,v) + S(obs,v))*ll_d2(obs);
						}
						clusterDeriv(index, v) = -value/sum_lld2[index];
					}

					// update of S
					for(k=0 ; k<K ; k++){
						index = start_cluster[q] + k;
						for(i=start_cluster_coef[index] ; i<end_cluster_coef[index] ; i++){
							obs = obsCluster(i, q);
							S(obs, v) += clusterDeriv(index, v);
						}
					}

					// Rprintf("q=%i - clusterDeriv(0, %i)=%f\n", q, v, clusterDeriv(0, v));
					// Rprintf("q=%i - clusterDeriv(4393, %i)=%f\n", q, v, clusterDeriv(4393, v));
					// Rprintf("q=%i - S(0,%i)=%f\n", q, v, S(0, v));


					// control
					if(ok == false){
						for(k=0 ; k<K ; k++){
							index = start_cluster[q] + k;
							if(cpppar_abs(clusterDeriv(index, v)) > epsDeriv){
								ok = true;
								break;
							}
						}
					}
				}
			}

			// Rprintf("K=%i, nb iter=%i\n", v, iter);
			if(iter == iterMax){
				// on calcule la max diff pour info
				max_diff = 0;
				for(index=0 ; index<sum_cases ; index++){
					if( cpppar_abs(clusterDeriv(index, v)) > max_diff ){
						max_diff = cpppar_abs(clusterDeriv(index, v));
					}
				}

				Rprintf("[Getting cluster deriv.] Max iterations reached (%i) (max diff.: %0.6f)\n", iterMax, max_diff);
			}

		}
	}

	return(S);
}

// [[Rcpp::export]]
NumericVector cpppar_exp(NumericVector x){
	// parallel exponentiation using omp

	int n = x.length();
	NumericVector res(n);

	#pragma omp parallel
	{
		#pragma omp for
		for(int i = 0 ; i < n ; i++) {
			res[i] = exp(x[i]);
		}
	}

	return(res);
}

// [[Rcpp::export]]
NumericVector cpppar_log(NumericVector x){
	// parallel exponentiation using omp

	int n = x.length();
	NumericVector res(n);

	#pragma omp parallel
	{
		#pragma omp for
		for(int i = 0 ; i < n ; i++) {
			res[i] = log(x[i]);
		}
	}

	return(res);
}

// [[Rcpp::export]]
NumericVector cpppar_log_a_exp(double a, NumericVector mu, NumericVector exp_mu){
	// faster this way

	int n = mu.length();
	NumericVector res(n);

	#pragma omp parallel
	{
		#pragma omp for
		for(int i=0 ; i<n ; i++) {
			if(mu[i] < 200){
				res[i] = log(a + exp_mu[i]);
			} else {
				res[i] = mu[i];
			}
		}
	}

	return(res);
}

// [[Rcpp::export]]
NumericVector cpppar_lgamma(NumericVector x){
	// parallel lgamma using omp

	int n = x.length();
	NumericVector res(n);

	#pragma omp parallel
	{
		#pragma omp for
		for(int i = 0 ; i < n ; i++) {
			res[i] = lgamma(x[i]);
		}
	}

	return(res);
}

// [[Rcpp::export]]
NumericVector cpppar_digamma(NumericVector x){
	// parallel digamma using omp

	int n = x.length();
	NumericVector res(n);

	#pragma omp parallel
	{
		#pragma omp for
		for(int i = 0 ; i < n ; i++) {
			res[i] = R::digamma(x[i]);
		}
	}

	return(res);
}

// [[Rcpp::export]]
NumericVector cpppar_trigamma(NumericVector x){
	// parallel trigamma using omp

	int n = x.length();
	NumericVector res(n);

	#pragma omp parallel
	{
		#pragma omp for
		for(int i = 0 ; i < n ; i++) {
			res[i] = R::trigamma(x[i]);
		}
	}

	return(res);
}







































































