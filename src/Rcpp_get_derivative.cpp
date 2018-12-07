#include <Rcpp.h>
#include <stdio.h>

using namespace Rcpp;

// This file contains the algorithm to get the partial derivative wrt dummies

// Different types:
// int, IntegerVector
// NumericVector, NumericMatrix
// bool
// Rcpp::stop("problem here");
// R_CheckUserInterrupt(); // to allow the user to stop the algo if the process is too long

// Some mathematical expressions that might be problematic
// exp log
// Attention: pour les fonctions log ou exp: il faut absolument du numeric vector (et pas du integer ou ca cree du overloading)

double cpp_abs(double x){
	//simple function to compute the absolute value
	if(x >= 0){
		return(x);
	} else {
		return(-x);
	}
}

// [[Rcpp::export]]
NumericVector cpp_lgamma(NumericVector x){
	// simple function to compute lgamma of a vector

	int n = x.length();
	NumericVector res(n);

	for(int i=0 ; i<n ; i++){
		res[i] = lgamma(x[i]);
	}

	return(res);
}

// [[Rcpp::export]]
NumericVector cpp_log_a_exp(double a, NumericVector mu, NumericVector exp_mu){
	// faster this way

	int n = mu.length();
	NumericVector res(n);

	for(int i=0 ; i<n ; i++){
		if(mu[i] < 200){
			res[i] = log(a + exp_mu[i]);
		} else {
			res[i] = mu[i];
		}
	}

	return(res);
}


// [[Rcpp::export]]
NumericMatrix RcppPartialDerivative(int iterMax, int Q, int N, int K, double epsDeriv, NumericVector ll_d2,	NumericMatrix F, NumericVector init, IntegerMatrix dumMat, IntegerVector nbCluster){
	// takes in:
	// dumMat: the matrix of dummies (n X c) each obs => cluster
	// init: the initialisation of the sum of derivatives vector
	// ll_d2: the second derivative
	// F: the jacobian matrix

	int iter;

	int i, q, k, c;
	int index;
	int sum_cases=0;
	bool ok;
	double new_value;
	IntegerVector start(Q), end(Q);

	for(q=0 ; q<Q ; q++){
		sum_cases += nbCluster(q);
		if(q == 0){
			start(q) = 0;
			end(q) = nbCluster(q);
		} else {
			start(q) = start(q-1) + nbCluster(q-1);
			end(q) = end(q-1) + nbCluster(q);
		}
	}

	NumericMatrix clusterDeriv(sum_cases, K); // the output
	NumericVector sum_lld2(sum_cases);

	// Creation of the sum_lld2
	for(i=0 ; i<N ; i++){
		for(q=0 ; q<Q ; q++){
			// index = start[q] + dumMat(i, q) - 1;
			index = start[q] + dumMat(i, q);
			sum_lld2[index] += ll_d2(i);
		}
	}

	// the result to return
	NumericMatrix S(N, K);
	for(i=0 ; i<N ; i++){
		for(k=0 ; k<K ; k++){
			S(i,k) = init(i,k);
		}
	}

	for(k=0 ; k<K ; k++){

		ok = true;
		iter = 0;
		while( ok & (iter<iterMax) ){
			iter++;
			ok = false;

			for(q=0 ; q<Q ; q++){
				R_CheckUserInterrupt();

				// init of the vector
				for(c=start[q] ; c<end[q] ; c++){
					clusterDeriv(c, k) = 0;
				}

				for(i=0 ; i<N ; i++){
					// index = start[q] + dumMat(i, q) - 1;
					index = start[q] + dumMat(i, q);
					clusterDeriv(index, k) += (F(i,k) + S(i,k))*ll_d2(i);
				}

				// on finit de calculer clusterDeriv + controle
				for(c=start[q] ; c<end[q] ; c++){
					new_value = -clusterDeriv(c, k)/sum_lld2[c];
					clusterDeriv(c, k) = new_value;
					if(cpp_abs(new_value) > epsDeriv){
						ok = true;
					}
				}

				// on ajoute la derivee a S:
				for(i=0 ; i<N ; i++){
					// index = start[q] + dumMat(i, q) - 1;
					index = start[q] + dumMat(i, q);
					S(i,k) += clusterDeriv(index, k);
				}

			}
		}

		// Rprintf("K=%i, nb iter=%i\n", k, iter);
		if(iter == iterMax){
			// on clalcule la max diff pour info
			double max_diff = 0;
			for(index=0 ; index<sum_cases ; index++){
				if( cpp_abs(clusterDeriv(index, k)) > max_diff ){
					max_diff = cpp_abs(clusterDeriv(index, k));
				}
			}

			Rprintf("[(cpp loop) Getting cluster deriv.] Max iterations reached (%i) (max diff.: %0.6f)\n", iterMax, max_diff);
		}

	}

	return(S);
}

// [[Rcpp::export]]
NumericMatrix RcppPartialDerivative_gaussian(int iterMax, int Q, int N, int K, double epsDeriv, NumericMatrix F, NumericVector init, IntegerMatrix dumMat, IntegerVector nbCluster){
	// takes in:
	// dumMat: the matrix of dummies (n X c) each obs => cluster
	// init: the initialisation of the sum of derivatives vector
	// ll_d2: the second derivative
	// F: the jacobian matrix

	int iter;

	int i, q, k, c;
	int index;
	int sum_cases=0;
	bool ok;
	double new_value;
	IntegerVector start(Q), end(Q);

	for(q=0 ; q<Q ; q++){
		sum_cases += nbCluster(q);
		if(q == 0){
			start(q) = 0;
			end(q) = nbCluster(q);
		} else {
			start(q) = start(q-1) + nbCluster(q-1);
			end(q) = end(q-1) + nbCluster(q);
		}
	}

	NumericMatrix clusterDeriv(sum_cases, K); // the output

	NumericVector n_group(sum_cases);

	// Creation of the sum_lld2
	for(i=0 ; i<N ; i++){
		for(q=0 ; q<Q ; q++){
			// index = start[q] + dumMat(i, q) - 1;
			index = start[q] + dumMat(i, q);
			n_group[index] ++;
		}
	}

	// the result to return
	NumericMatrix S(N, K);
	for(i=0 ; i<N ; i++){
		for(k=0 ; k<K ; k++){
			S(i,k) = init(i,k);
		}
	}

	for(k=0 ; k<K ; k++){

		ok = true;
		iter = 0;
		while( ok & (iter<iterMax) ){
			iter++;
			ok = false;

			for(q=0 ; q<Q ; q++){
				R_CheckUserInterrupt();

				// init of the vector
				for(c=start[q] ; c<end[q] ; c++){
					clusterDeriv(c, k) = 0;
				}

				for(i=0 ; i<N ; i++){
					// index = start[q] + dumMat(i, q) - 1;
					index = start[q] + dumMat(i, q);
					clusterDeriv(index, k) += F(i,k) + S(i,k);
				}

				// on finit de calculer clusterDeriv + controle
				for(c=start[q] ; c<end[q] ; c++){
					new_value = -clusterDeriv(c, k)/n_group[c];
					clusterDeriv(c, k) = new_value;
					if(cpp_abs(new_value) > epsDeriv/Q){
						ok = true;
					}
				}

				// on ajoute la derivee a S:
				for(i=0 ; i<N ; i++){
					// index = start[q] + dumMat(i, q) - 1;
					index = start[q] + dumMat(i, q);
					S(i,k) += clusterDeriv(index, k);
				}

			}
		}

		// Rprintf("K=%i, nb iter=%i\n", k, iter);
		if(iter == iterMax){
			// on clalcule la max diff pour info
			double max_diff = 0;
			for(index=0 ; index<sum_cases ; index++){
				if( cpp_abs(clusterDeriv(index, k)) > max_diff ){
					max_diff = cpp_abs(clusterDeriv(index, k));
				}
			}

			Rprintf("[Getting cluster deriv.] Max iterations reached (%i) (max diff.: %0.6f)\n", iterMax, max_diff);
		}

	}

	return(S);
}

// [[Rcpp::export]]
NumericVector RcppPartialDerivative_other(int iterMax, int Q, int N, double epsDeriv, NumericVector ll_d2,	NumericVector dx_dother, NumericVector init, IntegerMatrix dumMat, IntegerVector nbCluster){
	// takes in:
	// dumMat: the matrix of dummies (n X c) each obs => cluster // must be in cpp index!!!
	// init: the initialisation of the sum of derivatives vector
	// ll_d2: the second derivative
	// dx_dother: the vector of dx_dother

	int iter;

	int i, q, c;
	int index;
	int sum_cases=0;
	bool ok;
	double new_value;
	IntegerVector start(Q), end(Q);

	for(q=0 ; q<Q ; q++){
		// the total number of clusters (eg if man/woman and 10 countries: total of 12 cases)
		sum_cases += nbCluster(q);
		if(q == 0){
			start(q) = 0;
			end(q) = nbCluster(q);
		} else {
			start(q) = start(q-1) + nbCluster(q-1);
			end(q) = end(q-1) + nbCluster(q);
		}
	}

	NumericVector clusterDeriv(sum_cases); // the derivatives per cluster, stacked in one vector
	NumericVector sum_lld2(sum_cases);

	// Creation of the sum_lld2
	for(i=0 ; i<N ; i++){
		for(q=0 ; q<Q ; q++){
			index = start[q] + dumMat(i, q);
			sum_lld2[index] += ll_d2(i);
		}
	}

	// the result to return
	NumericVector S(N);
	for(i=0 ; i<N ; i++){
		S[i] = init(i);
	}

	ok = true;
	iter = 0;
	while( ok & (iter<iterMax) ){
		iter++;
		ok = false;

		for(q=0 ; q<Q ; q++){
			R_CheckUserInterrupt();

			// init of the vector
			for(c=start[q] ; c<end[q] ; c++){
				clusterDeriv(c) = 0;
			}

			for(i=0 ; i<N ; i++){
				index = start[q] + dumMat(i, q);
				clusterDeriv(index) += dx_dother(i) + S(i)*ll_d2(i);
			}

			// on finit de calculer clusterDeriv + controle
			for(c=start[q] ; c<end[q] ; c++){
				new_value = -clusterDeriv(c)/sum_lld2[c];
				clusterDeriv(c) = new_value;
				if(cpp_abs(new_value) > epsDeriv){
					ok = true;
				}
			}

			// on ajoute la derivee a S:
			for(i=0 ; i<N ; i++){
				index = start[q] + dumMat(i, q);
				S[i] += clusterDeriv(index);
			}

		}
	}

	// Rprintf("other, nb iter=%i\n", iter);
	if(iter == iterMax){
		Rprintf("[Getting cluster deriv. other] Max iterations reached (%i)\n", iterMax);
	}

	return(S);
}

// [[Rcpp::export]]
List RcppGetFE(int Q, int N, NumericVector S, IntegerMatrix dumMat, IntegerVector nbCluster, IntegerVector obsCluster){
	// This function returns a list of the cluster coefficients for each cluster
	// dumMat: the matrix of cluster ID for each observation, with cpp index style
	// Q, N: nber of clusters / obs
	// nbCluster: vector of the number of cases per cluster
	// obsCluster: the integer vector that is equal to order(dum[[g]])
	// RETURN:
	// a list for each cluster of the cluster coefficient value
	// + the last element is the number of clusters that have been set as references (nb_ref)


	int iter=0, iterMax=10000;
	int iter_loop=0, iterMax_loop=10000;


	// Creation of the indices to put all the cluster values into a single vector
	int sum_cases=0;
	IntegerVector start(Q), end(Q), nb_ref(Q); // nb_ref takes the nb of elements set as ref
	int q;

	for(q=0 ; q<Q ; q++){
		// the total number of clusters (eg if c1: man/woman and c2: 10 countries: total of 12 cases)
		sum_cases += nbCluster(q);
		if(q == 0){
			start(q) = 0;
			end(q) = nbCluster(q);
		} else {
			start(q) = start(q-1) + nbCluster(q-1);
			end(q) = end(q-1) + nbCluster(q);
		}
	}

	NumericVector cluster_values(sum_cases);

	// Now we create the vector of observations for each cluster
	// we need a strating and an end vector as well
	IntegerVector start_cluster(sum_cases), end_cluster(sum_cases);

	int i, k, index;

	for(q=0 ; q<Q ; q++){

		// table cluster: nber of elements for each cluster class
		IntegerVector tableCluster(nbCluster(q));
		for(i=0 ; i<N ; i++){
			k = dumMat(i, q);
			tableCluster(k) += 1; // the number of obs per case
		}

		// now creation of the start/end vectors
		for(k=0 ; k<nbCluster(q) ; k++){
			index = start(q) + k;
			if(k == 0){
				start_cluster[index] = 0;
				end_cluster[index] = tableCluster[k];
			} else {
				start_cluster[index] = end_cluster[index-1];
				end_cluster[index] = end_cluster[index-1] + tableCluster[k];
			}
		}
	}

	// matrix of the clusters that have been computed
	IntegerMatrix mat_done(N, Q);

	// vector of elements to loop over
	IntegerVector id2do(N);
	int nb2do = N;
	for(i=0 ; i<nb2do ; i++){
		id2do(i) = i;
	}

	// Other indices and variables
	int qui_max, obs;
	int rs, rs_max;
	int j, ind, id_cluster;
	double other_value;
	bool first, ok;

	//
	// THE MAIN LOOP
	//

	while(iter < iterMax){
		iter++;

		//
		// Finding the row where to put the 0s
		//


		if(iter == 1){
			// 1st iter, we select the first element
			qui_max = 0;
		} else {
			// we find the row that has the maximum of items done

			qui_max = 0;
			rs_max = 0;
			for(i=0 ; i<nb2do ; i++){
				obs = id2do[i];

				rs = 0;
				for(q=0 ; q<Q ; q++){
					rs += mat_done(obs, q);
				}

				if(rs == Q-2){
					// if rs == Q-2 => its the maximum possible, no need to go further
					qui_max = obs;
					break;
				} else if(rs<Q && rs>rs_max){
					// this is to handle complicated cases with more than 3+ clusters
					qui_max = obs;
					rs_max = rs;
				} else if(qui_max == 0 && rs == 0){
					// used to initialize qui_max
					qui_max = obs;
				}
			}
		}

		//
		// Putting the 0s, ie setting the references
		//

		// the first element is spared
		first = true;
		for(q=0 ; q<Q ; q++){
			if(mat_done(qui_max, q) == 0){
				if(first){
					// we spare the first element
					first = false;
				} else {
					// we set the cluster to 0
					// 1) we find the cluster
					id_cluster = dumMat(qui_max, q);
					// Rprintf("Cluster: %i\n", id_cluster + 1);
					// 2) we get the index of the cluster vector
					index = start[q] + id_cluster;
					// 3) we set the cluster value to 0
					cluster_values(index) = 0;
					// 4) we update the mat_done matrix for the elements of this cluster
					for(i=start_cluster[index] ; i<end_cluster[index] ; i++){
						obs = obsCluster(i, q);
						mat_done(obs, q) = 1;
					}
					// 5) => we save the information on which cluster was set as a reference
					nb_ref(q)++;
				}
			}
		}

		//
		// LOOP OF ALL OTHER UPDATES (CRITICAL)
		//

		iter_loop = 0;
		while(iter_loop < iterMax_loop){
			iter_loop++;

			R_CheckUserInterrupt();

			//
			// Selection of indexes (new way) to be updated
			//

			IntegerVector qui_indexes(sum_cases), qui_obs(sum_cases); // init at the max
			int nb_index = 0;

			for(i=0 ; i<nb2do ; i++){
				// we compute the rowsum of the obs that still needs to be done
				obs = id2do[i];

				rs = 0;
				for(q=0 ; q<Q ; q++){
					rs += mat_done(obs, q);
				}

				if(rs == Q-1){
					// means: needs to be updated
					for(q=0 ; mat_done(obs, q)!=0 ; q++){
						// selection of the q that is equal to 0
					}

					index = start(q) + dumMat(obs, q);

					ok = true;
					for(j=0 ; j<nb_index ; j++){
						if(index == qui_indexes[j]){
							ok = false;
							break;
						}
					}

					if(ok){
						// the index was not already there
						qui_indexes[nb_index] = index; // the 'cluster index'
						qui_obs[nb_index] = obs; // the observation to be updated
						nb_index++;
					}
				}
			}

			if(nb_index == 0) break;

			// Rprintf("nb = %i\n", nb_index);

			// => we obtain a unique list of index to be updated

			//
			// We update each index
			//

			for(ind=0 ; ind<nb_index ; ind++){

				int index_select = qui_indexes[ind];
				int q_select=0;

				other_value = 0;
				// Computing the sum of the other cluster values
				// and finding the cluster to be updated (q_select)
				for(q=0 ; q<Q ; q++){
					// we can loop over all q because cluster_values is initialized to 0
					obs = qui_obs[ind];
					index = start(q) + dumMat(obs, q);
					other_value += cluster_values(index);

					if(index == index_select){
						q_select = q;
					}

				}

				// the index to update
				cluster_values(index_select) = S(obs) - other_value;

				// Update of the mat_done and the id2do
				for(i=start_cluster[index_select] ; i<end_cluster[index_select] ; i++){
					obs = obsCluster(i, q_select);
					mat_done(obs, q_select) = 1;
				}
			}
		}

		// Check that everything is all right
		if(iter_loop == iterMax_loop){
			Rprintf("Problem getting FE, maximum iterations reached (2nd order loop).");
		}

		// now the control
		IntegerVector id2do_new(nb2do);
		int nb2do_new = 0;

		for(i=0 ; i<nb2do ; i++){
			obs = id2do[i];

			rs = 0;
			for(q=0 ; q<Q ; q++){
				rs += mat_done(obs, q);
			}

			if(rs < Q){
				id2do_new[nb2do_new] = obs;
				nb2do_new++;
			}
		}

		if(nb2do_new == 0) break;

		// update id2do
		IntegerVector id2do(nb2do_new);
		int nb2do = nb2do_new;

		for(i=0 ; i<nb2do ; i++){
			id2do[i] = id2do_new[i];
		}

	}

	if(iter == iterMax){
		Rprintf("Problem getting FE, maximum iterations reached (1st order loop).");
	}

	// final formatting and save
	List res(Q + 1);
	for(q=0 ; q<Q ; q++){
		NumericVector quoi(nbCluster(q));
		for(k=0 ; k<nbCluster(q) ; k++){
			index = start(q) + k;
			quoi(k) = cluster_values(index);
		}
		res(q) = quoi;
	}

	res(Q) = nb_ref;

	return(res);
}

// Function to get the evaluation of the dummy of the NB
// [[Rcpp::export]]
double cpp_NB_dum_fx(double theta, NumericVector lhs, NumericVector mu, double x1, IntegerVector obsCluster, int start, int end){

	double res=0;
	int i, obs;
	// double exp_mu;

	for(i=start ; i<end ; i++){

		obs = obsCluster(i); // obscluster must be in cpp index

		// exp_mu = exp(x1 + mu(obs));
		// res += lhs(obs) - (theta + lhs(obs)) * exp_mu / (exp_mu + theta);

		// I rewrite it to avoid problems of identification
		res += lhs(obs) - (theta + lhs(obs)) / (1 + theta*exp(-x1 - mu(obs)));
	}

	return(res);
}

// Function to get the evaluation of the dummy for the Logit
double cpp_LOGIT_dum_fx(NumericVector lhs, NumericVector mu, double x1, IntegerVector obsCluster, int start, int end){

	double res=0;
	int i, obs;
	// double exp_mu;

	for(i=start ; i<end ; i++){

		obs = obsCluster(i); // obscluster must be in cpp index

		// exp_mu = exp(x1 + mu(obs));
		// res += lhs(obs) - exp_mu / (1 + exp_mu);

		// I rewrite it to avoid problems of identification
		res += lhs(obs) - 1 / (1 + exp(-x1 - mu(obs)));
	}

	return(res);
}

// Function to get the evaluation of the dummy
double cpp_dum_fx(int family, double theta, NumericVector lhs, NumericVector mu, double x1, IntegerVector obsCluster, int start, int end){

	double res;

	if(family == 2){
		// => the negative binomial
		res = cpp_NB_dum_fx(theta, lhs, mu, x1, obsCluster, start, end);
	} else if(family == 3){
		// => Logit
		res = cpp_LOGIT_dum_fx(lhs, mu, x1, obsCluster, start, end);
	} else {
		stop("Unrecognized family");
	}

	return(res);
}

// [[Rcpp::export]]
double cpp_NB_dum_dfx(double theta, NumericVector lhs, NumericVector mu, double x1, IntegerVector obsCluster, int start, int end){

	double res=0;
	int i, obs;
	double exp_mu;
	// double new_mu;

	for(i=start ; i<end ; i++){

		obs = obsCluster(i); // obscluster must be in cpp index

		// exp_mu = exp(x1 + mu(obs));
		// res += - theta * (theta + lhs(obs)) * exp_mu / (theta + exp_mu) / (theta + exp_mu);

		// I rewrite it to avoid problems of identification
		exp_mu = exp(x1 + mu(obs));
		res += - theta * (theta + lhs(obs)) / ( (theta/exp_mu + 1) * (theta + exp_mu) );
	}

	return(res);
}

// Function to get the evaluation of the dummy for the Logit
double cpp_LOGIT_dum_dfx(NumericVector lhs, NumericVector mu, double x1, IntegerVector obsCluster, int start, int end){

	double res=0;
	int i, obs;
	double exp_mu;

	for(i=start ; i<end ; i++){

		obs = obsCluster(i); // obscluster must be in cpp index

		// exp_mu = exp(x1 + mu(obs));
		// res += - exp_mu / (1 + exp_mu) / (1 + exp_mu);

		// I rewrite it to avoid problems of identification
		exp_mu = exp(x1 + mu(obs));
		res += - 1 / ( (1/exp_mu + 1) * (1 + exp_mu) );
	}

	return(res);
}

// Function to get the evaluation of the dummy
double cpp_dum_dfx(int family, double theta, NumericVector lhs, NumericVector mu, double x1, IntegerVector obsCluster, int start, int end){

	double res;

	if(family == 2){
		// => the negative binomial
		res = cpp_NB_dum_dfx(theta, lhs, mu, x1, obsCluster, start, end);
	} else if(family == 3){
		// => Logit
		res = cpp_LOGIT_dum_dfx(lhs, mu, x1, obsCluster, start, end);
	} else {
		stop("Unrecognized family");
	}

	return(res);
}

// [[Rcpp::export]]
NumericMatrix RcppCreate_start_end_indexes(IntegerVector nbCluster, IntegerVector tableCluster_vect){
	// this function creates the start and end indexes for the clusters to be used in the dichotomy function
	// this is very specific

	int sum_cases = tableCluster_vect.size();

	NumericMatrix start_end(sum_cases, 2);

	for(int k=0 ; k<sum_cases ; k++){

		if(k == 0){
			start_end(k, 0) = 0; // start
			start_end(k, 1) = tableCluster_vect[k]; // end
		} else {
			start_end(k, 0) = start_end(k - 1, 1); // start, end
			start_end(k, 1) = start_end(k - 1, 1) + tableCluster_vect[k]; // end, end
		}
	}

	return(start_end);
}

// [[Rcpp::export]]
NumericVector cpp_DichotomyNR(int N, int K, int family, double theta, double epsDicho, NumericVector lhs, NumericVector mu, NumericVector borne_inf, NumericVector borne_sup, IntegerVector obsCluster, IntegerVector tableCluster){
	// ARGUMENTS:
	// N: number of observations
	// K: number of clusters classes
	// family:
	// 	- 1: negative binomial
	// 	- 2: logit
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

	int itermax = 100, iter, iterFullDicho = 10;
	double value, x0, derivee = 0; // evaluation of the dummy
	int k; // classical index
	int start, end; // the indices used to compute sums wrt clusters
	NumericVector res(K);
	double lower_bound, upper_bound;

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


	for(int k=0 ; k<K ; k++){
		// we loop over each cluster

		double x1 = 0; // we initialise the cluster coefficient at 0 (it should converge to 0)
		bool ok = true;
		iter = 0;

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
			R_CheckUserInterrupt();
			iter++;

			// 1st step: initialisation des bornes
			value = cpp_dum_fx(family, theta, lhs, mu, x1, obsCluster, start, end);

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
				derivee = cpp_dum_dfx(family, theta, lhs, mu, x1, obsCluster, start, end);
				x1 = x0 - value / derivee;
				// Rprintf("x1: %5.2f\n", x1);

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
				Rprintf("[Getting cluster coefficients nber %i] max iterations reached (%i).\n", k, itermax);
				Rprintf("Value Sum Deriv (NR) = %f. Difference = %f.\n", value, cpp_abs(x0-x1));
			}

			if(cpp_abs(x0-x1) < epsDicho){
				ok = false;
			}

		}
		// Rprintf("iter=%i.\n", iter);
		// res[k] = iter;
		// res[k] = value;
		res[k] = x1;
	}

	return(res);
}

NumericVector DichotomyNR_internal(int family, int K, int start, double theta, double epsDicho, NumericVector lhs, NumericVector mu, NumericVector borne_inf, NumericVector borne_sup, IntegerVector obsCluster_vect, IntegerVector start_vector, IntegerVector end_vector){
	// ARGUMENTS:
	// N: number of observations
	// K: number of clusters classes
	// family:
	// 	- 1: negative binomial
	// 	- 2: logit
	// theta: used only for the NB, but necessary
	// epsDicho: the stopping criterion
	// lhs: the dependent variable
	// mu: the vector of the value of mu
	// borne_inf, init, borne_sup => their length is equal to the number of classes
	// obsCluster: the integer vector that is equal to order(dum[[g]])

	// INFO:
	// this algorithm applies dichotomy to find out the cluster coefficients
	// this is for the logit and for the negative binomial

	int itermax = 1000, iter;
	double value, x0, derivee;; // evaluation of the dummy
	int k; // classical index
	int start_cluster, end_cluster; // the indices used to compute sums wrt clusters
	NumericVector res(K);
	double lower_bound, upper_bound;

	for(k=0 ; k<K ; k++){
		// we loop over each cluster

		double x1 = 0; // we initialise the cluster coefficient at 0 (it should converge to 0)
		bool ok = true;
		iter = 0;

		// the indices of the observations
		start_cluster = start_vector(start + k);
		end_cluster = end_vector(start + k);

		// the bounds
		lower_bound = borne_inf(k);
		upper_bound = borne_sup(k);

		// Update of the value if it goes out of the boundaries
		if(x1 >= upper_bound || x1 <= lower_bound){
			x1 = (lower_bound + upper_bound)/2;
		}

		while(ok){
			R_CheckUserInterrupt();
			iter++;

			// 1st step: initialisation des bornes
			value = cpp_dum_fx(family, theta, lhs, mu, x1, obsCluster_vect, start_cluster, end_cluster);

			// update of the bounds.
			if(value > 0){
				lower_bound = x1;
			} else {
				upper_bound = x1;
			}

			// 2nd step: NR iteration
			x0 = x1;
			derivee = cpp_dum_dfx(family, theta, lhs, mu, x1, obsCluster_vect, start_cluster, end_cluster);
			x1 = x0 - value / derivee;
			// Rprintf("x1: %5.2f\n", x1);

			// 3rd step: dichotomy (if necessary)
			// Update of the value if it goes out of the boundaries
			if(x1 >= upper_bound || x1 <= lower_bound){
				x1 = (lower_bound + upper_bound)/2;
			}

			// the stopping criteria
			if(iter == itermax){
				ok = false;
				Rprintf("[Getting cluster coefficients] max iterations reached (%i).\n", itermax);
				Rprintf("Value Sum Deriv = %f.\n", value);
			}
			if(cpp_abs(x0 - x1) < epsDicho){
				ok = false;
			}
		}
		// Rprintf("iter=%i.\n", iter);
		res[k] = x1;
		// Rprintf("%3.3f ", x1);
	}
	// Rprintf("\n");

	return(res);
}

// Function to get the conditional min, the max and the mean of a vector
// [[Rcpp::export]]
NumericMatrix cpp_conditional_minMaxMean(int K, int N, NumericVector mu, IntegerVector dum, IntegerVector tableCluster){
	// K: nber of classes
	// N: nber of observations
	// mu: a vector
	// dum: the N vector of clusters
	// tableCluster: the number of obs per cluster

	// returns a Q x 3 matrix

	NumericMatrix res(K,3);
	IntegerVector ok(K);
	int i, k;

	for(i=0 ; i<N ; i++){
		k = dum(i) - 1; // we take 1 off => different indexation in C

		if(ok[k] == 0){
			// means its the first time
			res(k, 0) = res(k, 1) = res(k, 2) = mu(i);
			ok[k] = 1;
		} else {
			// the min
			if(mu(i) < res(k, 0)){
				res(k, 0) = mu(i);
			}

			// the max
			if(mu(i) > res(k, 1)){
				res(k, 1) = mu(i);
			}

			// the mean
			res(k, 2) += mu(i);
		}
	}

	// the mean
	for(k=0 ; k<K ; k++) res(k, 2) = res(k, 2)/tableCluster(k);

	return(res);
}

// [[Rcpp::export]]
NumericMatrix cpp_conditional_minMax(int K, int N, NumericVector mu, IntegerVector dum, IntegerVector tableCluster){
	// Q: nber of classes
	// N: nber of observations
	// mu: a vector
	// dum: the N vector of clusters
	// tableCluster: the number of obs per cluster

	// returns a Q x 3 matrix

	NumericMatrix res(K,2);
	IntegerVector ok(K); // init a 0 par defaut
	int i, k;

	for(i=0 ; i<N ; i++){
		k = dum(i) - 1; // we take 1 off => different indexation in C

		if(ok[k] == 0){
			// means its the first time
			res(k, 0) = res(k, 1) = mu(i);
			ok[k] = 1;
		} else {
			// the min
			if(mu(i) < res(k, 0)){
				res(k, 0) = mu(i);
			}

			// the max
			if(mu(i) > res(k, 1)){
				res(k, 1) = mu(i);
			}
		}
	}

	return(res);
}

// Function to get the conditional sum of a matrix
// [[Rcpp::export]]
NumericMatrix cpp_tapply_sum(int Q, NumericMatrix x, IntegerVector dum){
	// Q: nber of classes
	// N: nber of observations
	// x: a matrix
	// dum: the N vector of clusters

	int N = x.nrow();
	int K = x.ncol();

	NumericMatrix res(Q, K);
	int i, q, k;

	for(i=0 ; i<N ; i++){
		q = dum(i) - 1; // we take 1 off => different indexation in C

		for(k=0 ; k<K ; k++){
			res(q, k) += x(i, k);
		}
	}

	return(res);
}

// Function to get the conditional sum of a vector
// [[Rcpp::export]]
NumericVector cpp_tapply_vsum(int Q, NumericVector x, IntegerVector dum){
	// Q: nber of classes
	// x: a matrix
	// dum: the N vector of clusters

	int N = x.size();

	NumericVector res(Q);
	int i, q;

	for(i=0 ; i<N ; i++){
		q = dum(i) - 1; // we take 1 off => different indexation in C
		res(q) += x(i);
	}

	return(res);
}

// similar a table but faster
// [[Rcpp::export]]
NumericVector cpp_table(int Q, IntegerVector dum){
	// Q: nber of classes
	// dum: the N vector of clusters

	int N = dum.size();

	NumericVector res(Q);
	int i, q;

	for(i=0 ; i<N ; i++){
		q = dum(i) - 1; // we take 1 off => different indexation in C
		res(q)++;
	}

	return(res);
}

// [[Rcpp::export]]
IntegerMatrix cpp_make_contrast(int N, int K, IntegerVector fact_num, bool addRef){
	// This function creates a NxK-1 matrix of contrasts
	// we take the first level as a reference

	if(addRef){
		IntegerMatrix res(N, K-1);

		for(int i=0 ; i<N ; i++){
			if(fact_num(i) != 1){
				res(i, fact_num(i) - 2) = 1;
			}
		}

		return(res);
	} else {
		IntegerMatrix res(N, K);

		for(int i=0 ; i<N ; i++){
			res(i, fact_num(i) - 1) = 1;
		}
		return(res);
	}

	// for the compiler to be happy
	IntegerMatrix res(1,1);
	return(res);
}

// [[Rcpp::export]]
NumericVector cpp_unik_id(IntegerVector x_sorted){
	// Q: the number of values of x_sorted
	// x_sorted: a sorted vector
	// => this function provides the indices of unique elements

	int N = x_sorted.size();

	int i, Q=1, q=0;

	// we get the number of elments of x_sorted: Q
	for(i=1 ; i<N ; i++){
		if(x_sorted(i-1) != x_sorted(i)){
			Q++;
		}
	}


	NumericVector res(Q);

	res[q] = 1; // the first element
	q++; // q est lindice courant

	for(i=1 ; i<N ; i++){
		if(x_sorted(i-1) != x_sorted(i)){
			res[q] = i+1;
			q++;
		}
	}

	return(res);
}


NumericVector cpp_poisson_closedForm(int q, int start, int K, NumericVector mu_in, NumericVector sum_y, IntegerVector dum_vect){
	// => computes the cluster coefs for the poisson

	int N = mu_in.size();

	int i, k;

	NumericVector cond_mu(K);

	for(i=0 ; i<N ; i++){
		k = dum_vect(q*N + i);
		cond_mu(k) += exp(mu_in(i));
	}

	NumericVector cluster_coef(K);

	for(k=0 ; k<K ; k++){
		cluster_coef(k) = log(sum_y(start + k)) - log(cond_mu(k));
	}

	return(cluster_coef);
}

NumericVector cpp_gaussian_closedForm(int q, int start, int K, NumericVector mu_in, NumericVector sum_y, IntegerVector dum_vect, NumericVector tableCluster_vect){
	// => computes the cluster coefs for the poisson
	// tableCluster_vec => as NumericVector to avoid problems with the olog function

	int N = mu_in.size();

	int i, k;
	NumericVector res(N);

	NumericVector cond_mu(K);

	for(i=0 ; i<N ; i++){
		k = dum_vect(q*N + i);
		cond_mu(k) += mu_in(i);
	}

	NumericVector cluster_coef(K);

	for(k=0 ; k<K ; k++){
		cluster_coef(k) = (sum_y(start + k) - cond_mu(k)) / tableCluster_vect(start + k);
	}

	return(cluster_coef);
}

NumericVector cpp_negbin_guessDummy(int K, int start, NumericVector sum_y, NumericVector tableCluster_vect, NumericVector mu_clust){
	// sum_y: sum of the lhs per cluster
	// tableCluster: nber of observation per cluster case => I put it as numeric and not integer to avoid overloading
	// mu_clust: value of mu for each cluster (usually min or max)

	NumericVector res(K);

	for(int k=0 ; k<K ; k++){
		res(k) = log(sum_y(start + k)) - log(tableCluster_vect(start + k)) - mu_clust(k);
	}

	return(res);
}

NumericVector cpp_negbin_clusterCoef(int family, int q, int K, int start, double theta, double epsDicho, NumericVector lhs, NumericVector sum_y, NumericVector mu, IntegerVector obsCluster_vect, IntegerVector start_vector, IntegerVector end_vector, IntegerVector dum_vect, NumericVector tableCluster_vect){

	// tablecluster_vec => put to numeric to avoid later problems

	int N = lhs.size();
	NumericVector cluster_coef(K);

	//
	// Obtaining min/max values pour les bornes
	//

	NumericVector mu_min(K), mu_max(K);
	IntegerVector ok(K); // init a 0 par defaut
	int i, k;

	for(i=0 ; i<N ; i++){
		k = dum_vect(q*N + i);

		if(ok[k] == 0){
			// means its the first time
			mu_min(k) = mu(i);
			mu_max(k) = mu(i);
			ok[k] = 1;
		} else {
			// the min
			if(mu(i) < mu_min(k)){
				mu_min(k) = mu(i);
			}

			// the max
			if(mu(i) > mu_max(k)){
				mu_max(k) = mu(i);
			}
		}
	}

	NumericVector borne_inf = cpp_negbin_guessDummy(K, start, sum_y, tableCluster_vect, mu_max);
	NumericVector borne_sup = cpp_negbin_guessDummy(K, start, sum_y, tableCluster_vect, mu_min);

	//
	// Computing the cluster coefficients
	//

	NumericVector clusterCoef = DichotomyNR_internal(family, K, start, theta, epsDicho, lhs, mu, borne_inf, borne_sup, obsCluster_vect, start_vector, end_vector);

	return(clusterCoef);
}

NumericVector cpp_logit_guessDummy(int K, int start, NumericVector sum_y, NumericVector tableCluster_vect, NumericVector mu_clust){
	// sum_y: sum of the lhs per cluster
	// tableCluster: nber of observation per cluster case => numeric to avoid overloading
	// mu_clust: value of mu for each cluster (usually min or max)
	// log(sum_y) - log(n_group-sum_y) - mu

	NumericVector res(K);

	for(int k=0 ; k<K ; k++){
		res(k) = log(sum_y(start + k)) - log(tableCluster_vect(start + k) - sum_y(start + k)) - mu_clust(k);
	}

	return(res);
}

NumericVector cpp_logit_clusterCoef(int family, int q, int K, int start, double theta, double epsDicho, NumericVector lhs, NumericVector sum_y, NumericVector mu, IntegerVector obsCluster_vect, IntegerVector start_vector, IntegerVector end_vector, IntegerVector dum_vect, NumericVector tableCluster_vect){

	// tableCluster_vec => put to NumericVector to avaoid later overloading problem with function log

	int N = lhs.size();
	NumericVector cluster_coef(K);

	//
	// Obtaining min/max values pour les bornes
	//

	NumericVector mu_min(K), mu_max(K);
	IntegerVector ok(K); // init a 0 par defaut
	int i, k;

	for(i=0 ; i<N ; i++){
		k = dum_vect(q*N + i);

		if(ok[k] == 0){
			// means its the first time
			mu_min(k) = mu(i);
			mu_max(k) = mu(i);
			ok[k] = 1;
		} else {
			// the min
			if(mu(i) < mu_min(k)){
				mu_min(k) = mu(i);
			}

			// the max
			if(mu(i) > mu_max(k)){
				mu_max(k) = mu(i);
			}
		}
	}

	NumericVector borne_inf = cpp_logit_guessDummy(K, start, sum_y, tableCluster_vect, mu_max);
	NumericVector borne_sup = cpp_logit_guessDummy(K, start, sum_y, tableCluster_vect, mu_min);

	//
	// Computing the cluster coefficients
	//

	NumericVector clusterCoef = DichotomyNR_internal(family, K, start, theta, epsDicho, lhs, mu, borne_inf, borne_sup, obsCluster_vect, start_vector, end_vector);

	return(clusterCoef);
}

// function that computes the clusters coefficents for a given cluster
NumericVector cpp_compute_single_cluster(int family, double theta, double epsDicho, int q, int start, int K, NumericVector mu_in, NumericVector lhs, NumericVector sum_y, NumericVector tableCluster_vect, IntegerVector dum_vect, IntegerVector obsCluster_vect, IntegerVector start_cluster_vect, IntegerVector end_cluster_vect){
	// tableCluster_vec put to NumericVector to avoid later overloading problem with the function log

	int N = lhs.size();
	NumericVector res(N), cluster_coef(K);

	if(family == 1){
		// Poisson
		cluster_coef = cpp_poisson_closedForm(q, start, K, mu_in, sum_y, dum_vect);
	} else if(family == 2){
		// negbin
		cluster_coef = cpp_negbin_clusterCoef(family, q, K, start, theta, epsDicho, lhs, sum_y, mu_in, obsCluster_vect, start_cluster_vect, end_cluster_vect, dum_vect, tableCluster_vect);
	} else if(family == 3){
		// Logit
		cluster_coef = cpp_logit_clusterCoef(family, q, K, start, theta, epsDicho, lhs, sum_y, mu_in, obsCluster_vect, start_cluster_vect, end_cluster_vect, dum_vect, tableCluster_vect);
	} else if(family == 4){
		// gaussian
		cluster_coef = cpp_gaussian_closedForm(q, start, K, mu_in, sum_y, dum_vect, tableCluster_vect);
	}

	return(cluster_coef);
}


// [[Rcpp::export]]
NumericVector Rcpp_compute_sum_clusters(int N, int Q, int family, double theta, double epsMu, NumericVector init, NumericVector lhs, NumericVector sum_y, NumericVector mu, IntegerVector dum_vect, IntegerVector nbCluster, NumericVector tableCluster_vect, IntegerVector obsCluster_vect, IntegerVector start_cluster_vect, IntegerVector end_cluster_vect){
	// computes the sum of clusters iteratively
	// sum_y: vector of length sum_cases
	// dum_vect: N*Q stacked vector of the clusters ID, is c(dum[[1]], ..., dum[[q]], ..., dum[[Q]])
	// tableCluster_vect: sum_cases stacked vector of the number of obs per cluster category // as numeric vector to avoid later problem with overloading in the function log
	// obsCluster_vect: N*Q stacked vector of c(order(dum[[1]]), ..., order(dum[[q]]), ..., order(dum[[Q]]))
	// start_cluster_vect (end_cluster_vect): the IDs of the observation of obsCluster_vect for each cluster case (used only in dichotomy, and obtained using function RcppCreate_start_end_indexes)

	int q, i, k, K;
	bool ok = true;
	int iter=0, iterMax = 100;
	double epsDicho = epsMu/10;
	NumericVector mu_in(N);

	// initialisation
	for(i=0 ; i<N ; i++){
		mu_in(i) = mu(i) + init(i);
	}

	//
	// Preliminary information
	//

	IntegerVector start(Q);

	for(q=0 ; q<Q ; q++){
		if(q == 0){
			start(q) = 0;
		} else {
			start(q) = start(q-1) + nbCluster(q-1);
		}
	}

	//
	// The main loop
	//

	while(ok && iter<iterMax){
		iter++;
		ok = false;
		// Rprintf("Iter: %i\n", iter);

		for(q=0 ; q<Q ; q++){
			// Rprintf("q=%i", q);

			K = nbCluster(q);
			NumericVector cluster_coef = cpp_compute_single_cluster(family, theta, epsDicho, q, start(q), K, mu_in, lhs, sum_y, tableCluster_vect, dum_vect, obsCluster_vect, start_cluster_vect, end_cluster_vect);

			// Rprintf(" ok\n");

			// control
			for(k=0 ; k<K ; k++){
				if(cpp_abs(cluster_coef(k)) > epsMu){
					ok = true;
				}
			}

			// now the full vector of cluster coefficients
			for(i=0 ; i<N ; i++){
				k = dum_vect(q*N + i);
				mu_in(i) += cluster_coef(k);
			}
		}
	}

	// Rprintf("iter=%i\n", iter);

	if(iter == iterMax){
		Rprintf("[Getting dummies] maximum iterations reached (%i)\n", iterMax);
	}

	return(mu_in);
}




// [[Rcpp::export]]
IntegerVector Rcpp_unclassFactor(NumericVector x){
	// x: a sorted integer vector

	int N = x.size();

	IntegerVector res(N);
	int k=1;
	res[0] = 1;

	for(int i=1 ; i<N ; i++){
		if(x(i-1)!=x(i)) k++;
		res[i] = k;
	}

	return(res);
}


// [[Rcpp::export]]
NumericVector Rcpp_Ax_Ais1(IntegerMatrix dumMat, IntegerVector nbCluster, IntegerVector tableCluster_all, NumericVector x){
	// This functions just computes A%*%x, but in fact it creates the (large) matrix A
	// by using the cluster IDs (dumMat)
	// takes in:
	// dumMat: the matrix of dummies (n X c) each obs => cluster
	// nbCluster: size of each cluster
	// x: the vector to be multiplied

	int Q = dumMat.ncol(), N = dumMat.nrow();

	int i, q, q_bis;
	int index, index_bis;
	int sum_cases=0;
	IntegerVector start(Q), end(Q);

	for(q=0 ; q<Q ; q++){
		sum_cases += nbCluster(q);
		if(q == 0){
			start(q) = 0;
			end(q) = nbCluster(q);
		} else {
			start(q) = start(q-1) + nbCluster(q-1);
			end(q) = end(q-1) + nbCluster(q);
		}
	}

	// the output
	NumericVector res(sum_cases);
	// Initialisation
	for(i=0 ; i<sum_cases ; i++){
		res[i] = x[i];
	}

	// The multiplication
	for(q=0 ; q<Q ; q++){
		for(i=0 ; i<N ; i++){
			// we get the index
			index = start[q] + dumMat(i, q);

			// we sum the elements
			for(q_bis=0 ; q_bis<Q ; q_bis++){
				if(q_bis != q){
					index_bis = start[q_bis] + dumMat(i, q_bis);
					res[index] += x(index_bis) / tableCluster_all(index);
				}
			}

		}
	}

	return(res);
}



//
// New version of the dichotomy: using the exp of the lcuster coefficient
//


// [[Rcpp::export]]
double new_cpp_NB_dum_fx(double theta, NumericVector lhs, NumericVector exp_mu, double x1, IntegerVector obsCluster, int start, int end){

	double res=0;
	int i, obs;
	double val;

	for(i=start ; i<end ; i++){

		obs = obsCluster(i); // obscluster must be in cpp index

		val = exp(x1)*exp_mu(obs);

		res += lhs(obs) - (theta + lhs(obs)) * val / (val + theta);
	}

	return(res);
}

// Function to get the evaluation of the dummy for the Logit
double new_cpp_LOGIT_dum_fx(NumericVector lhs, NumericVector exp_mu, double x1, IntegerVector obsCluster, int start, int end){

	double res=0;
	int i, obs;
	double val;

	for(i=start ; i<end ; i++){

		obs = obsCluster(i); // obscluster must be in cpp index

		val = exp(x1)*exp_mu(obs);

		res += lhs(obs) - val / (1 + val);
	}

	return(res);
}

// Function to get the evaluation of the dummy
double new_cpp_dum_fx(int family, double theta, NumericVector lhs, NumericVector exp_mu, double x1, IntegerVector obsCluster, int start, int end){

	double res;

	if(family == 2){
		// => the negative binomial
		res = new_cpp_NB_dum_fx(theta, lhs, exp_mu, x1, obsCluster, start, end);
	} else if(family == 3){
		// => Logit
		res = new_cpp_LOGIT_dum_fx(lhs, exp_mu, x1, obsCluster, start, end);
	} else {
		stop("Unrecognized family");
	}

	return(res);
}

// [[Rcpp::export]]
double new_cpp_NB_dum_dfx(double theta, NumericVector lhs, NumericVector exp_mu, double x1, IntegerVector obsCluster, int start, int end){

	double res=0;
	int i, obs;
	double val, num;

	for(i=start ; i<end ; i++){

		obs = obsCluster(i); // obscluster must be in cpp index

		num = exp(x1)*exp_mu(obs);
		val = theta + num;

		res += - theta * (theta + lhs(obs)) * num / val / val;
	}

	return(res);
}

// Function to get the evaluation of the dummy for the Logit
double new_cpp_LOGIT_dum_dfx(NumericVector lhs, NumericVector exp_mu, double x1, IntegerVector obsCluster, int start, int end){

	double res=0;
	int i, obs;
	double val, num;

	for(i=start ; i<end ; i++){

		obs = obsCluster(i); // obscluster must be in cpp index

		num = exp(x1)*exp_mu(obs);
		val = 1 + num;

		res += - num / val / val;
	}

	return(res);
}

// Function to get the evaluation of the dummy
double new_cpp_dum_dfx(int family, double theta, NumericVector lhs, NumericVector exp_mu, double x1, IntegerVector obsCluster, int start, int end){

	double res;

	if(family == 2){
		// => the negative binomial
		res = new_cpp_NB_dum_dfx(theta, lhs, exp_mu, x1, obsCluster, start, end);
	} else if(family == 3){
		// => Logit
		res = new_cpp_LOGIT_dum_dfx(lhs, exp_mu, x1, obsCluster, start, end);
	} else {
		stop("Unrecognized family");
	}

	return(res);
}



// [[Rcpp::export]]
NumericVector new_RcppDichotomyNR(int N, int K, int family, double theta, double epsDicho, NumericVector lhs, NumericVector exp_mu, NumericVector borne_inf, NumericVector borne_sup, IntegerVector obsCluster, IntegerVector tableCluster){
	// ARGUMENTS:
	// N: number of observations
	// K: number of clusters classes
	// family:
	// 	- 1: negative binomial
	// 	- 2: logit
	// theta: used only for the NB, but necessary
	// epsDicho: the stopping criterion
	// lhs: the dependent variable
	// exp_mu: the vector of the value of exp_mu
	// borne_inf, init, borne_sup => their length is equal to the number of classes
	// obsCluster: the integer vector that is equal to order(dum[[g]])
	// tableCluster: it gives for each classes the number of element

	// INFO:
	// this algorithm applies dichotomy to find out the cluster coefficients
	// this is for the logit and for the negative binomial

	int itermax = 200, iter;
	double value, x0, derivee;; // evaluation of the dummy
	int k; // classical index
	int start, end; // the indices used to compute sums wrt clusters
	NumericVector res(K);
	double lower_bound, upper_bound;

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


	for(int k=0 ; k<K ; k++){
		// we loop over each cluster

		double x1 = 0; // we initialise the cluster coefficient at 0 (it should converge to 0)
		bool ok = true;
		iter = 0;

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
			R_CheckUserInterrupt();
			iter++;

			// 1st step: initialisation des bornes
			value = new_cpp_dum_fx(family, theta, lhs, exp_mu, x1, obsCluster, start, end);

			// update of the bounds.
			if(value > 0){
				lower_bound = x1;
			} else {
				upper_bound = x1;
			}

			// 2nd step: NR iteration
			x0 = x1;
			derivee = new_cpp_dum_dfx(family, theta, lhs, exp_mu, x1, obsCluster, start, end);
			x1 = x0 - value / derivee;
			// Rprintf("x1: %5.2f\n", x1);

			// 3rd step: dichotomy (if necessary)
			// Update of the value if it goes out of the boundaries
			if(x1 >= upper_bound || x1 <= lower_bound){
				x1 = (lower_bound + upper_bound)/2;
			}

			// the stopping criteria
			if(iter == itermax){
				ok = false;
				Rprintf("[Getting cluster coefficients] max iterations reached (%i).\n", itermax);
				Rprintf("Value Sum Deriv = %f.\n", value);
			}
			if(cpp_abs(x0 - x1) < epsDicho){
				ok = false;
			}
		}
		// Rprintf("iter=%i.\n", iter);
		res[k] = x1;
	}

	return(res);
}


// [[Rcpp::export]]
List sum_double_index(int n_i, int n_j, IntegerVector index_i, IntegerVector index_j, NumericVector x){
	// ARGUMENTS:
	// n_i (= length(unique(index_i)))
	// index_i, index_j, two indexes that MUST be sorted before hand
	// x and the other mindexes must be of the same length
	// WHAT:
	// it sums all x per combinaison of indexes


	int n = x.length();

	// we first put the results into a list
	int i;
	int index_current = 0;
	double value = 0;

	NumericVector col_1(n), col_2(n), col_value(n);

	value = x[0];

	for(i=1 ; i<n ; i++){
		if(index_j[i] != index_j[i-1] || index_i[i] != index_i[i-1]){
			// save the value if we change index
			col_1(index_current) = index_i[i-1];
			col_2(index_current) = index_j[i-1];
			col_value(index_current) = value;

			// new row
			index_current++;

			// the new value
			value = x[i];
		} else {
			value += x[i];
		}
	}

	// last save
	col_1(index_current) = index_i[i-1];
	col_2(index_current) = index_j[i-1];
	col_value(index_current) = value;

	// we fill the matrix
	List res;

	int n_row = index_current + 1;

	IntegerVector index_i_out(n_row), index_j_out(n_row);
	NumericVector sum_x(n_row);

	for(i=0 ; i < n_row ; i++){
		index_i_out(i) = col_1(i) - 1;
		index_j_out(i) = col_2(i) - 1;
		sum_x(i) = col_value(i);
	}

	res["index_i"] = index_i_out;
	res["index_j"] = index_j_out;
	res["coefmat"] = sum_x;
	res["n_i"] = n_i;
	res["n_j"] = n_j;

	return(res);
}


// [[Rcpp::export]]
NumericVector matmult(IntegerVector index_i, IntegerVector index_j, NumericVector matcoef, NumericVector x){
	// suppprimer => use mmmult instead, this one is SLOW

	int n_obs = matcoef.length();
	int n = index_i[n_obs - 1] + 1;

	int obs;

	NumericVector res(n);

	for(obs=0 ; obs<n_obs ; obs++){
		res[index_i[obs]] += matcoef[obs] * x[index_j[obs]];
	}

	return(res);
}

// [[Rcpp::export]]
SEXP mmult(int n, SEXP index_i, SEXP index_j, SEXP coefmat, SEXP x){
	// I don't know why: lots of overheads in Rcpp (~20ms)
	// this way is MUCH faster
	// ARGS:
	// index_i, index_j: must be cpp indexes

	int n_obs = Rf_length(coefmat);
	int obs;

	// pointers
	double *pcoefmat, *px, *pres;
	int *pindex_i, *pindex_j;

	SEXP res = PROTECT(Rf_allocVector(REALSXP, n));

	pcoefmat = REAL(coefmat);
	px = REAL(x);
	pres = REAL(res);

	pindex_i = INTEGER(index_i);
	pindex_j = INTEGER(index_j);

	// double mysum = 0;

	// init
	for(obs=0 ; obs<n ; obs++){
		pres[obs] = 0;
	}

	for(obs=0 ; obs<n_obs ; obs++){
		pres[pindex_i[obs]] += pcoefmat[obs] * px[pindex_j[obs]];
		// if(pindex_i[obs] == 1){
		// 	int k = pindex_j[obs];
		// 	mysum += px[k] * pcoefmat[obs];
		// 	Rprintf("k = %i, coefmat = %3.2f, px = %3.2f, sum = %3.2f, res = %3.2f\n", k, pcoefmat[obs], px[k], mysum, pres[pindex_i[obs]]);
		// }
	}

	UNPROTECT(1);

	return(res);
}


////
//// New convergence functions ////
////


// [[Rcpp::export]]
List RcppPartialDerivative_new(int iterMax, int Q, int N, int K, double epsDeriv, NumericVector ll_d2,	NumericMatrix F, NumericVector init, IntegerMatrix dumMat, IntegerVector nbCluster){
	// takes in:
	// dumMat: the matrix of dummies (n X c) each obs => cluster
	// init: the initialisation of the sum of derivatives vector
	// ll_d2: the second derivative
	// F: the jacobian matrix

	int iter, max_iter = 1;

	int i, q, k, c;
	int index;
	int sum_cases = 0;
	bool ok;
	double new_value;
	IntegerVector start(Q), end(Q);

	for(q=0 ; q<Q ; q++){
		sum_cases += nbCluster(q);
		if(q == 0){
			start(q) = 0;
			end(q) = nbCluster(q);
		} else {
			start(q) = start(q-1) + nbCluster(q-1);
			end(q) = end(q-1) + nbCluster(q);
		}
	}

	NumericMatrix clusterDeriv(sum_cases, K); // the output
	NumericVector sum_lld2(sum_cases);

	// Creation of the sum_lld2
	for(i=0 ; i<N ; i++){
		for(q=0 ; q<Q ; q++){
			// index = start[q] + dumMat(i, q) - 1;
			index = start[q] + dumMat(i, q);
			sum_lld2[index] += ll_d2(i);
		}
	}

	// the result to return
	NumericMatrix S(N, K);
	for(i=0 ; i<N ; i++){
		for(k=0 ; k<K ; k++){
			S(i,k) = init(i,k);
		}
	}

	for(k=0 ; k<K ; k++){

		ok = true;
		iter = 0;
		while( ok & (iter<iterMax) ){
			iter++;
			ok = false;

			for(q=0 ; q<Q ; q++){
				R_CheckUserInterrupt();

				// init of the vector
				for(c=start[q] ; c<end[q] ; c++){
					clusterDeriv(c, k) = 0;
				}

				for(i=0 ; i<N ; i++){
					// index = start[q] + dumMat(i, q) - 1;
					index = start[q] + dumMat(i, q);
					clusterDeriv(index, k) += (F(i,k) + S(i,k))*ll_d2(i);
				}

				// on finit de calculer clusterDeriv + controle
				for(c=start[q] ; c<end[q] ; c++){
					new_value = -clusterDeriv(c, k)/sum_lld2[c];
					clusterDeriv(c, k) = new_value;
					if(cpp_abs(new_value) > epsDeriv){
						ok = true;
					}
				}

				// on ajoute la derivee a S:
				for(i=0 ; i<N ; i++){
					// index = start[q] + dumMat(i, q) - 1;
					index = start[q] + dumMat(i, q);
					S(i,k) += clusterDeriv(index, k);
				}

			}
		}

		if(iter > max_iter){
			max_iter = iter;
		}

	}

	List res;
	res["dxi_dbeta"] = S;
	res["iter"] = max_iter;

	return(res);
}

// [[Rcpp::export]]
List RcppPartialDerivative_gaussian_new(int iterMax, int Q, int N, int K, double epsDeriv, NumericMatrix F, NumericVector init, IntegerMatrix dumMat, IntegerVector nbCluster){
	// takes in:
	// dumMat: the matrix of dummies (n X c) each obs => cluster
	// init: the initialisation of the sum of derivatives vector
	// ll_d2: the second derivative
	// F: the jacobian matrix

	int iter, max_iter = 1;

	int i, q, k, c;
	int index;
	int sum_cases = 0;
	bool ok;
	double new_value;
	IntegerVector start(Q), end(Q);

	for(q=0 ; q<Q ; q++){
		sum_cases += nbCluster(q);
		if(q == 0){
			start(q) = 0;
			end(q) = nbCluster(q);
		} else {
			start(q) = start(q-1) + nbCluster(q-1);
			end(q) = end(q-1) + nbCluster(q);
		}
	}

	NumericMatrix clusterDeriv(sum_cases, K); // the output

	NumericVector n_group(sum_cases);

	// Creation of the sum_lld2
	for(i=0 ; i<N ; i++){
		for(q=0 ; q<Q ; q++){
			// index = start[q] + dumMat(i, q) - 1;
			index = start[q] + dumMat(i, q);
			n_group[index] ++;
		}
	}

	// the result to return
	NumericMatrix S(N, K);
	for(i=0 ; i<N ; i++){
		for(k=0 ; k<K ; k++){
			S(i,k) = init(i,k);
		}
	}

	for(k=0 ; k<K ; k++){

		ok = true;
		iter = 0;
		while( ok & (iter<iterMax) ){
			iter++;
			ok = false;

			for(q=0 ; q<Q ; q++){
				R_CheckUserInterrupt();

				// init of the vector
				for(c=start[q] ; c<end[q] ; c++){
					clusterDeriv(c, k) = 0;
				}

				for(i=0 ; i<N ; i++){
					// index = start[q] + dumMat(i, q) - 1;
					index = start[q] + dumMat(i, q);
					clusterDeriv(index, k) += F(i,k) + S(i,k);
				}

				// on finit de calculer clusterDeriv + controle
				for(c=start[q] ; c<end[q] ; c++){
					new_value = -clusterDeriv(c, k)/n_group[c];
					clusterDeriv(c, k) = new_value;
					if(cpp_abs(new_value) > epsDeriv/Q){
						ok = true;
					}
				}

				// on ajoute la derivee a S:
				for(i=0 ; i<N ; i++){
					// index = start[q] + dumMat(i, q) - 1;
					index = start[q] + dumMat(i, q);
					S(i,k) += clusterDeriv(index, k);
				}

			}
		}

		if(iter > max_iter){
			max_iter = iter;
		}

	}

	List res;
	res["dxi_dbeta"] = S;
	res["iter"] = max_iter;

	return(res);
}





