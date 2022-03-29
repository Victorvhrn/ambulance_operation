#include "../include/generator.h"

using namespace std;

Generator::Generator(string calls_path, string neighbors_path, string info_path){
	auto info_arq = ifstream(info_path, ios::in);
	info_arq >> T >> G >> R >> P >> nb_land_types >> nb_holidays_years; 
	slot_duration = 24 / T;
	daily_obs = std::vector<int>(G, 0);
	for(int g = 0; g < G; ++g){
		info_arq >> daily_obs[g];
	}
	info_arq.close();

	duration = vector<double>(T, slot_duration);
	xt::xarray<double>::shape_type shape = {T,G,R,P};
	nb_observations = xt::zeros<int>(shape);
	nb_calls = xt::zeros<int>(shape);
	lambda_teorico = xt::zeros<double>(shape);

	int max_obs = *max_element(daily_obs.begin(), daily_obs.end());
	
	shape = {T,G,R,P,static_cast<ulong>(max_obs)};
	sample_calls = xt::zeros<int>(shape);

	auto calls_arq = ifstream(calls_path, ios::in);
	std::string aux_str;
	do{
		std::getline(calls_arq, aux_str);
		if(aux_str == "END"){
			break;
		}
		std::istringstream ss(aux_str);
		int t, g, r, p, j, val, h;
		ss >> t >> g >> r >> p >> j >> val >> h;
		sample_calls(t,g,r,p,j) = val;
		nb_calls(t,g,r,p) += sample_calls(t,g,r,p,j);
		nb_observations(t,g,r,p) = daily_obs[g];
	}while(true);
	calls_arq.close();

	xt::xarray<double>::shape_type shape_gt = {G,T};
	xt::xarray<double>::shape_type shape_tgrp = {T,G,R,P};
	lambda_teorico_agg = xt::zeros<double>(shape_gt);
	calls_agg = xt::zeros<double>(shape_gt);
	ofstream lambda_out(fmt::format("{}/lambda_max_likelihood.txt", 
		g_params.generator_folder), ios::out);
	ofstream calls_out(fmt::format("{}/calls_agregados.txt", g_params.generator_folder), 
		ios::out);
	ofstream empirical_out(fmt::format("{}/empirical_estimation.txt", g_params.generator_folder),
		ios::out);
	for(int g = 0; g < G; ++g){
		lambda_out << g << " ";
		calls_out << g << " ";
		for(int t = 0; t < T; ++t){
			for(int r = 0; r < R; ++r){
				for(int p = 0; p < P; ++p){
					lambda_teorico(t,g,r,p) = static_cast<double>(nb_calls(t,g,r,p))/
						daily_obs[g];
					empirical_out << t << " " << g << " " << r << " " << p << " ";
					empirical_out << lambda_teorico(t,g,r,p) << "\n";
					lambda_teorico_agg(g,t) += lambda_teorico(t,g,r,p);
				}
			}
			calls_agg(g,t) = lambda_teorico_agg(g,t)*daily_obs[g];
			lambda_out << lambda_teorico_agg(g,t) << " ";
			calls_out << calls_agg(g,t) << " ";
		}
		lambda_out << "\n";
		calls_out << "\n";
	}

	lambda_out.close();
	calls_out.close();


	neighbors = std::vector<vector<int>>(R, std::vector<int>());
	shape = {R, R};
	regions = std::vector<Location>(R, null_location);
	distance = xt::ones<double>(shape);
	type = std::vector<int>(R,0);
	shape = {R, nb_land_types};
	regressors = xt::zeros<double>(shape);
	auto neighbors_arq = ifstream(neighbors_path, ios::in);
	while(true){
		int ind, terrain_type, s; 
		double lat, longi, dist;
		std::getline(neighbors_arq, aux_str);
		if(aux_str == "END"){
			break;
		}
		std::istringstream ss(aux_str);
		ss >> ind >> lat >> longi >> terrain_type;
		type[ind] = terrain_type;
		regions[ind] = make_pair(lat, longi);
		for(int j = 0; j < nb_land_types; ++j){
			ss >> regressors(ind,j);
		}
		while(ss >> s >> dist){
			distance(ind,s) = dist;
			neighbors[ind].push_back(s);
		}
	}
	neighbors_arq.close();
	std::cout << "Initialized\n";
}

double Generator::cross_validation(xt::xarray<double>& x, double sigma, double beta_bar, 
		double proportion){
	std::vector<double> alphas{0.01,0.1,0.5,1,10,100,1000, pow(10,7)};
	// std::vector<double> alphas{0.001,0.5,1};
	double max_likelihood = GRB_INFINITY;
	double best_alpha = -1;
	int nb_obs = static_cast<int>(sample_calls.shape(4));
	int nb_in_block = nb_obs*proportion;
	xt::xarray<double>::shape_type shape = {T,G,R,P};

	for(double alpha: alphas){
		double likelihood = 0;
		for(int ind_cross = 0; ind_cross < floor(1/proportion); ++ind_cross){
			nb_observations = xt::zeros<int>(shape);
			nb_calls = xt::zeros<int>(shape);
			for(int ind = ind_cross*nb_in_block; ind < (ind_cross+1)*nb_in_block; ++ind){
				for(int t = 0; t < T; ++t){
					for(int g = 0; g < G; ++g){
						for(int r = 0; r < R; ++r){
							for(int p = 0; p < P; ++p){
								if(ind <= daily_obs[g]){
									nb_observations(t,g,r,p) += 1;
									nb_calls(t,g,r,p) += sample_calls(t,g,r,p,ind);
								}
							}
						}
					}
				}
			}
			x = epsilon*xt::ones<double>(shape);
			auto f_val = projected_gradient_armijo_feasible(x, alpha, sigma);
			xt::xarray<int> nb_calls_remaining = xt::zeros<int>(shape);

			for(int ind = 0; ind < ind_cross*nb_in_block; ++ind){
				for(int t = 0; t < T; ++t){
					for(int g = 0; g < G; ++g){
						for(int r = 0; r < R; ++r){
							for(int p = 0; p < P; ++p){
								nb_calls_remaining(t,g,r,p) += sample_calls(t,g,r,p,ind);
							}
						}
					}
				}
			}

			for(int ind = (ind_cross+1)*nb_in_block; ind < nb_obs; ++ind){
				for(int t = 0; t < T; ++t){
					for(int g = 0; g < G; ++g){
						for(int r = 0; r < R; ++r){
							for(int p = 0; p < P; ++p){
								if(ind <= daily_obs[g]){
									nb_calls_remaining(t,g,r,p) += sample_calls(t,g,r,p,ind);
								}
							}
						}
					}
				}
			}

			double f = 0;
			for(int t = 0; t < T; ++t){
				for(int g = 0; g < G; ++g){
					for(int r = 0; r < R; ++r){
						for(int p = 0; p < P; ++p){
							double current_lambda = x(t,g,r,p);
							f += (daily_obs[g] - nb_observations(t,g,r,p))*current_lambda - 
								nb_calls_remaining(t,g,r,p)*log(current_lambda);
						}
					}
				}
			}
			likelihood += f;
		}
		fmt::print("Alpha = {}, likelihood = {}\n", alpha, likelihood);
		if(likelihood < max_likelihood){
			max_likelihood = likelihood;
			best_alpha = alpha;
		}
	}
	
	nb_observations = xt::zeros<int>(shape);
	nb_calls = xt::zeros<int>(shape);
	for(int t = 0; t < T; ++t){
		for(int g = 0; g < G; ++g){
			for(int r = 0; r < R; ++r){
				for(int p = 0; p < P; ++p){
					for(int ind = 0; ind < daily_obs[g]; ++ind){
						nb_observations(t,g,r,p) += 1;
						nb_calls(t,g,r,p) += sample_calls(t,g,r,p,ind);

					}
				}
			}
		}
	}
	auto f_val_best = projected_gradient_armijo_feasible(x, best_alpha, sigma);
	for(int t = 0; t < T; ++t){
		for(int g = 0; g < G; ++g){
			for(int r = 0; r < R; ++r){
				for(int p = 0; p < P; ++p){
					if(x(t,g,r,p) < epsilon){
						x(t,g,r,p) = 0;
					}
				}
			}
		}
	}

	return best_alpha;
}


void Generator::test(){
	epsilon = 0.01;
	xt::xarray<double>::shape_type shape = {T,G,R,P};
	sigma = 0.5;
	double beta_bar = 1;
	max_iter = 10;
	double alpha = 0;

	double beta_tilde = 1;
	double beta_hat = 2;

	xt::xarray<double> x = epsilon*xt::ones<double>(shape);
	auto f_val1 = projected_gradient_armijo_feasible(x, alpha, sigma);

	write_params(x, alpha);
	fmt::print("END Alpha = {}\n", alpha);

	double proportion = 0.2;
	x = epsilon*xt::ones<double>(shape);
	double best_alpha = cross_validation(x,sigma, beta_bar, proportion);
	write_params(x,best_alpha);
	fmt::print("END Alpha = {}\n", best_alpha);

	alpha = pow(10,7);

	x = epsilon*xt::ones<double>(shape);
	auto f_val2 = projected_gradient_armijo_feasible(x, alpha, sigma);
	write_params(x, alpha);
	fmt::print("END Alpha = {}\n", alpha);
}


void Generator::write_params(xt::xarray<double>& x, double alpha){
	ofstream out_file(fmt::format("{}/xNonRegressorT{}G{}I{}P{}_alpha{}.txt",
		g_params.generator_folder, T,G,R,P,alpha), ios::out);

	for(int t = 0; t < T; ++t){
		for(int g = 0; g < G; ++g){
			for(int r = 0; r < R; ++r){
				for(int p = 0; p < P; ++p){
					out_file << t << " " << g << " " << r << " " << p << " ";
					out_file << x(t,g,r,p) << "\n";
				}
			}
		}
	}
	out_file << "END";
	out_file.close();

	std::vector<std::vector<double>> sum_by_gt(G, std::vector<double>(T, 0));
	for(int g = 0; g < G; ++g){
		for(int t = 0; t < T; ++t){
			double sum = 0;
			for(int r = 0; r < R; ++r){
				for(int p = 0; p < P; ++p){
					sum += x(t,g,r,p);
				}
			}
			sum_by_gt[g][t] = sum;
		}
	}

	ofstream plot(fmt::format("{}/params_plot_{}.txt",g_params.generator_folder, alpha), ios::out);
	for(int g = 0; g < G; ++g){
		plot << g << " ";
		fmt::print("{} ",g);
		for(int t = 0; t < T; ++t){
			plot << sum_by_gt[g][t] << " ";
			fmt::print("{} ", sum_by_gt[g][t]);
		}
		plot << "\n";
		fmt::print("\n");
	}
	plot.close();
}


std::vector<double> Generator::projected_gradient_armijo_boundary(xt::xarray<double>& x, 
	double alpha, double sigma, double beta_bar){
	xt::xarray<double>::shape_type shape = {T,G,R,P};
	xt::xarray<double> z = xt::zeros<double>(shape);
	int k = 0;
	std::vector<double> f_val;
	f_val.reserve(max_iter);
	while(k < max_iter){
		double fold = oracle_objective_model1(x, alpha);
		auto gradient = oracle_gradient_model1(x, alpha);
		bool stop = false;
		int j = 0;
		double f = 0;
		do{
			xt::xarray<double> diff = x - (beta_bar/pow(2,j))*gradient;
			comp_wise_max(z, diff, epsilon);
			
			f = oracle_objective_model1(z, alpha);
			double rhs = xt::sum(gradient*(x-z))();
			if(f <= fold- sigma*rhs){
				stop = true;
			}else{
				++j;
			}
		}while(!stop);
		// fmt::print("k = {}, f = {}\n", k,f);
		f_val.push_back(f);
		x = z;
		++k;
	}

	return f_val;
}

std::vector<double> Generator::projected_gradient_armijo_feasible(xt::xarray<double>& x, 
		double alpha, double sigma){
	xt::xarray<double>::shape_type shape = {T,G,R,P};
	xt::xarray<double> z = xt::zeros<double>(shape);
	int k = 0;
	std::vector<double> f_val;
	f_val.reserve(max_iter);
	double b_param = 2;
	double beta_k = b_param;

	while(k < max_iter){
		double fold = oracle_objective_model1(x, alpha);
		xt::xarray<double> gradient = oracle_gradient_model1(x, alpha);
		xt::xarray<double> diff = x - beta_k*gradient;
		comp_wise_max(z, diff, epsilon);

		bool stop = false;
		int j = 0;
		double f = 0;
		double rhs = xt::sum(gradient*(x-z))();
		xt::xarray<double> z_aux = xt::zeros<double>(x.shape());
		do{
			z_aux = x + (1/pow(2,j))*(z-x);
			f = oracle_objective_model1(z_aux, alpha);
			if(f <= fold-(sigma/pow(2,j))*rhs){
				stop = true;
			}else{
				++j;
			}
		}while(!stop);
		f_val.push_back(f);
		x = z_aux;
		beta_k = b_param/pow(2,j);
		++k;
	}

	return f_val;
}

xt::xarray<double> Generator::oracle_gradient_model1(xt::xarray<double>& lambda, 
	double alpha){
	xt::xarray<double>::shape_type shape = {T,G,R,P};
	xt::xarray<double> gradient = xt::zeros<double>(shape);

	for(int t = 0; t < T; ++t){
		for(int g = 0; g < G; ++g){
			for(int r = 0; r < R; ++r){
				for(int p = 0; p < P; ++p){
					double current_lambda = lambda(t,g,r,p);
					double grad_component = nb_observations(t,g,r,p) - 
						(nb_calls(t,g,r,p) / current_lambda);

					for(auto s: neighbors[r]){
						if(type[r] == type[s]){
							grad_component += 2*alpha*(lambda(t,g,r,p) - 
								lambda(t,g,s,p)) / distance(r,s);
						}
					}
					gradient(t,g,r,p) = grad_component;
				}
			}
		}
	}

	return gradient;
}


double Generator::oracle_objective_model1(xt::xarray<double>& lambda,  double alpha){
	double f = 0;
	for(int t = 0; t < T; ++t){
		for(int g = 0; g < G; ++g){
			for(int r = 0; r < R; ++r){
				for(int p = 0; p < P; ++p){
					double current_lambda = lambda(t,g,r,p);
					f += nb_observations(t,g,r,p)*current_lambda-nb_calls(t,g,r,p)*
						log(current_lambda);
					for(auto s: neighbors[r]){
						if(type[r] == type[s]){
							f += (0.5*alpha)*(pow(current_lambda - 
								lambda(t,g,s,p), 2)) / distance(r,s);
						}
					}
				}
			}
		}
	}
	return f;
}


double Generator::average_difference(xt::xarray<double>& x){
	int count = T*G*R*P;
	double sum_diff = 0;
	double min = GRB_INFINITY;
	double max = -GRB_INFINITY;
	for(int t = 0; t < T; ++t){
		for(int g = 0; g < G; ++g){
			for(int r = 0; r < R; ++r){
				for(int p = 0; p < P; ++p){
					sum_diff += abs(lambda_teorico(t,g,r,p) - x(t,g,r,p));
					if(x(t,g,r,p) < min){
						min = x(t,g,r,p);
					}

					if(x(t,g,r,p) > max){
						max = x(t,g,r,p);
					}
				}	
			}
		}
	}
	fmt::print("Min x = {}, Max x = {}\n",min, max);
	return sum_diff / count;
}


double Generator::max_difference(xt::xarray<double>& x){
	int count = T*G*R*P;
	double max_diff = -GRB_INFINITY;
	for(int t = 0; t < T; ++t){
		for(int g = 0; g < G; ++g){
			for(int r = 0; r < R; ++r){
				for(int p = 0; p < P; ++p){
					double val = abs(lambda_teorico(t,g,r,p) - x(t,g,r,p));
					if(val > max_diff){
						max_diff = val;
					}
				}	
			}
		}
	}	

	return max_diff;
}

void Generator::comp_wise_max(xt::xarray<double>& z ,xt::xarray<double>& a, double eps){
	for(int t = 0; t < T; ++t){
		for(int g = 0; g < G; ++g){
			for(int r = 0; r < R; ++r){
				for(int p = 0; p < P; ++p){
					z(t,g,r,p) = max(a(t,g,r,p), eps);
				}
			}
		}
	}
}


bool Generator::is_neighbor(int r, int s){
	return r != s;
}



// std::vector<double> absc(n_x, 0);
// std::vector<double> ord(n_y, 0);

// for(int i = 0; i < n_x; ++i){
// 	absc[i] = x_max/(2*n_x) + (x_max/n_x)*i;
// }

// for(int j = 0; j < n_y; ++j){
// 	ord[j] = y_max/(2*n_y) + (y_max/n_y)*j;
// }

// type = std::vector<int>(R,-1);
// for(int i = 1; i <= R; ++i){
// 	type[i-1] = 1;
// 	int x_i = i % n_x;
// 	int y_i;
// 	if(x_i == 0){
// 		y_i = i / n_x;
// 	}else{
// 		y_i = 1 + (i-x_i) / n_x;
// 	}
// 	if(x_i == 0){
// 		x_i = n_x;
// 	}
// 	//setting neighbours
// 	// Possible neighbors are (Xi,Yi-1),(Xi,Yi+1),
//  	// (Xi-1,Yi),(Xi+1,Yi),
//  	// (Xi-1,Yi+1),(Xi-1,Yi-1),
//  	// (Xi+1,Yi+1), and (Xi+1,Yi-1),
// 	if(y_i - 1 >= 1){ //case 1
// 		neighbors[i-1].push_back((y_i-2)*n_x+x_i - 1);
// 	}
// 	if(y_i + 1 <= n_y){ //case 2
// 		neighbors[i-1].push_back(y_i*n_x+x_i - 1);
// 	}
// 	if(x_i-1 >= 1){ //case 3
// 		neighbors[i-1].push_back((y_i-1)*n_x+x_i - 2);
// 	}
// 	if(x_i+1 <= n_x){ // case 4
// 		neighbors[i-1].push_back((y_i-1)*n_x+x_i);
// 	}
// 	if((y_i+1<=n_y) && (x_i-1>=1)){ //c5
// 		neighbors[i-1].push_back(y_i*n_x+x_i-2);
// 	}
// 	if((y_i-1>=1)&&(x_i-1)>=1){ //c6
// 		neighbors[i-1].push_back((y_i-2)*n_x+x_i-2);
// 	}
// 	if((y_i + 1 <= n_y)&&(x_i + 1) <= n_x){//c7
// 		neighbors[i-1].push_back(y_i*n_x+x_i);
// 	}
// 	if((y_i - 1 >= 1) && (x_i+1)<=n_x){ //c8
// 		neighbors[i-1].push_back((y_i-2)*n_x+x_i);
// 	}

// 	for(int j = 1; j <= R; ++j){
// 		int x_j = j % n_x;
// 		int y_j = 0;
// 		if(x_j == 0){
// 			y_j = j / n_x;
// 		}else{
// 			y_j = 1 + (j-x_j)/n_x;
// 		}
// 		if(x_j == 0){
// 			x_j = n_x;
// 		}
// 		distance(i-1, j-1) = sqrt(pow(absc[x_i-1] - absc[x_j-1], 2) + pow(ord[y_i-1] - 
// 			ord[y_j-1], 2));
// 	}
// }