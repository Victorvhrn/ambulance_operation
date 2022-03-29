#include "../include/generator_regressor.h"

using namespace std;

GeneratorRegressor::GeneratorRegressor(GRBEnv& env, std::string calls_path, 
		std::string neighbors_path, std::string info_path): env(env){
	auto info_arq = ifstream(info_path, ios::in);
	info_arq >> T >> G >> R >> P >> nb_land_types >> nb_holidays_years; 
	slot_duration = 24 / T;
	daily_obs = std::vector<int>(G, 0);
	for(int g = 0; g < G; ++g){
		info_arq >> daily_obs[g];
	}
	info_arq.close();

	int max_obs = *max_element(daily_obs.begin(), daily_obs.end());

	xt::xarray<double>::shape_type shape_d = {T,G,R,P};
	xt::xarray<int>::shape_type shape_i = {T,G,R,P};
	nb_observations = xt::zeros<int>(shape_i);

	shape_i = {T,G,R,P,static_cast<ulong>(max_obs)};
	sample_calls = xt::zeros<int>(shape_i);

	shape_i = {T, nb_holidays_years, R, P};
	nb_observations_holidays = xt::zeros<int>(shape_i);
	shape_i = {T,nb_holidays_years,G,R,P};
	nb_calls_holidays = xt::zeros<int>(shape_i);
	shape_i = {T,G,R,P};
	nb_calls_no_holidays = xt::zeros<int>(shape_i);

	auto calls_arq = ifstream(calls_path, ios::in);
	std::string aux_str;
	is_holidays = std::vector<std::pair<bool, int>>(max_obs, make_pair(false,-1));
	do{
		std::getline(calls_arq, aux_str);
		if(aux_str == "END"){
			break;
		}
		std::istringstream ss(aux_str);
		int t, g, r, p, j, h, val;
		ss >> t >> g >> r >> p >> j >> val >> h;
		sample_calls(t,g,r,p,j) = val;
		if(h == -1){
			nb_calls_no_holidays(t,g,r,p) += sample_calls(t,g,r,p,j);
		}else{
			is_holidays[j] = make_pair(true, h);
			nb_observations_holidays(t,h,r,p) += 1;
			nb_calls_holidays(t,h,g,r,p) += sample_calls(t,g,r,p,j);
		}
		nb_observations(t,g,r,p) += 1;
	}while(true);

	calls_arq.close();


	neighbors = std::vector<vector<int>>(R, std::vector<int>());
	xt::xarray<double>::shape_type dist_shape = {R, R};
	regions = std::vector<Location>(R, null_location);
	distance = xt::xarray<double>(dist_shape);
	type = std::vector<int>(R,-1);
	xt::xarray<double>::shape_type reg_shape = {R, nb_land_types};
	regressors = xt::zeros<double>(reg_shape);
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

double GeneratorRegressor::cross_validation(xt::xarray<double>& x_beta, 
	xt::xarray<double>& x_delta, double sigma, double beta_bar, 
		double proportion){
	std::vector<double> alphas{0.001,0.01,0.05,0.1,0.5,1,2,5,10,50,100,1000};
	// std::vector<double> alphas{0.001,0.5,1};
	double max_likelihood = GRB_INFINITY;
	double best_alpha = -1;
	int nb_obs = static_cast<int>(sample_calls.shape(4));
	int nb_in_block = nb_obs*proportion;
	xt::xarray<double>::shape_type shape = {T,G,R,P};

	double beta_tilde = 1;
	double beta_hat = 2;

	for(double alpha: alphas){
		double likelihood = 0;
		for(int ind_cross = 0; ind_cross < floor(1/proportion); ++ind_cross){
			nb_observations = nb_in_block*xt::ones<int>(shape);
			nb_observations_holidays = xt::zeros<int>(nb_observations_holidays.shape());
			nb_calls_holidays = xt::zeros<int>(nb_calls_holidays.shape());
			nb_calls_no_holidays = xt::zeros<int>(nb_calls_no_holidays.shape());

			for(int ind = ind_cross*nb_in_block; ind < (ind_cross+1)*nb_in_block; ++ind){
				for(int t = 0; t < T; ++t){
					for(int g = 0; g < G; ++g){
						for(int r = 0; r < R; ++r){
							for(int p = 0; p < P; ++p){
								if(is_holidays[ind].first){
									int k = is_holidays[ind].second;
									nb_observations_holidays(t,k,r,p) += 1;
									nb_calls_holidays(t,k,g,r,p) += sample_calls(t,g,r,p,ind);
								}else{
									nb_calls_no_holidays(t,g,r,p) += sample_calls(t,g,r,p,ind);
								}
							}
						}
					}
				}
			}

			auto f_val = projected_gradient_armijo_feasible(x_beta, x_delta, alpha, 
				sigma, beta_tilde, beta_hat);
			xt::xarray<int> nb_observations_remaining = (nb_obs-nb_in_block)*xt::ones<int>(shape);
			xt::xarray<int> nb_observations_holidays_remaining = xt::zeros<int>(
				nb_observations_holidays.shape());
			xt::xarray<int> nb_calls_holidays_remaining = xt::zeros<int>(nb_calls_holidays.shape());
			xt::xarray<int> nb_calls_no_holidays_remaining = xt::zeros<int>(
				nb_calls_no_holidays.shape());

			for(int ind = 0; ind < ind_cross*nb_in_block; ++ind){
				for(int t = 0; t < T; ++t){
					for(int g = 0; g < G; ++g){
						for(int r = 0; r < R; ++r){
							for(int p = 0; p < P; ++p){
								if(is_holidays[ind].first){
									int k = is_holidays[ind].second;
									nb_observations_holidays_remaining(t,k,r,p) += 1;
									nb_calls_holidays_remaining(t,k,g,r,p) += sample_calls(t,g,r,p,ind);
								}else{
									nb_calls_no_holidays_remaining(t,g,r,p) += sample_calls(t,g,r,p,ind);
								}
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
								if(is_holidays[ind].first){
									int k = is_holidays[ind].second;
									nb_observations_holidays_remaining(t,k,r,p) += 1;
									nb_calls_holidays_remaining(t,k,g,r,p) += sample_calls(t,g,r,p,ind);
								}else{
									nb_calls_no_holidays_remaining(t,g,r,p) += sample_calls(t,g,r,p,ind);
								}
							}
						}
					}
				}
			}

			double f = 0;
			xt::xarray<double> rates = xt::zeros<double>(shape);
			for(int t = 0; t < T; ++t){
				for(int g = 0; g < G; ++g){
					for(int r = 0; r < R; ++r){
						for(int p = 0; p < P; ++p){
							for(int j = 0; j < nb_land_types; ++j){
								rates(t,g,r,p) += x_beta(t,g,p,j)*regressors(r,j);
							}
							f += nb_observations_remaining(t,g,r,p)*rates(t,g,r,p) - 
								nb_calls_no_holidays_remaining(t,g,r,p)*log(rates(t,g,r,p));
						}
					}
				}
			}
			for(int t = 0; t < T; ++t){
				for(int k = 0; k < nb_holidays_years; ++k){
					for(int r = 0; r < R; ++r){
						for(int p = 0; p < P; ++p){
							f += nb_observations_holidays_remaining(t,k,r,p)*x_delta(t,k,r,p);
							for(int g = 0; g < G; ++g){
								f -= nb_calls_holidays_remaining(t,k,g,r,p)*log(rates(t,g,r,p)*x_delta(t,k,r,p));
							}
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
	
	nb_observations = nb_obs*xt::ones<int>(shape);
	nb_observations_holidays = xt::zeros<int>(nb_observations_holidays.shape());
	nb_calls_holidays = xt::zeros<int>(nb_calls_holidays.shape());
	nb_calls_no_holidays = xt::zeros<int>(nb_calls_no_holidays.shape());
	for(int ind = 0; ind < nb_obs; ++ind){
		for(int t = 0; t < T; ++t){
			for(int g = 0; g < G; ++g){
				for(int r = 0; r < R; ++r){
					for(int p = 0; p < P; ++p){
						if(is_holidays[ind].first){
							int k = is_holidays[ind].second;
							nb_observations_holidays(t,k,r,p) += 1;
							nb_calls_holidays(t,k,g,r,p) += sample_calls(t,g,r,p,ind);
						}else{
							nb_calls_no_holidays(t,g,r,p) += sample_calls(t,g,r,p,ind);
						}
					}
				}
			}
		}
	}
	auto f_val_best = projected_gradient_armijo_feasible(x_beta, x_delta, best_alpha, 
		sigma, beta_tilde, beta_hat);

	return best_alpha;
}


void GeneratorRegressor::test(){
	double epsilon = 0.1;
	xt::xarray<double> x_delta = epsilon*xt::ones<double>(delta_teorico.shape());
	xt::xarray<double> x_beta = epsilon*xt::ones<double>(beta_teorico.shape());
	sigma = 0.5;
	double beta_bar = 1;
	max_iter = 10;
	double alpha = 0;

	double beta_tilde = 1;
	double beta_hat = 2;

	auto f_val1 = projected_gradient_armijo_boundary(x_beta, x_delta, alpha, sigma, 
		beta_bar);

	write_params(x_beta, x_delta, alpha);


	double proportion = 0.2;
	x_beta = epsilon*xt::ones<double>(x_beta.shape());
	x_delta = epsilon*xt::ones<double>(x_delta.shape());
	double best_alpha = cross_validation(x_beta, x_delta,sigma, beta_bar, proportion);
	write_params(x_beta, x_delta, best_alpha);


	alpha = pow(10,7);
	x_beta = epsilon*xt::ones<double>(x_beta.shape());
	x_delta = epsilon*xt::ones<double>(x_delta.shape());

	auto f_val2 = projected_gradient_armijo_feasible(x_beta, x_delta, alpha, sigma, 
		beta_tilde, beta_hat);

	write_params(x_beta, x_delta, alpha);

	// std::cout << "Avg diff Feasible = " << average_difference(x_beta, x_delta) << "\n";
	
	// double proportion = 0.2;

	// x_beta = epsilon*xt::ones<double>(x_beta.shape());
	// x_delta = epsilon*xt::ones<double>(x_delta.shape());
	// double best_alpha = cross_validation(x_beta, x_delta,sigma, beta_bar, proportion);
	// std::cout << "Best alpha: " << best_alpha << "\n";
	// std::cout << "Avg diff Cross = " << average_difference(x_beta, x_delta) << "\n";
}


std::vector<double> GeneratorRegressor::projected_gradient_armijo_boundary(
	xt::xarray<double>& x_beta, xt::xarray<double>& x_delta, double alpha, 
	double sigma, double beta_bar){

	xt::xarray<double> z_beta = xt::zeros<double>(x_beta.shape());
	xt::xarray<double> z_delta = xt::zeros<double>(x_delta.shape());
	int k = 0;
	double eps = 0.1;
	std::vector<double> f_val;
	f_val.reserve(max_iter);
	while(k < max_iter){
		double fold = oracle_objective_model2(x_beta, x_delta, alpha);
		xt::xarray<double> gradient_beta, gradient_delta;
		tie(gradient_beta, gradient_delta) = oracle_gradient_model2(x_beta, x_delta, 
			alpha);
		bool stop = false;
		int j = 0;
		double f = 0;
		do{
			xt::xarray<double> x_beta_aux = x_beta - (beta_bar/pow(2,j))*gradient_beta;
			xt::xarray<double> x_delta_aux = x_delta - (beta_bar/pow(2,j))*gradient_delta;
			xt::xarray<double> z_beta, z_delta;
			tie(z_beta, z_delta) = projection_regressors(x_beta_aux, x_delta_aux, eps);
			f = oracle_objective_model2(z_beta, z_delta, alpha);
			double rhs = fold - sigma*(xt::sum(gradient_beta*(x_beta-z_beta))() + 
				xt::sum(gradient_delta*(x_delta-z_delta))());
			fmt::print("k = {}, j = {}, f = {}, rhs = {}\n", k, j, f, rhs);
			if(f <= rhs){
				stop = true;
			}else{
				++j;
			}
		}while(!stop);
		f_val.push_back(f);
		x_beta = z_beta;
		x_delta = z_delta;
		++k;
	}

	return f_val;
}

std::vector<double> GeneratorRegressor::projected_gradient_armijo_feasible(
	xt::xarray<double>& x_beta, xt::xarray<double>& x_delta, double alpha, 
	double sigma, double beta_tilde, double beta_hat){

	xt::xarray<double> z_beta = xt::zeros<double>(x_beta.shape());
	xt::xarray<double> z_delta = xt::zeros<double>(x_delta.shape());
	int k = 0;
	double eps = 0.1;
	double b_param = 2;
	std::vector<double> f_val;
	f_val.reserve(max_iter);
	int j = 0;
	while(k < max_iter){
		double beta_k = 0;
		if(k == 0){
			beta_k = b_param;
		}else{
			beta_k = b_param/pow(2,j);
		}
		double fold = oracle_objective_model2(x_beta, x_delta, alpha);
		xt::xarray<double> gradient_beta, gradient_delta;
		tie(gradient_beta, gradient_delta) = oracle_gradient_model2(x_beta, 
			x_delta, alpha);
		xt::xarray<double> x_beta_aux = x_beta-beta_k*gradient_beta;
		xt::xarray<double> x_delta_aux = x_delta - beta_k*gradient_delta;
		xt::xarray<double> z_beta, z_delta;
		tie(z_beta, z_delta) = projection_regressors(x_beta_aux, x_delta_aux, eps);
		bool stop = false;
		j = 0;
		double f = 0;
		double rhs = xt::sum(gradient_beta*(x_beta-z_beta))() + xt::sum(gradient_delta*
			(x_delta - z_delta))();
		do{
			xt::xarray<double> z_aux_beta = x_beta + (1/pow(2,j))*(z_beta-x_beta);
			xt::xarray<double> z_aux_delta = x_delta + (1/pow(2,j))*(z_delta-x_delta);
			f = oracle_objective_model2(z_aux_beta, z_aux_delta, alpha);

			if(f <= fold-(sigma/pow(2,j))*rhs){
				stop = true;
			}else{
				++j;
			}
		}while(!stop);
		fmt::print("k = {}, f = {}, j = {}\n", k, f, j);
		f_val.push_back(f);
		x_beta += (1/pow(2,j))*(z_beta-x_beta);
		x_delta += (1/pow(2,j))*(z_delta-x_delta);
		++k;
	}
	return f_val;
}

std::pair<xt::xarray<double>,xt::xarray<double>> 
	GeneratorRegressor::oracle_gradient_model2(xt::xarray<double>& x_beta, 
		xt::xarray<double>& x_delta, double alpha){
	
	xt::xarray<double> gradient_beta = xt::zeros<double>(x_beta.shape());
	xt::xarray<double> gradient_delta = xt::zeros<double>(x_delta.shape());
	xt::xarray<double>::shape_type shape = {T,G,R,P};
	xt::xarray<double> rates = xt::zeros<double>(shape);

	for(int t = 0; t < T; ++t){
		for(int g = 0; g < G; ++g){
			for(int r = 0; r < R; ++r){
				for(int p = 0; p < P; ++p){
					for(int j = 0; j < nb_land_types; ++j){
						rates(t,g,r,p)+= x_beta(t,g,p,j)*regressors(r,j);
					}
					for(int j = 0; j < nb_land_types; ++j){
						gradient_beta(t,g,p,j) += nb_observations(t,g,r,p)*regressors(r,j) 
							- nb_calls_no_holidays(t,g,r,p)*regressors(r,j)/rates(t,g,r,p);
						for(int k = 0; k < nb_holidays_years; ++k){
							gradient_beta(t,g,p,j) -= nb_calls_holidays(t,k,g,r,p)*regressors(r,j) / 
								(rates(t,g,r,p)+x_delta(t,k,r,p));
						}
					}
				}
			}
		}
	}

	for(int t = 0; t < T; ++t){
		for(int k = 0; k < nb_holidays_years; ++k){
			for(int r = 0; r < R; ++r){
				for(int p = 0; p < P; ++p){
					gradient_delta(t,k,r,p) += nb_observations_holidays(t,k,r,p);
					for(int g = 0; g < G; ++g){
						gradient_delta(t,k,r,p) -= nb_calls_holidays(t,k,g,r,p) /
							(rates(t,g,r,p)+x_delta(t,k,r,p));
					}
					for(auto s: neighbors[r]){
						gradient_delta(t,k,r,p) += 2*alpha*(x_delta(t,k,r,p) - 
							x_delta(t,k,s,p))/distance(r,s);
					}
				}
			}
		}
	}
	return make_pair(gradient_beta, gradient_delta);
}


double GeneratorRegressor::oracle_objective_model2(xt::xarray<double>& x_beta, 
		xt::xarray<double>& x_delta, double alpha){
	double f = 0;
	xt::xarray<double>::shape_type shape = {T,G,R,P};
	xt::xarray<double> rates = xt::zeros<double>(shape);
	for(int t = 0; t < T; ++t){
		for(int g = 0; g < G; ++g){
			for(int r = 0; r < R; ++r){
				for(int p = 0; p < P; ++p){
					for(int j = 0; j < nb_land_types; ++j){
						rates(t,g,r,p)+= x_beta(t,g,p,j)*regressors(r,j);
					}
					f += nb_observations(t,g,r,p)*rates(t,g,r,p) - 
						nb_calls_no_holidays(t,g,r,p)*log(rates(t,g,r,p));
				}
			}
		}
	}
	for(int t = 0; t < T; ++t){
		for(int k = 0; k < nb_holidays_years; ++k){
			for(int r = 0; r < R; ++r){
				for(int p = 0; p < P; ++p){
					f += nb_observations_holidays(t,k,r,p)*x_delta(t,k,r,p);
					for(int g = 0; g < G; ++g){
						f -= nb_calls_holidays(t,k,g,r,p)*log(rates(t,g,r,p) + 
							x_delta(t,k,r,p));
					}
					for(auto s: neighbors[r]){
						f += (0.5*alpha)*pow(x_delta(t,k,r,p) - x_delta(t,k,s,p), 2) /
							distance(r,s);
					}
				}
			}
		}
	}
	return f;
}

std::pair<xt::xarray<double>,xt::xarray<double>> 
	GeneratorRegressor::projection_regressors(xt::xarray<double>& x_beta, 
		xt::xarray<double>& x_delta, double epsilon){

	xt::xarray<GRBVar>::shape_type delta_shape = {T,nb_holidays_years,R,P};
	xt::xarray<GRBVar> y_delta(delta_shape);
	xt::xarray<GRBVar>::shape_type beta_shape = {T,G,P,nb_land_types};
	xt::xarray<GRBVar> y_beta(beta_shape);

	GRBModel model(env);
	stringstream name;
	for(int t = 0; t < T; ++t){
		for(int k = 0; k < nb_holidays_years; ++k){
			for(int r = 0; r < R; ++r){
				for(int p = 0; p < P; ++p){
					name << "yd_" << t << "_" << k << "_" << r << "_" << p;
					y_delta(t,k,r,p) = model.addVar(-GRB_INFINITY, GRB_INFINITY,
						0,GRB_CONTINUOUS, name.str());
					name.str("");
				}
			}
		}
	}

	for(int t = 0; t < T; ++t){
		for(int g = 0; g < G; ++g){
			for(int p = 0; p < P; ++p){
				for(int j = 0; j < nb_land_types; ++j){
					name << "yb_" << t << "_" << g << "_" << p << "_" << j;
					y_beta(t,g,p,j) = model.addVar(-GRB_INFINITY, GRB_INFINITY,
						0,GRB_CONTINUOUS, name.str());
					name.str("");
				}
			}
		}
	}

	GRBQuadExpr obj;
	for(int t = 0; t < T; ++t){
		for(int k = 0; k < nb_holidays_years; ++k){
			for(int r = 0; r < R; ++r){
				for(int p = 0; p < P; ++p){
					obj += 0.5*y_delta(t,k,r,p)*y_delta(t,k,r,p) -
						x_delta(t,k,r,p)*y_delta(t,k,r,p);
				}
			}
		}
	}

	for(int t = 0; t < T; ++t){
		for(int g = 0; g < G; ++g){
			for(int p = 0; p < P; ++p){
				for(int j = 0; j < nb_land_types; ++j){
					obj += 0.5*y_beta(t,g,p,j)*y_beta(t,g,p,j) -
						x_beta(t,g,p,j)*y_beta(t,g,p,j);
				}
			}
		}
	}
	try{
		model.setObjective(obj, GRB_MINIMIZE);
	}catch(GRBException& ex){
		cout << ex.getMessage() << "\n";
	}

	GRBLinExpr con1;
	for(int t = 0; t < T; ++t){
		for(int g = 0; g < G; ++g){
			for(int p = 0; p < P; ++p){
				for(int r = 0; r < R; ++r){

					for(int j = 0; j < nb_land_types; ++j){
						con1 += y_beta(t,g,p,j)*regressors(r,j);
					}
					name << "con1_" << t << "_" << g << "_" << p;
					model.addConstr(con1, GRB_GREATER_EQUAL, epsilon, name.str());
					name.str("");
					con1 = 0;
				}
			}
		}
	}
	GRBLinExpr con2;
	for(int t = 0; t < T; ++t){
		for(int g = 0; g < G; ++g){
			for(int p = 0; p < P; ++p){
				for(int r = 0; r < R; ++r){
					for(int k = 0; k < nb_holidays_years; ++k){
						for(int j = 0; j < nb_land_types; ++j){
							con2 += y_beta(t,g,p,j)*regressors(r,j) + y_delta(t,k,r,p);
						}
						name << "con2_" << t << "_" << g << "_" << p << "_" << r;
						name << "_" << k;
						model.addConstr(con2, GRB_GREATER_EQUAL, epsilon, name.str());
						name.str("");
						con2 = 0;
					}
				}
			}
		}
	}

	model.set(GRB_IntParam_OutputFlag,0);
	model.optimize();

	xt::xarray<double> beta_val(y_beta.shape());
	xt::xarray<double> delta_val(y_delta.shape());
	for(int t = 0; t < T; ++t){
		for(int k = 0; k < nb_holidays_years; ++k){
			for(int r = 0; r < R; ++r){
				for(int p = 0; p < P; ++p){
					delta_val(t,k,r,p) = y_delta(t,k,r,p).get(GRB_DoubleAttr_X);
				}
			}
		}
	}

	for(int t = 0; t < T; ++t){
		for(int g = 0; g < G; ++g){
			for(int p = 0; p < P; ++p){
				for(int j = 0; j < nb_land_types; ++j){
					beta_val(t,g,p,j) = y_beta(t,g,p,j).get(GRB_DoubleAttr_X);
				}
			}
		}
	}

	return make_pair(beta_val, delta_val);

}


double GeneratorRegressor::average_difference(xt::xarray<double>& x_beta,
	xt::xarray<double>& x_delta){
	int count = 0;
	double sum_diff = 0;
	for(int t = 0; t < T; ++t){
		for(int g = 0; g < G; ++g){
			for(int p = 0; p < P; ++p){
				for(int j = 0; j < nb_land_types; ++j){
					sum_diff += abs(beta_teorico(t,g,p,j) - x_beta(t,g,p,j));
					++count;
				}	
			}
		}
	}

	for(int t = 0; t < T; ++t){
		for(int k = 0; k < nb_holidays_years; ++k){
			for(int r = 0; r < R; ++r){
				for(int p = 0; p < P; ++p){
					sum_diff += abs(delta_teorico(t,k,r,p) - x_delta(t,k,r,p));
					++count;
				}
			}
		}
	}


	return sum_diff / count;
}




void GeneratorRegressor::comp_wise_max(xt::xarray<double>& z ,xt::xarray<double>& a, double eps){
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


bool GeneratorRegressor::is_neighbor(int r, int s){
	return r != s;
}


void GeneratorRegressor::write_params(xt::xarray<double>& x_beta, 
	xt::xarray<double>& x_delta, double alpha){
	ofstream out_file(fmt::format("{}/xRegressorT{}G{}I{}P{}K{}J{}_alpha{}.txt",
		g_params.generator_folder,T,G,R,P, nb_holidays_years, nb_land_types,alpha), ios::out);

	for(int t = 0; t < T; ++t){
		for(int k = 0; k < nb_holidays_years; ++k){
			for(int r = 0; r < R; ++r){
				for(int p = 0; p < P; ++p){
					out_file << t << " " << k << " " << r << " " << p;
					out_file << " " << x_delta(t,k,r,p) << "\n";
				}
			}
		}
	}

	for(int t = 0; t < T; ++t){
		for(int g = 0; g < G; ++g){
			for(int p = 0; p < P; ++p){
				for(int j = 0; j < nb_land_types; ++j){
					out_file << t << " " << g << " " << p << " " << j;
					out_file << " " << x_beta(t,g,p,j) << "\n";
				}
			}
		}
	}

	out_file << "END";
	out_file.close();

	ofstream plot(fmt::format("{}/params_plot_reg_{}.txt", g_params.generator_folder, alpha), 
		ios::out);
	for(int k = 0; k < nb_holidays_years; ++k){
		plot << k << " ";
		for(int t = 0; t < T; ++t){
			double sum = 0;
			for(int r = 0; r < R; ++r){
				for(int p = 0; p < P; ++p){
					sum += x_delta(t,k,r,p);
				}
			}
			plot << sum << " ";
		}
		plot << "\n";
	}
	for(int g = 0; g < G; ++g){
		plot  << g << " ";
		for(int t = 0; t < T; ++t){
			double sum = 0;
			for(int p = 0; p < P; ++p){
				for(int j = 0; j < nb_land_types; ++j){
					sum += x_beta(t,g,p,j);
				}
			}
			plot << sum << " ";
		}
		plot << "\n";
	}
	plot.close();
}


	// for(int j = 0; j < nb_weeks*nb_years*G; ++j){
	// 	is_holidays.push_back(make_pair(false,0));
	// }

	// for(int year = 0; year < nb_years; ++year){
	// 	int day = (year+1) % 7;
	// 	if(day == 0){
	// 		day = 7;
	// 	}
	// 	is_holidays[((year)*nb_weeks)*G + day] = make_pair(true,0);
	// }

	// for(int year = 0; year < nb_years; ++year){
	// 	int day = (year+4) % 7;
	// 	if(day == 0){
	// 		day = 7;
	// 	}
	// 	is_holidays[(11 + (year)*nb_weeks)*G + day] = make_pair(true,1);
	// }

	// for(int year = 0; year < nb_years; ++year){
	// 	int day = (year+6) % 7;
	// 	if(day == 0){
	// 		day = 7;
	// 	}
	// 	is_holidays[(23 + (year)*nb_weeks)*G + day] = make_pair(true,2);
	// }

	// for(int year = 0; year < nb_years; ++year){
	// 	int day = (year+6) % 7;
	// 	if(day == 0){
	// 		day = 7;
	// 	}
	// 	is_holidays[(30 + (year)*nb_weeks)*G + day] = make_pair(true,3);
	// }


	// for(int year = 0; year < nb_years; ++year){
	// 	int day = (year+7) % 7;
	// 	if(day == 0){
	// 		day = 7;
	// 	}
	// 	is_holidays[(40 + (year)*nb_weeks)*G + day] = make_pair(true,4);
	// }

	// for(int year = 0; year < nb_years; ++year){
	// 	int day = (year+3) % 7;
	// 	if(day == 0){
	// 		day = 7;
	// 	}
	// 	is_holidays[(50 + (year)*nb_weeks)*G + day] = make_pair(true,5);
	// }