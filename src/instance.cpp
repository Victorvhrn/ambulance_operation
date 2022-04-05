#include "../include/instance.h"
#include "../include/call.h"
#include "../include/travel.h"
#include "../include/ambulance.h"
#include "../include/osrm_helper.h"


// Constructor for tests generating scenarios for all 7 days;
Instance::Instance(): path(""),travel(false){
	nb_times = 48;
	slot_duration = 0.5*3600;
	int nb_days = 7;
	nb_scenarios = nb_days*g_params.n_scenarios;
	nb_hospitals = g_params.n_hospitals;
	nb_bases = g_params.n_bases;
	nb_cleaning_bases = 4;
	nb_types_ambulance = 3;
	nb_ambulances = g_params.n_ambulances;
	nb_priorities = 3;
	nb_regions = 76;
	int nb_land_types = 5;

	xt::xarray<double>::shape_type shape = {(ulong)nb_times, (ulong) nb_days, (ulong)nb_regions, 
		(ulong)nb_priorities};
	xt::xarray<double> lambda = xt::zeros<double>(shape);

	auto lambda_arq = ifstream(fmt::format("{}/empirical_estimation.txt", 
		g_params.generator_folder), ios::in);
	int t,g,r,p;
	double val;
	for(int i = 0; i < nb_times*nb_days*nb_regions*nb_priorities; ++i){
		lambda_arq >> t >> g >> r >> p >> val;
		lambda(t,g,r,p) = val;
	}
	lambda_arq.close();

	fmt::print("Loaded Lambda\n");
	penalties = {4,2,1};

	A.push_back({0});
	A.push_back({0,1});
	A.push_back({0,1,2});

	auto neighbors_arq = ifstream(fmt::format("{}/neighbors.dat", g_params.generator_folder), 
		ios::in);
	centers = vector<Location>(nb_regions, null_location);

	string aux_str;
	while(true){
		int ind, terrain_type, s; 
		double lat, longi, dist;
		std::getline(neighbors_arq, aux_str);
		if(aux_str == "END"){
			break;
		}
		std::istringstream ss(aux_str);
		ss >> ind >> lat >> longi >> terrain_type;
		centers[ind] = make_pair(lat,longi);
		int regressor;
		for(int j = 0; j < nb_land_types; ++j){
			ss >> regressor;
		}
		while(ss >> s >> dist){
		}
	}

	// for(int r = 0; r < nb_regions; ++r){
	// 	fmt::print("{} {} {}\n",r,centers[r].first, centers[r].second);
	// }

	// std::cin.get();
	neighbors_arq.close();

	fmt::print("Loaded Regions\n");

	auto hosp_arq = ifstream("hospitals.txt", ios::in);
	int total_h;
	hosp_arq >> total_h;
	if(nb_hospitals > total_h){
		fmt::print("nb_hospitals {} > total hospitals {}\n", nb_hospitals, total_h);
		exit(1);
	}
	for(int h = 0; h < nb_hospitals; ++h){
		double lat, longi;
		hosp_arq >> lat >> longi;
		hospitals.push_back(make_pair(lat,longi));
	}
	hosp_arq.close();
	std::vector<int> nearest_hospital(nb_regions, -1);
	auto regions_hospitals_table = travel.table_in_out(centers, hospitals);
	for(int r = 0; r < nb_regions; ++r){
		double min_time = GRB_INFINITY;
		for(int h = 0; h < nb_hospitals; ++h){
			// fmt::print("{} {} to {} {} = {}\n",
			// 	centers[r].first, centers[r].second,
			// 	hospitals[h].first, hospitals[h].second, 
			// 	regions_hospitals_table[r][h]);
			if(regions_hospitals_table[r][h] < min_time){
				min_time = regions_hospitals_table[r][h];
				nearest_hospital[r] = h;
			}
		}
	}

	fmt::print("Loaded Region/Hospital Distances\n");

	auto base_arq = ifstream("bases.txt", ios::in);
	int total_b;
	base_arq >> total_b;
	if(nb_bases > total_b){
		fmt::print("nb_bases {} > total bases {}\n", nb_bases, total_b);
		exit(1);
	}
	for(int b = 0; b < nb_bases; ++b){
		double lat, longi;
		base_arq >> lat >> longi;
		bases.push_back(make_pair(lat,longi));
		cap_bases.push_back((nb_ambulances/nb_bases) + 4);
	}
	base_arq.close();
	std::default_random_engine gen;
	auto samples_arq = ifstream(fmt::format("{}/samples.dat", g_params.generator_folder), 
		ios::in);
	std::vector<std::vector<Location>> samples(nb_regions, std::vector<Location>(100, 
		null_location));
	for(int i = 0; i < nb_regions*100; ++i){
		int r, s;
		double lat, longi;
		samples_arq >> r >> s >> longi >> lat;
		samples[r][s].first = lat;
		samples[r][s].second = longi;
	}

	samples_arq.close();

	for(int g = 0; g < nb_days; ++g){
		for(int s = 0; s < g_params.n_scenarios; ++s){
			std::vector<Call> scenario;
			int id = 0;
			for(int t = 0; t < nb_times; ++t){
				for(int r = 0; r < nb_regions; ++r){
					for(int p = 0; p < nb_priorities; ++p){
						poisson_distribution<int> dist(lambda(t,g,r,p));
						int nb_calls = dist(gen);
						for(int k = 0; k < nb_calls; ++k){
							std::uniform_int_distribution<int> int_dist(t*slot_duration, 
								(t+1)*slot_duration);
							std::uniform_int_distribution<int> sample_rand(0,99);
							double time = int_dist(gen);
							Location location = samples[r][sample_rand(gen)];
							double time_on_scene = 1200;
							double time_at_hospital = 1200;
							double cleaning_time = 2400;
							int min_h = nearest_hospital[r];
							int cb = 0;
							scenario.push_back(Call(id++, time, location, min_h, cb, p,
								true, false, time_on_scene, time_at_hospital, cleaning_time));
							scenario.back().region = r;
							scenario.back().disc_time = t;
						}
					}
				}
			}
			std::sort(scenario.begin(),scenario.end());
			for(int i = 0; i < scenario.size(); ++i){
				scenario[i].id = i;
			}
			calls.push_back(scenario);
			nb_calls.push_back(scenario.size());
		}
	}
	
	fmt::print("Loaded Calls, {} scenarios {}\n", calls.size(),nb_days);

	// std::uniform_int_distribution<int> amb_dist(0, nb_bases-1);
	for(int a = 0; a < nb_ambulances; ++a){
		// int b = amb_dist(gen);
		Location base_location = bases[a % nb_bases];
		Location hospital_location = null_location;
		ambulances.push_back(Ambulance(a, base_location, hospital_location, -1, 0, 
			(a % nb_types_ambulance), 0.005));
	}


	nearest_base_to_region = std::vector<int>(nb_regions, -1);

	for(int r = 0; r < nb_regions; ++r){
		double min_d = GRB_INFINITY;
		for(int b = 0; b < nb_bases; ++b){
			double d = travel.norm(centers[r], bases[b]);
			if(d < min_d){
				min_d = d;
				nearest_base_to_region[r] = b;
			}
		}	
	}
}

// Instance::Instance(): path(g_params.instance), travel(false){
// 	nb_types_ambulance = 3;
// 	nb_priorities = 3;
// 	slot_duration = 3600;
// 	path = g_params.instance;
// 	nb_times = g_params.n_time_slots;
// 	nb_regions = g_params.n_regions;
// 	nb_hospitals = g_params.n_hospitals;
// 	nb_bases = g_params.n_bases;
// 	nb_cleaning_bases = g_params.n_cleaning_bases;
// 	nb_ambulances = g_params.n_ambulances;

// 	bool is_euclidian = path == simulated_rect_1 || path == simulated_rect_2;

// 	penalties = {4,2,1};

// 	A.push_back({0});
// 	A.push_back({0,1});
// 	A.push_back({0,1,2});


// 	if(is_euclidian){
// 		load_euclidian();
// 	}else{
// 		load_real_data();
// 	}


// 	for(int a = 0; a < nb_ambulances; ++a){
// 		// int b = amb_dist(gen);
// 		Location base_location = bases[a % nb_bases];
// 		Location hospital_location = null_location;
// 		const double speed = 0.078; //approximate speed needed to cross square in 30m
// 		ambulances.push_back(Ambulance(a, base_location, hospital_location, -1, 0, 
// 			(a % nb_types_ambulance), speed));
// 	}
// 	// fmt::print("Ambulances:\n");
// 	// for(auto& amb: ambulances){
// 	// 	std::cout << amb << " " << amb.base_location.first << " " << amb.base_location.second;
// 	// 	std::cout << "\n";
// 	// }
// }

void Instance::load_euclidian(){
	travel.euclidian = true;

	nb_hospitals = 4;
	nb_bases = 4;
	nb_cleaning_bases = 1;

	hospitals.push_back(make_pair(50,0));
	hospitals.push_back(make_pair(100,50));
	hospitals.push_back(make_pair(50,100));
	hospitals.push_back(make_pair(0,50));

	bases.push_back(make_pair(0,0));
	bases.push_back(make_pair(100,0));
	bases.push_back(make_pair(100,100));
	bases.push_back(make_pair(0,100));


	cleaning_bases.push_back(make_pair(50,50));

	default_random_engine gen;
	auto calls_arq = ifstream(path, ios::in);
	int scenario_size = 0;
	double time, time_on_scene, lat, longi, time_at_hospital, time_at_cleaning;
	int region, priority, day, cleaning_needed, hospital_needed, index_hospital, 
		index_cleaning;
	int k = 0;
	while(k++ < g_params.n_scenarios){
		calls_arq >> scenario_size;
		nb_calls.push_back(scenario_size);
		std::vector<Call> scenario;
		scenario.reserve(scenario_size);
		int id = 0;
		for(int i = 0; i < scenario_size; ++i){
			calls_arq >> time >> lat >> longi >> priority >> cleaning_needed;
			calls_arq >> time_at_cleaning >> index_hospital >> time_on_scene;
			calls_arq >> hospital_needed >> time_at_hospital;
			region = -1;
			priority -= 1;
			index_hospital -= 1;
			uniform_int_distribution<int> int_dist((time-1)*slot_duration, 
				time*slot_duration);

			Location call_location = make_pair(lat,longi);
			double min_d = GRB_INFINITY;
			for(int h = 0; h < nb_hospitals; ++h){
				double d = travel.norm(call_location, hospitals[h]);
				if(d < min_d){
					min_d = d;
					index_hospital = h;
				}
			}


			int index_cleaning = 0;

			if(priority == 0){
				std::uniform_int_distribution<int> cb_rand(0,4);
				if(cb_rand(gen) == 0){
					cleaning_needed = 1;
				}
			}

			scenario.push_back(Call(id++, time*slot_duration, call_location, index_hospital, 
				index_cleaning, priority, hospital_needed, cleaning_needed,
				0.333*slot_duration, 0.333*slot_duration, 
				0.666*slot_duration));
		}
		// fmt::print("scenario size {}\n", scenario.size());
		// scenario = std::vector<Call>(scenario.begin(), scenario.begin()+10);
		calls.push_back(scenario);
	}
	calls_arq.close();
	nb_scenarios = min(static_cast<int>(calls.size()), g_params.n_scenarios);
	calls = std::vector<std::vector<Call>>(calls.begin(), calls.begin() + nb_scenarios);
}


void Instance::load_real_data(){
	auto hosp_arq = ifstream("hospitals.txt", ios::in);
	int total_h;
	hosp_arq >> total_h;
	if(nb_hospitals > total_h){
		fmt::print("nb_hospitals {} > total hospitals {}\n", nb_hospitals, total_h);
		exit(1);
	}
	for(int h = 0; h < nb_hospitals; ++h){
		double lat, longi;
		hosp_arq >> lat >> longi;
		hospitals.push_back(make_pair(lat,longi));
	}
	hosp_arq.close();
	// fmt::print("Hospitals:\n");
	// for(auto h: hospitals){
	// 	fmt::print("{} {}\n",h.first, h.second);
	// }

	auto base_arq = ifstream("bases.txt", ios::in);
	int total_b;
	base_arq >> total_b;
	if(nb_bases > total_b){
		fmt::print("nb_bases {} > total bases {}\n", nb_bases, total_b);
		exit(1);
	}
	for(int b = 0; b < nb_bases; ++b){
		double lat, longi;
		base_arq >> lat >> longi;
		bases.push_back(make_pair(lat,longi));
		cap_bases.push_back((nb_ambulances/nb_bases) + 4);
	}
	base_arq.close();
	// fmt::print("Bases:\n");
	// for(auto h: bases){
	// 	fmt::print("{} {}\n",h.first, h.second);
	// }


	auto cleaning_arq = ifstream("cleaning.txt", ios::in);
	int total_cb;
	cleaning_arq >> total_cb;
	if(nb_cleaning_bases > total_cb){
		fmt::print("cleaning_bases {} > total cleaning bases {}\n", nb_cleaning_bases,
			total_cb);
		exit(1);
	}
	for(int c = 0; c < nb_cleaning_bases; ++c){
		double lat, longi;
		cleaning_arq >> lat >> longi;
		cleaning_bases.push_back(std::make_pair(lat,longi));
	}
	cleaning_arq.close();
	// fmt::print("Cleaning Bases:\n");
	// for(auto h: cleaning_bases){
	// 	fmt::print("{} {}\n",h.first, h.second);
	// }

	auto samples_arq = ifstream(fmt::format("{}/samples.dat", g_params.generator_folder),
		ios::in);
	std::vector<std::vector<Location>> samples(nb_regions, std::vector<Location>(100, 
		null_location));
	for(int i = 0; i < nb_regions*100; ++i){
		int r, s;
		double lat, longi;
		samples_arq >> r >> s >> longi >> lat;
		samples[r][s].first = lat;
		samples[r][s].second = longi;
	}
	samples_arq.close();



	std::default_random_engine gen;
	auto calls_arq = ifstream(path, ios::in);
	int scenario_size = 0;
	// Time - Region index - Priority - Day - Time on scene - Lat - Long - Time at hospital - 
	// TimeCleaningBase - Cleaning needed - Hospital needed - index hospital
	double time, time_on_scene, lat, longi, time_at_hospital, time_at_cleaning;
	int region, priority, day, cleaning_needed, hospital_needed, index_hospital, index_cleaning;
	while(calls_arq >> scenario_size){
		nb_calls.push_back(scenario_size);
		std::vector<Call> scenario;
		scenario.reserve(scenario_size);
		int id = 0;
		for(int i = 0; i < scenario_size; ++i){
			calls_arq >> time >> region >> priority >> day >> time_on_scene;
			calls_arq >> lat >> longi >> time_at_hospital >> time_at_cleaning;
			calls_arq >> cleaning_needed >> hospital_needed >> index_hospital;
			// fmt::print("line: {} {} {} {} {} {} {} {} {} {} {} {}\n",time,region,priority,
			// 	day,time_on_scene,lat,longi,time_at_hospital,time_at_cleaning, cleaning_needed,
			// 	hospital_needed, index_hospital);
			region -= 1;
			priority -= 1;
			std::uniform_int_distribution<int> int_dist((time-1)*slot_duration, time*slot_duration);
			std::uniform_int_distribution<int> sample_rand(0,99);
			double min_d = GRB_INFINITY;
			int sample_id = sample_rand(gen);
			auto call_location = samples[region][sample_id];
			index_cleaning = -1;
			Location previous_location = (hospital_needed) ? hospitals[index_hospital] : 
					call_location;
			for(int c = 0; c < nb_cleaning_bases; ++c){
				double travel_time = travel.travel_time(previous_location, cleaning_bases[c]);
				if(travel_time < min_d){
					min_d = travel_time;
					index_cleaning = c;
				}
			}
			min_d = GRB_INFINITY;
			index_hospital = -1;
			if(g_params.h_random_hospital){
				std::uniform_int_distribution<int> h_rand(0,nb_hospitals-1);
				index_hospital = h_rand(gen);
			}else{
				for(int h = 0; h < nb_hospitals; ++h){
					double travel_time = travel.travel_time(call_location,hospitals[h]);
					if(travel_time < min_d){
						min_d = travel_time;
						index_hospital = h;
					}
				}	
			}
			

			if(priority == 0){
				std::uniform_int_distribution<int> cb_rand(0,4);
				if(cb_rand(gen) == 0){
					cleaning_needed = 1;
				}
			}

			scenario.push_back(Call(id++, time*slot_duration, call_location, index_hospital, 
				index_cleaning, priority, hospital_needed, cleaning_needed,
				0.333*slot_duration, 0.333*slot_duration, 
				0.666*slot_duration));
			auto& call = scenario.back();
			call.region = region;
			// std::cout << call << " " << call.location.first << " " << call.location.second << "\n";
		}
		// fmt::print("=============================\n");
		calls.push_back(scenario);
	}
	calls_arq.close();
	nb_scenarios = min(static_cast<int>(calls.size()), g_params.n_scenarios);
	calls = std::vector<std::vector<Call>>(calls.begin(), calls.begin() + nb_scenarios);
}


Instance::Instance(std::string path): path(path), travel(true){
	std::ifstream arq(path, std::ios::in);
	std::cout << path << "\n";
	arq >> nb_scenarios;
	// arq >> nb_calls;
	arq >> nb_hospitals;
	arq >> nb_bases;
	arq >> nb_cleaning_bases;
	arq >> nb_types_ambulance;
	arq >> nb_ambulances;
	arq >> nb_priorities;
	arq >> x_min >> x_max >> y_min >> y_max;

	nb_calls = vector<int>(nb_scenarios, 0);
	std::cout << nb_scenarios << " ";
	std::cout << nb_hospitals << " " << nb_bases << " ";
	std::cout << nb_cleaning_bases << " " << nb_types_ambulance << " ";
	std::cout << nb_ambulances << " " << nb_priorities << "\n";
	std::cout << x_min << " " << x_max << " " << y_min << " " << y_max << "\n";

	A = vector<set<int>>(nb_priorities, set<int>());
	int aux = 0;
	std::string aux_str;
	// arq.ignore();
	// for(int c = 0; c < nb_priorities; ++c){
	// 	std::getline(arq, aux_str);
	// 	std::istringstream ss(aux_str);
	// 	while(ss >> aux){
	// 		A[c].insert(aux);
	// 	}
	// }

	// for(int c = 0; c < nb_priorities; ++c){
	// 	std::cout << c << ": ";
	// 	for(auto a: A[c]){
	// 		std::cout << a << " ";
	// 	}
	// 	std::cout << "\n";
	// }
	// std::cin.get();
	
	double time, lat, longi, time_on_scene, time_at_hospital, cleaning_time;
	int  priority, hospital, cleaning, hosp_needed, clean_needed;
	
	penalties.reserve(nb_priorities);
	for(int i = 0; i < nb_priorities; ++i){
		double aux;
		arq >> aux;
		penalties.push_back(aux);
		std::cout << penalties.back() << " ";
	}
	std::cout << "\n";

	penalties = {4,2,1};

	for(int k = 0; k < nb_scenarios; ++k){
		calls.push_back(vector<Call>());
		arq >> nb_calls[k];
		for(int i = 0; i < nb_calls[k]; ++i){
			arq >> time >> lat >> longi >> hospital >> cleaning >> priority;
			arq >> hosp_needed >> clean_needed >> time_on_scene;
			arq >> time_at_hospital >> cleaning_time;
			Location call_location(lat, longi);
			calls.back().push_back(Call(i,time, call_location, hospital, 
				cleaning, priority, hosp_needed, clean_needed, time_on_scene, 
				time_at_hospital, cleaning_time));
			std::cout << time << " " << lat << " " << longi << " " << hospital << " ";
			std::cout << cleaning << " " << priority << " " << hosp_needed << " ";
			std::cout << clean_needed << " " << time_on_scene << " " << time_at_hospital;
			std::cout << " " << cleaning_time << "\n"; 
		}
	}


	for(int h = 0; h < nb_hospitals; ++h){
		arq >> lat >> longi;
		hospitals.push_back(std::make_pair(lat,longi));
		std::cout << lat << " " << longi << "\n";
	}

	int cap;
	for(int b = 0; b < nb_bases; ++b){
		arq >> lat >> longi >> cap;
		bases.push_back(std::make_pair(lat,longi));
		cap_bases.push_back(cap);
		std::cout << lat << " " << longi << " " << cap << "\n";
	}

	for(int cb = 0; cb < nb_cleaning_bases; ++cb){
		arq >> lat >> longi;
		cleaning_bases.push_back(std::make_pair(lat, longi));
		std::cout << lat << " " << longi << "\n";
	}
	
	int type;
	double speed;
	for(int a = 0; a < nb_ambulances; ++a){
		arq >> lat >> longi >> type >> speed;
		Location base_location(lat, longi);
		Location hospital_location = null_location;
		ambulances.push_back(Ambulance(a, base_location, hospital_location, -1, 0, 
			type, speed));
		std::cout << lat << " " << longi << " " << type << " " << speed << "\n";
	}

	if(path == "instances/no_myopic7.3.txt"){
		ambulances[2].free_location = std::make_pair(5,10);
		double t_to_b0 = travel.travel_time(ambulances[2].free_location,
			bases[0]);
		double t_to_b1 = travel.travel_time(ambulances[2].free_location,
			bases[1]);
		if(t_to_b0 < t_to_b1){
			ambulances[2].base_location = bases[0];
		}else{
			ambulances[2].base_location = bases[1];
		}
		ambulances[2].arrival_time_at_f_last_trip = 1;
		ambulances[2].arrival_time_at_b_last_trip = 1 + travel.travel_time(ambulances[2].free_location,
			ambulances[2].base_location);
		std::cout << ambulances[2] << "\n";
		ambulances[2].set_new_point(0, hospitals[0], TripType::TO_HOSPITAL);
		ambulances[2].set_new_point(1, hospitals[0], TripType::AT_HOSPITAL);
		if(t_to_b0 < t_to_b1){
			ambulances[2].set_new_point(1 + t_to_b0, ambulances[2].base_location, 
				TripType::TO_BASE);
		}else{
			ambulances[2].set_new_point(1 + t_to_b1, ambulances[2].base_location, 
				TripType::TO_BASE);
		}


	}else if(path == "instances/no_myopic7.4.txt"){
		ambulances[1].free_location = hospitals[0];
		ambulances[1].base_location = bases[0];
		ambulances[1].arrival_time_at_f_last_trip = 3;
		double t_to_b0 = travel.travel_time(ambulances[1].free_location,
			ambulances[1].base_location);
		ambulances[1].arrival_time_at_b_last_trip = 3 + t_to_b0;

		ambulances[1].set_new_point(0, hospitals[0], TripType::TO_HOSPITAL);
		ambulances[1].set_new_point(3, hospitals[0], TripType::AT_HOSPITAL);
		ambulances[1].set_new_point(3 + t_to_b0, ambulances[1].base_location, 
			TripType::TO_BASE);
	}

	// exit(1);
}

Instance::~Instance(){}



	// auto scenarios_arq = ofstream(fmt::format("{}/scenarios.txt", g_params.generator_folder),
	// 	ios::out);

	// for(int s = 0; s < nb_scenarios; ++s){
	// 	int id = 0;
	// 	for(int t = 0; t < nb_times; ++t){
	// 		std::uniform_int_distribution<int> int_dist(t*slot_duration, (t+1)*slot_duration);
	// 		std::uniform_int_distribution<int> sample_rand(0,99);
	// 		for(int g = 0; g < 1; ++g){
	// 			for(int r = 0; r < nb_regions; ++r){
	// 				for(int p = 0; p < nb_priorities; ++p){
	// 					poisson_distribution<int> dist(lambda(t,g,r,p));
	// 					int this_nb_calls = dist(gen);
	// 					nb_calls[s] += this_nb_calls;
	// 					for(int c = 0; c < this_nb_calls; ++c){
	// 						double time = int_dist(gen);
	// 						// Location location = samples[r][sample_rand(gen)];
	// 						Location location = centers[r];
	// 						double time_on_scene = 0.1*slot_duration;
	// 						double time_at_hospital = 0.2*slot_duration;
	// 						double cleaning_time = 1*slot_duration;
	// 						double min_time = GRB_INFINITY;
	// 						int min_h = nearest_hospital[r];
	// 						int cb = 0;
	// 						calls[s].push_back(Call(id++, time, location, min_h, cb, p,
	// 							true, false, time_on_scene, time_at_hospital, cleaning_time));
	// 						calls[s].back().region = r;
	// 						calls[s].back().disc_time = t;
	// 					}
	// 				}
	// 			}
	// 		}
	// 	}
	// 	std::sort(calls[s].begin(), calls[s].end());
	// 	for(int i = 0; i < calls[s].size(); ++i){
	// 		calls[s][i].id = i;
	// 	}
	// 	scenarios_arq << calls[s].size() << "\n";
	// 	int this_g = 5;
	// 	for(auto& call: calls[s]){
	// 		scenarios_arq << call.disc_time << " ";
	// 		scenarios_arq << call.region << " " << call.priority << " " << this_g << " ";
	// 		scenarios_arq << call.time_on_scene/slot_duration << " " << call.location.first;
	// 		scenarios_arq << " " << call.location.second << " ";
	// 		scenarios_arq << call.time_at_hospital/slot_duration << " ";
	// 		scenarios_arq << call.cleaning_time/slot_duration << " ";
	// 		scenarios_arq << call.clean_needed << " " << call.hosp_needed << "\n";
	// 		// fmt::print("{} {} ({}, {}) {} {}\n", call.id, call.time, 
	// 		// 	call.location.first, call.location.second, call.hosp_needed, 
	// 		// 	call.clean_needed);
	// 	}
	// }
	// scenarios_arq.close();

