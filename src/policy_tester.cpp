#include "../include/osrm_helper.h"
#include "../include/travel.h"
#include "../include/solver.h"
#include "../include/policy_tester.h"


PolicyTester::PolicyTester(Instance& ins): ins(ins), travel(ins.travel),
	ambulances(ins.ambulances), time(0),
	waiting_on_scene(vector<vector<double>>(ins.nb_scenarios, vector<double>())),
	waiting_on_scene_penalized(vector<vector<double>>(ins.nb_scenarios, vector<double>())),
	waiting_to_hospital(vector<vector<double>>(ins.nb_scenarios, vector<double>())),
	calls_end(vector<vector<double>>(ins.nb_scenarios, vector<double>())),
	which_ambulance(vector<vector<int>>(ins.nb_scenarios, vector<int>())){
	for(int i = 0; i < ins.nb_scenarios; ++i){
		waiting_on_scene[i] = vector<double>(ins.nb_calls[i], GRB_INFINITY);
		waiting_on_scene_penalized[i] = vector<double>(ins.nb_calls[i], GRB_INFINITY);
		waiting_to_hospital[i] = vector<double>(ins.nb_calls[i], GRB_INFINITY);
		calls_end[i] = vector<double>(ins.nb_calls[i], GRB_INFINITY);
		which_ambulance[i] = vector<int>(ins.nb_calls[i], -1);
	}
}

shared_ptr<Solver> PolicyTester::get_solver(const string& policy, vector<Call>& calls, 
	vector<Ambulance>& ambulances, Travel& travel){

	if(policy == "forward"){
		return static_pointer_cast<Solver>(make_shared<ForwardSolver>(env,calls,ambulances,
			ins, travel));
	}else if(policy == "queue"){
		return static_pointer_cast<Solver>(make_shared<QueueSolver>(env,calls,ambulances,
			ins, travel));
	}else if(policy == "priorities" || policy == "priorities_a"){
		auto solver = static_pointer_cast<Solver>(make_shared<PrioritySolver>(env,calls,ambulances,
			ins, travel));
		solver->travel.set_forward(policy == "priorities");
		return solver;
	}else if(policy == "minmax"){
		return static_pointer_cast<Solver>(make_shared<MinMaxSolver>(env,calls,ambulances,
			ins, travel));
	}else if(policy == "gen_forward"){
		return static_pointer_cast<Solver>(make_shared<GenForwardSolver>(env,calls,ambulances,
			ins, travel));
	}else if(policy == "gen_minmax"){
		return static_pointer_cast<Solver>(make_shared<GenMinMaxSolver>(env,calls,ambulances,
			ins, travel));
	}else if(policy == "cg"){
		return static_pointer_cast<Solver>(make_shared<CGSolver>(env,calls,ambulances,
			ins, travel));
	}else if(policy == "model"){
		return static_pointer_cast<Solver>(make_shared<ModelSolver>(env,calls,ambulances,
			ins, travel));
	}else if(policy == "minmaxp" || policy == "minmaxp_a"){
		auto solver = static_pointer_cast<Solver>(make_shared<MinMaxPSolver>(env,calls,ambulances,
			ins, travel));
		solver->travel.set_forward(policy == "minmaxp");
		return solver;
	}else if(policy == "non_miopyc"){
		return static_pointer_cast<Solver>(make_shared<NonMiopycSolver>(env,calls,ambulances,
			ins, travel));
	}

	cout << "ERROR: Unknown Policy: " << policy << "\n";
	exit(1);
}


void PolicyTester::run(){
	one_stage_old();
	// one_stage();
	fmt::print("=============== END Simulate Policies ====================\n");
	// two_stage();
	// two_stage_tree();
	// fmt::print("=============== END Two Stage ====================\n");
}


void PolicyTester::one_stage_old(){
	for(int s = 0; s < ins.calls.size(); ++s){
		for(int i = 0; i < ins.calls[s].size(); ++i){
			for(int a = 0; a < ambulances.size(); ++a){
				ins.calls[s][i].ambulances.push_back(a);
			}
		}
	}
	int nb_scen = g_params.n_scenarios;

	vector<double> all_waiting_on_scene;
	vector<double> all_waiting_on_scene_penalized;
	vector<double> all_waiting_to_hospital;
	all_waiting_on_scene.reserve(ins.nb_scenarios*ins.calls[0].size());
	all_waiting_on_scene_penalized.reserve(ins.nb_scenarios*ins.calls[0].size());
	all_waiting_to_hospital.reserve(ins.nb_scenarios*ins.calls[0].size());
	for(auto& policy: policies){
		int s = 0;
		all_waiting_on_scene.clear();
		all_waiting_on_scene_penalized.clear();
		all_waiting_to_hospital.clear();
		for(int s = 0; s < ins.calls.size(); ++s){
			waiting_on_scene[s].clear();
			waiting_on_scene_penalized[s].clear();
		}
		run_times.clear();
		fmt::print("Policy {}\n", policy);
		for(auto& scenario: ins.calls){
			// fmt::print("Scenario {} ({} calls)\n", s, scenario.size());
			auto solver = get_solver(policy, scenario, ins.ambulances, travel);
			solver->run();
			// solver->print_results();
			run_times.insert(run_times.end(), solver->run_times.begin(), 
				solver->run_times.end());
			for(int i = 0; i < scenario.size(); ++i){
				waiting_on_scene[s].push_back(solver->waiting_on_scene[i]);
				waiting_on_scene_penalized[s].push_back(solver->waiting_on_scene_penalized[i]);
			}
			all_waiting_on_scene.insert(all_waiting_on_scene.end(), 
				solver->waiting_on_scene.begin(), solver->waiting_on_scene.end());
			all_waiting_on_scene_penalized.insert(all_waiting_on_scene_penalized.end(), 
				solver->waiting_on_scene_penalized.begin(), 
				solver->waiting_on_scene_penalized.end());
			all_waiting_to_hospital.insert(all_waiting_to_hospital.end(), 
				solver->waiting_to_hospital.begin(), solver->waiting_to_hospital.end());
			++s;
		}

		for(int g = 0; g < 7; ++g){
			const string policy_name = policies_names.find(policy)->second;
			ofstream daily(fmt::format("daily_results/g{}/{}/{}.txt", g, policy_name, 
				policy_name), ios::out);
			ofstream daily_pen(fmt::format("daily_results/g{}/{}/{}_pen.txt", g, policy_name, 
				policy_name), ios::out);
			for(int t = 0; t < ins.nb_times; ++t){
				for(int s = 0; s < nb_scen; ++s){
					int nb_t_calls = 0;
					int s_ind = g*nb_scen+s;
					for(int i = 0; i < ins.calls[s_ind].size(); ++i){
						if(ins.calls[s_ind][i].time > (t+1)*1800){
							break;
						}
						++nb_t_calls;
					}
					daily << t << " ";
					daily_pen << t << " ";
					for(int i = 0; i < nb_t_calls; ++i){
						daily << waiting_on_scene[s_ind][i] << " ";
						daily_pen << waiting_on_scene_penalized[s_ind][i] << " ";
					}
					daily << "\n";
					daily_pen << "\n";
				}
			}
			daily.close();
			daily_pen.close();
		}
	
		Stats stats(all_waiting_on_scene, all_waiting_on_scene_penalized, all_waiting_to_hospital);
		fmt::print("Policy {} Mean time = {:.1f} Mean pen = {:.1f}, Max pen = {:.1f} \
		Total pen = {:.1f}\n", policy, stats.mean_waiting_on_scene,
			stats.mean_waiting_on_scene_penalized, stats.max_waiting_on_scene_penalized,
			accumulate(all_waiting_on_scene_penalized.begin(),
				all_waiting_on_scene_penalized.end(), 0.0));
		double min_run_time = *min_element(run_times.begin(), run_times.end());
		double avg_run_time = accumulate(run_times.begin(),
			run_times.end(), 0.0)/run_times.size();
		double max_run_time = *max_element(run_times.begin(), run_times.end());
		fmt::print("Run times: {}\t{}\t{}\n",min_run_time,avg_run_time,max_run_time);
		// std::cin.get();
	}

}



void PolicyTester::one_stage(){

	for(int s = 0; s < ins.calls.size(); ++s){
		for(int i = 0; i < ins.calls[s].size(); ++i){
			for(int a = 0; a < ambulances.size(); ++a){
				ins.calls[s][i].ambulances.push_back(a);
			}
		}
	}

	int nb_scen_simu = 100;
	double t0 = 36*1800;
	int T = 4;
	int nb_realizations = 1;
	vector<vector<vector<Call>>> my_scenarios(T, 
		vector<vector<Call>>(nb_realizations, vector<Call>()));

	for(int scen = 0; scen < nb_realizations; ++scen){
		int j = 0;
		for(int t = 0; t < T; ++t){
			while(j < ins.calls[scen].size() && 
				ins.calls[scen][j].time <= t0 + (t+1)*1800){
				my_scenarios[t][scen].push_back(ins.calls[scen][j]);
				++j;
			}
		}
	}

	// for(int t = 0; t < T; ++t){
	// 	for(int scen = 0; scen < nb_realizations; ++scen){
	// 		my_scenarios[t][scen] = vector<Call>(my_scenarios[t][scen].begin(),
	// 			my_scenarios[t][scen].begin()+2);
	// 	}
	// }


	vector<double> all_waiting_on_scene;
	vector<double> all_waiting_on_scene_penalized;
	vector<double> all_waiting_to_hospital;
	all_waiting_on_scene.reserve(ins.nb_scenarios*ins.calls[0].size());
	all_waiting_on_scene_penalized.reserve(ins.nb_scenarios*ins.calls[0].size());
	all_waiting_to_hospital.reserve(ins.nb_scenarios*ins.calls[0].size());
	for(auto& policy: policies){
		// int s = 0;
		all_waiting_on_scene.clear();
		all_waiting_on_scene_penalized.clear();
		all_waiting_to_hospital.clear();
		run_times.clear();
		fmt::print("Policy {}\n", policy);
		// for(auto& scenario: ins.calls){
		for(int s = 0; s < 1; ++s){
			vector<Call> this_scenario;
			vector<int> this_scenario_indexes;
			default_random_engine gen;
			for(int t = 0; t < T; ++t){
				uniform_int_distribution<int> rand_scen(0,nb_realizations-1);
				int u = rand_scen(gen);
				// int u = 0;
				this_scenario_indexes.push_back(u);
				this_scenario.insert(this_scenario.end(), my_scenarios[t][u].begin(),
					my_scenarios[t][u].end());
			}

			for(int i = 0; i < this_scenario.size(); ++i){
				this_scenario[i].id = i;
			}

			waiting_on_scene[s] = vector<double>(this_scenario.size(), GRB_INFINITY);
			waiting_on_scene_penalized[s] = vector<double>(this_scenario.size(), GRB_INFINITY);
			waiting_to_hospital[s] = vector<double>(this_scenario.size(), GRB_INFINITY);
			calls_end[s] = vector<double>(this_scenario.size(), GRB_INFINITY);
			which_ambulance[s] = vector<int>(this_scenario.size(), -1);
			// fmt::print("Scenario {} ({} calls)\n", s, scenario.size());
			auto solver = get_solver(policy, this_scenario, ins.ambulances, travel);
			solver->run();
			// solver->print_results();
			run_times.insert(run_times.end(), solver->run_times.begin(), 
				solver->run_times.end());
			all_waiting_on_scene.insert(all_waiting_on_scene.end(), 
				solver->waiting_on_scene.begin(), solver->waiting_on_scene.end());
			all_waiting_on_scene_penalized.insert(all_waiting_on_scene_penalized.end(), 
				solver->waiting_on_scene_penalized.begin(), 
				solver->waiting_on_scene_penalized.end());
			all_waiting_to_hospital.insert(all_waiting_to_hospital.end(), 
				solver->waiting_to_hospital.begin(), solver->waiting_to_hospital.end());
			++s;
		}

		Stats stats(all_waiting_on_scene, all_waiting_on_scene_penalized, all_waiting_to_hospital);
		fmt::print("Policy {} Mean pen = {:.1f}, Max pen = {:.1f} Total pen = {:.1f}\n", policy,
			stats.mean_waiting_on_scene_penalized, stats.max_waiting_on_scene_penalized,
			accumulate(all_waiting_on_scene_penalized.begin(),
				all_waiting_on_scene_penalized.end(), 0.0));
		double min_run_time = *min_element(run_times.begin(), run_times.end());
		double avg_run_time = accumulate(run_times.begin(),
			run_times.end(), 0.0)/run_times.size();
		double max_run_time = *max_element(run_times.begin(), run_times.end());
		fmt::print("Run times: {}\t{}\t{}\n",min_run_time,avg_run_time,max_run_time);
		// std::cin.get();
	}
}


void PolicyTester::two_stage(){

	vector<double> all_waiting_on_scene;
	vector<double> all_waiting_on_scene_penalized;
	vector<double> all_waiting_to_hospital;

	bool debug = false;

	for(int s = 0; s < ins.calls.size(); ++s){
		for(int i = 0; i < ins.calls[s].size(); ++i){
			for(int a = 0; a < ambulances.size(); ++a){
				ins.calls[s][i].ambulances.push_back(a);
			}
		}
	}

	for(auto& policy: policies){ // each policy
		int sc = 0;
		fmt::print("Policy: {}\n", policy);
		all_waiting_on_scene.clear();
		all_waiting_on_scene_penalized.clear();
		all_waiting_to_hospital.clear();
		run_times.clear();
		for(int i = 0; i < ins.nb_scenarios; ++i){
			which_ambulance[i] = std::vector<int>(ins.nb_calls[i], -1);
			waiting_on_scene[i] = std::vector<double>(ins.nb_calls[i], GRB_INFINITY);
			waiting_on_scene_penalized[i] = std::vector<double>(ins.nb_calls[i], GRB_INFINITY);
			waiting_to_hospital[i] = std::vector<double>(ins.nb_calls[i], GRB_INFINITY);
		}

		// each scenario
		for(auto& scenario: ins.calls){
			auto nearest_base = get_nearest_base(scenario);
			ambulances = ins.ambulances;
			int calls_attended = 0;
			int event_call = 1;
			int index_call = 0;
			time = scenario[0].time;
			int amb_finish = -1;

			std::vector<int> queue; 
			queue.reserve(scenario.size());
			int iter_count = 0;
			while(calls_attended < scenario.size()){
				if(debug)
					fmt::print("calls_attended {} index_call {}\n",calls_attended, 
						index_call);
				vector<int> index_scenarios(ins.calls.size(), 0);
				//Find future calls in each scenario.
				for(int j = 0; j < ins.calls.size(); ++j){
					while(index_scenarios[j] < ins.calls[j].size() && 
						ins.calls[j][index_scenarios[j]].time <= time){
						++index_scenarios[j];
					}
				}

				if(event_call == 1){
					queue.push_back(index_call);
					std::vector<int> queue_aux;
					queue_aux.reserve(queue.size());
					bool increased_ca = false;
					for(int i = 0; i < queue.size(); ++i){
						auto t0 = std::chrono::high_resolution_clock::now();
						auto& call = scenario[queue[i]];
						double min_time = GRB_INFINITY;
						int index_amb = -1;
						double min_waiting_on_scene = GRB_INFINITY;
						double min_waiting_to_hospital = GRB_INFINITY;
						if(debug)
							cout << "Call " << call << "\n";
						for(int k = 0; k < ambulances.size(); ++k){
							auto& ambulance = ambulances[k];
							if(ambulance.type <= call.priority && 
								ambulance.arrival_time_at_f_last_trip <= time){
								auto result = get_waiting_time(ambulance.id, call, policy,
									index_scenarios, nearest_base[queue[i]],queue,sc, time);
								if(debug){
									cout << "\t" << ambulance << " " << result.mean_total_time;
									cout << "\n";
								}
								if(result.mean_total_time < min_time){
									index_amb = k;
									min_time = result.mean_total_time;
									min_waiting_on_scene = result.waiting_on_scene;
									min_waiting_to_hospital = result.waiting_to_hospital;
								}
							}
						}

						auto result = get_waiting_time(-1,call,policy,index_scenarios,
							nearest_base[queue[i]],queue,sc, time);

						if(debug)
							cout << "\tqueue " << result.mean_total_time << "\n";
						if(result.mean_total_time < min_time){
							index_amb = -1;
							min_time = result.mean_total_time;
							min_waiting_on_scene = result.waiting_on_scene;
							min_waiting_to_hospital = result.waiting_to_hospital;
						}

						if(index_amb >= 0){
							min_time = travel.get_response_time(ambulances[index_amb],
								call,time,true);
							double waiting_on_scene_i = time + min_time - call.time;
							double waiting_to_hospital_i = ambulances[index_amb].answer_call(
								call, travel, ins, time, min_time, nearest_base[queue[i]]);
							which_ambulance[sc][queue[i]] = index_amb;
							waiting_on_scene[sc][queue[i]] = waiting_on_scene_i;
							waiting_on_scene_penalized[sc][queue[i]] = waiting_on_scene_i*
								ins.penalties[call.priority];
							waiting_to_hospital[sc][queue[i]] = waiting_to_hospital_i;
							calls_attended++;
							increased_ca = true;
						}else{
							if(debug)
								cout << "queue\n";
							for(int a = 0; a < ambulances.size(); ++a){
								if(ambulances[a].arrival_time_at_f_last_trip <= time){
									call.ambulances.erase(remove(call.ambulances.begin(),
										call.ambulances.end(), a), call.ambulances.end());
								}
							}
							queue_aux.push_back(queue[i]);
						}
						auto dt = std::chrono::high_resolution_clock::now();
						run_times.push_back(std::chrono::duration_cast<chrono::nanoseconds>(dt - t0).count() 
							/ pow(10,9));
					}
					queue.clear();
					for(auto i: queue_aux){
						queue.push_back(i);
					}

					if(increased_ca){
						iter_count = 0;
					}else{
						iter_count++;
					}
				}else{
					//Event_call == 0, get ambulance that just finished.
					if(amb_finish == -1){
						fmt::print("ERROR! event_call = {}, but no ambulance is returning",
							event_call);
						cin.get();
					}
					auto& amb = ambulances[amb_finish];
					if(debug){
						cout << "amb_finish = " << amb_finish << "\n";
						cout << "return " << amb << "\n";
					}
					std::vector<int> queue_aux;
					queue_aux.reserve(queue.size());

					int ind_call  = -1;
					int ind_base = -1;
					double min_time = GRB_INFINITY;
					double min_waiting_on_scene = GRB_INFINITY;
					double min_waiting_to_hospital = GRB_INFINITY;
					for(int i = 0; i < queue.size(); ++i){
						auto& call = scenario[queue[i]];
						if(amb.type <= call.priority){
							auto result = get_waiting_time(amb_finish, call, policy,
								index_scenarios, nearest_base[queue[i]],queue,sc,time);
							if(debug)
								cout << "\tCall "<<call<< " " << result.mean_total_time << "\n";
							if(result.mean_total_time < min_time){
								ind_call = queue[i];
								min_time = result.mean_total_time;
								min_waiting_on_scene = result.waiting_on_scene;
								min_waiting_to_hospital = result.waiting_to_hospital;
							}
						}
					}

					// for(int b = 0; b < ins.bases.size(); ++b){
					// 	auto result = get_waiting_time_return(amb_finish, policy,
					// 				index_scenarios, time,queue,sc,b);

					// 	if(debug){
					// 		cout << "\tbase " << b << " " << result.mean_total_time << "\n";
					// 	}

					// 	if(result.mean_total_time < min_time){
					// 		ind_call = -1;
					// 		ind_base = b;
					// 		min_time = result.mean_total_time;
					// 		min_waiting_on_scene = result.waiting_on_scene;
					// 		min_waiting_to_hospital = result.waiting_to_hospital;
					// 	}

					// }


					if(ind_call >= 0){
						auto& call = scenario[ind_call];
						double min_time = travel.get_response_time(amb, call,
							time, true);
						double waiting_on_scene_i = time + min_time - call.time;
						calls_attended++;
						iter_count = 0;
						double waiting_to_hospital_i = amb.answer_call(call,
							travel, ins, time, min_time, nearest_base[ind_call]);
						// if(debug)
						// 	cout << "Amb after " << amb << " " << waiting_on_scene_i*
						// 	ins.penalties[call.priority];
						which_ambulance[sc][ind_call] = amb_finish;
						waiting_on_scene[sc][ind_call] = waiting_on_scene_i;
						waiting_on_scene_penalized[sc][ind_call] = waiting_on_scene_i*
							ins.penalties[call.priority];
						waiting_to_hospital[sc][ind_call] = waiting_to_hospital_i;
						queue.erase(remove(queue.begin(), queue.end(), ind_call),
							queue.end());
					}
					// else if(ind_base >= 0){
					// 	auto base = ins.bases[ind_base];
					// 	amb.base_location = base;
					// 	amb.arrival_time_at_b_last_trip = time + travel.travel_time(
					// 		amb.free_location,amb.base_location);

					// 	fmt::print("Chosen Base {}\n",ind_base);

					// 	for(auto i: queue){
					// 		scenario[i].ambulances.erase(remove(
					// 			scenario[i].ambulances.begin(),
					// 			scenario[i].ambulances.end(), amb_finish),
					// 			scenario[i].ambulances.end());
					// 	}

					// 	iter_count++;
					// }
				}
				amb_finish = set_next_event(scenario, event_call, index_call);
				// if(debug)
				// 	cin.get();
			}
		}

		int k = 0;
		for(int s = 0; s < ins.nb_scenarios; ++s){
			for(int i = 0; i < ins.calls[s].size(); ++i){
				// fmt::print("{} {} {}\n",s,i,waiting_on_scene_penalized[s][i]);
				all_waiting_on_scene.push_back(waiting_on_scene[s][i]);
				all_waiting_on_scene_penalized.push_back(waiting_on_scene_penalized[s][i]);
				all_waiting_to_hospital.push_back(waiting_to_hospital[s][i]);
			}
		}

		Stats stats(all_waiting_on_scene, all_waiting_on_scene_penalized, all_waiting_to_hospital);
		fmt::print("Policy {} Mean pen = {:.1f}, Max pen = {:.1f} Total pen = {:.1f}\n", policy,
			stats.mean_waiting_on_scene_penalized, stats.max_waiting_on_scene_penalized,
			accumulate(all_waiting_on_scene_penalized.begin(),
				all_waiting_on_scene_penalized.end(), 0.0));
		double min_run_time = *min_element(run_times.begin(), run_times.end());
		double avg_run_time = accumulate(run_times.begin(),
			run_times.end(), 0.0)/run_times.size();
		double max_run_time = *max_element(run_times.begin(), run_times.end());
		fmt::print("Run times: {}\t{}\t{}\n",min_run_time,avg_run_time,max_run_time);

		std::string base_filename = g_params.instance.substr(
			g_params.instance.find_last_of("/\\") + 1);
		std::string::size_type const p(base_filename.find_last_of('.'));
		std::string file_without_extension = base_filename.substr(0, p);
		std::ofstream results(fmt::format("results/{}_{}_{}.txt", 
			file_without_extension, policy, ins.nb_scenarios), std::ios::out);
		for(int i = 0; i < ins.nb_scenarios; ++i){
			results << ins.nb_calls[i] << " ";
			for(int j = 0; j < ins.nb_calls[i]; ++j){
				results << waiting_on_scene[i][j] << " ";
			}
			results << "\n";
		}
		results.close();
		// std::cin.get();
	}	
}


void PolicyTester::two_stage_tree(){
	vector<double> all_waiting_on_scene;
	vector<double> all_waiting_on_scene_penalized;
	vector<double> all_waiting_to_hospital;

	bool debug = false;

	for(int s = 0; s < ins.calls.size(); ++s){
		for(int i = 0; i < ins.calls[s].size(); ++i){
			for(int a = 0; a < ambulances.size(); ++a){
				ins.calls[s][i].ambulances.push_back(a);
			}
		}
	}


	int nb_scen_simu = 100;
	double t0 = 36*1800;
	int T = 4;
	int nb_realizations = 1;
	vector<vector<vector<Call>>> my_scenarios(T, 
		vector<vector<Call>>(nb_realizations, vector<Call>()));

	for(int scen = 0; scen < nb_realizations; ++scen){
		int j = 0;
		for(int t = 0; t < T; ++t){
			while(j < ins.calls[scen].size() && 
				ins.calls[scen][j].time <= t0 + (t+1)*1800){
				my_scenarios[t][scen].push_back(ins.calls[scen][j]);
				++j;
			}
		}
	}

	// for(int t = 0; t < T; ++t){
	// 	for(int scen = 0; scen < nb_realizations; ++scen){
	// 		my_scenarios[t][scen] = vector<Call>(my_scenarios[t][scen].begin(),
	// 			my_scenarios[t][scen].begin()+2);
	// 	}
	// }

	default_random_engine gen;
	uniform_int_distribution<int> rand_scen(0,nb_realizations-1);
	for(auto& policy: policies){
		fmt::print("Policy: {}\n", policy);
		run_times.clear();
		// for(int simu = 0; simu < g_params.n_scenarios; ++simu){
		for(int simu = 0; simu < 1; ++simu){
			vector<Call> this_scenario;
			vector<int> this_scenario_indexes;
			for(int t = 0; t < T; ++t){
				int u = rand_scen(gen);
				// int u = 0;
				this_scenario_indexes.push_back(u);
				this_scenario.insert(this_scenario.end(), my_scenarios[t][u].begin(),
					my_scenarios[t][u].end());
			}

			for(int i = 0; i < this_scenario.size(); ++i){
				this_scenario[i].id = i;
			}

			waiting_on_scene[simu] = vector<double>(this_scenario.size(), GRB_INFINITY);
			waiting_on_scene_penalized[simu] = vector<double>(this_scenario.size(), GRB_INFINITY);
			waiting_to_hospital[simu] = vector<double>(this_scenario.size(), GRB_INFINITY);
			calls_end[simu] = vector<double>(this_scenario.size(), GRB_INFINITY);
			which_ambulance[simu] = vector<int>(this_scenario.size(), -1);

			auto nearest_base = get_nearest_base(this_scenario);
			ambulances = ins.ambulances;
			int calls_attended = 0;
			int event_call = 1;
			int index_call = 0;
			time = this_scenario[0].time;
			int amb_finish = -1;

			std::vector<int> queue; 
			queue.reserve(this_scenario.size());
			int iter_count = 0;
			while(calls_attended < this_scenario.size()){
				fmt::print("Runtimes.size() = {}\n",run_times.size());
				double min_run_time = *min_element(run_times.begin(), run_times.end());
				double avg_run_time = accumulate(run_times.begin(),
					run_times.end(), 0.0)/run_times.size();
				double max_run_time = *max_element(run_times.begin(), run_times.end());
				fmt::print("Run times: {}\t{}\t{}\n",min_run_time,avg_run_time,
					max_run_time);

				vector<vector<Call>> future_scenarios = get_future_scenarios(
					time, this_scenario, nb_realizations, T, my_scenarios);
				if(debug){
					fmt::print("Time = {}, index_call = {}, calls_attended = {}\n", time,
						index_call, calls_attended);
				}
				if(event_call == 1){
					queue.push_back(index_call);
					std::vector<int> queue_aux;
					queue_aux.reserve(queue.size());
					for(int i = 0; i < queue.size(); ++i){
						auto t0 = std::chrono::high_resolution_clock::now();
						auto& call = this_scenario[queue[i]];
						double min_time = GRB_INFINITY;
						int index_amb = -1;
						double min_waiting_on_scene = GRB_INFINITY;
						double min_waiting_to_hospital = GRB_INFINITY;
						
						if(debug){
							cout << "Call " << call << "\n";
						}

						for(int k = 0; k < ambulances.size(); ++k){
							auto& ambulance = ambulances[k];
							if(can_answer(ambulance, call) && 
								ambulance.arrival_time_at_f_last_trip <= time){
								auto result = get_waiting_time_tree(ambulance.id, call.id,
								 	policy, future_scenarios, this_scenario, time, 
								 	queue, nearest_base[queue[i]]);
								if(debug){
									cout << "\t" << ambulance << " => " << result.mean_total_time;
									cout << "\n";
								}
								if(result.mean_total_time < min_time){
									index_amb = k;
									min_time = result.mean_total_time;
									min_waiting_on_scene = result.waiting_on_scene;
									min_waiting_to_hospital = result.waiting_to_hospital;
								}
							}
						}
						auto result = get_waiting_time_tree(-1,call.id,policy,
							future_scenarios, this_scenario, time, queue,
							nearest_base[queue[i]]);
						if(debug){
							cout << "\tQueue =>" << result.mean_total_time << "\n"; 
						}
						if(result.mean_total_time < min_time){
							index_amb = -1;
							min_time = result.mean_total_time;
							min_waiting_on_scene = result.waiting_on_scene;
							min_waiting_to_hospital = result.waiting_to_hospital;
						}

						if(index_amb >= 0){
							min_time = travel.get_response_time(ambulances[index_amb],
								call,time,true);
							double waiting_on_scene_i = time + min_time - call.time;
							double waiting_to_hospital_i = ambulances[index_amb].answer_call(
								call, travel, ins, time, min_time, nearest_base[queue[i]]);
							which_ambulance[simu][queue[i]] = index_amb;
							waiting_on_scene[simu][queue[i]] = waiting_on_scene_i;
							waiting_on_scene_penalized[simu][queue[i]] = 
								waiting_on_scene_i * ins.penalties[call.priority];
							if(debug){
								cout << "Chosen " << ambulances[index_amb] << " pt ";
								cout << waiting_on_scene_penalized[simu][queue[i]] << "\n";
							}
							waiting_to_hospital[simu][queue[i]] = waiting_to_hospital_i;
							calls_attended++;
						}else{
							queue_aux.push_back(queue[i]);
							if(debug){
								cout << "Chose Queue\n";
							}
						}
						auto dt = std::chrono::high_resolution_clock::now();
						run_times.push_back(std::chrono::duration_cast<chrono::nanoseconds>(dt - t0).count() 
							/ pow(10,9));
					}
					queue.clear();
					for(auto i: queue_aux){
						queue.push_back(i);
					}
				}else{
					//Event_call == 0, get ambulance that just finished.
					if(amb_finish == -1){
						fmt::print("WEIRD!!! event_call = 0 but no return\n");
						fmt::print("Runtimes.size() = {}\n",run_times.size());
						double min_run_time = *min_element(run_times.begin(), run_times.end());
						double avg_run_time = accumulate(run_times.begin(),
							run_times.end(), 0.0)/run_times.size();
						double max_run_time = *max_element(run_times.begin(), run_times.end());
						fmt::print("Run times: {}\t{}\t{}\n",min_run_time,avg_run_time,
							max_run_time);
						cin.get();
					}
					auto& amb = ambulances[amb_finish];
					std::vector<int> queue_aux;
					queue_aux.reserve(queue.size());
					if(debug){
						cout << "Return " << amb << "\n";
					}
					int ind_call  = -1;
					int ind_base = -1;
					double min_time = GRB_INFINITY;
					double min_waiting_on_scene = GRB_INFINITY;
					double min_waiting_to_hospital = GRB_INFINITY;
					for(int i = 0; i < queue.size(); ++i){
						auto t0 = std::chrono::high_resolution_clock::now();
						auto& call = this_scenario[queue[i]];
						if(can_answer(amb,call)){
							auto result = get_waiting_time_tree(amb_finish, call.id, policy,
								future_scenarios, this_scenario, time, queue, 
								nearest_base[queue[i]]);
							if(debug){
								cout << "\tCall " << call << " => " << result.mean_total_time;
								cout << "\n";
							}
							if(result.mean_total_time < min_time){
								ind_call = queue[i];
								min_time = result.mean_total_time;
								min_waiting_on_scene = result.waiting_on_scene;
								min_waiting_to_hospital = result.waiting_to_hospital;
							}
						}
						auto dt = std::chrono::high_resolution_clock::now();
						run_times.push_back(std::chrono::duration_cast<chrono::nanoseconds>(dt - t0).count() 
							/ pow(10,9));
					}
					
					for(int b = 0; b < ins.bases.size(); ++b){
						auto result = get_waiting_time_tree_return(amb_finish, b,
							policy, future_scenarios, this_scenario, time, queue);
						if(debug){
							cout << "\tBase " << b << " => " << result.mean_total_time << "\n";
						}
						if(result.mean_total_time < min_time){
							min_time = result.mean_total_time;
							ind_call = -1;
							ind_base = b;
						}
					}


					if(ind_call >= 0){
						auto& call = this_scenario[ind_call];
						double min_time = travel.get_response_time(amb, call,
							time, true);
						double waiting_on_scene_i = time + min_time - call.time;
						calls_attended++;
						iter_count = 0;
						double waiting_to_hospital_i = amb.answer_call(call,
							travel, ins, time, min_time, nearest_base[ind_call]);

						which_ambulance[simu][ind_call] = amb_finish;
						waiting_on_scene[simu][ind_call] = waiting_on_scene_i;
						waiting_on_scene_penalized[simu][ind_call] = waiting_on_scene_i*
							ins.penalties[call.priority];
						if(debug){
							cout << "Chosen call " << call << " ";
							cout << waiting_on_scene_penalized[simu][ind_call] << "\n";
						}
						waiting_to_hospital[simu][ind_call] = waiting_to_hospital_i;
						queue.erase(remove(queue.begin(), queue.end(), ind_call),
							queue.end());
					}else if(ind_base >= 0){
						// cout << "Chose Return\n";
						auto base = ins.bases[ind_base];
						amb.base_location = base;
						amb.arrival_time_at_b_last_trip = time + travel.travel_time(
							amb.free_location,amb.base_location);

						if(debug){
							fmt::print("Chosen Base {}\n",ind_base);
						}
						
						// If ambulance must return to base, it can't answer any calls in queue.
						for(auto i: queue){
							this_scenario[i].ambulances.erase(remove(
								this_scenario[i].ambulances.begin(),
								this_scenario[i].ambulances.end(), amb_finish),
								this_scenario[i].ambulances.end());
						}
					}

				}
				amb_finish = set_next_event(this_scenario, event_call, index_call);
				if(debug){
					cout << "==================================\n";
					cin.get();
				}
			}
			// cout << "Call\tWait On Scene(s)\tWait Penalized(s)\tAmb:\n";
			for(int i = 0; i < this_scenario.size(); ++i){
				// auto& call = this_scenario[i];
				// cout << call << "\t"<< waiting_on_scene[simu][i] << "\t";
				// cout << waiting_on_scene_penalized[simu][i] << "\t";
				// cout << which_ambulance[simu][i] << "\n";
				all_waiting_on_scene.push_back(waiting_on_scene[simu][i]);
				all_waiting_on_scene_penalized.push_back(waiting_on_scene_penalized[simu][i]);
				all_waiting_to_hospital.push_back(waiting_to_hospital[simu][i]);
			}
		}

		Stats stats(all_waiting_on_scene, all_waiting_on_scene_penalized, all_waiting_to_hospital);
		fmt::print("Policy {} Mean pen = {:.1f}, Max pen = {:.1f} Total pen = {:.1f}\n", policy,
			stats.mean_waiting_on_scene_penalized, stats.max_waiting_on_scene_penalized,
			accumulate(all_waiting_on_scene_penalized.begin(),
				all_waiting_on_scene_penalized.end(), 0.0));
		double min_run_time = *min_element(run_times.begin(), run_times.end());
		double avg_run_time = accumulate(run_times.begin(),
			run_times.end(), 0.0)/run_times.size();
		double max_run_time = *max_element(run_times.begin(), run_times.end());
		fmt::print("Run times: {}\t{}\t{}\n",min_run_time,avg_run_time,max_run_time);
		// std::cin.get();
	}
}

SecondStageWaitingTime PolicyTester::get_waiting_time_tree(int amb_id,  int call_id,
		const string& policy, vector<vector<Call>>& future_scenarios, 
		vector<Call>& this_scenario, double time, 
		std::vector<int>& queue, int base){
	Call* call = NULL;
	if(call_id >= 0){
		call = &this_scenario[call_id];
	}
	vector<Ambulance> temp_ambulances = ambulances;
	double min_time = 0;
	double waiting_time = 0;
	double waiting_to_hospital = 0;
	if(amb_id >= 0 && call_id >= 0){
		Ambulance& amb = temp_ambulances[amb_id];
		min_time = travel.get_response_time(amb, *call, time, true);
		waiting_to_hospital = amb.answer_call(*call, travel, ins, time, min_time, 
			base);
		waiting_time = time + min_time - call->time;
	}

	double sum_total = 0;
	int num_total = 0;

	vector<int> busy_ambulances;
	if(amb_id == -1){
		for(auto& amb: temp_ambulances){
			if(can_answer(amb,*call) && amb.arrival_time_at_f_last_trip > time){
				busy_ambulances.push_back(amb.id);
			}
		}

		if(busy_ambulances.size() == 0){
			return SecondStageWaitingTime{0, 0, GRB_INFINITY};
		}
	}


	for(int j = 0; j < future_scenarios.size(); ++j){
		auto future_calls = future_scenarios[j];
		vector<Ambulance> aux_ambulances = temp_ambulances;
	
		if(!queue.empty()){
			//if queue isn't empty, calls of queue are added to future scenarios.
			for(int i = queue.size()-1; i >= 0; --i){
				if(amb_id == -1 || queue[i] != call_id){ //if call in queue is not the current call
					future_calls.insert(future_calls.begin(), this_scenario[queue[i]]);
					if(queue[i] == call_id){
						// If the current call must go to queue, 
						// all free ambulances are excluded from attending it.
						for(auto& amb: aux_ambulances){
							if(amb.arrival_time_at_f_last_trip <= time){
								future_calls[0].ambulances.erase(
									remove(future_calls[0].ambulances.begin(),
										future_calls[0].ambulances.end(), amb.id),
										future_calls[0].ambulances.end());
							}
						}
					}
					// if(amb_id != -1 && call_id == -1){
					// 	// If current ambulance must return to base:
					// 	// current_ambulance can't be assigned to calls in queue
					// 	future_calls[0].ambulances.erase(remove(
					// 		future_calls[0].ambulances.begin(), future_calls[0].ambulances.end(), 
					// 		amb_id), future_calls[0].ambulances.end());
					// }
				}
			}
		}

		if(!future_calls.empty()){
			for(int i = 0; i < future_calls.size(); ++i){
				future_calls[i].id = i;
			}
			auto solver = get_solver(policy, future_calls, aux_ambulances, travel);
			solver->run();

			double sum_scenario = 0;
			int num_scenario = 0;

			if(policy != "cg"){
				sum_scenario = accumulate(solver->waiting_on_scene_penalized.begin(),
					solver->waiting_on_scene_penalized.end(), 0.0);
				sum_total += sum_scenario;
				// sum_total += accumulate(solver->waiting_to_hospital.begin(),
				// 	solver->waiting_to_hospital.end(), 0.0);
				num_scenario = solver->waiting_on_scene_penalized.size();
				num_total += num_scenario;
			}else{
				sum_total += solver->obj;
				++num_total;
			}
		}

		if(amb_id != -1 && call_id != -1){
			sum_total += waiting_time*ins.penalties[call->priority];
			num_total += 1;
		}
	}
	// cout << "call0 pt = " << waiting_time*ins.penalties[call->priority] << "\n";

	double result_total = (num_total == 0) ? 0 : sum_total / num_total;
	return SecondStageWaitingTime{0, 0, result_total};
}


SecondStageWaitingTime PolicyTester::get_waiting_time_tree_return(int amb_id, int base,
	const string& policy, vector<vector<Call>>& future_scenarios, 
	vector<Call>& this_scenario, double time, 
	std::vector<int>& queue){
	
	vector<Ambulance> temp_ambulances = ambulances;
	auto& amb = temp_ambulances[amb_id];
	amb.base_location = ins.bases[base];
	amb.arrival_time_at_b_last_trip = time + travel.travel_time(
		amb.free_location, amb.base_location);
	
	double min_time = 0;
	double waiting_to_hospital = 0;
	double sum_total = 0;
	int num_total = 0;

	for(int j = 0; j < future_scenarios.size(); ++j){
		auto future_calls = future_scenarios[j];
		
		vector<Ambulance> aux_ambulances = temp_ambulances;
		for(int i = queue.size()-1; i >= 0; --i){
			future_calls.insert(future_calls.begin(), this_scenario[queue[i]]);
			// By returning the current ambulance to base,
			// we forbid it to answer any calls in queue:
			future_calls[0].ambulances.erase(remove(future_calls[0].ambulances.begin(),
				future_calls[0].ambulances.end(), amb_id),
				future_calls[0].ambulances.end());
		}
		if(!future_calls.empty()){
			for(int i = 0; i < future_calls.size(); ++i){
				future_calls[i].id = i;
			}
			auto solver = get_solver(policy, future_calls, aux_ambulances, travel);
			solver->run();

			double sum_scenario = 0;
			int num_scenario = 0;
			if(policy != "cg"){
				sum_scenario = accumulate(solver->waiting_on_scene_penalized.begin(),
					solver->waiting_on_scene_penalized.end(), 0.0);
				sum_total += sum_scenario;
				// sum_total += accumulate(solver->waiting_to_hospital.begin(),
				// 	solver->waiting_to_hospital.end(), 0.0);
				num_scenario = solver->waiting_on_scene_penalized.size();
				num_total += num_scenario;
			}else{
				sum_total += solver->obj;
				++num_total;
			}
		}
	}

	double result_total = (num_total == 0) ? 0 : sum_total / num_total;
	return SecondStageWaitingTime{0, 0, result_total};
	
}

vector<vector<Call>> PolicyTester::get_future_scenarios(double time, vector<Call>& this_scenario,
	int nb_realizations, int T, vector<vector<vector<Call>>>& my_scenarios){
	int t_slot = 0;
	double t0 = 18*3600;
	double slot_duration = 1800;
	using namespace fmt;

	if(time <= t0 + slot_duration){
		t_slot = 0;
	}else if(time <= t0 + 2*slot_duration){
		t_slot = 1;
	}else if(time <= t0 + 3*slot_duration){
		t_slot = 2;
	}else{
		t_slot = 3;
	}

	int lb_scen = 0;
	while(lb_scen < this_scenario.size() && this_scenario[lb_scen].time <= time){
		++lb_scen;
	}
	int ub_scen = lb_scen;
	while(ub_scen < this_scenario.size()&& this_scenario[ub_scen].time <= 
		t0 + (t_slot+1)*slot_duration){
		++ub_scen;
	}

	int nb_scen = 0;
	if(time >= t0 + 3*slot_duration){ //Last time slot, no future scenarios.
		nb_scen = 1;
		vector<vector<Call>> future_scenarios(nb_scen, vector<Call>());
		int ind = 0;
		for(int i = lb_scen; i < ub_scen; ++i){
			future_scenarios[ind].push_back(this_scenario[i]);
		}
		return future_scenarios;
	}

	default_random_engine gen;
	uniform_int_distribution<int> rand_scen(0,nb_realizations-1);
	nb_realizations = 1;
	if(time < t0 + slot_duration){
		t_slot = 0;
		// nb_scen = 125;
		nb_scen = 100;
		vector<vector<Call>> future_scenarios(nb_scen, vector<Call>());
		int ind = 0;
		// for(int u = 0; u < nb_realizations; ++u){
		// 	for(int v = 0; v < nb_realizations; ++v){
		// 		for(int w = 0; w < nb_realizations; ++w){
		// 			for(int i = lb_scen; i < ub_scen; ++i){
		// 				future_scenarios[ind].push_back(this_scenario[i]);
		// 			}
		// 			future_scenarios[ind].insert(future_scenarios[ind].end(),
		// 			my_scenarios[1][u].begin(), my_scenarios[1][u].end());

		// 			future_scenarios[ind].insert(future_scenarios[ind].end(),
		// 			my_scenarios[2][v].begin(), my_scenarios[2][v].end());

		// 			future_scenarios[ind].insert(future_scenarios[ind].end(),
		// 			my_scenarios[3][v].begin(), my_scenarios[3][v].end());
		// 			++ind;
		// 		}
		// 	}
		// }

		for(int ind = 0; ind < nb_scen; ++ind){
			for(int i = lb_scen; i < ub_scen; ++i){
				future_scenarios[ind].push_back(this_scenario[i]);
			}
			
			for(int t = 1; t < T; ++t){
				int u = rand_scen(gen);
				future_scenarios[ind].insert(future_scenarios[ind].end(),
					my_scenarios[t][u].begin(), my_scenarios[t][u].end());
			}
		}

		for(int i = 0; i < nb_scen; ++i){
			future_scenarios[i] = future_scenarios[0];
			for(int j = 0; j < future_scenarios[i].size(); ++j){
				future_scenarios[i][j].id = j;
			}
		}
		return future_scenarios;	
	}else if(time < t0 + 2*slot_duration){
		t_slot = 1;
		nb_scen = 25;
		vector<vector<Call>> future_scenarios(nb_scen, vector<Call>());
		int ind = 0;
		for(int u = 0; u < nb_realizations; ++u){
			for(int v = 0; v < nb_realizations; ++v){
				for(int i = lb_scen; i < ub_scen; ++i){
					future_scenarios[ind].push_back(this_scenario[i]);
				}
				future_scenarios[ind].insert(future_scenarios[ind].end(),
				my_scenarios[2][u].begin(), my_scenarios[2][u].end());

				future_scenarios[ind].insert(future_scenarios[ind].end(),
				my_scenarios[3][v].begin(), my_scenarios[3][v].end());
				++ind;
			}
		}
		for(int i = 0; i < nb_scen; ++i){
			future_scenarios[i] = future_scenarios[0];
			for(int j = 0; j < future_scenarios[i].size(); ++j){
				future_scenarios[i][j].id = j;
			}
		}
		return future_scenarios;
	}else if(time < t0 + 3*slot_duration){
		t_slot = 2;
		nb_scen = 5;
		vector<vector<Call>> future_scenarios(nb_scen, vector<Call>());
		for(int u = 0; u < nb_realizations; ++u){
			for(int i = lb_scen; i < ub_scen; ++i){
				future_scenarios[u].push_back(this_scenario[i]);
			}
			future_scenarios[u].insert(future_scenarios[u].end(),
					my_scenarios[3][u].begin(), my_scenarios[3][u].end());
		}
		for(int i = 0; i < nb_scen; ++i){
			future_scenarios[i] = future_scenarios[0];
			for(int j = 0; j < future_scenarios[i].size(); ++j){
				future_scenarios[i][j].id = j;
			}
		}
		return future_scenarios;
	}else{
		return vector<vector<Call>>();
	}
}



SecondStageWaitingTime PolicyTester::get_waiting_time(int amb_id, Call& call, 
	const string& policy, vector<int>& index_scenarios, int nearest_base_id, 
	std::vector<int> queue, int sc, double time){
	vector<Ambulance> temp_ambulances = ambulances;
	double min_time = 0;
	double waiting_to_hospital = 0;
	if(amb_id >= 0){
		Ambulance& amb = temp_ambulances[amb_id];
		min_time = travel.get_response_time(amb, call, time, true);
		waiting_to_hospital = amb.answer_call(call, travel, ins, time, min_time, 
			nearest_base_id);	
	}

	double sum_total = 0;
	int num_total = 0;

	vector<int> busy_ambulances;
	if(amb_id == -1){
		for(auto& amb: temp_ambulances){
			if(amb.arrival_time_at_f_last_trip > time){
				busy_ambulances.push_back(amb.id);
			}
		}

		if(busy_ambulances.size() == 0){
			return SecondStageWaitingTime{0, 0, GRB_INFINITY};
		}
	}

	for(int j = 0; j < index_scenarios.size(); ++j){
		if(index_scenarios[j] < ins.calls[j].size() && 
			ins.calls[j][index_scenarios[j]].time >= call.time){
			vector<Ambulance> aux_ambulances = temp_ambulances;
			vector<Call> future_calls(ins.calls[j].begin() + index_scenarios[j],
				ins.calls[j].end());

			
			for(auto i: queue){
				if(amb_id == -1 ||  (i != call.id)){
					future_calls.insert(future_calls.begin(), ins.calls[sc][i]);
					if(i == call.id){
						// If the current call must go to queue, 
						// all free ambulances are excluded.
						for(auto& amb: aux_ambulances){
							if(amb.arrival_time_at_f_last_trip <= time){
								future_calls[0].ambulances.erase(
									remove(future_calls[0].ambulances.begin(),
										future_calls[0].ambulances.end(), amb_id),
										future_calls[0].ambulances.end());
							}
						}
					}
				}
			}

			for(int i = 0; i < future_calls.size(); ++i){
				future_calls[i].id = i;
			}
			
			auto solver = get_solver(policy, future_calls, aux_ambulances, travel);
			solver->run();

			// if(amb_id == -1){
			// 	solver->print_results();
			// 	cout << "curr call " << min_time*ins.penalties[call.priority] << "\n";
			// 	std::cin.get();
			// }


			if(policy != "cg"){
				sum_total += accumulate(solver->waiting_on_scene_penalized.begin(),
					solver->waiting_on_scene_penalized.end(), 0.0);
				// sum_total += accumulate(solver->waiting_to_hospital.begin(),
				// 	solver->waiting_to_hospital.end(), 0.0);
				num_total += solver->waiting_on_scene.size();
			}else{
				sum_total += solver->obj;
				++num_total;
			}
			// fmt::print("Sum total {}\n", sum_total);
			// std::cin.get();
		}
	}

	if(amb_id != -1){
		sum_total += min_time*ins.penalties[call.priority];
		num_total += 1;
	}

	double result_total = (num_total == 0) ? 0 : sum_total / num_total;
	return SecondStageWaitingTime{(time + min_time) - call.time,
		waiting_to_hospital, result_total};
}

SecondStageWaitingTime PolicyTester::get_waiting_time_return(int amb_id, const string& policy,
		vector<int>& index_scenarios, double time, std::vector<int> queue, int sc, int base){
	vector<Ambulance> temp_ambulances = ambulances;

	auto& amb = temp_ambulances[amb_id];
	amb.base_location = ins.bases[base];
	amb.arrival_time_at_b_last_trip = time + travel.travel_time(
		amb.free_location, amb.base_location);

	double min_time = 0;
	double waiting_to_hospital = 0;
	double sum_total = 0;
	int num_total = 0;
	for(int j = 0; j < index_scenarios.size(); ++j){
		if(index_scenarios[j] < ins.calls[j].size() && 
			ins.calls[j][index_scenarios[j]].time >= time){
			vector<Ambulance> aux_ambulances = temp_ambulances;
			vector<Call> future_calls(ins.calls[j].begin() + index_scenarios[j],
				ins.calls[j].end());

			for(auto i: queue){
				future_calls.insert(future_calls.begin(), ins.calls[sc][i]);
				// By returning the current ambulance to base,
				// we forbid it to answer any calls in queue:
				future_calls[0].ambulances.erase(remove(future_calls[0].ambulances.begin(),
					future_calls[0].ambulances.end(), amb_id),
					future_calls[0].ambulances.end());
			}

			for(int i = 0; i < future_calls.size(); ++i){
				future_calls[i].id = i;
			}
			
			auto solver = get_solver(policy, future_calls, aux_ambulances, travel);
			solver->run();
			// solver->print_results();
			// std::cin.get();

			if(policy != "cg"){
				sum_total += accumulate(solver->waiting_on_scene_penalized.begin(),
					solver->waiting_on_scene_penalized.end(), 0.0);
				sum_total += accumulate(solver->waiting_to_hospital.begin(),
					solver->waiting_to_hospital.end(), 0.0);
				num_total += solver->waiting_on_scene.size();
			}else{
				sum_total += solver->obj;
				++num_total;
			}
			// fmt::print("Sum total {}\n", sum_total);
			// std::cin.get();
		}
	}

	double result_total = (num_total == 0) ? 0 : sum_total / num_total;
	return SecondStageWaitingTime{(time + min_time) - time,
		waiting_to_hospital, result_total};
}


vector<int> PolicyTester::get_nearest_base(vector<Call>& calls){
	vector<int> nearest_base;
	for(int i = 0; i < calls.size(); ++i){
		auto& call = calls[i];
		double min_dist = GRB_INFINITY;
		int min_ind = -1;
		for(int b = 0; b < ins.nb_bases; ++b){
			double dist = GRB_INFINITY;
			auto& base = ins.bases[b];
			if(call.clean_needed){
				dist = travel.norm(ins.cleaning_bases[call.cleaning], base);
			}else if(call.hosp_needed){
				dist = travel.norm(ins.hospitals[call.hospital], base); 
			}else{
				dist = travel.norm(call.location, base);
			}

			if(dist < min_dist){
				min_dist = dist;
				min_ind = b;
			}
		}
		nearest_base.push_back(min_ind);
	}

	return nearest_base;
}


int PolicyTester::set_next_event(vector<Call>& this_scenario, int& event_call, 
	int& index_call){
	vector<pair<double, int>> future_arrival_times;
	for(int i = 0; i < ambulances.size(); ++i){
		auto& amb = ambulances[i];
		if(amb.arrival_time_at_f_last_trip > time + g_params.EPS){
			future_arrival_times.push_back(make_pair(amb.arrival_time_at_f_last_trip,
				i));
		}
	}

	int nb_calls = this_scenario.size();
	if(future_arrival_times.size() > 0){
		auto min_elem = *min_element(future_arrival_times.begin(),
			future_arrival_times.end());

		double min_arrival_time = min_elem.first;
		if(index_call < nb_calls-1 && this_scenario[index_call+1].time <= 
			min_arrival_time){
			event_call = 1;
			index_call += 1;
			time = this_scenario[index_call].time;
		}else{
			event_call = 0;
			time = min_arrival_time;
			return min_elem.second;
		}
	}else if(index_call < nb_calls - 1){
		event_call = 1;
		index_call += 1;
		time = this_scenario[index_call].time;
	}

	return -1;
}


Stats::Stats(vector<double> waiting_on_scene, vector<double> waiting_on_scene_penalized, 
	vector<double> waiting_to_hospital){
	assert(waiting_on_scene.size() == waiting_to_hospital.size() &&
		"Stats ERROR: Waiting times vectors must have the same size!");
	if(waiting_on_scene.size() == 0){
		mean_waiting_on_scene = GRB_INFINITY;
		mean_waiting_on_scene_penalized = GRB_INFINITY;
		mean_waiting_to_hospital = GRB_INFINITY;
		max_waiting_on_scene = GRB_INFINITY;
		max_waiting_on_scene_penalized = GRB_INFINITY;
		max_waiting_to_hospital = GRB_INFINITY;
		mean_total = GRB_INFINITY;
		max_total = GRB_INFINITY;
		waiting_on_scene_q90 = -1;
		waiting_to_hospital_q90 = -1;
		q90_total = -1;
	}else{
		mean_waiting_on_scene = accumulate(waiting_on_scene.begin(),
			waiting_on_scene.end(), 0.0)/waiting_on_scene.size(); 
		mean_waiting_on_scene_penalized = accumulate(waiting_on_scene_penalized.begin(),
			waiting_on_scene_penalized.end(), 0.0)/waiting_on_scene_penalized.size(); 
		mean_waiting_to_hospital = accumulate(waiting_to_hospital.begin(),
			waiting_to_hospital.end(), 0.0)/waiting_to_hospital.size();
		max_waiting_on_scene = *max_element(waiting_on_scene.begin(),
			waiting_on_scene.end());
		max_waiting_on_scene_penalized = *max_element(waiting_on_scene_penalized.begin(),
			waiting_on_scene_penalized.end());
		max_waiting_to_hospital = *max_element(waiting_to_hospital.begin(),
			waiting_to_hospital.end());
		vector<double> total = waiting_on_scene;
		transform(total.begin(), total.end(), waiting_to_hospital.begin(),
			total.begin(), plus<double>());
		mean_total = accumulate(total.begin(), total.end(), 0.0) / total.size();
		max_total = *max_element(total.begin(), total.end());
		sort(waiting_on_scene.begin(), waiting_on_scene.end());
		sort(waiting_to_hospital.begin(), waiting_to_hospital.end());
		sort(total.begin(), total.end());

		waiting_on_scene_q90 = waiting_on_scene[(int) floor(
			0.9*waiting_on_scene.size())];
		waiting_to_hospital_q90 = waiting_to_hospital[(int) floor(
			0.9*waiting_on_scene.size())];
		q90_total = total[(int) floor(0.9*waiting_on_scene.size())];
	}
}

bool PolicyTester::can_answer(Ambulance& amb, Call& call){
	bool amb_allowed = false;
	for(int amb_id: call.ambulances){
		if(amb_id == amb.id){
			amb_allowed = true;
			break;
		}
	}

	return amb_allowed && amb.type <= call.priority;
}


PolicyTester::~PolicyTester(){

}