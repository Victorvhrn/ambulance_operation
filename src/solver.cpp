#include "../include/travel.h"
#include "../include/instance.h"
#include "../include/call_model.h"
#include "../include/full_model_det.h"
#include "../include/data.h"
#include "../include/cg.h"
#include "../include/future.h"
#include "../include/solver.h"


Solver::Solver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
	Instance& ins, Travel& travel): env(env), calls(calls), ambulances(ambulances), 
	travel(travel), ins(ins), time(0), first_time(calls.front().time), 
	last_time(calls.back().time+7200), waiting_on_scene(calls.size(), GRB_INFINITY),
	waiting_on_scene_penalized(calls.size(),GRB_INFINITY),
	waiting_to_hospital(calls.size(), GRB_INFINITY),
	calls_end(calls.size(), GRB_INFINITY),
	which_ambulance(calls.size(), -1), obj(0){

	set_calls_nearest_bases();
}

void Solver::set_calls_nearest_bases(){
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
}


void Solver::set_next_event(int& event_call, int& index_call){
	vector<double> future_arrival_times;
	for(int i = 0; i < ambulances.size(); ++i){
		auto& amb = ambulances[i];
		if(amb.arrival_time_at_f_last_trip > time + g_params.EPS){
			future_arrival_times.push_back(amb.arrival_time_at_f_last_trip);
		}
	}

	int nb_calls = calls.size();
	if(future_arrival_times.size() > 0){
		double min_arrival_time = *min_element(future_arrival_times.begin(),
			future_arrival_times.end());

		if(index_call < nb_calls-1 && calls[index_call+1].time <= 
			min_arrival_time){
			event_call = 1;
			index_call += 1;
			time = calls[index_call].time;
		}else{
			event_call = 0;
			time = min_arrival_time;
		}
	}else if(index_call < nb_calls - 1){
		event_call = 1;
		index_call += 1;
		time = calls[index_call].time;
	}
}

bool Solver::can_answer(Ambulance& amb, Call& call){
	bool amb_allowed = false;
	for(int amb_id: call.ambulances){
		if(amb_id == amb.id){
			amb_allowed = true;
			break;
		}
	}

	return amb_allowed && amb.type <= call.priority;
}


void Solver::print_results(){
	double average_scene = 0;
	double average_hospital = 0;
	double average_pen = 0;
	cout << "Call\tWait On Scene(s)\tWait Penalized(s)\tAmb:\n";
	for(int i = 0; i < calls.size(); ++i){
		auto& call = calls[i];
		cout << call << "\t"<< waiting_on_scene[i] << "\t";
		cout << waiting_on_scene_penalized[i] << "\t" << which_ambulance[i] << "\n";

		average_scene += waiting_on_scene[i];
		average_pen += waiting_on_scene_penalized[i];
		average_hospital += waiting_to_hospital[i];
	}
	average_scene /= calls.size();
	average_hospital /= calls.size();
	average_pen /= calls.size();
	cout << "Average: Real " << average_scene << " Pen " << average_pen << "\n";
	cout << "Total: Real " << average_scene*calls.size() << "  Pen "; 
	cout << average_pen*calls.size() << "\n";
}

Solver::~Solver(){

}

QueueSolver::QueueSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
	Instance& ins, Travel& travel): Solver(env, calls,ambulances,ins,travel){
	travel.set_forward(false);
	queue.reserve(calls.size());
}


void QueueSolver::run(){
	int event_call = 1;
	time = calls[0].time;
	int index_call = 0;
	int calls_attended = 0;
	travel.set_forward(false);
	using fmt::print;
	while(calls_attended < calls.size()){

		if(event_call == 1){
			queue.push_back(index_call);
		}
		std::vector<int> queue_aux;
		for(int i = 0; i < queue.size(); ++i){
			auto t0 = std::chrono::high_resolution_clock::now();
			auto& call = calls[queue[i]];
			double min_time = GRB_INFINITY;
			int index_amb = -1;
			// print("Call {} {:.3f} ({:.3f},{:.3f}) {}:\n", call.id, call.time, 
			// 	call.location.first, call.location.second, call.priority);
			for(auto& amb: ambulances){
				if(can_answer(amb,call)){
					double time_amb = travel.get_response_time(amb,call,time);
					// print("\tAmb {} {} {}\n", amb.id, time_amb, amb.type);
					if(time_amb < min_time){
						min_time = time_amb;
						index_amb = amb.id;
					}
				}
			}
			if(index_amb >= 0){
				// auto& amb = ambulances[index_amb];
				// if(amb.arrival_time_at_b_last_trip <= time){
				// 	fmt::print("Ambulance {} at base {}\n",index_amb,amb.base_location);
				// }else{
				// 	Location current_location = travel.ambulance_position(amb,time);
				// 	fmt::print("Ambulance {} returning {}\n",index_amb, current_location);
				// }

				// cin.get();
				double waiting_on_scene_i = time + min_time - call.time;
				double waiting_to_hospital_i = ambulances[index_amb].answer_call(call, 
					travel, ins, time, min_time, nearest_base[queue[i]]);
				waiting_on_scene[queue[i]] = waiting_on_scene_i;
				waiting_on_scene_penalized[queue[i]] = waiting_on_scene_i*ins.penalties[call.priority];
				waiting_to_hospital[queue[i]] = waiting_to_hospital_i;
				which_ambulance[queue[i]] = index_amb;
				calls_end[queue[i]] = call.end;
				calls_attended++;
			}else{
				queue_aux.push_back(queue[i]);
			}
			auto dt = std::chrono::high_resolution_clock::now();
			run_times.push_back(std::chrono::duration_cast<chrono::microseconds>(dt - 
				t0).count());
		}
		queue.clear();
		for(auto i: queue_aux){
			queue.push_back(i);
		}
		set_next_event(event_call, index_call);
	}	
}


ForwardSolver::ForwardSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
	Instance& ins, Travel& travel): Solver(env, calls,ambulances, ins,travel){
	travel.set_forward(true);
}

void ForwardSolver::run(){
	travel.set_forward(true);
	for(int i = 0; i < calls.size(); ++i){
		auto t0 = std::chrono::high_resolution_clock::now();
		double min_time = GRB_INFINITY;
		int index_amb = -1;
		auto& call = calls[i];
		time = call.time;
		for(auto& amb: ambulances){
			if(can_answer(amb,call)){
				double time_amb = travel.get_response_time(amb, call, time);
				if(amb.id == 26){
					std::cout << amb << " " << time_amb << "\n";
					std::cout << "Dist = " << travel.lat_long_travel_distance(
						amb.base_location, call.location) << "\n";
					fmt::print("base {}\n", amb.base_location);
					fmt::print("call {}\n", call.location);
					std::cin.get();
				}
				if(time_amb < min_time){
					min_time = time_amb;
					index_amb = amb.id;
				}else if(abs(time_amb - min_time) <  g_params.EPS && index_amb != -1 &&
					amb.type > ambulances[index_amb].type){
					min_time = time_amb;
					index_amb = amb.id;
				}
			}
		}
		
		if(index_amb >= 0){
			double waiting_on_scene_i = time + min_time - call.time;
			double waiting_to_hospital_i = ambulances[index_amb].answer_call(call, travel, 
				ins, time, min_time, nearest_base[i]);

			waiting_on_scene[i] = waiting_on_scene_i;
			waiting_on_scene_penalized[i] = waiting_on_scene_i*ins.penalties[call.priority];
			waiting_to_hospital[i] = waiting_to_hospital_i;
			which_ambulance[i] = index_amb;
			calls_end[i] = call.end;
		}else{
			waiting_on_scene[i] = GRB_INFINITY;
			waiting_on_scene_penalized[i] = GRB_INFINITY;
			waiting_to_hospital[i] = GRB_INFINITY;
			calls_end[i] = GRB_INFINITY;
			which_ambulance[i] = -1;
		}
		auto dt = std::chrono::high_resolution_clock::now();
		run_times.push_back(std::chrono::duration_cast<chrono::microseconds>(dt - 
			t0).count());
	}
}



CGSolver::CGSolver(GRBEnv& env, vector<Call>& calls, 
	vector<Ambulance>& ambulances, Instance& ins, Travel& travel): 
	Solver(env, calls,ambulances,ins,travel){
	g_params.extended_model = false;
	travel.set_forward(false);
	queue.reserve(calls.size());
}


void CGSolver::run(){
	if(calls.size() == 0){
		fmt::print("No calls\n");
		return;
	}
	time = calls[0].time;
	Data data(this, &calls[0], NULL);
	CallModel1 cg(data,env,*this);
	cg.solve();
	obj = cg.model.get(GRB_DoubleAttr_ObjVal);
}


int CGSolver::get_return_base(CGCall& cg, Data& data, int t, int c, int a, int l1,
	int l, int h){
	int end_t = data.get_time_slot(time + t*data.quantum + data.tao[t][c][a][l1][l][h]);
	fmt::print("End t = {} {}\n", end_t, time + t*data.quantum + 
		data.tao[t][c][a][l1][l][h]);
	fmt::print("quantum = {}\n", data.quantum);
	for(int b = 0; b < data.num_bases; ++b){
		auto name = cg.yt_ahb[end_t][a][h][b].get(GRB_StringAttr_VarName);
		auto val = cg.yt_ahb[end_t][a][h][b].get(GRB_DoubleAttr_X);
		fmt::print("{} = {}\n", name, val);
		if(val > 0.5){
			return b;
		}
	}
	
	for(int c1 = 0; c1 < data.types_call; ++c1){
		if(a <= c1){
			for(int l2 = 0; l2 < data.num_locals; ++l2){
				for(auto h1: data.H[c1][l2]){
					auto name = cg.xt_cahlh[end_t][c1][a][h][l2][h1].get(
						GRB_StringAttr_VarName);
					auto val = cg.xt_cahlh[end_t][c1][a][h][l2][h1].get(GRB_DoubleAttr_X);
					if(val > g_params.EPS){
						fmt::print("{} {}\n",name, val);
						return ins.nearest_base_to_region[l2];
					}
				}
			}
		}
	}
	int ind_b = -1;
	double min_d = GRB_INFINITY;
	for(int b = 0; b < data.num_bases; ++b){
		if(travel.norm(ins.hospitals[h], ins.bases[b]) < min_d){
			double d = travel.norm(ins.hospitals[h], ins.bases[b]);
			if(d < min_d){
				min_d = d;
				ind_b = b;
			}
			
		}
	}
	
	return ind_b;
}

ModelSolver::ModelSolver(GRBEnv& env, vector<Call>& calls, 
	vector<Ambulance>& ambulances, Instance& ins, Travel& travel): 
	Solver(env, calls,ambulances,ins,travel){
	g_params.extended_model = false;
	travel.set_forward(false);
	queue.reserve(calls.size());
}

void ModelSolver::run(){
	int event_call = 1;
	time = calls[0].time;
	int index_call = 0;
	int calls_attended = 0;

	using fmt::print;
	std::vector<std::set<int>> H(calls.size(), std::set<int>());

	for(int i = 0; i < calls.size(); ++i){
		std::vector<std::pair<double, int>> sorted_hospitals;
		sorted_hospitals.reserve(ins.nb_hospitals);
		for(int h = 0; h < ins.nb_hospitals; ++h){
			double d = travel.norm(calls[i].location, ins.hospitals[h]);
			sorted_hospitals.push_back(std::make_pair(d,h));
		}
		std::sort(sorted_hospitals.begin(), sorted_hospitals.end());
		for(int k = 0; k < g_params.n_nearest_hospitals && k < ins.nb_hospitals; ++k){
			H[i].insert(sorted_hospitals[k].second);
		}
	}
	while(calls_attended < calls.size()){
		//Find future calls in each scenario.
		vector<int> index_scenarios(ins.calls.size(), 0);
		for(int j = 0; j < ins.calls.size(); ++j){
			while(index_scenarios[j] < ins.calls[j].size() && 
				ins.calls[j][index_scenarios[j]].time <= time){
				++index_scenarios[j];
			}
		}
		// fmt::print("Time {}, Queue Size {}\n",time, queue.size());
		if(event_call == 1){
			auto& call = calls[index_call];
			Data data(this, &call, NULL);
			double min_avg_cost = GRB_INFINITY;
			int best_amb = -1;
			int best_h = -1;
			//for each ambulance capable of answering call
			std::cout << "Call " << call << ":\n";
			for(int a = 0; a < ambulances.size(); ++a){
				auto& amb = ambulances[a];
				if(amb.arrival_time_at_f_last_trip <= time && can_answer(amb,call)){
					int l0 = -1, b0 = -1;

					double min_d = GRB_INFINITY;
					for(int b = 0; b < ins.nb_bases; ++b){
						double d = travel.norm(amb.base_location, ins.bases[b]);
						if(d < min_d){
							min_d = d;
							b0 = b;
						}
					}

					if(amb.arrival_time_at_b_last_trip > time){
						auto current_location = travel.ambulance_position(amb,time);
						l0 = data.get_location_index(current_location);	
					}

					for(auto h: H[call.id]){
						//Average over all scenarios
						double sum_cost = 0;
						for(int j = 0; j < ins.calls.size(); ++j){
							std::vector<Call> future_calls(ins.calls[j].begin() + 
								index_scenarios[j], ins.calls[j].end());
							FutureCall fc(data,env,*this, amb.type, l0,b0,h,future_calls);
							fc.solve();
							sum_cost += fc.obj;
						}

						double avg_cost = sum_cost / ins.calls.size();
						fmt::print("a{} l{} b{} h{} => {}\n",amb.type, l0,b0,h, avg_cost);

						if(avg_cost < min_avg_cost){
							min_avg_cost = avg_cost;
							best_amb = a;
							best_h = h;
						}
					}
				}
			}

			// Average cost of not answering over all scenarios
			double sum_cost = 0;
			for(int j = 0; j < ins.calls.size(); ++j){
				std::vector<Call> future_calls(ins.calls[j].begin() + 
					index_scenarios[j], ins.calls[j].end());
				FutureCall fc(data,env,*this, -1, -1,-1,-1,future_calls);
				fc.solve();
				sum_cost += fc.obj;
			}

			double avg_cost = sum_cost / ins.calls.size();
			fmt::print("queue => {}\n",avg_cost);
			if(avg_cost < min_avg_cost){
				min_avg_cost = avg_cost;
				best_amb = best_h = -1;
			}

			//Answer call with best average cost
			if(best_amb == -1){
				fmt::print("queue\n");
				queue.push_back(index_call);
			}else{
				Ambulance& amb = ambulances[best_amb];
				std::cout << "Attended by " << amb << "\n";
				call.hospital = best_h;
				double min_time = travel.get_response_time(amb,call,time);
				double waiting_on_scene_i = time + min_time - call.time;
				double waiting_to_hospital_i = amb.answer_call(call, travel, ins, time,
					min_time, 0); //all returns to base 0 temporarily
				waiting_on_scene[call.id] = waiting_on_scene_i;
				waiting_on_scene_penalized[call.id] = waiting_on_scene_i*ins.penalties[call.priority];
				waiting_to_hospital[call.id] = waiting_to_hospital_i;
				which_ambulance[call.id] = best_amb;
				calls_end[call.id] = call.end;
				calls_attended++;
			}

			fmt::print("==========================\n");
		}else{
			Ambulance* ambulance = NULL;
			for(int a = 0; a < ambulances.size(); ++a){
				if(ambulances[a].arrival_time_at_f_last_trip == time){
					ambulance = &ambulances[a];
				}
			}
			std::cout << "Ambulance " << *ambulance << ":\n";
			Data data(this, NULL, ambulance);
			// std::cout << "Returning " << *ambulance << "\n";
			//for each base
			double min_avg_cost = GRB_INFINITY;
			int best_b = -1;
			for(int b = 0; b < ins.nb_bases; ++b){
				int sum_b = 0;
				for(int a = 0; a < data.types_amb; ++a){
					sum_b += data.A0_ab[a][b];
				}

				if(sum_b < data.cap_bases[b]){
					double sum_cost = 0;
					for(int j = 0; j < ins.calls.size(); ++j){
						std::vector<Call> future_calls(ins.calls[j].begin() + 
							index_scenarios[j], ins.calls[j].end());
						FutureAmbulance fc(data,env,*this, NULL, -1, b, future_calls);
						fc.solve();
						// if(j == 1 && index_call == calls.size() - 1 && time > 8410){
						// 	GRBVar* vars = fc.model.getVars();
						// 	for(int i = 0; i < fc.model.get(GRB_IntAttr_NumVars); ++i){
						// 		auto name = vars[i].get(GRB_StringAttr_VarName); 
						// 		auto val = vars[i].get(GRB_DoubleAttr_X);
						// 		if(val > g_params.EPS){
						// 			fmt::print("{} = {}\n",name, val);
						// 		}
						// 	}
						// 	delete[] vars;

						// 	// fc.model.write(fmt::format("model_b{}.lp", b));
						// 	std::cin.get();
						// }
						sum_cost += fc.obj;
					}

					//average cost over all scenarios
					double avg_cost = sum_cost / ins.calls.size();
					fmt::print("Return to b => {}\n", avg_cost);
					if(avg_cost < min_avg_cost){
						min_avg_cost = avg_cost;
						best_b = b;
					}
				}
			}
			int call_ind = -1;
			int best_h = -1;
			for(auto i: queue){
				if(can_answer(*ambulance,calls[i])){
					for(auto h: H[calls[i].id]){
						double sum_cost = 0;
						for(int j = 0; j < ins.calls.size(); ++j){
							std::vector<Call> future_calls(ins.calls[j].begin() + 
								index_scenarios[j], ins.calls[j].end());
							FutureAmbulance fc(data,env,*this, &calls[i], h, -1, future_calls);
							fc.solve();
							// if(j == 1 && index_call == calls.size() - 1 && time > 8410){
							// 	GRBVar* vars = fc.model.getVars();
							// 	for(int i = 0; i < fc.model.get(GRB_IntAttr_NumVars); ++i){
							// 		auto name = vars[i].get(GRB_StringAttr_VarName); 
							// 		auto val = vars[i].get(GRB_DoubleAttr_X);
							// 		if(val > g_params.EPS){
							// 			fmt::print("{} = {}\n",name, val);
							// 		}
							// 	}
							// 	delete[] vars;

							// 	// fc.model.write(fmt::format("model_b{}.lp", b));
							// 	std::cin.get();
							// }
							sum_cost += fc.obj;
						}

						double avg_cost = sum_cost / ins.calls.size();
						fmt::print("Answer call {} => {}\n", i, avg_cost);
						if(avg_cost < min_avg_cost){
							min_avg_cost = avg_cost;
							best_b =  -1;
							call_ind = i;
							best_h = h;
						}
					}
				}
			}

			if(best_b == -1){
				//assign ambulance to call_ind and best_h
				auto& call = calls[call_ind];
				std::cout << "Attending " << call << "\n";
				double min_time = travel.get_response_time(*ambulance, call, time);
				double waiting_on_scene_i = time + min_time - call.time;
				double waiting_to_hospital_i = ambulance->answer_call(call, travel, ins, time,
					min_time, 0); //returns to base 0 temporarily

				waiting_on_scene[call.id] = waiting_on_scene_i;
				waiting_on_scene_penalized[call.id] = waiting_on_scene_i*ins.penalties[call.priority];
				waiting_to_hospital[call.id] = waiting_to_hospital_i;
				which_ambulance[call.id] = ambulance->id;
				calls_end[call.id] = call.end;
				//remove call_ind from queue
				queue.erase(std::remove(queue.begin(), queue.end(), call_ind), 
					queue.end());
				calls_attended++;
			}else{
				fmt::print("returning to {}\n", best_b);
				ambulance->base_location = ins.bases[best_b];
				ambulance->arrival_time_at_b_last_trip = time + travel.travel_time(
					ambulance->free_location, ambulance->base_location);
			}
		}
		fmt::print("==========================\n");
		set_next_event(event_call, index_call);
	}	
}

