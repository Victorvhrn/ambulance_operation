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

NonMiopycSolver::NonMiopycSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel): Solver(env, calls,ambulances, ins,travel){
	travel.set_forward(true);
}


void NonMiopycSolver::run(){
	std::vector<bool> is_allocated(calls.size(), false);
	std::vector<std::vector<double>> travel_times(calls.size(), 
			std::vector<double>(ambulances.size(), GRB_INFINITY));
	int max_index = -1;
	std::vector<double> min_times(calls.size(), GRB_INFINITY);
	std::vector<int> index_amb(calls.size(), -1);
	std::vector<int> amb_type(calls.size(), -1);
	travel.set_forward(true);

	int i = 0;

	while(i < calls.size()){
		if(!is_allocated[i]){
			auto t0 = std::chrono::high_resolution_clock::now();
			if(i == max_index + 1){
				for(int j = 0; j < ambulances.size(); ++j){
					auto & amb = ambulances[j];
					if(can_answer(amb, calls[i])){
						travel_times[i][j] = travel.get_response_time(amb,calls[i],
							calls[i].time);
					}
				}

				compute_min_times(i, min_times, index_amb, amb_type);
				max_index = i;
			}

			while(!is_allocated[i]){
				auto& amb = ambulances[index_amb[i]];

				if(amb.arrival_time_at_f_last_trip <= calls[i].time){
					is_allocated[i] = true;
					double waiting_on_scene_i = min_times[i];
					double waiting_to_hospital_i = amb.answer_call(calls[i],travel,ins,calls[i].time,
						travel_times[i][amb.id], nearest_base[i]);

					waiting_on_scene[i] = waiting_on_scene_i;
					waiting_on_scene_penalized[i] = waiting_on_scene_i*
						ins.penalties[calls[i].priority];
					waiting_to_hospital[i] = waiting_to_hospital_i;
					which_ambulance[i] = amb.id;
					calls_end[i] = calls[i].end;

					for(int k = i+1; k <= max_index; ++k){
						if(!is_allocated[k]){
							if(can_answer(amb, calls[k])){
								travel_times[k][index_amb[i]] = travel.get_response_time(
									amb,calls[k],calls[k].time);

								compute_min_times(k,min_times, index_amb, amb_type);
							}
						}
					}
				}else{
					int k = max_index + 1;
					bool cont = true;
					while(cont){
						if(k >= calls.size()){
							cont = false;
						}else{
							if(calls[k].time > amb.arrival_time_at_f_last_trip){
								cont = false;
							}else{
								if(!is_allocated[k]){
									for(int j = 0; j < ambulances.size(); ++j){
										auto& amb_j = ambulances[j];
										if(can_answer(amb_j,calls[k])){
											travel_times[k][j] = travel.get_response_time(amb_j,
												calls[k], calls[k].time);
										}
									}
									compute_min_times(k,min_times,index_amb,amb_type);
								}
								k += 1;
							}
						}
					}

					max_index = k - 1;
					double maxmin = -1;
					int which_k = -1;
					double t1 = ins.penalties[calls[i].priority]*min_times[i];
					k = i + 1;
					if(i < calls.size() - 1){
						auto& amb_i = ambulances[index_amb[i]];
						bool contd = calls[k].time <= amb_i.arrival_time_at_f_last_trip;

						while(contd){
							if(!is_allocated[k]){
								double t2 = ins.penalties[calls[k].priority]*min_times[k];

								if(abs(min_times[k] - travel_times[k][index_amb[i]]) < g_params.EPS && 
									t2 > maxmin){
									which_k = k;
                                	maxmin = t2;
								}
							}

							if(k+1 >= calls.size()){
								contd = false;
							}else if(calls[k+1].time > amb_i.arrival_time_at_f_last_trip){
								contd = false;
							}else{
								k += 1;
							}
						}
					}

					if(maxmin > t1){
						is_allocated[which_k] = true;
						auto& amb_i = ambulances[index_amb[i]];
						double waiting_on_scene_i = min_times[which_k];
						double waiting_to_hospital_i = amb_i.answer_call(calls[which_k],
							travel,ins, calls[which_k].time, travel_times[which_k][amb_i.id], 
							nearest_base[which_k]);
						waiting_on_scene[which_k] = waiting_on_scene_i;
						waiting_on_scene_penalized[which_k] = waiting_on_scene_i*
							ins.penalties[calls[which_k].priority];
						waiting_to_hospital[which_k] = waiting_to_hospital_i;
						which_ambulance[which_k] = amb_i.id;
						calls_end[which_k] = calls[which_k].end;

						for(int k = i; k <= max_index; ++k){
							if(!is_allocated[k]){
								if(can_answer(amb_i,calls[k])){
									travel_times[k][index_amb[i]] = travel.get_response_time(
										amb_i, calls[k], calls[k].time);
								}
								compute_min_times(k,min_times,index_amb,amb_type);
							}
						}
					}else{
						is_allocated[i] = true;
						auto& amb_i = ambulances[index_amb[i]];
						double waiting_on_scene_i = min_times[i];
						double waiting_to_hospital_i = amb_i.answer_call(calls[i],travel,ins,calls[i].time,
							travel_times[i][index_amb[i]], nearest_base[i]);

						waiting_on_scene[i] = waiting_on_scene_i;
						waiting_on_scene_penalized[i] = waiting_on_scene_i*
							ins.penalties[calls[i].priority];
						waiting_to_hospital[i] = waiting_to_hospital_i;
						which_ambulance[i] = amb_i.id;
						calls_end[i] = calls[i].end;

						for(int k = i; k <= max_index; ++k){
							if(!is_allocated[k]){
								if(can_answer(amb_i,calls[k])){
									travel_times[k][index_amb[i]] = travel.get_response_time(
										amb_i,calls[k],calls[k].time);
								}

								compute_min_times(k,min_times,index_amb,amb_type);
							}
						}
					}
				}
			}
			auto dt = std::chrono::high_resolution_clock::now();
			run_times.push_back(std::chrono::duration_cast<chrono::microseconds>(dt - 
				t0).count());
		}else{
			++i;
		}
	}
}

void NonMiopycSolver::compute_min_times(int k, vector<double>& min_times, 
	vector<int>& index_amb, vector<int> & amb_type){
	double min_time = GRB_INFINITY;
	int best_amb = -1;
	int best_type = -1;

	for(int j = 0; j < ambulances.size(); ++j){
		auto& amb = ambulances[j];
		double time_amb = GRB_INFINITY;
		if(can_answer(amb, calls[k])){
			if(amb.arrival_time_at_b_last_trip <= calls[k].time){
				time_amb = travel.travel_time(amb.base_location, calls[k].location);
			}else if(amb.arrival_time_at_f_last_trip <= calls[k].time){
				Location current_location = travel.ambulance_position(amb,calls[k].time);
				time_amb = travel.travel_time(current_location,calls[k].location);
			}else{
				double time_free = amb.arrival_time_at_f_last_trip-calls[k].time;
				double time_free_to_call = travel.travel_time(amb.free_location,
					calls[k].location);
				time_amb = time_free + time_free_to_call;
			}
			

			if((time_amb < min_time) || (abs(time_amb - min_time) < g_params.EPS &&
				(amb.type > best_type))) {
				best_amb = j;
				best_type = amb.type;
				min_time = time_amb;
			}
		}
	}

	min_times[k] = min_time;
	index_amb[k] = best_amb;
	amb_type[k]  = best_type;
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
				// std::cout << amb << " " << time_amb << "\n";
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
			// print("Call {} {} ({},{}) {}:\n", call.id, call.time, call.location.first,
			// 	call.location.second, call.priority);
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


PrioritySolver::PrioritySolver(GRBEnv& env, vector<Call>& calls, 
	vector<Ambulance>& ambulances, Instance& ins, Travel& travel): 
	Solver(env, calls, ambulances, ins, travel){
	// travel.set_forward(false);
	queue.reserve(calls.size());
}

void PrioritySolver::run(){
	int event_call = 1;
	time = calls[0].time;
	int index_call = 0;
	int calls_attended = 0;

	// travel.set_forward(false);

	while(calls_attended < calls.size()){
		if(event_call == 1){
			queue.push_back(index_call);
		}

		std::vector<std::pair<double, int>> sorted_queue;
		if(queue.size() > 1){
			sorted_queue.reserve(queue.size());
			for(int i = 0; i < queue.size(); ++i){
				auto& call = calls[queue[i]];
				sorted_queue.push_back(std::make_pair(ins.penalties[call.priority]*
					(time - call.time), queue[i]));
			}
			std::sort(sorted_queue.begin(), sorted_queue.end(),
				std::greater<std::pair<double, int>>());
		}else if(queue.size() == 1){
			auto& call = calls[queue[0]];
			sorted_queue.push_back(std::make_pair(ins.penalties[call.priority]*
				(time - call.time), queue[0]));
		}

		std::vector<int> queue_aux;
		for(int k = 0; k < sorted_queue.size(); ++k){
			auto t0 = std::chrono::high_resolution_clock::now();
			auto call_ind = sorted_queue[k].second;
			auto& call = calls[call_ind];
			double min_time = GRB_INFINITY;
			int index_amb = -1;

			for(auto& amb: ambulances){
				if(can_answer(amb,call)){
					double time_amb = travel.get_response_time(amb,call,time);

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

			if(index_amb >= 0 && ambulances[index_amb].arrival_time_at_f_last_trip <= time){
				double waiting_on_scene_i = time + min_time - call.time;
				double waiting_to_hospital_i = ambulances[index_amb].answer_call(call, 
					travel, ins, time, min_time, nearest_base[call_ind]);
				waiting_on_scene[call_ind] = waiting_on_scene_i;
				waiting_on_scene_penalized[call_ind] = waiting_on_scene_i*ins.penalties[call.priority];
				waiting_to_hospital[call_ind] = waiting_to_hospital_i;
				which_ambulance[call_ind] = index_amb;
				calls_end[call_ind] = call.end;
				calls_attended++;
			}else{
				queue_aux.push_back(sorted_queue[k].second);
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

MinMaxPSolver::MinMaxPSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel): Solver(env, calls,ambulances,ins,travel){
	// travel.set_forward(false);
	queue.reserve(calls.size());
}


void MinMaxPSolver::run(){
	int event_call = 1;
	time = calls[0].time;
	int index_call = 0;
	int calls_attended = 0;
	// travel.set_forward(false);
	while(calls_attended < calls.size()){
		if(event_call == 1){
			queue.push_back(index_call);
		}

		std::vector<double> min_times(queue.size(), GRB_INFINITY);
		std::vector<double> min_times_p(queue.size(), GRB_INFINITY);
		std::vector<std::vector<double>> travel_times(queue.size(), 
			std::vector<double>(ambulances.size(), GRB_INFINITY));
		std::vector<int> index_ambs(queue.size(), -1);

		// For every call in queue, determine the best ambulance and the best
		// response time
		for(int i = 0; i < queue.size(); ++i){
			auto& call = calls[queue[i]];
			for(auto& amb: ambulances){
				int j = amb.id;
				if(can_answer(amb,call)){
					travel_times[i][j] = travel.get_response_time(amb, call, time);

					if(travel_times[i][j] < min_times[i]){
						min_times[i] = travel_times[i][j];
						index_ambs[i] = j;
					}else if(abs(travel_times[i][j] - min_times[i]) <  g_params.EPS && 
						index_ambs[i] != -1 && amb.type > ambulances[index_ambs[i]].type){
						min_times[i] = travel_times[i][j];
						index_ambs[i] = j;
					}
				}
			}
			min_times[i] += time - call.time;
			min_times_p[i] = min_times[i]*ins.penalties[call.priority];
		}


		std::vector<int> remaining_indexes(queue.size(),0);
		for(int i = 0; i < queue.size(); ++i){
			remaining_indexes[i] = i;
		}

		std::vector<int> queue_aux;
		int nb_treated = 0;
		int total_queue = queue.size();

		while(nb_treated < total_queue){
			auto t0 = std::chrono::high_resolution_clock::now();
			//get time and index of the worst call, and the best ambulance to such call
			auto max_it = std::max_element(min_times_p.begin(), min_times_p.end());
			int current_call = std::distance(min_times_p.begin(), max_it);
			double max_time = min_times[current_call];
			int call_ind = queue[remaining_indexes[current_call]];
			auto& call = calls[call_ind];
			int index_amb = index_ambs[remaining_indexes[current_call]];
			auto& best_amb = ambulances[index_amb];

			//if best ambulance is free, it should answer the call
			if(index_amb >= 0 && best_amb.arrival_time_at_f_last_trip <= time){
				double waiting_on_scene_i = max_time;
				double waiting_to_hospital_i = ambulances[index_amb].answer_call(call, 
					travel, ins, time, max_time - (time - call.time), nearest_base[call_ind]);
				waiting_on_scene[call_ind] = waiting_on_scene_i;
				waiting_on_scene_penalized[call_ind] = waiting_on_scene_i*ins.penalties[call.priority];
				waiting_to_hospital[call_ind] = waiting_to_hospital_i;
				which_ambulance[call_ind] = index_amb;
				calls_end[call_ind] = call.end;
				
				calls_attended++;
				//update travel times
				for(int k = 0; k < total_queue-nb_treated; ++k){
					int ind = remaining_indexes[k];
					if(k != current_call && index_amb == index_ambs[ind]){
						double time_to_h = best_amb.arrival_time_at_f_last_trip -
							time;
						double time_from_h_to_c = travel.travel_time(best_amb.free_location,
							calls[queue[ind]].location);
						travel_times[ind][best_amb.id] = time_to_h + time_from_h_to_c;
						auto min_it = std::min_element(travel_times[ind].begin(),
							travel_times[ind].end());
						double min_val = *min_it;
						int min_ind = std::distance(travel_times[ind].begin(), min_it);
						min_times[k] = min_val + time - calls[queue[ind]].time;
						min_times_p[k] = min_times[k] *
							ins.penalties[calls[queue[ind]].priority];
						index_ambs[ind] = min_ind;
					}
				}
			}else{
				queue_aux.push_back(queue[remaining_indexes[current_call]]);
			}
			if(nb_treated+1 < total_queue){
				min_times_p.erase(max_it);
				remaining_indexes.erase(std::remove(remaining_indexes.begin(), 
					remaining_indexes.end(), remaining_indexes[current_call]), 
					remaining_indexes.end());
				min_times.erase(std::remove(min_times.begin(), min_times.end(),
					min_times[current_call]), min_times.end());

			}
			nb_treated++;
			auto dt = std::chrono::high_resolution_clock::now();
			run_times.push_back(std::chrono::duration_cast<chrono::microseconds>(dt - 
				t0).count()*queue.size());
		}

		queue.clear();
		for(auto i: queue_aux){
			queue.push_back(i);
		}
	
		set_next_event(event_call, index_call);
	}
}

GenForwardSolver::GenForwardSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
	Instance& ins, Travel& travel): Solver(env, calls,ambulances,ins,travel){
	g_params.extended_model = false;
	travel.set_forward(true);
	queue.reserve(calls.size());
}

void GenForwardSolver::run(){
	int index_call = 0;
	queue.clear();
	while(index_call < calls.size()){
		fmt::print("index_call {}\n", index_call);
		time = calls[index_call].time;
		queue.push_back(index_call);
		int current = index_call+1;
		while(current < calls.size() && calls[current].time-time < g_params.EPS){
			queue.push_back(current);
			current += 1;
		}
		fmt::print("queue size {}\n",queue.size());
		if(queue.size() > 1){
			//call model
			//call answer_call to each call in queue
			//set index_call += |answered calls|
			try{
				CallModel cm(env,*this);
				cm.solve();
				for(int i = 0; i < queue.size(); ++i){
					auto& best_amb = ambulances[cm.opt_which_amb[i]];
					double min_time = cm.arrival_times[i] - time;
					double waiting_to_hospital_i = best_amb.answer_call(calls[queue[i]], 
						travel, ins, time, min_time, nearest_base[queue[i]]);
					auto& call = calls[queue[i]];
					waiting_on_scene[queue[i]] = cm.arrival_times[i] - calls[queue[i]].time;
					waiting_on_scene_penalized[queue[i]] = waiting_on_scene[queue[i]]*ins.penalties[call.priority];
					waiting_to_hospital[queue[i]] = waiting_to_hospital_i;
					which_ambulance[queue[i]] = cm.opt_which_amb[i];
					calls_end[queue[i]] = calls[queue[i]].end;	
				}
				index_call = current;
				// print_results();
			}catch(GRBException& ex){
				std::cout << ex.getMessage() << "\n";
				exit(1);
			}
		}else{
			double min_time = GRB_INFINITY;
			int index_amb = -1;
			auto& call = calls[queue[0]];
			time = call.time;
			for(auto& amb: ambulances){
				if(can_answer(amb,call)){
					double time_amb = travel.get_response_time(amb, call, time);
					if(time_amb < min_time){
						min_time = time_amb;
						index_amb = amb.id;
					}
				}
			}
			if(index_amb >= 0){
				double waiting_to_hospital_i = ambulances[index_amb].answer_call(call, 
					travel, ins, time, min_time, nearest_base[queue[0]]);
				waiting_on_scene[queue[0]] = time - call.time + min_time;
				waiting_on_scene_penalized[queue[0]] = waiting_on_scene[queue[0]]*ins.penalties[call.priority];
				waiting_to_hospital[queue[0]] = waiting_to_hospital_i;
				which_ambulance[queue[0]] = index_amb;
				calls_end[queue[0]] = calls[queue[0]].end;
				index_call += 1;
			}else{
				waiting_on_scene[queue[0]] = GRB_INFINITY;
				waiting_on_scene_penalized[queue[0]] = GRB_INFINITY;
				waiting_to_hospital[queue[0]] = GRB_INFINITY;
				which_ambulance[queue[0]] = -1;
				calls_end[queue[0]] = GRB_INFINITY;
			}
		}
		queue.clear();
		fmt::print("==================\n");
	}
}

GenMinMaxSolver::GenMinMaxSolver(GRBEnv& env, vector<Call>& calls, 
	vector<Ambulance>& ambulances, Instance& ins, Travel& travel): 
	Solver(env, calls,ambulances,ins,travel){
	g_params.extended_model = false;
	travel.set_forward(false);
	queue.reserve(calls.size());
}

void GenMinMaxSolver::run(){
	int event_call = 1;
	time = calls[0].time;
	int index_call = 0;
	int calls_attended = 0;

	while(calls_attended < calls.size()){
		if(event_call == 1){
			queue.push_back(index_call);
		}

		int current = -1;
		if(event_call == 1){
			current = index_call+1;
			while(current < calls.size() && calls[current].time-time < g_params.EPS){
				queue.push_back(current);
				current += 1;
			}
		}

		if(queue.size() > 1){
			try{
				CallModel cm(env,*this);
				cm.solve();
				for(int i = 0; i < queue.size(); ++i){
					auto& best_amb = ambulances[cm.opt_which_amb[i]];
					double min_time = cm.arrival_times[i] - time;
					double waiting_to_hospital_i = best_amb.answer_call(calls[queue[i]], 
						travel, ins, time, min_time, nearest_base[queue[i]]);
					waiting_on_scene[queue[i]] = cm.arrival_times[i] - calls[queue[i]].time;
					auto& call = calls[queue[i]];
					waiting_on_scene_penalized[queue[i]] = waiting_on_scene[queue[i]]*ins.penalties[call.priority];
					waiting_to_hospital[queue[i]] = waiting_to_hospital_i;
					which_ambulance[queue[i]] = cm.opt_which_amb[i];
					calls_end[queue[i]] = calls[queue[i]].end;	
					calls_attended++;
				}
				// cm.print_results();
				if(event_call == 1)
					index_call = current-1;
				queue.clear();
				set_next_event(event_call, index_call);
				continue;
			}catch(GRBException& ex){
				std::cout << ex.getMessage() << "\n";
				exit(1);
			}
		}

		
		std::vector<double> min_times(queue.size(), GRB_INFINITY);
		std::vector<std::vector<double>> travel_times(queue.size(), 
			std::vector<double>(ambulances.size(), GRB_INFINITY));
		std::vector<int> index_ambs(queue.size(), -1);

		for(int i = 0; i < queue.size(); ++i){
			auto& call = calls[queue[i]];
			for(auto& amb: ambulances){
				int j = amb.id;
				if(can_answer(amb,call)){
					travel_times[i][j] = travel.get_response_time(amb, call, time);

					if(travel_times[i][j] < min_times[i]){
						min_times[i] = travel_times[i][j];
						index_ambs[i] = j;
					}
				}
			}
			min_times[i] += time - call.time;
		}

		std::vector<int> remaining_indexes(queue.size(),0);
		for(int i = 0; i < queue.size(); ++i){
			remaining_indexes[i] = i;
		}

		std::vector<int> queue_aux;
		int nb_treated = 0;
		int total_queue = queue.size();

		while(nb_treated < total_queue){
			auto max_it = std::max_element(min_times.begin(), min_times.end());
			double max_time = *max_it;
			int current_call = std::distance(min_times.begin(), max_it);
			auto& call = calls[queue[remaining_indexes[current_call]]];
			int index_amb = index_ambs[remaining_indexes[current_call]];
			auto& best_amb = ambulances[index_amb];

			if(index_amb >= 0 && best_amb.arrival_time_at_f_last_trip <= time){
				int call_ind = queue[remaining_indexes[current_call]];
				double min_time = max_time - (time - call.time);
				double waiting_to_hospital_i = best_amb.answer_call(call, travel, ins,
					time, min_time, nearest_base[call_ind]);

				waiting_on_scene[call_ind] = max_time;
				waiting_on_scene_penalized[call_ind] = waiting_on_scene[call_ind]*ins.penalties[call.priority];
				waiting_to_hospital[call_ind] = waiting_to_hospital_i;
				which_ambulance[call_ind] = index_amb;
				calls_end[call_ind] = call.end;


				calls_attended++;
				//update travel times
				for(int k = 0; k < total_queue-nb_treated; ++k){
					int ind = remaining_indexes[k];
					if(k != current_call && index_amb == index_ambs[ind]){
						double time_to_h = best_amb.arrival_time_at_f_last_trip -
							time;
						double time_from_h_to_c = travel.travel_time(best_amb.free_location,
							calls[queue[ind]].location);
						travel_times[ind][best_amb.id] = time_to_h + time_from_h_to_c;
						auto min_it = std::min_element(travel_times[ind].begin(),
							travel_times[ind].end());
						double min_val = *min_it;
						int min_ind = std::distance(travel_times[ind].begin(), min_it);
						min_times[k] = min_val + time - calls[queue[ind]].time;
						index_ambs[ind] = min_ind;
					}
				}
			}else{
				queue_aux.push_back(queue[remaining_indexes[current_call]]);
			}
			if(nb_treated+1 < total_queue){
				//erase current_call from min_times e remaining_indexes
				min_times.erase(max_it);
				remaining_indexes.erase(std::remove(remaining_indexes.begin(), 
					remaining_indexes.end(), remaining_indexes[current_call]), 
					remaining_indexes.end());
			}
			nb_treated++;
		}

		queue.clear();
		for(auto i: queue_aux){
			queue.push_back(i);
		}
		set_next_event(event_call, index_call);
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



MinMaxSolver::MinMaxSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
	Instance& ins, Travel& travel): Solver(env, calls,ambulances,ins,travel){
	travel.set_forward(false);
	queue.reserve(calls.size());
}

void MinMaxSolver::run(){
	int event_call = 1;
	time = calls[0].time;
	int index_call = 0;
	int calls_attended = 0;
	travel.set_forward(true);
	while(calls_attended < calls.size()){
		if(event_call == 1){
			queue.push_back(index_call);
		}

		std::vector<double> min_times(queue.size(), GRB_INFINITY);
		std::vector<std::vector<double>> travel_times(queue.size(), 
			std::vector<double>(ambulances.size(), GRB_INFINITY));
		std::vector<int> index_ambs(queue.size(), -1);

		// For every call in queue, determine the best ambulance and the best
		// response time
		for(int i = 0; i < queue.size(); ++i){
			auto& call = calls[queue[i]];
			for(auto& amb: ambulances){
				int j = amb.id;
				if(can_answer(amb,call)){
					travel_times[i][j] = travel.get_response_time(amb, call, time);

					if(travel_times[i][j] < min_times[i]){
						min_times[i] = travel_times[i][j];
						index_ambs[i] = j;
					}else if(abs(travel_times[i][j] - min_times[i]) <  g_params.EPS && 
						index_ambs[i] != -1 && amb.type > ambulances[index_ambs[i]].type){
						min_times[i] = travel_times[i][j];
						index_ambs[i] = j;
					}
				}
			}
			min_times[i] += time - call.time;
		}


		std::vector<int> remaining_indexes(queue.size(),0);
		for(int i = 0; i < queue.size(); ++i){
			remaining_indexes[i] = i;
		}

		std::vector<int> queue_aux;
		int nb_treated = 0;
		int total_queue = queue.size();

		while(nb_treated < total_queue){
			auto t0 = std::chrono::high_resolution_clock::now();
			//get time and index of the worst call, and the best ambulance to such call
			auto max_it = std::max_element(min_times.begin(), min_times.end());
			double max_time = *max_it;
			int current_call = std::distance(min_times.begin(), max_it);
			int call_ind = queue[remaining_indexes[current_call]];
			auto& call = calls[call_ind];
			int index_amb = index_ambs[remaining_indexes[current_call]];
			auto& best_amb = ambulances[index_amb];

			//if best ambulance is free, it should answer the call
			if(index_amb >= 0 && best_amb.arrival_time_at_f_last_trip <= time){
				double waiting_on_scene_i = max_time;
				double waiting_to_hospital_i = ambulances[index_amb].answer_call(call, 
					travel, ins, time, max_time - (time - call.time), nearest_base[call_ind]);
				waiting_on_scene[call_ind] = waiting_on_scene_i;
				waiting_on_scene_penalized[call_ind] = waiting_on_scene_i*ins.penalties[call.priority];
				waiting_to_hospital[call_ind] = waiting_to_hospital_i;
				which_ambulance[call_ind] = index_amb;
				calls_end[call_ind] = call.end;
				
				calls_attended++;
				//update travel times
				for(int k = 0; k < total_queue-nb_treated; ++k){
					int ind = remaining_indexes[k];
					if(k != current_call && index_amb == index_ambs[ind]){
						double time_to_h = best_amb.arrival_time_at_f_last_trip -
							time;
						double time_from_h_to_c = travel.travel_time(best_amb.free_location,
							calls[queue[ind]].location);
						travel_times[ind][best_amb.id] = time_to_h + time_from_h_to_c;
						auto min_it = std::min_element(travel_times[ind].begin(),
							travel_times[ind].end());
						double min_val = *min_it;
						int min_ind = std::distance(travel_times[ind].begin(), min_it);
						min_times[k] = min_val + time - calls[queue[ind]].time;
						index_ambs[ind] = min_ind;
					}
				}
			}else{
				queue_aux.push_back(queue[remaining_indexes[current_call]]);
			}
			if(nb_treated+1 < total_queue){
				min_times.erase(max_it);
				remaining_indexes.erase(std::remove(remaining_indexes.begin(), 
					remaining_indexes.end(), remaining_indexes[current_call]), remaining_indexes.end());
			}
			nb_treated++;
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
