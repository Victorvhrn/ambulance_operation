#include "../include/instance.h"
#include "../include/travel.h"
#include "../include/solver.h"
#include "../include/call.h"
#include "../include/ambulance.h"


Ambulance::Ambulance(int id, Location& base_location, Location& free_location,
		double arrival_time_at_f_last_trip, double arrival_time_at_b_last_trip, 
		int type, double speed): 
		call(NULL), id(id), base_location(base_location), 
		free_location(free_location), clean_location(null_location),
		arrival_time_at_f_last_trip(arrival_time_at_f_last_trip),
		arrival_time_at_b_last_trip(arrival_time_at_b_last_trip),
		type(type), speed(speed), busy(false), 
		last_origin_location(base_location), departure_time(0){
}

double Ambulance::answer_call(Call& call, Travel& travel, Instance& ins, double time, 
	double min_time, int nearest_base_id){
	this->call = &call;
	last_origin_location = call.location;
	departure_time = time;

	if(arrival_time_at_b_last_trip <= time){
		set_new_point(time, base_location, TripType::AT_BASE);
	}else if(arrival_time_at_f_last_trip < time){
		auto current_location = travel.ambulance_position(*this, time);
		times.back() = time;
		trips.back() = current_location;
	}else{
		times.erase(times.begin() + times.size() - 1);
		trips.erase(trips.begin() + trips.size() - 1);
		trip_types.erase(trip_types.begin() + trip_types.size()-2, 
			trip_types.end());
	}

	double time_arrival_scene = time + min_time;
	set_new_point(time_arrival_scene, call.location, TripType::TO_CALL);

	double time_leave_scene = time_arrival_scene + call.time_on_scene;
	set_new_point(time_leave_scene, call.location, TripType::AT_SCENE);

	double waiting_to_hospital_i = GRB_INFINITY;
	Location base;
	if(!g_params.h_use_fixed_bases)
		base = ins.bases[nearest_base_id];
	else
		base = base_location;

	if(call.clean_needed && call.hosp_needed){
		auto& hospital = ins.hospitals[call.hospital];
		auto& cb = ins.cleaning_bases[call.cleaning];
		double t_from_c_to_h = travel.travel_time(call.location, hospital);
		double t_from_h_to_cb = travel.travel_time(hospital, cb);
		double t_from_cb_to_b = travel.travel_time(cb,base);
		double time_leave_hospital = time_leave_scene + t_from_c_to_h +
			call.time_at_hospital;
		double time_arrival_cb = time_leave_hospital + t_from_h_to_cb;
		double time_leave_cb = time_arrival_cb + call.cleaning_time;
		double time_arrival_base = time_leave_cb + t_from_cb_to_b;

		call.end = time_leave_scene + t_from_c_to_h;
		waiting_to_hospital_i = call.time_on_scene + t_from_c_to_h;

		free_location = cb;
		base_location = base;
		arrival_time_at_f_last_trip = time_leave_cb;
		arrival_time_at_b_last_trip = time_arrival_base;

		set_new_point(time_leave_hospital - call.time_at_hospital, hospital,
			TripType::TO_HOSPITAL);
		set_new_point(time_leave_hospital, hospital, TripType::AT_HOSPITAL);

		set_new_point(time_arrival_cb, cb, TripType::TO_CB);
		set_new_point(time_leave_cb, cb, TripType::AT_CB);

		set_new_point(time_arrival_base, base_location, TripType::TO_BASE);
	}else if(call.clean_needed && !call.hosp_needed){
		auto& cb = ins.cleaning_bases[call.cleaning];
		double t_from_c_to_cb = travel.travel_time(call.location, cb);
		double t_from_cb_to_b = travel.travel_time(cb,base);
		double time_arrival_cb = time_leave_scene + t_from_c_to_cb;
		double time_leave_cb = time_arrival_cb + call.cleaning_time;
		double time_arrival_base = time_leave_cb + t_from_cb_to_b;

		call.end = time_leave_scene;
		waiting_to_hospital_i = call.time_on_scene;

		free_location = cb;
		base_location = base;
		arrival_time_at_f_last_trip = time_leave_cb;
		arrival_time_at_b_last_trip = time_arrival_base;


		set_new_point(time_arrival_cb, cb, TripType::TO_CB);
		set_new_point(time_leave_cb, cb, TripType::AT_CB);
		set_new_point(time_arrival_base, base_location, TripType::TO_BASE);

	}else if(!call.clean_needed && call.hosp_needed){
		auto& hospital = ins.hospitals[call.hospital];
		double t_from_c_to_h = travel.travel_time(call.location, hospital);
		double t_from_h_to_b = travel.travel_time(hospital, base);

		double time_leave_hospital = time_leave_scene + t_from_c_to_h +
			call.time_at_hospital;

		double time_arrival_base = time_leave_hospital + t_from_h_to_b;

		call.end = time_leave_scene + t_from_c_to_h;
		waiting_to_hospital_i = call.time_on_scene + t_from_c_to_h;

		free_location = hospital;
		base_location = base;
		arrival_time_at_f_last_trip = time_leave_hospital;
		arrival_time_at_b_last_trip = time_arrival_base;


		set_new_point(time_leave_hospital - call.time_at_hospital, hospital,
			TripType::TO_HOSPITAL);
		set_new_point(time_leave_hospital, hospital, TripType::AT_HOSPITAL);
		set_new_point(time_arrival_base, base_location, TripType::TO_BASE);

	}else if(!call.clean_needed && !call.hosp_needed){
		double t_from_c_to_b = travel.travel_time(call.location, base);

		double time_arrival_base = time_leave_scene + t_from_c_to_b;

		call.end = time_leave_scene;
		waiting_to_hospital_i = call.time_on_scene;

		free_location = call.location;
		base_location = base;
		arrival_time_at_f_last_trip = time_leave_scene;
		arrival_time_at_b_last_trip = time_arrival_base;

		set_new_point(time_arrival_base, base_location, TripType::TO_BASE);
	}

	return waiting_to_hospital_i;
}




void Ambulance::set_new_point(double time, Location& location, TripType trip_type){
	// fmt::print("time location {} {} {}\n", time, location.first, location.second);
	times.push_back(time);
	trips.push_back(location);
	trip_types.push_back(trip_type);
}





std::ostream& operator<<(std::ostream& out, const Ambulance& amb){
	out << "Amb " << amb.id << ", type " << amb.type << " | ";
	out << "Arr b " << amb.arrival_time_at_b_last_trip << ", Arr f ";
	out << amb.arrival_time_at_f_last_trip;
	// out << ", last call: ";
	// if(amb.call != NULL){
	// 	out << amb.call->id;
	// }else{
	// 	out << " - ";
	// }
	return out;
}

Ambulance::~Ambulance(){}