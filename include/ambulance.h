#ifndef _AMBULANCE_H
#define _AMBULANCE_H


#include "main.h"

class Simulator;
class Solver;
class Instance;
class Travel;
class Call;

using std::shared_ptr;

enum class TripType{
	AT_BASE, //matlab type 1
	TO_CALL, //matlab type 2
	AT_SCENE, //matlab type 3
	TO_HOSPITAL,//matlab type 4
	AT_HOSPITAL,//matlab type 5
	TO_CB, //matlab type 6
	AT_CB, //matlab type 7
	TO_BASE //matlab type 8
};

class Ambulance{
public:
	Ambulance();
	~Ambulance();

	Ambulance(int id, Location& base_location, Location& free_location,
			double arrival_time_at_f_last_trip,
			double arrival_time_at_b_last_trip, 
			int type, double speed);
	Call* call;
	Location base_location;
	Location free_location;
	Location clean_location;
	Location last_origin_location;
	int id;
	bool busy;
	double arrival_time_at_f_last_trip; //free
	double arrival_time_at_b_last_trip;
	double departure_time;

	std::vector<double> times;
	std::vector<Location> trips;
	std::vector<TripType> trip_types;
	int type;
	double speed;

	double answer_call(Call& call, Travel& travel, Instance& ins, double time, double min_time, 
		int nearest_base_id);

	void set_new_point(double time, Location& location, TripType trip_type);
	// bool is_at_hospital(double time, Travel& travel);
	// Location get_last_free_location();

	friend std::ostream& operator<<(std::ostream& out, const Ambulance& amb);
	
};


#endif