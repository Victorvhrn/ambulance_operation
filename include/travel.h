#ifndef _TRAVEL_H
#define _TRAVEL_H

#include "main.h"

class Instance;
class Call;
class Ambulance;

class Travel{

public:
	Travel(bool euclidian = true);

	bool euclidian;
	bool forward;

	void set_forward(bool a_forward);
	double get_response_time(Ambulance& amb, Call& call, double current_time,
		bool force_forward=false);
	double travel_time(Location& a, Location& b, double speed);
	double travel_time(Location& a, Location& b);
	Location ambulance_position(Ambulance& amb, double t);
	double travel_time_from_position(Ambulance& amb,  Call& call);
	double norm(Location& a, Location& b);

	std::vector<std::vector<double>> table_in_out(std::vector<Location>& in, 
        std::vector<Location>& out);

	double lat_long_travel_time(Location& a, Location& b);
	double lat_long_travel_distance(Location& a, Location& b);
	Location position_between_origin_destination(Location& a, Location& b, 
	double t0, double t);
};


#endif