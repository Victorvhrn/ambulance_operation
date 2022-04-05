#include "../include/osrm_helper.h"
#include "../include/call.h"
#include "../include/ambulance.h"
#include "../include/travel.h"


Travel::Travel(bool euclidian): osrm(), euclidian(euclidian), forward(false){}



void Travel::set_forward(bool a_forward){
	forward = a_forward;
}


double Travel::norm(Location& a, Location& b){
	if(euclidian){
		double x1,y1,x2,y2;
		std::tie(x1,y1) = a;
		std::tie(x2,y2) = b;
		return sqrt(pow(x2-x1,2) + pow(y2-y1,2));
	}else{
		// double result = osrm.get_distance(a,b);
		// return (result > g_params.EPS) ? result : lat_long_travel_distance(a,b);
		return lat_long_travel_distance(a,b);
	}
}


double Travel::travel_time(Location& a, Location& b, double speed){
	return norm(a,b)/speed;
}

double Travel::travel_time(Location& a, Location& b){
	// double result = osrm.get_duration(a,b);
	// return (result > g_params.EPS) ? result : lat_long_travel_time(a,b);
	return lat_long_travel_time(a,b);
}

double Travel::get_response_time(Ambulance& amb, Call& call, double current_time, 
	bool force_forward){

	double time_this_ambulance = GRB_INFINITY;

	if(amb.arrival_time_at_b_last_trip <= current_time){
		if(euclidian)
			time_this_ambulance = travel_time(amb.base_location,call.location, amb.speed);
		else
			time_this_ambulance = travel_time(amb.base_location,call.location);
	}else if(amb.arrival_time_at_f_last_trip <= current_time){
		Location current_location = ambulance_position(amb,current_time);
		if(euclidian)
			time_this_ambulance = travel_time(current_location, call.location, amb.speed);
		else
			time_this_ambulance = travel_time(current_location, call.location);
	}else if(forward || force_forward){
		double time_to_free = amb.arrival_time_at_f_last_trip-current_time;
		double time_free_to_call = GRB_INFINITY;
		if(euclidian){
			time_free_to_call = travel_time(amb.free_location, call.location, amb.speed);
		}else{
			time_free_to_call = travel_time(amb.free_location, call.location);
		}
		time_this_ambulance = time_to_free + time_free_to_call;
	}

	return time_this_ambulance;
}


Location Travel::ambulance_position(Ambulance& amb, double t){
	Location free_location = amb.free_location;
	double t0 = amb.arrival_time_at_f_last_trip;

	if(euclidian){
		double d = norm(free_location, amb.base_location);
		double ttravel = d/amb.speed;
		Location location = free_location;
		if(t < t0+ttravel){
			double x1,x2,y1,y2;
			std::tie(x1,x2) = free_location;
			std::tie(y1,y2) = amb.base_location;
			Location direction;
			direction.first = (y1-x1)/d;
			direction.second = (y2-x2)/d;
			location.first = x1 + amb.speed*direction.first*(t-t0);
			location.second = x2 + amb.speed*direction.second*(t-t0);
			return location;
		}else{
			return amb.base_location;
		}
	}else{
		const double R = 6371;
		const double radian = M_PI/180;
		const double v = 60;

		double p11 = R*cos(radian*free_location.first)*cos(radian*free_location.second);
		double p12 = R*cos(radian*free_location.first)*sin(radian*free_location.second);
		double p13 = R*sin(radian*free_location.first);

		double p21 = R*cos(radian*amb.base_location.first)*cos(radian*amb.base_location.second);
		double p22 = R*cos(radian*amb.base_location.first)*sin(radian*amb.base_location.second);
		double p23 = R*sin(radian*amb.base_location.first);

		double d = sqrt(pow(p11 - p21, 2) + pow(p12 - p22, 2) + pow(p13 - p23, 2));
		double alpha = 2*asin(d/(2*R));
		double dearth = R*alpha;
		double ttravel = dearth/v;
		double t0 = amb.arrival_time_at_f_last_trip;
		Location location = amb.free_location;
		if(t < t0 + ttravel){
			double alpha0 = (t-t0)*v/R;
			double beta = sin(alpha-alpha0)/sin(alpha);
			double gamma = cos(alpha-alpha0)-sin(alpha-alpha0)*cos(alpha)/sin(alpha);
			double curr_p1 = beta*p11 + gamma*p21;
			double curr_p2 = beta*p12 + gamma*p22;
			double curr_p3 = beta*p13 + gamma*p23;
			location.first = (asin(curr_p3/R)/M_PI)*180;
			if(curr_p2 > -g_params.EPS){
				location.second = (acos(curr_p1/sqrt(pow(R,2)-pow(curr_p3,2)))/M_PI)*180;
			}else{
				location.second = -(acos(curr_p1/sqrt(pow(R,2)-pow(curr_p3,2)))/M_PI)*180;
			}
			return location;
		}else{
			return amb.base_location;
		}
	}
}


double Travel::travel_time_from_position(Ambulance& amb,  Call& call){
	if(euclidian){
		double d = norm(amb.free_location, amb.base_location);
		double speed = amb.speed;
		double t =  call.time;
		double ttravel = d/speed;
		double t0 = amb.arrival_time_at_f_last_trip;
		Location location = amb.free_location;
		if(t < t0+ttravel){
			double x1,x2,y1,y2;
			std::tie(x1,x2) = amb.free_location;
			std::tie(y1,y2) = amb.base_location;
			Location direction;
			direction.first = (y1-x1)/d;
			direction.second = (y2-x2)/d;
			location.first = x1 + amb.speed*direction.first*(t-t0);
			location.second = x2 + amb.speed*direction.second*(t-t0);
			return travel_time(location,call.location, speed);
		}else{
			return travel_time(amb.base_location, call.location, speed);
		}
	}else{
		Location location = ambulance_position(amb, call.time);
		// double result = osrm.get_duration(location, call.location);
		// return (result < g_params.EPS) ? lat_long_travel_time(location, call.location) :
		// 	result;
		return travel_time(location, call.location);
	}
}


Location Travel::position_between_origin_destination(Location& a, Location& b, 
	double t0, double t){

	if(euclidian){
		double speed = 0.078;
		double d = norm(a, b);
		double ttravel = d/speed;
		Location location = a;
		if(t < t0+ttravel){
			double x1,x2,y1,y2;
			std::tie(x1,x2) = a;
			std::tie(y1,y2) = b;
			Location direction;
			direction.first = (y1-x1)/d;
			direction.second = (y2-x2)/d;
			location.first = x1 + speed*direction.first*(t-t0);
			location.second = x2 + speed*direction.second*(t-t0);
			return location;
		}else{
			return b;
		}
	}else{
		const double R = 6371;
		const double radian = M_PI/180;
		const double v = 60;

		double p11 = R*cos(radian*a.first)*cos(radian*a.second);
		double p12 = R*cos(radian*a.first)*sin(radian*a.second);
		double p13 = R*sin(radian*a.first);

		double p21 = R*cos(radian*b.first)*cos(radian*b.second);
		double p22 = R*cos(radian*b.first)*sin(radian*b.second);
		double p23 = R*sin(radian*b.first);

		double d = sqrt(pow(p11 - p21, 2) + pow(p12 - p22, 2) + pow(p13 - p23, 2));
		double alpha = 2*asin(d/(2*R));
		double dearth = R*alpha;
		double ttravel = dearth/v;
		Location location = a;
		if(t < t0 + ttravel){
			double alpha0 = (t-t0)*v/R;
			double beta = sin(alpha-alpha0)/sin(alpha);
			double gamma = cos(alpha-alpha0)-sin(alpha-alpha0)*cos(alpha)/sin(alpha);
			double curr_p1 = beta*p11 + gamma*p21;
			double curr_p2 = beta*p12 + gamma*p22;
			double curr_p3 = beta*p13 + gamma*p23;
			location.first = (asin(curr_p3/R)/M_PI)*180;
			if(curr_p2 > -g_params.EPS){
				location.second = (acos(curr_p1/sqrt(pow(R,2)-pow(curr_p3,2)))/M_PI)*180;
			}else{
				location.second = -(acos(curr_p1/sqrt(pow(R,2)-pow(curr_p3,2)))/M_PI)*180;
			}
			return location;
		}else{
			return b;
		}
	}
}


std::vector<std::vector<double>> Travel::table_in_out(std::vector<Location>& in, 
	std::vector<Location>& out){
	// auto result = osrm.table_in_out(in,out);
	auto result = std::vector<std::vector<double>>(in.size(), 
		std::vector<double>(out.size(), 0));
	for(int i = 0; i < in.size(); ++i){
		for(int j = 0; j < out.size(); ++j){
			result[i][j] = lat_long_travel_time(in[i], out[j]);
		}
	}

	return result;
}


double Travel::lat_long_travel_time(Location& a, Location& b){
	const double v = 60;
	if(euclidian){
		double speed = 0.078;
		return travel_time(a,b,speed);
	}else{
		return 3600*(lat_long_travel_distance(a,b)/v);
	}
}


double Travel::lat_long_travel_distance(Location& a, Location& b){
	if(euclidian){
		return norm(a,b);
	}

	const double R = 6371;
	const double radian = M_PI/180;

	double p11 = R*cos(radian*a.first)*cos(radian*a.second);
	double p12 = R*cos(radian*a.first)*sin(radian*a.second);
	double p13 = R*sin(radian*a.first);

	double p21 = R*cos(radian*b.first)*cos(radian*b.second);
	double p22 = R*cos(radian*b.first)*sin(radian*b.second);
	double p23 = R*sin(radian*b.first);

	double d = sqrt(pow(p11 - p21, 2) + pow(p12 - p22, 2) + pow(p13 - p23, 2));
	return 2*R*asin(d/(2*R));
}