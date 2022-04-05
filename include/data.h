#ifndef DATA_H
#define DATA_H

#include "main.h"
#include "ambulance.h"
#include "instance.h"
// #include "event_simulator.h"

class Ambulance;
class Call;
class Travel;
class Instance;
class Solver;
class CGSolver;


namespace bg = boost::geometry;

typedef bg::model::point<double, 2, bg::cs::cartesian> point_t;
typedef bg::model::polygon<point_t> polygon_t;

class Data{
public:
	// Data(const Data & data);
	Data(Solver* solver, Call* call, Ambulance* ambulance);
	~Data();
	std::string path;

	Solver* solver;
	Instance& ins;
	Call* call;
	Ambulance* ambulance;
	Travel travel;

	int num_hosps, num_bases, num_locals, num_times, types_call, 
		types_amb, num_ambs, total_locals, num_calls;
		
	int* hospitals; // Index of hospitals
	int* bases; // Index of bases
	int* cap_bases; //Capacity of bases

	double quantum;
	double max_cost;

	bool* is_hospital;
	bool* is_base;
	
	std::vector<polygon_t> regions;
	std::vector<Location> locations;

	std::vector<std::pair<int,int>> amb_nodes_map;

	std::vector<std::set<int>> A; // Types of calls attended by each ambulance.
	std::set<int>** H; // Hospitals that can answer call c at location l
	int** C; //Initial queue

	int** A0_ab; //Ambulances type a in base b.
	int** A0_ah; //Ambulances type a in hospital h.
	int*** A0_alb; //Ambulances type a at location l heading to b.
	int**** A0_calh; //Ambulances type a at location l 
	// answering call c heading to h.
	int***** A0_callh; //Ambulances type a at location l1 heading to
	// l to answer a call type c and take patients to hospital h.
	std::set<int>* G; // Graph
	std::vector<std::tuple<int,int,int>> calls;
	int*** lambda; //Forecasted calls
	// double** times; // travel times
	std::set<int>*** L_tab;
	double****** tao; //Forecasted time
	double**** tao_0; // Forecasted time
	std::map<std::string, bool> relax; //Relaxed variables on the model
	int local_call0, type_call0; //c0,l0
	bool debug; //Debug mode

	void get_times_matrix();
	void set_regions();
	void set_times();

	int L(int t, int a, int l1, int b);
	int get_time_slot(double time);
	int get_location_index(Location& location);
};

#endif