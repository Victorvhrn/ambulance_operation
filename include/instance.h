#ifndef _INSTANCE
#define _INSTANCE

#include "main.h"
#include "travel.h"

using namespace std;

class Call;
class Ambulance;
class Travel;
class OSRMHelper;

const string simulated_rect_1 = 
"calibration/Scenarios_For_Tests/Simulated_Data_Rectangle/simualtedRectangleUniform.txt";
const string simulated_rect_2 = 
"calibration/Scenarios_For_Tests/Simulated_Data_Rectangle/simualtedRectanglePoisson.txt";



class Instance{
public:
	Instance();
	explicit Instance(string path);
	~Instance();
	
	int nb_scenarios;
	vector<int> nb_calls;
	int nb_hospitals;
	int nb_bases;
	int nb_cleaning_bases;
	int nb_types_ambulance;
	int nb_ambulances;
	int nb_priorities;
	int nb_regions;
	int nb_times;

	string path;
	Travel travel;

	double x_min, x_max, y_min, y_max;
	double slot_duration;
	vector<set<int>> A; //for each priority c, set of ambulance types that can respond to c

	vector<Location> centers;
	vector<Location> bases;
	vector<Location> cleaning_bases;
	vector<Location> hospitals;
	vector<int> cap_bases;
	vector<double> penalties;
	vector<int> nearest_base_to_region;

	vector<vector<Call>> calls;
	vector<Ambulance> ambulances;
	vector<vector<vector<Call>>> scenarios_by_day;


	void load_euclidian();
	void load_real_data();
	

};

#endif