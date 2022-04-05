#ifndef _POLICY_TESTER_H
#define _POLICY_TESTER_H

#include "main.h"
#include "instance.h"

class Instance;
class OSRMHelper;
class Travel;

using namespace std;
struct SecondStageWaitingTime{
	double waiting_on_scene;
	double waiting_to_hospital;
	double mean_total_time;
};


struct Stats{
	double mean_waiting_on_scene;
	double mean_waiting_on_scene_penalized;
	double mean_waiting_to_hospital;
	double max_waiting_on_scene;
	double max_waiting_on_scene_penalized;
	double max_waiting_to_hospital;
	double waiting_on_scene_q90;
	double waiting_to_hospital_q90;
	double mean_total;
	double max_total;
	double q90_total;

	Stats(vector<double> waiting_on_scene, vector<double> waiting_on_scene_penalized, 
		vector<double> waiting_to_hospital);
};

class PolicyTester{
public:
	PolicyTester(Instance& ins);
	~PolicyTester();


	const string real_data1 = 
	"calibration/Scenarios_For_Tests/Real_Data_Continuous_Scenarios/baseScenario.txt";

	GRBEnv env;
	Instance& ins;
	Travel& travel;
	vector<Ambulance> ambulances;
	double time;

	vector<vector<double>> waiting_on_scene;
	vector<vector<double>> waiting_on_scene_penalized;
	vector<vector<double>> waiting_to_hospital;
	vector<vector<double>> calls_end;
	vector<vector<int>> which_ambulance;


	vector<double> run_times;


	// const vector<string> policies{"queue", "forward", "priorities", "minmax",
	// 	"gen_forward", "gen_minmax","non_miopyc"};

	const vector<string> policies{"queue", "forward", "priorities", "minmax", "non_miopyc"};
	const map<string,string> policies_names{{"queue", "CA"}, {"forward", "BM"}, {"priorities", 
		"GHP1"}, {"minmax", "GHP2"}, {"non_miopyc", "NM"}};

	shared_ptr<Solver> get_solver(const string& policy, vector<Call>& calls, 
	vector<Ambulance>& ambulances, Travel& travel);

	void run();

	void one_stage_old();
	void one_stage();
	void two_stage();
	void two_stage_tree();


	int set_next_event(vector<Call>& this_scenario, int& event_call, int& index_call);

	SecondStageWaitingTime get_waiting_time(int amb_id, Call& call, const string& policy,
		vector<int>& index_scenarios, int nearest_base_id, std::vector<int> queue, int sc,
		double time);

	vector<vector<Call>> get_future_scenarios(double time, vector<Call>& this_scenario,
	int nb_realizations, int T, vector<vector<vector<Call>>>& my_scenarios);

	SecondStageWaitingTime get_waiting_time_tree(int amb_id,  int call_id,
		const string& policy, vector<vector<Call>>& future_scenarios, 
		vector<Call>& this_scenario, double time, 
		std::vector<int>& queue, int base);

	SecondStageWaitingTime get_waiting_time_tree_return(int amb_id, int base,
		const string& policy, vector<vector<Call>>& future_scenarios, 
		vector<Call>& this_scenario, double time, 
		std::vector<int>& queue);

	SecondStageWaitingTime get_waiting_time_return(int amb_id, const string& policy,
		vector<int>& index_scenarios, double time, std::vector<int> queue, int sc, int base);

	vector<int> get_nearest_base(vector<Call>& calls);


	bool can_answer(Ambulance& amb, Call& call);
};


#endif