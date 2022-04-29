#ifndef _SOLVER_H
#define _SOLVER_H

using namespace std;

#include "main.h"
#include "call.h"
#include "ambulance.h"
#include "travel.h"

class Travel;
class CGCall;

class Solver{
protected:
	Solver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel);
public:
	virtual ~Solver();

	GRBEnv& env;
	vector<Call> calls;
	vector<Ambulance> ambulances;
	Travel travel;
	Instance& ins;

	double time;
	double first_time;
	double last_time;

	vector<double> waiting_on_scene;
	vector<double> waiting_on_scene_penalized;
	vector<double> waiting_to_hospital;
	vector<double> calls_end;
	vector<int> which_ambulance;

	vector<int> queue;
	vector<int> nearest_base;
	double obj;

	std::vector<double> run_times;

	virtual void run() = 0;
	void set_next_event(int& event_call, int& index_call);
	void set_calls_nearest_bases();
	void print_results();
	bool can_answer(Ambulance& amb, Call& call);
};


class ForwardSolver: public Solver{
public:
	ForwardSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel);
	
	virtual void run();
};

class QueueSolver: public Solver{
public:
	QueueSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel);
	
	virtual void run();
};

class PrioritySolver: public Solver{
public:
	PrioritySolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel);

	virtual void run();
};

class MinMaxSolver: public Solver{
public:
	MinMaxSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel);

	virtual void run();	
};


class MinMaxPSolver: public Solver{
public:
	MinMaxPSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel);

	virtual void run();
};

class NonMiopycSolver: public Solver{
public:
	NonMiopycSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel);

	virtual void run();
	void compute_min_times(int k, vector<double>& min_times, vector<int>& index_amb, 
		vector<int> & amb_type);
};

class GenForwardSolver: public Solver{
public:
	GenForwardSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel);
	virtual void run();
};


class GenMinMaxSolver: public Solver{
public:
	GenMinMaxSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel);

	virtual void run();
};


class DeterministicSolver: public Solver{
public:
	DeterministicSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel);

	virtual void run();
};

class CGSolver: public Solver{
public:
	CGSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel);

	int get_return_base(CGCall& cg, Data& data, int t, int c, int a, int l1, 
		int l, int h);
	virtual void run();
};


class ModelSolver: public Solver{
public:
	ModelSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel);

	virtual void run();
};
 

#endif