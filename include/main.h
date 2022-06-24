#ifndef _MAIN_H
#define _MAIN_H


#include "param.h"
#include <iostream>
#include <chrono>
#include <vector>
#include "gurobi_c++.h"
#include <cstdio>
#include <cstdlib>
#include <set>
#include <map>
#include <float.h>
#include <iomanip>
#include <fstream>
#include <pthread.h>
#include <algorithm>
#include <memory>
#include <numeric>
#include <random>
#include <cmath>
#include <queue>
#include <tuple>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <string>
#include <unordered_map>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/multi_array.hpp>
#include <cassert>
#include <sstream>
#include <climits>
#include <sys/wait.h>


#include <osrm/match_parameters.hpp>
#include <osrm/nearest_parameters.hpp>
#include <osrm/route_parameters.hpp>
#include <osrm/table_parameters.hpp>
#include <osrm/trip_parameters.hpp>

#include <osrm/coordinate.hpp>
#include <osrm/engine_config.hpp>
#include <osrm/json_container.hpp>

#include <osrm/osrm.hpp>
#include <osrm/status.hpp>

#include <xtensor/xarray.hpp>
#include <xtensor/xbuilder.hpp>
#include <xtensor/xstrided_view.hpp>

#define BASE 0
#define HOSP 1
#define LOC 2

class Data;
// class Simulator;

void version();
Param load_params(int argc, char** argv);

typedef std::pair<double, double> Location;
const Location null_location(-100,-100);

class Ambulance;
class Data;

typedef struct {
	int thread_id;
	Data* data_obj;
	// GRBEnv* env;
	double* costs;
	Ambulance& amb;
	int*** lambda;
	int type_call0;
	int local_call0;
	int hospital;
	int base_ret;
	bool is_call;
} thread_data;

class Node{
public:
	int id;
	int val;
	Node(int id, int val): id(id), val(val){}
	~Node(){}
};

class Compare{
public:
	bool operator()(const Node & u, const Node& v){
		return u.val > v.val;
	}
};

double f(int t, int c, int a, int b_h1, int l, int h, Data& data,
	int factor);
double f(int t, int c, int a, int l1, int b, int l, int h,
	Data& data, int factor);
double g_tcl(int t, int c, int l, double factor);
double g_tab(int t, int a, int b, double factor);
double g_talb(int t, int a, int l1, int b, Data& data, double factor);



//Parameter global object
extern Param g_params;

#endif