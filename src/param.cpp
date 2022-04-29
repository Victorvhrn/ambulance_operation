#include "../include/param.h"

Param::Param(){}

Param::Param(const Param &p): instance(p.instance), 
	solver(p.solver),
	generator_folder(p.generator_folder),
	debug(p.debug),
	h_use_fixed_bases(p.h_use_fixed_bases),
	h_forward(p.h_forward),
	h_discard(p.h_discard),
	h_order_priorities(p.h_order_priorities),
	h_order_time(p.h_order_time),
	n_nearest_hospitals(p.n_nearest_hospitals),
	n_nearest_ambs(p.n_nearest_ambs),
	n_queue_calls_eval(p.n_queue_calls_eval),
	n_scenarios(p.n_scenarios),
	n_time_slots(p.n_time_slots),
	n_regions(p.n_regions),
	n_hospitals(p.n_hospitals),
	n_bases(p.n_bases),
	n_cleaning_bases(p.n_cleaning_bases),
	n_ambulances(p.n_ambulances),
	osrm_map_path(p.osrm_map_path),
	extended_model(p.extended_model),
	h_random_hospital(p.h_random_hospital){}

Param::Param(boost::program_options::variables_map vm): 
	instance{vm["instance"].as<std::string>()},
	solver{vm["solver"].as<std::string>()},
	generator_folder{vm["generator_folder"].as<std::string>()},
	debug{vm["debug"].as<bool>()},
	h_use_fixed_bases{vm["h_use_fixed_bases"].as<bool>()},
	h_forward{vm["h_forward"].as<bool>()},
	h_discard{vm["h_discard"].as<bool>()},
	h_order_priorities{vm["h_order_priorities"].as<bool>()},
	h_order_time{vm["h_order_time"].as<bool>()},
	n_nearest_hospitals{vm["n_nearest_hospitals"].as<int>()},
	n_nearest_ambs{vm["n_nearest_ambs"].as<int>()},
	n_queue_calls_eval{vm["n_queue_calls_eval"].as<int>()},
	n_scenarios{vm["n_scenarios"].as<int>()},
	n_time_slots{vm["n_time_slots"].as<int>()},
	n_regions{vm["n_regions"].as<int>()},
	n_hospitals{vm["n_hospitals"].as<int>()},
	n_bases{vm["n_bases"].as<int>()},
	n_cleaning_bases{vm["n_cleaning_bases"].as<int>()},
	n_ambulances{vm["n_ambulances"].as<int>()},
	osrm_map_path{vm["osrm_map_path"].as<std::string>()},
	extended_model{vm["extended_model"].as<bool>()},
	h_random_hospital{vm["h_random_hospital"].as<bool>()}{
}

Param& Param::operator=(const Param& p){
	if(this == &p)
		return *this;

	instance = p.instance;
	solver = p.solver;
	generator_folder = p.generator_folder;
	debug = p.debug;
	h_use_fixed_bases = p.h_use_fixed_bases;
	h_forward = p.h_forward;
	h_discard = p.h_discard;
	h_order_priorities = p.h_order_priorities;
	h_order_time = p.h_order_time;
	n_nearest_hospitals = p.n_nearest_hospitals;
	n_nearest_ambs = p.n_nearest_ambs;
	n_queue_calls_eval = p.n_queue_calls_eval;
	n_scenarios = p.n_scenarios;
	n_time_slots = p.n_time_slots;
	n_regions = p.n_regions;
	n_hospitals = p.n_hospitals;
	n_bases = p.n_bases;
	n_cleaning_bases = p.n_cleaning_bases;
	n_ambulances = p.n_ambulances;
	osrm_map_path = p.osrm_map_path;
	extended_model = p.extended_model;
	h_random_hospital = p.h_random_hospital;
	return *this;
}

std::ostream& operator<<(std::ostream& out, const Param &p){
	out << "Parameters:\n";
	if(p.debug){
		out << "DEBUG MODE\n";
	}
	out << "Instance: " << p.instance << "\n";
	out << "Solver: " << p.solver << "\n";
	out << "Generator File: " << p.generator_folder << "\n";
	if(p.solver == "heuristic"){
		if(p.h_use_fixed_bases)
			out << "Using fixed bases\n";

		if(p.h_discard)
			out << "Discarding calls that can't be attended\n";
	}
	return out;
}

Param::~Param(){
}