#ifndef _CG_H
#define _CG_H


#include "main.h"
#include "data.h"
#include "solver.h"

class EventSimulator;
class Solver;

int get_time(int time);

class CGCall{
public:
	CGCall(Data& data, GRBEnv& env, Solver& solver);
	~CGCall();

	Data& data;
	GRBEnv& env;
	GRBModel model;
	Solver& solver;


	int src_type;
	int src_location;
	int amb_type;
	int loc_base;
	int dst_hosp;

	int type_call0;
	int local_call0;
	int iter;
	int columns_added;
	int*** lambda;
	
	int** A0_ab;
	int** A0_ah;
	int*** A0_alb;
	int**** A0_calh;
	int***** A0_callh;
	int** C;

	std::set<int>*** L_tab;

	std::vector<std::vector<int>> closest_l1_to_l;

	int**** vh_tall;
	int** vc_tl;
	int*** vb_tal;


	GRBVar*** x0_abh;
	GRBVar*** x0_ahh;
	GRBVar**** x0_albh;
	GRBVar*** y0_ahb;

	GRBVar****** xt_cablh;
	GRBVar****** xt_cahlh;
	GRBVar******* xt_calblh;
	GRBVar**** yt_ahb;
	GRBVar*** Ct_cl;
	GRBVar*** At_ab;
	GRBVar**** At_alb;


	GRBConstr** con_base_t0;
	GRBConstr** con_hospital_t0;
	GRBConstr*** con_location_t0;
	GRBConstr** con_queue_t0;


	GRBConstr*** con_base;
	GRBConstr*** con_hospital;
	GRBConstr**** con_location;
	GRBConstr*** con_queue;
	GRBConstr**** con_amb_location;
	GRBConstr*** con_amb_base;
	GRBConstr** con_cap_base;


	double** beta0_ab;
	double** alpha0_ah;
	double*** psi0_alb;
	double** phi0_cl;
	double*** beta_tab;
	double*** alpha_tah;
	double**** psi_talb;
	double*** phi_tcl;
	double** ups_tb;
	double**** theta_talb;
	double*** gamma_tab;



	void add_bases_constraints_t0();
	void add_hospitals_constraints_t0();
	void add_locations_constraints_t0();
	void add_queues_constraints_t0();

	void add_bases_constraints();
	void add_hospitals_constraints();
	void add_locations_constraints();
	void add_queues_constraints();
	void add_ambs_location_constraints();
	void add_base_cap_constraints();


	int c_tl(int t, int l);
	int b_tal(int t, int a, int l);
	int h_tall(int t, int a, int l1, int l);

	void solve();

	void set_dual();
	void add_column(std::vector<int>& column);
	double sub_problem(int t, int a, int l, int l1, int c, int h, int b);
	bool pricing();

	// int is_at_hospital(Ambulance & amb);
	// int is_at_base(Ambulance & amb);
	// std::pair<int,int> is_at_location_base(Ambulance & amb);
	
};


class CGAmbulance{
public:
	CGAmbulance(Data & data, GRBEnv & env, Solver& solver);
	~CGAmbulance();

	Data& data;
	GRBEnv& env;
	GRBModel model;
	Solver& solver;

	int amb_type;
	int src_hosp;
	int base_return;
	
	int call_type;
	int call_location;
	int dst_hosp;

	int iter;
	int columns_added;
	
	int** A0_ab;
	int** A0_ah;
	int*** A0_alb;
	int**** A0_calh;
	int***** A0_callh;
	int** C;

	std::vector<std::vector<int>> closest_l1_to_l;
	int**** vh_tall;
	int** vc_tl;
	int*** vb_tal;

	GRBVar*** y0_ahb;
	GRBVar***** y0_ahclh;


	GRBVar****** xt_cablh;
	GRBVar****** xt_cahlh;
	GRBVar******* xt_calblh;
	GRBVar**** yt_ahb;
	GRBVar*** Ct_cl;
	GRBVar*** At_ab;
	GRBVar**** At_alb;

	GRBConstr** con_base_t0;
	GRBConstr** con_hospital_t0;
	GRBConstr*** con_location_t0;
	GRBConstr** con_queue_t0;

	GRBConstr*** con_base;
	GRBConstr*** con_hospital;
	GRBConstr**** con_location;
	GRBConstr*** con_queue;
	GRBConstr**** con_amb_location;
	GRBConstr*** con_amb_base;
	GRBConstr** con_cap_base;


	double** beta0_ab;
	double** alpha0_ah;
	double*** psi0_alb;
	double** phi0_cl;
	double*** beta_tab;
	double*** alpha_tah;
	double**** psi_talb;
	double*** phi_tcl;
	double** ups_tb;
	double**** theta_talb;
	double*** gamma_tab;


	void add_bases_constraints_t0();
	void add_hospitals_constraints_t0();
	void add_locations_constraints_t0();
	void add_queues_constraints_t0();


	void add_bases_constraints();
	void add_hospitals_constraints();
	void add_locations_constraints();
	void add_queues_constraints();
	void add_ambs_location_constraints();
	void add_base_cap_constraints();

	int c_tl(int t, int l);
	int b_tal(int t, int a, int l);
	int h_tall(int t, int a, int l1, int l);
	
	void solve();

	void set_dual();
	void add_column(std::vector<int>& column);
	double sub_problem(int t, int a, int l, int l1, int c, int h, int b);
	bool pricing();

	// int is_at_hospital(Ambulance & amb);
	// int is_at_base(Ambulance & amb);
	// std::pair<int,int> is_at_location_base(Ambulance & amb);
};


#endif