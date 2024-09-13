/// This defines the structure, some constants and functions for sampling

enum PriorType { UNIFORM_PRIOR, GAMMA_PRIOR, TRUNC_GAMMA_PRIOR, EXP_PRIOR, BETA_PRIOR, LOG_BETA_PRIOR};// Defines different types of prior

enum EventType { EXPOSE, INFECT, RECOVERY, REMOVAL, TERMINATE};  // Defines different types of event

enum PropType{ ADD_UNOBS, REM_UNOBS, ADD_IM, REM_IM, RESAMPLE_OBS, RESAMPLE_UNOBS, RESAMPLE_IM, ADD_HISTIMM, REM_HISTIMM};

enum SwapType{ NOSWAP, SWAP_OBS, SWAP_UNOBS};

enum MBPType{ MU_L_MBP, SH_L_MBP, MU_I_MBP, SH_I_MBP};

struct TransParam{         // Transition parameters
	double mu_L_value;       // The values for the transition parameters
	double sh_L_value;
	double mu_I_value;
	double sh_I_value;
	double mu_L_sample;      // The values used for sampling (this account for when phiG is not one)
	double sh_L_sample;
	double mu_I_sample;
	double sh_I_sample;	
};

struct Individual{         // Information about an individuals
	bool inf;                // Set to true if the individual is infected
	bool dead;               // Set to true if the individual died
	bool found;				 // Set to true if the corpse was recovered
    bool histimm;             // Set to true if the individual has immunity acquired from a previous outbreak
	double Exp_t;            // The exposure time
	double Inf_t;            // The infection time
	double Rec_t;            // The recovery time
	double latentP;          // The latent period
	double infectiousP;      // The infectious period
	
	double Li_gamma;         // The likelihood of the gamma distributions

	double Li_gamma_store;   

	double likelihood_gamma(); // Calculates the likelihood
	
	void sample(const double expt, const TransParam &tp);
	void sample_sim(const double expt, const TransParam &tp);
	void sample_backwards(const TransParam &tp);
	double sample_probability(const TransParam &tp);
};

struct Herd{
	int index;               // A number index for the herd
	
	vector <Individual> ind; // A list of individuals in the herds
	double T_cull;           // The time for culling the herd
	double T0;		         // The time of the first exposure(s)
	
	int beta1_param;	     // References the parameter for the expernal force of infection
	int R0_param;            // References the parameter for R0
	
	double nprop;
	vector <double> inf_sample; // This is used to sample the infection time
	
	
	double Li;               // The likelihood
	
	double Pr_index;         // The prior on the number of index cases
	
	int nmulti_propose, nmulti_accept;// Used to calculate proposal acceptance probabilities 
	int nT0_propose, nT0_accept;        

	double T0_joint_jump;    // Size of joint T0 event change
	int nT0_joint_propose, nT0_joint_accept;    

	void initialise(int index_);
	
	void simulate(const int ninfected);
	void load_data();
	void sample_exp_inf();
	void remove_exp_inf();
	double likelihood();
	double prior_index();
	void set_T0();
	void beta1_proposal();
	void R0_proposal();
	void p_rem_proposal();
	void p_histimm_proposal();
	void p_nd_proposal();
	void multi_proposal();
	void swap_proposal();
	void joint_T0_event_proposal();
	void T0_proposal();
	void add_unobs();
};


struct Param{
	string name;             // The name of the parameter
	
	double value;            // The parameter value
	PriorType prior_type;    // The type of the prior
	double prior_val1;       // Values used to specify the prior
	double prior_val2;
	double max_val;			 // Upper limit
	double upper_bound;
	
	double jump;             // The size of random walk jumps
	int npropose;            // The number of MH proposals and accpetance 
	int naccept;
	
	double jump_MBP;
	int npropose_MBP;        // The number of MH proposals and accpetance 
	int naccept_MBP;
	int nfail_MBP; 
};

struct Event{
	double t;                // The time of an event
	EventType type;          // The type of an event
};

const auto inf_sampler_bin = 100u;  // The number of bins used for the infection sampler

const auto LARGE = 1000000.0;       // A token large number
const auto TINY = 0.000000001;      // A token small number

const auto UNSET = 99999999.0;      // A token unset number

const auto map_out_prior = false;   // Set to true to just map out the prior