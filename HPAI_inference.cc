#include <iostream>
#include <sstream>
#include <math.h>
#include <map>
#include <algorithm>
#include <vector>
#include <iterator>
#include <fstream>
#include <iomanip>
#include <string>
#include <random>

using namespace std;

#include "struct.hh"

enum Model{ MULTI_INIT, EXT_INF};

enum Transmission{ FD, DD};

Model model = EXT_INF;						 // Chosen EXT_INF here which switches beta_1 on.  If MULTI_INIT is chosen with index_max=1 then this would give a single index case and no external transmission.
	
Transmission transmission = FD;				 // Can choose DD (for density-dependent transmission) or FD (for frequency-dependent transmission)

const auto phi_L = 1.0;                      // Inverse temperature multipling the infection process
const auto phi_G = phi_L;                    // Inverse temperature multipling the transition process

const auto ncolony = 1u;                       // The number of colonies

const auto alpha = 0.5;                      // Determines an exponentially distributed prior on index cases - irrelavant for ASF data as we will fix to have only one index case
const auto index_max = 1;                    // The maximum limit on the number of index cases - the above line only matters if this is > 1

const auto nupd = 40;			             // The multi_proposals on each iteration

const auto nsamp = 3000001u;                  // The number of MCMC samples
const auto nburnin = nsamp/10;                // The number step for burnin 
unsigned int s;                              // The sample number

const int NO_INF = -10000;                   // The log likelihoFod of getting exposed when no infected

bool simulated = false;				// Use simulated or observed removal times
//const vector <int> colony_sizes = { 1061};										// Colony size
const vector <int> colony_sizes = { 2071};										// Colony size
const vector <double> pct_seen = { 1.0};										// Proportion of the outbreak that is oberved (NOT RELEVANT FOR HPAI ANALYSES IN THE PAPER - LEAVE AT 1)

// Carcass collection times for each colony
//std::vector<std::vector<double>> remtimes = {{ 1.0, 2.0, 3.0, 4.0, 6.0, 6.0, 7.0, 8.0, 9.0, 9.0, 10.0, 11.0, 11.0, 11.0, 12.0, 12.0, 13.0, 13.0, 13.0, 14.0, 14.0, 14.0, 14.0, 14.0, 15.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 17.0, 17.0, 17.0, 18.0, 18.0, 19.0, 19.0, 19.0, 19.0, 19.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 27.0, 27.0, 27.0, 27.0, 27.0, 28.0, 28.0, 28.0, 28.0, 29.0, 29.0, 29.0, 30.0, 30.0, 31.0, 31.0, 32.0, 32.0, 32.0, 34.0, 34.0, 38.0}};
std::vector<std::vector<double>> remtimes = {{ 1.0, 2.0, 5.0, 5.0, 7.0, 8.0, 8.0, 9.0, 9.0, 10.0, 10.0, 10.0, 10.0, 10.0, 11.0, 11.0, 11.0, 11.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 21.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 22.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 26.0, 27.0, 27.0, 27.0, 27.0, 27.0, 27.0, 27.0, 27.0, 27.0, 27.0, 27.0, 27.0, 27.0, 27.0, 27.0, 27.0, 27.0, 27.0, 27.0, 27.0, 27.0, 27.0, 27.0, 27.0, 27.0, 27.0, 28.0, 28.0, 28.0, 28.0, 28.0, 28.0, 28.0, 28.0, 28.0, 28.0, 28.0, 28.0, 28.0, 28.0, 28.0, 28.0, 28.0, 28.0, 28.0, 28.0, 28.0, 28.0, 28.0, 28.0, 28.0, 28.0, 29.0, 29.0, 29.0, 29.0, 29.0, 29.0, 29.0, 29.0, 29.0, 29.0, 29.0, 29.0, 29.0, 29.0, 29.0, 29.0, 29.0, 29.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 31.0, 31.0, 31.0, 31.0, 31.0, 31.0, 31.0, 31.0, 31.0, 32.0, 32.0, 32.0, 32.0, 32.0, 32.0, 32.0, 32.0, 32.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 34.0, 34.0, 34.0, 34.0, 34.0, 34.0, 35.0, 35.0, 35.0, 35.0, 35.0, 35.0, 35.0, 35.0, 36.0, 37.0, 37.0, 37.0, 37.0, 38.0, 39.0, 40.0, 40.0, 41.0, 43.0, 43.0, 46.0, 46.0, 53.0}};

vector <Colony> colony(ncolony); // Stores information 
vector <Colony> colony_ppc(ncolony);	

vector <Param> param;      // A list of all the model parameters

int mu_L;                  // The parameter giving the mean of the latent period
int sh_L;                  // The parameter giving the standard deviation of the latent period
int mu_I;                  // The parameter giving the mean of the infectious period
int sh_I;                  // The parameter giving the standard deviation of the infectious period
int p_rem;                 // The parameter giving the probability of removal (rather than recovery)
int p_nd;                  // The parameter giving the probability that a dead bird is not detected
int p_histimm;              // The parameter giving the probability that a bird in the colony has immunity from a previous outbreak

double Pr;                 // The log of the prior

double time_total;
double time_likelihood = 0;

bool is_ppc = false;		// Switched between true and false for generating posterior predictive checks rather than proposing new times.

int slurm_id;				// Parameters below are for running jobs in parallel on HPC using slurm

std::string filename;
std::string filename2;
std::string filename3;

ofstream trace;
ofstream ppc;
ofstream colony_data;

default_random_engine generator;

bool EV_ord (Event ev1, Event ev2){ return (ev1.t < ev2.t); };  

int add_param(string name, PriorType prior_type, double prior_val1, double prior_val2, double max_val);
void mcmc_initialise();
void param_prior_sample();
void param_prior_init();
double prior();
void randomwalk_proposal();
void MBP_proposal(MBPType type);
double choose(int N, int n);
double choose_mult(int N, int n1, int n2);
void mcmc_diagnostics();
void check(int num);
		
void trace_init();
void colony_init();
void colony_write();
void ppc_init();
void trace_plot();
void ppc_plot();

TransParam set_trans_param();
double gamma_sample(const double mu, const double sh);
double trunc_gamma_sample(const double mu, const double sh, const string nm);
double beta_sample(const double a, const double b, const double c);
double gamma_probability(const double x, const double mu, const double sh);
double trunc_gamma_probability(const double x, const double mu, const double sh, const string nm);
double beta_probability(const double x, const double a, const double b);
double exp_sample(const double rate);
double exp_probability(const double x, const double mu);
double normal_sample(const double mu, const double sd);
double normal_probability(const double x, const double mean, const double var);
double uniform_sample(const double lwr, const double upr);
void emsg(const string& msg);
void choose_init();

int main(int argc, char** argv)
{

	if(argc != 2) emsg("Must provide a seed");
	
	// Initialise all the model parameters
	
	std::stringstream ss;					// Define filename based on job number
	std::stringstream ss2;
	std::stringstream ss3;

	// Initialise all the model parameters

	ss << "trace_hpai_obs_year1_fd_" << argv[1] << ".txt";			// Open output files
	filename = ss.str();
	ss2 << "colony_data_hpai_obs_year1_fd_" << argv[1] << ".txt";
	filename2 = ss2.str();
	ss3 << "ppc_hpai_obs_year1_fd_" << argv[1] << ".txt";
	filename3 = ss3.str();
	trace.open(filename);
	colony_data.open(filename2);
	ppc.open(filename3);
	
	slurm_id = std::stoi(argv[1]);			// Set the seed so that results can be replicated and each chain is different
	
	mu_L = add_param("mu_L",GAMMA_PRIOR,2.0,10.0,UNSET);	// Set up some parameters
	sh_L = add_param("sh_L",GAMMA_PRIOR,5.0,5.0,UNSET);
	
	mu_I = add_param("mu_I",GAMMA_PRIOR,5.0,10.0,UNSET);
	sh_I = add_param("sh_I",GAMMA_PRIOR,5.0,5.0,UNSET);

    p_rem = add_param("p_rem",BETA_PRIOR,2.0,6.0,0.5);
	p_nd = add_param("p_nd",BETA_PRIOR,4,200,0.1);
    p_histimm = add_param("p_histimm",BETA_PRIOR,11.0,40.0,0.5);
	
	generator.seed(14);
	
	for(auto h = 0u; h < ncolony; h++){		// Set up each colony
		colony[h].initialise(h);
		is_ppc = true;
		colony_ppc[h].initialise(h);
		is_ppc = false;
	}
	
	colony_init();							// Sets up an output file for data on removal times in each colony			
	
	if(simulated){                                                    // Simulates some data and writes to a file
		param_prior_init();            								  // Sets the parameter values
		for(auto &co : colony) co.simulate(int(pct_seen[co.index]*colony_sizes[co.index]));
		//for(auto &co : colony) co.simulate(int(1*colony_sizes[co.index]));
	} 
	else {                                                            // Loads some data
		for(auto &co : colony) co.load_data();
	}
	
	colony_write();
	generator.seed(slurm_id+10);
	
	param_prior_sample();											  // Samples values for each parameter from the priors  
	for(auto &co : colony){        									  // Reconstructs other events from recovery times
		if(simulated == true){
			co.remove_exp_inf();
		} else {		
			co.sample_exp_inf();                                          // Sets exposed and infection times for observed
		}
		co.set_T0();                                                  // Sets initial infection time
	}

	mcmc_initialise();												  // For outputing information about the chains

	for(auto &co : colony){											  // Outputs some info on the initial state straight away
		auto ninf = 0u, ndead = 0u, nfound = 0u;
		for(const auto &in : co.ind){
			if(in.inf == true) ninf++;
			if(in.dead == true) ndead++;
			if(in.found == true) nfound++;
		}
			
		cout << "Simulation of colony " << co.index << "     ";
		cout << "  Infected: " << ninf << "      ";
		cout << "  Dead: " << ndead << "      ";
		cout << "  Found: " << nfound << "      ";
		cout << "  End time: " << co.T_end << endl;
	
		if(true){
			for(const auto &in : co.ind){
				if(in.inf == true){
					cout << "  Exposed: " << in.Exp_t << " ";
					cout << "Infected: " << in.Inf_t << " ";
					cout << "Recovered: " << in.Rec_t << " ";
					cout << "LatentP: " << in.latentP << " ";
					cout << "InfectiousP: " << in.infectiousP << " ";
					if(in.dead == true) cout << "Dead";
					if(in.found == true) cout << "Found";
					cout << endl;
				}
			}
		}
		cout << endl;
	}

	check(0);

	// Perform mcmc

	trace_init();
	ppc_init();
	
	time_total = -clock();
	for(s = 0; s < nsamp; s++){
		if(s%1000 == 0){
			cout << "Sample " << s << endl;
			/* auto historical_immunity_check = 0;
			auto N_susceptible = 0;
			for(auto &co : colony){
				for(const auto &in : co.ind){
					if(in.inf == true){
						cout << "  Exposed: " << in.Exp_t << " ";
						cout << "Infected: " << in.Inf_t << " ";
						cout << "Recovered: " << in.Rec_t << " ";
						cout << "LatentP: " << in.latentP << " ";
						cout << "InfectiousP: " << in.infectiousP << " ";
						if(in.dead == true) cout << "Dead";
						if(in.found == true) cout << "Found";
						cout << endl;
					}
					if(in.histimm == true) historical_immunity_check++;
					if(in.histimm == false && in.inf == false) N_susceptible++;	
				}
			}
			cout << "N_hist_imm: " << historical_immunity_check << " ";
			cout << "N_sus: " << N_susceptible << " ";
			cout << endl; */
		}
		
		if(s%10000 == 0){
			for(auto &co : colony_ppc) co.simulate(int(1.0*colony_sizes[co.index]));
			ppc_plot();
		}
		
		if(s%100 == 0) trace_plot();
	
		for(auto &co : colony){
			if(model == EXT_INF) co.beta1_proposal();              // Proposals to change the external infection rate

			co.R0_proposal();                   // Proposals to change transmission rate

			co.p_rem_proposal();				// Proposals to change probability of dying

			co.p_histimm_proposal();			// Proposals to change probability of having historically acquired immunity

			co.p_nd_proposal();					// Proposals to change probability of mortality not being detected

			co.multi_proposal();                // Proposals which change exposure / infection times
	
			co.swap_proposal();                 // Proposals which change which individuals are indexes
			
			co.T0_proposal();                   // Proposals which change initial infection time
			
			co.joint_T0_event_proposal();       // Proposals which change T0 and other events
		}
		
		randomwalk_proposal();                                   // Proposals to change parameters for gamma distributed transitions
		
		check(1);                                                // Checks algorithm is performing correctly
	}
	time_total += clock();
	
	if(true){
		for(auto co : colony){
			for(auto i = 0u; i < inf_sampler_bin; i++){
				cout << co.inf_sample[i] << ",";
			}
			cout << "\n";
		}
	}
	
	mcmc_diagnostics();
}	


/// Initialisises mcmc
void mcmc_initialise()
{
	choose_init();
	
	Pr = prior();
	for(auto &co : colony){
		co.Li = co.likelihood();
		
		for(auto &in : co.ind) in.Li_gamma = in.likelihood_gamma();
		
		co.Pr_index = co.prior_index();
	
		co.inf_sample.resize(inf_sampler_bin);
		for(auto &bin : co.inf_sample) bin = 100;
		
		co.nprop = 10;
			
		co.nmulti_propose = 0; co.nmulti_accept = 0; 
		co.nT0_propose = 0; co.nT0_accept = 0;
		
		co.T0_joint_jump = 1;
		co.nT0_joint_propose = 0; co.nT0_joint_accept = 0;
	}
}

/// Initialises a colony with a transmission rate and number of individuals 
void Colony::initialise(int index_)
{
	index = index_;
	
	if(!is_ppc){	// There are duplicate colonies for simulating the PPCs.
		if(model == EXT_INF) beta1_param = add_param("beta1_"+to_string(index),UNIFORM_PRIOR,-15.0,0.0,UNSET);  // Lower limit of -15 because lower cases -inf values and numerical issues
		else beta1_param = UNSET;
	
		R0_param = add_param("R0_"+to_string(index),UNIFORM_PRIOR,1.0,30.0,UNSET);
	} else {
		beta1_param = param.size()-2;
		R0_param = param.size()-1;
	}
}

/// Simulates from a colony until a certain number of individuals recovers (or epidemic dies out)
void Colony::simulate(const int nrecover)
{
	ind.resize(colony_sizes[index]);

	for(auto &in : ind){
		in.inf = false;
		in.dead = false;
		in.found = false;
        in.histimm = false;
	}
	
	auto beta1 = 0.0; if(model == EXT_INF) beta1 = exp(param[beta1_param].value);
	auto beta2 = param[R0_param].value/param[mu_I].value; if(transmission == DD) beta2 = beta2/colony_sizes[index];

	auto tp = set_trans_param();

	// For simplicity, assume the epidemic starts at time T0
	T0 = 100.0;
	
	switch(model){
		case EXT_INF:
			ind[0].inf = true;
			ind[0].dead = true;
			ind[0].found = true;
			ind[0].sample_sim(T0,tp);
			break;
			
		case MULTI_INIT:
			auto N = ind.size();                      // Samples the number of index cases
			vector <double> prob_sum(N);
			prob_sum[0] = 0;
			auto sum = 0.0;
			for(auto n = 1u; n < N; n++){
				if(n <= index_max) sum += exp(-alpha*n);
				prob_sum[n] = sum;
			}
			
			auto z = uniform_sample(0.0,1.0)*sum;
			auto n = 0u; while(n < N && z > prob_sum[n]) n++;
			
			cout << "Colony " << index << " number of index cases: " << n << endl;
			
			for(auto j = 0u; j < n; j++){
				ind[j].inf = true;
				ind[j].found = true;
				ind[j].sample_sim(T0,tp);
				ind[j].dead = true;
			}
			break;
	}
	
	auto t = T0;
	
	vector <int> sus_ind;     // A list of susceptible individuals
	vector <int> inf_ind;     // A list of infected individuals  
	vector <int> histimm_ind;  // A list of historically immune individuals
		
	for(auto i = 0u; i < ind.size(); i++){
		if(ind[i].inf == false){
            if(uniform_sample(0.0,1.0) < param[p_histimm].value){
                histimm_ind.push_back(i);
            } else {
                sus_ind.push_back(i);
            }
        } else {
            inf_ind.push_back(i);
        }
	}
		
	auto I = 0u;                       // The number of infectious individuals
	auto N = ind.size();               // The number of alive individuals

	auto N_min = N - nrecover;         
	while(N > N_min){
		auto tnext = LARGE;              // This works out when the next event is
		int tnexti = UNSET;
		for(auto i : inf_ind){
			if(ind[i].inf == true){ 
				auto tt = ind[i].Inf_t;
				if(tt > t && tt < tnext){ tnext = tt; tnexti = i;}
			
				tt = ind[i].Rec_t;
				if(tt > t && tt < tnext){ tnext = tt; tnexti = i;}
			}
		}

		if(I == 0 && tnext == LARGE) break;
		
        auto rate = 0.0;
		if(transmission == DD){
			rate = sus_ind.size()*(beta1 + I*beta2);
		} else {
			rate = sus_ind.size()*(beta1 + I*beta2/N);
		}

		auto tnew = t + exp_sample(rate);
		
		if(tnew < tnext){                 // Expose a new individual
			t = tnew;
			
			auto j = int(uniform_sample(0.0,1.0)*sus_ind.size());
			auto i = sus_ind[j];
			sus_ind.erase(sus_ind.begin()+j);
			inf_ind.push_back(i);
			
			ind[i].inf = true;
			auto tempno1 = uniform_sample(0.0,1.0);
			auto tempno2 = uniform_sample(0.0,1.0);
			if(tempno1 > param[p_rem].value){
                ind[i].dead = false;
				ind[i].found = false;
			} else if(tempno2 < param[p_nd].value){
				ind[i].dead = true;
				ind[i].found = false;
			} else {
				ind[i].dead = true;
				ind[i].found = true;
			}
			
			ind[i].sample_sim(t,tp);
		}
		else{
			t = tnext;
			
			if(tnext == ind[tnexti].Inf_t){ // An individuals becomes infectious
				I++;	
			}
			else{                           // An individual recovers
				if(tnext != ind[tnexti].Rec_t) emsg("Problem 1");
				if(ind[tnexti].dead == true){
					I--; N--;
				} else {
					I--;
				}
			}
		}
	}
	
	T_end = t+0.00001;
}


/// This loads up the data into the colony
void Colony::load_data()
{
	ind.resize(colony_sizes[index]);
	
	for(auto &i : ind){
		i.inf = false;
		i.dead = false;
		i.found = false;
        i.histimm = false;
	}
	
	cout << remtimes[index].size() << " size\n";
	
	T_end = -LARGE;
	for(auto ind_index = 0u; ind_index < remtimes[index].size(); ind_index++){
		auto Rt = remtimes[index][ind_index]+100.0;//-uniform_sample(0.0,1.0);
		
		ind[ind_index].Rec_t = Rt;
		ind[ind_index].inf = true;
		ind[ind_index].dead = true;
		ind[ind_index].found = true;
        ind[ind_index].histimm = false;
		
		if(Rt > T_end) T_end = Rt; 
	}
	
	T_end += TINY;                                 // End happens just after the last recovery
}

/// Remove exposure and infection times and unobserved mortalities if data is simulated
void Colony::remove_exp_inf()
{
	auto tp = set_trans_param();
	vector <int> dead_ind;     // A list of dead individuals  
	for(auto i = 0u; i < ind.size(); i++){
		auto &in = ind[i]; 
		if(in.Exp_t == T0) in.sample_backwards(tp);
		dead_ind.push_back(i);
	}
	for(auto i = 0u; i < ind.size(); i++){
		auto &in = ind[i]; 
		if(in.dead == true && in.found == true){
			if(in.Exp_t != T0){
				auto infector = int(uniform_sample(0.0,1.0)*dead_ind.size());
				infector = dead_ind[infector];
				while(ind[infector].Rec_t > in.Rec_t || ind[infector].inf != true){
					infector = int(uniform_sample(0.0,1.0)*dead_ind.size());
					infector = dead_ind[infector];
				}
				in.Exp_t = (ind[infector].Inf_t + ind[infector].Rec_t)/2.0;
				in.Inf_t = (in.Exp_t + in.Rec_t)/2.0;
				in.latentP = in.Inf_t - in.Exp_t;
				in.infectiousP = in.Rec_t - in.Inf_t;
				dead_ind.push_back(i);
			}
		} else {
			in.inf = false;
			in.dead = false;
			in.found = false;
			//if(in.Exp_t < 0.0) emsg("Problem");
		}
	}
}

/// Adds exposure and infection times to the known removal times
void Colony::sample_exp_inf()
{
	auto tp = set_trans_param();
	if(model == EXT_INF){
		for(auto &in : ind){
			if(in.dead == true && in.found == true){
				in.sample_backwards(tp);
				if(in.Exp_t < 0.0) emsg("Problem");
			}
		}
	} else {
		for(auto i = 0u; i < ind.size(); i++){
			auto &in = ind[i]; 
			if(in.dead == true && in.found == true){
				if(i == 0){
					in.sample_backwards(tp);
				} else if(i == 1) {
					in.Exp_t = (ind[0].Inf_t + ind[0].Rec_t)/2.0;
					in.Inf_t = (in.Exp_t + in.Rec_t)/2.0;
					in.latentP = in.Inf_t - in.Exp_t;
					in.infectiousP = in.Rec_t - in.Inf_t;
				} else {
					int infector = rand() % (i-1);
					while(ind[infector].Rec_t > in.Rec_t || ind[infector].inf != true){
						infector = rand() % ind.size();
					}
					in.Exp_t = (ind[infector].Inf_t + ind[infector].Rec_t)/2.0;
					in.Inf_t = (in.Exp_t + in.Rec_t)/2.0;
					in.latentP = in.Inf_t - in.Exp_t;
					in.infectiousP = in.Rec_t - in.Inf_t;
				}
			}
		}
		for(auto &in : ind){
			if(in.dead == false || in.found == false){
				in.inf = false;
				in.dead = false;
				in.found = false;
				//if(in.Exp_t < 0.0) emsg("Problem");
			}
		}
	}
}

/// Sets the time of the initial infection(s)
void Colony::set_T0()
{
	T0 = LARGE;
	for(const auto &in : ind){
		if(in.inf == true){
			if(in.Exp_t < T0) T0 = in.Exp_t;
		}
	}
}

 
/// Performs a random walk MH proposal for beta1
void Colony::beta1_proposal()
{
	auto &par = param[beta1_param];

	auto param_store = par.value;
	
	par.value += normal_sample(0,par.jump);
	
	auto L_prop = likelihood();                 // Calculate the new likelihood and prior
	auto Pr_prop = prior();
	
	auto al = exp(phi_L*(L_prop  - Li) + Pr_prop - Pr);  // Calculates the Metropolis-Hastings acceptance probability

	par.npropose++;
	if(uniform_sample(0.0,1.0) < al){
		par.naccept++;
		Li = L_prop;
		Pr = Pr_prop;
		if(s < nburnin) par.jump *= 1.001;          // The dynamically adapts the jumping size
	}
	else{
		par.value = param_store;
		if(s < nburnin) par.jump *= 0.9995;
	}
}

/// Performs a random walk MH proposal for R0 parameter
void Colony::R0_proposal()
{
	auto &par = param[R0_param];

	auto param_store = par.value;
	
	par.value += normal_sample(0,par.jump);
	
	auto L_prop = likelihood();                 // Calculate the new likelihood and prior
	auto Pr_prop = prior();
	
	if(par.value < 0) L_prop = -LARGE;
	
	auto al = exp(phi_L*(L_prop - Li) + Pr_prop - Pr);  // Calculates the Metropolis-Hastings acceptance probability

	par.npropose++;
	if(uniform_sample(0.0,1.0) < al){
		par.naccept++;
		Li = L_prop;
		Pr = Pr_prop;
		if(s < nburnin) par.jump *= 1.001;          // The dynamically adapts the jumping size
	}
	else{
		par.value = param_store;
		if(s < nburnin) par.jump *= 0.9995;
	}
}

void Colony::p_rem_proposal()
{
	auto &par = param[p_rem];

	auto param_store = par.value;
	
	par.value += normal_sample(0,par.jump);
	
	auto L_prop = likelihood();                 // Calculate the new likelihood and prior
	auto Pr_prop = prior();
	
	if(par.value < 0) L_prop = -LARGE;
	if(par.value > par.max_val) L_prop = -LARGE;

	
	auto al = exp(phi_L*(L_prop - Li) + Pr_prop - Pr);  // Calculates the Metropolis-Hastings acceptance probability

	par.npropose++;
	if(uniform_sample(0.0,1.0) < al){
		par.naccept++;
		Li = L_prop;
		Pr = Pr_prop;
		if(s < nburnin) par.jump *= 1.001;          // The dynamically adapts the jumping size
	}
	else{
		par.value = param_store;
		if(s < nburnin) par.jump *= 0.9995;
	}
}

void Colony::p_histimm_proposal()
{
	auto &par = param[p_histimm];

	auto param_store = par.value;
	
	par.value += normal_sample(0,par.jump);
	
	auto L_prop = likelihood();                 // Calculate the new likelihood and prior
	auto Pr_prop = prior();
	
	if(par.value < 0) L_prop = -LARGE;
	if(par.value > par.max_val) L_prop = -LARGE;
	
	auto al = exp(phi_L*(L_prop - Li) + Pr_prop - Pr);  // Calculates the Metropolis-Hastings acceptance probability

	par.npropose++;
	if(uniform_sample(0.0,1.0) < al){
		par.naccept++;
		Li = L_prop;
		Pr = Pr_prop;
		if(s < nburnin) par.jump *= 1.001;          // The dynamically adapts the jumping size
	}
	else{
		par.value = param_store;
		if(s < nburnin) par.jump *= 0.9995;
	}
}

void Colony::p_nd_proposal()
{
	auto &par = param[p_nd];

	auto param_store = par.value;
	
	par.value += normal_sample(0,par.jump);
	
	auto L_prop = likelihood();                 // Calculate the new likelihood and prior
	auto Pr_prop = prior();
	
	if(par.value < 0) L_prop = -LARGE;
	if(par.value > par.max_val) L_prop = -LARGE;
	
	auto al = exp(phi_L*(L_prop - Li) + Pr_prop - Pr);  // Calculates the Metropolis-Hastings acceptance probability

	par.npropose++;
	if(uniform_sample(0.0,1.0) < al){
		par.naccept++;
		Li = L_prop;
		Pr = Pr_prop;
		if(s < nburnin) par.jump *= 1.001;          // The dynamically adapts the jumping size
	}
	else{
		par.value = param_store;
		if(s < nburnin) par.jump *= 0.9995;
	}
}

/// This performs random walk proposals for the four parameters associates with latent and infectious times
void randomwalk_proposal()
{
	for(auto loop = 0u; loop < 4; loop++){
		int p;
		switch(loop%4){
			case 0: p = mu_L; break;
			case 1: p = mu_I; break;
			case 2: p = sh_L; break;
			case 3: p = sh_I; break;
		}

		auto &par = param[p];

		auto colony_store = colony;
		auto param_store = par.value;
		
		par.value += normal_sample(0,par.jump);
		
		auto dLi = 0.0, dLi_gamma = 0.0;
		for(auto &co : colony){
			if(p == mu_I){
				auto Li_prop = co.likelihood();
				dLi += Li_prop-co.Li;
				co.Li = Li_prop;
			}
			
		    for(auto &in : co.ind){
			    auto Li_gamma_prop = in.likelihood_gamma();
			    dLi_gamma += Li_gamma_prop - in.Li_gamma;
			    in.Li_gamma = Li_gamma_prop;
            }
		}
		
		auto Pr_prop = prior();
			
		auto al = exp(phi_L*dLi + phi_G*dLi_gamma + Pr_prop - Pr); 
	
		par.npropose++;
		if(uniform_sample(0.0,1.0) < al){
			par.naccept++;
			Pr = Pr_prop;
			if(s < nburnin) par.jump *= 1.001;                // The dynamically adapts the jumping size
		}
		else{
			par.value = param_store;
			colony = colony_store;
			if(s < nburnin) par.jump *= 0.9995;
		}
	}
}



/// Performs proposals which change T0 
void Colony::T0_proposal()
{		
	auto tp = set_trans_param();

	vector <int> list_index;               

	auto time_max = LARGE;                                   // Finds the maximum time for T0
	auto histI = 0;
	for(auto i = 0u; i < ind.size(); i++){
		const auto &in = ind[i]; 
		if(in.inf == true){
			if(in.Exp_t == T0){
				list_index.push_back(i);
				if(in.Inf_t < time_max) time_max = in.Inf_t;
			} else {
				if(in.Exp_t < time_max) time_max = in.Exp_t;
			}
		}
		if(in.histimm == true) histI++;
	}

	auto var = tp.mu_L_sample*tp.mu_L_sample/(list_index.size()*tp.sh_L_sample);

	auto T0_store = T0;

	T0 += normal_sample(0,sqrt(var)); 

	if(T0 >= time_max){
		T0 = T0_store;
	}
	else{
		auto beta1 = 0.0; if(model == EXT_INF) beta1 = exp(param[beta1_param].value);
		auto dLi = beta1*(T0-T0_store)*(ind.size()-histI-list_index.size()); // DE edit 05 Nov to make this affect everything other than index individuals 
		
		auto dLi_gamma = 0.0;
		for(auto i : list_index){
			auto &in = ind[i];
			in.Exp_t = T0;	
			in.latentP = in.Inf_t-in.Exp_t;
			
			in.Li_gamma_store = in.Li_gamma;
			in.Li_gamma = in.likelihood_gamma(); 
			dLi_gamma += in.Li_gamma - in.Li_gamma_store;
		}
		
		
		auto al = exp(phi_L*dLi + phi_G*dLi_gamma);  // Calculates the Metropolis-Hastings acceptance probability

		nT0_propose++;
		if(uniform_sample(0.0,1.0) < al){
			Li += dLi;
			nT0_accept++;
		}
		else{
			T0 = T0_store;
			for(auto i : list_index){
				auto &in = ind[i];
				in.Exp_t = T0;
				in.latentP = in.Inf_t-in.Exp_t;
				in.Li_gamma = in.Li_gamma_store;
			}
		}
	}
}


/// The performs proposals which swap individuals between being index / non-index cases
void Colony::swap_proposal()
{
	// First we make a list 
	vector <Event> event;
	auto num_index = 0u;
    auto histI = 0;
	for(const auto &in : ind){
		if(in.inf == true){
			if(in.Exp_t == T0) num_index++;
					
			if(in.Inf_t <= T_end){
				Event ev; ev.type = INFECT; ev.t = in.Inf_t;
				event.push_back(ev);
			}
			
			if(in.Rec_t <= T_end && in.dead == false){
				Event ev; ev.type = RECOVERY; ev.t = in.Rec_t;
				event.push_back(ev);
			}

			if(in.Rec_t <= T_end && in.dead == true){
				Event ev; ev.type = REMOVAL; ev.t = in.Rec_t;
				event.push_back(ev);
			}
		}
        if(in.histimm == true){
            histI++;
        }
	}
	
	sort(event.begin(),event.end(),EV_ord);
	Event ev; ev.type = TERMINATE; ev.t = T_end;
	event.push_back(ev);
	
	auto nevent = event.size();

	auto I = 0;             // The number of infected
	auto S = ind.size() - histI;    // The number of susceptible
	auto N = S + histI;             // The number alive
	vector <int> Nstore(nevent);
	vector <int> Istore(nevent);
	
	for(auto e = 0u; e < nevent; e++){
		Nstore[e] = N; Istore[e] = I;
		switch(event[e].type){
		case EXPOSE: break;
		case INFECT: I++; break;
		case RECOVERY: I--; break;
		case REMOVAL: I--; N--; break;
		case TERMINATE: break;
		}
	}
	
	auto tp = set_trans_param();
	
	auto beta1 = 0.0; if(model == EXT_INF) beta1 = exp(param[beta1_param].value);
	auto beta2 = param[R0_param].value/param[mu_I].value; if(transmission == DD) beta2 = beta2/ind.size();
	
	auto Exp_t_max = T0 + 2*tp.mu_L_value;     // This sets the limit for the exposure time
	
	for(auto loop = 0u; loop < 3; loop++){
		for(auto &in : ind){
			if(in.inf == true){
				if(in.Exp_t > T0){    // Add an index case
					if(in.Exp_t < Exp_t_max && num_index < index_max){
						auto t = in.Exp_t;
					
						auto dLi = 0.0;
					
						auto e = 0u; 
						auto tt = T0;
                        auto rate = 0.0;
						while(e < nevent && t > event[e].t){
							auto tnew = event[e].t;
							if(transmission == DD){
								auto rate = beta1 + beta2*Istore[e];
							} else {
								auto rate = beta1 + beta2*Istore[e]/Nstore[e];
							}
							dLi += (tnew-tt)*rate; tt = tnew;
							e++;
						}
						if(e == nevent) emsg("prob");
						if(transmission == DD){
							auto rate = beta1 + beta2*Istore[e];
						} else {
							auto rate = beta1 + beta2*Istore[e]/Nstore[e];
						}
						dLi += (t-tt)*rate - choose(ind.size()-histI,num_index+1) + choose(ind.size()-histI,num_index);
						if(rate > 0) dLi -= log(rate); else dLi -= NO_INF;
						
						auto dPr = -alpha;
						
						auto probif = 0;
						auto probfi = gamma_probability(in.latentP,tp.mu_L_sample,tp.sh_L_sample);
						
						auto latP_st = in.latentP;
						in.latentP = in.Inf_t-T0;
						auto Li_gamma_prop = in.likelihood_gamma();
						
						auto al = exp(phi_L*dLi + phi_G*(Li_gamma_prop - in.Li_gamma) + dPr + probfi - probif);
						if(uniform_sample(0.0,1.0) < al){
							Li += dLi;
							in.Li_gamma = Li_gamma_prop;
							Pr_index += dPr;
							in.Exp_t = T0;
							num_index++;
						}
						else{
							in.latentP = latP_st;
						}
					}
				} else {                 // Remove an index
					if(num_index > 1){
						auto latP = gamma_sample(tp.mu_L_sample,tp.sh_L_sample);
						auto t = in.Inf_t-latP;
						if(t > T0 && t < T_end && t < Exp_t_max){
							auto dLi = 0.0;
					
							auto e = 0u; 
							auto tt = T0;
                            auto rate = 0.0;
							while(e < nevent && t > event[e].t){
								auto tnew = event[e].t;
								if(transmission == DD){
									auto rate = beta1 + beta2*Istore[e]; 
								} else {
									auto rate = beta1 + beta2*Istore[e]/Nstore[e];
								}
								dLi -= (tnew-tt)*rate; tt = tnew;
								e++;
							}
							if(e == nevent) emsg("prob");
							if(transmission == DD){
								auto rate = beta1 + beta2*Istore[e];
							} else {
								auto rate = beta1 + beta2*Istore[e]/Nstore[e];
							}
							dLi += -(t-tt)*rate   - choose(ind.size()-histI,num_index-1) + choose(ind.size()-histI,num_index);
							if(rate > 0) dLi += log(rate); else dLi += NO_INF;
							
							
							auto dPr = alpha;

							auto latP_st = in.latentP;
							in.latentP = latP;
							auto Li_gamma_prop = in.likelihood_gamma();
							
							auto probif = gamma_probability(latP,tp.mu_L_sample,tp.sh_L_sample);
							auto probfi = 0;
							
							auto al = exp(phi_L*dLi + phi_G*(Li_gamma_prop-in.Li_gamma) + dPr + probfi - probif);
							if(uniform_sample(0.0,1.0) < al){
								Li += dLi;
								in.Li_gamma = Li_gamma_prop;
								Pr_index += dPr;
								
								in.Exp_t = t;		
								num_index--;
							}
							else{
								in.latentP = latP_st;
							}
						}	
					}
				}
			}
		}
	}
}

/// Makes multiple changes to events (the initial infection time is not changed)
void Colony::multi_proposal()
{
	const double prob_add_rem_unobs = 0.05;                  // This gives the probability of doing an add/rem proposal on unobserved mortality
	const double prob_add_rem_im = 0.3;                     // This gives the probability of doing an add/rem proposal on immune individual  
    const double prob_add_rem_histim = 0.3;
	const double prob_obs = 0.3;							// Prob of changing events on observed individual
	const double prob_unobs = 1-prob_add_rem_unobs-prob_add_rem_im-prob_add_rem_histim-prob_obs; 		// Prob of changing events on unobserved individual
	
	if(prob_unobs < 0) emsg("Samping error");
	
	// This generates a sampler for adding new exposure times
	vector <double> prob(inf_sampler_bin), prob_sum(inf_sampler_bin);

	auto sum = 0.0; for(auto i = 0u; i < inf_sampler_bin; i++) sum += inf_sample[i];
	
	auto sum2 = 0.0;
	for(auto i = 0u; i < inf_sampler_bin; i++){
		prob[i] = inf_sample[i]/sum;
		sum2 += prob[i];
		prob_sum[i] = sum2;
	}

	auto tp = set_trans_param();
	
	for(auto loop = 0u; loop < nupd; loop++){
		auto probif = 0.0, probfi = 0.0;
		
		// Makes lists of observed, infected (unobseved), susceptible and historically immune 
		vector <int> list_obs, list_unobs, list_im, list_sus, list_histimm;
		auto num_index = 0u;
		for(auto i = 0u; i < ind.size(); i++){
			if(ind[i].dead == false){
				if(ind[i].inf == false){
                    if(ind[i].histimm == true){
                        list_histimm.push_back(i);
                    } else {
                    	list_sus.push_back(i);
					}
                } else {
                    list_im.push_back(i);
                }
			} else if(ind[i].found == false){
				list_unobs.push_back(i);
			} else {
				list_obs.push_back(i);
			}
			
			if(ind[i].inf == true && ind[i].Exp_t == T0) num_index++;
		}
		
		auto nlist_im = list_im.size();
		auto nlist_sus = list_sus.size();
		auto nlist_obs = list_obs.size();
		auto nlist_unobs = list_unobs.size();
        auto nlist_histimm = list_histimm.size();
	
		vector <Individual> ind_store = ind;
	
		auto jmax = (unsigned int)(nprop+0.5); if(jmax == 0) jmax = 1;

		auto dLi_gamma = 0.0;
		for(auto j = 0u; j < jmax; j++){
			
			PropType type;
			
			auto z = uniform_sample(0.0,1.0);                                // Randomly select what change to make
			if(z < prob_add_rem_unobs){                          // Adds or removes an individual
				if(uniform_sample(0.0,1.0) < 0.5) type = ADD_UNOBS;
				else type = REM_UNOBS;
			} else if(z < (prob_add_rem_unobs+prob_add_rem_im)){
				if(uniform_sample(0.0,1.0) < 0.5) type = ADD_IM;
				else type = REM_IM;
			} else if(z < (prob_add_rem_unobs+prob_add_rem_im+prob_obs)){
				type = RESAMPLE_OBS;        // Resamples times for observed
            } else if(z < (prob_add_rem_unobs+prob_add_rem_im+prob_obs+prob_add_rem_histim)){
                if(uniform_sample(0.0,1.0) < 0.5) type = ADD_HISTIMM;
				else type = REM_HISTIMM;
			} else if(z < (prob_add_rem_unobs+prob_add_rem_im+prob_obs+prob_add_rem_histim+prob_unobs/2.0)){ 
				type = RESAMPLE_UNOBS;                  // Resamples  times for unobserved
			} else {
				type = RESAMPLE_IM;
			}
		
			switch(type){
				case ADD_UNOBS:                              
					if(nlist_sus > 0){
						auto sus_sel = int(uniform_sample(0.0,1.0)*nlist_sus);     // Randomly selects a susceptible individual
						auto i = list_sus[sus_sel];         
						
						auto z = uniform_sample(0.0,1.0);                          // Selects which bin to sample from
						auto bin = 0u; while(bin < inf_sampler_bin && z > prob_sum[bin]) bin++;
						if(bin == inf_sampler_bin) emsg("Problem 2");
						
						auto t_exp = (bin+uniform_sample(0.0,1.0))*T_end/inf_sampler_bin; // Samples exposure time 
						if(t_exp > T0){
							auto &in = ind[i];
							in.sample(t_exp,tp);
							
							if(in.Rec_t < T_end){               // Only accept if recovery time before end
								in.inf = true;                     
								in.dead = true;
								in.found = false;

								probif += in.sample_probability(tp);
								probif += log((1.0/nlist_sus)*prob[bin]*(inf_sampler_bin/T_end));
								probfi += log((1.0/(nlist_unobs+1)));
								
								list_unobs.push_back(i); nlist_unobs++;	   // Updates the lists
								list_sus.erase(list_sus.begin()+sus_sel); nlist_sus--;	
								
								auto Li_gamma_prop = in.likelihood_gamma();
								dLi_gamma += Li_gamma_prop - in.Li_gamma;
								in.Li_gamma = Li_gamma_prop;
							}	
						}
					}
					break;
					
				case REM_UNOBS:                              
					if(nlist_unobs > 0){
						auto inf_sel = int(uniform_sample(0.0,1.0)*nlist_unobs);     // Randomly selects a infected individual
						auto i = list_unobs[inf_sel];       
						
						auto t_exp = ind[i].Exp_t;
                        if(t_exp > T0){                           // Only accept if not an index case
							auto &in = ind[i];

							in.inf = false;
							in.dead = false;
							in.found = false;

							auto bin = int(inf_sampler_bin*t_exp/T_end);
														
							probif += log((1.0/nlist_unobs));
							probfi += log((1.0/(nlist_sus+1))*prob[bin]*(inf_sampler_bin/T_end));
							probfi += in.sample_probability(tp);
								
							list_sus.push_back(i); nlist_sus++;	
							list_unobs.erase(list_unobs.begin()+inf_sel); nlist_unobs--;	
							
							auto Li_gamma_prop = in.likelihood_gamma();
							dLi_gamma += Li_gamma_prop - in.Li_gamma;
							in.Li_gamma = Li_gamma_prop;
						}
					}
					break;

                case ADD_HISTIMM:                             
					if(nlist_sus > 0){
						auto sus_sel = int(uniform_sample(0.0,1.0)*nlist_sus);     // Randomly selects a susceptible individual
						auto i = list_sus[sus_sel];         
						
						auto &in = ind[i];
						
                        in.inf = false;                     
						in.dead = false;
						in.found = false;
                        in.histimm = true;

						probif += log(1.0/nlist_sus);
						probfi += log(1.0/(nlist_histimm+1));
								
						list_histimm.push_back(i); nlist_histimm++;	   // Updates the lists
						list_sus.erase(list_sus.begin()+sus_sel); nlist_sus--;	

						auto Li_gamma_prop = in.likelihood_gamma();
						dLi_gamma += Li_gamma_prop - in.Li_gamma;
						in.Li_gamma = Li_gamma_prop;
					}
					break;  

                case REM_HISTIMM:                           
					if(nlist_histimm > 0){
						auto histimm_sel = int(uniform_sample(0.0,1.0)*nlist_histimm);     // Randomly selects a infected individual
						auto i = list_histimm[histimm_sel];       
						
						auto &in = ind[i];

						in.inf = false;
						in.dead = false;
						in.found = false;
                        in.histimm = false;

						probif += log(1.0/nlist_histimm);
						probfi += log(1.0/(nlist_sus+1));
								
						list_sus.push_back(i); nlist_sus++;	
						list_histimm.erase(list_histimm.begin()+histimm_sel); nlist_histimm--;	

						auto Li_gamma_prop = in.likelihood_gamma();
						dLi_gamma += Li_gamma_prop - in.Li_gamma;
						in.Li_gamma = Li_gamma_prop;
					}
					break;  
					
				case ADD_IM:                          
					if(nlist_sus > 0){
						auto sus_sel = int(uniform_sample(0.0,1.0)*nlist_sus);     // Randomly selects a susceptible individual
						auto i = list_sus[sus_sel];         
						
						auto z = uniform_sample(0.0,1.0);                          // Selects which bin to sample from
						auto bin = 0u; while(bin < inf_sampler_bin && z > prob_sum[bin]) bin++;
						if(bin == inf_sampler_bin) emsg("Problem 2");
						
						auto t_exp = (bin+uniform_sample(0.0,1.0))*T_end/inf_sampler_bin; // Samples exposure time 
						if(t_exp > T0){
							auto &in = ind[i];
							in.sample(t_exp,tp);
							
							if(in.Rec_t < T_end){               // Only accept if recovery time before end
								in.inf = true;                     
								in.dead = false;
								in.found = false;

								probif += in.sample_probability(tp);
								probif += log((1.0/nlist_sus)*prob[bin]*(inf_sampler_bin/T_end));
								probfi += log((1.0/(nlist_im+1)));
								
								list_im.push_back(i); nlist_im++;	   // Updates the lists
								list_sus.erase(list_sus.begin()+sus_sel); nlist_sus--;	
								
								auto Li_gamma_prop = in.likelihood_gamma();
								dLi_gamma += Li_gamma_prop - in.Li_gamma;
								in.Li_gamma = Li_gamma_prop;
							}	
						}
					}
					break;
					
				case REM_IM:                     
					if(nlist_im > 0){
						auto inf_sel = int(uniform_sample(0.0,1.0)*nlist_im);     // Randomly selects a infected individual
						auto i = list_im[inf_sel];       
						
						auto t_exp = ind[i].Exp_t;
                        if(t_exp > T0){                           // Only accept if not an index case
							auto &in = ind[i];

							in.inf = false;
							in.dead = false;
							in.found = false;

							auto bin = int(inf_sampler_bin*t_exp/T_end);
														
							probif += log((1.0/nlist_im));
							probfi += log((1.0/(nlist_sus+1))*prob[bin]*(inf_sampler_bin/T_end));
							probfi += in.sample_probability(tp);
								
							list_sus.push_back(i); nlist_sus++;	
							list_im.erase(list_im.begin()+inf_sel); nlist_im--;	
							
							auto Li_gamma_prop = in.likelihood_gamma();
							dLi_gamma += Li_gamma_prop - in.Li_gamma;
							in.Li_gamma = Li_gamma_prop;
						}
					}
					break;	
				
				case RESAMPLE_OBS:                        
					if(nlist_obs > 0){
						auto &in = ind[list_obs[int(uniform_sample(0.0,1.0)*nlist_obs)]];
						if(in.Exp_t == T0){                     // Resamples infection time for index case
							auto mu1 = in.Exp_t + tp.mu_L_sample;
							auto var1 = tp.mu_L_sample*tp.mu_L_sample/tp.sh_L_sample;
							
							auto mu2 = in.Rec_t - tp.mu_I_sample;
							auto var2 = tp.mu_I_sample*tp.mu_I_sample/tp.sh_I_sample;
							
							auto mu = (mu1*var2 + mu2*var1)/(var1+var2);
							auto var = var1*var2/(var1+var2);
						
							auto Inf_t = normal_sample(mu,sqrt(var));
							if(Inf_t > in.Exp_t && Inf_t < in.Rec_t){
								probfi += normal_probability(in.Inf_t,mu,var);
								probif += normal_probability(Inf_t,mu,var);
								in.Inf_t = Inf_t;
								in.latentP = in.Inf_t - in.Exp_t;
								in.infectiousP = in.Rec_t - in.Inf_t;	
								auto Li_gamma_prop = in.likelihood_gamma();
								dLi_gamma += Li_gamma_prop - in.Li_gamma;
								in.Li_gamma = Li_gamma_prop;								
							}
						}
						else{
							Individual ind_store = in;
							auto dprobfi = in.sample_probability(tp);
							in.sample_backwards(tp);
							if(in.Exp_t <= T0) in = ind_store;		
							else{
								probfi += dprobfi;
								probif += in.sample_probability(tp); 
								auto Li_gamma_prop = in.likelihood_gamma();
								dLi_gamma += Li_gamma_prop - in.Li_gamma;
								in.Li_gamma = Li_gamma_prop;
							}
						}
					}
					break;
					
				case RESAMPLE_UNOBS:                     
					if(nlist_unobs > 0){
						auto &in = ind[list_unobs[int(uniform_sample(0.0,1.0)*nlist_unobs)]];
						auto dprobfi = in.sample_probability(tp);
							
						Individual ind_store = in;
						if(uniform_sample(0.0,1.0) < 0.5){
							if(in.Exp_t != T0){	
								in.sample_backwards(tp);
								if(in.Exp_t >= T_end || in.Exp_t <= T0) in = ind_store; 
								else{
									probfi += dprobfi;
									probif += in.sample_probability(tp); 
									auto Li_gamma_prop = in.likelihood_gamma();
									dLi_gamma += Li_gamma_prop - in.Li_gamma;
									in.Li_gamma = Li_gamma_prop;
								}
							}
						}
						else{
							in.sample(in.Exp_t,tp);
							if(in.Rec_t > T_end) in = ind_store;
							else{
								probfi += dprobfi;
								probif += in.sample_probability(tp); 
								auto Li_gamma_prop = in.likelihood_gamma();
								dLi_gamma += Li_gamma_prop - in.Li_gamma;
								in.Li_gamma = Li_gamma_prop;
							}
						}
					}
					break;
					
				case RESAMPLE_IM:                     
					if(nlist_im > 0){
						auto &in = ind[list_im[int(uniform_sample(0.0,1.0)*nlist_im)]];
						auto dprobfi = in.sample_probability(tp);
							
						Individual ind_store = in;
						if(uniform_sample(0.0,1.0) < 0.5){
							if(in.Exp_t != T0){	
								in.sample_backwards(tp);
								if(in.Exp_t >= T_end || in.Exp_t <= T0) in = ind_store; 
								else{
									probfi += dprobfi;
									probif += in.sample_probability(tp); 
									auto Li_gamma_prop = in.likelihood_gamma();
									dLi_gamma += Li_gamma_prop - in.Li_gamma;
									in.Li_gamma = Li_gamma_prop;
								}
							}
						}
						else{
							in.sample(in.Exp_t,tp);
							if(in.Rec_t > T_end) in = ind_store;
							else{
								probfi += dprobfi;
								probif += in.sample_probability(tp); 
								auto Li_gamma_prop = in.likelihood_gamma();
								dLi_gamma += Li_gamma_prop - in.Li_gamma;
								in.Li_gamma = Li_gamma_prop;
							}
						}
					}
					break;
			}
		}

		auto Li_prop = likelihood();  
		
		// Calculates the Metropolis-Hastings acceptance probability
		auto al = exp(phi_L*(Li_prop - Li) + phi_G*dLi_gamma + probfi - probif);  

		nmulti_propose++;
		if(uniform_sample(0.0,1.0) < al){
			nmulti_accept++;
			Li = Li_prop;
			if(s < nburnin){ nprop *= 1.001; if(nprop > 100) nprop = 100;}
		}
		else{
			ind = ind_store;
			if(s < nburnin){ nprop *= 0.999; if(nprop < 1) nprop = 1;}
		}
	}	
	
	if(s < nburnin){                                    // This dynamically adapts the exposure time sampler
		for(auto i = 0u; i < ind.size(); i++){
			auto t_exp = ind[i].Exp_t;
			if((ind[i].inf == true && ind[i].dead == false) || (ind[i].inf == true && ind[i].dead == true && ind[i].found == false) && t_exp > 0){
				auto bin = int(inf_sampler_bin*t_exp/T_end);
				if(bin < 0 || bin >= inf_sampler_bin) emsg("Problem 3");
				inf_sample[bin]++;
			}
		}
	}
}

/// Sets the parameter values
void param_prior_init()
{
	for(auto &par : param){
		if(par.name == "mu_I"){
			par.value = 5.0;
		} else if (par.name == "mu_L"){
			par.value = 2.0;
		} else if (par.name == "sh_I"){
			par.value = 5.0;
		} else if (par.name == "sh_L"){
			par.value = 5.0;
		} else if (par.name == "beta1_0"){
			par.value = log(0.01/colony_sizes[0]);
		} else if (par.name == "R0_0"){
			par.value = 5.0;
		} else if (par.name == "p_rem"){
            par.value = 0.15;
        } else if (par.name == "p_nd"){
			par.value = 0.02;
		} else if (par.name == "p_histimm"){
            par.value = 0.1;
        }
		cout << "Sampled value for parameter " << par.name << ": " << par.value << "\n";			
	}
}

/// Samples the parameters from the model prior
void param_prior_sample()
{
	for(auto &par : param){
		switch(par.prior_type){
		case UNIFORM_PRIOR:
			par.value = par.prior_val1 + uniform_sample(0.0,1.0)*(par.prior_val2-par.prior_val1);
			break;
			
		case BETA_PRIOR:
			par.value = beta_sample(par.prior_val1,par.prior_val2,par.max_val);
			break;
			
		case GAMMA_PRIOR: 
			par.value = gamma_sample(par.prior_val1,par.prior_val2);
			break;
			
		case TRUNC_GAMMA_PRIOR:
			par.value = trunc_gamma_sample(par.prior_val1,par.prior_val2,par.name);
			break;
			
		case EXP_PRIOR:
			par.value = exp_sample(par.prior_val1);
			break;
		}
		cout << "Sampled value for parameter " << par.name << ": " << par.value << "\n";			
	}
}

/// Calculates the prior
double prior()
{
	auto Pr = 0.0;
	for(auto &par : param){
		switch(par.prior_type){
			case UNIFORM_PRIOR:
				if(par.value < par.prior_val1 || par.value > par.prior_val2) Pr -= LARGE; 
				else Pr += log(1.0/(par.prior_val2-par.prior_val1));
				break;
				
			case GAMMA_PRIOR:
				Pr += gamma_probability(par.value, par.prior_val1, par.prior_val2); 
				break;
				
			case TRUNC_GAMMA_PRIOR:
				Pr += trunc_gamma_probability(par.value, par.prior_val1, par.prior_val2, par.name);
				break;
				
			case EXP_PRIOR: 
				if(par.value < 0.0 || par.value > par.max_val) Pr -= LARGE;
				Pr += exp_probability(par.value, par.prior_val1); 
				break;	
				
			case BETA_PRIOR: 
				if(par.value < 0.0 || par.value > par.max_val) Pr -= LARGE;
				else Pr += beta_probability(par.value, par.prior_val1, par.prior_val2);
				break;
		}
	}
	
	return Pr;
}
 

/// Calculate the log likelihood associated with exposure events
double Colony::likelihood()
{
	if(map_out_prior == true) return 0;
	time_likelihood -= clock();
	
	// Gererates a sorted list of event for the entire colony
	
	vector <Event> event;
    auto histI = 0;
	auto Nexp = 0;
	auto Ninf = 0;
	auto Nrem = 0;
	auto Nrec = 0;
	for(const auto &in : ind){
		if(in.inf == true){
			Event ev; ev.type = EXPOSE; ev.t = in.Exp_t;
			event.push_back(ev);
			Nexp++;
			
			if(in.Inf_t <= T_end){
				Event ev; ev.type = INFECT; ev.t = in.Inf_t;
				event.push_back(ev);
				Ninf++;
			}
			
			if((in.Rec_t <= T_end) & (in.dead == false)){
				Event ev; ev.type = RECOVERY; ev.t = in.Rec_t;
				event.push_back(ev);
				Nrec++;
			}

            if((in.Rec_t <= T_end) & (in.dead == true)){
				Event ev; ev.type = REMOVAL; ev.t = in.Rec_t;
				event.push_back(ev);
				Nrem++;
			}
		}
        if(in.histimm == true){
            histI++;
        }
	}
	
	Event ev; ev.type = TERMINATE; ev.t = T_end;
	event.push_back(ev);
	
	sort(event.begin(),event.end(),EV_ord);

	if(T0 != event[0].t){ cout << T0 << " " << event[0].t << " j\n"; emsg(" Initial event times do not agree");}
	
	// Calculates the log likelihood
	
	auto beta1 = 0.0; if(model == EXT_INF) beta1 = exp(param[beta1_param].value);
	auto beta2 = param[R0_param].value/param[mu_I].value; if(transmission == DD) beta2 = beta2/ind.size();
		
	auto L = 0.0;           // The log of the likelihood
	
	auto t = 0.0;           // The start time
	
	auto I = 0;             // The number of infected
	auto S = ind.size() - histI;    // The number of susceptible
	auto N = S + histI;             // The number alive
    auto Rec = 0;           // The number recovered
    auto Rem = 0;           // The number removed
	
	auto num_index = 0u;  
	for(auto i = 0u; i < event.size(); i++){
		if(N == 0) break;
		
		auto tt = event[i].t;

		auto rate = 0.0;
		if(i != 0){
			if (transmission == DD){
				rate = beta1 + beta2*I;
			} else {
				rate = beta1 + beta2*I/N;
			}
		}

		L -= (tt-t)*rate*S;
			
		if(tt == T0) num_index++;
		else{
			if(event[i].type == EXPOSE){
				if(rate == 0.0) L += NO_INF; // Make an exposure with no infected individuals very unlikely
				else L += log(rate);
			}					
		}

		switch(event[i].type){
			case EXPOSE: S--; break;
			case INFECT: I++; break;
            case RECOVERY: I--; Rec++; break;
			case REMOVAL: 
					I--; N--; Rem++;  
					if(I < 0){
						ppc_plot();
						trace_plot();
						cout << "Sinit = " << ind.size() - histI << "  ";
						cout << "Ninit = " << ind.size() << "  ";
						cout << "Nevents = " << event.size() << "  ";
						cout << "This event = " << i << "  ";
						cout << "Nexp = " << Nexp << "  ";
						cout << "Ninf = " << Ninf << "  ";
						cout << "Nrec = " << Nrec << "  ";
						cout << "Nrem = " << Nrem << "  ";
						cout << "N = " << N << "  ";
						cout << "S = " << S << "  ";
						cout << "I = " << I << "  ";
						cout << "Rec = " << Rec << "  ";
						cout << "Rem = " << Rem << "  ";
						cout << "histImm = " << histI << "  ";
						cout << "T0 = " << T0 << "  ";
						auto historical_immunity_check = 0;
						auto N_susceptible = 0;
						for(const auto &in : ind){
							if(in.inf == true){
								cout << "  Exposed: " << in.Exp_t << " ";
								cout << "Infected: " << in.Inf_t << " ";
								cout << "Recovered: " << in.Rec_t << " ";
								cout << "LatentP: " << in.latentP << " ";
								cout << "InfectiousP: " << in.infectiousP << " ";
								if(in.dead == true) cout << "Dead";
								if(in.found == true) cout << "Found";
								cout << endl;
							}
							if(in.histimm == true) historical_immunity_check++;
							if(in.histimm == false && in.inf == false) N_susceptible++;	
						}
						cout << "N_hist_imm: " << historical_immunity_check << " ";
						cout << "N_sus: " << N_susceptible << " ";
						cout << endl;
						emsg("Negative 1");
					}
					if(N < 0) emsg("Negative 2"); break;
			case TERMINATE: break;
		}
		t = tt;
	}

	auto Found = 0;
	for(const auto &in : ind){
		if(in.dead == true && in.found == true){
			Found += 1;
		}
	}

    //L += Rec * log(1-param[p_rem].value-param[p_nd].value) + Found * log(param[p_rem].value) + (Rem-Found) * log(param[p_nd].value);
    L += Rec * log(1-param[p_rem].value) + Rem * log(param[p_rem].value) + (Rem-Found) * log(param[p_nd].value) + Found * log(1-param[p_nd].value) + histI * log(param[p_histimm].value) + (ind.size() - histI) * log(1-param[p_histimm].value);
	//L += Rec * log(1-param[p_rem].value) + (Rem-Found) * log(param[p_rem].value*param[p_nd].value) + Found * log(param[p_rem].value*(1-param[p_nd].value));

	if(index_max > 1){
		L -= choose(ind.size(),num_index);   // This implements the likelihood of distribution in the index cases at T0
	}

	time_likelihood += clock();
	return L;
}

/// Calculate the likelihood for the 
double Individual::likelihood_gamma()
{
	if(inf == false) return 0;

	return gamma_probability(latentP,param[mu_L].value,param[sh_L].value) +
       	 gamma_probability(infectiousP,param[mu_I].value,param[sh_I].value);
}

/// Calculates the prior for the number of index cases
double Colony::prior_index()
{
	auto num_index = 0u;                   
	for(const auto &in : ind){
		if(in.Exp_t == T0) num_index++;
	}
	return -alpha*num_index;
}

/// Sets transition parameters
TransParam set_trans_param()
{
	TransParam tp;
	
	tp.mu_L_value = param[mu_L].value;
	tp.sh_L_value = param[sh_L].value;
	tp.mu_I_value = param[mu_I].value;
	tp.sh_I_value = param[sh_I].value;
	
	tp.sh_L_sample = 1+phi_G*(tp.sh_L_value-1);
	tp.mu_L_sample = tp.mu_L_value*tp.sh_L_sample/(phi_G*tp.sh_L_value);
	tp.sh_I_sample = 1+phi_G*(tp.sh_I_value-1);
	tp.mu_I_sample = tp.mu_I_value*tp.sh_I_sample/(phi_G*tp.sh_I_value);

	return tp;
}



/// Based on a specified exposure time this samples the individual's time line
void Individual::sample(const double expt, const TransParam &tp)
{
	latentP = gamma_sample(tp.mu_L_sample, tp.sh_L_sample);
	infectiousP = gamma_sample(tp.mu_I_sample,tp.sh_I_sample);
	
	Exp_t = expt; 
	Inf_t = Exp_t + latentP;
	Rec_t = Inf_t + infectiousP;
}

/// Very slightly modified version of the above which allows me to pass different values
void Individual::sample_sim(const double expt, const TransParam &tp)
{
	latentP = gamma_sample(tp.mu_L_value, tp.sh_L_value);
	infectiousP = gamma_sample(tp.mu_I_value,tp.sh_I_value);
	
	Exp_t = expt; 
	Inf_t = Exp_t + latentP;
	Rec_t = Inf_t + infectiousP;
}

/// Gets the probability of sampling a certain set of latent and infectious periods
double Individual::sample_probability(const TransParam &tp)
{
	return gamma_probability(latentP,tp.mu_L_sample,tp.sh_L_sample) +
       	 gamma_probability(infectiousP,tp.mu_I_sample,tp.sh_I_sample);
}

/// Based on a recovery time the exposure time and infecttous times are sampled
void Individual::sample_backwards(const TransParam &tp)
{
	latentP = gamma_sample(tp.mu_L_sample, tp.sh_L_sample);
	infectiousP = gamma_sample(tp.mu_I_sample,tp.sh_I_sample);

	Inf_t = Rec_t - infectiousP;
	Exp_t = Inf_t - latentP;
}


/// Adds a parameter to the list of parameters
int add_param(string name, PriorType prior_type, double prior_val1, double prior_val2, double max_val)
{
	auto p = param.size();

	double mean = UNSET;
	switch(prior_type){
		case UNIFORM_PRIOR: mean = (prior_val1+prior_val2)/2; break;
		case BETA_PRIOR: mean = prior_val1/(prior_val1+prior_val2); break;
		case GAMMA_PRIOR: mean = prior_val1; break;
		case TRUNC_GAMMA_PRIOR: mean = prior_val1; break;
		case EXP_PRIOR: mean = 1/prior_val1; break;
	}
		
	Param par;
	par.name = name;
	par.value = UNSET;
	par.prior_type = prior_type; par.prior_val1 = prior_val1; par.prior_val2 = prior_val2; par.max_val = max_val;
	if(mean >= 0){
		par.jump = mean/10; 
	} else {
		par.jump = -mean/50;
	}
	par.npropose = 0; par.naccept = 0;
	par.jump_MBP = mean/50;
	par.npropose_MBP = 0; par.naccept_MBP = 0; par.nfail_MBP = 0;
	param.push_back(par);
	
	return p;
}

/// Initialises trace plots
void trace_init()
{
	trace << "State"; 
	for(auto &pa : param) trace << "\t" << pa.name; 
	for(auto &co : colony) trace << "\tNinf Colony" << co.index;
    for(auto &co : colony) trace << "\tHisotically Immune" << co.index;
	for(auto &co : colony) trace << "\tTinit Colony" << co.index;
	for(auto &co : colony) trace << "\tInit.infected" << co.index;
	for(auto &co : colony) trace << "\tExt.infected" << co.index;	
	for(auto &co : colony) trace << "\tLi Colony" << co.index;
	for(auto &co : colony) trace << "\tPr index Colony" << co.index;
	trace << "\tPrior";
	trace << endl;	
}

/// Initialises ppc plots
void ppc_init()
{
	ppc << "State";
	ppc << "\tColony";
	ppc << "\tExposure";
	ppc << "\tInfection";
	ppc << "\tRemoval";
	ppc << "\tDead";
	ppc << "\tFound";
    ppc << "\tHistorically Immune";
	ppc << endl;	
}


/// Plots useful quantities
void trace_plot()
{
	trace << s;
	for(auto &pa : param) trace << "\t" << pa.value;
	
	for(auto &co : colony){
		auto ninf = 0u; for(const auto &in: co.ind){ if(in.inf == true) ninf++;}
		trace << "\t" << ninf;
	}

    for(auto &co : colony){
        auto nhistimm = 0u; for(const auto &in: co.ind){ if(in.histimm == true) nhistimm++;}
        trace << "\t" << nhistimm;
    } 
	
	for(auto &co : colony) trace << "\t" << co.T0;
	
	for(auto &co : colony){
		auto num = 0u;  for(const auto &in: co.ind){ if(in.inf == true && in.Exp_t == co.T0) num++;}
		trace << "\t" << num;
	}
	
	for(auto &co : colony){
		auto num = 0u;
		auto min_inf = LARGE;
		for(const auto &in: co.ind){
			if(in.inf == true && in.Inf_t < min_inf) min_inf = in.Inf_t;
		}
		for(const auto &in: co.ind){
			if(in.inf == true && in.Exp_t < min_inf) num++;
		}
		trace << "\t" << num;
	}
	
	for(auto &co : colony) trace << "\t" << co.Li;
	
	for(auto &co : colony) trace << "\t" << co.Pr_index;
		
	trace << "\t" << Pr;
	
	trace << endl;
}

/// Plots useful quantities
void ppc_plot()
{
	for(auto &co : colony_ppc){
		for(const auto &in: co.ind){
			if(in.inf == true){
				ppc << s << "\t" << co.index << "\t" << in.Exp_t << "\t" << in.Inf_t << "\t" << in.Rec_t << "\t" << in.dead << "\t" << in.found << "\t" << in.histimm;
				ppc << endl;
			}
		}
	}
}

/// Plots useful quantities
void colony_write()
{
	for(auto &co : colony){
		for(const auto &in: co.ind){
			if(in.inf == true){
				colony_data << in.Rec_t << "\t" << in.dead << "\t" << in.found << "\t" << co.index;
				colony_data << endl;
			}
		}
	}
}

/// Initialises trace plots
void colony_init()
{
	colony_data << "Removal Time";
	colony_data << "Dead";
	colony_data << "Found";
	colony_data << "\tColony";
	colony_data << endl;	
}

/// Outputs information about how well the chain is mixing
void mcmc_diagnostics()
{
	for(auto &pa : param){
		cout << pa.name << " ";
		cout << int(100*pa.naccept/(pa.npropose+0.001)) << "% random walk acceptace probability   ";
		cout << "Jump size: " << pa.jump << "    "; 
		if(pa.npropose_MBP){
			cout << int(100*pa.naccept_MBP/(pa.npropose_MBP+0.001)) << "% MBP   ";
			cout << int(100*pa.nfail_MBP/(pa.npropose_MBP+0.001)) << "% fail MBP   ";
			cout << "MBP Jump size: " << pa.jump_MBP << "    "; 
		}
		
		cout << endl;
	}
	
	for(auto &co : colony){
		cout << "Colony " << co.index << ": ";
		cout << int(100*co.nmulti_accept/(co.nmulti_propose+0.001)) << "% multi AP       ";
		cout << int(100*co.nT0_accept/(co.nT0_propose+0.001)) << "% T0 AP        ";
		cout << int(100*co.nT0_joint_accept/(co.nT0_joint_propose+0.001)) << "% T0 joint AP   size:" << co.T0_joint_jump << "      ";
		cout << co.nprop << " Number of proposals          ";
		cout << endl;
	}
	
	cout << time_total/(60*CLOCKS_PER_SEC) << " Time in minutes" << endl;
	cout << int(100*time_likelihood/(time_total)) << " Percentage in likelihood" << endl;
}

/// Draws a sample from the uniform distribution
double uniform_sample(const double lwr, const double upr)
{
	uniform_real_distribution<double> distribution(lwr,upr);
	return distribution(generator);
}

/// Draws a sample from the gamma distribution with mean and shape
double gamma_sample(const double mu, const double sh)
{
	gamma_distribution<double> distribution(sh,mu/sh);
	return distribution(generator);
}

/// Draws a sample from the gamma distribution with mean and shape.  The name argument allows for truncation in different places dependent on the parameter
double trunc_gamma_sample(const double mu, const double sh, const string nm)
{
	gamma_distribution<double> distribution(sh,mu/sh);
	auto valid = false;
	auto value = 0.0;
	while(valid == false){
		value = distribution(generator);
		if(nm != "sh_L" && nm != "sh_I" && nm != "mu_L" && nm != "mu_I" && value < 20.0) valid = true;
		if((nm == "sh_L" || nm == "sh_I") && value > 5.0) valid = true;
		if(nm == "mu_L" && value < 10.0) valid = true;
		if(nm == "mu_I" && value < 15.0) valid = true;
	}
	return value;
}

/// Draws a sample from the beta distribution 
double beta_sample(const double a, const double b, const double c)
{
	auto proceed = true;
	auto X1 = 0.0;
	auto X2 = 0.0;
	auto sampled_value = 0.0;
	while(proceed == true){
		X1 = gamma_sample(a,a);
		X2 = gamma_sample(b,b);
		sampled_value = X1/(X1+X2);
		if(sampled_value < c) proceed = false;
	}
	return sampled_value;
}

/// The log of the probability from the gamma distribution
double gamma_probability(const double x, const double mu, const double sh)
{
	auto b = sh/mu;
	return (sh-1)*log(x) - b*x + sh*log(b) - lgamma(sh);
}

/// The log of the probability from the gamma distribution.  The nm parameter specifies how to deal with the truncation for each parameter
double trunc_gamma_probability(const double x, const double mu, const double sh, const string nm)
{
	auto b = sh/mu;
	auto prob = -LARGE;
	if(x > 5.0 && nm == "sh_I") prob = ((sh-1)*log(x) - b*x + sh*log(b) - lgamma(sh))/0.9940157;
	if(x > 5.0 && nm == "sh_L") prob = ((sh-1)*log(x) - b*x + sh*log(b) - lgamma(sh))/0.9896794;
	if(x < 15.0 && nm == "mu_L") prob = ((sh-1)*log(x) - b*x + sh*log(b) - lgamma(sh))/0.9995746;
	if(x < 15.0 && nm == "mu_I") prob = ((sh-1)*log(x) - b*x + sh*log(b) - lgamma(sh))/0.9653511;	
	if(x < 30.0 && nm != "sh_L" && nm != "sh_I" && nm != "mu_L" && nm != "mu_I" && nm != "T0") prob = ((sh-1)*log(x) - b*x + sh*log(b) - lgamma(sh))/0.9972577;
	if(x < 30.0 && nm == "T0") prob = ((sh-1)*log(x) - b*x + sh*log(b) - lgamma(sh))/0.9658502;
	
	return prob;
}

/// The log of the probability from the beta distribution
double beta_probability(const double x, const double a, const double b)
{
  return (a-1)*log(x) + (b-1)*log(1-x) - lgamma(a) - lgamma(b) + lgamma(a+b);
}

/// Samples from the exponential distribution using a rate
double exp_sample(const double rate)
{
	exponential_distribution<double> distribution(rate);
	return distribution(generator);
}

/// The log of the probability from the exponential distribution
double exp_probability(const double x, const double mu)
{
	return log(mu) - x*mu;
}

/// Draws a normally distributed number with mean mu and standard deviation sd
double normal_sample(const double mu, const double sd)
{
	normal_distribution<double> distribution(mu,sd);
	return distribution(generator);
}

/// The log of the probability from the normal distribution
double normal_probability(const double x, const double mean, const double var)
{
  if(var <= 0) emsg("Negative");
  return -0.5*log(2*M_PI*var) - (x-mean)*(x-mean)/(2*var);
}

vector <double> logsum;

/// Initialises calculation of choose function
void choose_init()                 
{
	auto Nmax = 0u;
	for(const auto &co : colony){
		if(co.ind.size() > Nmax) Nmax = co.ind.size();
	}
	Nmax++;
	logsum.resize(Nmax);
	
	auto sum = 0.0;
	for(auto n = 1u; n < Nmax; n++){
		sum += log(n);
		logsum[n] = sum;
	}
}

/// Calculates N choose n
double choose(int N, int n)
{
	if(N < 1 || n < 1 || n > N) emsg("choose problem"); 
	return logsum[N] - logsum[N-n] - logsum[n];
}

/// Calculates N choose n
double choose_mult(int N, int n1, int n2)
{
	if(N < 1 || (n1+n2) < 1 || n1 > N || n2 > N) emsg("choose problem"); 
	return logsum[N] - logsum[N-n1-n2] - logsum[n1] - logsum[n2];
}
	
/// Displays an error message
void emsg(const string& msg)
{
	cout << msg << endl;
	exit (EXIT_FAILURE);
}

/// Performs proposals which changes T0 and alters other event time 
void Colony::joint_T0_event_proposal()
{	
	auto T0_store = T0;
	
	auto d_T0 = normal_sample(0,T0_joint_jump);
	
	T0 += d_T0;
		
	auto T_end = T_end;                                  // This is the maximum time over which scaling is performed

	auto ind_store = ind;

	auto illegal = false;

	auto dLi_gamma = 0.0;

	auto fac = 0.0;  // This is used to account for changes in phase space under the proposal
	for(auto &in : ind){
		if(in.inf == true){
			auto t = in.Rec_t;
			if(t > T_end) t = T_end;
			auto T = t-T0_store;
			
			if(in.Exp_t == T0_store){
				in.Exp_t = T0;
				if(T0 > t){ illegal = true; break;}
			}
			else{
				auto frac = (t-in.Exp_t)/T;
				if(frac < 0) illegal = true;
				fac += log(1-(d_T0/T));
				in.Exp_t += d_T0*frac; 
				if(in.Exp_t > T_end){ illegal = true; break;}
			}
					
			auto frac = (t-in.Inf_t)/T;
			if(frac > 0){
				fac += log(1-(d_T0/T));
				in.Inf_t += d_T0*frac; 
			}
			
			in.latentP = in.Inf_t - in.Exp_t;
			in.infectiousP = in.Rec_t - in.Inf_t;
		
			auto Li_gamma_prop = in.likelihood_gamma();
			dLi_gamma += Li_gamma_prop - in.Li_gamma;
			in.Li_gamma = Li_gamma_prop;
			
			if(in.latentP <= 0 || in.infectiousP <= 0){ illegal = true; break;}
		}
	}
		
	double al, Li_prop;
	if(illegal == true) al = 0;
	else{	
		Li_prop = likelihood();   
		
		al = exp(phi_L*(Li_prop - Li) + phi_G*dLi_gamma + fac); 
	}

	nT0_joint_propose++;
	if(uniform_sample(0.0,1.0) < al){
		nT0_joint_accept++;
		Li = Li_prop;
		if(s < nburnin) T0_joint_jump *= 1.001;                // The dynamically adapts the jumping size
	}
	else{
		T0 = T0_store;
		ind = ind_store;
		if(s < nburnin) T0_joint_jump *= 0.9995;
	}
}


/// This checks that all quantities are correctly specified
void check(int num)
{
	double dd;
	
	for(auto &co : colony){
		for(const auto &in : co.ind){                                         // Checks timings are correctly specified
			if(in.latentP < 0){
				cout << "  Exposed: " << in.Exp_t << "      ";
				cout << "  Infected: " << in.Inf_t << "      ";
				cout << "  Removed: " << in.Rec_t << "      " << endl;
				emsg("Problem latentP");
			}
			if(in.infectiousP < 0) emsg("Problem infectiousP");
			
			dd = in.latentP -	(in.Inf_t - in.Exp_t);
			if(dd*dd > TINY) emsg("Problem latentP2");
		
			dd = in.infectiousP -	(in.Rec_t - in.Inf_t);
			if(dd*dd > TINY)	emsg("Problem infectiousP2");
		}
		
		auto T = LARGE;                                                   // CHecks T0 correctly specified
		for(const auto &in : co.ind){
			if(in.inf == true){
				if(in.Exp_t < T) T = in.Exp_t;
				
				if(in.Exp_t > co.T_end) emsg("Infection must be before T_end");
				
				if(in.Rec_t > co.T_end) emsg("Must recover before T_end");
			}
		}
		if(T != co.T0) emsg ("T0 wrong");
		
		dd = co.likelihood() - co.Li;	
		if(dd*dd > TINY){
			emsg("Problem likelihood"+to_string(num));
		}
		
		for(auto &in : co.ind){     
			dd = in.likelihood_gamma() - in.Li_gamma;
			if(dd*dd > TINY){
				emsg("Problem likelihood_gamma");
			}
		}	
		
		dd = co.prior_index() - co.Pr_index;	
		if(dd*dd > TINY) emsg("Problem Pr_index");
	}
	
	dd = prior() - Pr;	
	if(dd*dd > TINY) emsg("Problem Pr");
}

