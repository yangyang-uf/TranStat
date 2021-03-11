/*datastruct.h*/

// CONTACT is a structure to record information for a single contact.
// 'contact_mode' is the type of the contact, 'infective_prob' is the probability of
// the case being infectious. 'covariates' is an array recording the covariate values
// associate with this contact being infectious. 
// 'size' is the frequency of contacts with the same values for contact_mode, infective prob
// and covariates on that day. For analytic purpose, we do not bother to record each single contact, 
// but only record the frequency of contacts of the same type.
// 'ignore' is designed to allow more flexibility in analysis. contacts with 'ignore'=1 will
// be excluded from analysis. 'next' is a pointer, we use CONTACT for chain operations, because
// we do not know the number of meaningful contacts for each individual each day, and the number
// may be small or may be very large in some circumstances (such as school).

typedef struct contact_struct
{
   int contact_mode;
   int contact_id;
   double offset;
   struct contact_struct *next;
   struct contact_struct *pair;
} CONTACT;

//CONTACT_HISTORY is a structure to record all daily contacts.
//an array of this structure is needed for each individual with the length
//equal to the epidemic duration, each array element recording contacts of that day.
typedef struct
{
   CONTACT *c2p_contact; //pointer to head node of c2p contacts
   CONTACT *c2p_contact_rear; //pointer to rear node of c2p contacts
   CONTACT *p2p_contact; //pointer to head node of p2p contacts
   CONTACT *p2p_contact_rear; //pointer to rear node of p2p contacts
} CONTACT_HISTORY;

// RISK is a structure similar to CONTACT. It records the information
// about the risk of transmission (or infection), or equivalently,
// the exposure information. This structure is useful in the following sense.
// Instead of infering exposure by going through all contacts one by one,
// for estimation purpose we only care how many (captured by "size") infective sources a susceptible
// person is exposed to on a given day for a given type of source. The type of source
// is defined by the covariates of the source, his/her infectivity level (captured
// by "infective_prob", the symptom status (captured by "symptom"), and the contact mode.
// "offset" is an auxilliary variable, e.g., for introducing time-varying c2p exposure.
typedef struct risk_struct
{
   int contact_mode;
   int symptom;
   double infective_prob;
   double *covariate;
   double offset;
   int size;
   struct risk_struct *next;
} RISK;

typedef struct individual_risk_struct
{
   int *infective_subject;
   int *contact_mode;
   double *infective_prob;
   int size;
} RISK_INDIVIDUAL;

//RISK_HISTORY is a structure to record all daily risks of infection.
//An array of this structure is needed for each individual with the length
//equal to the epidemic duration, each array element recording risks of that day.
typedef struct
{
   RISK *c2p_risk; //pointer to head node of c2p risk profile
   RISK *c2p_risk_rear; //pointer to rear node of c2p risk profile
   RISK *p2p_risk; //pointer to head node of p2p risk profile
   RISK *p2p_risk_rear; //pointer to rear node of p2p risk profile
} RISK_HISTORY;

//See below in the definition of CFG_PARS structure for explanation of the PAR_EQUICLASS structure.
typedef struct
{
   int size;
   int *member;
} PAR_EQUICLASS;


typedef struct int_chain
{
   int id;
   struct int_chain *next;
} INTEGER_CHAIN;

typedef struct risk_group
{
   INTEGER_CHAIN *member;
   INTEGER_CHAIN *member_rear;
   int size;
   struct risk_group *next;
   RISK_HISTORY *risk_history;
} RISK_CLASS;

// STATUS means the final outcome of a subject. We consider four types of outcome:
// 1. symptomatic infection; 2. asymptomatic infection; 3. escape any infection; 4. pre-season immune.
// The first 3 types assume the subject is susceptible at the beginning of the epidemic.
// An outcome is defined by two variables: infection or not and the infectiousness onset day.
// For types 1, 3 and 4, the two variables are both known. Type 1 means infection=yes and 
// infectiousness onset day is the observed symptom onset day (by our assumption that incubation period
// and latent period overlap). Types 3 and 4 both have infection=no, and infectiousness onset day is defined
// as INFINITY for the former and -INFINITY for the latter. For type 2, infection=yes, but 
// the infectiousness onset day is not observable, and therefore type 2 corresponds to many possible
// outcomes, each defined by a possible infectiousness onset day. 
// STATUS is a structure built for a community. The array "status" stores a possible combination of status
// of people with uncertain status in the community. In each of these people, there is an 2 by M array of possible
// outcomes, M being the number of possible outcomes. Each element of the array "status" is an integer
// refering to the location of an assigned outcome for the corresponding person.
// "log_L_ini" and "log_L_cur" are log-likelihood of the community evaluated at the initial and current
// estimates of parameters. NOTE, initial estimates here are not the same as the initial parameter values
// input by the user in the configuration file. These initial estimates are obtained using MCMC sampling
// and are used for importance sampling.  "weight" is the weight of importance sampling,
// used to perform the E-step in the MCEM algorithm (similar to numerical integration).

typedef struct my_state
{
   double log_L_ini;
   double log_L_cur;
   double log_L_pre;
   double ratio;
   int *states;
   struct my_state *next;
} STATE;

/*struture of individual subject*/
typedef struct {
        int id; 
        int idx; /* 1: index case, 0: contact*/
        int community; 
        int pre_immune;
        int infection;
        int symptom;
        int day_ill; // Infectiousness onset day. For symptomatic infeciton, this is also illness onset day
        int any_exposure; //an indicator used for checking model admissibility
        int day_infection;
        int day_infection_lower; //lower bound of possible infection days
        int day_infection_upper; //upper bound of possible infection days
        int day_infective_lower; //lower bound of days being infectious
        int day_infective_upper; //upper bound of days being infectious
        int final_risk_day;
        int exit;
        int day_exit;
        int u_mode; //the pathogenicity group
        int q_mode; //the preseason immunity group
        int ignore; //people with ignore=1 will be excluded form analysis
        double *time_ind_covariate; //time-independent covariates
        double **time_dep_covariate; //time-dependent covariates
        CONTACT_HISTORY *contact_history; //contact history over the epidemic duration
        RISK_HISTORY *risk_history; //contact history over the epidemic duration
        RISK_CLASS  *risk_class;
        double *imm_covariate; //covariates affecting preseason immunity
        double **pat_covariate; //covariates affecting pathogenicity
        int **possible_states;
        int size_possible_states;
        int current_state;
        double weight;
}PEOPLE;

//struture of community. A community defined in TranStat is a closed setting
//where people make contacts with each other and henceforth transmisison is possible.
//In other words, people in the same community are connected, but people in different communities
//are disconnected. A community does not have to be a physical community in common sense.
//For example, in a household-based study, if all households are considered independent, 
//Each household is considered as a community in TranStat.
typedef struct {
        int id; 
        int day_epi_start; //starting day of the epidemic to be analyzed
        int day_epi_stop; //stopping day of the epidemic to be analyzed
        int epi_duration; //duration of the epidemic, must be day_epi_stop - day_epi_start + 1
        int index_case; //the index case. This field is only useful for simulation.
        int idx_day_ill; //ILI onset day of index case. This field is only useful for simulation.
        int *member; //array storing all members of the community
        int size; //number of individuals in the community
        int *idx; //array storing all members of the community
        int size_idx; //number of individuals in the community
        int earliest_idx_day_ill; //ILI onset day of the earliest index case
        int latest_idx_day_ill; //ILI onset day of the latest index case
        int *member_impute; //array storing all members with uncertain status
        int size_impute; //number of individuals with uncertain status
        int size_possible_states;
        int counter; //An auxilliray variable for constructing the community
        int ignore; // for case ascertained design, families without index case will be ignored
        CONTACT_HISTORY *contact_history; //contact history over the epidemic duration, used if all people share the same contact history

        RISK_CLASS *risk_class;
        RISK_CLASS *risk_class_rear;

        STATE *sample_states; //array of potential status for members with uncertain status, 
                               //used for importance sampling in the Monte Carlo EM (MCEM) algorithm 
        STATE *sample_states_rear; //pointer to the last node in the array of sample status. 
                                    //It is needed for augmentation of the array in the MCEM algorithm.
        STATE *list_states; //array of potential status for members with uncertain status, 
                               //used for importance sampling in the Monte Carlo EM (MCEM) algorithm 
        STATE *list_states_rear; //pointer to the last node in the array of sample status. 
                                    //It is needed for augmentation of the array in the MCEM algorithm.
}COMMUNITY;

/* configuration parameters for the disease and the trial.
   some parameters are used only for estimation, 
   some are used only for simulations, 
   while others are used for both*/
typedef struct {
    //=============================================================
    // configuration parameters shared by simulation and estimation
    //=============================================================
    char path_in[100];
    char path_out[100];
    
    int n_b_mode; // number of c2p contact modes
    int n_p_mode; // number of p2p contact modes
    int n_u_mode; // number of pathogenicity groups
    int n_q_mode; // number of pre-immunity groups
    int n_c2p_covariate; //number of covariates affecting c2p transmission
    int n_sus_p2p_covariate; //number of covariates affecting susceptibility for p2p transmission
    int n_inf_p2p_covariate; //number of covariates affecting infectiousness for p2p transmission
    int n_int_p2p_covariate; //number of interactions for p2p transmission
    int n_p2p_covariate; //number of covariates for p2p transmission, must be
                         //n_sus_p2p_covariate +  n_inf_p2p_covariate + n_int_p2p_covariate
    int n_imm_covariate; //number of covariates for pre-immunity, which should be <= n_time_ind_covariate
    int n_pat_covariate; //number of covariates for pre-pathogenicity
    int n_par; //total number of parameters to be estimated.
               // must be n_b_mode + n_p_mode + n_c2p_covariate + n_p2p_covariate
    int n_time_ind_covariate; //number of time-independent covariates
    int n_time_dep_covariate; //number of time-dependent covariates
    int n_covariate; //total number of covariates, must be n_time_ind_covariate + n_time_dep_covariate
    
    int illness_as_covariate; // whether the illness period is used as a time dependent covariate.
    int illness_covariate_id; // This is useful if you want to compare infectivity level during incubation period vs during illness period.
    
    int min_incubation;  //the minimum duration of incubation period.
    int max_incubation;  //the maximum duration of incubation period.
    double * prob_incubation;  //empirical probability for incubation period.
    
    int lower_infectious; //the lower and upper bounds of days that a case is infective, given the ILI sonset day.
    int upper_infectious; // a typical value is 0 (ILI onset day) or 1 (one day after ILI onset).
    double *prob_infectious;  // empirical probability of being infectious for a given day in the infectious period

    int n_inc;
    int *min_inc;  //the minimum duration of incubation period.
    int *max_inc;  //the maximum duration of incubation period.
    double **prob_inc;  //empirical probability for incubation period.
    
    int n_inf;
    int *lower_inf; //the lower and upper bounds of days that a case is infective, given the ILI sonset day.
    int *upper_inf; // a typical value is 0 (ILI onset day) or 1 (one day after ILI onset).
    double **prob_inf;  // empirical probability of being infectious for a given day in the infectious period
    
    int *c2p_covariate; //indicating which covariates modify c2p susceptibility
    int *sus_p2p_covariate;//indicating which covariates modify p2p susceptibility
    int *inf_p2p_covariate;//indicating which covariates modify p2p infectiousness
    int **interaction; //indicating which covariates in the list sus_p2p_covariate and the list
                       //inf_p2p_covariate form interaction that affect p2p transmission.
    int *imm_covariate; //indicating which covariates modify pre-immunity
    int *pat_covariate; //indicating which covariates modify pathogenicity

    //User may specify which parameters share the same value, i.e., force them to be equal.
    //We group parameters into so-called equivalence classes, all parameters in the same class being
    //equal. For the convenience of programing, we require user to specify equivalence classes even
    //if no parameters are considered equal; that is, each parameter alone is considered an equivalence class.
    //Note: parameters to be fixed as constants (see below) should not be listed in any equivalence class;
    //In other words, only parameters to be estimated are specified in equivalence classes.
    int n_par_equiclass;
    PAR_EQUICLASS *par_equiclass;
    
    //User may specify which parameters should be hold fixed instead of to be estimated.
    //This is very useful for sensitivity analysis or profile likelihood.
    //Fixed parameters should not be listed in any equivalence class defined above.
    int n_par_fixed;
    int *par_fixed_id;
    double *par_fixed_value;

    //============================================
    // configuration parameters used by simulation
    //============================================
    // Simulation is a useful tool for learning or validating TranStat and for investigating
    // intervention strategies.
    int simulation; // Do you want to simulate epidemics?
    int n_simulation; //Number of simulations to be performed
    //The following two arrays specify parameter values used for simulation.
    // "sim_par_effective" refer to parameters you want to estimate for the simulated epidemic,
    // and "sim_par_fixed" refer to parameters you want to fix at given values 
    // during estimation for the simulated epidemic,
    double *sim_par_effective; //parameter values used for simulation. If you have parameters
                               // that need be fixed during estimation, use above "n_par_fixed",
                               // "par_fixed_id" and "par_fixed_value" to input the information.
    double prop_mix_imm_esc; //proportion of subjects with ambiguity about pre-immunity and escape.
                             // People without symptom onset and not confirmed as cases are either preimmune
                             // or escapes of infection. We may not the exact status for some of them but not for all.

    double asym_effect_sim; //relative infectiouness of asymptomatic infection compared to symptomatic infection, used for simulation.

    //============================================
    // configuration parameters used by estimation
    //============================================
    //If people in the same community contact with each other every day 
    //and there is only one c2p contact mode and one p2p contact mode, user can
    //ask TranStat to generate c2p contact and p2p contact files
    int generate_c2p_contact;
    int generate_p2p_contact;
    int c2p_offset;
    int p2p_offset;
    // Do all people of all communities share the same c2p contact history?
    // This flag is to simplify the c2p cotact file.
    //int common_c2p_contact_history_across_community;
    // Do people in the same community share the same contact history?
    int common_contact_history_within_community;
    
    int optimization_choice; //0: Newton-Raphson, 1:Newton-Raphson + Nelder-Mead,  2:Nelder-Mead
    int search_bound_provided; //Are search bounds of the Nelder-Mead approach provided? 1:yes, 0:no
    double *lower_search_bound;
    double *upper_search_bound;
    int ini_par_provided; // Are initial parameter values provided for optimization?
    int n_ini; //number of sets of initial estimates
    double **ini_par_effective; //Values of initial estimates
    int converge_criteria_provided; //Is convergence criteria provided?
    double *converge_criteria;//threshholds for judging convergence in parameter estimates

    int EM; //Will EM algorithm be performed? 1:yes, 0:no
    int min_size_MCEM; //minimum number of possible infection status combinations in a household
                       //to consider MCEM. For lower numbers, we simply use the EM, 
                       //averaging over all possible status combinations
    int community_specific_weighting; // in traditional MCEM algorithms, imputation (and sampling) is done at the population level,
                       // so the importance weights are also calculated at the population level, i.e., all communities share the same
                       // numbers and values of weights. In our setting, communities are independent of each other, and sampling is
                       // done independently for each community, so we can use community-level likelihoods to calculate community-specific
                       // weights, and the number of importance samples can also vary from community to community.
                       // 1: turn on community-level weighting, 0: keep population-level weighting.
    int n_base_sampling; //number of importance samples to start the EM algorithm.
                         //This base number will be augmented when necessary during the algorithm.
    int n_burnin_sampling; //number of burn-in samples in MCMC sampling to generate importance samples.
    int n_burnin_iter; //Before the official runs of the MCEM algorithm, 
                       //we need a number of burn-in iterations of MCMC sampling to get sensible initial estimates
                       // based on which importance samples are generated.
    double asym_effect_est; //relative infectiouness of asymptomatic infection compared to symptomatic infection, used for estimating other parameters.
                            //This quantity is not to be estimated, as the information about it is generally scarce.
                            //We suggest user perform sensitivity analysis by adjusting its value.

    // MCE stands for Monte Carlo error, the extra variation due to sampling in the MCEM algorithm.
    // n_sampling_for_mce gives the number of samplings to adjust for MCE. If it's 0, then no adjustment is done.
    // use_bootstrap_for_mce indicates whether bootstrap from old importance samples or draw new samples via MCMC.
    // It's much faster to bootstrap, but bootstrap likely requires relatively large number of old importance samples.    
    int n_sampling_for_mce;
    int use_bootstrap_for_mce;
    int skip_Evar_for_mce;
    int preset_index; //if preset_index=0, then index cases may change when onset times of some infections are resampled; otherwise, index cases are fixed.
    int simplify_output; //if simulation or multiple imputation is performed, it may be desired to output results without any comments.
    int simplify_output_SAR;
    int simplify_output_R0;
    int check_missingness;
    int check_runtime; // whether output the running times for the importance sampling step and the EM iteration steps or not.
    int check_mixing; // whether output importance samples, basically sampled infectiouseness onset times for asymptomatic infections
                      // to check how well they mix.

    int goodness_of_fit; //Should files for assessing goodness-of-fit be generated?
                         //It
    int adjust_for_left_truncation;//indicate whether to perform adjustment for left truncation.
                                   // A more appropriate term than "left truncation" is "selection bias".
                                   // It is necessary to adjust for "selection bias" when all communities
                                   // are enrolled upon identification of index cases, which is called
                                   // the case-ascertained design. See Yang et al. (2006, appl stat) for details.
    int adjust_for_right_censoring;//Indicate whether to perform adjustment for right censoring of ILI.
                                   //People who have not shown symptoms by the end of an epidemic may have
                                   // already been infected with a incubation period so long that the symptom onset
                                   // is not observed by the last day of follow-up. If the follow-up period
                                   // is long enough, the impact of right-censoring is very small.
    int CPI_duration; //duration of the epidemic for calculating CPI
    // the following are covariate values used for calculating SAR
    int SAR_covariate_provided;
    double *SAR_sus_time_ind_covariate;
    double *SAR_inf_time_ind_covariate;
    int SAR_time_dep_lower;
    int SAR_time_dep_upper;
    MATRIX SAR_sus_time_dep_covariate;
    MATRIX SAR_inf_time_dep_covariate;
    //The following two quantities are used for calculating R0 
    int n_R0_multiplier; // number of sets of multipliers
    double *R0_multiplier; // R0_multiplier * SAR is the estimated R0 accounting for transmission in and out of household. 
                          // One possible such multiplier is given in Yang et al (2009)
    double *R0_multiplier_var; // variance of R0_multiplier, used to calculate variance of R0
    
    // R0_by_time indicates whether divide the epidemic period into two sub-periods by a series of cut-points so that one can plot how R0 estimates varies over time.
    // For each cut point, transmission probability p is assumed different before and after the cut point and remain constant within each of the two sub-period.
    // The cut points are set by R0_divide_start and R0_divide_stop. The cut points should be between day_epi_start and day_epi_stop of all communities (not even equal).
    int R0_divide_by_time; 
    int R0_divide_start;
    int R0_divide_stop;
    int R0_window_size;
    int output_R0_only; // whether only output R0 estimates and ignore other parameter estimates.
    
    //The following two are effective lower  and upper bound of infectious period for calculating SAR and R.
    //In some settings such as school, a case may withdraw to home soon after symptom onset,
    //so it is unrealistic to use the whole biological infectious period to calculate school-based SAR.
    int effective_bounds_provided;
    int *effective_lower_infectious; 
    int *effective_upper_infectious; 
                                     
    int stat_test; //Should statistical test for the existence of human-to-human transmission
                   // be performed? This function is yet to be built in.

    int skip_variance;
    int skip_output;
    int silent_run; //When TranStat is running, some  messages are shown. These messages are likely not wanted
                    // in simulations. User can surpress the messages using this switch.
    int write_error_log;                
}CFG_PARS;

//We define global variables here. These variables are used throughout the running of TranStat.
  double *b, *p, *u, *q, *lb, *lp, *lu, *lq, *CPI, *SAR, *R0;
  double *coeff_c2p, *coeff_p2p, *coeff_pat, *coeff_imm, *OR_c2p, *OR_p2p, *OR_pat, *OR_imm;
  double *se_b, *se_p, *se_u, *se_q, *se_coeff_c2p, *se_coeff_p2p, *se_coeff_pat, *se_coeff_imm;
  double *se_lb, *se_lp, *se_lu, *se_lq, *se_OR_c2p, *se_OR_p2p, *se_OR_pat, *se_OR_imm;
  double *se_CPI, *se_SAR, *se_R0;
  double *lower_b, *lower_p, *lower_u, *lower_q;
  double *lower_OR_c2p, *lower_OR_p2p, *lower_OR_pat, *lower_OR_imm;
  double *lower_CPI, *lower_SAR, *lower_R0;
  double *upper_b, *upper_p, *upper_u, *upper_q;
  double *upper_OR_c2p, *upper_OR_p2p, *upper_OR_pat, *upper_OR_imm;
  double *upper_CPI, *upper_SAR, *upper_R0;
  
  double *ee, *cum_ee;
  double *log_f_lb, *log_f_c2p, *log_f_lp, *log_f_p2p;
  double *log_f_lb_lb, *log_f_lp_lp;
  double **log_f_lb_c2p, **log_f_lp_p2p, **log_f_c2p_c2p, **log_f_p2p_p2p;

  double *log_e_lb, *log_e_c2p, *log_e_lp, *log_e_p2p;
  double *log_e_lb_lb, *log_e_lp_lp;
  double **log_e_lb_c2p, **log_e_lp_p2p, **log_e_c2p_c2p, **log_e_p2p_p2p;

  double **log_ee_lb, **log_ee_c2p, **log_ee_lp, **log_ee_p2p;
  double **log_ee_lb_lb, **log_ee_lp_lp;
  double ***log_ee_lb_c2p, ***log_ee_c2p_c2p, ***log_ee_lp_p2p, ***log_ee_p2p_p2p;

  double *cum_log_ee, **cum_log_ee_lb, **cum_log_ee_c2p, **cum_log_ee_lp, **cum_log_ee_p2p;
  double **cum_log_ee_lb_lb, **cum_log_ee_lp_lp;
  double ***cum_log_ee_lb_c2p, ***cum_log_ee_c2p_c2p, ***cum_log_ee_lp_p2p, ***cum_log_ee_p2p_p2p;

  double *cum_log_e_lb, *cum_log_e_c2p, *cum_log_e_lp, *cum_log_e_p2p;
  double *cum_log_e_lb_lb, *cum_log_e_lp_lp;
  double **cum_log_e_lb_c2p, **cum_log_e_c2p_c2p, **cum_log_e_lp_p2p, **cum_log_e_p2p_p2p;

  double *temp_lb, *temp_c2p, *temp_lp, *temp_p2p;
  double *temp_lb_lb, *temp_lp_lp;
  double **temp_lb_c2p, **temp_c2p_c2p, **temp_lp_p2p, **temp_p2p_p2p;

  double *log_day_L_lb, *log_day_L_c2p, *log_day_L_lp, *log_day_L_p2p;
  double **log_day_L_lb_lb, **log_day_L_lb_c2p, **log_day_L_lb_lp, **log_day_L_lb_p2p;
  double **log_day_L_c2p_c2p, **log_day_L_c2p_p2p, **log_day_L_lp_lp, **log_day_L_lp_p2p, **log_day_L_lp_c2p, **log_day_L_p2p_p2p;

  double *L_lb,  *L_lp;
  double *L_c2p, *L_p2p;
  double **L_lb_lb, **L_lb_lp, **L_lb_c2p, **L_lb_p2p;
  double **L_lp_lp, **L_lp_c2p, **L_lp_p2p;
  double **L_c2p_c2p, **L_c2p_p2p, **L_p2p_p2p;
  double L_lu, L_lq, L_lu_lu, L_lu_lq, L_lq_lq;
  double *L_lb_lu, *L_lp_lu, *L_lu_c2p, *L_lu_p2p, *L_lu_imm, *L_lu_pat;
  double *L_lb_lq, *L_lp_lq, *L_lq_c2p, *L_lq_p2p, *L_lq_imm, *L_lq_pat;
  double *L_pat, *L_imm;
  double **L_lb_pat, **L_lp_pat, **L_c2p_pat, **L_p2p_pat;
  double **L_lb_imm, **L_lp_imm, **L_c2p_imm, **L_p2p_imm;
  double **L_pat_pat, **L_pat_imm, **L_imm_imm;

  double *LQ_lb,  *LQ_lp;
  double *LQ_c2p, *LQ_p2p;
  double **LQ_lb_lb, **LQ_lb_lp, **LQ_lb_c2p, **LQ_lb_p2p;
  double **LQ_lp_lp, **LQ_lp_c2p, **LQ_lp_p2p;
  double **LQ_c2p_c2p, **LQ_c2p_p2p, **LQ_p2p_p2p;
  double LQ_lu, LQ_lq, LQ_lu_lu, LQ_lu_lq, LQ_lq_lq;
  double *LQ_lb_lu, *LQ_lp_lu, *LQ_lu_c2p, *LQ_lu_p2p, *LQ_lu_imm, *LQ_lu_pat;
  double *LQ_lb_lq, *LQ_lp_lq, *LQ_lq_c2p, *LQ_lq_p2p, *LQ_lq_imm, *LQ_lq_pat;
  double *LQ_pat, *LQ_imm;
  double **LQ_lb_pat, **LQ_lp_pat, **LQ_c2p_pat, **LQ_p2p_pat;
  double **LQ_lb_imm, **LQ_lp_imm, **LQ_c2p_imm, **LQ_p2p_imm;
  double **LQ_pat_pat, **LQ_pat_imm, **LQ_imm_imm;

  double *log_L_lb,  *log_L_lp;
  double *log_L_c2p, *log_L_p2p;
  double **log_L_lb_lb, **log_L_lb_lp, **log_L_lb_c2p, **log_L_lb_p2p;
  double **log_L_lp_lp, **log_L_lp_c2p, **log_L_lp_p2p;
  double **log_L_c2p_c2p, **log_L_c2p_p2p, **log_L_p2p_p2p;
  double log_L_lu, log_L_lq, log_L_lu_lu, log_L_lu_lq, log_L_lq_lq;
  double *log_L_lb_lu, *log_L_lp_lu, *log_L_lu_c2p, *log_L_lu_p2p, *log_L_lu_imm, *log_L_lu_pat;
  double *log_L_lb_lq, *log_L_lp_lq, *log_L_lq_c2p, *log_L_lq_p2p, *log_L_lq_imm, *log_L_lq_pat;
  double *log_L_pat, *log_L_imm;
  double **log_L_lb_pat, **log_L_lp_pat, **log_L_c2p_pat, **log_L_p2p_pat;
  double **log_L_lb_imm, **log_L_lp_imm, **log_L_c2p_imm, **log_L_p2p_imm;
  double **log_L_pat_pat, **log_L_pat_imm, **log_L_imm_imm;

  double *my_log_L_lb,  *my_log_L_lp;
  double *my_log_L_c2p, *my_log_L_p2p;
  double **my_log_L_lb_lb, **my_log_L_lb_lp, **my_log_L_lb_c2p, **my_log_L_lb_p2p;
  double **my_log_L_lp_lp, **my_log_L_lp_c2p, **my_log_L_lp_p2p;
  double **my_log_L_c2p_c2p, **my_log_L_c2p_p2p, **my_log_L_p2p_p2p;
  double my_log_L_lu, my_log_L_lq, my_log_L_lu_lu, my_log_L_lu_lq, my_log_L_lq_lq;
  double *my_log_L_lb_lu, *my_log_L_lp_lu, *my_log_L_lu_c2p, *my_log_L_lu_p2p, *my_log_L_lu_imm, *my_log_L_lu_pat;
  double *my_log_L_lb_lq, *my_log_L_lp_lq, *my_log_L_lq_c2p, *my_log_L_lq_p2p, *my_log_L_lq_imm, *my_log_L_lq_pat;
  double *my_log_L_pat, *my_log_L_imm;
  double **my_log_L_lb_pat, **my_log_L_lp_pat, **my_log_L_c2p_pat, **my_log_L_p2p_pat;
  double **my_log_L_lb_imm, **my_log_L_lp_imm, **my_log_L_c2p_imm, **my_log_L_p2p_imm;
  double **my_log_L_pat_pat, **my_log_L_pat_imm, **my_log_L_imm_imm;

  double Q_lq, *Q_imm, Q_lq_lq, *Q_lq_imm, **Q_imm_imm;
  double U_lu, *U_pat, U_lu_lu, *U_lu_pat, **U_pat_pat;
  double pr_lu, *pr_pat, pr_lu_lu, *pr_lu_pat, **pr_pat_pat;

  double *score_lb, *score_lp, *score_lu, *score_lq, *score_c2p, *score_p2p, *score_pat, *score_imm;
  double **info_lb_lb, **info_lb_lp, **info_lb_c2p, **info_lb_p2p;
  double **info_lp_lp, **info_lp_c2p, **info_lp_p2p;
  double **info_c2p_c2p, **info_c2p_p2p, **info_p2p_p2p;
  double **info_lu_lu, **info_lq_lq, **info_lu_lq;
  double **info_lb_lu, **info_lp_lu, **info_lu_c2p, **info_lu_p2p, **info_lu_pat, **info_lu_imm;
  double **info_lb_lq, **info_lp_lq, **info_lq_c2p, **info_lq_p2p, **info_lq_pat, **info_lq_imm;
  double **info_lb_pat, **info_lp_pat, **info_c2p_pat, **info_p2p_pat;
  double **info_lb_imm, **info_lp_imm, **info_c2p_imm, **info_p2p_imm;
  double **info_pat_pat, **info_pat_imm, **info_imm_imm;

       
