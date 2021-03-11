/*config.h*/
/*this function reads in configuration parameters about trial, disease and intervention*/

void RemoveSpaces(char* source)
{
  char* i = source;
  char* j = source;
  while(*j != 0)
  {
    *i = *j++;
    if(*i != ' ')
      i++;
  }
  *i = 0;
}

#include <string.h>
CFG_PARS get_cfg_pars(FILE *file)
{
  char c, *cc;
  int i, j, k, l, m, n, t, len, n_covariate;
  int start_day, stop_day;
  double *value;
  char st[100];
  char string[100];
  CFG_PARS cfg_pars;

  cfg_pars.path_in[0] = 0;
  cfg_pars.path_out[0] = 0;
  cfg_pars.min_incubation = cfg_pars.max_incubation = cfg_pars.lower_infectious = cfg_pars.upper_infectious =
  cfg_pars.n_b_mode = cfg_pars.n_p_mode = cfg_pars.n_u_mode = cfg_pars.n_q_mode =
  cfg_pars.n_covariate = cfg_pars.n_time_ind_covariate = cfg_pars.n_time_dep_covariate =
  cfg_pars.n_sus_p2p_covariate = cfg_pars.n_inf_p2p_covariate = cfg_pars.n_int_p2p_covariate =
  cfg_pars.n_par = cfg_pars.n_c2p_covariate = cfg_pars.n_p2p_covariate = 
  cfg_pars.n_pat_covariate = cfg_pars.n_imm_covariate =
  cfg_pars.n_par_equiclass = cfg_pars.n_par_fixed = 0;
  cfg_pars.generate_c2p_contact = cfg_pars.generate_p2p_contact = 0;
  cfg_pars.c2p_offset = cfg_pars.p2p_offset = 1;
  cfg_pars.check_missingness = cfg_pars.check_runtime = cfg_pars.check_mixing = 0;
  cfg_pars.simulation = cfg_pars.n_simulation = cfg_pars.EM = 0;
  cfg_pars.SAR_time_dep_lower = cfg_pars.SAR_time_dep_upper = 0;
  cfg_pars.n_R0_multiplier = cfg_pars.R0_divide_by_time = cfg_pars.R0_divide_start = cfg_pars.R0_divide_stop = cfg_pars.R0_window_size = 0;
  
  cfg_pars.effective_lower_infectious = cfg_pars.effective_upper_infectious = NULL;
  cfg_pars.prob_incubation = cfg_pars.prob_infectious = cfg_pars.converge_criteria = NULL; 
  cfg_pars.c2p_covariate = cfg_pars.pat_covariate = cfg_pars.imm_covariate = 
  cfg_pars.sus_p2p_covariate =  cfg_pars.inf_p2p_covariate = NULL; 
  cfg_pars.interaction = NULL; 
  cfg_pars.lower_search_bound = cfg_pars.upper_search_bound = NULL;
  cfg_pars.ini_par_effective = NULL;
  cfg_pars.sim_par_effective = NULL;
  cfg_pars.converge_criteria = NULL;
  cfg_pars.SAR_sus_time_ind_covariate = cfg_pars.SAR_inf_time_ind_covariate = NULL;
  initialize_matrix(&cfg_pars.SAR_sus_time_dep_covariate);
  initialize_matrix(&cfg_pars.SAR_inf_time_dep_covariate);
  cfg_pars.R0_multiplier = cfg_pars.R0_multiplier_var = NULL;
  cfg_pars.par_equiclass = NULL;
  cfg_pars.par_fixed_id = NULL;
  cfg_pars.par_fixed_value = NULL;
  
  cfg_pars.preset_index = 1; //For ebola analysis.
  cfg_pars.simplify_output = cfg_pars.simplify_output_SAR = cfg_pars.simplify_output_R0 = 0;
  /* The following code gets the input probabilites from config.file */

  rewind(file);

  
  while( ! feof(file))
  {
     c=fgetc(file);
     if(c == '#')
     {
        fgets(st, 100, file);
        cc = strchr(st, '\n');
        if(cc!=NULL) (*cc) = 0;
        cc = strchr(st, '\r');
        if(cc!=NULL) (*cc) = 0;
        RemoveSpaces(st);
        //==================================================================
        // configuration parameters shared by both simulation and estimation
        //==================================================================
        //printf("%s\n", st);
        if(strcmp(st, "input-path") == 0)  
        {
           fgets(cfg_pars.path_in, 100, file);
           cc = strchr(cfg_pars.path_in, '\n');
           if(cc!=NULL) (*cc) = 0;
           cc = strchr(cfg_pars.path_in, '\r');
           if(cc!=NULL) (*cc) = 0;
        }
        if(strcmp(st, "output-path") == 0)  
        {
           fgets(cfg_pars.path_out, 100, file);
           cc = strchr(cfg_pars.path_out, '\n');
           if(cc!=NULL) (*cc) = 0;
           cc = strchr(cfg_pars.path_out, '\r');
           if(cc!=NULL) (*cc) = 0;
        }
        if(strcmp(st, "min-max-days-and-probs-of-incubation-period") == 0)
        {
           fscanf(file,"%d", &cfg_pars.n_inc);
           cfg_pars.min_inc = (int *) malloc((size_t) (cfg_pars.n_inc * sizeof(int)));
           cfg_pars.max_inc = (int *) malloc((size_t) (cfg_pars.n_inc * sizeof(int)));
           cfg_pars.prob_inc = (double **) malloc((size_t) (cfg_pars.n_inc * sizeof(double)));
           
           for(i=0; i<cfg_pars.n_inc; i++)
           {
              fscanf(file, "%d  %d", &cfg_pars.min_inc[i], &cfg_pars.max_inc[i]);
              //printf("%d  %d\n", cfg_pars.min_inc[i], cfg_pars.max_inc[i]);
              len = cfg_pars.max_inc[i] - cfg_pars.min_inc[i] + 1;
              cfg_pars.prob_inc[i] = (double *)malloc((size_t) (len * sizeof(double)));
              for(j = 0 ; j < len ; j++)
              {
                 fscanf(file, "%lf", cfg_pars.prob_inc[i] + j);
                 //printf("%lf  ", cfg_pars.prob_inc[i][j]);
              }
              //printf("\n");
           }


              //fscanf(file,"%d  %d", &cfg_pars.min_incubation, &cfg_pars.max_incubation);
              //len = cfg_pars.max_incubation - cfg_pars.min_incubation + 1;
              //cfg_pars.prob_incubation = (double *)malloc((size_t) (len*sizeof(double)));
              //for(i = 0 ; i < len ; i++)
              //   fscanf(file,"%lf", &cfg_pars.prob_incubation[i]);
        }

        if(strcmp(st, "primary-lower-upper-bounds-and-probs-of-infectious-days-relative-to-symptom-onset-day") == 0)
        {
           fscanf(file,"%d", &cfg_pars.n_inf);
           cfg_pars.lower_inf = (int *) malloc((size_t) (cfg_pars.n_inf * sizeof(int)));
           cfg_pars.upper_inf = (int *) malloc((size_t) (cfg_pars.n_inf * sizeof(int)));
           cfg_pars.prob_inf = (double **) malloc((size_t) (cfg_pars.n_inf * sizeof(double)));
           
           for(i=0; i<cfg_pars.n_inf; i++)
           {
              fscanf(file, "%d  %d", &cfg_pars.lower_inf[i], &cfg_pars.upper_inf[i]);
              //printf("%d  %d\n", cfg_pars.lower_inf[i], cfg_pars.upper_inf[i]);
              len = cfg_pars.upper_inf[i] - cfg_pars.lower_inf[i] + 1;
              cfg_pars.prob_inf[i] = (double *) malloc((size_t) (len * sizeof(double)));
              for(j = 0 ; j < len ; j++)
              {
                 fscanf(file, "%lf", cfg_pars.prob_inf[i] + j);
                 //printf("%lf  ", cfg_pars.prob_inf[i][j]);
              }
              //printf("\n");
           }
              //fscanf(file,"%d  %d", &cfg_pars.lower_infectious, &cfg_pars.upper_infectious);
              //len = cfg_pars.upper_infectious - cfg_pars.lower_infectious + 1;
              //cfg_pars.prob_infectious = (double *)malloc((size_t) (len*sizeof(double)));
              //for(i = 0 ; i < len ; i++)
              //   fscanf(file,"%lf", &cfg_pars.prob_infectious[i]);
        }

        // These are the actual lower and upper bounds of the infectious period to be used in analysis, because
        // in reality contact may not last the whole infectious period but exact contact history is not clear.
        // These effective bounds are not used in calculaing SAR and R

        if(strcmp(st, "number-of-c2p-transmission-probabilities") == 0)
           fscanf(file,"%d", &cfg_pars.n_b_mode);

        if(strcmp(st, "number-of-p2p-transmission-probabilities") == 0)
           fscanf(file,"%d", &cfg_pars.n_p_mode);

        if(strcmp(st, "number-of-pathogenicity-groups") == 0)
           fscanf(file,"%d", &cfg_pars.n_u_mode);

        if(strcmp(st, "number-of-preseason-immunity-groups") == 0)
           fscanf(file,"%d", &cfg_pars.n_q_mode);

        if(strcmp(st, "number-of-time-independent-covariates") == 0)
           fscanf(file,"%d", &cfg_pars.n_time_ind_covariate);

        if(strcmp(st, "number-of-time-dependent-covariates") == 0)
           fscanf(file,"%d", &cfg_pars.n_time_dep_covariate);
           
        if(cfg_pars.n_time_dep_covariate > 0)
        {
           if(strcmp(st, "illness-period-as-a-time-dependent-covariate-for-infectivity") == 0)
           {
              fscanf(file,"%d:", &cfg_pars.illness_as_covariate);
              if(cfg_pars.illness_as_covariate == 1)
              { 
                 fscanf(file,"%d", &cfg_pars.illness_covariate_id);
                 if(cfg_pars.illness_covariate_id > cfg_pars.n_time_dep_covariate)
                 {
                    printf("illness_as_covariate=%d illness_covariate_id=%d\n",  cfg_pars.illness_as_covariate, cfg_pars.illness_covariate_id);
                    printf("Please input a valid time dependent covariate for illness indicator\n");
                    exit(0);
                 }
              }   
           }      
        }
           
        if(strcmp(st, "covariates-affecting-susceptibility-for-c2p-transmission") == 0)
        {
           fscanf(file,"%d:", &cfg_pars.n_c2p_covariate);
           if(cfg_pars.n_c2p_covariate > 0)
           {
              cfg_pars.c2p_covariate = (int *)malloc((size_t) (cfg_pars.n_c2p_covariate * sizeof(int)));
              for(i=0; i<cfg_pars.n_c2p_covariate; i++)
 	         fscanf(file,"%d", &cfg_pars.c2p_covariate[i]);
           }
        }

        if(strcmp(st, "covariates-affecting-susceptibility-for-p2p-transmission") == 0)
        {
           fscanf(file,"%d:", &cfg_pars.n_sus_p2p_covariate);
           if(cfg_pars.n_sus_p2p_covariate > 0)
           {
              cfg_pars.sus_p2p_covariate = (int *)malloc((size_t) (cfg_pars.n_sus_p2p_covariate * sizeof(int)));
              for(i=0; i<cfg_pars.n_sus_p2p_covariate; i++)
                 fscanf(file,"%d", &cfg_pars.sus_p2p_covariate[i]);
           }
        }


        if(strcmp(st, "covariates-affecting-infectiousness-for-p2p-transmission") == 0)
        {
           fscanf(file,"%d:", &cfg_pars.n_inf_p2p_covariate);
           if(cfg_pars.n_inf_p2p_covariate > 0)
           {
              cfg_pars.inf_p2p_covariate = (int *)malloc((size_t) (cfg_pars.n_inf_p2p_covariate * sizeof(int)));
              for(i=0; i<cfg_pars.n_inf_p2p_covariate; i++)
                 fscanf(file,"%d", &cfg_pars.inf_p2p_covariate[i]);
           }
        }

        if(strcmp(st, "interactions-for-p2p-transmission") == 0)
        {
           fscanf(file,"%d:", &cfg_pars.n_int_p2p_covariate);
           if(cfg_pars.n_int_p2p_covariate > 0)
           {
              make_2d_array_int(&cfg_pars.interaction, cfg_pars.n_int_p2p_covariate, 2, 0);
              for(i=0; i<cfg_pars.n_int_p2p_covariate; i++)
                 fscanf(file,"%d  %d", &cfg_pars.interaction[i][0], &cfg_pars.interaction[i][1]);
           }
        }

        if(strcmp(st, "covariates-affecting-pathogenicity") == 0)
        {
           fscanf(file, "%d:", &cfg_pars.n_pat_covariate);
           if(cfg_pars.n_pat_covariate > 0)
           {
              cfg_pars.pat_covariate = (int *)malloc((size_t) (cfg_pars.n_pat_covariate * sizeof(int)));
              for(i=0; i<cfg_pars.n_pat_covariate; i++)
                 fscanf(file,"%d", &cfg_pars.pat_covariate[i]);
           }
        }

        if(strcmp(st, "covariates-affecting-preseason-immunity") == 0)
        {
           fscanf(file, "%d:", &cfg_pars.n_imm_covariate);
           if(cfg_pars.n_imm_covariate > 0)
           {
              cfg_pars.imm_covariate = (int *)malloc((size_t) (cfg_pars.n_imm_covariate * sizeof(int)));
              for(i=0; i<cfg_pars.n_imm_covariate; i++)
                 fscanf(file,"%d", &cfg_pars.imm_covariate[i]);
           }
        }

        cfg_pars.n_covariate = cfg_pars.n_time_ind_covariate + cfg_pars.n_time_dep_covariate;
        cfg_pars.n_p2p_covariate = cfg_pars.n_sus_p2p_covariate + cfg_pars.n_inf_p2p_covariate + cfg_pars.n_int_p2p_covariate;
        cfg_pars.n_par = cfg_pars.n_b_mode + cfg_pars.n_p_mode + cfg_pars.n_u_mode + cfg_pars.n_q_mode +
                         cfg_pars.n_c2p_covariate + cfg_pars.n_p2p_covariate + 
                         cfg_pars.n_pat_covariate + cfg_pars.n_imm_covariate;

        if(strcmp(st, "equal-parameters") == 0)
        {
           fscanf(file, "%d:", &cfg_pars.n_par_equiclass);
           if(cfg_pars.n_par_equiclass > 0)
           {
              cfg_pars.par_equiclass = (PAR_EQUICLASS *)malloc((size_t) (cfg_pars.n_par_equiclass * sizeof(PAR_EQUICLASS)));
              for(i=0; i<cfg_pars.n_par_equiclass; i++)
              {
                 cfg_pars.par_equiclass[i].size = 0;
                 cfg_pars.par_equiclass[i].member = NULL;

                 fscanf(file, "%d:", &cfg_pars.par_equiclass[i].size);
                 cfg_pars.par_equiclass[i].member = (int *) malloc((size_t) (cfg_pars.par_equiclass[i].size * sizeof(int)));
                 for(j=0; j<cfg_pars.par_equiclass[i].size; j++)
                 {
                    fscanf(file,"%d", &cfg_pars.par_equiclass[i].member[j]);
                 }
              }
           }
        }

        if(strcmp(st, "fixed-parameters") == 0)
        {
           fscanf(file,"%d:", &cfg_pars.n_par_fixed);
           if(cfg_pars.n_par_fixed > 0)
           {
              cfg_pars.par_fixed_id = (int *)malloc((size_t) (cfg_pars.n_par_fixed * sizeof(int)));
              cfg_pars.par_fixed_value = (double *)malloc((size_t) (cfg_pars.n_par_fixed * sizeof(double)));
              for(i=0; i<cfg_pars.n_par_fixed; i++)
              {
                 fscanf(file,"%d:  %lf", cfg_pars.par_fixed_id + i, cfg_pars.par_fixed_value + i);
                 k = cfg_pars.par_fixed_id[i] - 1;
                 if(k < cfg_pars.n_b_mode + cfg_pars.n_p_mode + cfg_pars.n_u_mode + cfg_pars.n_q_mode)
                    cfg_pars.par_fixed_value[i] = logit(cfg_pars.par_fixed_value[i]);
                 else
                    cfg_pars.par_fixed_value[i] = log(cfg_pars.par_fixed_value[i]);
              }
           }
        }

        //============================================
        // configuration parameters used by simulation
        //============================================

        if(strcmp(st, "perform-simulation") == 0)
           fscanf(file,"%d", &cfg_pars.simulation);

        if(strcmp(st, "proportion-with-ambiguity-about-preimmunity-and-escape-status") == 0)
           fscanf(file,"%lf", &cfg_pars.prop_mix_imm_esc);

        if(cfg_pars.simulation == 1)
        {
           if(strcmp(st, "parameters-for-simulation") == 0)
           {
              cfg_pars.sim_par_effective = (double *)malloc((size_t) (cfg_pars.n_par_equiclass * sizeof(double)));
              for(i=0; i<cfg_pars.n_par_equiclass; i++)
              {
                 fscanf(file,"%s  %lf", string, &cfg_pars.sim_par_effective[i]);
                 k = cfg_pars.par_equiclass[i].member[0] - 1;
                 if(k < cfg_pars.n_b_mode + cfg_pars.n_p_mode + cfg_pars.n_u_mode + cfg_pars.n_q_mode)
                 {
                    cfg_pars.sim_par_effective[i] = logit(cfg_pars.sim_par_effective[i]);
                 }
                 else
                 {
                    cfg_pars.sim_par_effective[i] = log(cfg_pars.sim_par_effective[i]);
                 }
              }
           }   
           if(strcmp(st, "number-of-simulations") == 0)
              fscanf(file,"%d", &cfg_pars.n_simulation);

        }
        if(strcmp(st, "relative-infectivity-of-asymptomatic-case-for-simulation") == 0)
           fscanf(file,"%lf", &cfg_pars.asym_effect_sim);
         

        //============================================
        // configuration parameters used by estimation
        //============================================

        if(strcmp(st, "optimization-choice") == 0)
           fscanf(file,"%d", &cfg_pars.optimization_choice);
         
        //convergence criteria are based on transformed parameters
        if(strcmp(st, "converge-criteria") == 0)
        {
           fscanf(file,"%d:\n", &cfg_pars.converge_criteria_provided);
           cfg_pars.converge_criteria = (double *)malloc((size_t) (cfg_pars.n_par_equiclass * sizeof(double)));
           if(cfg_pars.converge_criteria_provided == 1)
           {
              for(i=0; i<cfg_pars.n_par_equiclass; i++)
              {
                 fscanf(file,"%s  %lf\n", string, &cfg_pars.converge_criteria[i]);
              }
           }
        }   
        // initial estimates are provided in original parameters, need to be transformed
        if(strcmp(st, "initial-estimates") == 0)
        {
           fscanf(file,"%d:%d\n", &cfg_pars.n_ini, &cfg_pars.ini_par_provided);
           make_2d_array_double(&cfg_pars.ini_par_effective, cfg_pars.n_ini, cfg_pars.n_par_equiclass, 0.0);
           if(cfg_pars.ini_par_provided == 1)
           {
              for(i=0; i<cfg_pars.n_ini; i++)
              {
                 for(j=0; j<cfg_pars.n_par_equiclass; j++)
                 {
                    fscanf(file,"%s  %lf", string, &cfg_pars.ini_par_effective[i][j]);
                    k = cfg_pars.par_equiclass[j].member[0] - 1;
                    if(k < cfg_pars.n_b_mode + cfg_pars.n_p_mode + cfg_pars.n_u_mode + cfg_pars.n_q_mode)
                       cfg_pars.ini_par_effective[i][j] = logit(cfg_pars.ini_par_effective[i][j]);
                    else
                       cfg_pars.ini_par_effective[i][j] = log(cfg_pars.ini_par_effective[i][j]);
                 }
              }
           }
        }   

        if(strcmp(st, "search-bounds") == 0)
        {
           fscanf(file,"%d:\n", &cfg_pars.search_bound_provided);
           make_1d_array_double(&cfg_pars.lower_search_bound, cfg_pars.n_par_equiclass, 0.0);
           make_1d_array_double(&cfg_pars.upper_search_bound, cfg_pars.n_par_equiclass, 0.0);
           if(cfg_pars.search_bound_provided == 1)
           {
              for(i=0; i<cfg_pars.n_par_equiclass; i++)
              {
                 fscanf(file,"%s  %lf  %lf", string, &cfg_pars.lower_search_bound[i], &cfg_pars.upper_search_bound[i]);
                 k = cfg_pars.par_equiclass[i].member[0] - 1;
                 if(k < cfg_pars.n_b_mode + cfg_pars.n_p_mode + cfg_pars.n_u_mode + cfg_pars.n_q_mode)
                 {
                    cfg_pars.lower_search_bound[i] = logit(cfg_pars.lower_search_bound[i]);
                    cfg_pars.upper_search_bound[i] = logit(cfg_pars.upper_search_bound[i]);
                 }
                 else
                 {
                    cfg_pars.lower_search_bound[i] = log(cfg_pars.lower_search_bound[i]);
                    cfg_pars.upper_search_bound[i] = log(cfg_pars.upper_search_bound[i]);
                 }
              }
           }
        } 

        if(strcmp(st, "perform-EM-algorithm") == 0)
           fscanf(file,"%d", &cfg_pars.EM);
        
        if(cfg_pars.EM == 1)
        {
           if(strcmp(st, "min-number-of-possible-status-to-use-mcem") == 0)
              fscanf(file,"%d", &cfg_pars.min_size_MCEM);

           if(strcmp(st, "use-community-specific-weighting") == 0)
              fscanf(file,"%d", &cfg_pars.community_specific_weighting);

           if(strcmp(st, "number-of-base-mcmc-samples") == 0)
              fscanf(file,"%d", &cfg_pars.n_base_sampling);
           
           if(strcmp(st, "number-of-burnin-mcmc-samples") == 0)
              fscanf(file,"%d", &cfg_pars.n_burnin_sampling);

           if(strcmp(st, "number-of-burnin-mcmc-iterations") == 0)
              fscanf(file,"%d", &cfg_pars.n_burnin_iter);

           if(strcmp(st, "number-of-samplings-for-mc-error") == 0)
              fscanf(file,"%d", &cfg_pars.n_sampling_for_mce);

           if(strcmp(st, "use-bootstrap-for-mc-error") == 0)
              fscanf(file,"%d", &cfg_pars.use_bootstrap_for_mce);

           if(strcmp(st, "do-not-calculate-average-variance-for-mc-error") == 0)
              fscanf(file,"%d", &cfg_pars.skip_Evar_for_mce);

           if(strcmp(st, "check-missingness") == 0)
              fscanf(file,"%d", &cfg_pars.check_missingness);
        
           if(strcmp(st, "check-mixing") == 0)
              fscanf(file,"%d", &cfg_pars.check_mixing);
        
           if(strcmp(st, "check-runtime") == 0)
              fscanf(file,"%d", &cfg_pars.check_runtime);
        
        }
           
        if(strcmp(st, "relative-infectivity-of-asymptomatic-case-for-estimation") == 0)
           fscanf(file,"%lf", &cfg_pars.asym_effect_est);
         
        //if(strcmp(st, "members-share-common-c2p-contact-history-across-communities") == 0)
        //   fscanf(file,"%d", &cfg_pars.common_c2p_contact_history_across_community);
        
        if(strcmp(st, "members-share-common-contact-history-within-communities") == 0)
           fscanf(file,"%d", &cfg_pars.common_contact_history_within_community);
        
        if(strcmp(st, "automatically-generate-c2p-contact-file") == 0)
           fscanf(file,"%d", &cfg_pars.generate_c2p_contact);

        if(strcmp(st, "automatically-generate-p2p-contact-file") == 0)
           fscanf(file,"%d", &cfg_pars.generate_p2p_contact);
        
        if(strcmp(st, "use-c2p-offset") == 0)
           fscanf(file,"%d", &cfg_pars.c2p_offset);
        
        if(strcmp(st, "use-p2p-offset") == 0)
           fscanf(file,"%d", &cfg_pars.p2p_offset);
        
        if(strcmp(st, "adjust-for-selection-bias") == 0)
           fscanf(file,"%d", &cfg_pars.adjust_for_left_truncation);
         
        if(strcmp(st, "adjust-for-right-censoring") == 0)
           fscanf(file,"%d", &cfg_pars.adjust_for_right_censoring);
         
        if(strcmp(st, "prefix-index-cases") == 0)
           fscanf(file,"%d", &cfg_pars.preset_index);
           
        if(strcmp(st, "epidemic-duration-for-calculating-CPI") == 0)
           fscanf(file,"%d", &cfg_pars.CPI_duration);
         
        if(strcmp(st, "effective-lower-upper-bounds-of-infectious-days-relative-to-symptom-onset-day") == 0)
        {
           fscanf(file, "%d:", &cfg_pars.effective_bounds_provided);
           cfg_pars.effective_lower_infectious = (int *)malloc((size_t) (cfg_pars.n_p_mode*sizeof(int)));
           cfg_pars.effective_upper_infectious = (int *)malloc((size_t) (cfg_pars.n_p_mode*sizeof(int)));
           if(cfg_pars.n_p_mode > 0)
           {
              if(cfg_pars.effective_bounds_provided == 1)
              {
                 for(i=0; i<cfg_pars.n_p_mode; i++)
                 {
                    fscanf(file,"%d %d", cfg_pars.effective_lower_infectious + i, 
                                       cfg_pars.effective_upper_infectious + i);
                 }
              }
              else if(cfg_pars.effective_bounds_provided == 0)
              {
                 for(i=0; i<cfg_pars.n_p_mode; i++)
                 {
                    cfg_pars.effective_lower_infectious[i] =  cfg_pars.lower_infectious;
                    cfg_pars.effective_upper_infectious[i] =  cfg_pars.upper_infectious;
                 }
              }
           }
        }

        if(strcmp(st, "covariates-for-calculating-SAR-provided") == 0)
        {
           if(cfg_pars.n_p_mode > 0)
           {
              fscanf(file,"%d:", &cfg_pars.SAR_covariate_provided);
              if(cfg_pars.SAR_covariate_provided > 0)
              {
                 if(cfg_pars.n_time_ind_covariate > 0)  
                 {
                    // for time-independent covariates, the format looks like
                    // sus-ind  0.3  0  1  3.2
                    // inf-ind  0.1  1  0  6.7
                    cfg_pars.SAR_sus_time_ind_covariate = (double *)malloc((size_t) (cfg_pars.n_time_ind_covariate * sizeof(double)));
                    cfg_pars.SAR_inf_time_ind_covariate = (double *)malloc((size_t) (cfg_pars.n_time_ind_covariate * sizeof(double)));
                    fscanf(file,"%s", string);
                    for(i=0; i<cfg_pars.n_time_ind_covariate; i++)
 	                    fscanf(file,"%lf", &cfg_pars.SAR_sus_time_ind_covariate[i]);
                    fscanf(file,"%s", string);
                    for(i=0; i<cfg_pars.n_time_ind_covariate; i++)
 	                    fscanf(file,"%lf", &cfg_pars.SAR_inf_time_ind_covariate[i]);
 	              }      
                 if(cfg_pars.n_time_dep_covariate > 0)  
                 {
                    // for time-dependent covariates, suppose infectious period is day -5 to day 13, with day 0 for symptom onset.
                    // suppose only 1 such covariate, and takes value 0 for the whole period for susceptible, and take value 0 for days -5 to 0
                    // and 1 for days 1 to 13. the format looks like
                    // -5  13
                    // sus-dep 1
                    // -5  13  0
                    // inf-dep 2 
                    // -5  0  0
                    // 1  13  1
                    value =  (double *)malloc((size_t) (cfg_pars.n_time_dep_covariate * sizeof(double)));
                    fscanf(file,"%d  %d", &cfg_pars.SAR_time_dep_lower, &cfg_pars.SAR_time_dep_upper);
                    len = cfg_pars.SAR_time_dep_upper - cfg_pars.SAR_time_dep_lower + 1;
                    initialize_matrix(&cfg_pars.SAR_sus_time_dep_covariate);
                    inflate_matrix(&cfg_pars.SAR_sus_time_dep_covariate, len, cfg_pars.n_time_dep_covariate, 0.0); 
                    initialize_matrix(&cfg_pars.SAR_inf_time_dep_covariate);
                    inflate_matrix(&cfg_pars.SAR_inf_time_dep_covariate, len, cfg_pars.n_time_dep_covariate, 0.0);                
                    // read in time-dependent covariates for susceptible 
                    fscanf(file, "%s  %d", string, &m); 
                    //printf("%s  %d\n", string, m);
                    for(i=0; i<m; i++)
                    {
                       fscanf(file,"%d  %d", &start_day, &stop_day);
                       //printf("start=%d  stop=%d:\n", start_day, stop_day);
                       for(j=0; j<cfg_pars.n_time_dep_covariate; j++)
                       {
 	                       fscanf(file,"%lf", &value[j]);
 	                       //printf("%e ", value[j]);
 	                    }   
 	                    //printf("\n"); 
                       for(t=start_day; t<=stop_day; t++)
                       {
                          k = t - cfg_pars.SAR_time_dep_lower;
 	                       for(l=0; l<cfg_pars.n_time_dep_covariate; l++) 
                             cfg_pars.SAR_sus_time_dep_covariate.data[k][l] = value[l];
                       } 
                    }   
                    // read in time-dependent covariates for infectious     
                    fscanf(file, "%s  %d", string, &m); 
                    //printf("%s  %d\n", string, m);
                    for(i=0; i<m; i++)
                    {
                       fscanf(file,"%d  %d", &start_day, &stop_day);
                       //printf("start=%d  stop=%d:\n", start_day, stop_day);
                       for(j=0; j<cfg_pars.n_time_dep_covariate; j++)
                       {
 	                       fscanf(file,"%lf", &value[j]);
 	                       //printf("%e ", value[j]);
 	                    }
                       //printf("\n");   
                       for(t=start_day; t<=stop_day; t++)
                       {
                          k = t - cfg_pars.SAR_time_dep_lower;
 	                       for(l=0; l<cfg_pars.n_time_dep_covariate; l++) 
                             cfg_pars.SAR_inf_time_dep_covariate.data[k][l] = value[l];
                       } 
                    } 
                    free(value);       
 	              }   
              }      
           }
        }
                 
        if(strcmp(st, "multiplier-for-calculating-R0") == 0)
        {
           if(cfg_pars.n_p_mode > 0)
           {
              fscanf(file,"%d:", &cfg_pars.n_R0_multiplier);
              if(cfg_pars.n_R0_multiplier > 0)
              {
                 n = cfg_pars.n_p_mode * cfg_pars.n_R0_multiplier;
                 cfg_pars.R0_multiplier = (double *)malloc((size_t) (n * sizeof(double)));
                 cfg_pars.R0_multiplier_var = (double *)malloc((size_t) (n * sizeof(double)));
                 for(i=0; i<n; i++)
                 {
                    fscanf(file,"%lf %lf", &cfg_pars.R0_multiplier[i], &cfg_pars.R0_multiplier_var[i]);
                    //printf("R0 multiplier=%lf var=%f\n", cfg_pars.R0_multiplier[i], cfg_pars.R0_multiplier_var[i]);
                 }
              }
           }
        }   
        
        if(strcmp(st, "serial-division-of-epidemic-for-calculating-time-varying-R0") == 0)
        {
           if(cfg_pars.n_p_mode > 0)
           {
              fscanf(file,"%d:", &cfg_pars.R0_divide_by_time);
              if(cfg_pars.R0_divide_by_time == 1)
              {
                  if(cfg_pars.n_p_mode == 2)
                     fscanf(file,"%d %d", &cfg_pars.R0_divide_start, &cfg_pars.R0_divide_stop);
                  else
                     fscanf(file,"%d %d %d", &cfg_pars.R0_divide_start, &cfg_pars.R0_divide_stop, &cfg_pars.R0_window_size);   
                  //printf("R0 multiplier=%lf var=%f\n", cfg_pars.R0_multiplier[i], cfg_pars.R0_multiplier_var[i]);
              }
           }
        }   

        if(strcmp(st, "goodness-of-fit") == 0)
           fscanf(file,"%d", &cfg_pars.goodness_of_fit);
        
        if(strcmp(st, "perform-statistical-test") == 0)
           fscanf(file,"%d", &cfg_pars.stat_test);
     
        if(strcmp(st, "do-not-estimate-variance") == 0)
           fscanf(file,"%d", &cfg_pars.skip_variance);

        if(strcmp(st, "do-not-output-estimates") == 0)
           fscanf(file,"%d", &cfg_pars.skip_output);

        if(strcmp(st, "simplify-output") == 0)
           fscanf(file,"%d", &cfg_pars.simplify_output);

        if(strcmp(st, "simplify-output-SAR") == 0)
           fscanf(file,"%d", &cfg_pars.simplify_output_SAR);
        
        if(strcmp(st, "simplify-output-R0") == 0)
           fscanf(file,"%d", &cfg_pars.simplify_output_R0);
         
        if(strcmp(st, "run-transtat-silently") == 0)
           fscanf(file,"%d", &cfg_pars.silent_run);

        if(strcmp(st, "write-error-log") == 0)
           fscanf(file,"%d", &cfg_pars.write_error_log);

     }
  }
  return(cfg_pars);
}


int write_cfg_pars(FILE *file, CFG_PARS *cfg_pars)
{
  int i, j, k, l, m, n, len, n_covariate;


  fprintf(file, "# min-max-days-and-probs-of-incubation-period\n");
  fprintf(file, "%d  %d\n", cfg_pars->min_incubation, cfg_pars->max_incubation);
  len = cfg_pars->max_incubation - cfg_pars->min_incubation + 1;
  for(i = 0 ; i < len ; i++)
     fprintf(file, "%5.3lf  ", cfg_pars->prob_incubation[i]);
  
  fprintf(file, "\n\n");
  fprintf(file, "# primary-lower-upper-bounds-and-probs-of-infectious-days-relative-to-symptom-onset-day\n");
  fprintf(file, "%d  %d\n", cfg_pars->lower_infectious, cfg_pars->upper_infectious);

  len = cfg_pars->upper_infectious - cfg_pars->lower_infectious + 1;
  for(i = 0 ; i < len ; i++)
     fprintf(file, "%5.3lf ", cfg_pars->prob_infectious[i]);


  fprintf(file, "\n\n");
  fprintf(file, "# number-of-c2p-transmission-probabilities\n");
  fprintf(file, "%d\n", cfg_pars->n_b_mode);

  fprintf(file, "\n\n");
  fprintf(file, "# number-of-p2p-transmission-probabilities\n");
  fprintf(file, "%d\n", cfg_pars->n_p_mode);

  fprintf(file, "\n\n");
  fprintf(file, "# number-of-pathogenicity-groups\n");
  fprintf(file, "%d\n", cfg_pars->n_u_mode);

  fprintf(file, "\n\n");
  fprintf(file, "# number-of-preseason-immunity-groups\n");
  fprintf(file, "%d\n", cfg_pars->n_q_mode);

  fprintf(file, "\n\n");
  fprintf(file, "# number-of-time-independent-covariates\n");
  fprintf(file, "%d\n", cfg_pars->n_time_ind_covariate);

  fprintf(file, "\n\n");
  fprintf(file, "# number-of-time-dependent-covariates\n");
  fprintf(file, "%d\n", cfg_pars->n_time_dep_covariate);

  fprintf(file, "\n\n");
  fprintf(file, "# covariates-affecting-susceptibility-for-c2p-transmission\n");
  fprintf(file, "%d:", cfg_pars->n_c2p_covariate);
  if(cfg_pars->n_c2p_covariate > 0)
  {
     for(i=0; i<cfg_pars->n_c2p_covariate; i++)
        fprintf(file, "%d  ", cfg_pars->c2p_covariate[i]);
  }

  fprintf(file, "\n\n");
  fprintf(file, "# covariates-affecting-susceptibility-for-p2p-transmission\n");
  fprintf(file, "%d:", cfg_pars->n_sus_p2p_covariate);
  if(cfg_pars->n_sus_p2p_covariate > 0)
  {
     for(i=0; i<cfg_pars->n_sus_p2p_covariate; i++)
        fprintf(file, "%d  ", cfg_pars->sus_p2p_covariate[i]);
  }

  fprintf(file, "\n\n");
  fprintf(file, "# covariates-affecting-infectiousness-for-p2p-transmission\n");
  fprintf(file, "%d:", cfg_pars->n_inf_p2p_covariate);
  if(cfg_pars->n_inf_p2p_covariate > 0)
  {
     for(i=0; i<cfg_pars->n_inf_p2p_covariate; i++)
        fprintf(file, "%d  ", cfg_pars->inf_p2p_covariate[i]);
  }

  fprintf(file, "\n\n");
  fprintf(file, "# interactions-for-p2p-transmission\n");
  fprintf(file, "%d:", cfg_pars->n_int_p2p_covariate);
  if(cfg_pars->n_int_p2p_covariate > 0)
  {
     for(i=0; i<cfg_pars->n_int_p2p_covariate; i++)
        fprintf(file, "%d  ", cfg_pars->interaction[i]);
  }

  fprintf(file, "\n\n");
  fprintf(file, "# covariates-affecting-pathogenicity\n");
  fprintf(file, "%d:", cfg_pars->n_pat_covariate);
  if(cfg_pars->n_pat_covariate > 0)
  {
     for(i=0; i<cfg_pars->n_pat_covariate; i++)
        fprintf(file, "%d  ", cfg_pars->pat_covariate[i]);
  }

  fprintf(file, "\n\n");
  fprintf(file, "# covariates-affecting-preseason-immunity\n");
  fprintf(file, "%d:", cfg_pars->n_imm_covariate);
  if(cfg_pars->n_imm_covariate > 0)
  {
     for(i=0; i<cfg_pars->n_imm_covariate; i++)
        fprintf(file, "%d  ", cfg_pars->imm_covariate[i]);
  }


  fprintf(file, "\n\n");
  fprintf(file, "# equal-parameters\n");
  fprintf(file, "%d:\n", cfg_pars->n_par_equiclass);
  if(cfg_pars->n_par_equiclass > 0)
  {
     for(i=0; i<cfg_pars->n_par_equiclass; i++)
     {
        fprintf(file, "%d:", cfg_pars->par_equiclass[i].size);
        for(j=0; j<cfg_pars->par_equiclass[i].size; j++)
        {
           fprintf(file, "%d  ", cfg_pars->par_equiclass[i].member[j]);
        }
        fprintf(file, "\n"); 
     }
  }

  fprintf(file, "\n");
  fprintf(file, "# fixed-parameters\n");
  fprintf(file, "%d:\n", cfg_pars->n_par_fixed);
  if(cfg_pars->n_par_fixed > 0)
  {
     for(i=0; i<cfg_pars->n_par_fixed; i++)
     {
        k = cfg_pars->par_fixed_id[i] - 1;
        if(k < cfg_pars->n_b_mode + cfg_pars->n_p_mode + cfg_pars->n_u_mode + cfg_pars->n_q_mode)
           fprintf(file, "%d:%e\n", cfg_pars->par_fixed_id[i], inv_logit(cfg_pars->par_fixed_value[i]));
        else
           fprintf(file, "%d:%e\n", cfg_pars->par_fixed_id[i], exp(cfg_pars->par_fixed_value[i]));
     }
  }

  //============================================
  // configuration parameters used by simulation
  //============================================

  fprintf(file, "\n");
  fprintf(file, "# perform-simulation\n");
  fprintf(file, "%d\n", cfg_pars->simulation);

  fprintf(file, "\n");
  fprintf(file, "# proportion-with-ambiguity-about-preimmunity-and-escape-status\n");
  fprintf(file, "%f\n", cfg_pars->prop_mix_imm_esc);

  if(cfg_pars->simulation == 1)
  {
     fprintf(file, "\n");
     fprintf(file, "# parameters-for-simulation\n");
     for(i=0; i<cfg_pars->n_par_equiclass; i++)
     {
        k = cfg_pars->par_equiclass[i].member[0] - 1;
        if(k < cfg_pars->n_b_mode + cfg_pars->n_p_mode + cfg_pars->n_u_mode + cfg_pars->n_q_mode)
        {
           fprintf(file, "class%d  %e\n", i, inv_logit(cfg_pars->sim_par_effective[i]));
        }
        else
        {
           fprintf(file, "class%d  %e\n", i, exp(cfg_pars->sim_par_effective[i]));
        }
     }
  

     fprintf(file, "\n");
     fprintf(file, "# number-of-simulations\n");
     fprintf(file, "%d\n", cfg_pars->n_simulation);

     fprintf(file, "\n");
     fprintf(file, "# relative-infectivity-of-asymptomatic-case-for-simulation\n");
     fprintf(file, "%f\n", cfg_pars->asym_effect_sim);
         
  }

  //============================================
  // configuration parameters used by estimation
  //============================================

  fprintf(file, "\n");
  fprintf(file, "# optimization-choice\n");
  fprintf(file, "%d\n", cfg_pars->optimization_choice);

  fprintf(file, "\n");
  fprintf(file, "# converge-criteria\n");
  fprintf(file, "%d:\n", cfg_pars->converge_criteria_provided);
  if(cfg_pars->converge_criteria_provided == 1)
  {
     for(i=0; i<cfg_pars->n_par_equiclass; i++)
     {
        fprintf(file,"class%d  %e\n", i, cfg_pars->converge_criteria[i]);
     }
  }

  fprintf(file, "\n");
  fprintf(file, "# initial-estimates\n");
  fprintf(file, "%d:%d\n", cfg_pars->n_ini, cfg_pars->ini_par_provided);
  if(cfg_pars->ini_par_provided == 1)
  {
     for(i=0; i<cfg_pars->n_ini; i++)
     {
        for(j=0; j<cfg_pars->n_par_equiclass; j++)
        {
           k = cfg_pars->par_equiclass[j].member[0] - 1;
           if(k < cfg_pars->n_b_mode + cfg_pars->n_p_mode + cfg_pars->n_u_mode + cfg_pars->n_q_mode)
              fprintf(file, "class%d  %e\n", i, inv_logit(cfg_pars->ini_par_effective[i][j]));
           else
              fprintf(file, "class%d  %e\n", i, exp(cfg_pars->ini_par_effective[i][j]));
        }
     }
  }

  fprintf(file, "\n");
  fprintf(file, "# search-bounds\n");
  fprintf(file, "%d:\n", cfg_pars->search_bound_provided);
  if(cfg_pars->converge_criteria_provided == 1)
  {
     for(i=0; i<cfg_pars->n_par_equiclass; i++)
     {
        k = cfg_pars->par_equiclass[i].member[0] - 1;
        if(k < cfg_pars->n_b_mode + cfg_pars->n_p_mode + cfg_pars->n_u_mode + cfg_pars->n_q_mode)
           fprintf(file,"class%d  %e  %e\n", i, inv_logit(cfg_pars->lower_search_bound[i]), inv_logit(cfg_pars->upper_search_bound[i]));
        else
           fprintf(file,"class%d  %e  %e\n", i, exp(cfg_pars->lower_search_bound[i]), exp(cfg_pars->upper_search_bound[i]));
     }
  }

  fprintf(file, "\n");
  fprintf(file, "# perform-EM-algorithm\n");
  fprintf(file, "%d\n", cfg_pars->EM);
        
  if(cfg_pars->EM == 1)
  {
     fprintf(file, "\n");
     fprintf(file, "# min-number-of-possible-status-to-use-mcem\n");
     fprintf(file, "%d\n", cfg_pars->min_size_MCEM);

     fprintf(file, "\n");
     fprintf(file, "# use-community-specific-weighting\n");
     fprintf(file, "%d\n", cfg_pars->community_specific_weighting);

     fprintf(file, "\n");
     fprintf(file, "# number-of-base-mcmc-samples\n");
     fprintf(file, "%d\n", cfg_pars->n_base_sampling);

     fprintf(file, "\n");
     fprintf(file, "# number-of-burnin-mcmc-samples\n");
     fprintf(file, "%d\n", cfg_pars->n_burnin_sampling);

     fprintf(file, "\n");
     fprintf(file, "# number-of-burnin-mcmc-iterations\n");
     fprintf(file, "%d\n", cfg_pars->n_burnin_iter);

     fprintf(file, "\n");
     fprintf(file, "# number-of-samplings-for-mc-error\n");
     fprintf(file, "%d\n", cfg_pars->n_sampling_for_mce);

     fprintf(file, "\n");
     fprintf(file, "# use-bootstrap-for-mc-error\n");
     fprintf(file, "%d\n", cfg_pars->use_bootstrap_for_mce);

     fprintf(file, "\n");
     fprintf(file, "# do-not-calculate-average-variance-for-mc-error\n");
     fprintf(file, "%d\n", cfg_pars->skip_Evar_for_mce);

     fprintf(file, "\n");
     fprintf(file, "# check-missingness\n");
     fprintf(file, "%d\n", cfg_pars->check_missingness);
        
     fprintf(file, "\n");
     fprintf(file, "# check-mixing\n");
     fprintf(file, "%d\n", cfg_pars->check_mixing);
        
     fprintf(file, "\n");
     fprintf(file, "# check_runtime\n");
     fprintf(file, "%d\n", cfg_pars->check_runtime);
        
  }
           
  fprintf(file, "\n");
  fprintf(file, "# relative-infectivity-of-asymptomatic-case-for-estimation\n");
  fprintf(file, "%d\n", cfg_pars->asym_effect_est);
         
  //fprintf(file, "\n");
  //fprintf(file, "# members-share-common-c2p-contact-history-across-communities\n");
  //fprintf(file, "%d\n", cfg_pars->common_c2p_contact_history_across_community);
        
  fprintf(file, "\n");
  fprintf(file, "# members-share-common-contact-history-within-communities\n");
  fprintf(file, "%d\n", cfg_pars->common_contact_history_within_community);
        
  fprintf(file, "\n");
  fprintf(file, "# automatically-generate-c2p-contact-file\n");
  fprintf(file, "%d\n", cfg_pars->generate_c2p_contact);

  fprintf(file, "\n");
  fprintf(file, "# automatically-generate-p2p-contact-file\n");
  fprintf(file, "%d\n", cfg_pars->generate_p2p_contact);
        
  fprintf(file, "\n");
  fprintf(file, "# use-c2p-offset\n");
  fprintf(file, "%d\n", cfg_pars->c2p_offset);
        
  fprintf(file, "\n");
  fprintf(file, "# use-p2p-offset\n");
  fprintf(file, "%d\n", cfg_pars->p2p_offset);
        
  fprintf(file, "\n");
  fprintf(file, "# adjust-for-selection-bias\n");
  fprintf(file, "%d\n", cfg_pars->adjust_for_left_truncation);
         
  fprintf(file, "\n");
  fprintf(file, "# adjust-for-right-censoring\n");
  fprintf(file, "%d\n", cfg_pars->adjust_for_right_censoring);
         
  fprintf(file, "\n");
  fprintf(file, "# epidemic-duration-for-calculating-CPI\n");
  fprintf(file, "%d\n", cfg_pars->CPI_duration);
         
  fprintf(file, "\n");
  fprintf(file, "# effective-lower-upper-bounds-of-infectious-days-relative-to-symptom-onset-day\n");
  fprintf(file, "%d:\n", cfg_pars->effective_bounds_provided);
  if(cfg_pars->n_p_mode > 0)
  {
     if(cfg_pars->effective_bounds_provided == 1)
     {
        for(i=0; i<cfg_pars->n_p_mode; i++)
        {
           fprintf(file, "%d %d\n", cfg_pars->effective_lower_infectious[i], 
                                    cfg_pars->effective_upper_infectious[i]);
        }
     }
  }
  fprintf(file, "\n");
  fprintf(file, "# multiplier-for-calculating-R0\n");
  fprintf(file, "%d:\n", cfg_pars->n_R0_multiplier);
  if(cfg_pars->n_p_mode > 0)
  {
     if(cfg_pars->n_R0_multiplier > 0)
     {
        n = cfg_pars->n_p_mode * cfg_pars->n_R0_multiplier;
        for(i=0; i<n; i++)
        {
           fprintf(file, "%e %e\n", cfg_pars->R0_multiplier[i], cfg_pars->R0_multiplier_var[i]);
        }
     }
  }   
        

  fprintf(file, "\n");
  fprintf(file, "# goodness-of-fit\n");
  fprintf(file, "%d\n", cfg_pars->goodness_of_fit);
        
  fprintf(file, "\n");
  fprintf(file, "# perform-statistical-test\n");
  fprintf(file, "%d\n", cfg_pars->stat_test);
     
  fprintf(file, "\n");
  fprintf(file, "# do-not-estimate-variance\n");
  fprintf(file, "%d\n", cfg_pars->skip_variance);

  fprintf(file, "\n");
  fprintf(file, "# do-not-output-estimates\n");
  fprintf(file, "%d\n", cfg_pars->skip_output);

  fprintf(file, "\n");
  fprintf(file, "# run-transtat-silently\n");
  fprintf(file, "%d\n", cfg_pars->silent_run);

  

  return(0);
}
