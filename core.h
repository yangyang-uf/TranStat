////// Author: Yang Yang
////// Department of Biostatistics and Emerging Pathogens Institute
////// University of Florida


int core(int id_inc, int id_inf, int id_time)
{
  FILE *file;    
  PEOPLE *person, *member;
  int h, i, j, k, l, n, m, r, t;
  int len_inc, n_iter, n_sym_inf, n_asym_inf;
  int ignore, found, contact_mode; 
  int estimation_error_type, test_error_type;
  int n_par, n_b_mode, n_p_mode, n_u_mode,  n_q_mode, n_c2p_covariate, n_p2p_covariate;
  int n_pat_covariate, n_imm_covariate;
  int n_covariate, n_time_ind_covariate, n_time_dep_covariate;
  int n_sus_p2p_covariate, n_inf_p2p_covariate, n_int_p2p_covariate;                     
  int n_record, possible_pre_immune, possible_escape;
  int possible_sym, possible_sym_start, possible_sym_stop;
  int possible_asym, possible_asym_start, possible_asym_stop;
  int start_day, stop_day;
  int n_idx, n_sym, n_asym, n_esc;
  char file_name[100];

  double se_logit_SAR, var_SAR, var_R0;
  double pr, mean_sym, infective_prob, p_value, p_value_OR;
  double s, log_likelihood, pvalue, offset, multiplier, factor1;
  double covariate_effect, logit_f, f, ff;
  double *est, *p2p_covariate, *der, *der1, *der2, *value;
  double p_esc, log_p_esc_der;

  CONTACT *ptr_contact, *sus_ptr_contact, *inf_ptr_contact, *ptr1_contact, *ptr2_contact;
  STATE *ptr_stat, *ptr2_stat;
  RISK *ptr_risk, *ptr2_risk;
  RISK_CLASS *ptr_class, *ptr2_class;
  INTEGER_CHAIN *ptr_integer, *ptr2_integer;

  MATRIX var, var_logit, der_mat;
  time_t t1, t2;
  //printf("argc=%d  argv[0]=%d\n", argc, argv[0]);
  t1 = time(NULL);

  // Do some initialization
  n_b_mode = cfg_pars.n_b_mode;//number of c2p modes
  n_p_mode = cfg_pars.n_p_mode;//number of p2p modes
  n_u_mode = cfg_pars.n_u_mode;//number of pre-immunity groups
  n_q_mode = cfg_pars.n_q_mode;//number of pre-immunity groups
  n_c2p_covariate = cfg_pars.n_c2p_covariate;//number of c2p covariates
  n_sus_p2p_covariate = cfg_pars.n_sus_p2p_covariate;//number of susceptible covariates for p2p
  n_inf_p2p_covariate = cfg_pars.n_inf_p2p_covariate;//number of infective covariates for p2p
  n_int_p2p_covariate = cfg_pars.n_int_p2p_covariate;//number of interaction terms for p2p 
  n_p2p_covariate = cfg_pars.n_p2p_covariate;//total number of covariates for p2p
  n_imm_covariate = cfg_pars.n_imm_covariate;//total number of covariates for pre-immunity
  n_pat_covariate = cfg_pars.n_pat_covariate;//total number of covariates for pre-immunity
  n_par = cfg_pars.n_par;//total number of parameters
  n_time_ind_covariate = cfg_pars.n_time_ind_covariate;//number of time-independent covariates
  n_time_dep_covariate = cfg_pars.n_time_dep_covariate;//number of time-dependent covariates
  
  initialize_matrix(&var);
  initialize_matrix(&var_logit);
  inflate_matrix(&var, n_par, n_par, 0.0);
  inflate_matrix(&var_logit, n_par, n_par, 0.0);
  initialize_matrix(&der_mat);
  inflate_matrix(&der_mat, n_p_mode, n_par, 0.0);

  importance_weight = NULL;
  pdf_incubation = sdf_incubation = NULL;

  l = cfg_pars.max_incubation - cfg_pars.min_incubation + 1;
  pdf_incubation = (double *)malloc((size_t) (l * sizeof(double)));
  sdf_incubation = (double *)malloc((size_t) (l * sizeof(double)));

  // probability density function of incubation period,
  // pdf_incubation[i] = Pr(incubation duration = i + min_incubation)
  for(t=cfg_pars.min_incubation; t<=cfg_pars.max_incubation; t++)
  {
     r = t - cfg_pars.min_incubation;
     pdf_incubation[r] = cfg_pars.prob_incubation[r];
  }
  // Survival function of incubation period,
  // sdf_incubation[i] = Pr(incubation duration > i + min_incubation)
  for(t=cfg_pars.min_incubation; t<=cfg_pars.max_incubation; t++)
  {
     r = t - cfg_pars.min_incubation;
     if(t == cfg_pars.min_incubation)  sdf_incubation[r] = 1.0 - pdf_incubation[r];
     else if(t == cfg_pars.max_incubation)  sdf_incubation[r] = 0.0;
     else sdf_incubation[r] = sdf_incubation[r-1] - pdf_incubation[r];
     //printf("t=%d  r=%d  pdf=%e sdf=%e\n", t, r, pdf_incubation[r], sdf_incubation[r]);
  }
  // As we allow infectiousness before symptom onset, for each day t before symptom onset, 
  // we need to account for the probability that incubation period is at least day_ill - t
  for(t=cfg_pars.lower_infectious; t<=cfg_pars.upper_infectious; t++)
  {
     r = t - cfg_pars.lower_infectious;
     if(t < 0)
     {
        len_inc = -t;
        if(len_inc >= cfg_pars.min_incubation)
        {
           if(len_inc <= cfg_pars.max_incubation)
              cfg_pars.prob_infectious[r] *= (sdf_incubation[len_inc - cfg_pars.min_incubation] + pdf_incubation[len_inc - cfg_pars.min_incubation]); 
           else  
              cfg_pars.prob_infectious[r] = 0;
        }
     }
  }
  /*for(t=cfg_pars.lower_infectious; t<=cfg_pars.upper_infectious; t++)
  {
     r = t - cfg_pars.lower_infectious;
     printf("day relative to onset=%d  prob_inf=%f\n", t, cfg_pars.prob_infectious[r]);
  }
  goto end;*/
  people = NULL;
  community = NULL;

  est = der = der1 = der2 = value = p2p_covariate = NULL;
  /*printf("input path=%s\n", cfg_pars.path_in);
  printf("output path=%s\n", cfg_pars.path_out);
  printf("n_p_mode=%d\n", cfg_pars.n_p_mode);
  printf("n_q_mode=%d\n", cfg_pars.n_q_mode);
  printf("n_par=%d\n", cfg_pars.n_par);
  printf("n_par_effective=%d\n", cfg_pars.n_par_equiclass);
  printf("simulation=%d\n", cfg_pars.simulation);
  if(cfg_pars.simulation == 1)
  {
     printf("n_simulation=%d\n", cfg_pars.n_simulation);
     printf("asym effect=%f(sim) %f(est)\n", cfg_pars.asym_effect_sim, cfg_pars.asym_effect_est);
     for(j=0; j<cfg_pars.n_par_equiclass; j++)
        printf("%f  ", cfg_pars.sim_par_effective[j]);
     printf("\n");
 }   
  printf("EM=%d\n", cfg_pars.EM);
  printf("n_base_sampling=%d\n", cfg_pars.n_base_sampling);
  printf("n_burnin_sampling=%d\n", cfg_pars.n_burnin_sampling);
  printf("n_burnin_iter=%d\n", cfg_pars.n_burnin_iter);
  
  printf("optimization_choice=%d\n", cfg_pars.optimization_choice);
  printf("n_ini=%d\n", cfg_pars.n_ini);
  for(i=0; i<cfg_pars.n_ini; i++)
  {
     for(j=0; j<cfg_pars.n_par_equiclass; j++)
        printf("%f  ", cfg_pars.ini_par_effective[i][j]);
     printf("\n");
  }
  printf("search_bound_provided=%d\n", cfg_pars.search_bound_provided);
  for(i=0; i<cfg_pars.n_par_equiclass; i++)
  {
     printf("%f  %f\n ", cfg_pars.lower_search_bound[i], cfg_pars.upper_search_bound[i]);
  }
  printf("%d  %d\n ", cfg_pars.effective_lower_infectious[0], cfg_pars.effective_upper_infectious[0]);
  printf("skip_variance=%d  skip_output=%d\n", cfg_pars.skip_variance, cfg_pars.skip_output);
  
  printf("SAR calculation: time range for time-dependent covariates input: lower=%d, upper=%d\n", cfg_pars.SAR_time_dep_lower, cfg_pars.SAR_time_dep_upper);
  printf("SAR calculation: time-dependent covariates of susceptible:\n");
  fmprintf(&cfg_pars.SAR_sus_time_dep_covariate);
  printf("SAR calculation: time-dependent covariates of infectious:\n");
  fmprintf(&cfg_pars.SAR_inf_time_dep_covariate);
  goto end;
  */
  
  if(n_par > 0)
  {
     est = (double *) malloc((size_t) n_par * sizeof(double));
     der = (double *) malloc((size_t) n_par * sizeof(double));
     der1 = (double *) malloc((size_t) n_par * sizeof(double));
     der2 = (double *) malloc((size_t) n_par * sizeof(double));
  }

  if(n_time_ind_covariate > 0 || n_time_dep_covariate > 0)
  {
     k = max(n_time_ind_covariate, n_time_dep_covariate);
     value = (double *) malloc((size_t) k * sizeof(double));
  }

  if(n_p2p_covariate > 0)  
  {
     p2p_covariate = (double *)malloc((size_t) (n_p2p_covariate * sizeof(double)));
  }
  people = NULL;  community = NULL;
  /************************************************************************************************************/
  /************************************************************************************************************/
  /*** construct population structure ***/
  /************************************************************************************************************/
  /************************************************************************************************************/
  if(cfg_pars.silent_run == 0)
  {
     printf("###################################################################\n");
     printf("Incubation Period: %d, Infectious Period: %d, Division Time: %d\n", id_inc, id_inf, id_time);
     printf("###################################################################\n");
  }   
  /* This opens the population input file */
  sprintf(file_name, "%spop.dat", cfg_pars.path_in);
  if((file = fopen(file_name, "r")) == NULL)
  {
     printf("No valid population file\n");
     goto end;
  } 
    
  /* determine population file size */
  p_size = 0;
  p_size = get_size(file);
  if(p_size <= 0)
  {
     printf("The number of population should be larger than 0\n");
     fclose(file);
     goto end;
  }   
  else
  {
     rewind(file); 

     /* Allocate space for the people array. */
     people = (PEOPLE *)malloc((size_t) p_size * sizeof(PEOPLE));
     n_idx = n_sym_inf = n_asym_inf = 0;
     for(i=0,person=people; i < p_size; i++,person++)
     {
        fscanf(file, "%d", &person->id); 
        fscanf(file, "%d", &person->community); 
        fscanf(file, "%d", &person->pre_immune); 
        fscanf(file, "%d", &person->infection); 
        fscanf(file, "%d", &person->symptom); 
        fscanf(file, "%d", &person->day_ill); 
        fscanf(file, "%d", &person->exit); 
        fscanf(file, "%d", &person->day_exit); 
        fscanf(file, "%d", &person->idx); 
        fscanf(file, "%d", &person->u_mode); 
        fscanf(file, "%d", &person->q_mode); 
        fscanf(file, "%lf", &person->weight); 
        fscanf(file, "%d", &person->ignore); 
        //printf("%d  %d  %d\n", person->id, person->community, person->day_ill); 
        //if(cfg_pars.simulation == 1)
        //{
        //   person->infection = person->symptom = person->idx = person->exit = person->ignore = 0;
        //}  
        if(cfg_pars.EM == 0)
        {
           if(person->infection == 1 && person->symptom == 0)
           {
              if(person->day_ill == MISSING)
              {
                 printf("Asymptomatic person %d has no illness (infectiousness) onset day\n", i);
                 printf("Either add illness(infectiousness) onset day for all asymptomatic infections, or turn on the EM option.\n");
                 goto end;
              }
              else
                 person->symptom = 1;
           }  
        }  
           
        person->contact_history = NULL;
        person->risk_history = NULL;
        person->risk_class = NULL;
        person->time_ind_covariate = NULL;
        person->time_dep_covariate = NULL;
        person->imm_covariate = NULL;
        person->pat_covariate = NULL;
        person->possible_states = NULL;
        people[i].day_infection_lower = people[i].day_infection_upper = 
        people[i].day_infective_lower = people[i].day_infective_upper = MISSING;
        if(person->infection == 0)  person->day_ill = MISSING; 
        if(person->exit == 0)  person->day_exit = INFINITY_INTEGER; 
        person->size_possible_states = 0;
        person->final_risk_day = MISSING;
        if(person->idx == 1) n_idx ++; 
        if(person->idx != 1 && person->infection == 1)
        {
           if(person->symptom == 1)  n_sym_inf ++;
           else  n_asym_inf++;
        }
     }
     fclose(file); 
  }
  if(cfg_pars.silent_run == 0)  
     printf("Population size=%d,  number of index cases=%d,  number of non-index infections: symptomatic=%d asymptomatic=%d\n", p_size, n_idx, n_sym_inf, n_asym_inf);
  if(cfg_pars.silent_run == 0 && cfg_pars.EM == 0 && n_asym_inf > 0)
  {
     printf("Asymptomatic infection is found and EM-MCEM is not on. Asymotmatic infections will be analyzed as symptomatic cases.\n");
  }

  sprintf(file_name, "%scommunity.dat", cfg_pars.path_in);
  if((file = fopen(file_name, "r")) == NULL)
  {
     printf("No valid community file\n");
     goto end;
  }   
  n_community = 0;
  n_community = get_size(file);
  if(n_community <= 0)
  {
     printf("The number of community should be larger than 0\n");
     fclose(file);
     goto end;
  }   
  else
  {
     rewind(file); 

     max_epi_duration = 0;
     community = (COMMUNITY *)malloc((size_t) n_community * sizeof(COMMUNITY));
     for(h=0; h < n_community; h++)
     {
        community[h].idx = NULL;
        community[h].size_idx = 0;
        community[h].earliest_idx_day_ill = MISSING;
        community[h].latest_idx_day_ill = MISSING;
        community[h].day_epi_start = MISSING;
        community[h].day_epi_stop = MISSING;
        community[h].epi_duration = 0;
        community[h].member = NULL;
        community[h].size = 0;
        community[h].risk_class = NULL;
        community[h].member_impute = NULL;
        community[h].size_impute = 0;
        community[h].size_possible_states = 1;
        community[h].sample_states = NULL;
        community[h].sample_states_rear = NULL;
        community[h].list_states = NULL;
        community[h].list_states_rear = NULL;
        community[h].contact_history = NULL;
        community[h].ignore = 0;

        fscanf(file, "%d", &community[h].id);  
        fscanf(file, "%d", &community[h].day_epi_start);  //start day of epidemic 
        fscanf(file, "%d", &community[h].day_epi_stop); //stop day of epidemic
        
        community[h].epi_duration = community[h].day_epi_stop - community[h].day_epi_start + 1;

        max_epi_duration = max(max_epi_duration, community[h].epi_duration);
        //0 < day_epi_start <= day_epi_stop must hold 
        if(community[h].day_epi_stop < community[h].day_epi_start)
        {
           printf("Community %d does not have correctly ordered key days: ", h);
           printf("epidemic start=%d  epidemic stop=%d\n",
                   community[h].day_epi_start, community[h].day_epi_stop);    
           fclose(file);
           goto end;
        }
     }
     fclose(file); 
  }
  if(cfg_pars.silent_run == 0)  printf("Number of communities (contact groups) = %d\n", n_community); 
  
  for(i=0; i<p_size; i++)
  if(people[i].ignore == 0)
  {
     h = people[i].community;
     community[h].size++;
  }
  for(h=0; h < n_community; h++)
  {
      if(community[h].size == 0) 
         printf("community %d has no individual!\n", community[h].id); 
      else
         community[h].member = (int *) malloc((size_t) community[h].size * sizeof(int));
  }       

  //add non-ignorable people to the communities they belong to
  for(h=0; h < n_community; h++) community[h].counter = 0;
  for(i=0; i<p_size; i++)
  if(people[i].ignore == 0)
  {
     h = people[i].community;
     community[h].member[community[h].counter] = i;
     community[h].counter++;
     // if infection status are fixed for all individuals, assign the bounds for infectious period.
     // this is convenient for constructing contact history, because we only need to consider contact with infected people. 
     //day_epi_start is considered the first day for any possible exposure.
     //Cases without any possible infection day >= day_epi_start, 
     //they will contribute to exposure of other members, but their own disease outcome
     //should not contribute to the likelihood for tranmission or pathogenicity.
     //day_epi_stop is considered the last observation day. 
     //for case with ILI onset after day_epi_stop, set the cases as none-cases,
     //and user can specify adjustment for right-censoring to counter-balance
     //possible bias caused by this misclassification. 
     if(cfg_pars.EM == 0)
     {
        if(people[i].day_ill > community[h].day_epi_stop)
        {
           people[i].infection = 0;
           people[i].symptom = 0;
           people[i].day_ill = MISSING;
        }
        if(people[i].infection == 1)
        { 
           people[i].day_infection_lower = max(people[i].day_ill - cfg_pars.max_incubation, community[h].day_epi_start);
           people[i].day_infection_upper  = max(people[i].day_ill - cfg_pars.min_incubation, community[h].day_epi_start);
           people[i].day_infective_lower = people[i].day_ill + cfg_pars.lower_infectious;
           people[i].day_infective_upper  = people[i].day_ill + cfg_pars.upper_infectious;
        }
     }   
  }
  /*
  for(i=0; i<5; i++)   
  {
     printf("%d: hh=%d  idx=%d  inf=%d  sym=%d  d_ill=%d  d_exit=%d\n", i, people[i].community, people[i].idx, people[i].infection, 
                                               people[i].symptom, people[i].day_ill, people[i].day_exit);
  }
  
  for(h=260511; h<=26052; h++)   
  {
     printf("comm %d (%d): ", h, community[h].size);
     for(i=0; i<community[h].size; i++)
        printf("%d  ", community[h].member[i]);
     printf("\n");
  }
  exit(0);*/
     
  /*********************************************************
    record covariates for each subject, including both
    time-dependent and time-independent covariates.
   *********************************************************/
   
  for(i=0; i<p_size; i++)
  if(people[i].ignore == 0)
  {
     h = people[i].community;
     if(n_time_ind_covariate > 0)
	       make_1d_array_double(&people[i].time_ind_covariate, n_time_ind_covariate, 0.0);
     if(n_time_dep_covariate > 0)
	       make_2d_array_double(&people[i].time_dep_covariate, community[h].epi_duration, n_time_dep_covariate, 0.0);
  }   
  
  if(n_time_ind_covariate > 0)
  {
     sprintf(file_name, "%stime_ind_covariate.dat", cfg_pars.path_in);
     if((file = fopen(file_name, "r")) == NULL)
     {
        printf("No valid time-independent covariates file\n");
        goto end;
     }   
     n_record = get_size(file);
     rewind(file); 
     printf("n_time_ind_covariate=%d,  n_record=%d\n", n_time_ind_covariate, n_record);
     for(m=0; m<n_record; m++)
     {
        fscanf(file, "%d", &i); 
        for(j=0; j<n_time_ind_covariate; j++)
           fscanf(file, "%lf", &value[j]); 
        if(people[i].ignore == 0)
        {
           for(j=0; j<n_time_ind_covariate; j++)
              people[i].time_ind_covariate[j] = value[j];
        }
     }
     fclose(file);
  }
  
  if(n_time_dep_covariate > 0)
  {
     sprintf(file_name, "%stime_dep_covariate.dat", cfg_pars.path_in);
     if((file = fopen(file_name, "r")) == NULL)
     {
        printf("No valid time-dependent covariates file\n");
        goto end;
     }   
     n_record = get_size(file);
     rewind(file); 
     printf("n_time_dep_covariate=%d,  n_record=%d\n", n_time_dep_covariate, n_record);
     for(m=0; m<n_record; m++)
     {
        /*time-dependent covariates are stored as subject, covariate index, day, value per line*/
        fscanf(file, "%d  %d  %d", &i, &start_day, &stop_day); 
        for(k=0; k<n_time_dep_covariate; k++)
           fscanf(file, "%lf", &value[k]); 
        if(people[i].ignore == 0)
        {
           h = people[i].community;
           for(t=start_day; t<=stop_day; t++)
           if(t >= community[h].day_epi_start && t <= community[h].day_epi_stop)
           {
              r = t - community[h].day_epi_start;
              for(k=0; k<n_time_dep_covariate; k++)
                 people[i].time_dep_covariate[r][k] = value[k];
           }
        } 
     }   
     fclose(file);
  }
  
  for(h=0; h<n_community; h++)
  if(community[h].size > 0)
  {
     if(cfg_pars.common_contact_history_within_community == 1)
     {
        community[h].contact_history = (CONTACT_HISTORY *)malloc((size_t) (community[h].epi_duration * sizeof(CONTACT_HISTORY)));
        for(t=community[h].day_epi_start; t<=community[h].day_epi_stop; t++)
        {
           r = t - community[h].day_epi_start;
           community[h].contact_history[r].c2p_contact = NULL;
           community[h].contact_history[r].c2p_contact_rear = NULL;
           community[h].contact_history[r].p2p_contact = NULL;
           community[h].contact_history[r].p2p_contact_rear = NULL;
           // risk history array is created in each risk class in the create_risk_class() function
        }
     }   
  }

  for(i=0; i<p_size; i++)
  if(people[i].ignore == 0)
  {
     h = people[i].community;
     if(cfg_pars.common_contact_history_within_community == 0)
     {
        people[i].contact_history = (CONTACT_HISTORY *)malloc((size_t) (community[h].epi_duration * sizeof(CONTACT_HISTORY)));
        people[i].risk_history = (RISK_HISTORY *)malloc((size_t) (community[h].epi_duration * sizeof(RISK_HISTORY)));
        for(t=community[h].day_epi_start; t<=community[h].day_epi_stop; t++)
        {
           r = t - community[h].day_epi_start;
           people[i].contact_history[r].c2p_contact = NULL;
           people[i].contact_history[r].c2p_contact_rear = NULL;
           people[i].contact_history[r].p2p_contact = NULL;
           people[i].contact_history[r].p2p_contact_rear = NULL;
           people[i].risk_history[r].c2p_risk = NULL;
           people[i].risk_history[r].c2p_risk_rear = NULL;
           people[i].risk_history[r].p2p_risk = NULL;
           people[i].risk_history[r].p2p_risk_rear = NULL;
        }   
     }
     
     if(n_pat_covariate > 0)  
     {
        make_2d_array_double(&people[i].pat_covariate, n_pat_covariate, community[h].epi_duration, 0.0);
        organize_pat_covariate(people + i);
     }

     if(n_imm_covariate > 0)  
     {
        make_1d_array_double(&people[i].imm_covariate, n_imm_covariate, 0.0);
        organize_imm_covariate(people + i);
     }
  }

  // ad-hoc creation of c2p-contact file
  if(cfg_pars.generate_c2p_contact == 1)
  {
  	  if(n_b_mode != 1)
  	  {
  	     printf("To generate c2p contact history automatically, need n_b_mode=1\n");
  	     exit(0);
     }
  	  printf("Auto-generate c2p contact history\n");
  	  if(cfg_pars.common_contact_history_within_community == 1)
  	  {
        for(h=0; h<n_community; h++)
        {
           if(community[h].size > 0)
              add_c2p_contact_history_to_community(h, community[h].day_epi_start, community[h].day_epi_stop, 0, 0.0);
        }
     }
     else
     {
        for(i=0; i<p_size; i++)
        if(people[i].ignore == 0)
        {
           h = people[i].community;
           for(t = community[h].day_epi_start; t <= community[h].day_epi_stop; t++)
           {
              if(t <= people[i].day_exit)  add_c2p_contact_history(t, people + i, 0, 0.0);
           }
        }
     }
  }
  
  // ad-hoc creation of p2p-contact file
  if(cfg_pars.generate_p2p_contact == 1)
  {
     if(cfg_pars.common_contact_history_within_community != 1 || n_p_mode != 1)
  	  {
  	     printf("To generate p2p contact history automatically, need common_contact_history_within_community= 1 AND n_p_mode=1\n");
  	     exit(0);
     }
  	  printf("Auto-generate p2p contact history\n");
     for(h=0; h<n_community; h++)
     if(community[h].size > 0)
        add_p2p_contact_history_to_community(h, community[h].day_epi_start, community[h].day_epi_stop, 0, 0.0);
  }

  /************************************************************* 
  If susceptible subject i is exposed to common source of some type, 
  record this exposure in the structure people[i].contact_history[r].c2p_contact. 
  **************************************************************/
  /* This opens the special exposure input file */
  if(n_b_mode > 0 && cfg_pars.generate_c2p_contact == 0)
  {
     if(cfg_pars.silent_run == 0)  printf("Reading c2p contact information\n");
     sprintf(file_name, "%sc2p_contact.dat", cfg_pars.path_in);
     if((file = fopen(file_name, "r")) == NULL)
     {
        printf("No valid c2p_contact information file\n");
        goto end;
     }   
  
     n_record = get_size(file);
     rewind(file);
     // when all people, no matter which comunity, share the same contact history,
     // the c2p contact input file has the simplest format.
     /*if(cfg_pars.common_c2p_contact_history_across_community == 1)
     {
        for(m=0; m<n_record; m++)
        {
           fscanf(file, "%d  %d  %d  %lf  %d\n", &start_day, &stop_day, &contact_mode, &offset, &ignore);
           if(cfg_pars.c2p_offset == 0)  offset = 0.0;
           if(ignore != 1)
           {
              for(h=0; h<n_community; h++)
              {
                 if(community[h].size > 0)
                    add_c2p_contact_history_to_community(h, start_day, stop_day, contact_mode, offset);
              }
           }
        } 
     }*/
     if(cfg_pars.common_contact_history_within_community == 1)
     {
        //when people in the same community share the same contact history,
        //we only need to know which type of unknown source on which days for which community. 
        for(m=0; m<n_record; m++)
        {
           fscanf(file, "%d  %d  %d  %d  %lf  %d\n", &h, &start_day, &stop_day, &contact_mode, &offset, &ignore);
           if(cfg_pars.c2p_offset == 0)  offset = 0.0;
           if(start_day > community[h].day_epi_stop || stop_day < community[h].day_epi_start)
           {
              printf("c2p contact dates are all outside the epidemic period of community %d\n", h);
              goto end;
           }
           if(ignore != 1 && community[h].size > 0)
              add_c2p_contact_history_to_community(h, start_day, stop_day, contact_mode, offset);
        } 
     }
     else
     {
        //Individualized contact history
        for(m=0; m<n_record; m++)
        {
           fscanf(file, "%d  %d  %d  %d  %lf  %d\n", &i, &start_day, &stop_day, &contact_mode, &offset, &ignore); 
           if(cfg_pars.c2p_offset == 0)  offset = 0.0;
           if(ignore != 1 && people[i].ignore != 1)
           {
              h = people[i].community;
              if(start_day > community[h].day_epi_stop || stop_day < community[h].day_epi_start)
              {
                 printf("c2p contact dates of person %d are all outside the epidemic period of community %d\n", i, h);
                 goto end;
              }
              for(t = start_day; t <= stop_day; t++)
              {
                 if(t >= community[h].day_epi_start && t <= community[h].day_epi_stop &&  t <= people[i].day_exit)
                 {
                    add_c2p_contact_history(t, people + i, contact_mode, offset);
                 }
              }
           }
        }
     }   
     fclose(file);
  }
  
  
  /*h = 2;
  for(t=community[h].day_epi_start; t<=community[h].day_epi_stop; t++)
  {
     r = t - community[h].day_epi_start;
     ptr_contact = community[h].contact_history[r].c2p_contact;
     if(ptr_contact == NULL)  printf("day=%d: no c2p contact\n", t);
     while(ptr_contact != NULL)
     {
        printf("day=%d  contact_mode:%d  offset:%f \n",
                t, ptr_contact->contact_mode, ptr_contact->offset);
        ptr_contact = ptr_contact->next;
     }
  }*/
  
    
  /*for(i=0; i<p_size; i++)
  if(people[i].id == 1)
  {
     h = people[i].community;
     printf("community %d: day_epi_start=%d,  day_epi_stop=%d\n", h, community[h].day_epi_start, community[h].day_epi_stop);
     for(t=community[h].day_epi_start; t<=community[h].day_epi_stop; t++)
     {
        r = t - community[h].day_epi_start;
        ptr_contact = people[i].contact_history[r].c2p_contact;
        while(ptr_contact != NULL)
        {
           printf("id=%d  day=%d  contact_mode:%d  offset:%f\n",
                  people[i].id, t, ptr_contact->contact_mode,
                  ptr_contact->offset);
           printf("\n");
           ptr_contact = ptr_contact->next;
        }
     }
  }
  exit(0);*/
  
  /************************************************************* 
  If subject i made contact with subject j when i is susceptible and
  j is infectious, record this contact in the structure
  people[i].contact_history[r].p2p_contact. 
  **************************************************************/

  
  if(n_p_mode > 0 && cfg_pars.generate_p2p_contact == 0)
  {
     if(cfg_pars.silent_run == 0)  printf("Reading p2p contact information\n");
     sprintf(file_name, "%sp2p_contact.dat", cfg_pars.path_in);
     if((file = fopen(file_name, "r")) == NULL)
     {
        printf("No valid p2p_contact information file\n");
        goto end;
     }   
  
     n_record = get_size(file);
     rewind(file); 
  
     if(cfg_pars.common_contact_history_within_community == 1)
     {
        //when people in the community share the same contact history,
        //we only need to know which contact for which days. 
        for(m=0; m<n_record; m++)
        {
           fscanf(file, "%d  %d  %d  %d  %lf  %d\n", &h, &start_day, &stop_day, &contact_mode, &offset, &ignore);
           if(cfg_pars.p2p_offset == 0)  offset = 0.0;
           if(start_day > community[h].day_epi_stop || stop_day < community[h].day_epi_start)
           {
              printf("p2p contact dates are all outside the epidemic period of community %d\n", h);
              goto end;
           }
           if(ignore != 1 && community[h].size > 0)
              add_p2p_contact_history_to_community(h, start_day, stop_day, contact_mode, offset);
        }
     }
     else
     {
        //Individualized contact history
        for(m=0; m<n_record; m++)
        {
           fscanf(file, "%d  %d  %d  %d  %d  %lf  %d\n", &start_day, &stop_day, &i, &j, &contact_mode, &offset, &ignore);
           if(cfg_pars.p2p_offset == 0)  offset = 0.0;
           if(ignore != 1 && people[i].ignore == 0 && people[j].ignore == 0 
              && people[i].community == people[j].community && people[i].id != people[j].id)
           {
              h = people[i].community;
              if(start_day > community[h].day_epi_stop || stop_day < community[h].day_epi_start)
              {
                 printf("p2p contact dates between person %d and person %d are all outside the epidemic period of community %d\n", i, j, h);
                 goto end;
              }
              for(t = start_day; t <= stop_day; t++)
              {
                 if(t >= community[h].day_epi_start && t <= community[h].day_epi_stop &&  t <= people[i].day_exit &&  t <= people[j].day_exit)
                 {
                    ptr1_contact = add_p2p_contact_history(t, people + i, people + j, contact_mode, offset);
                    ptr2_contact = add_p2p_contact_history(t, people + j, people + i, contact_mode, offset);
                    ptr1_contact->pair = ptr2_contact;
                    ptr2_contact->pair = ptr1_contact;
                 } 
              }     
           }  
        }
     }
     fclose(file); 
  }
  
  /*h = 2;
  for(t=community[h].day_epi_start; t<=community[h].day_epi_stop; t++)
  {
     r = t - community[h].day_epi_start;
     ptr_contact = community[h].contact_history[r].p2p_contact;
     if(ptr_contact == NULL)  printf("day=%d  no contact\n", t);
     while(ptr_contact != NULL)
     {
        i = ptr_contact->contact_id;
        printf("day=%d  j=%d  onset=%d  lower=%d  upper=%d  contact_mode:%d  offset:%f  covariate:",
               t, i, people[i].day_ill, people[i].day_infective_lower, people[i].day_infective_upper, ptr_contact->contact_mode,
               ptr_contact->offset);
        printf("\n");
        ptr_contact = ptr_contact->next;
     }
  }
  goto end;*/
  
  /*for(i=0; i<p_size; i++)
  if(people[i].id ==994)
  {
     h = people[i].community;
     for(t=community[h].day_epi_start; t<=community[h].day_epi_stop; t++)
     {
        r = t - community[h].day_epi_start;
        ptr_contact = people[i].contact_history[r].p2p_contact;
        while(ptr_contact != NULL)
        {
           printf("i=%d  day=%d  j=%d  contact_mode:%d  offset:%f\n",
                  people[i].id, t, ptr_contact->contact_id, ptr_contact->contact_mode,
                  ptr_contact->offset);
           printf("\n");
           ptr_contact = ptr_contact->next;
        }
     }
  }
  exit(0);*/
  
  create_arrays();
  // Create risk classes for each community.
  // People in the same class share contact history and covariates affecting risk.
  // Such classes expedite the analysis quite a lot.
  create_risk_class();

  if(cfg_pars.EM == 0)
  {
     //if(admissibility() == 1)  goto end;
  }  
     
  //begin iterations for simulation. If this is not a simulation, there will be only one loop.
  n_iter = (cfg_pars.simulation == 1)? cfg_pars.n_simulation : 1;
  printf("n_iter=%d\n", n_iter);
  for(ITER=0; ITER<n_iter; ITER++)
  {
     if(cfg_pars.simulation == 1 && cfg_pars.silent_run >= 0)  printf("\n SIMULATION: %d \n\n", ITER);
     if(cfg_pars.simulation == 1)  simulate(0, NULL);
     if(cfg_pars.output_simulation_data == 1) //output simulated data (only with time-independent covariates) for regression analysis
     {
        sprintf(file_name, "%ssim_pop_%d.txt", cfg_pars.path_out, ITER);
        file = fopen(file_name, "w");        
        for(i=0,person=people; i < p_size; i++,person++)
        {
           fprintf(file, "%d  %d  %d  %d  %d", person->id, person->community, person->idx, person->infection, person->day_ill); 
           for(j=0; j<n_time_ind_covariate; j++)
              fprintf(file, "  %e",  person->time_ind_covariate[j]);
           fprintf(file, "\n");
        }
        fclose(file);      
     }
     // From now on an epidemic has been generated, either from user input files or the simulator

     /*
     //we search the earliest user specified index case in the community, and treat that as
     //the base day for left-truncation adjustment. 
     //Attention: if a case has symptom onset way earlier than other index cases, it would be better
     //not to denote it as an index case in pop.dat, but let TranStat to set it as an index case,
     //so that TranStat does not set that case's symptom onset day as community[h].idx_day_ill.
     for(h=0; h < n_community; h++)
     {
        for(j=0; j<community[h].size; j++)
        {
           i = community[h].member[j];
           if(people[i].idx == 1 && people[i].infection == 1)
           {
              if(community[h].idx_day_ill == MISSING)
                 community[h].idx_day_ill = people[i].day_ill;
              else
                 community[h].idx_day_ill = min(community[h].idx_day_ill, people[i].day_ill);
           }
        }
     }
     


     //People with earlier symptom onset days than community[h].idx_day_ill
     //but not user-specified as index cases will be regarded as index cases
     //automatically, but they do not affect the base day for left-truncation adjustment,
     //because they are not the reason the community is selected. We mark them
     // as index cases simply to exclude them from the likelihood. Like user-specified
     //index cases, they do contribute to the exposure history of secondary cases though.

     for(h=0; h < n_community; h++)
     {
        for(j=0; j<community[h].size; j++)
        {
           i = community[h].member[j];
           if(people[i].idx == 0 && people[i].infection == 1 && people[i].symptom == 1 && community[h].idx_day_ill != MISSING 
              && people[i].day_ill <= community[h].idx_day_ill)
              people[i].idx = 1;
        }
     }   

     */
     
     // In above, we were trying to identify co-primary cases based on illness onset days.
     // For some diseases such as Cholera, infectiousness (shedding) onset can be earlier than illness onset.
     // Before TranStat can accomodate both types of onset, we can use illness onset days as infectiousness onset days,
     // which does not use the true illness onset day or ascertainment day.
     // In such situation, infectiousness onset day could be missing for index cases, and it is impossible 
     // to idenitfy co-primary cases based on uncertain infectiousness onset days of index cases.
     // Therefore, we disable above code, and do not try to identify co-primary cases. If user want some subjects
     // to be treated as coprimary cases in the analyses, these subjects should be denoted as index cases, i.e., set idx=1 in the input file "pop.dat".

     // For cases-ascertained design, determine who are the index cases, and delete communities without index cases.
     // Some infected subjects may have their symptom onsets sampled dynamically.
     // If you don't want to change who are the index cases, set cfg_pars.preset_index=1.
     // If you want to reassign who are index cases according to newly sampled onset times, set cfg_pars.preset_index=0. 
    
     if(cfg_pars.adjust_for_left_truncation == 1)
     {
        for(h=0; h<n_community; h++)   
        {
                 
           //printf("h=%d\n", h);
           //if cfg_pars.preset_index=1, set_index_cases() will only create an index case tracking array and identify the earliest 
           //and latest illness onset dates of index cases. if cfg_pars.preset_index=0, index case status will be assigned according to the sequqence of onsets.  
           set_index_cases(community+h);
           // For simulation with case-ascertained design, families without index cases or with only index cases are excluded from analysis.
           // For real data analysis, if EM-MCEM is not to be used, these families are also exlcuded; but if EM-MCEM is to be used,
           // we cannot exclude these families because index cases could be generated when we sample possible states in these families.
           // Example: a family may have only an asymptomatic infection, and the day_ill will be set to MISSING when pop.dat is read in.
           // TranStat will not detect this peron as an index case unless user set this infection as index case and set preset_index=1,
           // and thus size_idx=0. When EM-MCEM is running, this person will be set as the index case.    
           
           if(!(cfg_pars.preset_index == 0 && cfg_pars.EM == 1))
           {
              if(community[h].size_idx == 0 || community[h].size_idx == community[h].size)
              {
                 printf("community to be ignored:%d  size_idx=%d  size=%d\n", h, community[h].size_idx, community[h].size);
                 for(j=0; j<community[h].size; j++)
                 {
                    i = community[h].member[j];
                    printf("member %d: idx=%d  infection=%d  symptom=%d  day_ill=%d\n", i, people[i].idx, people[i].infection, people[i].symptom, people[i].day_ill);
                 }  
                 exit(0); 
                 community[h].ignore = 1;
                 if(community[h].idx != NULL)  
                 {
                    free(community[h].idx);
                    community[h].idx = NULL;
                 }   
                 for(j=0; j<community[h].size; j++)
                 {
                    i = community[h].member[j];
                    people[i].ignore = 1;
                 }
              }  
           }
        }   
     }
     estimation_error_type = test_error_type = -1;
     //count_cases();
     /*file = fopen("check_sim_pop.dat", "w");
     for(i=0; i<p_size; i++)   
     {
        fprintf(file, "%3d: hh=%3d  idx=%1d  inf=%1d  sym=%1d  d_inf=%5d  d_ill=%5d  exit=%1d  d_exit=%5d\n", 
                i, people[i].community, people[i].idx, people[i].infection, 
                people[i].symptom, people[i].day_infection, people[i].day_ill, people[i].exit, people[i].day_exit);
     }
     fclose(file);*/
     /*file = fopen("check_sim_community.dat", "w");
     for(h=0; h<n_community; h++)   
     {
        fprintf(file, "%3d: size_index=%3d  earliest_idx_day_ill=%5d  latest_idx_day_ill=%5d  day_ill=", h, community[h].size_idx, community[h].earliest_idx_day_ill, community[h].latest_idx_day_ill);
        for(j=0; j<community[h].size_idx; j++)  
        {
           i = community[h].idx[j];
           fprintf(file, "%5d  ", people[i].day_ill);
        }
        fprintf(file, "\n");
     }
     fclose(file);*/
     //printf("statistical test=%d\n", cfg_pars.stat_test);
     if(cfg_pars.stat_test == 1 && cfg_pars.EM == 0)
     {
        p_value = stat_test_cs_hh();
     } 
     
     /*file = fopen("check_risk_class.dat", "w");
     for(h=0; h<n_community; h++)
     {
        fprintf(file, "%d:\n", h);
        ptr_class = community[h].risk_class;
        while(ptr_class != NULL)
        {
           ptr_integer = ptr_class->member;
           while(ptr_integer != NULL)
           {
              fprintf(file, "%d  ", ptr_integer->id);
              ptr_integer = ptr_integer->next;
           }
           fprintf(file, "\n");
           ptr_class = ptr_class->next;
        }
     }   
     fclose(file);
     goto end;*/
        
     n_need_OEM = n_need_MCEM = 0;
     if(cfg_pars.EM == 1)
     {
        if(cfg_pars.simulation == 0)
        {
           sprintf(file_name, "%simpute.dat", cfg_pars.path_in);
           if((file = fopen(file_name, "r")) == NULL)
           {
              printf("No valid imputation information file\n");
              goto end;
           }   
        
           n_record = get_size(file);
           rewind(file); 
           for(m=0; m<n_record; m++)
           {
              fscanf(file, "%d", &i); 
              fscanf(file, "%d", &possible_pre_immune); 
              fscanf(file, "%d", &possible_escape); 
              fscanf(file, "%d  %d  %d", &possible_sym, &possible_sym_start, &possible_sym_stop); 
              fscanf(file, "%d  %d  %d", &possible_asym, &possible_asym_start, &possible_asym_stop); 
              if(people[i].ignore != 1)
              {
                 people[i].size_possible_states = 0;
                 if(possible_pre_immune == 1)  people[i].size_possible_states++;
                 if(possible_escape == 1)  people[i].size_possible_states++;
                 if(possible_sym == 1)  people[i].size_possible_states += (possible_sym_stop - possible_sym_start + 1);
                 if(possible_asym == 1)  people[i].size_possible_states += (possible_asym_stop - possible_asym_start + 1);
                 if(people[i].size_possible_states > 0)
                 {
                    make_2d_array_int(&people[i].possible_states, 2, people[i].size_possible_states, 0);
                    k = 0;
                    if(possible_pre_immune == 1)
                    {
                       people[i].possible_states[0][k] = -INFINITY_INTEGER;
                       k++;
                    }
                    if(possible_escape == 1)
                    {
                       people[i].possible_states[0][k] = INFINITY_INTEGER;
                       k++;
                    }
                    if(possible_sym == 1)
                    {
                       for(j=possible_sym_start; j<=possible_sym_stop; j++)
                       {
                          people[i].possible_states[0][k] = j;
                          people[i].possible_states[1][k] = 1;
                          k++;
                       }
                    }
                    if(possible_asym == 1)
                    {
                       for(j=possible_asym_start; j<=possible_asym_stop; j++)
                       {
                          people[i].possible_states[0][k] = j;
                          people[i].possible_states[1][k] = 0;
                          k++;
                       }
                    }
                 }
                 //if(people[i].id == 1764)  
                 //{    
                 //    print_2d_int(people[i].possible_states, 2, people[i].size_possible_states);
                 //    exit(0);
                 //}    
              }
           }
           fclose(file);
        }
        else
        {
           for(i=0; i<p_size; i++)
           {
              if(people[i].ignore != 1 && people[i].infection == 1 && people[i].symptom == 0)
              {
                 h = people[i].community;
                 possible_asym_start = community[h].day_epi_start + cfg_pars.min_incubation;
                 possible_asym_stop = community[h].day_epi_stop;
                 people[i].size_possible_states = (possible_asym_stop - possible_asym_start + 1);
                 make_2d_array_int(&people[i].possible_states, 2, people[i].size_possible_states, 0);
                 k = 0;
                 for(j=possible_asym_start; j<=possible_asym_stop; j++)
                 {
                    people[i].possible_states[0][k] = j;
                    people[i].possible_states[1][k] = 0;
                    k++;
                 }
              }
              // for 50% of preimmune or escaped people, we assume pre-immunity and escape are not distinguishable
              
              if(people[i].ignore != 1 && people[i].infection == 0 && people[i].symptom == 0 && cfg_pars.n_q_mode > 0)
              {
                 if(runiform(&seed) < cfg_pars.prop_mix_imm_esc)
                 {
                    people[i].size_possible_states = 2;
                    make_2d_array_int(&people[i].possible_states, 2, people[i].size_possible_states, 0);
                    people[i].possible_states[0][0] = -INFINITY_INTEGER;
                    people[i].possible_states[0][1] = INFINITY_INTEGER;
                 }
              }
              //if(people[i].community == 1 | people[i].community == 48)
              //   printf("fam %d  person %d: preimmune=%d  infect=%d  symptom=%d  size_possible_states=%d\n", 
              //           people[i].community, i, people[i].pre_immune, people[i].infection, people[i].symptom, people[i].size_possible_states);

           }
        } //if(cfg_pars.simulation == 0)
        
        //any subject need imputation in this community
        for(h=0; h < n_community; h++)
        {
           for(j=0; j<community[h].size; j++)
           {
              i = community[h].member[j];
              if(people[i].size_possible_states > 0)  
              {
                 community[h].size_impute ++;
                 community[h].size_possible_states *= people[i].size_possible_states;
              } 
           }
           //printf("%d: number of possible stauts=%d\n", h, community[h].size_possible_states);
           if(community[h].size_impute > 0)
           {
              community[h].member_impute = (int *) malloc((size_t) community[h].size_impute * sizeof(int));
              k = 0;
              for(j=0; j<community[h].size; j++)
              {
                 i = community[h].member[j];
                 if(people[i].size_possible_states > 0)  
                 {
                    community[h].member_impute[k++] = i;
                 }
              }
           }
           // if size_impute is <= 1 in all communities, we can do straightforward EM rather than MCEM.
           if(community[h].size_possible_states > 1 && 
              community[h].size_possible_states < cfg_pars.min_size_MCEM)  n_need_OEM++;
           if(community[h].size_possible_states >= cfg_pars.min_size_MCEM)  n_need_MCEM++;
        }
        if(cfg_pars.silent_run == 0)  printf("number of communities need OEM (%d) or MCEM (%d)\n", n_need_OEM, n_need_MCEM);
        
        if(cfg_pars.check_missingness == 1)
        {
           file = fopen("check_missingness.dat", "w");
           for(h=0; h<n_community; h++)
           fprintf(file, "%d  %d\n", h, community[h].size_possible_states);
           fclose(file);
           restore_susceptibility();
           free(importance_weight);
           goto end;
        }

        //day_epi_start is considered the first day for any possible exposure.
        //Cases without any possible infection day >= day_epi_start, 
        //they will contribute to exposure of other members, but their own disease outcome
        //should not contribute to the likelihood for tranmission or pathogenicity.
        //day_epi_stop is considered the last observation day. 
        //for case with ILI onset after day_epi_stop, set the cases as none-cases,
        //and user can specify adjustment for right-censoring to counter-balance
        //possible bias caused by this misclassification. 
        for(i=0; i<p_size; i++)
        if(people[i].ignore == 0 && people[i].infection == 1 && people[i].size_possible_states == 0)
        {
           h = people[i].community;
           if(h < 0 || h >= n_community)
           {
              printf("person %d has an invalid community id %d\n", people[i].id, h);
              goto end;
           }
           // Do not set ignore=1 as this person may still contribute exposure information for susceptibles
           //if((people[i].day_ill - cfg_pars.min_incubation) < community[h].day_epi_start)
           //   people[i].ignore = 1;
           if(people[i].day_ill > community[h].day_epi_stop)
           {
              people[i].infection = 0;
              people[i].symptom = 0;
              people[i].day_ill = MISSING;
           }
           people[i].day_infection_lower = max(people[i].day_ill - cfg_pars.max_incubation, community[h].day_epi_start);
           people[i].day_infection_upper  = max(people[i].day_ill - cfg_pars.min_incubation, community[h].day_epi_start);

           people[i].day_infective_lower = people[i].day_ill + cfg_pars.lower_infectious;
           people[i].day_infective_upper = people[i].day_ill + cfg_pars.upper_infectious;
        }
     
     
        /*file = fopen("check_impute.dat", "w");
        for(i=0; i<p_size; i++)
        if(people[i].size_possible_states > 0)
        {
           fprintf(file, "%d: size=%d  infection=%d  symptom=%d  day_ill=%d\n", 
                   i, people[i].size_possible_states, people[i].infection, people[i].symptom, people[i].day_ill);
           for(j=0; j<people[i].size_possible_states; j++)
              fprintf(file, "%d  ", people[i].possible_states[0][j]);
           fprintf(file, "\n");
           for(j=0; j<people[i].size_possible_states; j++)
              fprintf(file, "%d  ", people[i].possible_states[1][j]);
           fprintf(file, "\n\n");
        }
        fclose(file);*/
     

     
        // Before create risk history, we need to initialize infection status for people with uncertain status.
        //final_risk_day is used to track the last exposure day for people for whom imputation is not neccesary.
        //When imputing infection status of a subject, we need to update the risk history of his contacts;
         //but if we are sure a contact is not susceptible anymore, such updating can be skipped, which
        // will save some computation. 
        for(h=0; h<n_community; h++)
        if(community[h].size > 0)
        {
           for(m=0; m<community[h].size; m++)
           {
              person = people + community[h].member[m];
              if(person->size_possible_states > 0)
              {
                 person->current_state = floor(person->size_possible_states * runiform(&seed));
                 l = person->current_state;
                 if(person->possible_states[0][l] == -INFINITY_INTEGER) //pre-immune
                 {
                    person->pre_immune = 1;
                    person->infection = 0;
                    person->symptom = 0;
                    person->day_ill = MISSING;
                 }
                 else if(person->possible_states[0][l] == INFINITY_INTEGER) //escape
                 {
                    person->pre_immune = 0;
                    person->infection = 0;
                    person->symptom = 0;
                    person->day_ill = MISSING;
                 }
                 else // symptomatic or asymptomatic infection
                 {
                    person->pre_immune = 0;
                    person->infection = 1;
                    person->symptom = person->possible_states[1][l];
                    person->day_ill = person->possible_states[0][l];
                    person->day_infective_lower = person->day_ill + cfg_pars.lower_infectious;
                    person->day_infective_upper = person->day_ill + cfg_pars.upper_infectious;
                    person->day_infection_lower = max(person->day_ill - cfg_pars.max_incubation, community[h].day_epi_start);
                    person->day_infection_upper  = max(person->day_ill - cfg_pars.min_incubation, community[h].day_epi_start);
                 }
              }
              else
              {
                 if(person->pre_immune == 1)  person->final_risk_day = -INFINITY_INTEGER;
                 else
                 {
                    if(person->infection == 1)  person->final_risk_day = person->day_ill - cfg_pars.min_incubation;
                    else  person->final_risk_day = community[h].day_epi_stop;
                 }
              }
           }
        }
     
        // As missing status have been sampled, now set index cases.
        if(cfg_pars.adjust_for_left_truncation == 1 && cfg_pars.preset_index == 0)
        {
           for(h=0; h < n_community; h++)   set_index_cases(community + h);
        }
        // output the imputed population
     
        /*file = fopen("check_pop.dat", "w");
        for(i=0; i<p_size; i++)
        {
           fprintf(file, "%d  %d  preimmune=%d  infection=%d  symptom=%d  day_ill=%d  idx=%d  infection:(lower=%d  upper=%d)  infective:(lower=%d  upper=%d)  ignore=%d\n", 
                i, people[i].community,people[i].pre_immune,  people[i].infection, 
                people[i].symptom, people[i].day_ill, people[i].idx, 
                people[i].day_infection_lower, people[i].day_infection_upper, 
                people[i].day_infective_lower, people[i].day_infective_upper, 
                people[i].ignore);
        }
        fclose(file);*/
     } //if(cfg_pars.EM == 1)
     
     
     // now we are ready to build the risk history based on the contact history.
     // If individuals in the same community share the same contact history, risk history is built at the risk class level;
     // otherwise, risk history is built at the individual level.
     if(cfg_pars.common_contact_history_within_community == 1)
     {
        for(h=0; h<n_community; h++)
        if(community[h].ignore == 0 && community[h].size > 0)
        {
           for(t=community[h].day_epi_start; t<=community[h].day_epi_stop; t++)
           {
              r = t - community[h].day_epi_start;
              ptr_contact = community[h].contact_history[r].c2p_contact;
              while(ptr_contact != NULL)
              {
                 add_c2p_risk_history_to_risk_class(t, h, ptr_contact->contact_mode, ptr_contact->offset);
                 ptr_contact = ptr_contact->next;
              }
           }
        }     
        for(h=0; h<n_community; h++)
        if(community[h].ignore == 0 && community[h].size > 0)
        {
           for(t=community[h].day_epi_start; t<=community[h].day_epi_stop; t++)
           {
              r = t - community[h].day_epi_start;
              ptr_contact = community[h].contact_history[r].p2p_contact;
              while(ptr_contact != NULL)
              {
                  member = people + ptr_contact->contact_id;
                  if(member->infection == 1 && 
                     member->day_infective_lower != MISSING && member->day_infective_upper != MISSING &&
                     t>=member->day_infective_lower && t<=member->day_infective_upper && t <= member->day_exit)
                  {   
                     infective_prob = cfg_pars.prob_infectious[t - member->day_infective_lower];
                     //if(h == 2)  printf("id=%d  t=%d  infective_prob=%f\n", member->id, t, infective_prob);
                     add_p2p_risk_history_to_risk_class(t, member, ptr_contact->contact_mode, ptr_contact->offset, infective_prob, member->symptom);
                  }    
                  ptr_contact = ptr_contact->next;
                 
              }
           }
        } 
                    
        // print out community-level risk history for an arbitrary community to check if we get this part correctly
        if(0)
        {
           h = 2;
           printf("c2p risk history of community %d:\n", h);
           show_c2p_risk_in_risk_class(h, community[h].day_epi_start, community[h].day_epi_stop);
           printf("p2p risk history of community %d:\n",h);
           show_p2p_risk_in_risk_class(h, community[h].day_epi_start, community[h].day_epi_stop);
           goto end;
        }
     }
     else // build individual-level risk history
     {
        for(i=0; i<p_size; i++)
        {
           //printf("person %d\n", i);
           person = people + i;
           if(person->ignore == 0)
           {
              h = person->community;
              //build c2p risk history
              //if(i == 994)  printf("build c2p risk history\n");
              for(t=community[h].day_epi_start; t<=community[h].day_epi_stop; t++)
              if(!(person->final_risk_day != MISSING && t > person->final_risk_day))
              {
                 r = t - community[h].day_epi_start;
                 ptr_contact = person->contact_history[r].c2p_contact;
                 while(ptr_contact != NULL)
                 {
                    add_c2p_risk_history(t, person, ptr_contact->contact_mode, ptr_contact->offset);
                    ptr_contact = ptr_contact->next;
                 }   
              }
              //build p2p risk history
              //if(i == 994) printf("build p2p risk history\n");
              if(person->infection == 1 && 
                 person->day_infective_lower != MISSING && person->day_infective_upper != MISSING)
              {
                 for(t=person->day_infective_lower; t<=person->day_infective_upper; t++)
                 if(t >= community[h].day_epi_start && t <= community[h].day_epi_stop)
                 {
                    r = t - community[h].day_epi_start;
                    infective_prob = cfg_pars.prob_infectious[t - person->day_infective_lower];
                    //if(i == 994)   printf("infective=%d  t=%d  infective_prob=%f\n", person->id, t, infective_prob);
                    inf_ptr_contact = person->contact_history[r].p2p_contact;
                    while(inf_ptr_contact != NULL)
                    {
                       member = people + inf_ptr_contact->contact_id;
                       //if(i == 994)   printf("contact_id=%d  member=%d  final risk day=%d\n", inf_ptr_contact->contact_id, member->id, member->final_risk_day);
                       if(!(member->final_risk_day != MISSING && t > member->final_risk_day) &&  t <= member->day_exit)
                       {
                          sus_ptr_contact = inf_ptr_contact->pair;
                          add_p2p_risk_history(t, member, person, sus_ptr_contact->contact_mode, sus_ptr_contact->offset, infective_prob, person->symptom);
                       }
                       inf_ptr_contact = inf_ptr_contact->next;
                    }
                 }
              }
           }   
        }
        
        // print out individual-level risk history for an arbitrary person to check if we get this part correctly
        if(0)
        {
           for(i=0; i<p_size; i++)
           if(i == 992)
           {
              h = people[i].community;
              printf("c2p risk history of person %d:\n", i);
              show_c2p_risk(people+i, community[h].day_epi_start, community[h].day_epi_stop);
              printf("\n p2p risk history of person %d from %d to %d:\n", i, community[h].day_epi_start, community[h].day_epi_stop);
              show_p2p_risk(people+i, community[h].day_epi_start, community[h].day_epi_stop);
           }
           exit(0); 
        } 
     } 
     //if(cfg_pars.simulation == 0)
     //{
     //   if(cfg_pars.silent_run == 0)  printf("Checking admissibility now\n");
     //   if(admissibility() == 1)  goto end;
     //}
     estimation_error_type = test_error_type = -1;
     /*********************************************************************************
     Estimation starts here.
     estimation_error_type takes 3 values:
     0: no error
     1: The maximum number of function evaluation in Nelder_mead searching is reached
     2: error in variance calculating, which means the final estimates may not be valid
     ***********************************************************************************/
     estimation_error_type = estimation(id_inc, id_inf, id_time, est, &log_likelihood, &var, &var_logit);
     
     if(cfg_pars.write_error_log == 1)
     {      
        sprintf(file_name, "%serror.txt", cfg_pars.path_out);
        file = fopen(file_name, "a");
        fprintf(file, "\n#1.3 Estimation error: \n");
        fprintf(file, "%d\n", estimation_error_type);
        fclose(file);
     }
     if(cfg_pars.silent_run == 0)  printf("estimation error type=%d\n", estimation_error_type);
     
     if(estimation_error_type == 0 && cfg_pars.skip_output == 0)
     {
        //printf("check0, serial number=%d\n", serial_number);
        for(k=0; k<n_b_mode; k++)  
        {  
           b[k] = est[k];
           lb[k] = logit(b[k]);
           se_lb[k] = sqrt(var_logit.data[k][k]);   
           se_b[k] = sqrt(var.data[k][k]);
        }   
        for(k=0; k<n_p_mode; k++)  
        {
           m = n_b_mode + k;           
           p[k] = est[m];
           lp[k] = logit(p[k]);
           se_lp[k] = sqrt(var_logit.data[m][m]);   
           se_p[k] = sqrt(var.data[m][m]);   
        }   
        for(k=0; k<n_u_mode; k++)  
        {
           m = n_b_mode + n_p_mode + k;           
           u[k] = est[m];
           lu[k] = logit(u[k]);
           se_lu[k] = sqrt(var_logit.data[m][m]);   
           se_u[k] = sqrt(var.data[m][m]);   
        }   
        for(k=0; k<n_q_mode; k++)  
        {
           m = n_b_mode + n_p_mode + n_u_mode + k;           
           q[k] = est[m];
           lq[k] = logit(q[k]);
           se_lq[k] = sqrt(var_logit.data[m][m]);   
           se_q[k] = sqrt(var.data[m][m]);   
        }   
        for(k=0; k<n_c2p_covariate; k++)  
        {
           m = n_b_mode + n_p_mode + n_u_mode + n_q_mode + k;           
           OR_c2p[k] = est[m];
           coeff_c2p[k] = log(OR_c2p[k]);
           se_coeff_c2p[k] = sqrt(var_logit.data[m][m]);   
           se_OR_c2p[k] = sqrt(var.data[m][m]);   
        }   
        for(k=0; k<n_p2p_covariate; k++)  
        {
           m = n_b_mode + n_p_mode + n_u_mode + n_q_mode + n_c2p_covariate + k;
           OR_p2p[k] = est[m];
           coeff_p2p[k] = log(OR_p2p[k]);
           se_coeff_p2p[k] = sqrt(var_logit.data[m][m]);   
           se_OR_p2p[k] = sqrt(var.data[m][m]);   
        }
        for(k=0; k<n_pat_covariate; k++)  
        {
           m = n_b_mode + n_p_mode + n_u_mode + n_q_mode + n_c2p_covariate + n_p2p_covariate + k;
           OR_pat[k] = est[m];
           coeff_pat[k] = log(OR_pat[k]);
           se_coeff_pat[k] = sqrt(var_logit.data[m][m]);   
           se_OR_pat[k] = sqrt(var.data[m][m]);   
        }
        for(k=0; k<n_imm_covariate; k++)  
        {
           m = n_b_mode + n_p_mode + n_u_mode + n_q_mode + n_c2p_covariate + n_p2p_covariate + n_pat_covariate + k;
           OR_imm[k] = est[m];
           coeff_imm[k] = log(OR_imm[k]);
           se_coeff_imm[k] = sqrt(var_logit.data[m][m]);   
           se_OR_imm[k] = sqrt(var.data[m][m]);   
        }

        for(k=0; k<n_b_mode; k++)
        {
           lower_b[k] = inv_logit(lb[k] - 1.96 * se_lb[k]);
           upper_b[k] = inv_logit(lb[k] + 1.96 * se_lb[k]);

           CPI[k] = 1.0 - pow((1.0 - b[k]), cfg_pars.CPI_duration);
           der[0] = cfg_pars.CPI_duration * pow(1.0 - b[k], cfg_pars.CPI_duration - 1.0);
           se_CPI[k] = fabs(der[0]) * se_b[k];
           lower_CPI[k] = 1.0 - pow((1.0 - inv_logit(lb[k] - 1.96 * se_lb[k])), cfg_pars.CPI_duration);
           upper_CPI[k] = 1.0 - pow((1.0 - inv_logit(lb[k] + 1.96 * se_lb[k])), cfg_pars.CPI_duration);
        }   
        reset(&der_mat, 0.0);   
        for(k=0; k<n_p_mode; k++)
        {
           lower_p[k] = inv_logit(lp[k] - 1.96 * se_lp[k]);
           upper_p[k] = inv_logit(lp[k] + 1.96 * se_lp[k]);
              
           // Calculate unadjusted SAR
           p_esc = 1.0;
           reset_1d_array_double(der, n_par, 0.0);
           reset_1d_array_double(log_f_lp, n_p_mode, 0.0);
           for(t=cfg_pars.effective_lower_infectious[k]; t<=cfg_pars.effective_upper_infectious[k]; t++)  
           {
              logit_f = lp[k];
              l = t - cfg_pars.lower_infectious;
              s = cfg_pars.prob_infectious[l];
              ff = inv_logit(logit_f) * s;
              f = 1 - ff;
              factor1 = ff / s;
              log_f_lp[k] = - ff * (1 - factor1) / f;
              p_esc *= f;
              //der is d(log p_esc) / d (p[k]) or d(log p_esc) / d coeff_p2p
              der[n_b_mode + k] += log_f_lp[k];
           }
           SAR0[k] = 1.0 - p_esc;
           for(j=0; j<n_par; j++) der_mat.data[k][j] = - p_esc * der[j];  //der1 is d(SAR) / d (lp[k]) or d(SAR) / d coeff_p2p
           var_SAR = 0.0;
           for(j=0; j<n_par; j++)
           {
              for(l=0; l<n_par; l++)
              {
                 var_SAR += der_mat.data[k][j] * der_mat.data[k][l] * var_logit.data[j][l];
              }
           }  
           se_SAR0[k] = sqrt(var_SAR);
           se_logit_SAR = se_SAR0[k] / (SAR0[k] * (1 - SAR0[k]));
           lower_SAR0[k] = inv_logit(logit(SAR0[k]) - 1.96 * se_logit_SAR);
           upper_SAR0[k] = inv_logit(logit(SAR0[k]) + 1.96 * se_logit_SAR);

           // Calculate SAR adjusted for covariates
           for(m=0; m<cfg_pars.SAR_n_covariate_sets; m++)
           {
               p_esc = 1.0;
               reset_1d_array_double(der, n_par, 0.0);
               reset_1d_array_double(log_f_lp, n_p_mode, 0.0);
               if(n_p2p_covariate > 0)  reset_1d_array_double(log_f_p2p, n_p2p_covariate, 0.0);
               for(t=cfg_pars.effective_lower_infectious[k]; t<=cfg_pars.effective_upper_infectious[k]; t++)  
               {
                  if(n_p2p_covariate > 0)
                  {
                     r = t - cfg_pars.SAR_time_dep_lower;
                     if(cfg_pars.SAR_n_covariate_sets > 0)  
                        organize_p2p_covariate_4SAR(m, r, p2p_covariate);
                     else
                        reset_1d_array_double(p2p_covariate, n_p2p_covariate, 0.0);  
                  }       
                  covariate_effect = 0.0;
                  for(j=0; j<n_p2p_covariate; j++)
                      covariate_effect += p2p_covariate[j] * coeff_p2p[j];
                  logit_f = lp[k] + covariate_effect;
                  l = t - cfg_pars.lower_infectious;
                  s = cfg_pars.prob_infectious[l];
                  ff = inv_logit(logit_f) * s;
                  f = 1 - ff;
                  factor1 = ff / s;
                  log_f_lp[k] = - ff * (1 - factor1) / f;
                  for(j=0; j<n_p2p_covariate; j++)
                     log_f_p2p[j] = log_f_lp[k] * p2p_covariate[j];
                  p_esc *= f;
                  //der is d(log p_esc) / d (p[k]) or d(log p_esc) / d coeff_p2p
                  der[n_b_mode + k] += log_f_lp[k];
                  for(j=0; j<n_p2p_covariate; j++)
                     der[n_b_mode + n_p_mode + n_u_mode + n_q_mode + n_c2p_covariate + j] += log_f_p2p[j];
                  //printf("t=%d  s=%e\n", t, s);
                  //printf("   p2p covariates:");
                  //for(j=0; j<n_p2p_covariate; j++)  printf("%e  ", p2p_covariate[j]);
                  //printf("\n");
                  //printf("   covariate_effect=%e  f=%e  p_esc=%e\n", covariate_effect, f, p_esc);
               }
               SAR[m][k] = 1.0 - p_esc;
               for(j=0; j<n_par; j++) der_mat.data[k][j] = - p_esc * der[j];  //der1 is d(SAR) / d (lp[k]) or d(SAR) / d coeff_p2p
               var_SAR = 0.0;
               for(j=0; j<n_par; j++)
               {
                  for(l=0; l<n_par; l++)
                  {
                     var_SAR += der_mat.data[k][j] * der_mat.data[k][l] * var_logit.data[j][l];
                  }
               }  
               se_SAR[m][k] = sqrt(var_SAR);
               se_logit_SAR = se_SAR[m][k] / (SAR[m][k] * (1 - SAR[m][k]));
               lower_SAR[m][k] = inv_logit(logit(SAR[m][k]) - 1.96 * se_logit_SAR);
               upper_SAR[m][k] = inv_logit(logit(SAR[m][k]) + 1.96 * se_logit_SAR);
           }
        } 
        //printf("Derivative matrix:\n");  
        //fmprintf(&der_mat);


        if(cfg_pars.R0_multiplier_provided > 0 && n_p_mode > 0)   
        {
           //Calculate unadjusted R0
           R0 = 0.0;
           reset_1d_array_double(der1, n_par, 0.0);
           reset_1d_array_double(der2, n_par, 0.0);
           for(k=0; k<n_p_mode; k++)
           {
              R0 += cfg_pars.R0_multiplier[k] * SAR0[k];
              for(l=0; l<n_par; l++)  der1[l] += cfg_pars.R0_multiplier[k] * der_mat.data[k][l];  //der[k] is d(R0) / d (lp[k])
              der2[k] += SAR0[k];
           }
           var_R0 = 0.0;
           for(k=0; k<n_par; k++)
           {
              for(l=0; l<n_par; l++)
              {
                 var_R0 += der1[k] * der1[l] * var_logit.data[k][l];
              }
           }  
           for(k=0; k<n_p_mode; k++)  var_R0 += der2[k] * der2[k] * cfg_pars.R0_multiplier_var[k];

           se_R0 = sqrt(var_R0);
           lower_R0 = exp(log(R0) - 1.96 * se_R0 / R0);
           upper_R0 = exp(log(R0) + 1.96 * se_R0 / R0);
           
           // Calculate R adjusted for covariates. Sometimes R0 also needs adjustment for covariates, e.g., for SARS-CoV-2,
           // if you assume infectivity differs before and after symptom onset and estimate the difference via 
           // a time-dependent covariate (basically an indicator for before vs. after symptom onset),
           // then you need to adjust for that covariate to cacculate a meaningful R0.
           for(m=0; m<cfg_pars.SAR_n_covariate_sets; m++)
           {
               R0_adj[m] = 0.0;
               reset_1d_array_double(der1, n_par, 0.0);
               reset_1d_array_double(der2, n_par, 0.0);
               for(k=0; k<n_p_mode; k++)
               {
                  R0_adj[m] += cfg_pars.R0_multiplier[k] * SAR[m][k];
                  for(l=0; l<n_par; l++)  der1[l] += cfg_pars.R0_multiplier[k] * der_mat.data[k][l];  //der[k] is d(R0) / d (lp[k])
                  der2[k] += SAR[m][k];
               }
               var_R0 = 0.0;
               for(k=0; k<n_par; k++)
               {
                  for(l=0; l<n_par; l++)
                  {
                     var_R0 += der1[k] * der1[l] * var_logit.data[k][l];
                  }
               }  
               for(k=0; k<n_p_mode; k++)  var_R0 += der2[k] * der2[k] * cfg_pars.R0_multiplier_var[k];

               se_R0_adj[m] = sqrt(var_R0);
               lower_R0_adj[m] = exp( log(R0_adj[m]) - 1.96 * se_R0_adj[m] / R0_adj[m]);
               upper_R0_adj[m] = exp( log(R0_adj[m]) + 1.96 * se_R0_adj[m] / R0_adj[m]);
           }   
        }
           
        for(k=0; k<n_u_mode; k++)
        {
           lower_u[k] = inv_logit(lu[k] - 1.96 * se_lu[k]);
           upper_u[k] = inv_logit(lu[k] + 1.96 * se_lu[k]);
        }   
              
        for(k=0; k<n_q_mode; k++)
        {
           lower_q[k] = inv_logit(lq[k] - 1.96 * se_lq[k]);
           upper_q[k] = inv_logit(lq[k] + 1.96 * se_lq[k]);
        }   
              
        for(k=0; k<n_c2p_covariate; k++)
        {
           lower_OR_c2p[k] = exp(coeff_c2p[k] - 1.96 * se_coeff_c2p[k]);
           upper_OR_c2p[k] = exp(coeff_c2p[k] + 1.96 * se_coeff_c2p[k]);
        }   
           
        for(k=0; k<n_p2p_covariate; k++)
        {
           lower_OR_p2p[k] = exp(coeff_p2p[k] - 1.96 * se_coeff_p2p[k]);
           upper_OR_p2p[k] = exp(coeff_p2p[k] + 1.96 * se_coeff_p2p[k]);
        }   

        for(k=0; k<n_pat_covariate; k++)
        {
           lower_OR_pat[k] = exp(coeff_pat[k] - 1.96 * se_coeff_pat[k]);
           upper_OR_pat[k] = exp(coeff_pat[k] + 1.96 * se_coeff_pat[k]);
        }   
        for(k=0; k<n_imm_covariate; k++)
        {
           lower_OR_imm[k] = exp(coeff_imm[k] - 1.96 * se_coeff_imm[k]);
           upper_OR_imm[k] = exp(coeff_imm[k] + 1.96 * se_coeff_imm[k]);
        }   
        //printf("b=%e  p=%e  OR=%e\n", b[0], p[0], OR_p2p[0]);  

        sprintf(file_name, "%soutput.txt", cfg_pars.path_out);
        if((file = fopen(file_name, "a")) == NULL)
            file = fopen(file_name, "w");
        if(cfg_pars.simplify_output == 1)
        {
           fprintf(file, "%6d  %3d  %3d  ", serial_number, id_inc, id_inf);   
           if(cfg_pars.R0_divide_by_time == 1)  fprintf(file, "%6d  ", id_time);
           for(k=0; k<n_b_mode; k++)
           {
              if(cfg_pars.skip_variance == 0)
                 fprintf(file, "%e  %e  %e  %e  ", b[k], se_b[k], lower_b[k], upper_b[k]);
              else
                 fprintf(file, "%e  ", b[k]);
           }              
           for(k=0; k<n_p_mode; k++)
           {
              if(cfg_pars.skip_variance == 0)
                 fprintf(file, "%e  %e  %e  %e  ", p[k], se_p[k], lower_p[k], upper_p[k]);
              else
                 fprintf(file, "%e  ", p[k]);
           }              
           for(k=0; k<n_u_mode; k++)
           {
              if(cfg_pars.skip_variance == 0)
                 fprintf(file, "%e  %e  %e  %e  ", u[k], se_u[k], lower_u[k], upper_u[k]);
              else
                 fprintf(file, "%e  ", u[k]);
           }   
           for(k=0; k<n_q_mode; k++)
           {
              if(cfg_pars.skip_variance == 0)
                 fprintf(file, "%e  %e  %e  %e  ", q[k], se_q[k], lower_q[k], upper_q[k]);
              else
                 fprintf(file, "%e  ", q[k]);
           }   
              
           for(k=0; k<n_c2p_covariate; k++)
           {
              if(cfg_pars.skip_variance == 0)
                 fprintf(file, "%7.3f  %7.3f  %7.3f  %7.3f  ", OR_c2p[k], se_OR_c2p[k], lower_OR_c2p[k], upper_OR_c2p[k]);
              else
                 fprintf(file, "%e  ", OR_c2p[k]);
           }   
           
           for(k=0; k<n_p2p_covariate; k++)
           {
              if(cfg_pars.skip_variance == 0)
                 fprintf(file, "%7.3f  %7.3f  %7.3f  %7.3f  ", OR_p2p[k], se_OR_p2p[k], lower_OR_p2p[k], upper_OR_p2p[k]);
              else
                 fprintf(file, "%e  ", OR_p2p[k]);
           }   

           for(k=0; k<n_pat_covariate; k++)
           {
              if(cfg_pars.skip_variance == 0)
                 fprintf(file, "%7.3f  %7.3f  %7.3f  %7.3f", OR_pat[k], se_OR_pat[k], lower_OR_pat[k], upper_OR_pat[k]);
              else
                 fprintf(file, "%e  ", OR_pat[k]);
           }   
              
           for(k=0; k<n_imm_covariate; k++)
           {
              if(cfg_pars.skip_variance == 0)
                 fprintf(file, "%7.3f  %7.3f  %7.3f  %7.3f", OR_imm[k], se_OR_imm[k], lower_OR_imm[k], upper_OR_imm[k]);
              else
                 fprintf(file, "%e  ", OR_imm[k]);
           }   
             
           if(cfg_pars.skip_variance == 0)   
           for(m=0; m<n_par; m++)
           {
              for(n=m; n<n_par; n++)
              {
                 fprintf(file, "%e  ", var_logit.data[m][n]);
              }
           }
           fprintf(file, "\n");
        }
        else
        {
           if(cfg_pars.stat_test == 1)
              fprintf(file, "\n# P-value for testing p=%f\n", p_value);
     
           if(n_b_mode > 0)
           {
              fprintf(file, "\n# Estimates of b\n");
              for(k=0; k<n_b_mode; k++)
                 fprintf(file, "%e,  %e,  %e,  %e\n", b[k], se_b[k], lower_b[k], upper_b[k]);
           }
           if(n_p_mode > 0)
           {
              fprintf(file, "\n# Estimates of p\n");
              for(k=0; k<n_p_mode; k++)
                 fprintf(file, "%e,  %e,  %e,  %e\n", p[k], se_p[k], lower_p[k], upper_p[k]);
           }
           if(n_u_mode > 0)
           {
              fprintf(file, "\n# Estimates of u\n");
              for(k=0; k<n_u_mode; k++)
                 fprintf(file, "%e,  %e,  %e,  %e\n", u[k], se_u[k], lower_u[k], upper_u[k]);
           }
           if(n_q_mode > 0)
           {
              fprintf(file, "\n# Estimates of q\n");
              for(k=0; k<n_q_mode; k++)
                 fprintf(file, "%e,  %e,  %e,  %e\n", q[k], se_q[k], lower_q[k], upper_q[k]);
           }
           if(n_c2p_covariate > 0)
           {    
              fprintf(file, "\n# Estimates of odds ratios for c2p exposure\n");
              for(k=0; k<n_c2p_covariate; k++)
              {
                 p_value_OR = 2.0 * (1.0 - pnorm(fabs(coeff_c2p[k] / se_coeff_c2p[k]), 0, 1));
                 fprintf(file, "%20.8f,  %20.8f,  %20.8f,  %20.8f,  %20.8f\n", OR_c2p[k], se_OR_c2p[k], lower_OR_c2p[k], upper_OR_c2p[k], p_value_OR);
              }   
           }
           if(n_p2p_covariate > 0)
           {    
              fprintf(file, "\n# Estimates of odds ratios for p2p exposure\n");
              for(k=0; k<n_p2p_covariate; k++)
              {
                 p_value_OR = 2.0 * (1.0 - pnorm(fabs(coeff_p2p[k] / se_coeff_p2p[k]), 0, 1));
                 fprintf(file, "%20.8f,  %20.8f,  %20.8f,  %20.8f,  %20.8f\n", OR_p2p[k], se_OR_p2p[k], lower_OR_p2p[k], upper_OR_p2p[k], p_value_OR);
              }   
           }
           if(n_pat_covariate > 0)
           {    
              fprintf(file, "\n# Estimates of odds ratios for pathogenicity\n");
              for(k=0; k<n_pat_covariate; k++)
              {
                 p_value_OR = 2.0 * (1.0 - pnorm(fabs(coeff_pat[k] / se_coeff_pat[k]), 0, 1));
                 fprintf(file, "%20.8f,  %20.8f,  %20.8f,  %20.8f,  %20.8f\n", 
                         OR_pat[k], se_OR_pat[k], lower_OR_pat[k], upper_OR_pat[k], p_value_OR);
              }   
           }

           if(n_imm_covariate > 0)
           {    
              fprintf(file, "\n# Estimates of odds ratios for pre-immunity level\n");
              for(k=0; k<n_imm_covariate; k++)
              {
                 p_value_OR = 2.0 * (1.0 - pnorm(fabs(coeff_imm[k] / se_coeff_imm[k]), 0, 1));
                 fprintf(file, "%20.8f,  %20.8f,  %20.8f,  %20.8f,  %20.8f\n", 
                         OR_imm[k], se_OR_imm[k], lower_OR_imm[k], upper_OR_imm[k], p_value_OR);
              }   
           }
           if(n_b_mode > 0 && cfg_pars.CPI_duration > 0)
           {
              fprintf(file, "\n# Estimates of CPI\n");
              for(k=0; k<n_b_mode; k++)
              {
                 fprintf(file, "%e,  %e,  %e,  %e\n", CPI[k], se_CPI[k], lower_CPI[k], upper_CPI[k]);
              }   
           }
           if(n_p_mode > 0)
           {
              fprintf(file, "\n# Estimates of unadjusted SAR\n");
              for(k=0; k<n_p_mode; k++)
              {
                 fprintf(file, "%e,  %e,  %e,  %e\n", SAR0[k], se_SAR0[k], lower_SAR0[k], upper_SAR0[k]);
              }
              if(cfg_pars.SAR_n_covariate_sets > 0)
              {
                 fprintf(file, "\n# Estimates of SAR adjusted for covariates\n");
                 for(m=0; m < cfg_pars.SAR_n_covariate_sets; m++)
                 {
                    fprintf(file, "\n# Covariate set %d\n", m);
                    for(k=0; k<n_p_mode; k++)
                    {
                       fprintf(file, "%e,  %e,  %e,  %e\n", SAR[m][k], se_SAR[m][k], lower_SAR[m][k], upper_SAR[m][k]);
                    }   
                 }
              }
              if(cfg_pars.R0_multiplier_provided > 0)
              {
                 fprintf(file, "\n# Estimates of unadjusted R0\n");
                 fprintf(file, "%20.8f,  %20.8f,  %20.8f,  %20.8f\n", R0, se_R0, lower_R0, upper_R0);
                 if(cfg_pars.SAR_n_covariate_sets > 0)
                 {
                    fprintf(file, "\n# Estimates of R0 adjusted for covariates\n");
                    for(m=0; m<cfg_pars.SAR_n_covariate_sets; m++)
                    {
                       fprintf(file, "\n# Covariate set %d\n", m);
                       fprintf(file, "%20.8f,  %20.8f,  %20.8f,  %20.8f\n", R0_adj[m], se_R0_adj[m], lower_R0_adj[m], upper_R0_adj[m]);
                    }
                 }
              }   
           }

           fprintf(file, "\n#Log-likelihood\n");
           fprintf(file, "%20.8f\n", log_likelihood);
           
           if(cfg_pars.skip_variance == 0)   
           {
              fprintf(file, "\n#Covariance Matrix Estimates\n");
              fmfprintf(file, &var, 1);
     
              fprintf(file, "\n#Covariance Matrix Estimates for logit(b), logit(p) and log(OR)\n");
              fmfprintf(file, &var_logit, 1);
           }   
        }
        fclose(file);

        // output unadjusted SAR and R0 to separate files.
        if(cfg_pars.simplify_output_SAR == 1 && n_p_mode > 0)
        {
           sprintf(file_name, "%soutput_SAR.txt", cfg_pars.path_out);
           if((file = fopen(file_name, "a")) == NULL)
              file = fopen(file_name, "w");
           fprintf(file, "%6d  %3d  %3d  ", serial_number, id_inc, id_inf); 
           if(cfg_pars.R0_divide_by_time == 1)  fprintf(file, "%d  ", id_time);
           for(k=0; k<n_p_mode; k++)
              fprintf(file, "  %e  %e  %e  %e", SAR0[k], se_SAR0[k], lower_SAR0[k], upper_SAR0[k]);
           fprintf(file, "\n");
           fclose(file);   
        }
        if(cfg_pars.simplify_output_R0 == 1 && n_p_mode > 0)
        {
           sprintf(file_name, "%soutput_R0.txt", cfg_pars.path_out);
           if((file = fopen(file_name, "a")) == NULL)
              file = fopen(file_name, "w");
           fprintf(file, "%6d  %3d  %3d  %6d", serial_number, id_inc, id_inf, id_time);   
           fprintf(file, "  %e  %e  %e  %e\n", R0, se_R0, lower_R0, upper_R0);
           fclose(file);   
        }
     } //if(estimation_error_type == 0 && skip_output == 0)
     //clean structures that are used for each specific epidemic.
     //These structures are created during the estimation
loop_end:  
     restore_susceptibility();
     free(importance_weight);
     size_sample_states = 0;
  }
end:
  free_arrays();
  deflate_matrix(&var);
  deflate_matrix(&var_logit);
  deflate_matrix(&der_mat);
  free(est);
  free(der); free(der1); free(der2);
  free(p2p_covariate);
  free(value);
  free(pdf_incubation);
  free(sdf_incubation);



  if(people != NULL)
  {
     for(i=0; i<p_size; i++)
     {
        h = people[i].community;
        free(people[i].time_ind_covariate);
        free_2d_array_double(people[i].time_dep_covariate);

        if(people[i].pat_covariate != NULL) free_2d_array_double(people[i].pat_covariate);
        if(people[i].imm_covariate != NULL) free(people[i].imm_covariate);
        if(people[i].possible_states != NULL) free_2d_array_int(people[i].possible_states);

        if(people[i].contact_history != NULL)
        {
           for(r=0; r<community[h].epi_duration; r++)
           {
              ptr_contact = people[i].contact_history[r].c2p_contact;
              while(ptr_contact != NULL)
              {
                 ptr2_contact = ptr_contact->next;
                 //if(ptr_contact->covariate != NULL)  free(ptr_contact->covariate);
                 free(ptr_contact);
                 ptr_contact = ptr2_contact;
              }
              ptr_contact = people[i].contact_history[r].p2p_contact;
              while(ptr_contact != NULL)
              {
                 ptr2_contact = ptr_contact->next;
                 //if(ptr_contact->covariate != NULL)  free(ptr_contact->covariate);
                 free(ptr_contact);
                 ptr_contact = ptr2_contact;
              }
              ptr_risk = people[i].risk_history[r].c2p_risk;
              while(ptr_risk != NULL)
              {
                 ptr2_risk = ptr_risk->next;
                 if(ptr_risk->covariate != NULL)  free(ptr_risk->covariate);
                 free(ptr_risk);
                 ptr_risk = ptr2_risk;
              }
              ptr_risk = people[i].risk_history[r].p2p_risk;
              while(ptr_risk != NULL)
              {
                 ptr2_risk = ptr_risk->next;
                 if(ptr_risk->covariate != NULL)  free(ptr_risk->covariate);
                 free(ptr_risk);
                 ptr_risk = ptr2_risk;
              }
           }
           free(people[i].contact_history);
           free(people[i].risk_history);
        }   
     }
     free(people);
  }   
  
  if(community != NULL)
  {
     for(h=0; h < n_community; h++)
     {
        if(community[h].member != NULL)  free(community[h].member);
        if(community[h].idx != NULL)  free(community[h].idx);
        if(community[h].member_impute != NULL)  free(community[h].member_impute);
        if(community[h].sample_states != NULL)  free_state_chain(community[h].sample_states, community[h].sample_states_rear);
        if(community[h].list_states != NULL)  free_state_chain(community[h].list_states, community[h].list_states_rear);
        ptr_class = community[h].risk_class;
        while(ptr_class != NULL) 
        {
           ptr_integer = ptr_class->member;
           while(ptr_integer != NULL)
           {
              ptr2_integer = ptr_integer->next;
              free(ptr_integer);
              ptr_integer = ptr2_integer;
           }
           if(ptr_class->risk_history != NULL)  free(ptr_class->risk_history); // c2p and p2p risk history at each node of risk_history array has been released in restore_susceptibility()
           ptr2_class = ptr_class->next;
           free(ptr_class);
           ptr_class = ptr2_class;
        }    
        if(community[h].contact_history != NULL)
        {
           for(r=0; r<community[h].epi_duration; r++)
           {
              ptr_contact = community[h].contact_history[r].c2p_contact;
              while(ptr_contact != NULL)
              {
                 ptr2_contact = ptr_contact->next;
                 //if(ptr_contact->covariate != NULL)  free(ptr_contact->covariate);
                 free(ptr_contact);
                 ptr_contact = ptr2_contact;
              }
              ptr_contact = community[h].contact_history[r].p2p_contact;
              while(ptr_contact != NULL)
              {
                 ptr2_contact = ptr_contact->next;
                 //if(ptr_contact->covariate != NULL)  free(ptr_contact->covariate);
                 free(ptr_contact);
                 ptr_contact = ptr2_contact;
              }
           }
           free(community[h].contact_history);
        }   
     }
     free(community);
  }

  t2 = time(NULL);
  if(cfg_pars.silent_run == 0)  printf("\n use time: %f seconds\n", difftime(t2, t1));
  return(0);
}


