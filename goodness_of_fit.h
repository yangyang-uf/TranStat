
// Calculate observed and expected ILI frequencies as well as the deriavatives of the expected frequencies for each given community
int goodness_of_fit_community(int id_community, double *est, int *at_risk, double *observed_INF_freq, double *expected_INF_freq, MATRIX *first)
{
  int i, j, h, k, l, m, n, r, t, rr, rrr, tt, length, len;
  int start_day, stop_day, max_stop_day, likelihood_start_day;
  int inf1, inf2, inf3, inf4, shift;
  int found, b_mode, p_mode, q_mode, u_mode;
  int n_par, n_b_mode, n_p_mode, n_q_mode, n_u_mode;
  int n_c2p_covariate, n_p2p_covariate, n_imm_covariate, n_pat_covariate;
  int n_covariate, n_time_ind_covariate, n_time_dep_covariate;
  int n_sus_p2p_covariate, n_inf_p2p_covariate, n_int_p2p_covariate; 
  int n_par_equiclass;
  
  double factor1, factor2, factor3;
  double e, f, ff, s, p_inf, U, L, temp, mean_inc;
  double covariate_effect, cum_e, day_L, sum_obs_freq, sum_exp_freq;
  
  double sdf, pdf, pr, *der;
  double logit_f;
  double asym_effect;
  double prob_escape, prob_ILI, prob_noILI;
  double **score_lb, **score_lp, **score_lu, **score_lq; 
  double **score_c2p, **score_p2p, **score_pat, **score_imm;
  double *weight, sum_weight;
  double max_likelihood, max_weight;

  RISK *ptr_risk;
  PEOPLE *person, *member;
  RISK_CLASS *ptr_class;
  INTEGER_CHAIN *ptr_integer;
  STATE *ptr_stat, *ptr_stat_max_likelihood, *ptr_stat_max_weight;


  /***********************************************************************************************************
   * n_b_mode: number of types/modes of community-to-person contact. 
   * n_p_mode: number of types/modes of person-to-person contact.
   * n_time_ind_covariate: number of time-independent covariates
   * n_time_dep_covariate: number of time-dependent covariates
   * n_covariate: total number of covariates
   * n_c2p_covariate: number of covariates of the susceptible person that can modify community-to-person risk
   * n_sus_p2p_covariate: number of covariates of the susceptible person that can modify person-to-person risk
   * n_inf_p2p_covariate: number of covariates of the infective person that can modify person-to-person risk
   * n_inf_p2p_covariate: number of interactions between covariates of the susceptible person
   *                      and the infective person that can modify person-to-person risk
   * n_p2p_covariate: total number of covariates that can modify person-to-person risk, which is the sum
   *                  of above three.
   * The user need to supply  n_time_ind_covariate, n_time_dep_covariate, n_c2p_covariate,
   * n_sus_p2p_covariate, n_inf_p2p_covariate and n_int_p2p_covariate.
   ***********************************************************************************************************/
  n_b_mode = cfg_pars.n_b_mode;
  n_p_mode = cfg_pars.n_p_mode;
  n_u_mode = cfg_pars.n_u_mode;
  n_q_mode = cfg_pars.n_q_mode;
  n_time_ind_covariate = cfg_pars.n_time_ind_covariate;
  n_time_dep_covariate = cfg_pars.n_time_dep_covariate;
  n_covariate = cfg_pars.n_covariate;
  n_c2p_covariate = cfg_pars.n_c2p_covariate;
  n_sus_p2p_covariate = cfg_pars.n_sus_p2p_covariate;
  n_inf_p2p_covariate = cfg_pars.n_inf_p2p_covariate;
  n_int_p2p_covariate = cfg_pars.n_int_p2p_covariate;
  n_p2p_covariate = cfg_pars.n_p2p_covariate;
  n_pat_covariate = cfg_pars.n_pat_covariate;
  n_imm_covariate = cfg_pars.n_imm_covariate;
  n_par = cfg_pars.n_par;
  n_par_equiclass = cfg_pars.n_par_equiclass;
  asym_effect = cfg_pars.asym_effect_est;
 

  /************************************************************************************************************
   * f is person-to-person escape probability. e is the daily escape probability for a susceptible.
   * log_f_lb is the 1st derivative of log(f) with respect to lb.
   * log_f_lb_lb is the second derivative of log(f) with respect to lb.
   * Other derivative terms are similarly defined.
   *
   * e is the daily escape probability for a susceptible from all contacts.
   * cum_e and cum_log_e_* are cumulatives of e and log_e_* over a short period, 
   * from day_ill-max_latent to t, where day_ill-max_latent<=t<=day_ill-min_latent.
   * 
   * ee[t] is the escape probability for day t. 
   * cum_log_ee[t] = log(ee[1] * ee[2] * ... * ee[t])
   * 
   * temp_* variables are just auxiliary variables to help get cum_log_ee_* variables.
   * 
   * day_L is the likelihood for each infected person for the period 
   * from day_ill-max_latent+1 to t, where day_ill-max_latent<=t<=day_ill-min_latent.
   * log_day_L_* are derivatives of log(day_L). 
   *
   * L is the overall likelihood, and L_* and log_L_* are derivatives of L and log(L) respectively.
   * *********************************************************************************************************/
  
  for(k=0; k<n_b_mode; k++)  
  {  
     b[k] = est[k];
     lb[k] = logit(b[k]);
  }   
  m = n_b_mode;           
  for(k=0; k<n_p_mode; k++)  
  {
     p[k] = est[m+k];
     lp[k] = logit(p[k]);
  }   
  m = n_b_mode + n_p_mode + n_u_mode + n_q_mode;           
  for(k=0; k<n_c2p_covariate; k++)  
  {
     OR_c2p[k] = est[m+k];
     coeff_c2p[k] = log(OR_c2p[k]);
  }   
  m = n_b_mode + n_p_mode + n_u_mode + n_q_mode + n_c2p_covariate;           
  for(k=0; k<n_p2p_covariate; k++)  
  {
     OR_p2p[k] = est[m+k];
     coeff_p2p[k] = log(OR_p2p[k]);
  }

  h = id_community;
  
  if(community[h].size > 0 && community[h].ignore == 0)
  {
     length = community[h].day_epi_stop - community[h].day_epi_start + 1;
     if(n_b_mode > 0)  make_2d_array_double(&score_lb, length, n_b_mode, 0.0); 
     if(n_b_mode > 0)  make_2d_array_double(&score_lp, length, n_p_mode, 0.0); 
     if(n_c2p_covariate > 0)  make_2d_array_double(&score_c2p, length, n_c2p_covariate, 0.0); 
     if(n_p2p_covariate > 0)  make_2d_array_double(&score_p2p, length, n_p2p_covariate, 0.0); 
     
     // index cases are now included in the risk classes.
     ptr_class = community[h].risk_class;
     while(ptr_class != NULL)
     {
        if(cfg_pars.common_contact_history_within_community == 1)
        {
           // find the person with the largest stop_day
           max_stop_day = -INFINITY_INTEGER;
           ptr_integer = ptr_class->member;
           while(ptr_integer != NULL)
           {
              member = people + ptr_integer->id;
              if(member->idx == 0 && member->pre_immune == 0)
              {
                 stop_day = (member->infection == 1)? member->day_infection_upper:min(community[h].day_epi_stop, member->day_exit);
                 if(stop_day > max_stop_day)
                 {
                    max_stop_day = stop_day;
                    person = member;
                 }
              }
              //if(community[i].id == 6)  printf("%d  infection=%d idx=%d stop_day=%d max_stop_day=%d\n", member->id, member->infection, member->idx, stop_day, max_stop_day); 
              ptr_integer = ptr_integer->next;
           }
           
           //printf("person %d is chosen\n", person->id);
           if(cfg_pars.adjust_for_left_truncation == 1 && community[h].earliest_idx_day_ill != MISSING)
              start_day = max(community[h].day_epi_start, community[h].earliest_idx_day_ill - cfg_pars.max_incubation + 1);
           else  start_day = community[h].day_epi_start;
           stop_day = max_stop_day;
           //printf("start_day=%d  stop_day=%d  day_epi_stop=%d\n", start_day, stop_day, community[h].day_epi_stop);
           
           for(t=start_day; t<=stop_day; t++)
           {
              r = t - start_day;
              rr = t - community[h].day_epi_start; //the time reference for risk history is community->day_epi_start, not start_day.
              e = 1.0;
              for(j=0; j<n_b_mode; j++) log_e_lb[j] = 0.0;
              for(j=0; j<n_p_mode; j++) log_e_lp[j] = 0.0;
              for(k=0; k<n_c2p_covariate; k++) log_e_c2p[k] = 0.0;
              for(k=0; k<n_p2p_covariate; k++) log_e_p2p[k] = 0.0;
              //printf("person %d: start_day=%d  t=%d  r=%d  rr=%d\n", person->id, start_day, t, r, rr);
              //calculate log_f_* variables related to the probability of escaping risk from common source.
              //update log_e_* variables accordingly
              
              ptr_risk = ptr_class->risk_history[rr].c2p_risk;
              while(ptr_risk != NULL)
              {
                 covariate_effect = 0.0;
                 for(k=0; k<n_c2p_covariate; k++)
                    covariate_effect += ptr_risk->covariate[k] * coeff_c2p[k];
                 covariate_effect += ptr_risk->offset;
                 
                 b_mode = ptr_risk->contact_mode;
                 logit_f = lb[b_mode] + covariate_effect;
                 
                 ff = inv_logit(logit_f);
                 f = 1.0 - ff;
                 log_f_lb[b_mode] = - ff;
                 
                 e *= ipow(f, ptr_risk->size);
                 //printf("contact_mode=%d  lb=%e  covariate effect=%e  logit_f=%e  ff=%e  f=%e  e=%e\n", b_mode, lb[b_mode], covariate_effect, logit_f, ff, f, e);
                 log_e_lb[b_mode] += ptr_risk->size * log_f_lb[b_mode];
                 for(k=0; k<n_c2p_covariate; k++)
                    log_e_c2p[k] += ptr_risk->size * log_f_c2p[k];
                 ptr_risk = ptr_risk->next;
              }
              
              //calculate log_f_* variables related to the probability of escaping risk from infective people.
              //update log_e_* variables accordingly
              //The difference in infectiousness level between symptomatic and asymptomatic cases
              //is adjusted by one additional covariate that denote the symptom status of the infective person.
              ptr_risk = ptr_class->risk_history[rr].p2p_risk;
              while(ptr_risk != NULL)
              {
                 covariate_effect = 0.0;
                 for(k=0; k<n_p2p_covariate; k++)
                    covariate_effect += ptr_risk->covariate[k] * coeff_p2p[k];
                 covariate_effect += ptr_risk->offset;
                 
                 p_mode = ptr_risk->contact_mode;
                 logit_f = lp[p_mode] + covariate_effect;
                 s = ptr_risk->infective_prob;
                 //printf("id=%d  t=%d  s=%e  size=%d  p0=%e  p1=%e  ", 
                 //        person->id, t, s, ptr_risk->size, inv_logit(lp[p_mode]), inv_logit(logit_f));
                 //for(k=0; k<n_p2p_covariate; k++)  printf("x=%e  beta=%e  ", ptr_risk->covariate[k], coeff_p2p[k]);
                 //printf("\n");
                 if( s > 0)
                 {
                    ff = inv_logit(logit_f) * s * bipow(asym_effect, 1 - ptr_risk->symptom);
                    f = 1 - ff;
                    factor1 = ff / s;
                    log_f_lp[p_mode] = - ff * (1 - factor1) / f;
                    
                    for(k=0; k<n_p2p_covariate; k++)
                       log_f_p2p[k] = log_f_lp[p_mode] * ptr_risk->covariate[k];
                    
                    if(f <= 0.0 || f > 1.0)
                    {
                       printf("Goodness-of-fit: f=%e\n", f);
                       free(observed_INF_freq);
                       free(expected_INF_freq);
                       return(1);
                    }
                    e *= ipow(f, ptr_risk->size);
                    log_e_lp[p_mode] += ptr_risk->size * log_f_lp[p_mode];
                    for(k=0; k<n_p2p_covariate; k++)
                       log_e_p2p[k] += ptr_risk->size * log_f_p2p[k];
                 }
                 ptr_risk = ptr_risk->next;
              }
              //if(e > 0 && e < 1) printf("t=%d  e=%e\n", t, e);
              ee[r] = e;
              for(j=0; j<n_b_mode; j++) log_ee_lb[r][j] = log_e_lb[j];
              for(j=0; j<n_p_mode; j++) log_ee_lp[r][j] = log_e_lp[j];
              for(k=0; k<n_c2p_covariate; k++) log_ee_c2p[r][k] = log_e_c2p[k];
              for(k=0; k<n_p2p_covariate; k++) log_ee_p2p[r][k] = log_e_p2p[k];
           } /* end of day t */
        }
        
        ptr_integer = ptr_class->member;
        while(ptr_integer != NULL)
        {
           person = people + ptr_integer->id;
           if(person->idx == 1)  
           {
              ptr_integer = ptr_integer->next;
              continue;
           }   
           if(n_u_mode > 0)  u_mode = person->u_mode;

           // if people in the comunity do not share the same exposure/risk history,
           // calculate individual-level exposure/risk history.
           if(cfg_pars.common_contact_history_within_community == 0 && person->pre_immune == 0)
           {
              if(cfg_pars.adjust_for_left_truncation == 1 && community[h].earliest_idx_day_ill != MISSING)
                 start_day = max(community[h].day_epi_start, community[h].earliest_idx_day_ill - cfg_pars.max_incubation + 1);
              else  start_day = community[h].day_epi_start;
              stop_day = (person->infection == 1)? person->day_infection_upper:min(community[h].day_epi_stop, person->day_exit);
              for(t=start_day; t<=stop_day; t++)
              {
                 r = t - start_day;
                 rr = t - community[h].day_epi_start; //the time reference for risk history is community[h].day_epi_start, not start_day.

                 e = 1.0;
                 for(j=0; j<n_b_mode; j++) log_e_lb[j] = 0.0;
                 for(j=0; j<n_p_mode; j++) log_e_lp[j] = 0.0;
                 for(k=0; k<n_c2p_covariate; k++) log_e_c2p[k] = 0.0;
                 for(k=0; k<n_p2p_covariate; k++) log_e_p2p[k] = 0.0;
                 //calculate log_f_* variables related to the probability of escaping risk from community per day.
                 //update log_e_* variables accordingly
                 ptr_risk = person->risk_history[rr].c2p_risk;
                 while(ptr_risk != NULL)
                 {
                    covariate_effect = 0.0;
                    for(k=0; k<n_c2p_covariate; k++)
                       covariate_effect += ptr_risk->covariate[k] * coeff_c2p[k];
                    covariate_effect += ptr_risk->offset;
                    
                    b_mode = ptr_risk->contact_mode;
                    logit_f = lb[b_mode] + covariate_effect;

                    ff = inv_logit(logit_f);
                    f = 1.0 - ff;
                    log_f_lb[b_mode] = - ff;
                    
                    for(k=0; k<n_c2p_covariate; k++)
                       log_f_c2p[k] = -ptr_risk->covariate[k] * ff;

                    e *= ipow(f, ptr_risk->size);
                    log_e_lb[b_mode] += ptr_risk->size * log_f_lb[b_mode];
                    for(k=0; k<n_c2p_covariate; k++)
                       log_e_c2p[k] += ptr_risk->size * log_f_c2p[k];
                    ptr_risk = ptr_risk->next;
                 }
                 //calculate log_f_* variables related to the probability of escaping risk from infective people.
                 //update log_e_* variables accordingly.
                 //The difference in infectiousness level between symptomatic and asymptomatic cases
                 //is adjusted by one additional covariate that denote the symptom status of the infective person.

                 ptr_risk = person->risk_history[rr].p2p_risk;
                 while(ptr_risk != NULL)
                 {
                    covariate_effect = 0.0;
                    for(k=0; k<n_p2p_covariate; k++)
                       covariate_effect += ptr_risk->covariate[k] * coeff_p2p[k];
                    covariate_effect += ptr_risk->offset;
                                
                    p_mode = ptr_risk->contact_mode;
                    logit_f = lp[p_mode] + covariate_effect;
                    s = ptr_risk->infective_prob;
                    if( s > 0)
                    {
                       ff = inv_logit(logit_f) * s * bipow(asym_effect, 1 - ptr_risk->symptom);
                       f = 1 - ff;
                       factor1 = ff / s;
                       log_f_lp[p_mode] = - ff * (1 - factor1) / f;
                       
                       for(k=0; k<n_p2p_covariate; k++)\
                          log_f_p2p[k] = log_f_lp[p_mode] * ptr_risk->covariate[k];

                       if(f <= 0.0 || f > 1.0)
                       {
                          printf("Goodness-of-fit: f=%e\n", f);
                          free(observed_INF_freq);
                          free(expected_INF_freq);
                          return(1);
                       }
                       e *= ipow(f, ptr_risk->size);
                       log_e_lp[p_mode] += ptr_risk->size * log_f_lp[p_mode];
                       for(k=0; k<n_p2p_covariate; k++)
                          log_e_p2p[k] += ptr_risk->size * log_f_p2p[k];
                    }      
                    ptr_risk = ptr_risk->next;
                 }    /* if(person->p2p_contact_history[r].size > 0) */
                 ee[r] = e;
                 for(j=0; j<n_b_mode; j++) log_ee_lb[r][j] = log_e_lb[j];
                 for(j=0; j<n_p_mode; j++) log_ee_lp[r][j] = log_e_lp[j];
                 for(k=0; k<n_c2p_covariate; k++) log_ee_c2p[r][k] = log_e_c2p[k];
                 for(k=0; k<n_p2p_covariate; k++) log_ee_p2p[r][k] = log_e_p2p[k];
              } /* end of day t */
           }
           if( person->pre_immune == 0)
           {
              if(cfg_pars.adjust_for_left_truncation == 1 && community[h].earliest_idx_day_ill != MISSING)
              {
                 //start_day = max(community[h].day_epi_start, community[h].earliest_idx_day_ill + 1);
                 likelihood_start_day = max(community[h].day_epi_start, community[h].earliest_idx_day_ill - cfg_pars.max_incubation + 1);
                 start_day = likelihood_start_day;
              }   
              else  
              {
                 likelihood_start_day = community[h].day_epi_start;
                 start_day = likelihood_start_day;
              }  
              //if(h == 25)  printf("person id=%d  infection=%d  day ill=%d\n", person->id, person->infection, person->day_ill); 
              if(person->infection == 1)
              {
                 //determine if this person's community is OEM or MCEM community
                 /*
                 if(cfg_pars.EM == 1 && person->size_possible_states > 0)
                 {
                    if(community[h].size_possible_states > 1 && community[h].size_possible_states < cfg_pars.min_size_MCEM)
                    {
                       ptr_stat = community[h].list_states;
                       
                       ptr_stat_max_likelihood = ptr_stat;
                       max_likelihood = ptr_stat->log_L_cur;
                       while(ptr_stat != NULL)
                       {
                          if(ptr_stat->log_L_cur > max_likelihood)
                          {
                             max_likelihood = ptr_stat->log_L_cur;
                             ptr_stat_max_likelihood = ptr_stat;
                          }
                          ptr_stat = ptr_stat->next;
                       }
                       for(i=0; i<community[h].size_impute; i++)
                       {
                          if(community[h].member_impute[i] == person->id)
                             set_state(person, ptr_stat_max_likelihood->states[i]);
                       }      
                    }
                    else if(community[h].size_possible_states >= cfg_pars.min_size_MCEM)
                    {
                       // here we find community-specific maximum weight, not the global maximum weight for combining all communities,
                       // although for inference we used the global weights.
                       ptr_stat = community[h].sample_states;
                       ptr_stat_max_weight = ptr_stat;
                       max_weight = exp(ptr_stat->log_L_cur-ptr_stat->log_L_pre);
                       while(ptr_stat != NULL)
                       {
                          if(exp(ptr_stat->log_L_cur - ptr_stat->log_L_pre) > max_weight)
                          {
                             max_weight = exp(ptr_stat->log_L_cur - ptr_stat->log_L_pre);
                             ptr_stat_max_weight = ptr_stat;
                          }
                          ptr_stat = ptr_stat->next;
                       }
                       for(i=0; i<community[h].size_impute; i++)
                       {
                          if(community[h].member_impute[i] == person->id)
                             set_state(person, ptr_stat_max_weight->states[i]);
                       }      
                    }         
                 }            
                 */
                 inf1 = person->day_ill - cfg_pars.max_incubation;
                 inf2 = person->day_ill - cfg_pars.min_incubation;
                 if(inf2 >= likelihood_start_day)
                 {
                    inf1 = max(inf1, likelihood_start_day);
                    if(inf1 > likelihood_start_day)
                    for(t=likelihood_start_day; t<inf1; t++)
                    {
                       r = t - likelihood_start_day;
                       if(at_risk != NULL)  at_risk[r] += person->weight;
                       rr = t - community[h].day_epi_start;
                       expected_INF_freq[rr] += person->weight * (-log(ee[r]));
                       for(j=0; j<n_b_mode; j++) score_lb[rr][j] += person->weight * (-log_ee_lb[r][j]); 
                       for(j=0; j<n_p_mode; j++) score_lp[rr][j] += person->weight * (-log_ee_lp[r][j]); 
                       for(j=0; j<n_c2p_covariate; j++) score_c2p[rr][j] += person->weight * (-log_ee_c2p[r][j]); 
                       for(j=0; j<n_p2p_covariate; j++) score_p2p[rr][j] += person->weight * (-log_ee_p2p[r][j]); 
                    }   
                    len = inf2 - inf1 + 1;
                    weight = (double *) malloc((size_t) (len * sizeof(double)));
                    cum_e = 1.0;
                    sum_weight = 0.0;
                    for(t=inf1; t<=inf2; t++)
                    {
                       r = t - likelihood_start_day;
                       //p_inf = (fabs(log(ee[r]))<1e-7)? (-log(ee[r])) : (1.0 - ee[r]);
                       p_inf = 1.0 - ee[r];
                       weight[t - inf1] = cum_e * p_inf * pdf_incubation[inf2 - t];
                       cum_e *= ee[r];
                       sum_weight += weight[t - inf1];
                    }
                    for(l=0; l<len; l++)  weight[l] = weight[l]/sum_weight;
                    sum_weight = 0.0;
                    for(t=inf1; t<=inf2; t++)
                    {
                       r = t - likelihood_start_day;
                       if(at_risk != NULL)  at_risk[r] += person->weight;
                       rr = t - community[h].day_epi_start;
                       l = t - inf1;
                       //factor1 = ((1.0 - sum_weight - weight[l]) * 1 + weight[l] * 0.5);
                       factor1 = 1.0 - sum_weight;// prob of at-risk
                       expected_INF_freq[rr] += person->weight * (-log(ee[r])) * factor1;
                       observed_INF_freq[rr] += person->weight * weight[l];
                       for(j=0; j<n_b_mode; j++) score_lb[rr][j] += person->weight * (-log_ee_lb[r][j]) * factor1; 
                       for(j=0; j<n_p_mode; j++) score_lp[rr][j] += person->weight * (-log_ee_lp[r][j]) * factor1; 
                       for(j=0; j<n_c2p_covariate; j++) score_c2p[rr][j] += person->weight * (-log_ee_c2p[r][j]) * factor1; 
                       for(j=0; j<n_p2p_covariate; j++) score_p2p[rr][j] += person->weight * (-log_ee_p2p[r][j]) * factor1; 
                       sum_weight += weight[l];
                    }   
                    free(weight);
                 }
              } 
              else //escape
              {
                 stop_day = min(community[h].day_epi_stop, person->day_exit);
                 if(cfg_pars.adjust_for_right_censoring == 1)
                 {
                    inf1 = stop_day - cfg_pars.max_incubation + 1;
                    inf2 = stop_day - cfg_pars.min_incubation;
                    //if(h == 25)  printf("max_incubation=%d  min_incubation=%d\n", cfg_pars.max_incubation, cfg_pars.min_incubation);
                    //if(h == 25)  printf("start=%d likelihood start=%d  inf1=%d  inf2=%d  stop_day=%d\n", community[h].day_epi_start, likelihood_start_day, inf1, inf2, stop_day);
                    if(inf1 > likelihood_start_day)
                    for(t=likelihood_start_day; t<inf1; t++)
                    {
                       r = t - likelihood_start_day;
                       //if(h == 25)  printf("r=%d\n", r);
                       if(at_risk != NULL)  at_risk[r] += person->weight;
                       rr = t - community[h].day_epi_start;
                       //if(h==25)  printf("t=%d  r=%d  rr=%d\n", t, r, rr);
                       expected_INF_freq[rr] += person->weight *  (-log(ee[r]));
                       for(j=0; j<n_b_mode; j++) score_lb[rr][j] += person->weight * (-log_ee_lb[r][j]); 
                       for(j=0; j<n_p_mode; j++) score_lp[rr][j] += person->weight * (-log_ee_lp[r][j]); 
                       for(j=0; j<n_c2p_covariate; j++) score_c2p[rr][j] += person->weight * (-log_ee_c2p[r][j]); 
                       for(j=0; j<n_p2p_covariate; j++) score_p2p[rr][j] += person->weight * (-log_ee_p2p[r][j]); 
                    }  
                    //if(h==25)  printf("check2\n"); 
                    len = stop_day - inf1 + 2;
                    weight = (double *) malloc((size_t) (len * sizeof(double)));
                    cum_e = 1.0;
                    sum_weight = 0.0;
                    // possible infections with incubation period extending beyond stop_day
                    for(t=inf1; t<=stop_day; t++)
                    {
                       r = t - likelihood_start_day;
                       rr = t - community[h].day_epi_start;
                       //if(h==25)  printf("t=%d  r=%d  rr=%d\n", t, r, rr);
                       //p_inf = (fabs(log(ee[r]))<1e-7)? (-log(ee[r])) : (1.0 - ee[r]);
                       p_inf = 1.0 - ee[r];
                       if(t <= inf2)
                       {
                          weight[t - inf1] = cum_e * p_inf * sdf_incubation[stop_day - t-1];
                          //printf("stop_day-t=%d  sdf=%e\n", stop_day-t, sdf_incubation[stop_day - t-1]);
                       }   
                       else
                          weight[t - inf1] = cum_e * p_inf;   
                       //observed_INF_freq[rr] += person->weight * weight[t - inf1];
                       cum_e *= ee[r];
                       sum_weight += weight[t - inf1];
                    }
                    //if(h==25)  printf("check3\n"); 
                    // full escape
                    weight[len-1] = cum_e;
                    sum_weight += weight[len-1];
                    
                    for(l=0; l<len; l++)  weight[l] = weight[l]/sum_weight;
                    sum_weight = 0.0;
                    for(t=inf1; t<=stop_day; t++)
                    {
                       r = t - likelihood_start_day;
                       if(at_risk != NULL)  at_risk[r] += person->weight;
                       rr = t - community[h].day_epi_start;
                       //if(h==25)  printf("t=%d  r=%d  rr=%d\n", t, r, rr);
                       l = t - inf1;
                       //factor1 = ((1.0 - sum_weight - weight[l]) * 1 + weight[l] * 0.5);
                       factor1 = 1.0 - sum_weight; //prob of at-risk
                       expected_INF_freq[rr] += person->weight * (-log(ee[r])) * factor1;
                       observed_INF_freq[rr] += person->weight * weight[l];
                       for(j=0; j<n_b_mode; j++) score_lb[rr][j] += person->weight * (-log_ee_lb[r][j]) * factor1; 
                       for(j=0; j<n_p_mode; j++) score_lp[rr][j] += person->weight * (-log_ee_lp[r][j]) * factor1; 
                       for(j=0; j<n_c2p_covariate; j++) score_c2p[rr][j] += person->weight * (-log_ee_c2p[r][j]) * factor1; 
                       for(j=0; j<n_p2p_covariate; j++) score_p2p[rr][j] += person->weight * (-log_ee_p2p[r][j]) * factor1; 
                       sum_weight += weight[l];
                    } 
                    //if(h==25)  printf("check4\n"); 
                    // when t=stop_day, 1.0 - sum_weight already accounts for the weight of full escape.
                    free(weight);
                 }  
                 else
                 {
                    for(t=likelihood_start_day; t<=stop_day; t++)  
                    {
                       r = t - likelihood_start_day;
                       if(at_risk != NULL)  at_risk[r] += person->weight;
                       rr = t - community[h].day_epi_start;
                       expected_INF_freq[rr] += person->weight * (-log(ee[r]));
                       for(j=0; j<n_b_mode; j++) score_lb[rr][j] += person->weight * (-log_ee_lb[r][j]); 
                       for(j=0; j<n_p_mode; j++) score_lp[rr][j] += person->weight * (-log_ee_lp[r][j]); 
                       for(j=0; j<n_c2p_covariate; j++) score_c2p[rr][j] += person->weight * (-log_ee_c2p[r][j]); 
                       for(j=0; j<n_p2p_covariate; j++) score_p2p[rr][j] += person->weight * (-log_ee_p2p[r][j]); 
                    }
                 }
              }      
           } //if( person->pre_immune == 0)
           //printf("person %d: my_log_L=%e\n", person->id, my_log_L);
           ptr_integer = ptr_integer->next;
        } //while(ptr_integer != NULL)
        ptr_class = ptr_class->next;
     } //while(ptr_class != NULL)
     for(j=0; j<length; j++)
     {
        shift = 0;
        for(i=0; i<n_b_mode; i++)  first->data[j][shift + i] = score_lb[j][i];
        shift = n_b_mode;
        for(i=0; i<n_p_mode; i++)  first->data[j][shift + i] = score_lp[j][i];
        shift = n_b_mode + n_p_mode + n_u_mode + n_q_mode;
        for(i=0; i<n_c2p_covariate; i++)  first->data[j][shift + i] = score_c2p[j][i];
        shift = n_b_mode + n_p_mode + n_u_mode + n_q_mode + n_c2p_covariate;
        for(i=0; i<n_p2p_covariate; i++)  first->data[j][shift + i] = score_p2p[j][i];
        shift = n_b_mode + n_p_mode + n_u_mode + n_q_mode + n_c2p_covariate + n_p2p_covariate;
     }
     
     if(n_b_mode > 0)  free_2d_array_double(score_lb); 
     if(n_p_mode > 0)  free_2d_array_double(score_lp); 
     if(n_c2p_covariate > 0)  free_2d_array_double(score_c2p); 
     if(n_p2p_covariate > 0)  free_2d_array_double(score_p2p); 
  }
  return(0);
}

// For goodness of fit, we assume all communities are aligned, starting from the same time. If it's case-ascertained design,
// then the onset times of index cases are algined. Although this not necessarily the reality, it appears the alignment is necessary.
// For example, if each community is a housheold, and the housheolds have outbreak times all different from each other,
// then the housheold-specific observed onset frequencies and model-fitted onset frquencies will be very differnt from each other, because the sample size is very small.
// For goodness of fit, we care about the overall performance of the model, and it is reasonable to align all the households 
// and assess the sum of observed and fitted fruqencies over all households using the same timeline, as these households are independent.
// aux_community should not be empty. We assume all importance weights are 1 in this program, and therefore importance samples
// need to be generated using the final parameter estimates.
int goodness_of_fit(int id_inc, int id_inf, int id_time, double *est, COMMUNITY *aux_community, MATRIX *var_logit)
{
   int i, j, h, k, l, m, n, r, t, count, max_epi_duration;
   int n_b_mode, n_p_mode, n_q_mode, n_u_mode;
   int n_c2p_covariate, n_p2p_covariate, n_imm_covariate, n_pat_covariate;
   int n_covariate, n_time_ind_covariate, n_time_dep_covariate;
   int n_sus_p2p_covariate, n_inf_p2p_covariate, n_int_p2p_covariate;
   int n_par, n_par_equiclass, n_sim,  *at_risk; 
   int error = 0;
   double var_exp_INF_freq, se_log, mean_inc;
   double *obs_INF_freq, *exp_INF_freq, *lower, *upper, *sim_lower, *sim_upper;
   double *sum_obs_INF_freq, *sum_exp_INF_freq, *com_obs_INF_freq, *com_exp_INF_freq;
   double *par, *par_effective, *new_par, *new_par_effective, *new_est;
   double total_obs_INF_freq, total_exp_INF_freq;
   double tmp1, tmp2;
   MATRIX first, score, sim_INF_freq, subvar_logit;
   STATE *ptr_stat;
   PEOPLE *person;
   char file_name[100];
   FILE *file, *file1;

   n_b_mode = cfg_pars.n_b_mode;
   n_p_mode = cfg_pars.n_p_mode;
   n_u_mode = cfg_pars.n_u_mode;
   n_q_mode = cfg_pars.n_q_mode;
   n_time_ind_covariate = cfg_pars.n_time_ind_covariate;
   n_time_dep_covariate = cfg_pars.n_time_dep_covariate;
   n_covariate = cfg_pars.n_covariate;
   n_c2p_covariate = cfg_pars.n_c2p_covariate;
   n_sus_p2p_covariate = cfg_pars.n_sus_p2p_covariate;
   n_inf_p2p_covariate = cfg_pars.n_inf_p2p_covariate;
   n_int_p2p_covariate = cfg_pars.n_int_p2p_covariate;
   n_p2p_covariate = cfg_pars.n_p2p_covariate;
   n_pat_covariate = cfg_pars.n_pat_covariate;
   n_imm_covariate = cfg_pars.n_imm_covariate;
   n_par = cfg_pars.n_par;
   n_par_equiclass = cfg_pars.n_par_equiclass;
   
   max_epi_duration = community[0].epi_duration;
   for(h=1; h<n_community; h++)  max_epi_duration = max(community[h].epi_duration, max_epi_duration);
   //printf("max_epi_duration=%d\n", max_epi_duration);
   make_1d_array_double(&obs_INF_freq, max_epi_duration, 0.0);
   make_1d_array_double(&exp_INF_freq, max_epi_duration, 0.0);
   make_1d_array_double(&sum_obs_INF_freq, max_epi_duration, 0.0);
   make_1d_array_double(&sum_exp_INF_freq, max_epi_duration, 0.0);
   make_1d_array_double(&com_obs_INF_freq, max_epi_duration, 0.0);
   make_1d_array_double(&com_exp_INF_freq, max_epi_duration, 0.0);
   make_1d_array_double(&lower, max_epi_duration, 0.0);
   make_1d_array_double(&upper, max_epi_duration, 0.0);
   make_1d_array_double(&sim_lower, max_epi_duration, 0.0);
   make_1d_array_double(&sim_upper, max_epi_duration, 0.0);
   //at_risk only used for Rt calculation using epidemic curve data, wich typically has a single community and does not involve EM-MCEM
   make_1d_array_int(&at_risk, max_epi_duration, 0); 
   make_1d_array_double(&par, n_par, 0.0);
   make_1d_array_double(&par_effective, n_par_equiclass, 0.0);
   make_1d_array_double(&new_par, n_par, 0.0);
   make_1d_array_double(&new_par_effective, n_par_equiclass, 0.0);
   make_1d_array_double(&new_est, n_par, 0.0);
   
   initialize_matrix(&subvar_logit);
   initialize_matrix(&first);
   initialize_matrix(&score);
   inflate_matrix(&subvar_logit, n_par_equiclass, n_par_equiclass, 0);
   inflate_matrix(&first, max_epi_duration, n_par, 0);
   inflate_matrix(&score, max_epi_duration, n_par, 0);
   
   
   if(cfg_pars.EM == 0)  size_sample_states = 1;
   else
   {
      if(aux_community == NULL)
      {
         printf("Goodness-of-fit for EM-MCEM: new importance samples are always needed. Auxilliary community cannot be empty\n");
         error = 1;
         goto end;
      }
   }
   //print_2d_int(people[295].possible_states, 2, people[295].size_possible_states);
   //printf("size_sample_states=%d\n", size_sample_states);
   //sprintf(file_name, "%scommunity_obs_%d.txt", cfg_pars.path_out, serial_number);
   //file = fopen(file_name, "w");
   //sprintf(file_name, "%scommunity_fitted_%d.txt", cfg_pars.path_out, serial_number);
   //file1 = fopen(file_name, "w");
   for(h=0; h<n_community; h++)
   {
      if(community[h].size > 0 && community[h].ignore == 0)
      {
         reset_1d_array_double(com_obs_INF_freq, max_epi_duration, 0.0);
         reset_1d_array_double(com_exp_INF_freq, max_epi_duration, 0.0);
         //printf("h=%d  size=%d  size_impute=%d\n", h, community[h].size, community[h].size_impute);
         if(community[h].size_impute == 0) //when not EM, every community should have size_impute=0
         {
            reset_1d_array_double(obs_INF_freq, max_epi_duration, 0.0);
            reset_1d_array_double(exp_INF_freq, max_epi_duration, 0.0);
            reset(&first, 0.0);
            if(goodness_of_fit_community(h, est, at_risk, obs_INF_freq, exp_INF_freq, &first) == 1)  error = 1;
            addition(&score, 1.0, &first, (double) size_sample_states, &score);
            for(k=0; k<max_epi_duration; k++) 
            {   
               sum_obs_INF_freq[k] += size_sample_states * obs_INF_freq[k];
               sum_exp_INF_freq[k] += size_sample_states * exp_INF_freq[k];
               com_obs_INF_freq[k] += size_sample_states * obs_INF_freq[k];
               com_exp_INF_freq[k] += size_sample_states * exp_INF_freq[k];
            }
         }
         else
         {
            ptr_stat = aux_community[h].sample_states;
            while(ptr_stat != NULL)
            {
               // go over members with uncertain state, and assign sample state
               //if(h==93) printf("log_L_cur=%e: ", ptr_stat->log_L_cur);
               for(i=0; i<community[h].size_impute; i++)
               {
                  person = people + community[h].member_impute[i];
                  set_state(person, ptr_stat->states[i]);
                  //if(h==93) printf("  %d:%d", person->id, ptr_stat->states[i]);
               }     
               //if(h==93) printf("\n");
               if(cfg_pars.adjust_for_left_truncation == 1 && cfg_pars.preset_index == 0)   set_index_cases(community + h);
               reset_1d_array_double(obs_INF_freq, max_epi_duration, 0.0);
               reset_1d_array_double(exp_INF_freq, max_epi_duration, 0.0);
               reset(&first, 0.0);
            
               if(goodness_of_fit_community(h, est, NULL, obs_INF_freq, exp_INF_freq, &first) == 1)  error = 1;
            
               addition(&score, 1.0, &first, 1.0, &score);
               for(k=0; k<max_epi_duration; k++) 
               {    
                  sum_obs_INF_freq[k] += obs_INF_freq[k];
                  sum_exp_INF_freq[k] += exp_INF_freq[k];
                  com_obs_INF_freq[k] += obs_INF_freq[k];
                  com_exp_INF_freq[k] += exp_INF_freq[k];
               }
               ptr_stat = ptr_stat->next;
            }
         }
         //fprintf(file, "%3d\n", h);
         //for(k=0; k<max_epi_duration; k++)
         //   fprintf(file, "  %10e(%3d)", com_obs_INF_freq[k], k+1);
         //fprintf(file, "\n");
         //fprintf(file1, "%3d\n", h);
         //for(k=0; k<max_epi_duration; k++)
         //   fprintf(file1, "  %12.3e(%3d)", com_exp_INF_freq[k], k+1);
         //fprintf(file1, "\n");
      }   
      /*if(h == 154)
      {
         printf("h=%3d  size=%3d  size_impute=%3d", h, community[h].size, community[h].size_impute);  
         printf("  sum_obs=%e  sum_exp=%e\n", sum_double(obs_INF_freq, max_epi_duration)/size_sample_states, sum_double(exp_INF_freq, max_epi_duration)/size_sample_states);
         printf("  cum_sum_obs=%e  cum_sum_exp=%e\n", sum_double(sum_obs_INF_freq, max_epi_duration)/size_sample_states, sum_double(sum_exp_INF_freq, max_epi_duration)/size_sample_states);
      } */  
   }
   //fclose(file);
   //fclose(file1);
   scalar_addition(&score, 0.0, 1.0/size_sample_states, &score);
   total_obs_INF_freq = total_exp_INF_freq = 0;
   for(k=0; k<max_epi_duration; k++) 
   { 
      //printf("k=%d: obs=%e,  exp=%e, sum_obs=%e,  sum_exp=%e\n", k, obs_INF_freq[k], exp_INF_freq[k], sum_obs_INF_freq[k], sum_exp_INF_freq[k]);
      sum_obs_INF_freq[k] /= size_sample_states;
      sum_exp_INF_freq[k] /= size_sample_states;
      total_obs_INF_freq += sum_obs_INF_freq[k];
      total_exp_INF_freq += sum_exp_INF_freq[k];
   }
   if(cfg_pars.silent_run == 0)  printf("total_obs_INF_freq=%e,  total_exp_INF_freq=%e\n", total_obs_INF_freq, total_exp_INF_freq);
   if(cfg_pars.skip_variance == 0)
   {
      for(k=0; k<max_epi_duration; k++) 
      {
         if(sum_exp_INF_freq[k] > 0)
         {
            var_exp_INF_freq = 0;
            for(m=0; m<n_par; m++)
            {
               for(n=0; n<n_par; n++)
               {
                  var_exp_INF_freq += score.data[k][m] * var_logit->data[m][n] * score.data[k][n];
               }
            }
            se_log = sqrt(var_exp_INF_freq) / sum_exp_INF_freq[k];
            lower[k] = exp(log(sum_exp_INF_freq[k]) - 1.96 * se_log);     
            upper[k] = exp(log(sum_exp_INF_freq[k]) + 1.96 * se_log);     
         }  
         else  lower[k] = upper[k] = -1;   
      }
   }

   // use estimated parameter to simulate epidemics to see whether they cover observed data
   
   for(j=0; j<n_par; j++)
   {
      if(j < n_b_mode + n_p_mode + n_u_mode + n_q_mode)  par[j] = logit(est[j]);
      else par[j] = log(est[j]);
   }
   for(j=0; j<n_par_equiclass; j++)
   {
      m = cfg_pars.par_equiclass[j].member[0] - 1;
      par_effective[j] = par[m] ;
   }
   for(i=0; i<n_par_equiclass; i++)
   {
      for(j=i; j<n_par_equiclass; j++)
      {
         m = cfg_pars.par_equiclass[i].member[0] - 1;
         n = cfg_pars.par_equiclass[j].member[0] - 1;
         subvar_logit.data[i][j] = var_logit->data[m][n];
         subvar_logit.data[j][i] = var_logit->data[m][n];
      }
   }
   //count_cases();
   n_sim = 0;
   initialize_matrix(&sim_INF_freq);
   if(n_sim > 0)   inflate_matrix(&sim_INF_freq, max_epi_duration, n_sim, 0);
   for(k=0; k<n_sim; k++)
   {
      if(cfg_pars.silent_run == 0) printf("%d\n", k);
      // generate new parameters from multinormal with mean as the point estimates and covariance matrix as the covariance estimates
      rmultinorm (new_par_effective, n_par_equiclass, par_effective, &subvar_logit, &seed);
      // convert effective parameters to the full parameters
      for(j=0; j<n_par_equiclass; j++)
      {
         for(l=0; l<cfg_pars.par_equiclass[j].size; l++)
         {
            m = cfg_pars.par_equiclass[j].member[l] - 1;
            new_par[m] = new_par_effective[j];
         }
      }

      if(cfg_pars.n_par_fixed > 0)
      {
         for(i=0; i<cfg_pars.n_par_fixed; i++)
         {
            j = cfg_pars.par_fixed_id[i] - 1;
            new_par[j] =  cfg_pars.par_fixed_value[i];
         }
      }
      for(j=0; j<n_par; j++)
      {
         if(j < n_b_mode + n_p_mode + n_u_mode + n_q_mode)  new_est[j] = inv_logit(new_par[j]);
         else new_est[j] = exp(new_par[j]);
      }
      
      simulate(1, new_est);
      for(h=0; h<n_community; h++)
      if(community[h].size > 0 && community[h].ignore == 0)
      {          
         for(m=0; m<community[h].size; m++)
         {
            i = community[h].member[m];
            if(people[i].ignore == 0 && people[i].infection == 1)
            {
               r = people[i].day_infection - community[h].day_epi_start;
               if(cfg_pars.adjust_for_left_truncation == 1)
               {
                  if(people[i].idx != 1)  sim_INF_freq.data[r][k] ++;
               }   
               else sim_INF_freq.data[r][k] ++;
            }
         }
      }
   }
   if(n_sim > 1)
   {
      for(r=0; r<max_epi_duration; r++)
      {
         sim_lower[r] = quantile_double(sim_INF_freq.data[r], n_sim, 0.025);
         sim_upper[r] = quantile_double(sim_INF_freq.data[r], n_sim, 0.975);
      }   
   }   
   //sprintf(file_name, "%soutput_fitted_sim_%d.txt", cfg_pars.path_out, serial_number);
   //if((file = fopen(file_name, "a")) == NULL)
   //   file = fopen(file_name, "w");
   //for(k=0; k<max_epi_duration; k++)
   //{
   //   fprintf(file, "%3d  %3d   %5d", id_inc, id_inf, k);
   //   for(i=0; i<n_sim; i++)  fprintf(file, "  %18.3e", sim_INF_freq.data[k][i]);
   //   fprintf(file, "\n");
   //}   
   //fclose(file);                
   sprintf(file_name, "%soutput_goodness.txt", cfg_pars.path_out);
   if((file = fopen(file_name, "a")) == NULL)
      file = fopen(file_name, "w");
   if(cfg_pars.R0_divide_by_time == 0)
   {
      for(k=0; k<max_epi_duration; k++)
      {
         fprintf(file, "%6d  %3d  %3d  %5d  %15e  %18.3e", serial_number, id_inc, id_inf, k, sum_obs_INF_freq[k], sum_exp_INF_freq[k]);
         //if(cfg_pars.skip_variance == 0)
         //   fprintf(file, "  %18.3e  %18.3e", lower[k], upper[k]);
         if(n_sim > 1)
            fprintf(file, "  %18.3e  %18.3e", sim_lower[k], sim_upper[k]);
         fprintf(file, "\n");   
      }
   }
   else
   {
      mean_inc = 0.0;
      for(j=cfg_pars.min_incubation; j<=cfg_pars.max_incubation; j++)
         mean_inc += cfg_pars.prob_incubation[j - cfg_pars.min_incubation] * j;
      t = id_time - floor(cfg_pars.R0_window_size/2) + lround(mean_inc);
      k = t - community[0].day_epi_start;
      fprintf(file, "%6d  %3d  %3d  %6d  %15d  %18.3e  %18.3e\n", serial_number, id_inc, id_inf, t, at_risk[k], sum_obs_INF_freq[k], sum_exp_INF_freq[k]);
   }
   fclose(file);
   

end:
   deflate_matrix(&subvar_logit);
   deflate_matrix(&first);
   deflate_matrix(&score);
   deflate_matrix(&sim_INF_freq);
   free(obs_INF_freq);
   free(exp_INF_freq);
   free(sum_obs_INF_freq);
   free(sum_exp_INF_freq);
   free(com_obs_INF_freq);
   free(com_exp_INF_freq);
   free(lower);
   free(upper);
   free(sim_lower);
   free(sim_upper);
   free(at_risk);
   free(par);
   free(par_effective);
   free(new_par);
   free(new_par_effective);
   free(new_est);
   return(error);
}


