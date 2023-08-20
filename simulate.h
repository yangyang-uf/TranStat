
int simulate(int goodness_of_fit, double *est) // est is vector of parameter estimates on their original scale, i.e., not log- or logit-transformed
{
  FILE *file;     /* file pointer, needed to handle i/o to files */
  int loop, id, h, i, j, k, l, m, n, r, t, rr, tt;
  int resetable, length, count, start_day, stop_day;
  int n_exp_c, n_inf_c, n_exp_a, n_inf_a;
  int n_imm, n_sym_idx, n_asym_idx, n_sym_sec, n_asym_sec, n_esc;
  int b_mode, p_mode, q_mode, u_mode;
  int n_b_mode, n_p_mode, n_q_mode, n_u_mode;
  int n_c2p_covariate, n_p2p_covariate, n_imm_covariate, n_pat_covariate;
  int n_covariate, n_time_ind_covariate, n_time_dep_covariate;
  int n_sus_p2p_covariate, n_inf_p2p_covariate, n_int_p2p_covariate;
  int n_par, n_par_equiclass;

  double prob_esc, prob_imm, prob_sym, s;
  double covariate_effect, temp, *par; 
  double *c2p_covariate, *p2p_covariate;
  CONTACT *ptr_contact;

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

  if(goodness_of_fit == 0)
  {
     make_1d_array_double(&par, n_par, 0.0);
     for(i=0; i<n_par_equiclass; i++)
     {
        for(j=0; j<cfg_pars.par_equiclass[i].size; j++)
        {
           m = cfg_pars.par_equiclass[i].member[j] - 1;
           par[m] = cfg_pars.sim_par_effective[i];
        }
     }
  
     //if fixed parameters exist, need to set fixed parameters to prespecified values
     if(cfg_pars.n_par_fixed > 0)
     {
        for(i=0; i<cfg_pars.n_par_fixed; i++)
        {
           j = cfg_pars.par_fixed_id[i] - 1;
           par[j] =  cfg_pars.par_fixed_value[i];
        }
     }

     for(i=0; i<n_b_mode; i++)  {lb[i] = par[i]; b[i] = inv_logit(lb[i]);}
     for(i=0; i<n_p_mode; i++)  {lp[i] = par[n_b_mode + i]; p[i] = inv_logit(lp[i]);}
     for(i=0; i<n_u_mode; i++)  {lu[i] = par[n_b_mode + n_p_mode + i]; u[i] = inv_logit(lu[i]);}
     for(i=0; i<n_q_mode; i++)  {lq[i] = par[n_b_mode + n_p_mode + n_u_mode + i]; q[i] = inv_logit(lq[i]);}
     for(i=0; i<n_c2p_covariate; i++)  coeff_c2p[i] = par[n_b_mode + n_p_mode + n_u_mode + n_q_mode + i];
     for(i=0; i<n_p2p_covariate; i++)  coeff_p2p[i] = par[n_b_mode + n_p_mode + n_u_mode + n_q_mode
                                                    + n_c2p_covariate + i];
     for(i=0; i<n_pat_covariate; i++)  coeff_pat[i] = par[n_b_mode + n_p_mode + n_u_mode + n_q_mode
                                                    + n_c2p_covariate + n_p2p_covariate + i];
     for(i=0; i<n_imm_covariate; i++)  coeff_imm[i] = par[n_b_mode + n_p_mode + n_u_mode + n_q_mode
                                                    + n_c2p_covariate + n_p2p_covariate + n_pat_covariate + i];
     free(par);
  }
  else if(est != NULL)
  {
     for(i=0; i<n_b_mode; i++)  {b[i] = est[i]; lb[i] = logit(b[i]); }
     for(i=0; i<n_p_mode; i++)  {p[i] = est[n_b_mode + i]; lp[i] = logit(p[i]); }
     for(i=0; i<n_u_mode; i++)  {u[i] = est[n_b_mode + n_p_mode + i]; lu[i] = logit(u[i]); } 
     for(i=0; i<n_q_mode; i++)  {q[i] = est[n_b_mode + n_p_mode + n_u_mode + i]; lq[i] = logit(q[i]); }
     for(i=0; i<n_c2p_covariate; i++)  coeff_c2p[i] = log(est[n_b_mode + n_p_mode + n_u_mode + n_q_mode + i]);
     for(i=0; i<n_p2p_covariate; i++)  coeff_p2p[i] = log(est[n_b_mode + n_p_mode + n_u_mode + n_q_mode
                                                    + n_c2p_covariate + i]);
     for(i=0; i<n_pat_covariate; i++)  coeff_pat[i] = log(est[n_b_mode + n_p_mode + n_u_mode + n_q_mode
                                                    + n_c2p_covariate + n_p2p_covariate + i]);
     for(i=0; i<n_imm_covariate; i++)  coeff_imm[i] = log(est[n_b_mode + n_p_mode + n_u_mode + n_q_mode
                                                    + n_c2p_covariate + n_p2p_covariate + n_pat_covariate + i]);
  }    

  c2p_covariate = NULL;
  if(n_c2p_covariate > 0)  
      c2p_covariate = (double *)malloc((size_t) (n_c2p_covariate * sizeof(double)));

  p2p_covariate = NULL;
  if(n_p2p_covariate > 0)  
      p2p_covariate = (double *)malloc((size_t) (n_p2p_covariate * sizeof(double)));
 
  for(i=0; i < p_size; i++)
  if(people[i].ignore == 0)
  {
     h = people[i].community;
     resetable = 1;
     if(cfg_pars.adjust_for_left_truncation == 1 && people[i].idx == 1)  resetable = 0;
     if(resetable == 1)
     {
        people[i].pre_immune = people[i].infection = people[i].symptom = 0;
        people[i].day_infection = people[i].day_ill = MISSING;
        people[i].day_infection_lower = people[i].day_infection_upper = MISSING;
        people[i].day_infective_lower = people[i].day_infective_upper = MISSING;
        if(cfg_pars.n_time_dep_covariate > 0 && cfg_pars.illness_as_covariate == 1)
        {
           j = cfg_pars.illness_covariate_id - 1;
           for(k=0; k<community[h].epi_duration; k++)
              people[i].time_dep_covariate[k][j] = 0.0;
        }      
        //Reset time-independent indicators for intervals of infectivity for infected individuals
        //Only applies to the Wuhan housheold study
        /*
        people[i].time_ind_covariate[10-1] = 0;
        people[i].time_ind_covariate[11-1] = 0;
        people[i].time_ind_covariate[12-1] = 0;
        people[i].time_ind_covariate[13-1] = 0;
        people[i].time_ind_covariate[14-1] = 0;
        people[i].time_ind_covariate[15-1] = 0;
        */
        if(n_q_mode > 0)
        {
           q_mode = people[i].q_mode;
           covariate_effect = 0.0;
           for(k=0; k<n_imm_covariate; k++)
              covariate_effect += people[i].imm_covariate[k] * coeff_imm[k];
           prob_imm = inv_logit(lq[q_mode] + covariate_effect);
           if(runiform(&seed) <= prob_imm)  people[i].pre_immune = 1;
        }
     }
  }
  
  /************************************************************
   at this point, all is setup to begin the dynamics.
   start the do loop that will run until the epidemic is over
   ************************************************************/
  //printf("b=%e  p=%e  lb=%e  lp=%e  coeff_c2p=%e  coeff_p2p=%e\n", b[0], p[0], lb[0], lp[0], coeff_c2p[0], coeff_p2p[0]);
  for(h=0; h<n_community; h++)
  if(community[h].ignore == 0)
  {
     if(goodness_of_fit == 1)
     {
        if(cfg_pars.adjust_for_left_truncation == 1 && community[h].earliest_idx_day_ill != MISSING)
           start_day = max(community[h].day_epi_start, community[h].earliest_idx_day_ill - cfg_pars.max_incubation + 1);
        else 
           start_day = community[h].day_epi_start;  
     }       
     else  start_day = community[h].day_epi_start;
     stop_day = community[h].day_epi_stop;
     for(t=start_day; t <= stop_day; t++)
     {
        for(m=0; m<community[h].size; m++)
        {
           i = community[h].member[m];
           //printf("h=%d  t=%d  i=%d\n", h, t, i);
           if(n_u_mode > 0)  u_mode = people[i].u_mode;
           if(n_q_mode > 0)  q_mode = people[i].q_mode;
           if(people[i].ignore == 0 && people[i].day_infection == MISSING && people[i].pre_immune == 0 && people[i].idx != 1 && t <= people[i].day_exit)
           {
              prob_esc = 1.0;
              r = t - community[h].day_epi_start;
              if(n_c2p_covariate > 0)  
                 organize_c2p_covariate(t, people+i, c2p_covariate);
              if(cfg_pars.common_contact_history_within_community == 1)
                 ptr_contact = community[h].contact_history[r].c2p_contact;
              else
                 ptr_contact = people[i].contact_history[r].c2p_contact;
              while(ptr_contact != NULL)
              {
                 covariate_effect = 0.0;
                 for(k=0; k<n_c2p_covariate; k++)
                    covariate_effect += c2p_covariate[k] * coeff_c2p[k];
                 covariate_effect += ptr_contact->offset;

                 b_mode = ptr_contact->contact_mode;
                 prob_esc *= 1.0 - inv_logit(lb[b_mode] + covariate_effect);
                 ptr_contact = ptr_contact->next;
              }
              //printf("b_mode=%d  b=%e  lb=%e  inv_logit=%e  covariate_effect=%e  prob_esc=%e\n", 
              //        b_mode, b[b_mode], lb[b_mode], inv_logit(lb[b_mode]), covariate_effect, prob_esc);
              if(cfg_pars.common_contact_history_within_community == 1)
                 ptr_contact = community[h].contact_history[r].p2p_contact;
              else
                 ptr_contact = people[i].contact_history[r].p2p_contact;
              while(ptr_contact != NULL)
              {
                 j = ptr_contact->contact_id;
                 if(people[j].day_infection != MISSING && t >= people[j].day_infective_lower 
                                                       && t <= people[j].day_infective_upper )
                 {
                    s = cfg_pars.prob_infectious[t - people[j].day_infective_lower];
                    if( s > 0)
                    {
                       organize_p2p_covariate(t, people+i, people+j, p2p_covariate);
                       covariate_effect = 0.0;
                       for(k=0; k<n_p2p_covariate; k++)
                          covariate_effect += p2p_covariate[k] * coeff_p2p[k];
                       //if(people[j].id == 0)
                       //   printf("sus=%d  t=%d  sus_x=%f  contact_x=%f\n", people[i].id, t, people[i].time_ind_covariate[0], ptr_contact->covariate[0]); 
                       covariate_effect += ptr_contact->offset;

                       p_mode = ptr_contact->contact_mode;
                       temp = 1.0 - s * bipow(cfg_pars.asym_effect_sim, 1 - people[j].symptom)
                                 * inv_logit(lp[p_mode] + covariate_effect);
                       prob_esc *= temp;
                       //printf("i=%d  j=%d  symptom=%d  asym_effect_sim=%e  p=%e  lp=%e  covariate_effect=%e  p_effect=%e  s=%e p_esc=%e\n", 
                       //        people[i].id, people[j].id, people[j].symptom, cfg_pars.asym_effect_sim, p[p_mode], lp[p_mode], 
                       //        covariate_effect, inv_logit(lp[p_mode]+covariate_effect), s, temp);          
                    }
                 }
                 ptr_contact = ptr_contact->next;
              }    /* if(person->p2p_contact_history[r].size > 0) */

              //printf("h=%d t=%d i=%d x=%f prob_inf=%e day_infection=%d\n", 
              //        people[i].community, t, i, people[i].time_ind_covariate[0], 1 - prob_esc, people[i].day_infection);
              
              temp = runiform(&seed);
              if(temp < 1 - prob_esc)
              {
                 people[i].infection = 1;
                 people[i].day_infection = t;
                 if(n_u_mode > 0)
                 {
                    if(u[u_mode] == 0 || u[u_mode] == 1)  prob_sym = u[u_mode];
                    else
                    {
                       covariate_effect = 0.0;
                       for(k=0; k<n_pat_covariate; k++)
                          covariate_effect += people[i].pat_covariate[k][r] * coeff_pat[k];
                       prob_sym = inv_logit(lu[u_mode] + covariate_effect);
                    }
                 }
                 else  prob_sym = 1.0;
                 if(runiform(&seed) < prob_sym)  people[i].symptom = 1;
                 length = cfg_pars.max_incubation - cfg_pars.min_incubation + 1;
                 people[i].day_infection = t;
                 people[i].day_ill = t + cfg_pars.min_incubation + discrete(length, cfg_pars.prob_incubation, &seed);
                 people[i].day_infective_lower = people[i].day_ill + cfg_pars.lower_infectious;
                 people[i].day_infective_upper = people[i].day_ill + cfg_pars.upper_infectious;
                 people[i].day_infection_lower = max(people[i].day_ill - cfg_pars.max_incubation, community[h].day_epi_start);
                 people[i].day_infection_upper  = max(people[i].day_ill - cfg_pars.min_incubation, community[h].day_epi_start);

                 if(people[i].symptom == 1) //special scenario for COVID-19, remake time-dependent covariates
                 {
                    if(cfg_pars.n_time_dep_covariate > 0 && cfg_pars.illness_as_covariate == 1)
                    {
                       j = cfg_pars.illness_covariate_id - 1;
                       for(tt = people[i].day_ill; tt <= people[i].day_infective_upper; tt++)
                       if(tt >= start_day && tt <= stop_day)
                       {
                          rr = tt - community[h].day_epi_start;
                          people[i].time_dep_covariate[rr][j] = 1.0;
                       }   
                    }
                 }       
                 //the following is only for the Wuhan household analysis project
                 /*
                 if(people[i].day_ill > community[h].day_epi_stop - 88 && people[i].day_ill <= community[h].day_epi_stop)
                 {
                    if(people[i].day_ill <= community[h].day_epi_stop - 70) 
                       people[i].time_ind_covariate[11-1] = 1;
                    else  
                       people[i].time_ind_covariate[12-1] = 1;
                 }
                 if(people[i].symptom == 0 && people[i].day_ill <= community[h].day_epi_stop) //special scenario for COVID-19, remake time-dependent covariates
                 {
                    if(people[i].day_ill <= community[h].day_epi_stop - 58)
                       people[i].time_ind_covariate[10-1] = 1;
                    else if(people[i].day_ill <= community[h].day_epi_stop - 37) 
                       people[i].time_ind_covariate[13-1] = 1;
                    else if(people[i].day_ill <= community[h].day_epi_stop - 16) 
                       people[i].time_ind_covariate[14-1] = 1;
                    else
                       people[i].time_ind_covariate[15-1] = 1;
                 }   
                 */   
                 //printf("h=%d t=%d i=%d x=%f prob_inf=%e symptom=%d day_ill=%d\n", 
                 //        people[i].community, t, i, people[i].time_ind_covariate[0], 1 - prob_esc, people[i].symptom, people[i].day_ill);
              }
           }
        }
     }
  } 
  /* The folowing code is used when extracting case-ascertained housheolds form prespectively simulated outbreak
  if(cfg_pars.adjust_for_left_truncation == 1)
  {
  
     for(h=0; h < n_community; h++)
     {
        community[h].earliest_idx_day_ill = MISSING;
        community[h].latest_idx_day_ill = MISSING;
        community[h].size_idx = 0;
        for(j=0; j<community[h].size; j++)
        {
           i = community[h].member[j];
           if(people[i].symptom == 1) 
           {
              if(community[h].earliest_idx_day_ill == MISSING) 
                 community[h].earliest_idx_day_ill = people[i].day_ill;
              else
                 community[h].earliest_idx_day_ill = min(community[h].earliest_idx_day_ill, people[i].day_ill);    
              if(community[h].latest_idx_day_ill == MISSING) 
                 community[h].latest_idx_day_ill = people[i].day_ill;
              else
                 community[h].latest_idx_day_ill = max(community[h].latest_idx_day_ill, people[i].day_ill);    
           }
        }   
        for(j=0; j<community[h].size; j++)
        {
           i = community[h].member[j];
           if(people[i].symptom == 1 && people[i].day_ill == community[h].earliest_idx_day_ill) 
           {
              community[h].size_idx ++;
              people[i].idx = 1;              
           }
        }   
        if(community[h].size_idx == 0)
        {
           community[h].ignore = 1;
           for(j=0; j<community[h].size; j++)
           {
              i = community[h].member[j];
              people[i].ignore = 1;
           }
        }
        else
        {
           community[h].idx = (int *) malloc((size_t) community[h].size_idx * sizeof(int));
           community[h].counter = 0;
           for(j=0; j<community[h].size; j++)
           {
              i = community[h].member[j];
              if(people[i].idx == 1)  
              {
                 community[h].idx[community[h].counter] = i;
                 community[h].counter++;
              }
           }
        }
     }
  }   */   
  n_sym_idx = n_asym_idx = n_sym_sec = n_asym_sec = n_esc = n_imm = 0;
  for(i=0; i < p_size; i++)
  if(people[i].ignore == 0)
  {
     if(people[i].pre_immune == 1)  n_imm++;
     else if(people[i].infection == 0)  n_esc++;
     else
     {
        if(people[i].idx == 1)
        {
           if(people[i].symptom == 1)  n_sym_idx++;
           else  n_asym_idx++;
        }
        else
        {
           if(people[i].symptom == 1)  n_sym_sec++;
           else  n_asym_sec++;
        }
     }
  }
  if(cfg_pars.silent_run == 0)  
  {
     printf("before right-censoring: n_imm=%d  n_esc=%d  n_sym=%d (idx=%d sec=%d)  n_asym=%d (idx=%d  sec=%d)\n", 
             n_imm, n_esc, n_sym_idx+n_sym_sec, n_sym_idx, n_sym_sec, n_asym_idx+n_asym_sec, n_asym_idx, n_asym_sec);
     printf("preimmune: %d/%d=%e\n", n_imm, p_size, (double) n_imm/p_size);
     printf("attack rate: %d/%d=%e\n", n_sym_sec+n_asym_sec, n_sym_sec+n_asym_sec+n_esc, (double) (n_sym_sec+n_asym_sec) / (n_sym_sec+n_asym_sec+n_esc));
     printf("pathogenicity: %d/%d=%e\n", n_sym_sec, (n_sym_sec+n_asym_sec), (double) n_sym_sec/(n_sym_sec+n_asym_sec));
  }
  /*
  n_exp_c = n_inf_c = n_exp_a = n_inf_a = 0;
  for(i=0; i < p_size; i++)
  if(people[i].pre_immune == 0 && people[i].idx == 0)
  {
     if(people[i].time_ind_covariate[0] == 0)
     {
        n_exp_c++;
        if(people[i].infection == 1)  n_inf_c++;
     }
     else
     {
        n_exp_a++;
        if(people[i].infection == 1)  n_inf_a++;
     }
  }
  printf("children: %d/%d=%e\n", n_inf_c, n_exp_c, (double) n_inf_c/n_exp_c);
  printf("adult: %d/%d=%e\n", n_inf_a, n_exp_a, (double) n_inf_a/n_exp_a);
  printf("Crude risk ratio: %f\n", (double) n_inf_a/n_exp_a / ((double) n_inf_c/n_exp_c));  
  */
  //Some cases (symptomatic or asymptomatic) are infected right before the last day but infectiousness onsets extend over the observation period.
  //We set these people to escapes, and this right-censoring of symptoms can be adjusted for in analysis.
  if(goodness_of_fit >= 0)
  {
     count = 0;
     for(i=0; i < p_size; i++)
     {
        h = people[i].community;
        if(people[i].infection == 1 && people[i].day_ill > community[h].day_epi_stop)
        {
           people[i].infection = 0;
           people[i].day_ill = MISSING;
           people[i].symptom = 0;
        }  
        if(people[i].day_ill != MISSING) count++;
     }
  }
  n_sym_idx = n_asym_idx = n_sym_sec = n_asym_sec = n_esc = n_imm = 0;
  for(i=0; i < p_size; i++)
  if(people[i].ignore == 0)
  {
     if(people[i].pre_immune == 1)  n_imm++;
     else if(people[i].infection == 0)  n_esc++;
     else
     {
        if(people[i].idx == 1)
        {
           if(people[i].symptom == 1)  n_sym_idx++;
           else  n_asym_idx++;
        }
        else
        {
           if(people[i].symptom == 1)  n_sym_sec++;
           else  n_asym_sec++;
        }
     }
  }
  if(cfg_pars.silent_run == 0)  
  {
     printf("After right-censoring: n_imm=%d  n_esc=%d  n_sym=%d (idx=%d sec=%d)  n_asym=%d (idx=%d  sec=%d)\n", 
             n_imm, n_esc, n_sym_idx+n_sym_sec, n_sym_idx, n_sym_sec, n_asym_idx+n_asym_sec, n_asym_idx, n_asym_sec);
     printf("preimmune: %d/%d=%e\n", n_imm, p_size, (double) n_imm/p_size);
     printf("Secondary attack rate: %d/%d=%e\n", n_sym_sec+n_asym_sec, n_sym_sec+n_asym_sec+n_esc, (double) (n_sym_sec+n_asym_sec) / (n_sym_sec+n_asym_sec+n_esc));
     printf("pathogenicity: %d/%d=%e\n", n_sym_sec, (n_sym_sec+n_asym_sec), (double) n_sym_sec/(n_sym_sec+n_asym_sec));
  }

  free(c2p_covariate); 
  free(p2p_covariate);
  return(0);
}




