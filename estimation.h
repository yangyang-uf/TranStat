////// Author: Yang Yang
////// Department of Biostatistics and Emerging Pathogens Institute
////// University of Florida
/*
void brutal_set_state(PEOPLE *person, int new_state)
{
   int h, i, j, k, l, r, t, old_state;
   double infective_prob;
   CONTACT * ptr_contact;
   PEOPLE *subject, *member;
   RISK *ptr_risk, *ptr2_risk;

   h = person->community;
   old_state = person->current_state;
   person->current_state = new_state;
   l = person->current_state;

   person->day_infection_lower = person->day_infection_upper = MISSING;
   person->day_infective_lower = person->day_infective_upper = MISSING;
   //printf("new state: %d  %d\n", person->possible_state[0][l], person->possible_state[1][l]);
   if(person->possible_state[0][l] == -INFINITY) //pre-immune
   {
      person->pre_immune = 1;
      person->infection = 0;
      person->symptom = 0;
      person->day_ill = MISSING;
   }
   else if(person->possible_state[0][l] == INFINITY) //escape
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
      person->symptom = person->possible_state[1][l];
      person->day_ill = person->possible_state[0][l];

      person->day_infection_lower = max(person->day_ill - cfg_pars.max_incubation, community[h].day_epi_start);
      person->day_infection_upper  = person->day_ill - cfg_pars.min_incubation;
      if(person->day_ill + cfg_pars.lower_infectious <= community[h].day_epi_stop)
      {
         person->day_infective_lower = person->day_ill + cfg_pars.lower_infectious;
         person->day_infective_upper = min(person->day_ill + cfg_pars.upper_infectious, community[h].day_epi_stop);
      }
   }   

   k = old_state;
   if((person->possible_state[0][l] > -INFINITY && person->possible_state[0][l] < INFINITY) ||
      (person->possible_state[0][k] > -INFINITY && person->possible_state[0][k] < INFINITY))
   {
      // p2p risk history
      for(j=0; j<community[h].size; j++)
      {
         i = community[h].member[j];
         subject = people + i;
         for(r=0; r<community[h].epi_duration; r++)
         {
            ptr_risk = subject->risk_history[r].p2p_risk;
            while(ptr_risk != NULL)
            {
               ptr2_risk = ptr_risk->next;
               free(ptr_risk->covariate);
               free(ptr_risk);
               ptr_risk = ptr2_risk;
            }
            subject->risk_history[r].p2p_risk = NULL;
         }
      }
      //rebuild p2p risk history
      for(j=0; j<community[h].size; j++)
      {
         i = community[h].member[j];
         subject = people + i;
         if(subject->infection == 1 &&
            subject->day_infective_lower != MISSING && subject->day_infective_upper != MISSING)
         {
            for(t=subject->day_infective_lower; t<=subject->day_infective_upper; t++)
            if(t >= community[h].day_epi_start && t <= community[h].day_epi_stop)
            {
               r = t - community[h].day_epi_start;
               infective_prob = cfg_pars.prob_infectious[t - subject->day_infective_lower];
               ptr_contact = subject->contact_history[r].p2p_contact;
               while(ptr_contact != NULL)
               {
                  member = people + ptr_contact->contact_id;
                  if(!(member->final_risk_day != MISSING && t > member->final_risk_day))
                     add_p2p_risk_history(t, member, ptr_contact, infective_prob, subject->symptom);
                  ptr_contact = ptr_contact->next;
               }
            }
         }
      }
   }
}
*/
void set_state(PEOPLE *person, int new_state)
{
   int h, i, j, k, l, r, t;
   double infective_prob;
   CONTACT *ptr_contact, *inf_ptr_contact, *sus_ptr_contact;
   PEOPLE *member;

   h = person->community;
   l = person->current_state;
   //if current state is asymptomatic infection, remove the risk imposed by this person
   //to other community members
   //if(person->id >= 0)  printf("commuity=%d  person=%d, old state: %d  %d  %d\n", h, person->id, l, person->possible_states[0][l], person->possible_states[1][l]);
   if(person->possible_states[0][l] > -INFINITY_INTEGER && person->possible_states[0][l] < INFINITY_INTEGER)
   {
      //printf("day_infective_lower=%d  day_infective_upper=%d\n", person->day_infective_lower, person->day_infective_upper);
      for(t=person->day_infective_lower; t<=person->day_infective_upper; t++)
      if(t >= community[h].day_epi_start && t <= community[h].day_epi_stop)
      {
         r = t - community[h].day_epi_start;
         infective_prob = cfg_pars.prob_infectious[t - person->day_infective_lower];
         if(cfg_pars.common_contact_history_within_community == 1)
         {
         	 //printf("Remove old risk contribution\n");
             delete_p2p_risk_history_in_risk_class(t, person, infective_prob, person->symptom);
         }
         else
         {   
            inf_ptr_contact = person->contact_history[r].p2p_contact;
            //there could be multiple p2p contact types with the same person.  
            //Do not jump out of searching after one deletion.
            while(inf_ptr_contact != NULL)
            {
               member = people + inf_ptr_contact->contact_id;
               //if(DEBUG == 1) printf("t=%d  contact=%d  final_risk_day=%d\n", t, member->id, member->final_risk_day);
               if(!(member->final_risk_day != MISSING && t > member->final_risk_day))
               {
                  //if(DEBUG == 1)  printf("Before deletion\n");
                  //if(DEBUG == 1)  show_p2p_risk(member, t, t);
                  sus_ptr_contact = inf_ptr_contact->pair;
                  delete_p2p_risk_history(t, member, sus_ptr_contact, infective_prob, person->symptom);
                  //if(DEBUG == 1)  printf("After deletion\n");
                  //if(DEBUG == 1)  show_p2p_risk(member, t, t);
               }
               inf_ptr_contact = inf_ptr_contact->next;
            }
         }
      }   
   }
   person->current_state = new_state;
   l = person->current_state;
   person->day_infection_lower = person->day_infection_upper = MISSING;
   person->day_infective_lower = person->day_infective_upper = MISSING;
   //if(person->id >= 0)  printf("new state: %d  %d  %d\n", l, person->possible_states[0][l], person->possible_states[1][l]);
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
      
      //for asymptomatic, we still use day_ill to denote the infectiousness onset day.
      person->day_ill = person->possible_states[0][l];
      if(person->day_ill + cfg_pars.lower_infectious <= community[h].day_epi_stop)
      {
         person->day_infective_lower = person->day_ill + cfg_pars.lower_infectious;
         person->day_infective_upper = min(person->day_ill + cfg_pars.upper_infectious, community[h].day_epi_stop);
         //printf("day_infective_lower=%d  day_infective_upper=%d\n", 
         //        person->day_infective_lower, person->day_infective_upper);
         for(t=person->day_infective_lower; t<=person->day_infective_upper; t++)
         if(t >= community[h].day_epi_start && t <= community[h].day_epi_stop)
         {
            r = t - community[h].day_epi_start;
            infective_prob = cfg_pars.prob_infectious[t - person->day_infective_lower];
            if(cfg_pars.common_contact_history_within_community == 1)
            {
               //printf("Add new risk contribution for day t=%d\n", t);
               ptr_contact = community[h].contact_history[r].p2p_contact;
               while(ptr_contact != NULL)
               {
                  if(ptr_contact->contact_id == person->id)
                     add_p2p_risk_history_to_risk_class(t, person, ptr_contact->contact_mode, ptr_contact->offset, infective_prob, person->symptom);
                  ptr_contact = ptr_contact->next;
               }
            }
            else
            {   
               inf_ptr_contact = person->contact_history[r].p2p_contact;
               while(inf_ptr_contact != NULL)
               {
                  member = people + inf_ptr_contact->contact_id;
                  //printf("t=%d  contact=%d  final_risk_day=%d\n", t, member->id, member->final_risk_day);
                  if(!(member->final_risk_day != MISSING && t > member->final_risk_day))
                  {
                     //printf("Before deletion\n");
                     //show_p2p_risk(member, t, t);
                     sus_ptr_contact = inf_ptr_contact->pair;
                     add_p2p_risk_history(t, member, person, sus_ptr_contact->contact_mode, sus_ptr_contact->offset, infective_prob, person->symptom);
                     //printf("After deletion\n");
                     //show_p2p_risk(member, t, t);
                  }
                  inf_ptr_contact = inf_ptr_contact->next;
               }
            }
         }
      }
      person->day_infection_lower = max(person->day_ill - cfg_pars.max_incubation, community[h].day_epi_start);
      person->day_infection_upper  = person->day_ill - cfg_pars.min_incubation;
      //printf("day_infection_lower=%d  day_infection_upper=%d\n", 
      //        person->day_infection_lower, person->day_infection_upper);
   }
}


// draw MCMC samples of infection states for a given community.
// pointer "head" gives the chain in which the samples will be stored,
// it does not have to be the chain in the pointer "community".
// The reason for such design is the need for new MCMC samples for all communities
// in the variance calculation, in particular if MC error is to be accounted for.
// to evaluate MC error, we repeat the estimation process from new sets of importance samples obtained
// via bootstrap or new MCMC sampling. When calculating variance for each new process,
// it's not necessary to generate MCMC samples for all communities for each process.
// They can share the same MCMC samples for variance calculation. 

void generate_sample_states(COMMUNITY *community, int n_new_sample, int n_burnin_sample, 
                           double *par_effective, STATE *head) 
{
   int i, j, k, l, m, n, max_n_states;
   int n_sample;
   double *log_prob, log_L_ini, temp;
   PEOPLE *person;
   STATE *ptr_stat;

   n_sample = n_new_sample + n_burnin_sample;
   max_n_states = 2 * max_epi_duration + 2;
   make_1d_array_double(&log_prob, max_n_states, 0.0);
   ptr_stat = head; 
   n = 0;
   while(n < n_sample & ptr_stat != NULL)
   {
      for(m=0; m<community->size_impute; m++)
      {
         person = people + community->member_impute[m];
         //go through possible values
         // first set all probabilities to be close to 0
         reset_1d_array_double(log_prob, max_n_states, -1e200);
         for(l=0; l<person->size_possible_states; l++)
         {
            set_state(person, l);
            if(cfg_pars.adjust_for_left_truncation == 1)   
            {
               if(cfg_pars.preset_index == 0)  set_index_cases(community);
               if(!(community->size_idx == 0 || community->size_idx == community->size))
                  log_prob[l] = community_log_likelihood(community, par_effective);
               else log_prob[l] = -1e200;
            }
            else  log_prob[l] = community_log_likelihood(community, par_effective);      
         }
         //sample an infectiousness onset day
         k = centered_discrete(person->size_possible_states, log_prob);
         if(person->current_state != k)
         {
            set_state(person, k);
            if(cfg_pars.adjust_for_left_truncation == 1 && cfg_pars.preset_index == 0)   set_index_cases(community);
         }   
         log_L_ini = log_prob[k];
      }
      if(n >= n_burnin_sample)
      {
         for(m=0; m<community->size_impute; m++)
         {
            person = people + community->member_impute[m];
            ptr_stat->states[m] = person->current_state;
         }
         ptr_stat->log_L_ini = log_L_ini;
         ptr_stat->log_L_cur = log_L_ini;
         ptr_stat->log_L_pre = log_L_ini;
         ptr_stat->ratio = 1.0;
         ptr_stat = ptr_stat->next;
      }
      n++;
   }
   free(log_prob);
}


void generate_list_states(COMMUNITY *community, double *par_effective)
{
   int i, j, k, l, m, n, count;
   int *width, *location;
   PEOPLE *person;
   STATE *ptr_stat, *ptr2_stat, *ptr3_stat;

   width = (int *) malloc((size_t) community->size_impute * sizeof(int));
   location = (int *) malloc((size_t) community->size_impute * sizeof(int));

   k = 1;
   for(m=community->size_impute-1; m>=0; m--)
   {
      person = people + community->member_impute[m];
      k *= person->size_possible_states;
      width[m] = k;
   }
  
   ptr_stat = community->list_states; 
   for(count=0; count < community->size_possible_states; count++)
   {
      //printf("Community possible state %d:\n", count);
      //printf("===========================================================\n");
      for(m=0; m<community->size_impute; m++)
      {
         n = count % width[m];
         if( m < community->size_impute - 1)
         {
            location[m] = n / width[m+1];
         }
         else  location[m] = n;
      }
      for(m=0; m<community->size_impute; m++)
      {
         person = people + community->member_impute[m];
         //printf("Set person %d's state as %d\n", person->id, location[m]);
         //printf("---------------------------------------------------------\n");
         set_state(person, location[m]);
         //printf("---------------------------------------------------------\n");
         ptr_stat->states[m] = location[m];
      }
      if(cfg_pars.adjust_for_left_truncation == 1)   
      {
         if(cfg_pars.preset_index == 0)  set_index_cases(community);
         //for case-ascertained design, we need to check whether the possible state leads to
         // valid index case; if not, these states 
         if(!(community->size_idx == 0 || community->size_idx == community->size))
         {
            ptr_stat->log_L_ini = community_log_likelihood(community, par_effective);
            ptr_stat->log_L_cur = ptr_stat->log_L_ini;
            ptr_stat->log_L_pre = ptr_stat->log_L_ini;
            ptr2_stat = ptr_stat;
            ptr_stat = ptr_stat->next;
         }   
      }
      else  
      {
         ptr_stat->log_L_ini = community_log_likelihood(community, par_effective);
         ptr_stat->log_L_cur = ptr_stat->log_L_ini;
         ptr_stat->log_L_pre = ptr_stat->log_L_ini;
         ptr2_stat = ptr_stat;
         ptr_stat = ptr_stat->next;
      }   
   }
   if(count < community->size_possible_states && ptr_stat != NULL)
   {
      ptr3_stat = community->list_states_rear;
      community->list_states_rear = ptr2_stat;
      community->list_states_rear->next = NULL;
      free_state_chain(ptr_stat, ptr3_stat);
      community->size_possible_states = count;
   }   
   free(width);
   free(location);
   //exit(0);
}

double likelihood(double *par_effective, double *extra_par)
{

  int i, j, h, k, l, m, n;
  double log_L, log_L_grp, weight, sum_weight, weight_offset;
  STATE *ptr_stat;
  PEOPLE *person;
  
  log_L = 0.0;
  for(h=0; h<n_community; h++)
  if(community[h].size > 0)
  {
     if(community[h].size_impute == 0)
     {
        log_L += community_log_likelihood(community + h, par_effective);
     }
     else if(community[h].size_possible_states > 1 && community[h].size_possible_states < cfg_pars.min_size_MCEM)
     {
        log_L_grp = 0.0;
        sum_weight = 0.0;
        ptr_stat = community[h].list_states;
        weight_offset = ptr_stat->log_L_pre;
        while(ptr_stat != NULL)
        {
           for(i=0; i<community[h].size_impute; i++)
           {
              person = people + community[h].member_impute[i];
              set_state(person, ptr_stat->states[i]);
           }  
           if(cfg_pars.adjust_for_left_truncation == 1 && cfg_pars.preset_index == 0)   set_index_cases(community + h); 
           ptr_stat->log_L_cur = community_log_likelihood(community + h, par_effective);
           weight = exp(ptr_stat->log_L_pre - weight_offset);
           log_L_grp += weight * ptr_stat->log_L_cur;
           sum_weight += weight;
           ptr_stat = ptr_stat->next;
        }
        log_L += log_L_grp / sum_weight;
     }
     else if(community[h].size_possible_states >= cfg_pars.min_size_MCEM)
     {
        log_L_grp = 0.0;
        sum_weight = 0.0;
        ptr_stat = community[h].sample_states;
        k = 0;
        while(ptr_stat != NULL)
        {
           // go over members with uncertain state, and assign sample state
           for(i=0; i<community[h].size_impute; i++)
           {
              person = people + community[h].member_impute[i];
              set_state(person, ptr_stat->states[i]);
           }   
           if(cfg_pars.adjust_for_left_truncation == 1 && cfg_pars.preset_index == 0)   set_index_cases(community + h); 
           ptr_stat->log_L_cur = community_log_likelihood(community + h, par_effective);
           //if(cfg_pars.community_specific_weighting == 1)  
           //   weight = ptr_stat->ratio;// do not use log_L_cur. Its value changes during the maximization.
           //else  weight = importance_weight[k];
           weight = importance_weight[k];
           log_L_grp += weight * ptr_stat->log_L_cur;
    
           sum_weight += weight;
           ptr_stat = ptr_stat->next;
           k++;
        }
        log_L += log_L_grp / sum_weight;
     }
  }

  if(!(log_L <= 0))  log_L = -1e200;
  return(log_L);
}

double likelihood2(double *par_effective, int hid1, int hid2)
{

  int i, j, h, k, l, m, n;
  double log_L, log_L_grp, weight, sum_weight;
  STATE *ptr_stat;
  PEOPLE *person;
  
  log_L = 0.0;
  for(h=hid1; h<=hid2; h++)
  if(community[h].size > 0)
  {
     log_L += community_log_likelihood(community + h, par_effective);
     printf("HH %d: log_L=%e cum_log_L=%e\n", h - hid1, community_log_likelihood(community + h, par_effective), log_L);
  }
  if(!(log_L <= 0))  log_L = -1e200;
  return(log_L);
}

void update_weight(double *par_effective)
{

  int i, j, h, k, l, m, n;
  double *log_L_cur, *log_L_ini;
  STATE *ptr_stat;
  PEOPLE *person;
  

  make_1d_array_double(&log_L_cur, size_sample_states, 0.0);
  make_1d_array_double(&log_L_ini, size_sample_states, 0.0);
  if(importance_weight != NULL)  free(importance_weight);
  importance_weight = (double *) malloc((size_t) size_sample_states * sizeof(double));
  for(h=0; h<n_community; h++)
  {
     // we update both log_L_ini and log_L_cur for communities with possible state >1 but < min_size_MCEM.
     // These communities need EM but not MCEM. log_L_ini is used as the weights for the EM and should be
     // updated before each maximization step.
     
     if(community[h].size_possible_states > 1 && 
        community[h].size_possible_states < cfg_pars.min_size_MCEM)
     {
        ptr_stat = community[h].list_states;
        while(ptr_stat != NULL)
        {
           // go over members with uncertain state, and assign sample state
           if(par_effective != NULL)
           {
              for(i=0; i<community[h].size_impute; i++)
              {
                 person = people + community[h].member_impute[i];
                 set_state(person, ptr_stat->states[i]);
              }   
              if(cfg_pars.adjust_for_left_truncation == 1 && cfg_pars.preset_index == 0)   set_index_cases(community + h); 
              ptr_stat->log_L_cur = community_log_likelihood(community + h, par_effective);
              ptr_stat->log_L_pre = ptr_stat->log_L_cur; //store the likelihood based on parameter estimates from previous iteration,
                                                         // so that it can be used for the next iteration.
           }
           ptr_stat = ptr_stat->next;
        }
     }

     // update log_L_cur for communities with possible state > min_size_MCEM.
     // Do not change log_L_ini, as they are fixed at the parameters based on which
     // importance samples are obtained.
     if(community[h].size_possible_states >= cfg_pars.min_size_MCEM)
     {
        k = 0;
        ptr_stat = community[h].sample_states;
        while(ptr_stat != NULL)
        {
           // go over members with uncertain state, and assign sample state
           if(par_effective != NULL)
           {
              for(i=0; i<community[h].size_impute; i++)
              {
                 person = people + community[h].member_impute[i];
                 set_state(person, ptr_stat->states[i]);
              }   
              if(cfg_pars.adjust_for_left_truncation == 1 && cfg_pars.preset_index == 0)   set_index_cases(community + h); 
              ptr_stat->log_L_cur = community_log_likelihood(community + h, par_effective);
              ptr_stat->log_L_pre = ptr_stat->log_L_cur; //store the likelihood based on parameter estimates from previous iteration,
                                                         // so that it can be used for the next iteration.
              ptr_stat->ratio = exp(ptr_stat->log_L_pre - ptr_stat->log_L_ini);
           }
           log_L_cur[k] += ptr_stat->log_L_cur;
           log_L_ini[k] += ptr_stat->log_L_ini;
           ptr_stat = ptr_stat->next;
           k++;
        }
     }
     // update weights
     for(k=0; k<size_sample_states; k++)
        importance_weight[k] = exp(log_L_cur[k] - log_L_ini[k]);
  }

  free(log_L_cur);
  free(log_L_ini);
}
// calculate importance weights for importance sampling.
// If importance samples have been augmented in this iteration,
// we have log-likelihoods for new samples at current parameter values, 
// and need to calculate:
// (1) log-likelihoods for old samples at current parameter values
// (2) log-likelihoods for new samples at initial parameter values
// If there is no augmentation, we only need to calculate (1).
int E_step_importance(double *par_effective_cur, double *log_L, MATRIX *score, MATRIX *info)
{
   int i, j, h, k, l, m, n, count;
   int n_par; 
   int error = 0;
   double log_L_grp, weight, sum_weight, weight_offset;
   MATRIX first, second, score_grp, info_grp;
   PEOPLE *person;
   STATE *ptr_stat;
   FILE *file;
   
   n_par = cfg_pars.n_par;

   initialize_matrix(&first);
   initialize_matrix(&second);
   initialize_matrix(&score_grp);
   initialize_matrix(&info_grp);

   inflate_matrix(&first, n_par, 1, 0);
   inflate_matrix(&second, n_par, n_par, 0);
   inflate_matrix(&info_grp, n_par, n_par, 0);
   inflate_matrix(&score_grp, n_par, 1, 0);

   (* log_L) = 0.0;
   reset(score, 0.0);
   reset(info, 0.0);
   //file = fopen("check_derivatives.dat", "w");
   for(h=0; h<n_community; h++)
   if(community[h].size > 0 && community[h].ignore == 0)
   {
      //printf("%d: size_impute=%d  size_possible_states=%d  min_size_MCEM=%d\n", h, community[h].size_impute, community[h].size_possible_states, cfg_pars.min_size_MCEM);
      if(community[h].size_impute == 0)
      {
         error = community_derivatives(community+h, par_effective_cur, &log_L_grp, &score_grp, &info_grp);
         if(error == 1)  goto end;
         (* log_L) += log_L_grp;
         addition(score, 1.0, &score_grp, 1.0, score);
         addition(info, 1.0, &info_grp, -1.0, info);
		   //fprintf(file, "%d: %e  %e  %e\n", h, score_grp.data[0][0], score_grp.data[1][0], score_grp.data[2][0]);
		   //printf("Scores: %e  %e\n", score_grp.data[0][0], score_grp.data[1][0]);
      }
      // when size_impute>1 and <min_size.MCEM, we apply the EM algorithm, i.e., average derivatives over all possible state
      // with the likelihood based on parameter estimates from last EM iteration  as the weight.
      else if(community[h].size_possible_states > 1 && community[h].size_possible_states < cfg_pars.min_size_MCEM)
      {
         log_L_grp = 0.0;
         reset(&score_grp, 0.0);
         reset(&info_grp, 0.0);
         ptr_stat = community[h].list_states;
         weight_offset = ptr_stat->log_L_pre; // the weight is the likelihood, which may bt too close to 0 and will lead to numerical problem. We use an offset to solve this problem
         sum_weight = 0.0;
         count = 0;
         while(ptr_stat != NULL)
         {
            for(i=0; i<community[h].size_impute; i++)
            {
               person = people + community[h].member_impute[i];
               set_state(person, ptr_stat->states[i]);
            }   
            if(cfg_pars.adjust_for_left_truncation == 1 && cfg_pars.preset_index == 0)   set_index_cases(community + h); 
            error = community_derivatives(community+h, par_effective_cur, &ptr_stat->log_L_cur, &first, &second);
            if(error == 0)
            {
               weight = exp(ptr_stat->log_L_pre - weight_offset);// do not use log_L_cur. Its value changes during the maximization.
               sum_weight += weight;
               log_L_grp += weight * ptr_stat->log_L_cur;
               addition(&score_grp, 1.0, &first, weight, &score_grp);
               addition(&info_grp, 1.0, &second, -weight, &info_grp);
               //printf("State %d:  weight=%e  sum_weight=%e  log_L_pre=%e  log_L_cur=%e  log_L_grp=%e\n", count, weight, sum_weight, ptr_stat->log_L_pre, ptr_stat->log_L_cur, log_L_grp);
            }
            ptr_stat = ptr_stat->next;
            count++;
         }
            
         (* log_L) += log_L_grp / sum_weight;
         addition(score, 1.0, &score_grp, 1.0/sum_weight, score);
         addition(info, 1.0, &info_grp, 1.0/sum_weight, info);

      }
      // when size_impute>1, it is impossible to enumerate all possible state. 
      // We average derivatives over sampled state with weights just like in importance sampling.
      else if(community[h].size_possible_states >= cfg_pars.min_size_MCEM)
      {
         log_L_grp = 0.0;
         reset(&score_grp, 0.0);
         reset(&info_grp, 0.0);
         sum_weight = 0.0;
         ptr_stat = community[h].sample_states;
         k = 0;
         while(ptr_stat != NULL)
         {
            // go over members with uncertain state, and assign sample state
            for(i=0; i<community[h].size_impute; i++)
            {
               //if(h == 8)  printf("%d\n", ptr_stat->state[i]);
               person = people + community[h].member_impute[i];
               set_state(person, ptr_stat->states[i]);
               //if(h == 8)
               //   printf("%d: pre_immune=%d  infection=%d  symptom=%d  day_ill=%d  n_state=%d\n", 
               //           person->id, person->pre_immune, person->infection, person->symptom, 
               //           person->day_ill, person->size_possible_state);
            }   
            //if(h == 8)  printf("\n");
            if(cfg_pars.adjust_for_left_truncation == 1 && cfg_pars.preset_index == 0)   set_index_cases(community + h); 
            error = community_derivatives(community+h, par_effective_cur, &ptr_stat->log_L_cur, &first, &second);
            //if(h == 8)  printf("%d  error=%d\n", k++, error);
            if(error == 0)
            {
               //if(cfg_pars.community_specific_weighting == 1)  
               //   weight = ptr_stat->ratio;// do not use log_L_cur. Its value changes during the maximization.
               //else  weight = importance_weight[k];
               weight = importance_weight[k];
               sum_weight += weight;
               log_L_grp += weight * ptr_stat->log_L_cur;
               addition(&score_grp, 1.0, &first, weight, &score_grp);
               addition(&info_grp, 1.0, &second, -weight, &info_grp);
            }
            ptr_stat = ptr_stat->next;
            k++;
         }
         (* log_L) += log_L_grp / sum_weight;
         addition(score, 1.0, &score_grp, 1.0/sum_weight, score);
         addition(info, 1.0, &info_grp, 1.0/sum_weight, info);
      }
   }
   //fclose(file);
   //exit(0);
   if(!((* log_L) <= 0))  (* log_L) = -1e200;
   
end:
   deflate_matrix(&first);
   deflate_matrix(&second);
   deflate_matrix(&score_grp);
   deflate_matrix(&info_grp);

   return(error);
}


int M_step(double *par_effective, double *max_log_L, int optimization_method)
{
   int i, j, h, k, l, m, n, loop;
   int error = 0;
   int max_iter = 500;
   int n_par, n_par_equiclass;
   int converge, positive, n_iter, num_func_eval, stop;
   int n_b_mode, n_p_mode, n_q_mode, n_u_mode;

   double log_L, log_L_grp, log_L_new, log_L_try, ftol, gap;
   double *step;
   double *old_par_effective, *par_effective_new, *par_effective_try;

   MATRIX info, inv, diff, score;
   MATRIX sub_info, sub_inv, sub_diff, sub_score;
   
   n_b_mode = cfg_pars.n_b_mode;
   n_p_mode = cfg_pars.n_p_mode;
   n_u_mode = cfg_pars.n_u_mode;
   n_q_mode = cfg_pars.n_q_mode;

   old_par_effective = par_effective_new = par_effective_try = NULL;

   initialize_matrix(&score);
   initialize_matrix(&info);
   initialize_matrix(&inv);
   initialize_matrix(&diff);

   initialize_matrix(&sub_score);
   initialize_matrix(&sub_info);
   initialize_matrix(&sub_inv);
   initialize_matrix(&sub_diff);

   n_par = cfg_pars.n_par;
   n_par_equiclass = cfg_pars.n_par_equiclass;
   
   inflate_matrix(&info, n_par, n_par, 0);
   inflate_matrix(&inv, n_par, n_par, 0);
   inflate_matrix(&score, n_par, 1, 0);
   inflate_matrix(&diff,n_par,1,0);
  
   inflate_matrix(&sub_info, n_par_equiclass, n_par_equiclass, 0);
   inflate_matrix(&sub_inv, n_par_equiclass, n_par_equiclass, 0);
   inflate_matrix(&sub_score, n_par_equiclass, 1, 0);
   inflate_matrix(&sub_diff,n_par_equiclass,1,0);
  
   make_1d_array_double(&old_par_effective, n_par_equiclass, 0.0);
   make_1d_array_double(&par_effective_new, n_par_equiclass, 0.0);
   make_1d_array_double(&par_effective_try, n_par_equiclass, 0.0);
   make_1d_array_double(&step, n_par_equiclass, 0.0);

   // optimization_method=1 for Newton-Raphson
   // optimization_method=2 for Nelder-Mead
   if(optimization_method == 2)
   {
   
      for(j=0; j<n_par_equiclass; j++)  step[j] = 0.5;
      ftol = 1e-8;

      down_hill_simplex(par_effective, cfg_pars.lower_search_bound, cfg_pars.upper_search_bound,
         step, &log_L, n_par_equiclass, 0, ftol, likelihood, NULL, &num_func_eval);
      if(num_func_eval >= SIMPLEX_MAX)  error = 10;
   }
   else
   {
      loop = 0;
      converge = 0;
      while( converge == 0 && loop < max_iter)
      {
         //printf("\n M_step loop=%d\n", loop);
         loop ++;
         E_step_importance(par_effective, &log_L, &score, &info);
         //printf("log_L=%e\n", log_L);
         reset(&diff, 0.0);
         reset(&inv, 0.0);
         // if user specify any parameters should be equal,
         // we group the parameters into equi-classes.
         for(i=0; i<n_par_equiclass; i++)   
         {
            sub_score.data[i][0] = 0.0;
            for(j=i; j<n_par_equiclass; j++)
               sub_info.data[i][j] = sub_info.data[j][i] = 0.0;
         }

         for(i=0; i<n_par_equiclass; i++)
         {
            for(j=0; j<cfg_pars.par_equiclass[i].size; j++)
            {
               m = cfg_pars.par_equiclass[i].member[j] - 1;
               sub_score.data[i][0] += score.data[m][0];
            }
         }

         for(i=0; i<n_par_equiclass; i++)
         {
            for(j=i; j<n_par_equiclass; j++)
            {
               for(k=0; k<cfg_pars.par_equiclass[i].size; k++)
               {
                   
                  m = cfg_pars.par_equiclass[i].member[k] - 1;
                  for(l=0; l<cfg_pars.par_equiclass[j].size; l++)
                  {
                     n = cfg_pars.par_equiclass[j].member[l] - 1;
                     sub_info.data[i][j] += info.data[m][n];
                  }
               }
               sub_info.data[j][i] = sub_info.data[i][j];
            }
         }


         if(inverse(&sub_info,&sub_inv) > 0)
         {
            printf("information matrix is not invertable\n");
            fmprintf(&sub_score);
            fmprintf(&sub_info);
            exit(0);
            error = 1;
            goto end;
         }
         product(&sub_inv,&sub_score,&sub_diff);

         //parse the changes to original parameters
         for(i=0; i<n_par_equiclass; i++)
         {
            for(j=0; j<cfg_pars.par_equiclass[i].size; j++)
            {
               m = cfg_pars.par_equiclass[i].member[j] - 1;
               diff.data[m][0] = sub_diff.data[i][0];
            }
         }

         //parse the Fisher information to original parameters
         for(i=0; i<n_par_equiclass; i++)
         {
            for(j=i; j<n_par_equiclass; j++)
            {
               for(k=0; k<cfg_pars.par_equiclass[i].size; k++)
               {

                  m = cfg_pars.par_equiclass[i].member[k] - 1;
                  for(l=0; l<cfg_pars.par_equiclass[j].size; l++)
                  {
                     n = cfg_pars.par_equiclass[j].member[l] - 1;
                     inv.data[m][n] = sub_inv.data[i][j];
                     inv.data[n][m] = inv.data[m][n];
                  }
               }
            }
         }

         for(i=0; i<n_par_equiclass; i++)  old_par_effective[i] = par_effective[i];
         //for(i=0; i<n_par_equiclass; i++)  par_effective[i] = par_effective[i] + sub_diff.data[i][0];
         //show_par_effective(par_effective, NULL);
         
         positive = 1;
         for(i=0; i<n_par_equiclass; i++)
         {
            if(sub_inv.data[i][i] < 0)  positive=0;
         }
         //printf("positive=%d\n", positive);
         if(positive == 0)
         {
            for(i=0; i<n_par_equiclass; i++)
               par_effective_new[i] = par_effective[i] + sub_score.data[i][0] / fabs(sub_info.data[i][i]);
         }
         else
         {
            for(i=0; i<n_par_equiclass; i++)  par_effective_new[i] = par_effective[i] + sub_diff.data[i][0];
         }
         
         log_L_new = likelihood(par_effective_new, NULL);
         //printf("Before further searching: log_L=%20.12lf  log_L_new=%20.12lf\n", log_L, log_L_new);
         //show_par_effective(par_effective);
         //show_par_effective(par_effective_new);
         if(log_L_new < log_L)
         {
            //printf("contractive searching:\n");
            stop = 0;
            n_iter = 0;
            while(stop == 0)
            {
               n_iter++;
               for(i=0; i<n_par_equiclass; i++)  par_effective_try[i] = 0.5 * (par_effective[i] + par_effective_new[i]);
               log_L_try = likelihood(par_effective_try, NULL);
               for(i=0; i<n_par_equiclass; i++) par_effective_new[i] = par_effective_try[i];
               log_L_new = log_L_try;
               if(log_L_try > log_L || n_iter > 20) stop = 1;
            }
            //If contractive searching does not identify a better point, perform a contractive searching
            //in the opposite direction.
            
            if(log_L_new <= log_L)
            {
               //printf("contractive searching in opposite direction: \n");
               for(i=0; i<n_par_equiclass; i++)  par_effective_new[i] = par_effective[i] - sub_diff.data[i][0];
               stop = 0;
               n_iter = 0;
               while(stop == 0)
               {
                  n_iter++;
                  for(i=0; i<n_par_equiclass; i++)  par_effective_try[i] = 0.5 * (par_effective[i] + par_effective_new[i]);
                  log_L_try = likelihood(par_effective_try, NULL);
                  for(i=0; i<n_par_equiclass; i++) par_effective_new[i] = par_effective_try[i];
                  log_L_new = log_L_try;
                  if(log_L_try > log_L || n_iter > 20) stop = 1;
               }
               //if(log_L_try > log_L_all)  printf("succeeded\n");
               //else  printf("failed\n");
            }  
         }
         else if(positive == 0)
         {
            stop = 0;
            n_iter = 0;
            while(stop == 0)
            {
               n_iter++;
               for(i=0; i<n_par_equiclass; i++)
                  par_effective_try[i] = par_effective_new[i] + (par_effective_new[i] - par_effective[i]);
               log_L_try = likelihood(par_effective_try, NULL);
               if(log_L_try < log_L_new || n_iter > 20)  stop = 1;
               else
               {
                  for(i=0; i<n_par_equiclass; i++) par_effective_new[i] = par_effective_try[i];
                  log_L_new = log_L_try;
               }
            }
         }
         //printf("after further searching: log_L=%20.12lf  log_L_new=%20.12lf\n", log_L, log_L_new);
         //show_par_effective(par_effective_new, NULL);
         //fmprintf(&info);
         //printf("\n\n");
         //exit(0);
         if(log_L_new < log_L)
         {
            //printf("Can not find any point better than the current\n");
            break;
         }
         else
         {
            //if successfully found a better point, update par_effective and par
            for(i=0; i<n_par_equiclass; i++)
               par_effective[i] = par_effective_new[i];
            log_L = log_L_new;
         }
         
         
         converge = 1;
         for(i=0; i<n_par_equiclass; i++)
         {
            m = cfg_pars.par_equiclass[i].member[0] - 1;
            if(m < n_b_mode + n_p_mode + n_u_mode + n_q_mode)
               gap = fabs(inv_logit(old_par_effective[i]) - inv_logit(par_effective[i])) / fabs(inv_logit(old_par_effective[i]));
            else
               gap = fabs(exp(old_par_effective[i]) - exp(par_effective[i])) / fabs(exp(old_par_effective[i]));
            if(gap > cfg_pars.converge_criteria[i])
               converge = 0;
         }
      }   
      //exit(0);
      //if(loop >= max_iter)
      //{
      //   error = 1;
      //   printf("max number of iteration exceeded for the M-step=%d\n", loop);
      //   goto end;
      //}
   }
   
   if(cfg_pars.silent_run == 0)  printf("number of M-step iterations=%d\n", loop); 
 end:
   (* max_log_L) = log_L;
   deflate_matrix(&diff);
   deflate_matrix(&info);
   deflate_matrix(&score);
   deflate_matrix(&inv);
   deflate_matrix(&sub_diff);
   deflate_matrix(&sub_info);
   deflate_matrix(&sub_score);
   deflate_matrix(&sub_inv);
   free(old_par_effective);
   free(par_effective_new);
   free(par_effective_try);
   free(step);
   
   return(error);
}
/*
int select_sub_sample(int m, int *index)
{
   int count, i, j, k;
   double v, d, gamma;
   count = 0;
   k = 0;
   v = 1.0;
   d = 0.5;
   while(count < m)
   {
      gamma = v * pow((k+1), d);
      j = rpoisson(gamma, &seed) + 1;
      count += j;
      index[k] = count - 1;
      //printf("%d  %d  %d\n", k, j, index[k]);
      k++;
   }
   return(k);
}

int augmentation(double *ini_par_effective, double *cur_par_effective, int n_ind_sample, int *index, double *Q_stat, int Q_stat_only)
{
   int i, j, h, k, l, m, n;
   int i1, i2;
   int n_par, MC_swamp, n_new_sample;
   int error;
   double log_L, sum_weight, sub_sum_weight, sub_sum_weight_sq, factor, sd, lower, upper, average_weight;
   MATRIX first, score, score_by_sample, score2, sub_score;
   PEOPLE *person;
   STATE *ptr_stat;

   n_par = cfg_pars.n_par;

   initialize_matrix(&first);
   initialize_matrix(&score);
   initialize_matrix(&score2);
   initialize_matrix(&sub_score);
   initialize_matrix(&score_by_sample);

   inflate_matrix(&first, n_par, 1, 0);
   inflate_matrix(&score, n_par, 1, 0);
   inflate_matrix(&score2, n_par, 1, 0);
   inflate_matrix(&sub_score, n_par, 1, 0);
   inflate_matrix(&score_by_sample, n_par, size_sample_state, 0);


   reset(&score_by_sample, 0.0);
   for(h=0; h<n_community; h++)
   if(community[h].size > 0)
   {
      if(community[h].size_impute > 0)
      {
         // get the indices of a subset of samples that can be considered independent
         k = 0;
         ptr_stat = community[h].sample_state;
         while(ptr_stat != NULL)
         {
            // go over members with uncertain state, and assign sample state
            for(i=0; i<community[h].size_impute; i++)
            {
               person = people + community[h].member_impute[i];
               set_state(person, ptr_stat->state[i]);
            }   
            error = community_score(community+h, cur_par_effective, &log_L, &first);
            if(error == 0)
            {
               for(j=0; j<n_par; j++)  
                  score_by_sample.data[j][k] += first.data[j][0];

            }
            ptr_stat = ptr_stat->next;
            k++;
         }
      }
      else
      {
         error = community_score(community+h, cur_par_effective, &log_L, &first);
         if(error == 0)
         {
            for(k=0; k<size_sample_state; k++)  
               for(j=0; j<n_par; j++)  
                  score_by_sample.data[j][k] += first.data[j][0];
         }
      }
   }
   reset(&score, 0.0);
   reset(&score2, 0.0);
   reset(&sub_score, 0.0);
   sum_weight = sub_sum_weight = 0.0;
   j = 0;
   for(i=0; i<size_sample_state; i++)
   {
      sum_weight += importance_weight[i];
      for(k=0; k<n_par; k++)  
      {
         score.data[k][0] += importance_weight[i] * score_by_sample.data[k][i];
         score2.data[k][0] += importance_weight[i] * score_by_sample.data[k][i] * score_by_sample.data[k][i];
      }
      
      if(j < n_ind_sample)
      {
         if(i == index[j])
         {
            sub_sum_weight += importance_weight[i];
            for(k=0; k<n_par; k++)  sub_score.data[k][0] += importance_weight[i] * score_by_sample.data[k][i];
            j++;
         }
      }
   }
   for(i=0; i<n_par; i++)  
   {
      score.data[i][0] /= sum_weight;
      score2.data[i][0] /= sum_weight;
      sub_score.data[i][0] /= sub_sum_weight;
   }
   // no augmentation is done at iteration 0 as no Q statistics are available
   if(Q_stat_only == 0)
   {
      sub_sum_weight = sub_sum_weight_sq = 0.0;
      j = 0;
      for(i=0; i<size_sample_state; i++)
      if(j < n_ind_sample)
      {
         if(i == index[j])
         {
            sub_sum_weight += importance_weight[i];
            sub_sum_weight_sq += importance_weight[i] * importance_weight[i];
            j++;
         }
      }
      factor = sqrt(sub_sum_weight_sq) / sub_sum_weight; 
      printf("variance shrinking factor=%e\n", factor);
      MC_swamp = 0;
      for(k=0; k<n_par; k++)
      {
         sd = factor * sqrt(score2.data[k][0] - score.data[k][0] * score.data[k][0]);
         //1.15 is the 87.5% quantile og standard normal
         lower = score.data[k][0] - 1.15 * sd;
         upper = score.data[k][0] + 1.15 * sd;
         if(Q_stat[k] > lower && Q_stat[k] < upper)  MC_swamp = 1;
         printf("par %d: score=%e  sub_score=%e  lower=%e  upper=%e  Q=%e\n", 
                 k, score.data[k][0], sub_score.data[k][0], lower, upper, Q_stat[k]);
      }

      n_new_sample = (int) floor(size_sample_state / 3);
      // if the EM estimator is swamped by the MC errors, use gibbs sampler to
      // augment the set of sample state
      if(MC_swamp == 1)
      {
         for(h=0; h<n_community; h++)
            generate_mcmc_samples(community + h, ini_par_effective, cur_par_effective, n_new_sample, cfg_pars.n_burnin_sampling);
         size_sample_state += n_new_sample; 
      }
      printf("swamp=%d  n_new_sample=%d  size_sample_state=%d\n", MC_swamp, n_new_sample, size_sample_state);
   }

   // compute the Q statistic for MC error evaluation next iteration.
   for(k=0; k<n_par; k++) Q_stat[k] = sub_score.data[k][0]; 
   //for(k=0; k<n_par; k++) printf("%e  ", Q_stat[k]); 
   deflate_matrix(&first);
   deflate_matrix(&score);
   deflate_matrix(&score2);
   deflate_matrix(&sub_score);
   deflate_matrix(&score_by_sample);

   return(MC_swamp);
}
*/
/* this program is for estimating MLE for the model with
within-household person-to-person transmission probability p1
, across-household transmission probability p2, and common source probability of infection (CSPI) b. */


int covariance(double *par_effective, COMMUNITY *aux_community, MATRIX *var_logit, int set_equal_weights)
{
   int i, j, h, k, l, m, n, count;
   int n_par, n_par_equiclass, n_b_mode, n_p_mode, n_q_mode, n_u_mode; 
   int error = 0;
   double *der, log_L_grp, weight, sum_weight;
   double *log_L_cur, *log_L_ini;
   MATRIX first, second, score_grp, info_grp;
   MATRIX score, score_by_sample, score2, info;
   MATRIX em_info, sub_info, sub_inv;
   STATE *ptr_stat, *ptr2_stat;
   PEOPLE *person;

   n_b_mode = cfg_pars.n_b_mode;
   n_p_mode = cfg_pars.n_p_mode;
   n_u_mode = cfg_pars.n_u_mode;
   n_q_mode = cfg_pars.n_q_mode;
   n_par = cfg_pars.n_par;
   n_par_equiclass = cfg_pars.n_par_equiclass;

   initialize_matrix(&first);
   initialize_matrix(&second);
   initialize_matrix(&score_grp);
   initialize_matrix(&info_grp);
   initialize_matrix(&score);
   initialize_matrix(&score2);
   initialize_matrix(&score_by_sample);
   initialize_matrix(&info);
   initialize_matrix(&em_info);
   initialize_matrix(&sub_info);
   initialize_matrix(&sub_inv);

   inflate_matrix(&first, n_par, 1, 0);
   inflate_matrix(&second, n_par, n_par, 0);
   inflate_matrix(&score_grp, n_par, 1, 0);
   inflate_matrix(&info_grp, n_par, n_par, 0);
   inflate_matrix(&score, n_par, 1, 0);
   inflate_matrix(&score2, n_par, n_par, 0);
   inflate_matrix(&info, n_par, n_par, 0);
   inflate_matrix(&em_info, n_par, n_par, 0);
   inflate_matrix(&sub_info, n_par_equiclass, n_par_equiclass, 0);
   inflate_matrix(&sub_inv, n_par_equiclass, n_par_equiclass, 0);
    
   reset(&score, 0.0);
   reset(&score2, 0.0);
   reset(&info, 0.0);
   reset(var_logit, 0.0);
   if(size_sample_states > 0)  
   {
      inflate_matrix(&score_by_sample, n_par, size_sample_states, 0);
      reset(&score_by_sample, 0.0);
   }
   if(aux_community != NULL)
   {
      //communities  may have different number of possible state, which makes the expectation of the
      //the square of score function extrememly difficult to calculate. Our solution is to
      //regenerate importance samples via mcmc for everybody with > 1 possible state,
      //and all have equal number of imprtance samples. These samples are generated based on the final
      //parameter estimates, so all imprtance weightes are set to 1, as any further iterative update
      //shouldn't change the parameter estimates much.
      if(set_equal_weights == 1)
      {
         for(k=0; k<size_sample_states; k++)
            importance_weight[k] = 1.0;
      }
      else
      {
         make_1d_array_double(&log_L_cur, size_sample_states, 0.0);
         make_1d_array_double(&log_L_ini, size_sample_states, 0.0);
         for(h=0; h<n_community; h++)
         if(community[h].size_possible_states > 1)
         {
            k = 0;
            ptr_stat = aux_community[h].sample_states;
            while(ptr_stat != NULL)
            {
               // go over members with uncertain state, and assign sample state
               for(i=0; i<community[h].size_impute; i++)
               {
                  person = people + community[h].member_impute[i];
                  set_state(person, ptr_stat->states[i]);
               }   
               if(cfg_pars.adjust_for_left_truncation == 1 && cfg_pars.preset_index == 0)   set_index_cases(community + h);
               ptr_stat->log_L_cur = community_log_likelihood(community + h, par_effective);
               //printf("log_L_cur=%e  log_L_ini=%e\n", ptr_stat->log_L_cur, ptr_stat->log_L_ini);
               log_L_cur[k] += ptr_stat->log_L_cur;
               log_L_ini[k] += ptr_stat->log_L_ini;
               ptr_stat = ptr_stat->next;
               k++;
            }
         }
         // update weights
         for(k=0; k<size_sample_states; k++)
            importance_weight[k] = exp(log_L_cur[k] - log_L_ini[k]);
         free(log_L_cur);
         free(log_L_ini);
      }
   }
   
   for(h=0; h<n_community; h++)
   if(community[h].size > 0)
   {
      //printf("h=%d  size=%d  size_impute=%d\n", h, community[h].size, community[h].size_impute);  
      if(community[h].size_impute == 0)
      {
         error = community_derivatives(community+h, par_effective, &log_L_grp, &first, &second);
         if(error == 0)
         {
            addition(&info, 1.0, &second, -1.0, &info);
            if(cfg_pars.EM == 1)
            {
               for(k=0; k<size_sample_states; k++)  
                  for(j=0; j<n_par; j++)  
                     score_by_sample.data[j][k] += first.data[j][0];
            }
         }
      }
      else
      {
         reset(&info_grp, 0.0);
         ptr_stat = (aux_community == NULL)? community[h].sample_states : aux_community[h].sample_states;
         sum_weight = 0.0;
         k = 0;
         while(ptr_stat != NULL)
         {
            // go over members with uncertain state, and assign sample state
            for(i=0; i<community[h].size_impute; i++)
            {
               person = people + community[h].member_impute[i];
               set_state(person, ptr_stat->states[i]);
            }   
            if(cfg_pars.adjust_for_left_truncation == 1 && cfg_pars.preset_index == 0)   set_index_cases(community + h);
            error = community_derivatives(community+h, par_effective, &ptr_stat->log_L_cur, &first, &second);
            if(error == 0) 
            {
               sum_weight += importance_weight[k];
               addition(&info_grp, 1.0, &second, -importance_weight[k], &info_grp);
               for(j=0; j<n_par; j++)  
                  score_by_sample.data[j][k] += first.data[j][0];
            }
            ptr_stat = ptr_stat->next;
            k++;
         }
         addition(&info, 1.0, &info_grp, 1.0/sum_weight, &info);
      }
   }
   
   if(cfg_pars.EM == 1)
   {
      sum_weight = 0.0;
      for(k=0; k<size_sample_states; k++)  
      {
         for(i=0; i<n_par; i++)  
         {
            score.data[i][0] += importance_weight[k] * score_by_sample.data[i][k];
            for(j=0; j<n_par; j++)  
               score2.data[i][j] += importance_weight[k] * score_by_sample.data[i][k] * score_by_sample.data[j][k];
         }
         sum_weight += importance_weight[k];
      }
      for(i=0; i<n_par; i++)  
      {
         score.data[i][0] /= sum_weight;
         for(j=0; j<n_par; j++)  
            score2.data[i][j] /= sum_weight;
      }
   }
   
   /*printf("\nestimate for covariance:\n");
   show_par_effective(par_effective, NULL);
   printf("\nscore:\n");
   fmprintf(&score);
   printf("\nscore2:\n");
   fmprintf(&score2);
   printf("\ninfo:\n");
   fmprintf(&info);
   exit(0);*/
     
   reset(&em_info, 0.0);
   for(i=0; i<n_par; i++)
   {
      for(j=i; j<n_par; j++)
      {
         em_info.data[i][j] = info.data[i][j] - score2.data[i][j] + score.data[i][0] * score.data[j][0];
         em_info.data[j][i] = em_info.data[i][j];
      }
   }

   //printf("\nem_info:\n");
   //fmprintf(&em_info);
   
   reset(&sub_info, 0.0);
   for(i=0; i<n_par_equiclass; i++)
   {
      for(j=i; j<n_par_equiclass; j++)
      {
         for(k=0; k<cfg_pars.par_equiclass[i].size; k++)
         {

            m = cfg_pars.par_equiclass[i].member[k] - 1;
            for(l=0; l<cfg_pars.par_equiclass[j].size; l++)
            {
               n = cfg_pars.par_equiclass[j].member[l] - 1;
               sub_info.data[i][j] += em_info.data[m][n];
            }
         }
         sub_info.data[j][i] = sub_info.data[i][j];
      }
   }
   //printf("\nsub_info:\n");
   //fmprintf(&sub_info);

   if(inverse(&sub_info,&sub_inv) > 0)
   {
      printf("covariance calculation: information matrix is not invertable\n");
      error = 1;
      goto end;
   }

   //parse the Fisher information to original parameters
   for(i=0; i<n_par_equiclass; i++)
   {
      for(j=i; j<n_par_equiclass; j++)
      {
         for(k=0; k<cfg_pars.par_equiclass[i].size; k++)
         {

            m = cfg_pars.par_equiclass[i].member[k] - 1;
            for(l=0; l<cfg_pars.par_equiclass[j].size; l++)
            {
               n = cfg_pars.par_equiclass[j].member[l] - 1;
               var_logit->data[m][n] = sub_inv.data[i][j];
               var_logit->data[n][m] = var_logit->data[m][n];
            }
         }
      }
   }
   //it is likely that if this function is only called for variance calculation,
   //the variance matrix will have negative diagonal elements because
   //the provided parameter values are not truely MLEs.
   for(m=0; m<n_par; m++)
   {
      if(var_logit->data[m][m] < 0)
      {
         error = 1;
         printf("Serial number=%d: Covariance matrix var_logit is not positive definite\n", serial_number);
         fmprintf(var_logit);
         break;
      }
   }


end:
   deflate_matrix(&first);
   deflate_matrix(&second);
   deflate_matrix(&score_grp);
   deflate_matrix(&info_grp);
   deflate_matrix(&score);
   deflate_matrix(&score2);
   deflate_matrix(&info);
   deflate_matrix(&em_info);
   deflate_matrix(&sub_info);
   deflate_matrix(&sub_inv);
   deflate_matrix(&score_by_sample);
   return(error);
}

// factor Monte Carlo error (MCE) into the final variance estimation
int adjust_for_mce(double *ini_par_effective, int optimization_method, COMMUNITY *aux_community, MATRIX *final_var_logit)
{
   int h, i, j ,k, l, m, n;
   int n_valid_est, n_valid_var;
   int *index;
   int n_par, n_par_equiclass, n_b_mode, n_p_mode, n_q_mode, n_u_mode;
   int iter, loop, converge, error;
   int max_EM_iter = 200;
   double *par_effective, *old_par_effective;
   double *bootstrap_mean, **bootstrap_var, **bootstrap_var_logit, **store, **store_logit;
   double **log_L_ini;
   double gap, max_log_L, *der;
   STATE *ptr_stat, *ptr2_stat;
   MATRIX *sample_states, var, var_logit, E_var_logit, E_var;
   PEOPLE *person;
   FILE *file;
   time_t t1, t2, tt1, tt2;
   clock_t T1, T2, TT1, TT2;
   int TIME_SAMPLING_SEC, TIME_SAMPLING_CPU;
   int TIME_ITERATION_SEC, TIME_ITERATION_CPU;
   int TIME_VARIANCE_SEC, TIME_VARIANCE_CPU;
   int TIME_TOTAL_SEC, TIME_TOTAL_CPU;

   n_b_mode = cfg_pars.n_b_mode;
   n_p_mode = cfg_pars.n_p_mode;
   n_u_mode = cfg_pars.n_u_mode;
   n_q_mode = cfg_pars.n_q_mode;
   n_par_equiclass = cfg_pars.n_par_equiclass;
   n_par = cfg_pars.n_par;

   initialize_matrix(&var_logit);
   initialize_matrix(&var);
   initialize_matrix(&E_var_logit);
   initialize_matrix(&E_var);
   inflate_matrix(&var_logit, n_par, n_par, 0.0);
   inflate_matrix(&var, n_par, n_par, 0.0);
   inflate_matrix(&E_var_logit, n_par, n_par, 0.0);
   inflate_matrix(&E_var, n_par, n_par, 0.0);

   index = (int *) malloc((size_t) size_sample_states * sizeof(int));
   make_1d_array_double(&der, n_par, 0.0);
   par_effective = (double *) malloc((size_t) n_par_equiclass * sizeof(double));
   old_par_effective = (double *) malloc((size_t) n_par_equiclass * sizeof(double));
   bootstrap_mean = (double *) malloc((size_t) n_par_equiclass * sizeof(double));
   make_2d_array_double(&bootstrap_var, n_par_equiclass, n_par_equiclass, 0.0);
   make_2d_array_double(&bootstrap_var_logit, n_par_equiclass, n_par_equiclass, 0.0);
   make_2d_array_double(&store, n_par_equiclass, cfg_pars.n_sampling_for_mce, 0.0);
   make_2d_array_double(&store_logit, n_par_equiclass, cfg_pars.n_sampling_for_mce, 0.0);


   // The matrix array "sample_state" records the sample state structure of all communities. 
   // For OEM communities, this structure is lost after covariance calculation (because
   // we replace it with a MCMC chain in the covariance function). We need to restore this structure
   // from "sample_state"
   // For MCEM communities, we resample state from "sample_state".
   // "log_L_ini" records the initial log-likelihood of each community evaluated at
   // the initial paramer values.
   sample_states = (MATRIX *) malloc((size_t) n_community * sizeof(MATRIX));
   log_L_ini = (double **) malloc((size_t) n_community * sizeof(double *));
   for(h=0; h<n_community; h++)
   {
      initialize_matrix(sample_states + h);
      log_L_ini[h] = NULL;
      if(community[h].size_possible_states > 1)
      {
         if(community[h].size_possible_states < cfg_pars.min_size_MCEM)
         {
            inflate_matrix(sample_states + h, community[h].size_impute, community[h].size_possible_states, 0.0);
            log_L_ini[h] = (double *) malloc((size_t) community[h].size_possible_states * sizeof(double));
            ptr_stat = community[h].list_states;
         }
         else
         { 
            inflate_matrix(sample_states + h, community[h].size_impute, size_sample_states, 0.0);
            log_L_ini[h] = (double *) malloc((size_t) size_sample_states * sizeof(double));
            ptr_stat = community[h].sample_states;
         }
         k = 0;
         while(ptr_stat != NULL)
         {
            for(l=0; l<community[h].size_impute; l++)
               sample_states[h].data[l][k] = ptr_stat->states[l];
            log_L_ini[h][k] = ptr_stat->log_L_ini;
            ptr_stat = ptr_stat->next;
            k++;
         }
      }
   }

   TIME_SAMPLING_SEC = TIME_SAMPLING_CPU = 0;
   TIME_ITERATION_SEC = TIME_ITERATION_CPU = 0;
   TIME_VARIANCE_SEC = TIME_VARIANCE_CPU = 0;
   TIME_TOTAL_SEC = TIME_TOTAL_CPU = 0;
   n_valid_est = n_valid_var = 0;
   tt1 = time(NULL);
   TT1 = clock();
   for(iter=0; iter<cfg_pars.n_sampling_for_mce; iter++)
   {
      t1 = time(NULL);
      T1 = clock();
      for(j=0; j<n_par_equiclass; j++)  par_effective[j] = ini_par_effective[j];

      //initialize MCEM communities.
      bootstrap(size_sample_states, size_sample_states, index, &seed);

      if(cfg_pars.use_bootstrap_for_mce == 1)
      {
         for(h=0; h<n_community; h++)
         if(community[h].size_possible_states >= cfg_pars.min_size_MCEM)
         {
            ptr_stat = community[h].sample_states;
            k = 0;
            while(ptr_stat != NULL)
            {
               m = index[k];
               for(l=0; l<community[h].size_impute; l++)
                  ptr_stat->states[l] = sample_states[h].data[l][m];
               ptr_stat->log_L_ini = log_L_ini[h][m];
               ptr_stat->log_L_cur = log_L_ini[h][m];
               ptr_stat->log_L_pre = log_L_ini[h][m];
               ptr_stat->ratio = 1.0;
               ptr_stat = ptr_stat->next;
               k++;
            }
         }
      }
      else
      {
         seed = 123456 + iter * 100;
         mt_init(seed);

         for(h=0; h<n_community; h++)
         if(community[h].size_possible_states >= cfg_pars.min_size_MCEM)
         {
            generate_sample_states(community+h, size_sample_states, cfg_pars.n_burnin_sampling, par_effective, community[h].sample_states);
         }
      }
      t2 = time(NULL);
      T2 = clock();
      
      TIME_SAMPLING_SEC += difftime(t2,t1);
      TIME_SAMPLING_CPU += T2 - T1;

      // Initialize OEM communities
      for(h=0; h<n_community; h++)
      {
         if(community[h].size_possible_states > 1 &&
            community[h].size_possible_states < cfg_pars.min_size_MCEM)
         {
            k = 0;
            ptr_stat = community[h].list_states;
            while(ptr_stat != NULL)
            {
               ptr_stat->log_L_ini = log_L_ini[h][k];
               ptr_stat->log_L_cur = ptr_stat->log_L_ini;
               ptr_stat->log_L_pre = ptr_stat->log_L_ini;
               ptr_stat = ptr_stat->next;
               k++;
            }
         }
      }

      for(i=0; i<size_sample_states; i++)  importance_weight[i] = 1.0;

      loop = 0;
      converge = 0;
      
      
      t1 = time(NULL);
      T1 = clock();
      while( converge == 0 && loop < max_EM_iter)
      {
         //printf("%d:%d\n", iter, loop);

         for(i=0; i<n_par_equiclass; i++)  
            old_par_effective[i] = par_effective[i];
         
         error = 0;
         error = M_step(par_effective, &max_log_L, optimization_method);
         update_weight(par_effective);
         //show_par_effective(par_effective, NULL);
         if(error > 0)  break;
         converge = 1;
         for(i=0; i<n_par_equiclass; i++)
         {
            m = cfg_pars.par_equiclass[i].member[0] - 1;
            if(m < n_b_mode + n_p_mode + n_u_mode + n_q_mode)
               gap = fabs(inv_logit(old_par_effective[i]) - inv_logit(par_effective[i])) 
                   / fabs(inv_logit(old_par_effective[i]));
            else
               gap = fabs(exp(old_par_effective[i]) - exp(par_effective[i])) / fabs(exp(old_par_effective[i]));
            if(gap > cfg_pars.converge_criteria[i])
               converge = 0;
         }
         //printf("converge=%d\n", converge);
         loop ++;
      }
      /*
      for(h=0; h<n_community; h++)
      {
         if(community[h].size_possible_states > 1 &&
            community[h].size_possible_states < cfg_pars.min_size_MCEM)
         printf("%d:%d  log_L_ini=%e log_L_cur=%e\n", h, community[h].size_possible_states, community[h].list_states->log_L_ini, community[h].list_states->log_L_cur);
      }
      exit(0);*/
      
      t2 = time(NULL);
      T2 = clock();
      TIME_ITERATION_SEC += difftime(t2,t1);
      TIME_ITERATION_CPU += T2 - T1;
      
      t1 = time(NULL);
      T1 = clock();
      if( converge == 1 && error == 0)
      {
         n_valid_est++;
         for(j=0; j<n_par_equiclass; j++)
         {
            store_logit[j][iter] = par_effective[j];
            m = cfg_pars.par_equiclass[j].member[0] - 1;
            if(m < n_b_mode + n_p_mode + n_u_mode + n_q_mode)   
               store[j][iter] = inv_logit(par_effective[j]);
            else
               store[j][iter] = exp(par_effective[j]);
         }
         // variance calculation is very time-consuming.
         // if user is only interested in evaluating the additional variance
         // due to MC sampling, not in the final variance estimation, 
         // this step can be skipped to save time.
         // Note, As E_var and E_var_logit are zero matrices if this step is skipped,
         // The returned final_var_logit would contain variance informtion only due to MC sampling.
         if(cfg_pars.skip_Evar_for_mce == 0)
         {
            error = covariance(par_effective, aux_community, &var_logit, 0);
            if(error == 0)
            {
               n_valid_var ++;
               addition(&E_var_logit, 1.0, &var_logit, 1.0, &E_var_logit);
               for(j=0; j<n_par_equiclass; j++)
               {
                  for(k=0; k<cfg_pars.par_equiclass[j].size; k++)
                  {
                     m = cfg_pars.par_equiclass[j].member[k] - 1;
                     if(m < n_b_mode + n_p_mode + n_u_mode + n_q_mode)   
                        der[m] = inv_logit(par_effective[j]) * (1 - inv_logit(par_effective[j]));
                     else  der[m] = exp(par_effective[j]);
                  }
               }
               for(m=0; m<n_par; m++)
               {
                  for(n=m; n<n_par; n++)
                  {
                     var.data[m][n] = var.data[n][m] = der[m] * var_logit.data[m][n] * der[n];
                  }
               }
               addition(&E_var, 1.0, &var, 1.0, &E_var);
            }
         }
      }
      t2 = time(NULL);
      T2 = clock();
      TIME_VARIANCE_SEC += difftime(t2,t1);
      TIME_VARIANCE_CPU += T2 - T1;
      if(cfg_pars.silent_run == 0)  
         printf("\n\nBootstrap %d: n_loop=%d  converge=%d:\n", iter, loop, converge);
      
   } //for(iter=0; iter<cfg_pars.n_sampling_for_mce; iter++)
   tt2 = time(NULL);
   TT2 = clock();

   TIME_TOTAL_SEC = difftime(tt2,tt1);
   TIME_TOTAL_CPU = TT2 - TT1;

   if(cfg_pars.silent_run == 0)
   printf("\n\n Overall: sampling=%d(%d) iteration=%d(%d)  variance=%d(%d)  total=%d(%d)\n",
                TIME_SAMPLING_SEC, TIME_SAMPLING_CPU, TIME_ITERATION_SEC, TIME_ITERATION_CPU,
                TIME_VARIANCE_SEC, TIME_VARIANCE_CPU, TIME_TOTAL_SEC, TIME_TOTAL_CPU);

   if(n_valid_var > 0)
   {   
      scalar_addition(&E_var_logit, 0.0, 1.0/n_valid_var,  &E_var_logit);
      scalar_addition(&E_var, 0.0, 1.0/n_valid_var,  &E_var);
   }
   for(i=0; i<n_par_equiclass; i++)
   {
      bootstrap_mean[i] = mean_double(store[i], n_valid_est);

      for(j=i; j<n_par_equiclass; j++)
      {
         bootstrap_var_logit[j][i] = bootstrap_var_logit[i][j] 
                                   = covariance_double(store_logit[i], store_logit[j], n_valid_est);
         bootstrap_var[j][i] = bootstrap_var[i][j] 
                                   = covariance_double(store[i], store[j], n_valid_est);
         for(k=0; k<cfg_pars.par_equiclass[i].size; k++)
         {
            m = cfg_pars.par_equiclass[i].member[k] - 1;
            for(l=0; l<cfg_pars.par_equiclass[j].size; l++)
            {
               n = cfg_pars.par_equiclass[j].member[l] - 1;
               if(final_var_logit != NULL)  
                  final_var_logit->data[n][m] = final_var_logit->data[m][n] 
                                              = bootstrap_var_logit[i][j] + E_var_logit.data[m][n];
            }
         }
      }
   }
   


   if(cfg_pars.silent_run == 0)
   {
      printf("E_var:\n");
      fmprintf(&E_var);
      printf("E_var_logit:\n");
      fmprintf(&E_var_logit);
      printf("valid bootstrap samples:%d\n", n_valid_est);
      printf("Bootstrap mean:\n");
      show_par_effective(bootstrap_mean, NULL);
      printf("Bootstrap var:\n");
      for(i=0; i<n_par_equiclass; i++)
      {
         for(j=0; j<n_par_equiclass; j++)
            printf("%e  ", bootstrap_var[i][j]);
         printf("\n");
      }
      printf("Bootstrap var_logit:\n");
      for(i=0; i<n_par_equiclass; i++)
      {
         for(j=0; j<n_par_equiclass; j++)
            printf("%e  ", bootstrap_var_logit[i][j]);
         printf("\n");
      }
   }
   // restore original sample state in each community
   /*
   for(h=0; h<n_community; h++)
   {
      if(community[h].size_impute > 0)
      {
         ptr_stat = community[h].sample_state;
         k = 0;
         while(ptr_stat != NULL)
         {
            for(l=0; l<community[h].size_impute; l++)
               ptr_stat->state[l] = sample_state[h].data[l][k];
            ptr_stat->log_L_ini = log_L_ini[h][k];
            ptr_stat = ptr_stat->next;
            k++;
         }
      }
   }
   */
   for(h=0; h<n_community; h++) 
   {
      deflate_matrix(sample_states + h);
      free(log_L_ini[h]);
   } 
   free(sample_states);
   free(log_L_ini);
   free(index);
   free(der);
   free(par_effective);
   free_2d_array_double(store);
   free_2d_array_double(store_logit);
   free(bootstrap_mean);
   free_2d_array_double(bootstrap_var);
   free_2d_array_double(bootstrap_var_logit);
   deflate_matrix(&E_var);
   deflate_matrix(&E_var_logit);
   deflate_matrix(&var);
   deflate_matrix(&var_logit);

   return(n_valid_est);
}   



int control(double *ini_par_effective, double *par_effective, double *max_log_L, int optimization_method)
{
   int i, j, h, k, l, m, n, loop;
   int n_par_equiclass, n_b_mode, n_p_mode, n_q_mode, n_u_mode;
   int converge, error=0;
   int max_EM_iter = 500;
   int MC_swamp, n_ind_sample, *index;
   time_t t1, t2;
   clock_t T1, T2;


   char file_name[50];

   double *old_par_effective;
   double *Q_stat, gap, temp;
   FILE *file;
   STATE *ptr_stat, *ptr2_stat;

   old_par_effective = NULL;
   TIME_SAMPLING_SEC = TIME_SAMPLING_CPU =
   TIME_ITERATION_SEC = TIME_ITERATION_CPU = 0;

   n_b_mode = cfg_pars.n_b_mode;
   n_p_mode = cfg_pars.n_p_mode;
   n_u_mode = cfg_pars.n_u_mode;
   n_q_mode = cfg_pars.n_q_mode;
   n_par_equiclass = cfg_pars.n_par_equiclass;
   
   //initialize size_sample_state and set all weights to be 1
   if(cfg_pars.EM == 1)  
   {
       size_sample_states = cfg_pars.n_base_sampling;   
       make_1d_array_double(&importance_weight, size_sample_states, 1.0);
   }    
   else  
   {
       size_sample_states = 0;
       importance_weight = NULL;
   }    

   //show_par_effective(par_effective, NULL);
   make_1d_array_double(&old_par_effective, n_par_equiclass, 0.0);
   //make_1d_array_double(&Q_stat, n_par, 0.0);
   
   if(cfg_pars.EM == 0)
      error = M_step(par_effective, max_log_L, optimization_method);
   else
   {
      t1 = time(NULL);   
      T1 = clock();   

      if(cfg_pars.silent_run == 0)  
         printf("Generate all possible state for communities asigned to EM algorithm\n");
      
      // For communities with small number of possible state,
      // we apply original EM (OEM) 
      for(h=0; h<n_community; h++)
      {
         //if(community[h].size_possible_states > 1)  printf("h=%d  size_possible_states=%d\n", h, community[h].size_possible_states);
         if(community[h].size_possible_states > 1 && community[h].size_possible_states < cfg_pars.min_size_MCEM)
         {
            create_state_chain(community[h].size_possible_states, community[h].size_impute, &community[h].list_states, &community[h].list_states_rear); 
            generate_list_states(community+h, par_effective);
         }
         if(community[h].size_possible_states >= cfg_pars.min_size_MCEM)
            create_state_chain(size_sample_states, community[h].size_impute, &community[h].sample_states, &community[h].sample_states_rear); 
      }
      
      // burnin process to get sensable initial estimates based on which
      // importance samples will be drawn
      if(n_need_MCEM > 0)
      {
         if(cfg_pars.silent_run == 0)  
            printf("Begin burn-in iterations to generate MCMC importance samples for MCEM\n");
         
         for(i=0; i<cfg_pars.n_burnin_iter; i++)
         {
            // create MCMC samples for communities with too many possible state.
            // We apply MCEM to these communities
            for(h=0; h<n_community; h++)
            if(community[h].size_possible_states >= cfg_pars.min_size_MCEM)
            {
               generate_sample_states(community+h, size_sample_states, cfg_pars.n_burnin_sampling, par_effective, community[h].sample_states); 
            }
            
            //printf("Estimation for initial iterations\n");
            if(i < cfg_pars.n_burnin_iter )
            {
               error = M_step(par_effective, max_log_L, optimization_method);
               update_weight(par_effective);
            }
            if(cfg_pars.silent_run == 0)  
            {
               printf("burn-in iteration: %d\n", i);
               show_par_effective(par_effective, NULL);
            }
            if(error > 0)  goto end;
         }
         // set initial estimates, and calculate likelihoods for samples based on initial estimates, 
         // which are used for calculating weights. 
         for(i=0; i<n_par_equiclass; i++)  
            ini_par_effective[i] = par_effective[i];
         //update log_L_ini for OEM communities.
         for(h=0; h<n_community; h++)
         if(community[h].size_impute > 0 && community[h].size_possible_states > 1 && 
            community[h].size_possible_states < cfg_pars.min_size_MCEM)
         {
            ptr_stat = community[h].list_states;
            while(ptr_stat != NULL)
            {
               ptr_stat->log_L_ini = ptr_stat->log_L_cur;
               ptr_stat = ptr_stat->next;
            }
         }
      }
      t2 = time(NULL);   
      T2 = clock();   
      TIME_SAMPLING_SEC = difftime(t2,t1);
      TIME_SAMPLING_CPU = T2 - T1;
      if(cfg_pars.silent_run == 0)  
         printf("It takes %d seconds (%d cpu clocks) to generate MCMC importance samples for MCEM\n", TIME_SAMPLING_SEC, TIME_SAMPLING_CPU);
     
      if(cfg_pars.silent_run == 0)
         printf("\n\nBegin MCEM iterations with fixed sample size\n");
      
      t1 = time(NULL);   
      T1 = clock();   
      loop = 0;
      converge = 0;
      while( converge == 0 && loop < max_EM_iter)
      {
         for(i=0; i<n_par_equiclass; i++)  old_par_effective[i] = par_effective[i];
         error = M_step(par_effective, max_log_L, optimization_method);
         update_weight(par_effective);
         if(error > 0)  
         {
            printf("Error in M-step\n");
            goto end;
         }

         converge = 1;
         for(i=0; i<n_par_equiclass; i++)
         {
            m = cfg_pars.par_equiclass[i].member[0] - 1;
            if(m < n_b_mode + n_p_mode + n_u_mode + n_q_mode)
               gap = fabs(inv_logit(old_par_effective[i]) - inv_logit(par_effective[i])) / fabs(inv_logit(old_par_effective[i]));
            else
               gap = fabs(exp(old_par_effective[i]) - exp(par_effective[i])) / fabs(exp(old_par_effective[i]));
            if(gap > cfg_pars.converge_criteria[i])
               converge = 0;
         }
         if(cfg_pars.silent_run == 0)  
         {
            printf("\n\nMCEM iteration %d\n", loop);
            show_par_effective(par_effective, NULL);
            printf("max_log_L=%e  Converge=%d\n", (* max_log_L), converge);
         }
         if(n_need_OEM == 0 && n_need_MCEM == 0)
         {
            converge = 1;
            break;
         }
         loop ++;
      }
      
      t2 = time(NULL);   
      T2 = clock();   
      TIME_ITERATION_SEC = difftime(t2,t1);
      TIME_ITERATION_CPU = T2 - T1;
      if(cfg_pars.silent_run == 0)  
         printf("It takes %d seconds and %d iterations to finish the EM iterations\n", TIME_ITERATION_SEC, loop);
      
      // output sample infectiousness onset days of asymptomatic cases
      // so that we can check how well they mix
      if(cfg_pars.check_mixing == 1)
      {         
         sprintf(file_name, "check_mixing.dat");
         file = fopen(file_name, "w");
         for(h=0; h<n_community; h++)
         if(community[h].size_possible_states >= cfg_pars.min_size_MCEM)
         {
            ptr_stat = community[h].sample_states; 
            while(ptr_stat != NULL)
            {
               for(m=0; m<community[h].size_impute; m++)
               {
                  i = community[h].member_impute[m];
                  l = ptr_stat->states[m];
                  if(people[i].possible_states[0][l] != INFINITY_INTEGER && people[i].possible_states[0][l] != -INFINITY_INTEGER)
                  {
                     fprintf(file, "%d  %d  %d\n", h, i, people[i].possible_states[0][l]);
                  } 
               }
               ptr_stat = ptr_stat->next;
            }
         }
         fclose(file);
      }

      if( converge == 0)
      {
          printf("MCEM iterations with increasing sample size does not converge\n");
          error = 1;
          goto end;
      }
      //t2 = time(NULL);
      //printf("\n use time: %f seconds\n", difftime(t2, t1));
      //exit(0);
      //first select independent subsamples for checking the need for sample size increase
      // Note: in Levine and Casella's paper, it appears that the independent subsamples
      // are reselected at each iteration, which doesn't make much sense because such re-selection
      // adds another layer of randomness and will not reflect the difference of Q statistics between
      // iterations due to EM algorithm. Therefore, we fix the independent subsamples for each sample size,
      // and reselect them only when sample size is changed.
      /*
      index = (int *) malloc((size_t) size_sample_state * sizeof(int));
      n_ind_sample = select_sub_sample(size_sample_state, index);
      augmentation(ini_par_effective, par_effective, n_ind_sample, index, Q_stat, 1);
      
      printf("\n\nBegin MCEM iterations with increasing sample size\n");
      loop = 0;
      converge = 0;
      while( converge == 0 && loop < max_EM_iter)
      {
         printf("\n\niteration %d\n", loop);

         // the iterations with fixed sample size already reached convergence,
         // so we need to consider augmentation first before the E-M step and checking convergence.
         // In the first loop, the sample size will be almost surely increased because the parameters
         // already converge under the fixed-sample-size setting.
         printf("Check augmentation\n");
         if(size_sample_state < 1000)  
         {
            MC_swamp = augmentation(ini_par_effective, par_effective, n_ind_sample, index, Q_stat, 0);
            //if sample size is increased, reselect independent subsamples
            if(MC_swamp == 1)
            {
               free(index);
               index = (int *) malloc((size_t) size_sample_state * sizeof(int));
            }
            n_ind_sample = select_sub_sample(size_sample_state, index);
            printf("n_ind_sample=%d\n", n_ind_sample);
         }
         update_weight(NULL);
         
         for(i=0; i<n_par_equiclass; i++)  old_par_effective[i] = par_effective[i];
         error = M_step(par_effective, max_log_L, optimization_method);
         if(error > 0)  
         {
            printf("Error in M-step\n");
            goto end;
         }
         show_par_effective(par_effective);

         converge = 1;
         for(i=0; i<n_par_equiclass; i++)
         {
            m = cfg_pars.par_equiclass[i].member[0] - 1;
            if(m < n_b_mode + n_p_mode + n_u_mode + n_q_mode)
               gap = fabs(inv_logit(old_par_effective[i]) - inv_logit(par_effective[i])) / fabs(inv_logit(old_par_effective[i]));
            else
               gap = fabs(exp(old_par_effective[i]) - exp(par_effective[i])) / fabs(exp(old_par_effective[i]));
            if(gap > cfg_pars.converge_criteria[i])
               converge = 0;
         }
         printf("Converge=%d\n", converge);
         
         loop ++;
      }
      if(loop >= max_EM_iter)
      {
         error = 1;
         printf("max number of iteration exceeded for the EM procedure: %d\n", loop);
         goto end;
      }
      */
   }

end:
   free(old_par_effective);   
   // free(Q_stat);
   //free(index);
   return(error);
}
  

     
int estimation(int id_inc, int id_inf, int id_time, double *est, double *logL, MATRIX *var, MATRIX *var_logit)
{
  int h, i, j, k, l, m, n, count;
  int error_type=0;
  int error = 0;
  int n_ini, n_par, n_par_equiclass, n_b_mode, n_p_mode, n_q_mode, n_u_mode;
  time_t t1, t2;
  clock_t T1, T2;

  double temp_logL, NR_logL, NM_logL;
  double *logL_store, *der;
  double ini_min, ini_max;
  double *par, *par_effective, *ini_par_effective, **par_store;
  double temp;
  COMMUNITY *aux_community = NULL;
  char file_name[50];
  FILE *file;

  n_b_mode = cfg_pars.n_b_mode;
  n_p_mode = cfg_pars.n_p_mode;
  n_u_mode = cfg_pars.n_u_mode;
  n_q_mode = cfg_pars.n_q_mode;
  n_par = cfg_pars.n_par;
  n_par_equiclass = cfg_pars.n_par_equiclass;
  n_ini = cfg_pars.n_ini;

  make_1d_array_double(&par, n_par, 0.0);
  make_1d_array_double(&par_effective, n_par_equiclass, 0.0);
  make_1d_array_double(&ini_par_effective, n_par_equiclass, 0.0);
  make_2d_array_double(&par_store, n_ini, n_par_equiclass, 0.0);
  make_1d_array_double(&logL_store, n_ini, 0.0);
  make_1d_array_double(&der, n_par, 0.0);

  // if there is no parameter to estimate, we return only log-likelihood
  if(n_par_equiclass == 0)
  {
     (* logL) = likelihood(par_effective, NULL);
     //force fixed parameters to be prespecified values
     if(cfg_pars.n_par_fixed > 0)
     {
        for(i=0; i<cfg_pars.n_par_fixed; i++)
        {
           j = cfg_pars.par_fixed_id[i] - 1;
           par[j] =  cfg_pars.par_fixed_value[i];
        }
     }
     for(j=0; j<n_par; j++)
     {
        if(j < n_b_mode + n_p_mode + n_u_mode + n_q_mode)  est[j] = inv_logit(par[j]);
        else est[j] = exp(par[j]);
     }
     goto end;
 }  
  // first start the Newton-Raphson emthod from user-specified initial estimates,
  // and select optimal result from these. If user does not specify,
  // the initial estimates

  if(cfg_pars.ini_par_provided == 0)
  {
     for(i=0; i<n_ini; i++)
     {
        for(j=0; j<n_par_equiclass; j++)
        {
           m = cfg_pars.par_equiclass[j].member[0] - 1;
           if(m < n_b_mode)
           {
              ini_min = logit(1e-6);
              ini_max = logit(1e-2);
           }
           else if(m < n_b_mode + n_p_mode)
           {
              ini_min = logit(1e-5);
              ini_max = logit(1e-1);
           }
           else if(m < n_b_mode + n_p_mode + n_u_mode)
           {
              ini_min = logit(5e-2);
              ini_max = logit(5e-1);
           }
           else if(m < n_b_mode + n_p_mode + n_u_mode + n_q_mode)
           {
              ini_min = logit(5e-2);
              ini_max = logit(5e-1);
           }
           else
           {
              ini_min = -5.0;
              ini_max = 5.0;
           }

           cfg_pars.ini_par_effective[i][j] = ini_min + (ini_max - ini_min) * runiform(&seed);
        }
     }
  }

  if(cfg_pars.search_bound_provided == 0)
  {
     for(j=0; j<n_par_equiclass; j++)
     {
        m = cfg_pars.par_equiclass[j].member[0] - 1;
        if(m < n_b_mode + n_p_mode + n_u_mode + n_q_mode)
        {
           cfg_pars.lower_search_bound[j] = close_to_0;
           cfg_pars.upper_search_bound[j] = close_to_1;
        }
        else
        {
           cfg_pars.lower_search_bound[j] = -1e10;
           cfg_pars.upper_search_bound[j] = 1e10;
        }
     }
  }

  if(cfg_pars.converge_criteria_provided == 0)
  {
     for(j=0; j<n_par_equiclass; j++)
        cfg_pars.converge_criteria[j] = 1e-8;
  }
  
  //printf("config_pars initial estimates\n");
  //for(i=0; i<n_ini; i++)
  //{
  //   for(j=0; j<n_par_equiclass; j++) printf("%e  ", cfg_pars.ini_par_effective[i][j]);
  //   printf("\n");
  //}
  //show_par_effective(cfg_pars.ini_par_effective[0], NULL);
  //for(h=0; h<n_community; h++)
  //if(h<1)
  //{
     //h = 360;
  //   temp_logL = community_log_likelihood(community + h, cfg_pars.ini_par_effective[0]);
  //   printf("HH %d: temp_logL=%e\n", h, temp_logL);
  //}
  //printf("hh 216: log_l=%e\n", community_log_likelihood(community + 216, cfg_pars.ini_par_effective[0]));
  //printf("hh 0: log_l=%e\n", community_log_likelihood(community + 0, cfg_pars.ini_par_effective[0]));
  //likelihood2(cfg_pars.ini_par_effective[0], 0, 106);
  //likelihood2(cfg_pars.ini_par_effective[0], 125, 231);
  //printf("old=%e  new=%e  all=%e\n", likelihood2(cfg_pars.ini_par_effective[0], 0, 124),
  //                                   likelihood2(cfg_pars.ini_par_effective[0], 125, 231),
  //                                   likelihood2(cfg_pars.ini_par_effective[0], 0, 231));
  
  switch(cfg_pars.optimization_choice)
  {
      case 0:  // default choice, start Newton-Raphson algorithm from multiple
               // initial points, choose the best one.
           // count is the number of valid Newton-Raphson search
           count = 0;
           for(i=0; i<n_ini; i++)
           {
              for(j=0; j<n_par_equiclass; j++)  par_effective[j] = cfg_pars.ini_par_effective[i][j];
              error = control(ini_par_effective, par_effective, &temp_logL, 1);
              if(cfg_pars.silent_run == 0)  
              {
                 printf("\nchoice 0 run %d: logL=%e  error=%d  parameters:\n", i, temp_logL, error);
                 show_par_effective(par_effective, NULL);
              }
              if(error == 0)
              {
                 for(j=0; j<n_par_equiclass; j++)
                    par_store[count][j] = par_effective[j];
                 logL_store[count] = temp_logL;
                 count++;
              }
           }
           if(count > 0)
           {
              k = 0;
              NR_logL = logL_store[0];

              for(i=0; i<count; i++)
              {
                 if(NR_logL < logL_store[i])
                 {
                    NR_logL = logL_store[i];
                    k = i;
                 }
              }
              for(j=0; j<n_par_equiclass; j++)  par_effective[j] = par_store[k][j];
              (* logL) = NR_logL;
           }
           else
           {
              error_type = 14;
              goto end;
           }
           break;

      case 1: //start Newton-Raphson algorithm from multiple
               // auto-generated initial points, choose the best one, then
               // start the simplex method.
           count = 0;
           for(i=0; i<n_ini; i++)
           {
              for(j=0; j<n_par_equiclass; j++)  par_effective[j] = cfg_pars.ini_par_effective[i][j];
              error = control(ini_par_effective, par_effective, &temp_logL, 1);
              if(cfg_pars.silent_run == 0)  
              {
                 if(cfg_pars.simulation == 0)
                 {
                    printf("\nchoice 1 ini run %d: logL=%e  error=%d\n", i, temp_logL, error);
                    show_par_effective(par_effective, NULL);
                 }
              }
              if(error == 0)
              {
                 for(j=0; j<n_par_equiclass; j++)
                    par_store[count][j] = par_effective[j];
                 logL_store[count] = temp_logL;
                 count++;
              }
           }

           if(count > 0)
           {
              k = 0;
              NR_logL = logL_store[0];

              for(i=0; i<count; i++)
              {
                 if(NR_logL < logL_store[i])
                 {
                    NR_logL = logL_store[i];
                    k = i;
                 }
              }
              for(j=0; j<n_par_equiclass; j++)  par_effective[j] = par_store[k][j];
              //error = control(ini_par_effective, par_effective, &NM_logL, 1);
              NM_logL = NR_logL;
              if(cfg_pars.silent_run == 0)  
              {
                 if(cfg_pars.simulation == 0)
                 {
                    printf("\nchoice 1: logL=%e\n", NM_logL);
                    show_par_effective(par_effective, NULL);
                 }
              }

              if(NM_logL > NR_logL)
              {
                 (* logL) = NM_logL;
              }
              else
              {
                 for(j=0; j<n_par_equiclass; j++)  par_effective[j] = par_store[k][j];
                 (* logL) = NR_logL;
              }
           }
           else
           {
              error_type = 14;
              goto end;
           }
           break;

      case 2:

           count = 0;
           for(i=0; i<n_ini; i++)
           {
              for(j=0; j<n_par_equiclass; j++)  par_effective[j] = cfg_pars.ini_par_effective[i][j];

              error = control(ini_par_effective, par_effective, &temp_logL, 2);
              if(cfg_pars.silent_run == 0)  
              {
                 if(cfg_pars.simulation == 0)
                 {
                    printf("\nchoice 2 run %d: logL=%e\n", i, temp_logL);
                    show_par_effective(par_effective, file);
                 }
              }
              if(error == 0)
              {
                 for(j=0; j<n_par_equiclass; j++)
                    par_store[i][j] = par_effective[j];

                 logL_store[i] = temp_logL;
              }
           }
           if(count > 0)
           {
              k = 0;
              NM_logL = logL_store[0];
              for(i=0; i<n_ini; i++)
              {
                 if(NM_logL < logL_store[i])
                 {
                    NM_logL = logL_store[i];
                    k = i;
                 }
              }
              for(j=0; j<n_par_equiclass; j++)  par_effective[j] = par_store[k][j];
              (*logL) = logL_store[k];
           }
           else
           {
              error_type = 14;
              goto end;
           }
           break;
  }

  for(j=0; j<n_par_equiclass; j++)
  {
     for(l=0; l<cfg_pars.par_equiclass[j].size; l++)
     {
        m = cfg_pars.par_equiclass[j].member[l] - 1;
        par[m] = par_effective[j];
     }
  }

  //force fixed parameters to be prespecified values
  if(cfg_pars.n_par_fixed > 0)
  {
     for(i=0; i<cfg_pars.n_par_fixed; i++)
     {
        j = cfg_pars.par_fixed_id[i] - 1;
        par[j] =  cfg_pars.par_fixed_value[i];
     }
  }
  if(cfg_pars.silent_run >= 0) 
  {
     printf("Parameter estimates: ");
     for(j=0; j<n_par; j++)
     {
        if(j < n_b_mode + n_p_mode + n_u_mode + n_q_mode)  est[j] = inv_logit(par[j]);
        else est[j] = exp(par[j]);
        printf("%e  ", est[j]);
     }
     printf("\n");
  }
  TIME_VARIANCE_SEC = TIME_VARIANCE_CPU = 0;
  if(cfg_pars.skip_variance == 0)
  {
     if(cfg_pars.silent_run == 0)  printf("Calculating variance\n");
     t1 = time(NULL);   
     T1 = clock();   
     // In EM-MCEM algorithm, we need to generate new sample states for estimating final variance. 
     // If interest is not final variance, but only to evaluate Monte Carlo error, 
     // then no need to generate new sample states (by setting n_sampling_for_mce > 0 and skip_Evar_for_mce == 1) 
     // the structure aux_comminty only contains new samples of importance samples. 
     if(cfg_pars.EM == 0)  aux_community = NULL;  
     else if(!(cfg_pars.n_sampling_for_mce > 0 && cfg_pars.skip_Evar_for_mce == 1))  
     {
        //printf("Generate auxiliary communities\n");
        aux_community = (COMMUNITY *) malloc((size_t) n_community * sizeof(COMMUNITY));
        for(h=0; h<n_community; h++)
        {
           aux_community[h].sample_states = NULL;
           aux_community[h].sample_states_rear = NULL;
           //printf("h=%d  size_possible_states=%d\n", h, community[h].size_possible_states);
           if(community[h].size_possible_states > 1)
           { 
              create_state_chain(size_sample_states,  community[h].size_impute, &aux_community[h].sample_states, &aux_community[h].sample_states_rear); 
              generate_sample_states(community + h, size_sample_states, cfg_pars.n_burnin_sampling, par_effective, aux_community[h].sample_states);
           }
        }
     }
     error = 0;
     if(cfg_pars.EM == 1 && cfg_pars.n_sampling_for_mce > 0)
     {
        //printf("Adjust for Marte Carlo error\n");
        k = adjust_for_mce(ini_par_effective, 1, aux_community, var_logit);
        if(k == 0)  error = 1;
     }
     else
     {
        //printf("Call covariance function\n");
        if(cfg_pars.EM == 1 && aux_community == NULL)
        {
           printf("To calculate variance-covariance, new sample states are needed for all all OEM and MCEM communities\n");
           goto end;
        }
        error =  covariance(par_effective, aux_community, var_logit, 1);
     }
     if(error > 0)  error_type = 15;
     //printf("check2\n");
        
     t2 = time(NULL);   
     T2 = clock();   
     TIME_VARIANCE_SEC = difftime(t2,t1);
     TIME_VARIANCE_CPU = T2 - T1;

     if(cfg_pars.silent_run == 0)
     {
        printf("\nvar_logit:\n");
        fmprintf(var_logit);
     }
      
     for(j=0; j<n_par_equiclass; j++)
     {
        for(k=0; k<cfg_pars.par_equiclass[j].size; k++)
        {
           m = cfg_pars.par_equiclass[j].member[k] - 1;
           if(m < n_b_mode + n_p_mode + n_u_mode + n_q_mode)   
              der[m] = inv_logit(par_effective[j]) * (1 - inv_logit(par_effective[j]));
           else  der[m] = exp(par_effective[j]);
        }
     }
     for(m=0; m<n_par; m++)
     {
        for(n=m; n<n_par; n++)
        {
           var->data[m][n] = var->data[n][m] = der[m] * var_logit->data[m][n] * der[n];
        }
     }
     
     if(cfg_pars.silent_run == 0)
     {
        printf("\nvar:\n");
        fmprintf(var);
     }
  }
  // Goodness of fit
  // Warning: do not run goodness-of-fit if there is uncertainty about infection status, e.g., not sure if a person is an infection or an escape.
  // If uncertainty is only in infectiousness onset time of asymptomatic infection, or only between escape and pre-immunity, then it is fine.
  if(cfg_pars.goodness_of_fit == 1)
  {
     if(cfg_pars.silent_run == 0)   printf("Goodness of fit\n");
     if(cfg_pars.EM == 0)
     {
        goodness_of_fit(id_inc, id_inf, id_time, est, NULL, var_logit);
     }
     else
     {
        if(aux_community != NULL) 
           goodness_of_fit(id_inc, id_inf, id_time, est, aux_community, var_logit);
        else
        {
           printf("For assessing goodness-of-fit under EM, the auxiliary community structure containing new importance samples should not be NULL\n");
        }
     }
  }
  if(cfg_pars.check_runtime == 1)
  {
     sprintf(file_name, "check_runtime.dat");
     file = fopen(file_name, "w");
     fprintf(file, "%ld  %ld  %ld  %ld  %ld  %ld\n", TIME_SAMPLING_SEC, TIME_SAMPLING_CPU, 
                   TIME_ITERATION_SEC, TIME_ITERATION_CPU, TIME_VARIANCE_SEC, TIME_VARIANCE_CPU);
     fclose(file);
  }


end:
  free(par);
  free_2d_array_double(par_store);
  free(logL_store);
  free(par_effective);
  free(ini_par_effective);   
  free(der);
  if(aux_community != NULL)
  {
     for(h=0; h<n_community; h++)
     if(community[h].size_possible_states > 1)
        free_state_chain(aux_community[h].sample_states, aux_community[h].sample_states_rear); 
     free(aux_community);
  }

  return(error_type);
}

