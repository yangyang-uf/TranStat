/* this program is using resampling to estimate the null distribution of
   likelihood ratio statistics. Under the null hypothesis p=0 (no person-to-person
   transmission, we probably can move the symptom onset (or infectiousness onset)
   time around, while keeping the number of cases and the sum of time to symptom
   onset over all cases fixed. By this way, each resampled data will give about the same
   mle for the CSPI */

   
#include "test_estimation_cs.h"
#include "test_estimation_cs_hh.h"
   
unsigned long n_balls_to_m_boxes(int n, int m, int volume, MATRIX *store)
{
   int i, j ,k;
   unsigned long sum;
        
   if(n == 0)   
   {
      store->data[n][m] = 1;
      return(1);
   }
   if(m == 1)
   {
      if(n <= volume)
      {
         store->data[n][m] = 1;
         return(1);
      }
      else
      {
         store->data[n][m] = 0;
         return(0);
      }
   }
   if(store->data[n][m] != -1)  return((unsigned long) store->data[n][m]);
   
   sum = 0;
   k = min(n, volume);
   for(i=0; i<=k; i++)
   {
      sum +=  n_balls_to_m_boxes(n-i, m-1, volume, store);
   }
   store->data[n][m] = sum;
   
   return(sum);
}                      

int resampling_cs_hh(int n_resampling, double stat_obs, double *pvalue)
{
  int i, j, h, k, l, m, n, t, iter, count, width, baseline;
  int day_last_escape, lower, upper, escape_time, infection_time;
  int n_par_alternative = 2;
  int n_par_null = 1;
  int n_valid_sample;
  int n_ini=5;
  int n_print = 50;
  int error = 0;
  int volume, n_choice, choice;
  int available, inventory;
  int estimation_error_alternative, estimation_error_null;
  int n_eligible_subject, n_eligible_case, n_eligible_var_case, n_eligible_fix_case;
  int *eligible, *eligible_subjects, *eligible_cases, *eligible_var_cases, *eligible_fix_cases;
  int *day_ill_var, *day_exit_var, *day_ill_fix, *day_exit_fix, *offset_fix;
  unsigned int *sample, *sample_eligible_cases;
  int *allocate, *candidate;
  double *prob;
  double *est_alternative, *est_null;
  double s, temp, stat_null, logL_null, logL_alternative;
  PEOPLE *person;
  MATRIX store;
  FILE *file;

  eligible = eligible_subjects = eligible_cases = eligible_var_cases = eligible_fix_cases = NULL;
  day_ill_var = day_exit_var = day_ill_fix = day_exit_fix = offset_fix = NULL;
  sample = sample_eligible_cases = NULL;
  allocate = candidate = NULL;
  prob = est_alternative = est_null = NULL;
  
  make_1d_array_int(&eligible, p_size, 1);  
  
  // In case-ascertained design, non-index cases are called contact cases.
  // In prospective study, we do not distiguish between index cases and contact cases.
  // "eligible subjects" means all people who are not index cases
  // "eligible cases" are cases who are not index cases
  // "eligible variable cases" are eligible cases whose symptom onsets are sufficiently distant
  //  from the onset of index case, so that their symptom onset dates can be rearranged without any restriction.
  //  "eligible fixed cases" are eligible cases whose symptom onsets are not distant from index case's onset
  //  so that changing their onset date will change the null likelihood.
  
  n_eligible_subject = n_eligible_case = n_eligible_fix_case = n_eligible_var_case = 0;
  for(i=0, person=people; i<p_size; i++, person++)
  {
     h = person->community;

     eligible[i] = 1;
     if(cfg_pars.adjust_for_left_truncation == 1)
     {
        if(person->idx == 1)  eligible[i] = 0;
        else if(person->symptom == 1)
        {  
           day_last_escape = person->day_ill - cfg_pars.max_incubation - 1;
           lower = community[h].latest_idx_day_ill - cfg_pars.min_incubation;
           upper = community[h].day_epi_stop - cfg_pars.max_incubation - 1;
           //if(person->id==928) printf("day_ill=%d  day_last_escape=%d  lower=%d\n", person->day_ill, day_last_escape, lower);
           if(day_last_escape >= lower && day_last_escape <= upper)  
              eligible[i] = 2;
           else  eligible[i] = 3;
        }
     }
     else     
     {
        if(person->symptom == 1)
        { 
           day_last_escape = person->day_ill - cfg_pars.max_incubation - 1;
           lower = community[h].day_epi_start - 1;     
           upper = community[h].day_epi_stop - cfg_pars.max_incubation - 1;
           if(day_last_escape >= lower && day_last_escape <= upper)  
              eligible[i] = 2;
           else  eligible[i] = 3;
        }
     }
     if(eligible[i] > 0)
     {
        n_eligible_subject ++;     
        if(eligible[i] > 1)  n_eligible_case++;
        if(eligible[i] == 2)  n_eligible_var_case++;
        if(eligible[i] == 3)  n_eligible_fix_case++;
     }   
  }
  //printf("n_eligible_subject=%d  n_eligible_case=%d  n_eligible_var_case=%d  n_eligible_fix_case=%d\n", 
  //        n_eligible_subject, n_eligible_case, n_eligible_var_case, n_eligible_fix_case);
  if(n_eligible_case == 0) 
  {
     (* pvalue) = -1;
     error = 1;
     printf("No eligible case to perform the test\n");
     goto end;
  }
  
  make_1d_array_int(&eligible_subjects, n_eligible_subject, 0);  
  make_1d_array_int(&eligible_cases, n_eligible_case, 0);  
  if(n_eligible_fix_case > 0) 
  {
     make_1d_array_int(&eligible_fix_cases, n_eligible_fix_case, 0);  
     make_1d_array_int(&day_ill_fix, n_eligible_fix_case, 0);  
     make_1d_array_int(&day_exit_fix, n_eligible_fix_case, 0);  
     make_1d_array_int(&offset_fix, n_eligible_fix_case, 0);  
  }  
  if(n_eligible_var_case > 0) 
  {
     make_1d_array_int(&eligible_var_cases, n_eligible_var_case, 0);  
     make_1d_array_int(&day_ill_var, n_eligible_var_case, 0);  
     make_1d_array_int(&day_exit_var, n_eligible_var_case, 0);  
  }  
  //file = fopen("check_test.txt", "w");
  j = k = l = m = 0;
  for(i=0, person=people; i<p_size; i++, person++)
  {
     //fprintf(file, "id:%d  eligibility:%d\n", person->id, eligible[i]);
     h = person->community;
     if(eligible[i] > 0)
     {
        eligible_subjects[j++] = i;
        if(eligible[i] > 1)
        {        
           eligible_cases[k++] = i;
           if(eligible[i] == 2) 
           {
              eligible_var_cases[l] = i;
              day_ill_var[l] = person->day_ill; 
              day_exit_var[l] = person->day_exit;
              l++;
           }
           if(eligible[i] == 3) 
           {
              eligible_fix_cases[m] = i;
              day_ill_fix[m] = person->day_ill; 
              day_exit_fix[m] = person->day_exit;
              if(cfg_pars.adjust_for_left_truncation == 1)
                 offset_fix[m] = person->day_ill - community[h].latest_idx_day_ill;
              else
                 offset_fix[m] = person->day_ill - community[h].day_epi_start;
              m++;
           }
        }           
     }   
  }
  //fclose(file);
  //exit(0);
  initialize_matrix(&store);

  est_alternative = (double *) malloc((size_t) (n_par_alternative * sizeof(double)));  
  est_null = (double *) malloc((size_t) (n_par_null * sizeof(double)));  
 
  
  sample = (unsigned int *) malloc((size_t) (n_eligible_case * sizeof(unsigned int)));
  sample_eligible_cases = (unsigned int *) malloc((size_t) (n_eligible_case * sizeof(unsigned int)));
  
  if(n_eligible_var_case > 0)
  {
     make_1d_array_int(&allocate, n_eligible_var_case, 0);
  }

  make_1d_array_double(&prob, 1000, 0);
  make_1d_array_int(&candidate, 1000, 0);
  
  // instead of rearranging symptom onset days directly, we rearrange the earliest
  //   escape time, which is safer and more convenient. The re-arranged earliest escape time
  //   will plus the original gap (array max_latent) to get the new symptom onset time.
      
  // for now, volume has to be the same for all cases. This implies epi_duraiton should be equal for all communities.
  available = volume = 0;
  //printf("n_eligible_contact_case=%d\n", n_eligible_contact_case);
  for(i=0; i<n_eligible_var_case; i++)
  {
     person = people + eligible_var_cases[i];
     h = person->community;
     day_last_escape = person->day_ill - cfg_pars.max_incubation - 1;
     if(cfg_pars.adjust_for_left_truncation == 1)
     {
        lower = community[h].latest_idx_day_ill - cfg_pars.min_incubation;
        upper = community[h].day_epi_stop - cfg_pars.max_incubation - 1;
     }   
     else
     {     
        lower = community[h].day_epi_start - 1;     
        upper = community[h].day_epi_stop - cfg_pars.max_incubation - 1;
     }   
     width = day_last_escape - lower + 1;
     available += width;
     volume = max(volume, upper - lower + 1);
     //printf("person %d  day_ill=%d  idx_day_ill=%d  lower=%d  day_exit=%d  width=%d\n", 
     //        person->id, person->day_ill, community[h].latest_idx_day_ill, lower, person->day_exit, width);
  }
  inventory = available;      
  printf("volume=%d  available=%d\n", volume, available);
  
  if(available > 0)
  {
     inflate_matrix(&store, available+1, n_eligible_var_case+1, -1);
  
     n = available;
     m = n_eligible_var_case;
     k = n_balls_to_m_boxes(n, m, volume, &store);
  }
  //dmprintf(&store);      
  
  for(i=0; i<n_eligible_case; i++)
  {
     person = people + eligible_cases[i];
     person->symptom = 0;
     person->day_ill = INFINITY_INTEGER;
     person->exit = 0;
     person->day_exit = INFINITY_INTEGER;
  }
  
  count = 0; /*count records the number of null stat smaller than observed stat*/
  n_valid_sample = 0; /* num_rand is the number of successful resamplings*/
  for(iter=0; iter<n_resampling; iter++)
  {
     if((iter % n_print) == 0) printf("%d\n", iter);
     
     permute(n_eligible_subject, n_eligible_case, sample, &seed);
     for(i=0; i<n_eligible_case; i++)
     {
        j = sample[i];
        sample_eligible_cases[i] = eligible_subjects[j];
     }      
     if(n_eligible_var_case > 0)
     {
        available = inventory;
        for(i=0; i<n_eligible_var_case; i++)  allocate[i] = 0;
        if(available > 0)
        {
           t = 0;
           for(i=0; i<n_eligible_var_case; i++)
           {
              if(t < n_eligible_var_case - 1)
              {
                 k = min(available, volume);
                 n_choice = 0;
                 //calculate weights for choosing each possible allocation
                 for(j=0; j<=k; j++)
                 {
                    n = available - j;
                    m = n_eligible_var_case - t - 1;
                    if(store.data[n][m] > 0)
                    {
                       prob[n_choice] = store.data[n][m];
                       candidate[n_choice] = j;
                       n_choice++;
                    }
                 }
                 choice = discrete(n_choice, prob, &seed); 
                 allocate[i] = candidate[choice];
                 available -= allocate[i];     
              }
              else
              {
                 if(available > volume)
                 {
                    printf("This should ever happer\n");
                    error = 1;
                    goto end;;
                 }   
                 allocate[i] = available;
              }   
              t++;
           }
        }  
        for(i=0; i<n_eligible_var_case; i++)
        {
           j = sample_eligible_cases[i];
           person = people + j;
           h = person->community;
           // baseline is the earliest escape time
           if(cfg_pars.adjust_for_left_truncation == 1)
              baseline = community[h].latest_idx_day_ill - cfg_pars.min_incubation;
           else   
              baseline = community[h].day_epi_start - 1;
           escape_time = baseline + allocate[i];
           /* earliest infection time should be the allocated escape days plus 1*/
           infection_time = escape_time + 1;
           person->day_ill = infection_time + cfg_pars.max_incubation;  
           person->symptom = 1;
           if(day_exit_var[i] < INFINITY_INTEGER)   
           {     
              person->day_exit = person->day_ill+ (day_exit_var[i] - day_ill_var[i]);
              person->exit = 1;
           }
           else
              person->day_exit =  INFINITY_INTEGER;  
           //printf("person %d  day_ill=%d  idx_day_ill=%d  allocate=%d  day_exit=%d  hh size=%d\n", 
           //        person->id, person->day_ill, community[h].latest_idx_day_ill, 
           //        allocate[i], person->day_exit, community[h].size);
        }

     }   
     
     if(n_eligible_fix_case > 0)
     {
        for(i=0; i<n_eligible_fix_case; i++)
        {
           j = sample_eligible_cases[n_eligible_var_case + i];
           person = people + j;
           h = person->community;
           // baseline is just a reference time
           if(cfg_pars.adjust_for_left_truncation == 1)
              baseline = community[h].latest_idx_day_ill;
           else   
              baseline = community[h].day_epi_start;
           person->day_ill = baseline + offset_fix[i];  
           person->symptom = 1;
           if(day_exit_fix[i] < INFINITY_INTEGER)   
           {     
              person->day_exit = person->day_ill + (day_exit_fix[i] - day_ill_fix[i]);
              person->exit = 1;
           }
           else
              person->day_exit =  INFINITY_INTEGER;  
           //printf("person %d  day_ill=%d  idx_day_ill=%d  allocate=%d  day_exit=%d  hh size=%d\n", 
           //        person->id, person->day_ill, community[h].latest_idx_day_ill, 
           //        allocate[i], person->day_exit, community[h].size);
        }
     }
     estimation_error_alternative = estimation_cs_hh(est_alternative, &logL_alternative, n_ini, 0);
     estimation_error_null = estimation_cs(est_null, &logL_null, n_ini);
     //if((iter % n_print) == 0) printf("estimation error: null=%d  alternative=%d\n", estimation_error_null, estimation_error_alternative);
     if(estimation_error_alternative == 0 && estimation_error_null == 0)
     {
        if(est_alternative[1] <= close_to_0)
           stat_null = 0.0;
        else
           stat_null = -2.0 * (logL_null - logL_alternative);
        if(stat_null < 0.0)  stat_null = 0.0; 
        if(stat_null < stat_obs) count++;
        n_valid_sample ++;
        if((iter % n_print) == 0)
        {
           //printf("Null model: b=%e\n", est_null[0]);
           //printf("Full model: b=%e  p=%e\n", est_alternative[0], est_alternative[1]);
           //printf("log_L_null=%e  log_L_alternative=%e  stat_null=%e  count=%d  n_valid=%d\n", 
           //        logL_null, logL_alternative, stat_null, count, n_valid_sample);
        }
     }   
     
     for(i=0; i<n_eligible_case; i++)
     {
        person = people + sample_eligible_cases[i]; 
        person->symptom = 0;
        person->day_ill = INFINITY_INTEGER;
        person->exit = 0;
        person->day_exit = INFINITY_INTEGER;
     }
  }
  
  (*pvalue) = 1.0 - (double) count / (double) n_valid_sample;   
  for(i=0; i<n_eligible_var_case; i++)
  {
     person = people + eligible_var_cases[i];
     person->day_ill = day_ill_var[i];
     person->day_exit = day_exit_var[i];
     person->symptom = 1;
     if(person->day_exit < INFINITY_INTEGER)  person->exit = 1;
  }
  for(i=0; i<n_eligible_fix_case; i++)
  {
     person = people + eligible_fix_cases[i];
     person->day_ill = day_ill_fix[i];
     person->day_exit = day_exit_fix[i];
     person->symptom = 1;
     if(person->day_exit < INFINITY_INTEGER)  person->exit = 1;
  }

end:  
  if(eligible != NULL) free(eligible);
  if(eligible_subjects != NULL) free(eligible_subjects);
  if(eligible_cases != NULL) free(eligible_cases);
  if(eligible_fix_cases != NULL) free(eligible_fix_cases);
  if(eligible_var_cases != NULL) free(eligible_var_cases);
  if(day_ill_var != NULL) free(day_ill_var);
  if(day_exit_var != NULL) free(day_exit_var);  
  if(day_ill_fix != NULL) free(day_ill_fix);
  if(day_exit_fix != NULL) free(day_exit_fix);  
  if(offset_fix != NULL) free(offset_fix);  
  if(sample != NULL) free(sample);
  if(sample_eligible_cases != NULL) free(sample_eligible_cases);
  if(est_alternative != NULL) free(est_alternative);
  if(est_null != NULL) free(est_null);
   
  if(prob != NULL) free(prob);
  if(allocate != NULL) free(allocate);
  if(candidate != NULL) free(candidate);

  deflate_matrix(&store);

  
  return error;
}  
        
 


double stat_test_cs_hh(void)
{

   int i, j, k, l, count;
   int error, error_null, error_alternative;
   int n_asym = 10000;
   int n_resampling = 2000;
   int n_ini =3;
   double est_null[1], est_alternative[2];
   double temp, stat_obs, p_value_asym, p_value_resample;
   double log_L_null, log_L_alternative;
   double *chisq, *stat_null;
  
   // asymptotic null distribution
   chisq = (double *)malloc((size_t) (n_asym * sizeof(double)));
   stat_null = (double *)malloc((size_t) (n_asym * sizeof(double)));
   for(k=0; k<n_asym; k++)
   {
     chisq[k] = rchisq(1.0, &seed);
     temp = runiform(&seed);
     if( temp < 0.5)
        stat_null[k] = chisq[k];
     else
        stat_null[k] = 0;
   }

   error_null = estimation_cs(est_null, &log_L_null, n_ini);
   printf("Null model: b=%e\n", est_null[0]);
   error_alternative = estimation_cs_hh(est_alternative, &log_L_alternative, n_ini, 0);
   printf("Full model: b=%e  p=%e\n", est_alternative[0], est_alternative[1]);
   
   if(error_null == 1)
   {
      printf("Unable to find the null statistic\n");
      return(-1);
   }
   else if(error_alternative == 1)
   {
      printf("Unable to find the alternative statistic\n");
      return(-1);
   }
   else
   {
   
      stat_obs = -2.0 * (log_L_null - log_L_alternative);
      if(stat_obs < 0.0)  stat_obs = 0.0;
      
      count = 0;
      for(k=0; k<n_asym; k++)  if(stat_obs > stat_null[k])  count ++;
      p_value_asym = 1 - (double) count / (double) n_asym;
      printf("stat_obs=%e  asymptotic p-value=%e\n", stat_obs, p_value_asym);
      
      error = resampling_cs_hh(n_resampling, stat_obs, &p_value_resample);
      if(error == 1)  
      {
         printf("Error in resampling test\n");
         return (-1);
      }   
      printf("resampling p-value=%e\n", p_value_resample);
   }
   free(chisq);
   free(stat_null);
   return(p_value_resample);
   
}

