

/* this program is for estimating MLE for the model with
within-household person-to-person transmission probability p1
and common source probability of infection (CSPI) b. */

int newton_raphson_cs(double *est, double *log_likelihood)
{
  int i, j, h, k, l, m, n, r, t, error=0;
  int max_iter = 200;
  int start_day, stop_day;
  double factor1, factor2, factor3;
  int len, inf1, inf2;
  int loop, id, index, skip;
  double cutpoint, L, log_L, log_L_all, s;
  double e, f;
  double b;
  double b_old;
  double lb; /* logit tansformation of b, p and q*/
  double b_lb;
  double b_lb_lb;
  double f1, f2;
  double log_f_b;
  double log_f_b_b;
  double log_e_b;
  double log_e_b_b;
  double log_L_b;
  double log_L_b_b;
  double score_lb;
  double info_lb_lb;
  double *ini_ee, *ee, *log_ee_b;
  double *log_ee_b_b;
  double *cum_logee, *cum_logee_b;
  double *cum_logee_b_b;
  double cum_e, cum_loge_b;
  double cum_loge_b_b;
  double temp, temp_b;
  double temp_b_b;
  double L_b;
  double L_b_b;
  double day_L;
  double log_day_L_b;
  double log_day_L_b_b;
  double score_b;
  double info_b_b;
  double pr;
  double p_inf, p_esc;
  double info, inv, score, d, dd;
  PEOPLE *person, *member;
  FILE *file;
  
 
  max_epi_duration = 0;
  for(h=0; h<n_community; h++) max_epi_duration = max(max_epi_duration, community[h].epi_duration);
  
  ee = (double *)malloc((size_t) (max_epi_duration * sizeof(double)));
  log_ee_b = (double *)malloc((size_t) (max_epi_duration * sizeof(double)));
  log_ee_b_b = (double *)malloc((size_t) (max_epi_duration * sizeof(double)));
  cum_logee = (double *)malloc((size_t) (max_epi_duration * sizeof(double)));
  cum_logee_b = (double *)malloc((size_t) (max_epi_duration * sizeof(double)));
  cum_logee_b_b = (double *)malloc((size_t) (max_epi_duration * sizeof(double)));

  
  /* p1 is tran prob in house, b is tran prob from comon source*/

  b = est[0]; 
  cutpoint = 1e-10;
  loop = 0;
  dd = 1.0;
 
 
  while( (dd > cutpoint) && (loop < max_iter))
  {
     loop ++;
     score_b = 0.0;
     info_b_b = 0.0;
     log_L_all = 0.0;
     
     for(h=0; h<n_community; h++)
     {
        if(cfg_pars.adjust_for_left_truncation == 1)
            start_day = community[h].latest_idx_day_ill - cfg_pars.min_incubation + 1;
        else  start_day = community[h].day_epi_start;
        stop_day = community[h].day_epi_stop;
        for(t=start_day; t<=stop_day; t++)  
        {
           r = t - start_day;
           ee[r] = 1.0;
           log_ee_b[r] = 0.0;
           log_ee_b_b[r] = 0.0;
        
           f = 1.0 - b;
           if(f <= 0.0 || f > 1.0)
           {
              /*printf("Full model: f=%e  b=%e\n", f, b);*/
              error = 1;
              goto end;
           }
           log_f_b = -1.0 / f;
           log_f_b_b = -1.0 / ( f * f );

           ee[r] = f;
           
           log_ee_b[r] = log_f_b;
           log_ee_b_b[r] = log_f_b_b;
        }     
     
        for(i=0; i<community[h].size; i++)
        {
           person = people + community[h].member[i];
           if(person->idx == 0)
           {
              if(cfg_pars.adjust_for_left_truncation == 1)
                 start_day = community[h].latest_idx_day_ill - cfg_pars.min_incubation + 1;
              else  start_day = community[h].day_epi_start;
              stop_day = (person->symptom == 1)? (person->day_ill - cfg_pars.min_incubation):community[h].day_epi_stop;
              //printf("start_day=%d  stop_day=%d\n", start_day, stop_day);
              temp = 0.0;
              temp_b = temp_b_b = 0.0;

              L = 0.0;
              L_b = L_b_b = 0.0;

              log_L = 0.0;
              log_L_b = log_L_b_b = 0.0;

              for(t=start_day; t<=stop_day; t++)
              {
                 r = t - start_day;
                 temp += log(ee[r]);
                 cum_logee[r] = temp;

                 temp_b += log_ee_b[r];
                 cum_logee_b[r] = temp_b;

                 temp_b_b += log_ee_b_b[r];
                 cum_logee_b_b[r] = temp_b_b;
              }
           
              if(person->symptom == 1) /* if he was ill */
              {
                 inf1 = person->day_ill - cfg_pars.max_incubation;
                 inf2 = person->day_ill - cfg_pars.min_incubation;
                 if(inf1 < start_day)
                    inf1 = start_day;
                 if(inf2 < start_day)
                    inf2 = start_day;

                 cum_e = 1.0; 
                 cum_loge_b = cum_loge_b_b = 0.0; 
                 for(t=inf1; t<=inf2; t++)
                 {
                    r = t - start_day;
                    pr = pdf_incubation[inf2-t];
                    p_inf = (fabs(log(ee[r]))<1e-6)? (-log(ee[r])) : (1.0 - ee[r]);
                    if(ee[r] < 0.0 || ee[r] > 1.0)
                    {
                       /*printf("error in ee[nd], nd=%d ee[nd]=%f\n", nd, ee[nd]);*/
                       error = 1;
                       goto end;
                    }
                 
                    if(ee[r] >= 0.0 && ee[r] < 1.0)
                    {
                       factor1 = -ee[r] / p_inf;
                       factor2 = factor1 / p_inf;
                       
                       day_L = cum_e * p_inf;
                       log_day_L_b = cum_loge_b + factor1 * log_ee_b[r];
                       log_day_L_b_b = cum_loge_b_b + factor1 * log_ee_b_b[r] + factor2 * log_ee_b[r] * log_ee_b[r];

                       temp = day_L * pr; 
                       L += temp;
                       L_b += log_day_L_b * temp;
                       L_b_b += (log_day_L_b * log_day_L_b + log_day_L_b_b) * temp;
                    }
                           
                    cum_e *= ee[r];
                    cum_loge_b += log_ee_b[r];
                    cum_loge_b_b += log_ee_b_b[r];
                 }
                 if(L <= 0.0)
                 {
                    /*
                    printf("Likelihood must be positive!  L=%f\n", L);
                    printf("id=%d  start_day=%d  day_ill=%d  inf1=%d  inf2=%d\n", i, start_day, person->day_ill, inf1, inf2);
                    printf("b=%e  p1=%e\n", b, p1);
                    for(t=inf1;t<=inf2;t++)
                    {
                       nd = t - start_day;
                       printf("ee[nd]=%e 1-ee[nd]=%e\n", ee[nd], 1.0-ee[nd]);
                    }
                    */   
                    error = 1;
                    goto end;
                 }
                 temp = 1.0 / L;
                 log_L = log(L);
                 log_L_b = L_b * temp;
                 log_L_b_b = L_b_b * temp - log_L_b * log_L_b;

                 r = inf1 - start_day; 
                 if(r > 0)
                 {   
                    log_L += cum_logee[r-1];
                    log_L_b += cum_logee_b[r-1];
                    log_L_b_b += cum_logee_b_b[r-1];
                 }
              }
              else  /* if never been ill */
              {
                 if(cfg_pars.adjust_for_right_censoring == 0)
                 {
                    r = stop_day - start_day;
                    log_L = cum_logee[r];
                    log_L_b = cum_logee_b[r];
                    log_L_b_b = cum_logee_b_b[r];
                 }
                 else
                 {
                    inf1 = stop_day - cfg_pars.max_incubation + 1;
                    inf2 = stop_day - cfg_pars.min_incubation;
                    if(inf1 < start_day)
                       inf1 = start_day;
                    if(inf2 < start_day)
                       inf2 = start_day;
              
                    /* potential infection during days from inf1 to inf2*/
                    cum_e = 1.0; 
                    cum_loge_b = cum_loge_b_b = 0.0; 
                    for(t=inf1; t<=inf2; t++)
                    {
                       pr = sdf_incubation[inf2-t];
                       r = t - start_day;
                       if(ee[r] < 0.0 || ee[r] > 1.0)
                       {
                          /*printf("error in ee[nd], nd=%d ee[nd]=%f\n", nd, ee[nd]);*/
                          error = 1;
                          goto end;
                       }
                 
                       if(ee[r] >= 0.0 && ee[r] < 1.0)
                       {
                          p_inf = (fabs(log(ee[r]))<1e-6)? (-log(ee[r])) : (1.0 - ee[r]);
                          factor1 = -ee[r] / p_inf;
                          factor2 = factor1 / p_inf;
                       
                          day_L = cum_e * p_inf;
                          log_day_L_b = cum_loge_b + factor1 * log_ee_b[r];
                          log_day_L_b_b = cum_loge_b_b + factor1 * log_ee_b_b[r] + factor2 * log_ee_b[r] * log_ee_b[r];

                          temp = day_L * pr;
                          L += temp;
                          L_b += log_day_L_b * temp;
                          L_b_b += (log_day_L_b * log_day_L_b + log_day_L_b_b) * temp;
                       }
                          
                       cum_e *= ee[r];
                       cum_loge_b += log_ee_b[r];
                       cum_loge_b_b += log_ee_b_b[r];
                    }
              
                    /* escape during days from inf1 to inf2 */
                    day_L = cum_e;
                    log_day_L_b = cum_loge_b;
                    log_day_L_b_b = cum_loge_b_b;
                    
                    L += day_L;
                    L_b += log_day_L_b * day_L;
                    L_b_b += (log_day_L_b * log_day_L_b + log_day_L_b_b) * day_L;
                    if(L <= 0.0)
                    {
                       error = 1;
                       goto end;
                    }
                    temp = 1.0 / L;
                    log_L = log(L);
                    log_L_b = L_b * temp;
                    log_L_b_b = L_b_b * temp - log_L_b * log_L_b;

                    r = inf1 - start_day; 
                    if(r > 0)
                    {   
                       log_L += cum_logee[r-1];
                       log_L_b += cum_logee_b[r-1];
                       log_L_b_b += cum_logee_b_b[r-1];
                    }
                 }
              }
           
              /*************************************************************************/
           
              score_b += log_L_b;
              info_b_b += log_L_b_b;

              log_L_all += log_L;
              /*
              printf("person=%d  start=%d  stop=%d  day_ill=%d  symptom=%d  log_L=%f  score_b=%f  \n",
                      person->id, start_day, stop_day, person->day_ill, person->symptom, log_L, log_L_b);
              */        
           } /*end of if(person->idx == 0)*/
        }  /*end of loop over community members */
     } // end of loop over communities
     /* use logit transformation. newton-raphson algorithm is suitable for searching over the whole real line */
     b_lb = 0.0;
     b_lb_lb = 0.0;
     
     if(b > 0.0 && b < 1.0)
     {
        lb = logit(b);
        b_lb = b * (1.0 - b);
        b_lb_lb = (1.0 - 2.0 * b) * b_lb;
     }   
      
     score_lb = score_b * b_lb;
     info_lb_lb = info_b_b * b_lb * b_lb + score_b * b_lb_lb;


     info = -info_lb_lb;
     score = score_lb;  
     d = score / info;

     lb += d;
     b_old = b;
     
     b = inv_logit(lb);
     
     if(b < close_to_0)  b = close_to_0;
     if(b > close_to_1)  b = close_to_1;
     
     dd = fabs((b - b_old)/b_old);
     
     /*printf("b=%f\n", p);*/

  } /* end of while */
  
  if(loop >= max_iter)
  {
     error = 1;
     goto end;
  }
  
  (* log_likelihood) = log_L_all;
  est[0] = b;

end:

  free(ee);
  free(log_ee_b);
  free(log_ee_b_b);

  free(cum_logee);
  free(cum_logee_b);
  free(cum_logee_b_b);
        
  return(error);
}



  /********************************
  Get point estimates
  *********************************/
  
int estimation_cs(double *est, double *logL, int n_ini)  
{
  int i, j, k, l, count;
  int error = 0;
  int overall_error = 0;
  double b, lb;
  double temp_logL, NR_logL;
  double *logL_store;
  double min_b, max_b;
  double *est_store;
  
  make_1d_array_double(&est_store, n_ini, 0);
  make_1d_array_double(&logL_store, n_ini, 0);
  
  /* first start the Newton-Raphson emthod from randomly generated initial estimates,
     and select optimal result from these*/
  min_b = 1e-5; max_b = 0.2;
  
  count = 0;
  for(i=0; i<n_ini; i++)
  {
     
     lb = logit(min_b) + (logit(max_b) - logit(min_b)) * runiform(&seed);
  
     b = inv_logit(lb);
     est[0] = b;
     error = newton_raphson_cs(est, &temp_logL);
     if(error == 0)
     {
        est_store[count] = est[0];
        logL_store[count] = temp_logL;
        count++;
     }
  }
  
  /* count is the number of valid Newton-Raphson search */
    
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
     est[0] = est_store[k];
     /*
     printf("MLES found by Newton-Raphson algorithm\n");
     printf("Log likelihood = %e\n", NR_logL);
     printf("b = %e\n", est[0]);
     */
     (* logL) = NR_logL;
  }
  else  overall_error = 1;   
  
  
end:
  free(logL_store);
  free(est_store);  
  
  return(overall_error);
}  

