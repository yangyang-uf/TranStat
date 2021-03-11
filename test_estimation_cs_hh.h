

/* this program is for estimating MLE for the model with
within-household person-to-person transmission probability p1
and common source probability of infection (CSPI) b. */

int newton_raphson_cs_hh(RISK_INDIVIDUAL **risk_history, double *est, double *log_likelihood, int debug)
{
  int i, j, h, k, l, m, n, r, t, rr, error=0;
  int n_par = 2;
  int max_iter = 200;
  int start_day, stop_day;
  double factor1, factor2, factor3;
  int len, inf1, inf2;
  int loop, id, index, skip;
  int type_contact, found;
  double cutpoint1, cutpoint2, dd1, dd2, L, log_L, log_L_all, s;
  double e, f;
  double b, p;
  double b_old, p_old;
  double lb, lp; /* logit tansformation of b, p and q*/
  double b_lb, p_lp;
  double b_lb_lb, p_lp_lp;
  double f1, f2;
  double log_f_b, log_f_p;
  double log_f_b_b, log_f_p_p;
  double log_e_b, log_e_p;
  double log_e_b_b, log_e_p_p;
  double log_L_b, log_L_p;
  double log_L_b_b, log_L_b_p, log_L_p_p;
  double score_lb, score_lp;
  double info_lb_lb, info_lb_lp, info_lp_lp;
  double *ini_ee, *ee, *log_ee_b, *log_ee_p;
  double *log_ee_b_b, *log_ee_p_p;
  double *cum_logee, *cum_logee_b, *cum_logee_p;
  double *cum_logee_b_b, *cum_logee_p_p;
  double cum_e, cum_loge_b, cum_loge_p;
  double cum_loge_b_b, cum_loge_p_p;
  double temp, temp_b, temp_p;
  double temp_b_b, temp_b_p, temp_p_p;
  double L_b, L_p;
  double L_b_b, L_b_p, L_p_p;
  double day_L;
  double log_day_L_b, log_day_L_p;
  double log_day_L_b_b, log_day_L_b_p, log_day_L_p_p;
  double score_b, score_p;
  double info_b_b, info_b_p, info_p_p;
  double pr;
  double p_inf, p_esc;
  MATRIX info, inv, score, d;
  PEOPLE *person, *member;
  CONTACT *ptr_contact;
  FILE *file;
  
  initialize_matrix(&score);
  initialize_matrix(&info);
  initialize_matrix(&inv);
  initialize_matrix(&d);

  inflate_matrix(&info, n_par, n_par, 0);
  inflate_matrix(&inv, n_par, n_par, 0);
  inflate_matrix(&score, n_par, 1, 0);
  inflate_matrix(&d, n_par, 1, 0);
 
  max_epi_duration = 0;
  for(h=0; h<n_community; h++) max_epi_duration = max(max_epi_duration, community[h].epi_duration);
  
  ee = (double *)malloc((size_t) (max_epi_duration * sizeof(double)));
  log_ee_b = (double *)malloc((size_t) (max_epi_duration * sizeof(double)));
  log_ee_p = (double *)malloc((size_t) (max_epi_duration * sizeof(double)));
  log_ee_b_b = (double *)malloc((size_t) (max_epi_duration * sizeof(double)));
  log_ee_p_p = (double *)malloc((size_t) (max_epi_duration * sizeof(double)));
  cum_logee = (double *)malloc((size_t) (max_epi_duration * sizeof(double)));
  cum_logee_b = (double *)malloc((size_t) (max_epi_duration * sizeof(double)));
  cum_logee_p = (double *)malloc((size_t) (max_epi_duration * sizeof(double)));
  cum_logee_b_b = (double *)malloc((size_t) (max_epi_duration * sizeof(double)));
  cum_logee_p_p = (double *)malloc((size_t) (max_epi_duration * sizeof(double)));
           
  ini_ee = (double *)malloc((size_t) (max_epi_duration * sizeof(double)));
  /* p is tran prob in house, b is tran prob from comon source*/

  b = est[0]; p = est[1]; 
  cutpoint1 = cutpoint2 = 1e-10;
  loop = 0;
  dd1 = dd2 = 1.0;
  
  while( (dd1 > cutpoint1 || dd2 > cutpoint2) && (loop < max_iter))
  {
     loop ++;
     //printf("loop=%d\n", loop);
     score_b = score_p = 0.0;
     info_b_b = info_b_p = info_p_p = 0.0;
     log_L_all = 0.0;
     
     for(h=0; h<n_community; h++)
     {
        //printf("h=%d\n", h);
        if(cfg_pars.adjust_for_left_truncation == 1)
            start_day = community[h].latest_idx_day_ill - cfg_pars.min_incubation + 1;
        else  start_day = community[h].day_epi_start;
        stop_day = community[h].day_epi_stop;
        for(t=start_day; t<=stop_day; t++)  
        {
           r = t - start_day;
           ini_ee[r] = 1.0;
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

           ini_ee[r] = f;
           log_ee_b[r] = log_f_b;
           log_ee_b_b[r] = log_f_b_b;
        }     
     
        for(i=0; i<community[h].size; i++)
        {
           person = people + community[h].member[i];
           //printf("person=%d\n", person->id);
           if(person->idx == 0)
           {
              if(cfg_pars.adjust_for_left_truncation == 1)
                 start_day = community[h].latest_idx_day_ill - cfg_pars.min_incubation + 1;
              else  start_day = community[h].day_epi_start;
              stop_day = (person->symptom == 1)? (person->day_ill - cfg_pars.min_incubation):community[h].day_epi_stop;
              //if(person->id == 461)  printf("start_day=%d  stop_day=%d\n", start_day, stop_day);                               
              for(t=start_day; t<=stop_day; t++)
              {
                 e = 1.0;
                 log_e_p = log_e_p_p = 0.0;
                 r = t - start_day;
                 rr = t - community[h].day_epi_start; // both risk history and contact history start from day_epi_start of the community
                 //if(person->id == 461) printf("t=%d  risk_history size=%d\n", t, risk_history[h][rr].size);
                 if(risk_history[h][rr].size > 0)
                 for(j=0; j<risk_history[h][rr].size; j++)
                 {
                    member = people + risk_history[h][rr].infective_subject[j];
                    //printf("member %d  %d\n", risk_history[h][rr].infective_subject[j], member->id);
                    if(member != person)
                    {
                       /*decide the prob that member k is infectious today*/
                       s = risk_history[h][rr].infective_prob[j];
                    
                       if(s > 0.0)
                       {
                          /* find out whether there was a contact*/
                          ptr_contact = person->contact_history[rr].p2p_contact;
                          while(ptr_contact != NULL)
                          {
                             if(ptr_contact->contact_id == member->id)  break;
                             ptr_contact = ptr_contact->next;
                          }   
                          
                          if(ptr_contact != NULL) /* if it's a close contact */
                          {

                             f = 1.0;
                             log_f_p = 0.0;
                             log_f_p_p = 0.0;
                             f = 1.0 - p * s;
                             if(f <= 0.0 || f > 1.0)
                             {
                                if(debug == 1)  printf("Full model: f=%e  p=%e  s=%e\n", f, p, s);
                                error = 1;
                                goto end;
                             }
                             log_f_p = -s/f;
                             log_f_p_p = -log_f_p * log_f_p;

                             e *= f;
                             log_e_p += log_f_p;
                             log_e_p_p += log_f_p_p;
                          }
                       }
                    }    /* end of if (member != person && s > 0) */

                 } /* end loop over all infected subjects */
              
                 ee[r] = ini_ee[r] * e;
                 log_ee_p[r] = log_e_p;
                 log_ee_p_p[r] = log_e_p_p;

              } /* end of day r */
          
              temp = 0.0;
              temp_b = temp_p = temp_b_b = temp_b_p = temp_p_p = 0.0;

              L = 0.0;
              L_b = L_p = L_b_b = L_b_p = L_p_p = 0.0;

              log_L = 0.0;
              log_L_b = log_L_p = log_L_b_b = log_L_b_p = log_L_p_p = 0.0;

              for(t=start_day; t<=stop_day; t++)
              {
                 r = t - start_day;
                 temp += log(ee[r]);
                 cum_logee[r] = temp;

                 temp_b += log_ee_b[r];
                 cum_logee_b[r] = temp_b;

                 temp_p += log_ee_p[r];
                 cum_logee_p[r] = temp_p;

                 temp_b_b += log_ee_b_b[r];
                 cum_logee_b_b[r] = temp_b_b;

                 temp_p_p += log_ee_p_p[r];
                 cum_logee_p_p[r] = temp_p_p;
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
                 cum_loge_b = cum_loge_b_b =
                 cum_loge_p = cum_loge_p_p = 0.0; 
                 for(t=inf1; t<=inf2; t++)
                 {
                    r = t - start_day;
                    pr = pdf_incubation[inf2-t];
                    p_inf = (fabs(log(ee[r]))<1e-6)? (-log(ee[r])) : (1.0 - ee[r]);
                    if(ee[r] < 0.0 || ee[r] > 1.0)
                    {
                       if(debug == 1)  printf("first error in ee[r], t=%d   ee[r]=%f  id=%d  day_ill=%d  inf1=%d  inf2=%d\n", t, ee[r], person->id, person->day_ill, inf1, inf2);
                       error = 1;
                       goto end;
                    }
                 
                    if(ee[r] >= 0.0 && ee[r] < 1.0)
                    {
                       factor1 = -ee[r] / p_inf;
                       factor2 = factor1 / p_inf;
                       
                       day_L = cum_e * p_inf;
                       log_day_L_b = cum_loge_b + factor1 * log_ee_b[r];
                       log_day_L_p = cum_loge_p + factor1 * log_ee_p[r];
                       log_day_L_b_b = cum_loge_b_b + factor1 * log_ee_b_b[r] + factor2 * log_ee_b[r] * log_ee_b[r];
                       log_day_L_b_p = factor2 * log_ee_b[r] * log_ee_p[r];
                       log_day_L_p_p = cum_loge_p_p + factor1 * log_ee_p_p[r] + factor2 * log_ee_p[r] * log_ee_p[r];

                       temp = day_L * pr; 
                       L += temp;
                       L_b += log_day_L_b * temp;
                       L_p += log_day_L_p * temp;
                       L_b_b += (log_day_L_b * log_day_L_b + log_day_L_b_b) * temp;
                       L_b_p += (log_day_L_b * log_day_L_p + log_day_L_b_p) * temp;
                       L_p_p += (log_day_L_p * log_day_L_p + log_day_L_p_p) * temp;
                    }
                           
                    cum_e *= ee[r];
                    cum_loge_b += log_ee_b[r];
                    cum_loge_p += log_ee_p[r];
                    cum_loge_b_b += log_ee_b_b[r];
                    cum_loge_p_p += log_ee_p_p[r];
                 }
                 if(L <= 0.0)
                 {
                    if(debug == 1)
                    {
                       printf("Likelihood must be positive (1) !  L=%f\n", L);
                       printf("id=%d  start_day=%d  day_ill=%d  inf1=%d  inf2=%d  idx_day_ill=%d\n", person->id, start_day, person->day_ill, inf1, inf2, community[h].latest_idx_day_ill);
                      
                       printf("b=%e  p=%e\n", b, p);
                       for(t=inf1;t<=inf2;t++)
                       {
                          r = t - start_day;
                          printf("t=%d  ee[r]=%e 1-ee[r]=%e\n", t, ee[r], 1.0-ee[r]);
                          rr = t - community[h].day_epi_start;
                          for(j=0; j<risk_history[h][rr].size; j++)
                          {
                              member = people + risk_history[h][rr].infective_subject[j];
                              printf("inf_id=%d  inf_prob=%f\n", risk_history[h][rr].infective_subject[j], risk_history[h][rr].infective_prob[j]);
                          }
                       }
                    }   
                    error = 1;
                    goto end;
                 }
                 temp = 1.0 / L;
                 log_L = log(L);
                 log_L_b = L_b * temp;
                 log_L_p = L_p * temp;
                 log_L_b_b = L_b_b * temp - log_L_b * log_L_b;
                 log_L_b_p = L_b_p * temp - log_L_b * log_L_p;
                 log_L_p_p = L_p_p * temp - log_L_p * log_L_p;

                 r = inf1 - start_day; 
                 if(r > 0)
                 {   
                    log_L += cum_logee[r-1];
                    log_L_b += cum_logee_b[r-1];
                    log_L_p += cum_logee_p[r-1];
                    log_L_b_b += cum_logee_b_b[r-1];
                    log_L_p_p += cum_logee_p_p[r-1];
                 }
              }
              else  /* if never been ill */
              {
                 if(cfg_pars.adjust_for_right_censoring == 0)
                 {
                    r = stop_day - start_day;
              
                    log_L = cum_logee[r];
                    log_L_b = cum_logee_b[r];
                    log_L_p = cum_logee_p[r];
                    log_L_b_b = cum_logee_b_b[r];
                    log_L_p_p = cum_logee_p_p[r];
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
                    cum_loge_b = cum_loge_b_b =
                    cum_loge_p = cum_loge_p_p = 0.0; 
                    for(t=inf1; t<=inf2; t++)
                    {
                       pr = sdf_incubation[inf2-t];
                       r = t - start_day;
                       if(ee[r] < 0.0 || ee[r] > 1.0)
                       {
                          if(debug == 1)  printf("second error: t=%d ee[r]=%f id=%d  day_ill=%d  inf1=%d  inf2=%d\n", t, ee[r], person->id, person->day_ill, inf1, inf2);
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
                          log_day_L_p = cum_loge_p + factor1 * log_ee_p[r];
                          log_day_L_b_b = cum_loge_b_b + factor1 * log_ee_b_b[r] + factor2 * log_ee_b[r] * log_ee_b[r];
                          log_day_L_b_p = factor2 * log_ee_b[r] * log_ee_p[r];
                          log_day_L_p_p = cum_loge_p_p + factor1 * log_ee_p_p[r] + factor2 * log_ee_p[r] * log_ee_p[r];

                          temp = day_L * pr;
                          L += temp;
                          L_b += log_day_L_b * temp;
                          L_p += log_day_L_p * temp;
                          L_b_b += (log_day_L_b * log_day_L_b + log_day_L_b_b) * temp;
                          L_b_p += (log_day_L_b * log_day_L_p + log_day_L_b_p) * temp;
                          L_p_p += (log_day_L_p * log_day_L_p + log_day_L_p_p) * temp;
                       }
                          
                       cum_e *= ee[r];
                       cum_loge_b += log_ee_b[r];
                       cum_loge_p += log_ee_p[r];
                       cum_loge_b_b += log_ee_b_b[r];
                       cum_loge_p_p += log_ee_p_p[r];
                    }
              
                    /* escape during days from inf1 to inf2 */
                    day_L = cum_e;
                    log_day_L_b = cum_loge_b;
                    log_day_L_p = cum_loge_p;
                    log_day_L_b_b = cum_loge_b_b;
                    log_day_L_b_p = 0.0;
                    log_day_L_p_p = cum_loge_p_p;
                    
                    L += day_L;
                    L_b += log_day_L_b * day_L;
                    L_p += log_day_L_p * day_L;
                    L_b_b += (log_day_L_b * log_day_L_b + log_day_L_b_b) * day_L;
                    L_b_p += (log_day_L_b * log_day_L_p + log_day_L_b_p) * day_L;
                    L_p_p += (log_day_L_p * log_day_L_p + log_day_L_p_p) * day_L;
                    if(L <= 0.0)
                    {
                       if(debug == 1)
                       {
                          printf("Likelihood must be positive (2) !  L=%f\n", L);
                          printf("id=%d  start_day=%d  day_ill=%d  inf1=%d  inf2=%d\n", i, start_day, person->day_ill, inf1, inf2);
                          printf("b=%e  p=%e\n", b, p);
                          for(t=inf1;t<=inf2;t++)
                          {
                             r = t - start_day;
                             printf("ee[r]=%e 1-ee[r]=%e\n", ee[r], 1.0-ee[r]);
                          }
                       }   
                       error = 1;
                       goto end;
                    }
                    temp = 1.0 / L;
                    log_L = log(L);
                    log_L_b = L_b * temp;
                    log_L_p = L_p * temp;
                    log_L_b_b = L_b_b * temp - log_L_b * log_L_b;
                    log_L_b_p = L_b_p * temp - log_L_b * log_L_p;
                    log_L_p_p = L_p_p * temp - log_L_p * log_L_p;

                    r = inf1 - start_day; 
                    if(r > 0)
                    {   
                       log_L += cum_logee[r-1];
                       log_L_b += cum_logee_b[r-1];
                       log_L_p += cum_logee_p[r-1];
                       log_L_b_b += cum_logee_b_b[r-1];
                       log_L_p_p += cum_logee_p_p[r-1];
                    }
                 }
              }
              /*************************************************************************/
           
              score_b += log_L_b;
              score_p += log_L_p;
              info_b_b += log_L_b_b;
              info_b_p += log_L_b_p;
              info_p_p += log_L_p_p;

              log_L_all += log_L;
              /*
              printf("person=%d  start=%d  stop=%d  day_ill=%d  symptom=%d  log_L=%f  score_b=%f  score_p1=%f\n",
                      person->id, start_day, stop_day, person->day_ill, person->symptom, log_L, log_L_b, log_L_p1);
              */        
           } /*end of if(person->idx == 0)*/
        }  /*end of loop over community members */
     } // end of loop over communities
     /* use logit transformation. newton-raphson algorithm is suitable for searching over the whole real line */
     b_lb = p_lp = 0.0;
     b_lb_lb = p_lp_lp = 0.0;
     
     if(b > 0.0 && b < 1.0)
     {
        lb = logit(b);
        b_lb = b * (1.0 - b);
        b_lb_lb = (1.0 - 2.0 * b) * b_lb;
     }  
     else
     {
        printf("b can not be estimated in the test: %e\n", b);
        error = 1;
        goto end;
     } 
     if(p > 0.0 && p < 1.0)
     {
        lp = logit(p);
        p_lp = p * (1.0 - p);
        p_lp_lp = (1.0 - 2.0 * p) * p_lp;
     }
     else
     {
        printf("p can not be estimated in the test: %e\n", p);
        error = 1;
        goto end;
     } 
      
     score_lb = score_b * b_lb;
     info_lb_lb = info_b_b * b_lb * b_lb + score_b * b_lb_lb;

     score_lp = score_p * p_lp; 
     info_lp_lp = info_p_p * p_lp * p_lp + score_p * p_lp_lp;


     info_lb_lp = info_b_p * b_lb * p_lp;

     info.data[0][0] = -info_lb_lb;
     info.data[0][1] = -info_lb_lp;
     info.data[1][0] = info.data[0][1];
     info.data[1][1] = -info_lp_lp;
        
     score.data[0][0] = score_lb;  
     score.data[1][0] = score_lp;

     temp = info.data[0][0] * info.data[1][1] - info.data[0][1] * info.data[1][0];
     if(fabs(temp) < 1e-100)
     {
        if(debug == 1) printf("Can not invert information matrix, temp=%e\n", temp);
        error = 1;
        goto end;
     }           

     inv.data[0][0] = info.data[1][1] / temp;
     inv.data[0][1] = - info.data[0][1] / temp;
     inv.data[1][0] = - info.data[1][0] / temp;
     inv.data[1][1] = info.data[0][0] / temp;
     
     product(&inv,&score,&d);

     lb += d.data[0][0];
     lp += d.data[1][0];
     
     b_old = b;
     p_old = p;
     
     b = inv_logit(lb);
     p = inv_logit(lp);
     
     if(b < close_to_0)  b = close_to_0;
     if(b > close_to_1)  b = close_to_1;
     if(p < close_to_0)  p = close_to_0;
     if(p > close_to_1)  p = close_to_1;
     
     dd1 = fabs((b - b_old)/b_old);
     dd2 = fabs((p - p_old)/p_old);
     
  } /* end of while */
  
  if(debug == 1)  
  {
      fmprintf(&info);
      fmprintf(&inv);
      printf("score: %e  %e\n", score_lb, score_lp);
      printf("loop=%d  b=%f  p=%f\n", loop, b, p);
  }
  if(loop >= max_iter)
  {
     error = 1;
     goto end;
  }
  
  (* log_likelihood) = log_L_all;
  est[0] = b; est[1] = p;

end:

  free(ini_ee);
  free(ee);
  free(log_ee_b);
  free(log_ee_p);
  free(log_ee_b_b);
  free(log_ee_p_p);

  free(cum_logee);
  free(cum_logee_b);
  free(cum_logee_p);
  free(cum_logee_b_b);
  free(cum_logee_p_p);
        
  
  deflate_matrix(&d);
  deflate_matrix(&info);
  deflate_matrix(&score);
  deflate_matrix(&inv);
  
  return(error);
}



  /********************************
  Get point estimates
  *********************************/
  
int estimation_cs_hh(double *est, double *logL, int n_ini, int debug)  
{
  int h, i, j, k, l, r, t, count;
  int n_par = 2;
  int error = 0;
  int overall_error = 0;
  double b, p, s;
  double lb, lp;
  double temp_logL, NR_logL;
  double *logL_store;
  double min_b, min_p;
  double max_b, max_p;
  MATRIX est_store;
  RISK_INDIVIDUAL ** risk_history;
  PEOPLE *person;

  initialize_matrix(&est_store);
  inflate_matrix(&est_store, n_ini, n_par, 0);
  logL_store = (double *) malloc((size_t) n_ini * sizeof(double));

  risk_history = (RISK_INDIVIDUAL **) malloc((size_t) n_community * sizeof(RISK_INDIVIDUAL *));
  for(h=0; h<n_community; h++)
  {
     risk_history[h] = (RISK_INDIVIDUAL *) malloc((size_t) community[h].epi_duration * sizeof(RISK_INDIVIDUAL));
  
     for(i=0; i<community[h].epi_duration; i++)  
     {
        risk_history[h][i].size = 0;
        risk_history[h][i].infective_subject = (int *) malloc((size_t) community[h].size * sizeof(int));
        risk_history[h][i].infective_prob = (double *) malloc((size_t) community[h].size * sizeof(double));
     }   
     for(i=0; i<community[h].size; i++)
     {
        person = people + community[h].member[i];
        if(person->symptom == 1)
        for(j=cfg_pars.lower_infectious; j<=cfg_pars.upper_infectious; j++)
        {
           t = person->day_ill + j;
           if(t <= community[h].day_epi_stop && t <= person->day_exit)
           {
              r = t - community[h].day_epi_start;
              l = risk_history[h][r].size;
              risk_history[h][r].infective_subject[l] = community[h].member[i];
              s = cfg_pars.prob_infectious[j - cfg_pars.lower_infectious];
              risk_history[h][r].infective_prob[l] = s;
              risk_history[h][r].size++;
           }   
        }
     }              
  }
  
  /* first start the Newton-Raphson emthod from randomly generated initial estimates,
     and select optimal result from these*/
  min_b = 1e-5; max_b = 0.2;
  min_p = 1e-5; max_p = 0.2;
  
  count = 0;
  for(i=0; i<n_ini; i++)
  {
     lb = logit(min_b) + (logit(max_b) - logit(min_b)) * runiform(&seed);
     lp = logit(min_p) + (logit(max_p) - logit(min_p)) * runiform(&seed);
  
     b = inv_logit(lb);
     p = inv_logit(lp);
     est[0] = b; est[1] = p;
     if(debug == 1)  printf("\n\n initial:  b=%e  p=%e\n", est[0], est[1]);
     error = newton_raphson_cs_hh(risk_history, est, &temp_logL, debug);
     if(debug == 1) printf("outcome:  error=%d  b=%e  p=%e  logL=%e\n\n\n", error, est[0], est[1], temp_logL);
     if(error == 0)
     {
        for(j=0; j<n_par; j++)
           est_store.data[count][j] = est[j];
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
     for(j=0; j<n_par; j++)  est[j] = est_store.data[k][j];
     
     //printf("MLES found by Newton-Raphson algorithm\n");
     //printf("Log likelihood = %e\n", NR_logL);
     //printf("b = %e,  p = %e\n", est[0], est[1]);
     
     (* logL) = NR_logL;
  }
  else  overall_error = 1;
     
  
  
end:
  free(logL_store);
  deflate_matrix(&est_store);  
  for(h=0; h<n_community; h++)
  {
     for(i=0; i<community[h].epi_duration; i++)
     {
        free(risk_history[h][i].infective_subject);
        free(risk_history[h][i].infective_prob);
     } 
     free(risk_history[h]);
  }        
  free(risk_history);  
  risk_history = NULL; 
  
  return(overall_error);
}  

