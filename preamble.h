
#define INITIALIZE_LOG_E \
           e = 1.0;\
           for(j=0; j<n_b_mode; j++)\
           {\
              log_e_lb[j] = log_e_lb_lb[j] = 0.0;\
              for(k=0; k<n_c2p_covariate; k++)\
                 log_e_lb_c2p[j][k] = 0.0;\
           }\
           for(k=0; k<n_c2p_covariate; k++)\
           {\
              log_e_c2p[k] = 0.0;\
              for(l=k; l<n_c2p_covariate; l++)\
                 log_e_c2p_c2p[k][l] = 0.0;\
           }\
           for(j=0; j<n_p_mode; j++)\
           {\
              log_e_lp[j] = log_e_lp_lp[j] = 0.0;\
              for(k=0; k<n_p2p_covariate; k++)\
                 log_e_lp_p2p[j][k] = 0.0;\
           }\
           for(k=0; k<n_p2p_covariate; k++)\
           {\
              log_e_p2p[k] = 0.0;\
              for(l=k; l<n_p2p_covariate; l++)\
              log_e_p2p_p2p[k][l] = 0.0;\
           }\

#define UPDATE_LOG_E_C2P \
           if(cfg_pars.common_contact_history_within_community == 1)\
              ptr_risk = ptr_class->risk_history[rr].c2p_risk;\
           else\
              ptr_risk = person->risk_history[rr].c2p_risk;\
           while(ptr_risk != NULL)\
           {\
              covariate_effect = 0.0;\
              for(k=0; k<n_c2p_covariate; k++)\
                 covariate_effect += ptr_risk->covariate[k] * coeff_c2p[k];\
              covariate_effect += ptr_risk->offset;\
              \
              b_mode = ptr_risk->contact_mode;\
              logit_f = lb[b_mode] + covariate_effect;\
              \
              ff = inv_logit(logit_f);\
              f = 1.0 - ff;\
              log_f_lb[b_mode] = - ff;\
              log_f_lb_lb[b_mode] = -f * ff;\
              for(k=0; k<n_c2p_covariate; k++)\
                 log_f_lb_c2p[b_mode][k] = -ptr_risk->covariate[k] * f * ff;\
              \
              for(k=0; k<n_c2p_covariate; k++)\
              {\
                 log_f_c2p[k] = -ptr_risk->covariate[k] * ff;\
                 for(l=k; l<n_c2p_covariate; l++)\
                    log_f_c2p_c2p[k][l] = -ptr_risk->covariate[k] * ptr_risk->covariate[l] * f * ff;\
              }\
              e *= ipow(f, ptr_risk->size);\
              log_e_lb[b_mode] += ptr_risk->size * log_f_lb[b_mode];\
              log_e_lb_lb[b_mode] += ptr_risk->size * log_f_lb_lb[b_mode];\
              for(k=0; k<n_c2p_covariate; k++)\
              {\
                 log_e_c2p[k] += ptr_risk->size * log_f_c2p[k];\
                 log_e_lb_c2p[b_mode][k] += ptr_risk->size * log_f_lb_c2p[b_mode][k];\
                 for(l=k; l<n_c2p_covariate; l++)\
                    log_e_c2p_c2p[k][l] += ptr_risk->size * log_f_c2p_c2p[k][l];\
              }\
              ptr_risk = ptr_risk->next;\
           }

#define UPDATE_LOG_E_P2P \
           if(cfg_pars.common_contact_history_within_community == 1)\
              ptr_risk = ptr_class->risk_history[rr].p2p_risk;\
           else\
              ptr_risk = person->risk_history[rr].p2p_risk;\
           while(ptr_risk != NULL)\
           {\
              covariate_effect = 0.0;\
              for(k=0; k<n_p2p_covariate; k++)\
                 covariate_effect += ptr_risk->covariate[k] * coeff_p2p[k];\
              covariate_effect += ptr_risk->offset;\
              \
              p_mode = ptr_risk->contact_mode;\
              logit_f = lp[p_mode] + covariate_effect;\
              s = ptr_risk->infective_prob;\
              if( s > 0)\
              {\
                 ff = inv_logit(logit_f) * s * bipow(asym_effect, 1 - ptr_risk->symptom);\
                 f = 1 - ff;\
                 factor1 = ff / s;\
                 factor2 = (1 - (2 - ff) * factor1 ) / f;\
                 log_f_lp[p_mode] = - ff * (1 - factor1) / f;\
                 for(k=0; k<n_p2p_covariate; k++)\
                    log_f_p2p[k] = log_f_lp[p_mode] * ptr_risk->covariate[k];\
                 \
                 log_f_lp_lp[p_mode] = factor2 * log_f_lp[p_mode];\
                 \
                 for(k=0; k<n_p2p_covariate; k++)\
                 {\
                    log_f_lp_p2p[p_mode][k] =  factor2 * log_f_p2p[k];\
                    for(l=k; l<n_p2p_covariate; l++)\
                       log_f_p2p_p2p[k][l] = ptr_risk->covariate[k] * factor2 * log_f_p2p[l];\
                 }\
                 if(f <= 0.0 || f > 1.0)\
                 {\
                    printf("Full model: f=%e\n", f);\
                    error = 1;\
                    goto end;\
                 }\
                 e *= ipow(f, ptr_risk->size);\
                 log_e_lp[p_mode] += ptr_risk->size * log_f_lp[p_mode];\
                 log_e_lp_lp[p_mode] += ptr_risk->size * log_f_lp_lp[p_mode];\
                 for(k=0; k<n_p2p_covariate; k++)\
                 {\
                    log_e_p2p[k] += ptr_risk->size * log_f_p2p[k];\
                    log_e_lp_p2p[p_mode][k] += ptr_risk->size * log_f_lp_p2p[p_mode][k];\
                    for(l=k; l<n_p2p_covariate; l++)\
                       log_e_p2p_p2p[k][l] += ptr_risk->size * log_f_p2p_p2p[k][l];\
                 }\
              }\
              ptr_risk = ptr_risk->next;\
           }\

#define UPDATE_LOG_EE \
           ee[r] = e;\
           for(j=0; j<n_b_mode; j++)\
           {\
              log_ee_lb[r][j] = log_e_lb[j];\
              log_ee_lb_lb[r][j] = log_e_lb_lb[j];\
              for(k=0; k<n_c2p_covariate; k++)\
                 log_ee_lb_c2p[r][j][k] = log_e_lb_c2p[j][k];\
           }\
           \
           for(j=0; j<n_p_mode; j++)\
           {\
              log_ee_lp[r][j] = log_e_lp[j];\
              log_ee_lp_lp[r][j] = log_e_lp_lp[j];\
              for(k=0; k<n_p2p_covariate; k++)\
                 log_ee_lp_p2p[r][j][k] = log_e_lp_p2p[j][k];\
           }\
           \
           for(k=0; k<n_c2p_covariate; k++)\
           {\
              log_ee_c2p[r][k] = log_e_c2p[k];\
              for(l=k; l<n_c2p_covariate; l++)\
                 log_ee_c2p_c2p[r][k][l] = log_e_c2p_c2p[k][l];\
           }\
           \
           for(k=0; k<n_p2p_covariate; k++)\
           {\
              log_ee_p2p[r][k] = log_e_p2p[k];\
              for(l=k; l<n_p2p_covariate; l++)\
                 log_ee_p2p_p2p[r][k][l] = log_e_p2p_p2p[k][l];\
           }\


#define INITIALIZE_L \
        L = 0.0;\
        for(j=0; j<n_b_mode; j++)\
        {\
           L_lb[j] = 0.0;\
           for(k=j; k<n_b_mode; k++)          L_lb_lb[j][k] = 0.0;\
           for(k=0; k<n_c2p_covariate; k++)   L_lb_c2p[j][k] = 0.0;\
           for(k=0; k<n_p_mode; k++)          L_lb_lp[j][k] = 0.0;\
           for(k=0; k<n_p2p_covariate; k++)   L_lb_p2p[j][k] = 0.0;\
        }\
        for(j=0; j<n_p_mode; j++)\
        {\
           L_lp[j] = 0.0;\
           for(k=j; k<n_p_mode; k++)          L_lp_lp[j][k] = 0.0;\
           for(k=0; k<n_p2p_covariate; k++)   L_lp_p2p[j][k] = 0.0;\
           for(k=0; k<n_c2p_covariate; k++)   L_lp_c2p[j][k] = 0.0;\
        }\
        for(k=0; k<n_c2p_covariate; k++)\
        {\
           L_c2p[k] = 0.0;\
           for(l=k; l<n_c2p_covariate; l++)   L_c2p_c2p[k][l] = 0.0;\
           for(l=0; l<n_p2p_covariate; l++)   L_c2p_p2p[k][l] = 0.0;\
        }\
        for(k=0; k<n_p2p_covariate; k++)\
        {\
           L_p2p[k] = 0.0;\
           for(l=k; l<n_p2p_covariate; l++)   L_p2p_p2p[k][l] = 0.0;\
        }\
        if(n_u_mode > 0)\
        {\
           L_lu = L_lu_lu = 0.0;\
           for(j=0; j<n_b_mode; j++)  L_lb_lu[j] = 0.0;\
           for(j=0; j<n_p_mode; j++)  L_lp_lu[j] = 0.0;\
           for(k=0; k<n_c2p_covariate; k++) L_lu_c2p[k] = 0.0;\
           for(k=0; k<n_p2p_covariate; k++) L_lu_p2p[k] = 0.0;\
           \
           for(k=0; k<n_pat_covariate; k++)\
           {\
              L_pat[k] = 0.0;\
              L_lu_pat[k] = 0.0;\
              for(l=k; l<n_pat_covariate; l++)  L_pat_pat[k][l] = 0.0;\
              \
              for(j=0; j<n_b_mode; j++)  L_lb_pat[j][k] = 0.0;\
              for(j=0; j<n_p_mode; j++)  L_lp_pat[j][k] = 0.0;\
              for(j=0; j<n_c2p_covariate; j++)  L_c2p_pat[j][k] = 0.0;\
              for(j=0; j<n_p2p_covariate; j++)  L_p2p_pat[j][k] = 0.0;\
           }\
        }


              
#define INITIALIZE_LOG_L \
        log_L = 0.0;\
        for(j=0; j<n_b_mode; j++)\
        {\
           log_L_lb[j] = 0.0;\
           for(k=j; k<n_b_mode; k++)          log_L_lb_lb[j][k] = 0.0;\
           for(k=0; k<n_c2p_covariate; k++)   log_L_lb_c2p[j][k] = 0.0;\
           for(k=0; k<n_p_mode; k++)          log_L_lb_lp[j][k] = 0.0;\
           for(k=0; k<n_p2p_covariate; k++)   log_L_lb_p2p[j][k] = 0.0;\
        }\
        for(j=0; j<n_p_mode; j++)\
        {\
           log_L_lp[j] = 0.0;\
           for(k=j; k<n_p_mode; k++)          log_L_lp_lp[j][k] = 0.0;\
           for(k=0; k<n_c2p_covariate; k++)   log_L_lp_c2p[j][k] = 0.0;\
           for(k=0; k<n_p2p_covariate; k++)   log_L_lp_p2p[j][k] = 0.0;\
        }\
        for(k=0; k<n_c2p_covariate; k++)\
        {\
           log_L_c2p[k] = 0.0;\
           for(l=k; l<n_c2p_covariate; l++)   log_L_c2p_c2p[k][l] = 0.0;\
           for(l=0; l<n_p2p_covariate; l++)   log_L_c2p_p2p[k][l] = 0.0;\
        }\
        for(k=0; k<n_p2p_covariate; k++)\
        {\
           log_L_p2p[k] = 0.0;\
           for(l=k; l<n_p2p_covariate; l++)   log_L_p2p_p2p[k][l] = 0.0;\
        }\
        if(n_u_mode > 0)\
        {\
           log_L_lu = log_L_lu_lu = 0.0;\
           for(j=0; j<n_b_mode; j++)  log_L_lb_lu[j] = 0.0;\
           for(j=0; j<n_p_mode; j++)  log_L_lp_lu[j] = 0.0;\
           for(k=0; k<n_c2p_covariate; k++) log_L_lu_c2p[k] = 0.0;\
           for(k=0; k<n_p2p_covariate; k++) log_L_lu_p2p[k] = 0.0;\
           \
           for(k=0; k<n_pat_covariate; k++)\
           {\
              log_L_pat[k] = 0.0;\
              log_L_lu_pat[k] = 0.0;\
              for(l=k; l<n_pat_covariate; l++)  log_L_pat_pat[k][l] = 0.0;\
              \
              for(j=0; j<n_b_mode; j++)  log_L_lb_pat[j][k] = 0.0;\
              for(j=0; j<n_p_mode; j++)  log_L_lp_pat[j][k] = 0.0;\
              for(j=0; j<n_c2p_covariate; j++)  log_L_c2p_pat[j][k] = 0.0;\
              for(j=0; j<n_p2p_covariate; j++)  log_L_p2p_pat[j][k] = 0.0;\
           }\
        }\





#define INITIALIZE_CUM_LOG_E \
   cum_e = 1.0; \
   for(j=0; j<n_b_mode; j++) \
   { \
      cum_log_e_lb[j] = cum_log_e_lb_lb[j] = 0.0; \
      for(k=0; k<n_c2p_covariate; k++)   cum_log_e_lb_c2p[j][k] = 0.0; \
   } \
   for(j=0; j<n_p_mode; j++) \
   { \
      cum_log_e_lp[j] = cum_log_e_lp_lp[j] = 0.0; \
      for(k=0; k<n_p2p_covariate; k++)   cum_log_e_lp_p2p[j][k] = 0.0; \
   }  \
   for(k=0; k<n_c2p_covariate; k++) \
   { \
      cum_log_e_c2p[k] = 0.0; \
      for(l=k; l<n_c2p_covariate; l++)   cum_log_e_c2p_c2p[k][l] = 0.0; \
   } \
   for(k=0; k<n_p2p_covariate; k++) \
   { \
      cum_log_e_p2p[k] = 0.0; \
      for(l=k; l<n_p2p_covariate; l++)   cum_log_e_p2p_p2p[k][l] = 0.0; \
   } 

#define UPDATE_LOG_DAY_L_FOR_INFECTION \
      p_inf = (fabs(log(ee[r]))<1e-6)? (-log(ee[r])) : (1.0 - ee[r]);\
      factor1 = -ee[r] / p_inf; \
      factor2 = factor1 / p_inf; \
      day_L = cum_e * p_inf; \
      for(j=0; j<n_b_mode; j++) \
      { \
         log_day_L_lb[j] = cum_log_e_lb[j] + factor1 * log_ee_lb[r][j]; \
         for(k=j; k<n_b_mode; k++) \
         { \
            if(j == k) \
               log_day_L_lb_lb[j][k] = cum_log_e_lb_lb[j] + factor1 * log_ee_lb_lb[r][j] \
                                  + factor2 * log_ee_lb[r][j] * log_ee_lb[r][k]; \
            else \
               log_day_L_lb_lb[j][k] = factor2 * log_ee_lb[r][j] * log_ee_lb[r][k]; \
         } \
         for(k=0; k<n_p_mode; k++) \
            log_day_L_lb_lp[j][k] =  factor2 * log_ee_lb[r][j] * log_ee_lp[r][k]; \
         for(k=0; k<n_c2p_covariate; k++) \
            log_day_L_lb_c2p[j][k] = cum_log_e_lb_c2p[j][k] + factor1 * log_ee_lb_c2p[r][j][k] \
                                   + factor2 * log_ee_lb[r][j] * log_ee_c2p[r][k]; \
         for(k=0; k<n_p2p_covariate; k++) \
            log_day_L_lb_p2p[j][k] = factor2 * log_ee_lb[r][j] * log_ee_p2p[r][k]; \
      }\
      for(j=0; j<n_p_mode; j++) \
      { \
         log_day_L_lp[j] = cum_log_e_lp[j] + factor1 * log_ee_lp[r][j]; \
         for(k=j; k<n_p_mode; k++) \
         { \
            if(j == k) \
               log_day_L_lp_lp[j][k] = cum_log_e_lp_lp[j] + factor1 * log_ee_lp_lp[r][j]\
                                     + factor2 * log_ee_lp[r][j] * log_ee_lp[r][k]; \
            else \
               log_day_L_lp_lp[j][k] = factor2 * log_ee_lp[r][j] * log_ee_lp[r][k]; \
         }\
         for(k=0; k<n_c2p_covariate; k++) \
            log_day_L_lp_c2p[j][k] = factor2 * log_ee_lp[r][j] * log_ee_c2p[r][k]; \
         for(k=0; k<n_p2p_covariate; k++) \
            log_day_L_lp_p2p[j][k] = cum_log_e_lp_p2p[j][k] + factor1 * log_ee_lp_p2p[r][j][k] \
                                   + factor2 * log_ee_lp[r][j] * log_ee_p2p[r][k]; \
      } \
      for(k=0; k<n_c2p_covariate; k++) \
      { \
         log_day_L_c2p[k] = cum_log_e_c2p[k] + factor1 * log_ee_c2p[r][k]; \
         for(l=k; l<n_c2p_covariate; l++) \
            log_day_L_c2p_c2p[k][l] = cum_log_e_c2p_c2p[k][l] + factor1 * log_ee_c2p_c2p[r][k][l] \
                                    + factor2 * log_ee_c2p[r][k] * log_ee_c2p[r][l]; \
         for(l=0; l<n_p2p_covariate; l++) \
            log_day_L_c2p_p2p[k][l] = factor2 * log_ee_c2p[r][k] * log_ee_p2p[r][l]; \
      } \
      for(k=0; k<n_p2p_covariate; k++) \
      { \
         log_day_L_p2p[k] = cum_log_e_p2p[k] + factor1 * log_ee_p2p[r][k]; \
         for(l=k; l<n_p2p_covariate; l++) \
            log_day_L_p2p_p2p[k][l] = cum_log_e_p2p_p2p[k][l] + factor1 * log_ee_p2p_p2p[r][k][l] \
                                    + factor2 * log_ee_p2p[r][k] * log_ee_p2p[r][l]; \
      }
                       
#define UPDATE_L_FOR_INFECTION \
      temp = day_L * pr; \
      L += temp; \
      for(j=0; j<n_b_mode; j++) \
      { \
         L_lb[j] += log_day_L_lb[j] * temp; \
         for(k=j; k<n_b_mode; k++) \
            L_lb_lb[j][k] += (log_day_L_lb[j] * log_day_L_lb[k] + log_day_L_lb_lb[j][k]) * temp; \
         for(k=0; k<n_p_mode; k++) \
            L_lb_lp[j][k] += (log_day_L_lb[j] * log_day_L_lp[k] + log_day_L_lb_lp[j][k]) * temp; \
         for(k=0; k<n_c2p_covariate; k++) \
            L_lb_c2p[j][k] +=  (log_day_L_lb[j] * log_day_L_c2p[k] + log_day_L_lb_c2p[j][k]) * temp; \
         for(k=0; k<n_p2p_covariate; k++) \
            L_lb_p2p[j][k] +=  (log_day_L_lb[j] * log_day_L_p2p[k] + log_day_L_lb_p2p[j][k]) * temp; \
      } \
      for(j=0; j<n_p_mode; j++) \
      { \
         L_lp[j] += log_day_L_lp[j] * temp; \
         for(k=j; k<n_p_mode; k++) \
            L_lp_lp[j][k] += (log_day_L_lp[j] * log_day_L_lp[k] + log_day_L_lp_lp[j][k]) * temp; \
         for(k=0; k<n_c2p_covariate; k++) \
            L_lp_c2p[j][k] +=  (log_day_L_lp[j] * log_day_L_c2p[k] + log_day_L_lp_c2p[j][k]) * temp; \
         for(k=0; k<n_p2p_covariate; k++) \
            L_lp_p2p[j][k] +=  (log_day_L_lp[j] * log_day_L_p2p[k] + log_day_L_lp_p2p[j][k]) * temp; \
      } \
      for(k=0; k<n_c2p_covariate; k++) \
      { \
         L_c2p[k] += log_day_L_c2p[k] * temp; \
         for(l=k; l<n_c2p_covariate; l++) \
            L_c2p_c2p[k][l] += (log_day_L_c2p[k] * log_day_L_c2p[l] + log_day_L_c2p_c2p[k][l]) * temp; \
         for(l=0; l<n_p2p_covariate; l++) \
            L_c2p_p2p[k][l] += (log_day_L_c2p[k] * log_day_L_p2p[l] + log_day_L_c2p_p2p[k][l]) * temp; \
      } \
      for(k=0; k<n_p2p_covariate; k++) \
      { \
         L_p2p[k] += log_day_L_p2p[k] * temp; \
         for(l=k; l<n_p2p_covariate; l++) \
            L_p2p_p2p[k][l] += (log_day_L_p2p[k] * log_day_L_p2p[l] + log_day_L_p2p_p2p[k][l]) * temp; \
      }\
      if(n_u_mode > 0)\
      {\
         L_lu += day_L * pr_lu;\
         L_lu_lu += day_L * pr_lu_lu;\
         \
         for(k=0; k<n_b_mode; k++)\
            L_lb_lu[k] += log_day_L_lb[k] * day_L * pr_lu;\
         for(k=0; k<n_p_mode; k++)\
            L_lp_lu[k] += log_day_L_lp[k] * day_L * pr_lu;\
         for(k=0; k<n_c2p_covariate; k++)\
            L_lu_c2p[k] += log_day_L_c2p[k] * day_L * pr_lu;\
         for(k=0; k<n_p2p_covariate; k++)\
            L_lu_p2p[k] += log_day_L_p2p[k] * day_L * pr_lu;\
         \
         for(k=0; k<n_pat_covariate; k++)\
         {\
            L_pat[k] += day_L * pr_pat[k];\
            L_lu_pat[k] += day_L * pr_lu_pat[k];\
            for(l=k; l<n_pat_covariate; l++)\
               L_pat_pat[k][l] += day_L * pr_pat_pat[k][l];\
            for(l=0; l<n_b_mode; l++)\
               L_lb_pat[l][k] += log_day_L_lb[l] * day_L * pr_pat[k];\
            for(l=0; l<n_p_mode; l++)\
               L_lp_pat[l][k] += log_day_L_lp[l] * day_L * pr_pat[k];\
            for(l=0; l<n_c2p_covariate; l++)\
               L_c2p_pat[l][k] += log_day_L_c2p[l] * day_L * pr_pat[k];\
            for(l=0; l<n_p2p_covariate; l++)\
               L_p2p_pat[l][k] += log_day_L_p2p[l] * day_L * pr_pat[k];\
         }\
      }


#define UPDATE_CUM_LOG_E \
      cum_e *= ee[r]; \
      for(j=0; j<n_b_mode; j++) \
      { \
         cum_log_e_lb[j] += log_ee_lb[r][j]; \
         cum_log_e_lb_lb[j] += log_ee_lb_lb[r][j]; \
         for(k=0; k<n_c2p_covariate; k++) \
            cum_log_e_lb_c2p[j][k] += log_ee_lb_c2p[r][j][k]; \
      } \
      for(j=0; j<n_p_mode; j++) \
      { \
         cum_log_e_lp[j] += log_ee_lp[r][j]; \
         cum_log_e_lp_lp[j] += log_ee_lp_lp[r][j]; \
         for(k=0; k<n_p2p_covariate; k++) \
            cum_log_e_lp_p2p[j][k] += log_ee_lp_p2p[r][j][k]; \
      } \
      for(k=0; k<n_c2p_covariate; k++) \
      { \
         cum_log_e_c2p[k] += log_ee_c2p[r][k]; \
         for(l=k; l<n_c2p_covariate; l++) \
            cum_log_e_c2p_c2p[k][l] += log_ee_c2p_c2p[r][k][l]; \
      } \
      for(k=0; k<n_p2p_covariate; k++) \
      { \
         cum_log_e_p2p[k] += log_ee_p2p[r][k]; \
         for(l=k; l<n_p2p_covariate; l++) \
            cum_log_e_p2p_p2p[k][l] += log_ee_p2p_p2p[r][k][l]; \
      }      
                    
              
#define UPDATE_LOG_DAY_L_FOR_ESCAPE \
      day_L = cum_e; \
      for(j=0; j<n_b_mode; j++) \
      { \
         log_day_L_lb[j] = cum_log_e_lb[j]; \
         for(k=j; k<n_b_mode; k++) \
         { \
            if(j == k)  log_day_L_lb_lb[j][k] = cum_log_e_lb_lb[j]; \
            else  log_day_L_lb_lb[j][k] = 0.0; \
         } \
         for(k=0; k<n_c2p_covariate; k++) \
            log_day_L_lb_c2p[j][k] = cum_log_e_lb_c2p[j][k]; \
         for(k=0; k<n_p_mode; k++) \
            log_day_L_lb_lp[j][k] = 0.0; \
         for(k=0; k<n_p2p_covariate; k++) \
            log_day_L_lb_p2p[j][k] = 0.0; \
      } \
      for(j=0; j<n_p_mode; j++) \
      { \
         log_day_L_lp[j] = cum_log_e_lp[j]; \
         for(k=j; k<n_p_mode; k++) \
         { \
            if(j == k) log_day_L_lp_lp[j][k] = cum_log_e_lp_lp[j]; \
            else  log_day_L_lp_lp[j][k] = 0.0; \
         } \
         for(k=0; k<n_p2p_covariate; k++) \
            log_day_L_lp_p2p[j][k] = cum_log_e_lp_p2p[j][k]; \
         for(k=0; k<n_c2p_covariate; k++) \
            log_day_L_lp_c2p[j][k] = 0.0; \
      } \
      for(k=0; k<n_c2p_covariate; k++) \
      { \
         log_day_L_c2p[k] = cum_log_e_c2p[k]; \
         for(l=k; l<n_c2p_covariate; l++) \
            log_day_L_c2p_c2p[k][l] = cum_log_e_c2p_c2p[k][l]; \
         for(l=0; l<n_p2p_covariate; l++) \
            log_day_L_c2p_p2p[k][l] = 0.0; \
      } \
      for(k=0; k<n_p2p_covariate; k++) \
      { \
         log_day_L_p2p[k] = cum_log_e_p2p[k]; \
         for(l=k; l<n_p2p_covariate; l++) \
            log_day_L_p2p_p2p[k][l] = cum_log_e_p2p_p2p[k][l]; \
      }

#define UPDATE_L_FOR_ESCAPE \
                 L += day_L; \
                 for(j=0; j<n_b_mode; j++) \
                 { \
                    L_lb[j] += log_day_L_lb[j] * day_L; \
                    for(k=j; k<n_b_mode; k++) \
                       L_lb_lb[j][k] += (log_day_L_lb[j] * log_day_L_lb[k] + log_day_L_lb_lb[j][k]) * day_L; \
                    for(k=0; k<n_p_mode; k++) \
                       L_lb_lp[j][k] += (log_day_L_lb[j] * log_day_L_lp[k] + log_day_L_lb_lp[j][k]) * day_L; \
                    for(k=0; k<n_c2p_covariate; k++) \
                       L_lb_c2p[j][k] +=  (log_day_L_lb[j] * log_day_L_c2p[k] + log_day_L_lb_c2p[j][k]) * day_L; \
                    for(k=0; k<n_p2p_covariate; k++) \
                       L_lb_p2p[j][k] +=  (log_day_L_lb[j] * log_day_L_p2p[k] + log_day_L_lb_p2p[j][k]) * day_L; \
                 } \
                 for(j=0; j<n_p_mode; j++) \
                 { \
                    L_lp[j] += log_day_L_lp[j] * day_L; \
                    for(k=j; k<n_p_mode; k++) \
                       L_lp_lp[j][k] += (log_day_L_lp[j] * log_day_L_lp[k] + log_day_L_lp_lp[j][k]) * day_L; \
                    for(k=0; k<n_c2p_covariate; k++) \
                       L_lp_c2p[j][k] +=  (log_day_L_lp[j] * log_day_L_c2p[k] + log_day_L_lp_c2p[j][k]) * day_L; \
                    for(k=0; k<n_p2p_covariate; k++) \
                       L_lp_p2p[j][k] +=  (log_day_L_lp[j] * log_day_L_p2p[k] + log_day_L_lp_p2p[j][k]) * day_L; \
                 } \
                 for(k=0; k<n_c2p_covariate; k++) \
                 { \
                    L_c2p[k] += log_day_L_c2p[k] * day_L; \
                    for(l=k; l<n_c2p_covariate; l++) \
                       L_c2p_c2p[k][l] += (log_day_L_c2p[k] * log_day_L_c2p[l] + log_day_L_c2p_c2p[k][l]) * day_L; \
                    for(l=0; l<n_p2p_covariate; l++) \
                       L_c2p_p2p[k][l] += (log_day_L_c2p[k] * log_day_L_p2p[l] + log_day_L_c2p_p2p[k][l]) * day_L; \
                 } \
                 for(k=0; k<n_p2p_covariate; k++) \
                 { \
                    L_p2p[k] += log_day_L_p2p[k] * day_L; \
                    for(l=k; l<n_p2p_covariate; l++) \
                       L_p2p_p2p[k][l] += (log_day_L_p2p[k] * log_day_L_p2p[l] + log_day_L_p2p_p2p[k][l]) * day_L; \
                 }

              
#define UPDATE_LOG_L \
                 factor1 = 1.0 / L; \
                 factor2 = factor1 / L; \
                 log_L = log(L); \
                 for(j=0; j<n_b_mode; j++) \
                 { \
                    log_L_lb[j] = L_lb[j] * factor1; \
                    for(k=j; k<n_b_mode; k++) \
                       log_L_lb_lb[j][k] = L_lb_lb[j][k] * factor1 - L_lb[j] * L_lb[k] * factor2; \
                    for(k=0; k<n_p_mode; k++) \
                       log_L_lb_lp[j][k] = L_lb_lp[j][k] * factor1 - L_lb[j] * L_lp[k] * factor2; \
                    for(k=0; k<n_c2p_covariate; k++) \
                       log_L_lb_c2p[j][k] = L_lb_c2p[j][k] * factor1 - L_lb[j] * L_c2p[k] * factor2; \
                    for(k=0; k<n_p2p_covariate; k++) \
                       log_L_lb_p2p[j][k] = L_lb_p2p[j][k] * factor1 - L_lb[j] * L_p2p[k] * factor2; \
                 } \
                 for(j=0; j<n_p_mode; j++) \
                 { \
                    log_L_lp[j] = L_lp[j] * factor1; \
                    for(k=j; k<n_p_mode; k++) \
                       log_L_lp_lp[j][k] = L_lp_lp[j][k] * factor1 - L_lp[j] * L_lp[k] * factor2; \
                    for(k=0; k<n_c2p_covariate; k++) \
                       log_L_lp_c2p[j][k] = L_lp_c2p[j][k] * factor1 - L_lp[j] * L_c2p[k] * factor2; \
                    for(k=0; k<n_p2p_covariate; k++) \
                       log_L_lp_p2p[j][k] = L_lp_p2p[j][k] * factor1 - L_lp[j] * L_p2p[k] * factor2; \
                 } \
                 for(k=0; k<n_c2p_covariate; k++) \
                 { \
                    log_L_c2p[k] = L_c2p[k] * factor1; \
                    for(l=k; l<n_c2p_covariate; l++) \
                       log_L_c2p_c2p[k][l] = L_c2p_c2p[k][l] * factor1 - L_c2p[k] * L_c2p[l] * factor2; \
                    for(l=0; l<n_p2p_covariate; l++) \
                       log_L_c2p_p2p[k][l] = L_c2p_p2p[k][l] * factor1 - L_c2p[k] * L_p2p[l] * factor2; \
                 }  \
                 for(k=0; k<n_p2p_covariate; k++) \
                 { \
                    log_L_p2p[k] = L_p2p[k] * factor1; \
                    for(l=k; l<n_p2p_covariate; l++) \
                       log_L_p2p_p2p[k][l] = L_p2p_p2p[k][l] * factor1 - L_p2p[k] * L_p2p[l] * factor2; \
                 }\
                 if(n_u_mode > 0)\
                 {\
                    log_L_lu = L_lu * factor1;\
                    log_L_lu_lu = L_lu_lu * factor1 - L_lu * L_lu * factor2;\
                    \
                    for(k=0; k<n_b_mode; k++)\
                       log_L_lb_lu[k] = L_lb_lu[k] * factor1 - L_lb[k] * L_lu * factor2;\
                    for(k=0; k<n_p_mode; k++)\
                       log_L_lp_lu[k] = L_lp_lu[k] * factor1 - L_lp[k] * L_lu * factor2;\
                    for(k=0; k<n_c2p_covariate; k++)\
                       log_L_lu_c2p[k] = L_lu_c2p[k] * factor1 - L_lu * L_c2p[k] * factor2;\
                    for(k=0; k<n_p2p_covariate; k++)\
                       log_L_lu_p2p[k] = L_lu_p2p[k] * factor1 - L_lu * L_p2p[k] * factor2;\
                    \
                    for(k=0; k<n_pat_covariate; k++)\
                    {\
                       log_L_pat[k] = L_pat[k] * factor1;\
                       log_L_lu_pat[k] = L_lu_pat[k] * factor1 - L_lu * L_pat[k] * factor2;\
                       for(l=k; l<n_pat_covariate; l++)\
                          log_L_pat_pat[k][l] = L_pat_pat[k][l] * factor1 - L_pat[k] * L_pat[l] * factor2;\
                       \
                       for(l=0; l<n_b_mode; l++)\
                          log_L_lb_pat[l][k] = L_lb_pat[l][k] * factor1 - L_lb[l] * L_pat[k] * factor2;\
                       for(l=0; l<n_p_mode; l++)\
                          log_L_lp_pat[l][k] = L_lp_pat[l][k] * factor1 - L_lp[l] * L_pat[k] * factor2;\
                       for(l=0; l<n_c2p_covariate; l++)\
                          log_L_c2p_pat[l][k] = L_c2p_pat[l][k] * factor1 - L_c2p[l] * L_pat[k] * factor2;\
                       for(l=0; l<n_p2p_covariate; l++)\
                          log_L_p2p_pat[l][k] = L_p2p_pat[l][k] * factor1 - L_p2p[l] * L_pat[k] * factor2;\
                    }\
                 }


#define CONCATENATE_ESCAPE_HISTORY \
                 r = inf1 - start_day; \
                 if(r > 0) \
                 { \
                    log_L += cum_log_ee[r-1]; \
                    for(j=0; j<n_b_mode; j++) \
                    { \
                       log_L_lb[j] += cum_log_ee_lb[r-1][j]; \
                       log_L_lb_lb[j][j] += cum_log_ee_lb_lb[r-1][j]; \
                       for(k=0; k<n_c2p_covariate; k++) \
                          log_L_lb_c2p[j][k] += cum_log_ee_lb_c2p[r-1][j][k]; \
                    }  \
                    for(j=0; j<n_p_mode; j++) \
                    { \
                       log_L_lp[j] += cum_log_ee_lp[r-1][j]; \
                       log_L_lp_lp[j][j] += cum_log_ee_lp_lp[r-1][j]; \
                       for(k=0; k<n_p2p_covariate; k++) \
                          log_L_lp_p2p[j][k] += cum_log_ee_lp_p2p[r-1][j][k]; \
                    }  \
                    for(k=0; k<n_c2p_covariate; k++) \
                    {\
                       log_L_c2p[k] += cum_log_ee_c2p[r-1][k]; \
                       for(l=k; l<n_c2p_covariate; l++) \
                          log_L_c2p_c2p[k][l] += cum_log_ee_c2p_c2p[r-1][k][l]; \
                    }\
                    for(k=0; k<n_p2p_covariate; k++) \
                    {\
                       log_L_p2p[k] += cum_log_ee_p2p[r-1][k]; \
                       for(l=k; l<n_p2p_covariate; l++) \
                          log_L_p2p_p2p[k][l] += cum_log_ee_p2p_p2p[r-1][k][l]; \
                    }\
                 }

#define UPDATE_LQ \
                    L = exp(log_L);\
                    LQ = L * (1 - Q) + Q;\
                    factor = (1 - Q) * L;\
                    for(j=0; j<n_b_mode; j++)\
                    {\
                       LQ_lb[j] = factor * log_L_lb[j];\
                       for(k=j; k<n_b_mode; k++)\
                          LQ_lb_lb[j][k] = factor * (log_L_lb[j] * log_L_lb[k] + log_L_lb_lb[j][k]);\
                       for(k=0; k<n_p_mode; k++)\
                          LQ_lb_lp[j][k] = factor * (log_L_lb[j] * log_L_lp[k] + log_L_lb_lp[j][k]);\
                       for(k=0; k<n_c2p_covariate; k++)\
                          LQ_lb_c2p[j][k] = factor * (log_L_lb[j] * log_L_c2p[k] + log_L_lb_c2p[j][k]);\
                       for(k=0; k<n_p2p_covariate; k++)\
                          LQ_lb_p2p[j][k] = factor * (log_L_lb[j] * log_L_p2p[k] + log_L_lb_p2p[j][k]);\
                    }\
                    \
                    for(j=0; j<n_p_mode; j++)\
                    {\
                       LQ_lp[j] = factor * log_L_lp[j];\
                       for(k=j; k<n_p_mode; k++)\
                          LQ_lp_lp[j][k] = factor * (log_L_lp[j] * log_L_lp[k] + log_L_lp_lp[j][k]);\
                       for(k=0; k<n_c2p_covariate; k++)\
                          LQ_lp_c2p[j][k] = factor * (log_L_lp[j] * log_L_c2p[k] + log_L_lp_c2p[j][k]);\
                       for(k=0; k<n_p2p_covariate; k++)\
                          LQ_lp_p2p[j][k] = factor * (log_L_lp[j] * log_L_p2p[k] + log_L_lp_p2p[j][k]);\
                    }\
                    \
                    for(j=0; j<n_c2p_covariate; j++)\
                    {\
                       LQ_c2p[j] = factor * log_L_c2p[j];\
                       for(k=j; k<n_c2p_covariate; k++)\
                          LQ_c2p_c2p[j][k] = factor * (log_L_c2p[j] * log_L_c2p[k] + log_L_c2p_c2p[j][k]);\
                       for(k=0; k<n_p2p_covariate; k++)\
                          LQ_c2p_p2p[j][k] = factor * (log_L_c2p[j] * log_L_p2p[k] + log_L_c2p_p2p[j][k]);\
                    }\
                    \
                    for(j=0; j<n_p2p_covariate; j++)\
                    {\
                       LQ_p2p[j] = factor * log_L_p2p[j];\
                       for(k=j; k<n_p2p_covariate; k++)\
                          LQ_p2p_p2p[j][k] = factor * (log_L_p2p[j] * log_L_p2p[k] + log_L_p2p_p2p[j][k]);\
                    }\
                    \
                    LQ_lq = -Q_lq * L + Q_lq;\
                    LQ_lq_lq = -Q_lq_lq * L + Q_lq_lq;\
                    for(j=0; j<n_b_mode; j++)\
                       LQ_lb_lq[j] = L * log_L_lb[j] * (-Q_lq);\
                    for(j=0; j<n_p_mode; j++)\
                       LQ_lp_lq[j] = L * log_L_lp[j] * (-Q_lq);\
                    for(j=0; j<n_c2p_covariate; j++)\
                       LQ_lq_c2p[j] = L * log_L_c2p[j] * (-Q_lq);\
                    for(j=0; j<n_p2p_covariate; j++)\
                       LQ_lq_p2p[j] = L * log_L_p2p[j] * (-Q_lq);\
                    for(j=0; j<n_imm_covariate; j++)\
                       LQ_lq_imm[j] = - Q_lq_imm[j] * L + Q_lq_imm[j];\
                    \
                    for(k=0; k<n_imm_covariate; k++)\
                    {\
                       LQ_imm[k] = - Q_imm[k] * L;\
                       for(l=k; l<n_imm_covariate; l++)\
                          LQ_imm_imm[k][l] = - Q_imm_imm[k][l] * L + Q_imm_imm[k][l];\
                       \
                       for(j=0; j<n_b_mode; j++)\
                          LQ_lb_imm[j][k] = L * log_L_lb[j] * (-Q_imm[k]);\
                       for(j=0; j<n_p_mode; j++)\
                          LQ_lp_imm[j][k] = L * log_L_lp[j] * (-Q_imm[k]);\
                       for(j=0; j<n_c2p_covariate; j++)\
                          LQ_c2p_imm[j][k] = L * log_L_c2p[j] * (-Q_imm[k]);\
                       for(j=0; j<n_p2p_covariate; j++)\
                          LQ_p2p_imm[j][k] = L * log_L_p2p[j] * (-Q_imm[k]);\
                    }\
                    \
                    if(n_u_mode > 0)\
                    {\
                       LQ_lu = factor * log_L_lu;\
                       LQ_lu_lu = factor * (log_L_lu * log_L_lu + log_L_lu_lu);\
                       for(j=0; j<n_b_mode; j++)  LQ_lb_lu[j]= factor * (log_L_lb[j] * log_L_lu + log_L_lb_lu[j]);\
                       for(j=0; j<n_p_mode; j++)  LQ_lp_lu[j]= factor * (log_L_lp[j] * log_L_lu + log_L_lp_lu[j]);\
                       for(k=0; k<n_c2p_covariate; k++)\
                          LQ_lu_c2p[k] = factor * (log_L_lu * log_L_c2p[k] + log_L_lu_c2p[k]);\
                       for(k=0; k<n_p2p_covariate; k++)\
                          LQ_lu_p2p[k] = factor * (log_L_lu * log_L_p2p[k] + log_L_lu_p2p[k]);\
                       for(k=0; k<n_pat_covariate; k++)\
                          LQ_lu_pat[k] = factor * (log_L_lu * log_L_pat[k] + log_L_lu_pat[k]);\
                       \
                       for(k=0; k<n_pat_covariate; k++)\
                       {\
                          LQ_pat[k] = factor * log_L_pat[k];\
                          for(l=k; l<n_pat_covariate; l++)\
                             LQ_pat_pat[k][l] = factor * (log_L_pat[k] * log_L_pat[l] + log_L_pat_pat[k][l]);\
                          for(j=0; j<n_b_mode; j++)\
                             LQ_lb_pat[j][k] = factor * (log_L_lb[j] * log_L_pat[k] + log_L_lb_pat[j][k]);\
                          for(j=0; j<n_p_mode; j++)\
                             LQ_lp_pat[j][k] = factor * (log_L_lp[j] * log_L_pat[k] + log_L_lp_pat[j][k]);\
                          for(j=0; j<n_c2p_covariate; j++)\
                             LQ_c2p_pat[j][k] = factor * (log_L_c2p[j] * log_L_pat[k] + log_L_c2p_pat[j][k]);\
                          for(j=0; j<n_p2p_covariate; j++)\
                             LQ_p2p_pat[j][k] = factor * (log_L_p2p[j] * log_L_pat[k] + log_L_p2p_pat[j][k]);\
                       }\
                       \
                       LQ_lu_lq = L * log_L_lu * (-Q_lq);\
                       for(k=0; k<n_imm_covariate; k++)\
                          LQ_lu_imm[k] = L * log_L_lu * (-Q_imm[k]);\
                       for(k=0; k<n_pat_covariate; k++)\
                       {\
                          LQ_lq_pat[k] = L * log_L_pat[k] * (-Q_lq);\
                          for(l=0; l<n_imm_covariate; l++)\
                             LQ_pat_imm[k][l] = L * log_L_pat[k] * (-Q_imm[l]);\
                       }\
                    }


#define UPDATE_LOG_LQ \
                    log_L = log(LQ);\
                    factor1 = 1 / LQ;\
                    factor2 = -factor1/LQ;\
                    for(j=0; j<n_b_mode; j++)\
                    {\
                       log_L_lb[j] = factor1 * LQ_lb[j];\
                       for(k=j; k<n_b_mode; k++)\
                          log_L_lb_lb[j][k] = factor1 * LQ_lb_lb[j][k] + factor2 * LQ_lb[j] * LQ_lb[k];\
                       for(k=0; k<n_p_mode; k++)\
                          log_L_lb_lp[j][k] = factor1 * LQ_lb_lp[j][k] + factor2 * LQ_lb[j] * LQ_lp[k];\
                       for(k=0; k<n_c2p_covariate; k++)\
                          log_L_lb_c2p[j][k] = factor1 * LQ_lb_c2p[j][k] + factor2 * LQ_lb[j] * LQ_c2p[k];\
                       for(k=0; k<n_p2p_covariate; k++)\
                          log_L_lb_p2p[j][k] = factor1 * LQ_lb_p2p[j][k] + factor2 * LQ_lb[j] * LQ_p2p[k];\
                    }\
                    for(j=0; j<n_p_mode; j++)\
                    {\
                       log_L_lp[j] = factor1 * LQ_lp[j];\
                       for(k=j; k<n_p_mode; k++)\
                          log_L_lp_lp[j][k] = factor1 * LQ_lp_lp[j][k] + factor2 * LQ_lp[j] * LQ_lp[k];\
                       for(k=0; k<n_c2p_covariate; k++)\
                          log_L_lp_c2p[j][k] = factor1 * LQ_lp_c2p[j][k] + factor2 * LQ_lp[j] * LQ_c2p[k];\
                       for(k=0; k<n_p2p_covariate; k++)\
                          log_L_lp_p2p[j][k] = factor1 * LQ_lp_p2p[j][k] + factor2 * LQ_lp[j] * LQ_p2p[k];\
                    }\
                    \
                    for(j=0; j<n_c2p_covariate; j++)\
                    {\
                       log_L_c2p[j] = factor1 * LQ_c2p[j];\
                       for(k=j; k<n_c2p_covariate; k++)\
                          log_L_c2p_c2p[j][k] = factor1 * LQ_c2p_c2p[j][k] + factor2 * LQ_c2p[j] * LQ_c2p[k];\
                       for(k=0; k<n_p2p_covariate; k++)\
                          log_L_c2p_p2p[j][k] = factor1 * LQ_c2p_p2p[j][k] + factor2 * LQ_c2p[j] * LQ_p2p[k];\
                    }\
                    \
                    for(j=0; j<n_p2p_covariate; j++)\
                    {\
                       log_L_p2p[j] = factor1 * LQ_p2p[j];\
                       for(k=j; k<n_p2p_covariate; k++)\
                          log_L_p2p_p2p[j][k] = factor1 * LQ_p2p_p2p[j][k] + factor2 * LQ_p2p[j] * LQ_p2p[k];\
                    }\
                    \
                    log_L_lq = factor1 * LQ_lq;\
                    log_L_lq_lq = factor1 * LQ_lq_lq + factor2 * LQ_lq * LQ_lq;\
                    for(j=0; j<n_b_mode; j++)\
                       log_L_lb_lq[j] = factor1 * LQ_lb_lq[j] + factor2 * LQ_lb[j] * LQ_lq;\
                    for(j=0; j<n_p_mode; j++)\
                       log_L_lp_lq[j] = factor1 * LQ_lp_lq[j] + factor2 * LQ_lp[j] * LQ_lq;\
                    for(j=0; j<n_c2p_covariate; j++)\
                       log_L_lq_c2p[j] = factor1 * LQ_lq_c2p[j] + factor2 * LQ_c2p[j] * LQ_lq;\
                    for(j=0; j<n_p2p_covariate; j++)\
                       log_L_lq_p2p[j] = factor1 * LQ_lq_p2p[j] + factor2 * LQ_p2p[j] * LQ_lq;\
                    for(j=0; j<n_imm_covariate; j++)\
                       log_L_lq_imm[j] = factor1 * LQ_lq_imm[j] + factor2 * LQ_lq * LQ_imm[j];\
                    \
                    for(k=0; k<n_imm_covariate; k++)\
                    {\
                       log_L_imm[k] = factor1 * LQ_imm[k];\
                       for(l=k; l<n_imm_covariate; l++)\
                          log_L_imm_imm[k][l] = factor1 * LQ_imm_imm[k][l] + factor2 * LQ_imm[k] * LQ_imm[l];\
                       \
                       for(j=0; j<n_b_mode; j++)\
                          log_L_lb_imm[j][k] = factor1 * LQ_lb_imm[j][k] + factor2 * LQ_lb[j] * LQ_imm[k];\
                       for(j=0; j<n_p_mode; j++)\
                          log_L_lp_imm[j][k] = factor1 * LQ_lp_imm[j][k] + factor2 * LQ_lp[j] * LQ_imm[k];\
                       for(j=0; j<n_c2p_covariate; j++)\
                          log_L_c2p_imm[j][k] = factor1 * LQ_c2p_imm[j][k] + factor2 * LQ_c2p[j] * LQ_imm[k];\
                       for(j=0; j<n_p2p_covariate; j++)\
                          log_L_p2p_imm[j][k] = factor1 * LQ_p2p_imm[j][k] + factor2 * LQ_p2p[j] * LQ_imm[k];\
                    }\
                    if(n_u_mode > 0)\
                    {\
                       log_L_lu = factor1 * LQ_lu;\
                       log_L_lu_lu = factor1 * LQ_lu_lu + factor2 * LQ_lu * LQ_lu;\
                       for(j=0; j<n_b_mode; j++)\
                          log_L_lb_lu[j] = factor1 * LQ_lb_lu[j] + factor2 * LQ_lb[j] * LQ_lu;\
                       for(j=0; j<n_p_mode; j++)\
                          log_L_lp_lu[j] = factor1 * LQ_lp_lu[j] + factor2 * LQ_lp[j] * LQ_lu;\
                       for(k=0; k<n_c2p_covariate; k++)\
                          log_L_lu_c2p[k] = factor1 * LQ_lu_c2p[k] + factor2 * LQ_lu * LQ_c2p[k];\
                       for(k=0; k<n_p2p_covariate; k++)\
                          log_L_lu_p2p[k] = factor1 * LQ_lu_p2p[k] + factor2 * LQ_lu * LQ_p2p[k];\
                       for(k=0; k<n_pat_covariate; k++)\
                          log_L_lu_pat[k] = factor1 * LQ_lu_pat[k] + factor2 * LQ_lu * LQ_pat[k];\
                       \
                       for(k=0; k<n_pat_covariate; k++)\
                       {\
                          log_L_pat[k] = factor1 * LQ_pat[k];\
                          for(l=k; l<n_pat_covariate; l++)\
                             log_L_pat_pat[k][l] = factor1 * LQ_pat_pat[k][l] + factor2 * LQ_pat[k] * LQ_pat[l];\
                          \
                          for(j=0; j<n_b_mode; j++)\
                             log_L_lb_pat[j][k] = factor1 * LQ_lb_pat[j][k] + factor2 * LQ_lb[j] * LQ_pat[k];\
                          for(j=0; j<n_p_mode; j++)\
                             log_L_lp_pat[j][k] = factor1 * LQ_lp_pat[j][k] + factor2 * LQ_lp[j] * LQ_pat[k];\
                          for(j=0; j<n_c2p_covariate; j++)\
                             log_L_c2p_pat[j][k] = factor1 * LQ_c2p_pat[j][k] + factor2 * LQ_c2p[j] * LQ_pat[k];\
                          for(j=0; j<n_p2p_covariate; j++)\
                             log_L_p2p_pat[j][k] = factor1 * LQ_p2p_pat[j][k] + factor2 * LQ_p2p[j] * LQ_pat[k];\
                       }\
                       \
                       log_L_lu_lq = factor1 * LQ_lu_lq + factor2 * LQ_lu * LQ_lq;\
                       for(k=0; k<n_imm_covariate; k++)\
                          log_L_lu_imm[k] = factor1 * LQ_lu_imm[k] + factor2 * LQ_lu * LQ_imm[k];\
                       for(k=0; k<n_pat_covariate; k++)\
                       {\
                          log_L_lq_pat[k] = factor1 * LQ_lq_pat[k] + factor2 * LQ_lq * LQ_pat[k];\
                          for(l=0; l<n_imm_covariate; l++)\
                             log_L_pat_imm[k][l] = factor1 * LQ_pat_imm[k][l] + factor2 * LQ_pat[k] * LQ_imm[l];\
                       }\
                    }



