void print_community(int h)
{
   int i, j;
   printf("community %d  size_idx=%d  size=%d\n", h, community[h].size_idx, community[h].size);
   for(j=0; j<community[h].size; j++)
   {
     i = community[h].member[j];
     printf("member %d: idx=%d  infection=%d  symptom=%d  day_ill=%d\n", i, people[i].idx, people[i].infection, people[i].symptom, people[i].day_ill);
   }
   if(community[h].size_idx > 0)
   {
      printf("Index cases:");
      for(j=0; j<community[h].size_idx; j++)
      {
        i = community[h].idx[j];
        printf("%d ", i);
      }
      printf("\n");
   }
}

void count_cases(void)
{
   int i, n_ign, n_imm, n_idx, n_sym, n_asym, n_esc;
   n_ign = n_imm = n_idx = n_sym = n_asym = n_esc = 0;
   for(i=0; i < p_size; i++)
   {
      if(people[i].ignore == 1)  n_ign++;
      else
      {
         if(people[i].pre_immune == 0)
         {
            if(people[i].idx == 0)
            {
               if(people[i].infection == 1)
               {
                  if(people[i].symptom == 1)  n_sym++;
                  else  n_asym++;
               }
               else n_esc++;
            }
            else  n_idx++;
         }
         else  n_imm++;
      }
   }      
   printf("Ignore:%d  Preimm:%d  Index:%d, Secondary: n_sym=%d  n_asym=%d, Escape:%d\n", n_ign, n_imm, n_idx, n_sym, n_asym, n_esc);
}

void create_arrays(void)
{
  int n_par, n_b_mode, n_p_mode, n_q_mode, n_u_mode;
  int n_c2p_covariate, n_p2p_covariate, n_imm_covariate, n_pat_covariate;
  int n_covariate, n_time_ind_covariate, n_time_dep_covariate;
  int n_sus_p2p_covariate, n_inf_p2p_covariate, n_int_p2p_covariate;
  int n_par_equiclass;

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
  n_imm_covariate = cfg_pars.n_imm_covariate;
  n_pat_covariate = cfg_pars.n_pat_covariate;
  n_par = cfg_pars.n_par;
  n_par_equiclass = cfg_pars.n_par_equiclass;

  ee = cum_ee = cum_log_ee = NULL;
  lb = lp = lq = lu = coeff_c2p = coeff_p2p = coeff_imm = coeff_pat = NULL;
  b = p = q = u = OR_c2p = OR_p2p = OR_pat = OR_imm = CPI = SAR = R0 = NULL;
  
  se_lb = se_lp = se_lu = se_lq = se_coeff_c2p = se_coeff_p2p = se_coeff_pat = se_coeff_imm = NULL;
  se_b = se_p = se_u = se_q = se_OR_c2p = se_OR_p2p = se_OR_pat = se_OR_imm = se_CPI = se_SAR = se_R0 = NULL;

  lower_b = lower_p = lower_u = lower_q = lower_OR_c2p = lower_OR_p2p = lower_OR_pat = lower_OR_imm = lower_CPI = lower_SAR = lower_R0 = NULL;
  upper_b = upper_p = upper_u = upper_q = upper_OR_c2p = upper_OR_p2p = upper_OR_pat = upper_OR_imm = upper_CPI = upper_SAR = upper_R0 = NULL;

  log_f_lb = log_f_c2p = log_f_lp = log_f_p2p = NULL;
  log_f_lb_lb = log_f_lp_lp = NULL;
  log_f_lb_c2p = log_f_lp_p2p = log_f_c2p_c2p = log_f_p2p_p2p = NULL;

  log_e_lb = log_e_c2p = log_e_lp = log_e_p2p = NULL;
  log_e_lb_lb = log_e_lp_lp = NULL;
  log_e_lb_c2p = log_e_lp_p2p = log_e_c2p_c2p = log_e_p2p_p2p = NULL;

  log_ee_lb = log_ee_c2p = log_ee_lp = log_ee_p2p = NULL;
  log_ee_lb_lb = log_ee_lp_lp = NULL;
  log_ee_lb_c2p = log_ee_c2p_c2p = log_ee_lp_p2p = log_ee_p2p_p2p = NULL;

  cum_log_ee_lb = cum_log_ee_c2p = cum_log_ee_lp = cum_log_ee_p2p = NULL;
  cum_log_ee_lb_lb = cum_log_ee_lp_lp = NULL;
  cum_log_ee_lb_c2p = cum_log_ee_c2p_c2p = cum_log_ee_lp_p2p = cum_log_ee_p2p_p2p = NULL;

  cum_log_e_lb = cum_log_e_c2p = cum_log_e_lp = cum_log_e_p2p = NULL;
  cum_log_e_lb_lb = cum_log_e_lp_lp = NULL;
  cum_log_e_lb_c2p = cum_log_e_c2p_c2p = cum_log_e_lp_p2p = cum_log_e_p2p_p2p = NULL;

  temp_lb = temp_c2p = temp_lp = temp_p2p = NULL;
  temp_lb_lb = temp_lp_lp = NULL;
  temp_lb_c2p = temp_c2p_c2p = temp_lp_p2p = temp_p2p_p2p = NULL;

  log_day_L_lb = log_day_L_c2p = log_day_L_lp = log_day_L_p2p = NULL;
  log_day_L_lb_lb = log_day_L_lb_c2p = log_day_L_lb_lp = log_day_L_lb_p2p = NULL;
  log_day_L_c2p_c2p = log_day_L_c2p_p2p = log_day_L_lp_lp = NULL;
  log_day_L_lp_p2p = log_day_L_lp_c2p = log_day_L_p2p_p2p = NULL;

  L_lb = L_lp = L_c2p = L_p2p = NULL;
  L_lb_lb = L_lb_lp = L_lb_c2p = L_lb_p2p =
  L_lp_lp = L_lp_c2p = L_lp_p2p =
  L_c2p_c2p = L_c2p_p2p = L_p2p_p2p = NULL;
  //L_lu = L_lq = L_lu_lu = L_lu_lq = L_lq_lq = NULL;
  L_lb_lu = L_lp_lu = L_lu_c2p = L_lu_p2p = L_lu_imm = L_lu_pat = NULL;
  L_lb_lq = L_lp_lq = L_lq_c2p = L_lq_p2p = L_lq_imm = L_lq_pat = NULL;
  L_pat = L_imm = NULL;
  L_lb_pat = L_lp_pat = L_c2p_pat = L_p2p_pat =
  L_lb_imm = L_lp_imm = L_c2p_imm = L_p2p_imm =
  L_pat_pat = L_pat_imm = L_imm_imm = NULL;

  LQ_lb = LQ_lp = LQ_c2p = LQ_p2p = NULL;
  LQ_lb_lb = LQ_lb_lp = LQ_lb_c2p = LQ_lb_p2p =
  LQ_lp_lp = LQ_lp_c2p = LQ_lp_p2p =
  LQ_c2p_c2p = LQ_c2p_p2p = LQ_p2p_p2p = NULL;
  //LQ_lu = LQ_lq = LQ_lu_lu = LQ_lu_lq = LQ_lq_lq = NULL;
  LQ_lb_lu = LQ_lp_lu = LQ_lu_c2p = LQ_lu_p2p = LQ_lu_imm = LQ_lu_pat = NULL;
  LQ_lb_lq = LQ_lp_lq = LQ_lq_c2p = LQ_lq_p2p = LQ_lq_imm = LQ_lq_pat = NULL;
  LQ_pat = LQ_imm = NULL;
  LQ_lb_pat = LQ_lp_pat = LQ_c2p_pat = LQ_p2p_pat =
  LQ_lb_imm = LQ_lp_imm = LQ_c2p_imm = LQ_p2p_imm =
  LQ_pat_pat = LQ_pat_imm = LQ_imm_imm = NULL;

  log_L_lb = log_L_lp = log_L_c2p = log_L_p2p = NULL;
  log_L_lb_lb = log_L_lb_lp = log_L_lb_c2p = log_L_lb_p2p =
  log_L_lp_lp = log_L_lp_c2p = log_L_lp_p2p =
  log_L_c2p_c2p = log_L_c2p_p2p = log_L_p2p_p2p = NULL;
  //log_L_lu = log_L_lq = log_L_lu_lu = log_L_lu_lq = log_L_lq_lq = NULL;
  log_L_lb_lu = log_L_lp_lu = log_L_lu_c2p = log_L_lu_p2p = log_L_lu_imm = log_L_lu_pat = NULL;
  log_L_lb_lq = log_L_lp_lq = log_L_lq_c2p = log_L_lq_p2p = log_L_lq_imm = log_L_lq_pat = NULL;
  log_L_pat = log_L_imm = NULL;
  log_L_lb_pat = log_L_lp_pat = log_L_c2p_pat = log_L_p2p_pat =
  log_L_lb_imm = log_L_lp_imm = log_L_c2p_imm = log_L_p2p_imm =
  log_L_pat_pat = log_L_pat_imm = log_L_imm_imm = NULL;

  my_log_L_lb = my_log_L_lp = my_log_L_c2p = my_log_L_p2p = NULL;
  my_log_L_lb_lb = my_log_L_lb_lp = my_log_L_lb_c2p = my_log_L_lb_p2p =
  my_log_L_lp_lp = my_log_L_lp_c2p = my_log_L_lp_p2p =
  my_log_L_c2p_c2p = my_log_L_c2p_p2p = my_log_L_p2p_p2p = NULL;
  //my_log_L_lu = my_log_L_lq = my_log_L_lu_lu = my_log_L_lu_lq = my_log_L_lq_lq = NULL;
  my_log_L_lb_lu = my_log_L_lp_lu = my_log_L_lu_c2p = my_log_L_lu_p2p = my_log_L_lu_imm = my_log_L_lu_pat = NULL;
  my_log_L_lb_lq = my_log_L_lp_lq = my_log_L_lq_c2p = my_log_L_lq_p2p = my_log_L_lq_imm = my_log_L_lq_pat = NULL;
  my_log_L_pat = my_log_L_imm = NULL;
  my_log_L_lb_pat = my_log_L_lp_pat = my_log_L_c2p_pat = my_log_L_p2p_pat =
  my_log_L_lb_imm = my_log_L_lp_imm = my_log_L_c2p_imm = my_log_L_p2p_imm =
  my_log_L_pat_pat = my_log_L_pat_imm = my_log_L_imm_imm = NULL;

  pr_pat = pr_lu_pat = NULL; pr_pat_pat = NULL;
  U_pat = U_lu_pat = NULL; U_pat_pat = NULL;
  Q_imm = Q_lq_imm = NULL; Q_imm_imm = NULL;

  score_lb = score_lp = score_lu = score_lq = score_c2p = score_p2p = score_pat = score_imm = NULL;
  info_lb_lb = info_lb_lp = info_lb_c2p = info_lb_p2p =
  info_lp_lp = info_lp_c2p = info_lp_p2p =
  info_c2p_c2p = info_c2p_p2p = info_p2p_p2p =
  info_lu_lu = info_lq_lq = info_lu_lq =
  info_lb_lu = info_lp_lu = info_lu_c2p = info_lu_p2p = info_lu_pat = info_lu_imm =
  info_lb_lq = info_lp_lq = info_lq_c2p = info_lq_p2p = info_lq_pat = info_lq_imm =
  info_lb_pat = info_lp_pat = info_c2p_pat = info_p2p_pat =
  info_lb_imm = info_lp_imm = info_c2p_imm = info_p2p_imm =
  info_pat_pat = info_pat_imm = info_imm_imm = NULL;

  if(n_b_mode > 0)
  {
     make_1d_array_double(&b, n_b_mode, 0.0);
     make_1d_array_double(&lb, n_b_mode, 0.0);
     make_1d_array_double(&se_b, n_b_mode, 0.0);
     make_1d_array_double(&se_lb, n_b_mode, 0.0);
     make_1d_array_double(&lower_b, n_b_mode, 0.0);
     make_1d_array_double(&upper_b, n_b_mode, 0.0);
     
     make_1d_array_double(&CPI, n_b_mode, 0.0);
     make_1d_array_double(&se_CPI, n_b_mode, 0.0);
     make_1d_array_double(&lower_CPI, n_b_mode, 0.0);
     make_1d_array_double(&upper_CPI, n_b_mode, 0.0);
  }
  if(n_p_mode > 0)
  {
     make_1d_array_double(&p, n_p_mode, 0.0);
     make_1d_array_double(&lp, n_p_mode, 0.0);
     make_1d_array_double(&se_p, n_p_mode, 0.0);
     make_1d_array_double(&se_lp, n_p_mode, 0.0);
     make_1d_array_double(&lower_p, n_p_mode, 0.0);
     make_1d_array_double(&upper_p, n_p_mode, 0.0);

     make_1d_array_double(&SAR, n_p_mode, 0.0);
     make_1d_array_double(&se_SAR, n_p_mode, 0.0);
     make_1d_array_double(&lower_SAR, n_p_mode, 0.0);
     make_1d_array_double(&upper_SAR, n_p_mode, 0.0);
     
     if(cfg_pars.n_R0_multiplier > 0)
     {
        make_1d_array_double(&R0, cfg_pars.n_R0_multiplier, 0.0);
        make_1d_array_double(&se_R0, cfg_pars.n_R0_multiplier, 0.0);
        make_1d_array_double(&lower_R0, cfg_pars.n_R0_multiplier, 0.0);
        make_1d_array_double(&upper_R0, cfg_pars.n_R0_multiplier, 0.0);
     }    
  }
  if(n_u_mode > 0)
  {
     make_1d_array_double(&u, n_u_mode, 0.0);
     make_1d_array_double(&lu, n_u_mode, 0.0);
     make_1d_array_double(&se_u, n_u_mode, 0.0);
     make_1d_array_double(&se_lu, n_u_mode, 0.0);
     make_1d_array_double(&lower_u, n_u_mode, 0.0);
     make_1d_array_double(&upper_u, n_u_mode, 0.0);
  }
  if(n_q_mode > 0)
  {
     make_1d_array_double(&q, n_q_mode, 0.0);
     make_1d_array_double(&lq, n_q_mode, 0.0);
     make_1d_array_double(&se_q, n_q_mode, 0.0);
     make_1d_array_double(&se_lq, n_q_mode, 0.0);
     make_1d_array_double(&lower_q, n_q_mode, 0.0);
     make_1d_array_double(&upper_q, n_q_mode, 0.0);
  }
  if(n_c2p_covariate > 0)  
  {
     make_1d_array_double(&coeff_c2p, n_c2p_covariate, 0.0);
     make_1d_array_double(&OR_c2p, n_c2p_covariate, 0.0);
     make_1d_array_double(&se_coeff_c2p, n_c2p_covariate, 0.0);
     make_1d_array_double(&se_OR_c2p, n_c2p_covariate, 0.0);
     make_1d_array_double(&lower_OR_c2p, n_c2p_covariate, 0.0);
     make_1d_array_double(&upper_OR_c2p, n_c2p_covariate, 0.0);
  }
  if(n_p2p_covariate > 0)  
  {
     make_1d_array_double(&coeff_p2p, n_p2p_covariate, 0.0);
     make_1d_array_double(&OR_p2p, n_p2p_covariate, 0.0);
     make_1d_array_double(&se_coeff_p2p, n_p2p_covariate, 0.0);
     make_1d_array_double(&se_OR_p2p, n_p2p_covariate, 0.0);
     make_1d_array_double(&lower_OR_p2p, n_p2p_covariate, 0.0);
     make_1d_array_double(&upper_OR_p2p, n_p2p_covariate, 0.0);
  }
  if(n_pat_covariate > 0)  
  {
     make_1d_array_double(&coeff_pat, n_pat_covariate, 0.0);
     make_1d_array_double(&OR_pat, n_pat_covariate, 0.0);
     make_1d_array_double(&se_coeff_pat, n_pat_covariate, 0.0);
     make_1d_array_double(&se_OR_pat, n_pat_covariate, 0.0);
     make_1d_array_double(&lower_OR_pat, n_pat_covariate, 0.0);
     make_1d_array_double(&upper_OR_pat, n_pat_covariate, 0.0);
  }
  if(n_imm_covariate > 0)  
  {
     make_1d_array_double(&coeff_imm, n_imm_covariate, 0.0);
     make_1d_array_double(&OR_imm, n_imm_covariate, 0.0);
     make_1d_array_double(&se_coeff_imm, n_imm_covariate, 0.0);
     make_1d_array_double(&se_OR_imm, n_imm_covariate, 0.0);
     make_1d_array_double(&lower_OR_imm, n_imm_covariate, 0.0);
     make_1d_array_double(&upper_OR_imm, n_imm_covariate, 0.0);
  }

  make_1d_array_double(&ee, max_epi_duration, 1.0);
  make_1d_array_double(&cum_log_ee, max_epi_duration, 0.0);

  if(n_b_mode > 0)
  {
     make_1d_array_double(&log_f_lb, n_b_mode, 0.0);
     make_1d_array_double(&log_f_lb_lb, n_b_mode, 0.0);

     make_1d_array_double(&log_e_lb, n_b_mode, 0.0);
     make_1d_array_double(&log_e_lb_lb, n_b_mode, 0.0);

     make_1d_array_double(&cum_log_e_lb, n_b_mode, 0.0);
     make_1d_array_double(&cum_log_e_lb_lb, n_b_mode, 0.0);

     make_1d_array_double(&temp_lb, n_b_mode, 0.0);
     make_1d_array_double(&temp_lb_lb, n_b_mode, 0.0);

     make_2d_array_double(&log_ee_lb, max_epi_duration, n_b_mode, 0.0);
     make_2d_array_double(&log_ee_lb_lb, max_epi_duration, n_b_mode, 0.0);

     make_2d_array_double(&cum_log_ee_lb, max_epi_duration, n_b_mode, 0.0);
     make_2d_array_double(&cum_log_ee_lb_lb, max_epi_duration, n_b_mode, 0.0);

     make_1d_array_double(&log_day_L_lb, n_b_mode, 0.0);
     make_2d_array_double(&log_day_L_lb_lb, n_b_mode, n_b_mode, 0.0);

     make_1d_array_double(&L_lb, n_b_mode, 0.0);
     make_2d_array_double(&L_lb_lb, n_b_mode, n_b_mode, 0.0);

     make_1d_array_double(&LQ_lb, n_b_mode, 0.0);
     make_2d_array_double(&LQ_lb_lb, n_b_mode, n_b_mode, 0.0);

     make_1d_array_double(&log_L_lb, n_b_mode, 0.0);
     make_2d_array_double(&log_L_lb_lb, n_b_mode, n_b_mode, 0.0);

     make_1d_array_double(&my_log_L_lb, n_b_mode, 0.0);
     make_2d_array_double(&my_log_L_lb_lb, n_b_mode, n_b_mode, 0.0);

     make_1d_array_double(&score_lb, n_b_mode, 0.0);
     make_2d_array_double(&info_lb_lb, n_b_mode, n_b_mode, 0.0);
  }

  if(n_p_mode > 0)
  {
     make_1d_array_double(&log_f_lp, n_p_mode, 0.0);
     make_1d_array_double(&log_f_lp_lp, n_p_mode, 0.0);

     make_1d_array_double(&log_e_lp, n_p_mode, 0.0);
     make_1d_array_double(&log_e_lp_lp, n_p_mode, 0.0);

     make_1d_array_double(&cum_log_e_lp, n_p_mode, 0.0);
     make_1d_array_double(&cum_log_e_lp_lp, n_p_mode, 0.0);

     make_1d_array_double(&temp_lp, n_p_mode, 0.0);
     make_1d_array_double(&temp_lp_lp, n_p_mode, 0.0);

     make_2d_array_double(&log_ee_lp, max_epi_duration, n_p_mode, 0.0);
     make_2d_array_double(&log_ee_lp_lp, max_epi_duration, n_p_mode, 0.0);

     make_2d_array_double(&cum_log_ee_lp, max_epi_duration, n_p_mode, 0.0);
     make_2d_array_double(&cum_log_ee_lp_lp, max_epi_duration, n_p_mode, 0.0);

     make_1d_array_double(&log_day_L_lp, n_p_mode, 0.0);
     make_2d_array_double(&log_day_L_lp_lp, n_p_mode, n_p_mode, 0.0);

     make_1d_array_double(&L_lp, n_p_mode, 0.0);
     make_2d_array_double(&L_lp_lp, n_p_mode, n_p_mode, 0.0);

     make_1d_array_double(&LQ_lp, n_p_mode, 0.0);
     make_2d_array_double(&LQ_lp_lp, n_p_mode, n_p_mode, 0.0);

     make_1d_array_double(&log_L_lp, n_p_mode, 0.0);
     make_2d_array_double(&log_L_lp_lp, n_p_mode, n_p_mode, 0.0);

     make_1d_array_double(&my_log_L_lp, n_p_mode, 0.0);
     make_2d_array_double(&my_log_L_lp_lp, n_p_mode, n_p_mode, 0.0);

     make_1d_array_double(&score_lp, n_p_mode, 0.0);
     make_2d_array_double(&info_lp_lp, n_p_mode, n_p_mode, 0.0);
  }

  if(n_u_mode > 0)
  {
     make_1d_array_double(&score_lu, n_u_mode, 0.0);
     make_2d_array_double(&info_lu_lu, n_u_mode, n_u_mode, 0.0);
  }

  if(n_q_mode > 0)
  {
     make_1d_array_double(&score_lq, n_q_mode, 0.0);
     make_2d_array_double(&info_lq_lq, n_q_mode, n_q_mode, 0.0);
  }

  if(n_c2p_covariate > 0)
  {

     make_1d_array_double(&log_f_c2p, n_c2p_covariate, 0.0);
     make_2d_array_double(&log_f_c2p_c2p, n_c2p_covariate, n_c2p_covariate, 0.0);

     make_1d_array_double(&log_e_c2p, n_c2p_covariate, 0.0);
     make_2d_array_double(&log_e_c2p_c2p, n_c2p_covariate, n_c2p_covariate, 0.0);

     make_1d_array_double(&cum_log_e_c2p, n_c2p_covariate, 0.0);
     make_2d_array_double(&cum_log_e_c2p_c2p, n_c2p_covariate, n_c2p_covariate, 0.0);

     make_1d_array_double(&temp_c2p, n_c2p_covariate, 0.0);
     make_2d_array_double(&temp_c2p_c2p, n_c2p_covariate, n_c2p_covariate, 0.0);

     make_2d_array_double(&log_ee_c2p, max_epi_duration, n_c2p_covariate, 0.0);
     make_3d_array_double(&log_ee_c2p_c2p, max_epi_duration, n_c2p_covariate, n_c2p_covariate, 0.0);

     make_2d_array_double(&cum_log_ee_c2p, max_epi_duration, n_c2p_covariate, 0.0);
     make_3d_array_double(&cum_log_ee_c2p_c2p, max_epi_duration, n_c2p_covariate, n_c2p_covariate, 0.0);

     make_1d_array_double(&log_day_L_c2p, n_c2p_covariate, 0.0);
     make_2d_array_double(&log_day_L_c2p_c2p, n_c2p_covariate, n_c2p_covariate, 0.0);

     make_1d_array_double(&L_c2p, n_c2p_covariate, 0.0);
     make_2d_array_double(&L_c2p_c2p, n_c2p_covariate, n_c2p_covariate, 0.0);

     make_1d_array_double(&LQ_c2p, n_c2p_covariate, 0.0);
     make_2d_array_double(&LQ_c2p_c2p, n_c2p_covariate, n_c2p_covariate, 0.0);

     make_1d_array_double(&log_L_c2p, n_c2p_covariate, 0.0);
     make_2d_array_double(&log_L_c2p_c2p, n_c2p_covariate, n_c2p_covariate, 0.0);

     make_1d_array_double(&my_log_L_c2p, n_c2p_covariate, 0.0);
     make_2d_array_double(&my_log_L_c2p_c2p, n_c2p_covariate, n_c2p_covariate, 0.0);

     make_1d_array_double(&score_c2p, n_c2p_covariate, 0.0);
     make_2d_array_double(&info_c2p_c2p, n_c2p_covariate, n_c2p_covariate, 0.0);
  }

  if(n_p2p_covariate > 0)
  {

     make_1d_array_double(&log_f_p2p, n_p2p_covariate, 0.0);
     make_2d_array_double(&log_f_p2p_p2p, n_p2p_covariate, n_p2p_covariate, 0.0);

     make_1d_array_double(&log_e_p2p, n_p2p_covariate, 0.0);
     make_2d_array_double(&log_e_p2p_p2p, n_p2p_covariate, n_p2p_covariate, 0.0);

     make_1d_array_double(&cum_log_e_p2p, n_p2p_covariate, 0.0);
     make_2d_array_double(&cum_log_e_p2p_p2p, n_p2p_covariate, n_p2p_covariate, 0.0);

     make_1d_array_double(&temp_p2p, n_p2p_covariate, 0.0);
     make_2d_array_double(&temp_p2p_p2p, n_p2p_covariate, n_p2p_covariate, 0.0);

     make_2d_array_double(&log_ee_p2p, max_epi_duration, n_p2p_covariate, 0.0);
     make_3d_array_double(&log_ee_p2p_p2p, max_epi_duration, n_p2p_covariate, n_p2p_covariate, 0.0);

     make_2d_array_double(&cum_log_ee_p2p, max_epi_duration, n_p2p_covariate, 0.0);
     make_3d_array_double(&cum_log_ee_p2p_p2p, max_epi_duration, n_p2p_covariate, n_p2p_covariate, 0.0);

     make_1d_array_double(&log_day_L_p2p, n_p2p_covariate, 0.0);
     make_2d_array_double(&log_day_L_p2p_p2p, n_p2p_covariate, n_p2p_covariate, 0.0);

     make_1d_array_double(&L_p2p, n_p2p_covariate, 0.0);
     make_2d_array_double(&L_p2p_p2p, n_p2p_covariate, n_p2p_covariate, 0.0);

     make_1d_array_double(&LQ_p2p, n_p2p_covariate, 0.0);
     make_2d_array_double(&LQ_p2p_p2p, n_p2p_covariate, n_p2p_covariate, 0.0);

     make_1d_array_double(&log_L_p2p, n_p2p_covariate, 0.0);
     make_2d_array_double(&log_L_p2p_p2p, n_p2p_covariate, n_p2p_covariate, 0.0);

     make_1d_array_double(&my_log_L_p2p, n_p2p_covariate, 0.0);
     make_2d_array_double(&my_log_L_p2p_p2p, n_p2p_covariate, n_p2p_covariate, 0.0);

     make_1d_array_double(&score_p2p, n_p2p_covariate, 0.0);
     make_2d_array_double(&info_p2p_p2p, n_p2p_covariate, n_p2p_covariate, 0.0);
  }

  if(n_pat_covariate > 0)
  {
     make_1d_array_double(&pr_pat, n_pat_covariate, 0.0);
     make_2d_array_double(&pr_pat_pat, n_pat_covariate, n_pat_covariate, 0.0);

     make_1d_array_double(&U_pat, n_pat_covariate, 0.0);
     make_2d_array_double(&U_pat_pat, n_pat_covariate, n_pat_covariate, 0.0);

     make_1d_array_double(&L_pat, n_pat_covariate, 0.0);
     make_2d_array_double(&L_pat_pat, n_pat_covariate, n_pat_covariate, 0.0);

     make_1d_array_double(&LQ_pat, n_pat_covariate, 0.0);
     make_2d_array_double(&LQ_pat_pat, n_pat_covariate, n_pat_covariate, 0.0);

     make_1d_array_double(&log_L_pat, n_pat_covariate, 0.0);
     make_2d_array_double(&log_L_pat_pat, n_pat_covariate, n_pat_covariate, 0.0);

     make_1d_array_double(&my_log_L_pat, n_pat_covariate, 0.0);
     make_2d_array_double(&my_log_L_pat_pat, n_pat_covariate, n_pat_covariate, 0.0);

     make_1d_array_double(&score_pat, n_pat_covariate, 0.0);
     make_2d_array_double(&info_pat_pat, n_pat_covariate, n_pat_covariate, 0.0);
  }

  if(n_imm_covariate > 0)
  {
     make_1d_array_double(&Q_imm, n_imm_covariate, 0.0);
     make_2d_array_double(&Q_imm_imm, n_imm_covariate, n_imm_covariate, 0.0);

     make_1d_array_double(&L_imm, n_imm_covariate, 0.0);
     make_2d_array_double(&L_imm_imm, n_imm_covariate, n_imm_covariate, 0.0);

     make_1d_array_double(&LQ_imm, n_imm_covariate, 0.0);
     make_2d_array_double(&LQ_imm_imm, n_imm_covariate, n_imm_covariate, 0.0);

     make_1d_array_double(&log_L_imm, n_imm_covariate, 0.0);
     make_2d_array_double(&log_L_imm_imm, n_imm_covariate, n_imm_covariate, 0.0);

     make_1d_array_double(&my_log_L_imm, n_imm_covariate, 0.0);
     make_2d_array_double(&my_log_L_imm_imm, n_imm_covariate, n_imm_covariate, 0.0);

     make_1d_array_double(&score_imm, n_imm_covariate, 0.0);
     make_2d_array_double(&info_imm_imm, n_imm_covariate, n_imm_covariate, 0.0);
  }

  if(n_b_mode > 0 && n_p_mode > 0)
  {
     make_2d_array_double(&log_day_L_lb_lp, n_b_mode, n_p_mode, 0.0);
     make_2d_array_double(&L_lb_lp, n_b_mode, n_p_mode, 0.0);
     make_2d_array_double(&LQ_lb_lp, n_b_mode, n_p_mode, 0.0);
     make_2d_array_double(&log_L_lb_lp, n_b_mode, n_p_mode, 0.0);
     make_2d_array_double(&my_log_L_lb_lp, n_b_mode, n_p_mode, 0.0);
     make_2d_array_double(&info_lb_lp, n_b_mode, n_p_mode, 0.0);
  }
  if(n_b_mode > 0 && n_u_mode > 0)
  {
     make_1d_array_double(&L_lb_lu, n_b_mode, 0.0);
     make_1d_array_double(&LQ_lb_lu, n_b_mode, 0.0);
     make_1d_array_double(&log_L_lb_lu, n_b_mode, 0.0);
     make_1d_array_double(&my_log_L_lb_lu, n_b_mode, 0.0);
     make_2d_array_double(&info_lb_lu, n_b_mode, n_u_mode, 0.0);
  }
  if(n_b_mode > 0 && n_q_mode > 0)
  {
     make_1d_array_double(&L_lb_lq, n_b_mode, 0.0);
     make_1d_array_double(&LQ_lb_lq, n_b_mode, 0.0);
     make_1d_array_double(&log_L_lb_lq, n_b_mode, 0.0);
     make_1d_array_double(&my_log_L_lb_lq, n_b_mode, 0.0);
     make_2d_array_double(&info_lb_lq, n_b_mode, n_q_mode, 0.0);
  }
  if(n_b_mode > 0 && n_c2p_covariate > 0)
  {
     make_2d_array_double(&log_f_lb_c2p, n_b_mode, n_c2p_covariate, 0.0);
     make_2d_array_double(&log_e_lb_c2p, n_b_mode, n_c2p_covariate, 0.0);
     make_2d_array_double(&cum_log_e_lb_c2p, n_b_mode, n_c2p_covariate, 0.0);
     make_2d_array_double(&temp_lb_c2p, n_b_mode, n_c2p_covariate, 0.0);
     make_3d_array_double(&log_ee_lb_c2p, max_epi_duration, n_b_mode, n_c2p_covariate, 0.0);
     make_3d_array_double(&cum_log_ee_lb_c2p, max_epi_duration, n_b_mode, n_c2p_covariate, 0.0);
     make_2d_array_double(&log_day_L_lb_c2p, n_b_mode, n_c2p_covariate, 0.0);
     make_2d_array_double(&L_lb_c2p, n_b_mode, n_c2p_covariate, 0.0);
     make_2d_array_double(&LQ_lb_c2p, n_b_mode, n_c2p_covariate, 0.0);
     make_2d_array_double(&log_L_lb_c2p, n_b_mode, n_c2p_covariate, 0.0);
     make_2d_array_double(&my_log_L_lb_c2p, n_b_mode, n_c2p_covariate, 0.0);
     make_2d_array_double(&info_lb_c2p, n_b_mode, n_c2p_covariate, 0.0);
  }
  if(n_b_mode > 0 && n_p2p_covariate > 0)
  {
     make_2d_array_double(&log_day_L_lb_p2p, n_b_mode, n_p2p_covariate, 0.0);
     make_2d_array_double(&L_lb_p2p, n_b_mode, n_p2p_covariate, 0.0);
     make_2d_array_double(&LQ_lb_p2p, n_b_mode, n_p2p_covariate, 0.0);
     make_2d_array_double(&log_L_lb_p2p, n_b_mode, n_p2p_covariate, 0.0);
     make_2d_array_double(&my_log_L_lb_p2p, n_b_mode, n_p2p_covariate, 0.0);
     make_2d_array_double(&info_lb_p2p, n_b_mode, n_p2p_covariate, 0.0);
  }
  if(n_b_mode > 0 && n_pat_covariate > 0)
  {
     make_2d_array_double(&L_lb_pat, n_b_mode, n_pat_covariate, 0.0);
     make_2d_array_double(&LQ_lb_pat, n_b_mode, n_pat_covariate, 0.0);
     make_2d_array_double(&log_L_lb_pat, n_b_mode, n_pat_covariate, 0.0);
     make_2d_array_double(&my_log_L_lb_pat, n_b_mode, n_pat_covariate, 0.0);
     make_2d_array_double(&info_lb_pat, n_b_mode, n_pat_covariate, 0.0);
  }
  if(n_b_mode > 0 && n_imm_covariate > 0)
  {
     make_2d_array_double(&L_lb_imm, n_b_mode, n_imm_covariate, 0.0);
     make_2d_array_double(&LQ_lb_imm, n_b_mode, n_imm_covariate, 0.0);
     make_2d_array_double(&log_L_lb_imm, n_b_mode, n_imm_covariate, 0.0);
     make_2d_array_double(&my_log_L_lb_imm, n_b_mode, n_imm_covariate, 0.0);
     make_2d_array_double(&info_lb_imm, n_b_mode, n_imm_covariate, 0.0);
  }
  if(n_p_mode > 0 && n_u_mode > 0)
  {
     make_1d_array_double(&L_lp_lu, n_p_mode, 0.0);
     make_1d_array_double(&LQ_lp_lu, n_p_mode, 0.0);
     make_1d_array_double(&log_L_lp_lu, n_p_mode, 0.0);
     make_1d_array_double(&my_log_L_lp_lu, n_p_mode, 0.0);
     make_2d_array_double(&info_lp_lu, n_p_mode, n_u_mode, 0.0);
  }
  if(n_p_mode > 0 && n_q_mode > 0)
  {
     make_1d_array_double(&L_lp_lq, n_p_mode, 0.0);
     make_1d_array_double(&LQ_lp_lq, n_p_mode, 0.0);
     make_1d_array_double(&log_L_lp_lq, n_p_mode, 0.0);
     make_1d_array_double(&my_log_L_lp_lq, n_p_mode, 0.0);
     make_2d_array_double(&info_lp_lq, n_p_mode, n_q_mode, 0.0);
  }
  if(n_p_mode > 0 && n_c2p_covariate > 0)
  {
     make_2d_array_double(&log_day_L_lp_c2p, n_p_mode, n_c2p_covariate, 0.0);
     make_2d_array_double(&L_lp_c2p, n_p_mode, n_c2p_covariate, 0.0);
     make_2d_array_double(&LQ_lp_c2p, n_p_mode, n_c2p_covariate, 0.0);
     make_2d_array_double(&log_L_lp_c2p, n_p_mode, n_c2p_covariate, 0.0);
     make_2d_array_double(&my_log_L_lp_c2p, n_p_mode, n_c2p_covariate, 0.0);
     make_2d_array_double(&info_lp_c2p, n_p_mode, n_c2p_covariate, 0.0);
  }
  if(n_p_mode > 0 && n_p2p_covariate > 0)
  {
     make_2d_array_double(&log_f_lp_p2p, n_p_mode, n_p2p_covariate, 0.0);
     make_2d_array_double(&log_e_lp_p2p, n_p_mode, n_p2p_covariate, 0.0);
     make_2d_array_double(&cum_log_e_lp_p2p, n_p_mode, n_p2p_covariate, 0.0);
     make_2d_array_double(&temp_lp_p2p, n_p_mode, n_p2p_covariate, 0.0);
     make_3d_array_double(&log_ee_lp_p2p, max_epi_duration, n_p_mode, n_p2p_covariate, 0.0);
     make_3d_array_double(&cum_log_ee_lp_p2p, max_epi_duration, n_p_mode, n_p2p_covariate, 0.0);
     make_2d_array_double(&log_day_L_lp_p2p, n_p_mode, n_p2p_covariate, 0.0);
     make_2d_array_double(&L_lp_p2p, n_p_mode, n_p2p_covariate, 0.0);
     make_2d_array_double(&LQ_lp_p2p, n_p_mode, n_p2p_covariate, 0.0);
     make_2d_array_double(&log_L_lp_p2p, n_p_mode, n_p2p_covariate, 0.0);
     make_2d_array_double(&my_log_L_lp_p2p, n_p_mode, n_p2p_covariate, 0.0);
     make_2d_array_double(&info_lp_p2p, n_p_mode, n_p2p_covariate, 0.0);
  }
  if(n_p_mode > 0 && n_pat_covariate > 0)
  {
     make_2d_array_double(&L_lp_pat, n_p_mode, n_pat_covariate, 0.0);
     make_2d_array_double(&LQ_lp_pat, n_p_mode, n_pat_covariate, 0.0);
     make_2d_array_double(&log_L_lp_pat, n_p_mode, n_pat_covariate, 0.0);
     make_2d_array_double(&my_log_L_lp_pat, n_p_mode, n_pat_covariate, 0.0);
     make_2d_array_double(&info_lp_pat, n_p_mode, n_pat_covariate, 0.0);
  }
  if(n_p_mode > 0 && n_imm_covariate > 0)
  {
     make_2d_array_double(&L_lp_imm, n_p_mode, n_imm_covariate, 0.0);
     make_2d_array_double(&LQ_lp_imm, n_p_mode, n_imm_covariate, 0.0);
     make_2d_array_double(&log_L_lp_imm, n_p_mode, n_imm_covariate, 0.0);
     make_2d_array_double(&my_log_L_lp_imm, n_p_mode, n_imm_covariate, 0.0);
     make_2d_array_double(&info_lp_imm, n_p_mode, n_imm_covariate, 0.0);
  }
  if(n_u_mode > 0 && n_q_mode > 0)
  {
     make_2d_array_double(&info_lu_lq, n_u_mode, n_q_mode, 0.0);
  }
  if(n_u_mode > 0 && n_c2p_covariate > 0)
  {
     make_1d_array_double(&L_lu_c2p, n_c2p_covariate, 0.0);
     make_1d_array_double(&LQ_lu_c2p, n_c2p_covariate, 0.0);
     make_1d_array_double(&log_L_lu_c2p, n_c2p_covariate, 0.0);
     make_1d_array_double(&my_log_L_lu_c2p, n_c2p_covariate, 0.0);
     make_2d_array_double(&info_lu_c2p, n_u_mode, n_c2p_covariate, 0.0);
  }
  if(n_u_mode > 0 && n_p2p_covariate > 0)
  {
     make_1d_array_double(&L_lu_p2p, n_p2p_covariate, 0.0);
     make_1d_array_double(&LQ_lu_p2p, n_p2p_covariate, 0.0);
     make_1d_array_double(&log_L_lu_p2p, n_p2p_covariate, 0.0);
     make_1d_array_double(&my_log_L_lu_p2p, n_p2p_covariate, 0.0);
     make_2d_array_double(&info_lu_p2p, n_u_mode, n_p2p_covariate, 0.0);
  }
  if(n_u_mode > 0 && n_pat_covariate > 0)
  {
     make_1d_array_double(&pr_lu_pat, n_pat_covariate, 0.0);
     make_1d_array_double(&U_lu_pat, n_pat_covariate, 0.0);
     make_1d_array_double(&L_lu_pat, n_pat_covariate, 0.0);
     make_1d_array_double(&LQ_lu_pat, n_pat_covariate, 0.0);
     make_1d_array_double(&log_L_lu_pat, n_pat_covariate, 0.0);
     make_1d_array_double(&my_log_L_lu_pat, n_pat_covariate, 0.0);
     make_2d_array_double(&info_lu_pat, n_u_mode, n_pat_covariate, 0.0);
  }
  if(n_u_mode > 0 && n_imm_covariate > 0)
  {
     make_1d_array_double(&L_lu_imm, n_imm_covariate, 0.0);
     make_1d_array_double(&LQ_lu_imm, n_imm_covariate, 0.0);
     make_1d_array_double(&log_L_lu_imm, n_imm_covariate, 0.0);
     make_1d_array_double(&my_log_L_lu_imm, n_imm_covariate, 0.0);
     make_2d_array_double(&info_lu_imm, n_u_mode, n_imm_covariate, 0.0);
  }

  if(n_q_mode > 0 && n_c2p_covariate > 0)
  {
     make_1d_array_double(&L_lq_c2p, n_c2p_covariate, 0.0);
     make_1d_array_double(&LQ_lq_c2p, n_c2p_covariate, 0.0);
     make_1d_array_double(&log_L_lq_c2p, n_c2p_covariate, 0.0);
     make_1d_array_double(&my_log_L_lq_c2p, n_c2p_covariate, 0.0);
     make_2d_array_double(&info_lq_c2p, n_q_mode, n_c2p_covariate, 0.0);
  }
  if(n_q_mode > 0 && n_p2p_covariate > 0)
  {
     make_1d_array_double(&L_lq_p2p, n_p2p_covariate, 0.0);
     make_1d_array_double(&LQ_lq_p2p, n_p2p_covariate, 0.0);
     make_1d_array_double(&log_L_lq_p2p, n_p2p_covariate, 0.0);
     make_1d_array_double(&my_log_L_lq_p2p, n_p2p_covariate, 0.0);
     make_2d_array_double(&info_lq_p2p, n_q_mode, n_p2p_covariate, 0.0);
  }

  if(n_q_mode > 0 && n_pat_covariate > 0)
  {
     make_1d_array_double(&L_lq_pat, n_pat_covariate, 0.0);
     make_1d_array_double(&LQ_lq_pat, n_pat_covariate, 0.0);
     make_1d_array_double(&log_L_lq_pat, n_pat_covariate, 0.0);
     make_1d_array_double(&my_log_L_lq_pat, n_pat_covariate, 0.0);
     make_2d_array_double(&info_lq_pat, n_q_mode, n_pat_covariate, 0.0);
  }
  if(n_q_mode > 0 && n_imm_covariate > 0)
  {
     make_1d_array_double(&Q_lq_imm, n_imm_covariate, 0.0);
     make_1d_array_double(&L_lq_imm, n_imm_covariate, 0.0);
     make_1d_array_double(&LQ_lq_imm, n_imm_covariate, 0.0);
     make_1d_array_double(&log_L_lq_imm, n_imm_covariate, 0.0);
     make_1d_array_double(&my_log_L_lq_imm, n_imm_covariate, 0.0);
     make_2d_array_double(&info_lq_imm, n_q_mode, n_imm_covariate, 0.0);
  }
  if(n_c2p_covariate > 0 && n_p2p_covariate > 0)
  {
     make_2d_array_double(&log_day_L_c2p_p2p, n_c2p_covariate, n_p2p_covariate, 0.0);
     make_2d_array_double(&L_c2p_p2p, n_c2p_covariate, n_p2p_covariate, 0.0);
     make_2d_array_double(&log_L_c2p_p2p, n_c2p_covariate, n_p2p_covariate, 0.0);
     make_2d_array_double(&my_log_L_c2p_p2p, n_c2p_covariate, n_p2p_covariate, 0.0);
     make_2d_array_double(&info_c2p_p2p, n_c2p_covariate, n_p2p_covariate, 0.0);
  }

  if(n_c2p_covariate > 0 && n_pat_covariate > 0)
  {
     make_2d_array_double(&L_c2p_pat, n_c2p_covariate, n_pat_covariate, 0.0);
     make_2d_array_double(&LQ_c2p_pat, n_c2p_covariate, n_pat_covariate, 0.0);
     make_2d_array_double(&log_L_c2p_pat, n_c2p_covariate, n_pat_covariate, 0.0);
     make_2d_array_double(&my_log_L_c2p_pat, n_c2p_covariate, n_pat_covariate, 0.0);
     make_2d_array_double(&info_c2p_pat, n_c2p_covariate, n_pat_covariate, 0.0);
  }
  if(n_c2p_covariate > 0 && n_imm_covariate > 0)
  {
     make_2d_array_double(&L_c2p_imm, n_c2p_covariate, n_imm_covariate, 0.0);
     make_2d_array_double(&LQ_c2p_imm, n_c2p_covariate, n_imm_covariate, 0.0);
     make_2d_array_double(&log_L_c2p_imm, n_c2p_covariate, n_imm_covariate, 0.0);
     make_2d_array_double(&my_log_L_c2p_imm, n_c2p_covariate, n_imm_covariate, 0.0);
     make_2d_array_double(&info_c2p_imm, n_c2p_covariate, n_imm_covariate, 0.0);
  }
  if(n_p2p_covariate > 0 && n_pat_covariate > 0)
  {
     make_2d_array_double(&L_p2p_pat, n_p2p_covariate, n_pat_covariate, 0.0);
     make_2d_array_double(&LQ_p2p_pat, n_p2p_covariate, n_pat_covariate, 0.0);
     make_2d_array_double(&log_L_p2p_pat, n_p2p_covariate, n_pat_covariate, 0.0);
     make_2d_array_double(&my_log_L_p2p_pat, n_p2p_covariate, n_pat_covariate, 0.0);
     make_2d_array_double(&info_p2p_pat, n_p2p_covariate, n_pat_covariate, 0.0);
  }

  if(n_p2p_covariate > 0 && n_imm_covariate > 0)
  {
     make_2d_array_double(&L_p2p_imm, n_p2p_covariate, n_imm_covariate, 0.0);
     make_2d_array_double(&LQ_p2p_imm, n_p2p_covariate, n_imm_covariate, 0.0);
     make_2d_array_double(&log_L_p2p_imm, n_p2p_covariate, n_imm_covariate, 0.0);
     make_2d_array_double(&my_log_L_p2p_imm, n_p2p_covariate, n_imm_covariate, 0.0);
     make_2d_array_double(&info_p2p_imm, n_p2p_covariate, n_imm_covariate, 0.0);
  }

  if(n_pat_covariate > 0 && n_imm_covariate > 0)
  {
     make_2d_array_double(&L_pat_imm, n_pat_covariate, n_imm_covariate, 0.0);
     make_2d_array_double(&LQ_pat_imm, n_pat_covariate, n_imm_covariate, 0.0);
     make_2d_array_double(&log_L_pat_imm, n_pat_covariate, n_imm_covariate, 0.0);
     make_2d_array_double(&my_log_L_pat_imm, n_pat_covariate, n_imm_covariate, 0.0);
     make_2d_array_double(&info_pat_imm, n_pat_covariate, n_imm_covariate, 0.0);
  }
}


void free_arrays(void)
{
  int n_par, n_b_mode, n_p_mode, n_q_mode, n_u_mode;
  int n_c2p_covariate, n_p2p_covariate, n_imm_covariate, n_pat_covariate;
  int n_covariate, n_time_ind_covariate, n_time_dep_covariate;
  int n_sus_p2p_covariate, n_inf_p2p_covariate, n_int_p2p_covariate;
  int n_par_equiclass;

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
  n_imm_covariate = cfg_pars.n_imm_covariate;
  n_pat_covariate = cfg_pars.n_pat_covariate;
  n_par = cfg_pars.n_par;
  n_par_equiclass = cfg_pars.n_par_equiclass;

  if(ee != NULL)  free(ee);
  if(cum_log_ee != NULL)  free(cum_log_ee);

  if(n_b_mode > 0)
  {
     if(b != NULL)  free(b);
     if(lb != NULL)  free(lb);
     if(se_b != NULL)  free(se_b);
     if(se_lb != NULL)  free(se_lb);
     if(lower_b != NULL)  free(lower_b);
     if(upper_b != NULL)  free(upper_b);

     if(CPI != NULL)  free(CPI);
     if(se_CPI != NULL)  free(se_CPI);
     if(lower_CPI != NULL)  free(lower_CPI);
     if(upper_CPI != NULL)  free(upper_CPI);

     if(log_f_lb != NULL)  free(log_f_lb);
     if(log_f_lb_lb != NULL)  free(log_f_lb_lb);

     if(log_e_lb != NULL)  free(log_e_lb);
     if(log_e_lb_lb != NULL)  free(log_e_lb_lb);

     if(cum_log_e_lb != NULL)  free(cum_log_e_lb);
     if(cum_log_e_lb_lb != NULL)  free(cum_log_e_lb_lb);

     if(temp_lb != NULL)  free(temp_lb);
     if(temp_lb_lb != NULL)  free(temp_lb_lb);

     if(log_ee_lb != NULL)  free_2d_array_double(log_ee_lb);
     if(log_ee_lb_lb != NULL)  free_2d_array_double(log_ee_lb_lb);

     if(cum_log_ee_lb != NULL)  free_2d_array_double(cum_log_ee_lb);
     if(cum_log_ee_lb_lb != NULL)  free_2d_array_double(cum_log_ee_lb_lb);

     if(log_day_L_lb != NULL)  free(log_day_L_lb);
     if(log_day_L_lb_lb != NULL)  free_2d_array_double(log_day_L_lb_lb);

     if(L_lb != NULL)  free(L_lb);
     if(L_lb_lb != NULL)  free_2d_array_double(L_lb_lb);

     if(LQ_lb != NULL)  free(LQ_lb);
     if(LQ_lb_lb != NULL)  free_2d_array_double(LQ_lb_lb);

     if(log_L_lb != NULL)  free(log_L_lb);
     if(log_L_lb_lb != NULL)  free_2d_array_double(log_L_lb_lb);

     if(my_log_L_lb != NULL)  free(my_log_L_lb);
     if(my_log_L_lb_lb != NULL)  free_2d_array_double(my_log_L_lb_lb);

     if(score_lb != NULL)  free(score_lb);
     if(info_lb_lb != NULL)  free_2d_array_double(info_lb_lb);

  }

  if(n_p_mode > 0)
  {
     free(p);
     free(lp);
     free(se_p);
     free(se_lp);
     free(lower_p);
     free(upper_p);

     free(SAR);
     free(se_SAR);
     free(lower_SAR);
     free(upper_SAR);

     free(R0);
     free(se_R0);
     free(lower_R0);
     free(upper_R0);
     
     free(log_f_lp);
     free(log_f_lp_lp);

     free(log_e_lp);
     free(log_e_lp_lp);

     free(cum_log_e_lp);
     free(cum_log_e_lp_lp);

     free(temp_lp);
     free(temp_lp_lp);

     free_2d_array_double(log_ee_lp);
     free_2d_array_double(log_ee_lp_lp);

     free_2d_array_double(cum_log_ee_lp);
     free_2d_array_double(cum_log_ee_lp_lp);

     free(log_day_L_lp);
     free_2d_array_double(log_day_L_lp_lp);

     free(L_lp);
     free_2d_array_double(L_lp_lp);

     free(LQ_lp);
     free_2d_array_double(LQ_lp_lp);

     free(log_L_lp);
     free_2d_array_double(log_L_lp_lp);

     free(my_log_L_lp);
     free_2d_array_double(my_log_L_lp_lp);

     free(score_lp);
     free_2d_array_double(info_lp_lp);

  }

  if(n_u_mode > 0)
  {
     free(u);
     free(lu);
     free(se_u);
     free(se_lu);
     free(lower_u);
     free(upper_u);

     free(score_lu);
     free_2d_array_double(info_lu_lu);
  }

  if(n_q_mode > 0)
  {
     free(q);
     free(lq);
     free(se_q);
     free(se_lq);
     free(lower_q);
     free(upper_q);

     free(score_lq);
     free_2d_array_double(info_lq_lq);
  }

  if(n_c2p_covariate > 0)
  {
     free(coeff_c2p);
     free(se_coeff_c2p);
     free(OR_c2p);
     free(se_OR_c2p);
     free(lower_OR_c2p);
     free(upper_OR_c2p);
     free(log_f_c2p);
     free_2d_array_double(log_f_c2p_c2p);

     free(log_e_c2p);
     free_2d_array_double(log_e_c2p_c2p);

     free(cum_log_e_c2p);
     free_2d_array_double(cum_log_e_c2p_c2p);

     free(temp_c2p);
     free_2d_array_double(temp_c2p_c2p);

     free_2d_array_double(log_ee_c2p);
     free_3d_array_double(log_ee_c2p_c2p);

     free_2d_array_double(cum_log_ee_c2p);
     free_3d_array_double(cum_log_ee_c2p_c2p);

     free(log_day_L_c2p);
     free_2d_array_double(log_day_L_c2p_c2p);

     free(L_c2p);
     free_2d_array_double(L_c2p_c2p);

     free(LQ_c2p);
     free_2d_array_double(LQ_c2p_c2p);

     free(log_L_c2p);
     free_2d_array_double(log_L_c2p_c2p);

     free(my_log_L_c2p);
     free_2d_array_double(my_log_L_c2p_c2p);

     free(score_c2p);
     free_2d_array_double(info_c2p_c2p);
  }

  if(n_p2p_covariate > 0)
  {
     free(coeff_p2p);
     free(se_coeff_p2p);
     free(OR_p2p);
     free(se_OR_p2p);
     free(lower_OR_p2p);
     free(upper_OR_p2p);
     free(log_f_p2p);
     free_2d_array_double(log_f_p2p_p2p);

     free(log_e_p2p);
     free_2d_array_double(log_e_p2p_p2p);

     free(cum_log_e_p2p);
     free_2d_array_double(cum_log_e_p2p_p2p);

     free(temp_p2p);
     free_2d_array_double(temp_p2p_p2p);

     free_2d_array_double(log_ee_p2p);
     free_3d_array_double(log_ee_p2p_p2p);

     free_2d_array_double(cum_log_ee_p2p);
     free_3d_array_double(cum_log_ee_p2p_p2p);

     free(log_day_L_p2p);
     free_2d_array_double(log_day_L_p2p_p2p);

     free(L_p2p);
     free_2d_array_double(L_p2p_p2p);

     free(LQ_p2p);
     free_2d_array_double(LQ_p2p_p2p);

     free(log_L_p2p);
     free_2d_array_double(log_L_p2p_p2p);

     free(my_log_L_p2p);
     free_2d_array_double(my_log_L_p2p_p2p);

     free(score_p2p);
     free_2d_array_double(info_p2p_p2p);

  }
  if(n_pat_covariate > 0)
  {
     free(coeff_pat);
     free(se_coeff_pat);
     free(OR_pat);
     free(se_OR_pat);
     free(lower_OR_pat);
     free(upper_OR_pat);
     free(pr_pat);
     free_2d_array_double(pr_pat_pat);
     free(U_pat);
     free_2d_array_double(U_pat_pat);

     free(L_pat);
     free_2d_array_double(L_pat_pat);

     free(LQ_pat);
     free_2d_array_double(LQ_pat_pat);

     free(log_L_pat);
     free_2d_array_double(log_L_pat_pat);

     free(my_log_L_pat);
     free_2d_array_double(my_log_L_pat_pat);

     free(score_pat);
     free_2d_array_double(info_pat_pat);
  }

  if(n_imm_covariate > 0)
  {
     free(coeff_imm);
     free(se_coeff_imm);
     free(OR_imm);
     free(se_OR_imm);
     free(lower_OR_imm);
     free(upper_OR_imm);
     free(Q_imm);
     free_2d_array_double(Q_imm_imm);

     free(L_imm);
     free_2d_array_double(L_imm_imm);

     free(LQ_imm);
     free_2d_array_double(LQ_imm_imm);

     free(log_L_imm);
     free_2d_array_double(log_L_imm_imm);

     free(my_log_L_imm);
     free_2d_array_double(my_log_L_imm_imm);

     free(score_imm);
     free_2d_array_double(info_imm_imm);
  }

  if(n_b_mode > 0 && n_p_mode > 0)
  {
     free_2d_array_double(log_day_L_lb_lp);
     free_2d_array_double(L_lb_lp);
     free_2d_array_double(LQ_lb_lp);
     free_2d_array_double(log_L_lb_lp);
     free_2d_array_double(my_log_L_lb_lp);
     free_2d_array_double(info_lb_lp);
  }

  if(n_b_mode > 0 && n_u_mode > 0)
  {
     free(L_lb_lu);
     free(LQ_lb_lu);
     free(log_L_lb_lu);
     free(my_log_L_lb_lu);
     free_2d_array_double(info_lb_lu);
  }
  if(n_b_mode > 0 && n_q_mode > 0)
  {
     free(L_lb_lq);
     free(LQ_lb_lq);
     free(log_L_lb_lq);
     free(my_log_L_lb_lq);
     free_2d_array_double(info_lb_lq);
  }
  if(n_b_mode > 0 && n_c2p_covariate > 0)
  {
     free_2d_array_double(log_f_lb_c2p);
     free_2d_array_double(log_e_lb_c2p);
     free_2d_array_double(cum_log_e_lb_c2p);
     free_2d_array_double(temp_lb_c2p);
     free_3d_array_double(log_ee_lb_c2p);
     free_3d_array_double(cum_log_ee_lb_c2p);
     free_2d_array_double(log_day_L_lb_c2p);
     free_2d_array_double(L_lb_c2p);
     free_2d_array_double(LQ_lb_c2p);
     free_2d_array_double(log_L_lb_c2p);
     free_2d_array_double(my_log_L_lb_c2p);
     free_2d_array_double(info_lb_c2p);
  }

  if(n_b_mode > 0 && n_p2p_covariate > 0)
  {
     free_2d_array_double(log_day_L_lb_p2p);
     free_2d_array_double(L_lb_p2p);
     free_2d_array_double(LQ_lb_p2p);
     free_2d_array_double(log_L_lb_p2p);
     free_2d_array_double(my_log_L_lb_p2p);
     free_2d_array_double(info_lb_p2p);
  }
  if(n_b_mode > 0 && n_pat_covariate > 0)
  {
     free_2d_array_double(L_lb_pat);
     free_2d_array_double(LQ_lb_pat);
     free_2d_array_double(log_L_lb_pat);
     free_2d_array_double(my_log_L_lb_pat);
     free_2d_array_double(info_lb_pat);
  }
  if(n_b_mode > 0 && n_imm_covariate > 0)
  {
     if(L_lb_imm != NULL)  free_2d_array_double(L_lb_imm);
     if(LQ_lb_imm != NULL)  free_2d_array_double(LQ_lb_imm);
     if(log_L_lb_imm != NULL)  free_2d_array_double(log_L_lb_imm);
     if(my_log_L_lb_imm != NULL)  free_2d_array_double(my_log_L_lb_imm);
     if(info_lb_imm != NULL)  free_2d_array_double(info_lb_imm);
  }

  if(n_p_mode > 0 && n_u_mode > 0)
  {
     free(L_lp_lu);
     free(LQ_lp_lu);
     free(log_L_lp_lu);
     free(my_log_L_lp_lu);
     free_2d_array_double(info_lp_lu);
  }
  if(n_p_mode > 0 && n_q_mode > 0)
  {
     free(L_lp_lq);
     free(LQ_lp_lq);
     free(log_L_lp_lq);
     free(my_log_L_lp_lq);
     free_2d_array_double(info_lp_lq);
  }
  if(n_p_mode > 0 && n_c2p_covariate > 0)
  {
     free_2d_array_double(log_day_L_lp_c2p);
     free_2d_array_double(L_lp_c2p);
     free_2d_array_double(LQ_lp_c2p);
     free_2d_array_double(log_L_lp_c2p);
     free_2d_array_double(my_log_L_lp_c2p);
     free_2d_array_double(info_lp_c2p);
  }

  if(n_p_mode > 0 && n_p2p_covariate > 0)
  {
     free_2d_array_double(log_f_lp_p2p);
     free_2d_array_double(log_e_lp_p2p);
     free_2d_array_double(cum_log_e_lp_p2p);
     free_2d_array_double(temp_lp_p2p);
     free_3d_array_double(log_ee_lp_p2p);
     free_3d_array_double(cum_log_ee_lp_p2p);
     free_2d_array_double(log_day_L_lp_p2p);
     free_2d_array_double(L_lp_p2p);
     free_2d_array_double(LQ_lp_p2p);
     free_2d_array_double(log_L_lp_p2p);
     free_2d_array_double(my_log_L_lp_p2p);
     free_2d_array_double(info_lp_p2p);
  }
  if(n_p_mode > 0 && n_pat_covariate > 0)
  {
     free_2d_array_double(L_lp_pat);
     free_2d_array_double(LQ_lp_pat);
     free_2d_array_double(log_L_lp_pat);
     free_2d_array_double(my_log_L_lp_pat);
     free_2d_array_double(info_lp_pat);
  }
  if(n_p_mode > 0 && n_imm_covariate > 0)
  {
     if(L_lp_imm != NULL)  free_2d_array_double(L_lp_imm);
     if(LQ_lp_imm != NULL)  free_2d_array_double(LQ_lp_imm);
     if(log_L_lp_imm != NULL)  free_2d_array_double(log_L_lp_imm);
     if(my_log_L_lp_imm != NULL)  free_2d_array_double(my_log_L_lp_imm);
     if(info_lp_imm != NULL)  free_2d_array_double(info_lp_imm);
  }
  if(n_u_mode > 0 && n_q_mode > 0)
  {
     free_2d_array_double(info_lu_lq);
  }

  if(n_u_mode > 0 && n_c2p_covariate > 0)
  {
     free(L_lu_c2p);
     free(LQ_lu_c2p);
     free(log_L_lu_c2p);
     free(my_log_L_lu_c2p);
     free_2d_array_double(info_lu_c2p);
  }
  if(n_u_mode > 0 && n_p2p_covariate > 0)
  {
     free(L_lu_p2p);
     free(LQ_lu_p2p);
     free(log_L_lu_p2p);
     free(my_log_L_lu_p2p);
     free_2d_array_double(info_lu_p2p);
  }
  if(n_u_mode > 0 && n_pat_covariate > 0)
  {
     free(pr_lu_pat);
     free(U_lu_pat);
     free(L_lu_pat);
     free(LQ_lu_pat);
     free(log_L_lu_pat);
     free(my_log_L_lu_pat);
     free_2d_array_double(info_lu_pat);
  }
  if(n_u_mode > 0 && n_imm_covariate > 0)
  {
     if(L_lu_imm != NULL)  free(L_lu_imm);
     if(LQ_lu_imm != NULL)  free(LQ_lu_imm);
     if(log_L_lu_imm != NULL)  free(log_L_lu_imm);
     if(my_log_L_lu_imm != NULL)  free(my_log_L_lu_imm);
     if(info_lu_imm != NULL)  free_2d_array_double(info_lu_imm);
  }

  if(n_q_mode > 0 && n_c2p_covariate > 0)
  {
     free(L_lq_c2p);
     free(LQ_lq_c2p);
     free(log_L_lq_c2p);
     free(my_log_L_lq_c2p);
     free_2d_array_double(info_lq_c2p);
  }
  if(n_q_mode > 0 && n_p2p_covariate > 0)
  {
     free(L_lq_p2p);
     free(LQ_lq_p2p);
     free(log_L_lq_p2p);
     free(my_log_L_lq_p2p);
     free_2d_array_double(info_lq_p2p);
  }
  if(n_q_mode > 0 && n_pat_covariate > 0)
  {
     if(L_lq_pat != NULL)  free(L_lq_pat);
     if(LQ_lq_pat != NULL)  free(LQ_lq_pat);
     if(log_L_lq_pat != NULL)  free(log_L_lq_pat);
     if(my_log_L_lq_pat != NULL)  free(my_log_L_lq_pat);
     if(info_lq_pat != NULL)  free_2d_array_double(info_lq_pat);
  }
  if(n_q_mode > 0 && n_imm_covariate > 0)
  {
     if(Q_lq_imm != NULL)  free(Q_lq_imm);
     if(L_lq_imm != NULL)  free(L_lq_imm);
     if(LQ_lq_imm != NULL)  free(LQ_lq_imm);
     if(log_L_lq_imm != NULL)  free(log_L_lq_imm);
     if(my_log_L_lq_imm != NULL)  free(my_log_L_lq_imm);
     if(info_lq_imm != NULL)  free_2d_array_double(info_lq_imm);
  }

  if(n_c2p_covariate > 0 && n_p2p_covariate > 0)
  {
     free_2d_array_double(log_day_L_c2p_p2p);
     free_2d_array_double(L_c2p_p2p);
     free_2d_array_double(LQ_c2p_p2p);
     free_2d_array_double(log_L_c2p_p2p);
     free_2d_array_double(my_log_L_c2p_p2p);
     free_2d_array_double(info_c2p_p2p);
  }

  if(n_c2p_covariate > 0 && n_pat_covariate > 0)
  {
     if(L_c2p_pat != NULL)  free_2d_array_double(L_c2p_pat);
     if(LQ_c2p_pat != NULL)  free_2d_array_double(LQ_c2p_pat);
     if(log_L_c2p_pat != NULL)  free_2d_array_double(log_L_c2p_pat);
     if(my_log_L_c2p_pat != NULL)  free_2d_array_double(my_log_L_c2p_pat);
     if(info_c2p_pat != NULL)  free_2d_array_double(info_c2p_pat);
  }
  if(n_c2p_covariate > 0 && n_imm_covariate > 0)
  {
     if(L_c2p_imm != NULL)  free_2d_array_double(L_c2p_imm);
     if(LQ_c2p_imm != NULL)  free_2d_array_double(LQ_c2p_imm);
     if(log_L_c2p_imm != NULL)  free_2d_array_double(log_L_c2p_imm);
     if(my_log_L_c2p_imm != NULL)  free_2d_array_double(my_log_L_c2p_imm);
     if(info_c2p_imm != NULL)  free_2d_array_double(info_c2p_imm);
  }
  if(n_p2p_covariate > 0 && n_pat_covariate > 0)
  {
     if(L_p2p_pat != NULL)  free_2d_array_double(L_p2p_pat);
     if(LQ_p2p_pat != NULL)  free_2d_array_double(LQ_p2p_pat);
     if(log_L_p2p_pat != NULL)  free_2d_array_double(log_L_p2p_pat);
     if(my_log_L_p2p_pat != NULL)  free_2d_array_double(my_log_L_p2p_pat);
     if(info_p2p_pat != NULL)  free_2d_array_double(info_p2p_pat);
  }

  if(n_p2p_covariate > 0 && n_imm_covariate > 0)
  {
     if(L_p2p_imm != NULL)  free_2d_array_double(L_p2p_imm);
     if(LQ_p2p_imm != NULL)  free_2d_array_double(LQ_p2p_imm);
     if(log_L_p2p_imm != NULL)  free_2d_array_double(log_L_p2p_imm);
     if(my_log_L_p2p_imm != NULL)  free_2d_array_double(my_log_L_p2p_imm);
     if(info_p2p_imm != NULL)  free_2d_array_double(info_p2p_imm);
  }

  if(n_pat_covariate > 0 && n_imm_covariate > 0)
  {
     if(L_pat_imm != NULL)  free_2d_array_double(L_pat_imm);
     if(LQ_pat_imm != NULL)  free_2d_array_double(LQ_pat_imm);
     if(log_L_pat_imm != NULL)  free_2d_array_double(log_L_pat_imm);
     if(my_log_L_pat_imm != NULL)  free_2d_array_double(my_log_L_pat_imm);
     if(info_pat_imm != NULL)  free_2d_array_double(info_pat_imm);
  }
}
//this function extracts covariates affecting c2p transmission for model fitting from
//all time-independent and time-dependent covariates.
//the extracted covariates are stored in array c2p_covariate.
void organize_c2p_covariate(int t, PEOPLE *person, double *c2p_covariate)
{
   int h, i, j, k, l, m, n, r;
   int n_covariate, n_time_ind_covariate, n_time_dep_covariate;
   double *covariate;

   n_time_ind_covariate = cfg_pars.n_time_ind_covariate;
   n_time_dep_covariate = cfg_pars.n_time_dep_covariate;
   n_covariate = cfg_pars.n_covariate;
   covariate = NULL;

   if(n_covariate > 0)  
   {
      covariate = (double *)malloc((size_t) (n_covariate * sizeof(double)));

      h = person->community;
      r = t - community[h].day_epi_start;

      for(k=0; k<n_time_ind_covariate; k++)
         covariate[k] = person->time_ind_covariate[k];
      for(k=0; k<n_time_dep_covariate; k++)
         covariate[n_time_ind_covariate + k] = person->time_dep_covariate[r][k];

      l = 0;
      for(k=0; k<cfg_pars.n_c2p_covariate; k++)
      {
         m = cfg_pars.c2p_covariate[k] - 1;
         c2p_covariate[l++] = covariate[m];
      }
      free(covariate);
   }
}

//this function extracts covariates affecting p2p transmission for model fitting from
//all time-independent and time-dependent covariates.
//the extracted covariates are stored in array p2p_covariate.
void organize_p2p_covariate(int t, PEOPLE *sus_person, PEOPLE *inf_person, double *p2p_covariate)
{
   int h, i, j, k, l, m, n, r;
   int n_covariate, n_time_ind_covariate, n_time_dep_covariate;
   int n_sus_p2p_covariate, n_inf_p2p_covariate, n_int_p2p_covariate;
   double *sus_covariate, *inf_covariate;
   double *sus_p2p_covariate, *inf_p2p_covariate, *int_p2p_covariate;

   n_time_ind_covariate = cfg_pars.n_time_ind_covariate;
   n_time_dep_covariate = cfg_pars.n_time_dep_covariate;
   n_covariate = cfg_pars.n_covariate;
   n_sus_p2p_covariate = cfg_pars.n_sus_p2p_covariate;
   n_inf_p2p_covariate = cfg_pars.n_inf_p2p_covariate;
   n_int_p2p_covariate = cfg_pars.n_int_p2p_covariate;

   sus_covariate = inf_covariate = sus_p2p_covariate =
   inf_p2p_covariate = int_p2p_covariate = NULL;

   if(n_covariate > 0)  
   {
      sus_covariate = (double *)malloc((size_t) (n_covariate * sizeof(double)));
      inf_covariate = (double *)malloc((size_t) (n_covariate * sizeof(double)));
      if(n_sus_p2p_covariate > 0)  sus_p2p_covariate = (double *)malloc((size_t) (n_sus_p2p_covariate * sizeof(double)));
      if(n_inf_p2p_covariate > 0)  inf_p2p_covariate = (double *)malloc((size_t) (n_inf_p2p_covariate * sizeof(double)));
      if(n_int_p2p_covariate > 0)  int_p2p_covariate = (double *)malloc((size_t) (n_int_p2p_covariate * sizeof(double)));

      if(sus_person == NULL)
      {
         if(cfg_pars.SAR_covariate_provided == 0)
         {
            printf("Covariates for calculating SAR are not provided. The function organize_p2p_covariate should not be called\n");
            exit(0);
         }   
         r = t - cfg_pars.SAR_time_dep_lower;
         for(k=0; k<n_time_ind_covariate; k++)
            sus_covariate[k] = cfg_pars.SAR_sus_time_ind_covariate[k];
         for(k=0; k<n_time_dep_covariate; k++)
            sus_covariate[n_time_ind_covariate + k] = cfg_pars.SAR_sus_time_dep_covariate.data[r][k];

         for(k=0; k<n_time_ind_covariate; k++)
            inf_covariate[k] = cfg_pars.SAR_inf_time_ind_covariate[k];
         for(k=0; k<n_time_dep_covariate; k++)
            inf_covariate[n_time_ind_covariate + k] = cfg_pars.SAR_inf_time_dep_covariate.data[r][k];
      }      
      else
      {
         h = sus_person->community;
         r = t - community[h].day_epi_start;
         
         for(k=0; k<n_time_ind_covariate; k++)
            sus_covariate[k] = sus_person->time_ind_covariate[k];
         for(k=0; k<n_time_dep_covariate; k++)
            sus_covariate[n_time_ind_covariate + k] = sus_person->time_dep_covariate[r][k];

         for(k=0; k<n_time_ind_covariate; k++)
            inf_covariate[k] = inf_person->time_ind_covariate[k];
         for(k=0; k<n_time_dep_covariate; k++)
            inf_covariate[n_time_ind_covariate + k] = inf_person->time_dep_covariate[r][k];
      }
      l = 0;
      for(k=0; k<n_sus_p2p_covariate; k++)
      {
         m = cfg_pars.sus_p2p_covariate[k] - 1;
         sus_p2p_covariate[l++] = sus_covariate[m];
      }

      l = 0;
      for(k=0; k<n_inf_p2p_covariate; k++)
      {
         m = cfg_pars.inf_p2p_covariate[k] - 1;
         inf_p2p_covariate[l++] = inf_covariate[m];
      }

      l = 0;
      for(k=0; k<cfg_pars.n_int_p2p_covariate; k++)
      {
         m = cfg_pars.interaction[k][0] - 1;
         n = cfg_pars.interaction[k][1] - 1;
         int_p2p_covariate[l++] = sus_covariate[m] * inf_covariate[n];
      }

      l = 0;
      for(k=0; k<n_sus_p2p_covariate; k++)   p2p_covariate[l++] = sus_p2p_covariate[k];
      for(k=0; k<n_inf_p2p_covariate; k++)   p2p_covariate[l++] = inf_p2p_covariate[k];
      for(k=0; k<n_int_p2p_covariate; k++)   p2p_covariate[l++] = int_p2p_covariate[k];

      free(sus_covariate);
      free(inf_covariate);
      free(sus_p2p_covariate);
      free(inf_p2p_covariate);
      free(int_p2p_covariate);
   }
}


//this function extracts covariates affecting preseason immunity for model fitting from
//all time-independent covariates. It makes no sense to think of any time-dependent covariate
//affecting preseason immunity level. The extracted covariates are stored in array imm_covariate of each subject.
void organize_imm_covariate(PEOPLE *person)
{
   int h, i, j, k, l, m, n, r;
   int n_time_ind_covariate;
   double *covariate;

   n_time_ind_covariate = cfg_pars.n_time_ind_covariate;
   covariate = NULL;

   if(n_time_ind_covariate > 0)  
   {
      covariate = (double *)malloc((size_t) (n_time_ind_covariate * sizeof(double)));

      h = person->community;

      for(k=0; k<n_time_ind_covariate; k++)
         covariate[k] = person->time_ind_covariate[k];

      l = 0;
      for(k=0; k<cfg_pars.n_imm_covariate; k++)
      {
         m = cfg_pars.imm_covariate[k] - 1;
         person->imm_covariate[l++] = covariate[m];
      }
      free(covariate);
   }
}

//this function extracts covariates affecting pathpgenicity
//The extracted covariates are stored in 2-dim array pat_covariate of each subject.
void organize_pat_covariate(PEOPLE *person)
{
   int h, i, j, k, l, m, n, r, t;
   int n_time_ind_covariate, n_time_dep_covariate, n_covariate;
   double *covariate;

   n_time_ind_covariate = cfg_pars.n_time_ind_covariate;
   n_time_dep_covariate = cfg_pars.n_time_dep_covariate;
   n_covariate = cfg_pars.n_covariate;
   covariate = NULL;

   if(n_covariate > 0) 
   {
      covariate = (double *)malloc((size_t) (n_covariate * sizeof(double)));

      h = person->community;

      for(k=0; k<n_time_ind_covariate; k++)
         covariate[k] = person->time_ind_covariate[k];
      for(t=community[h].day_epi_start; t<=community[h].day_epi_stop; t++)
      {
         r = t - community[h].day_epi_start;
         for(k=0; k<n_time_dep_covariate; k++)
            covariate[n_time_ind_covariate + k] = person->time_dep_covariate[r][k];

         l = 0;
         for(k=0; k<cfg_pars.n_pat_covariate; k++)
         {
            m = cfg_pars.pat_covariate[k] - 1;
            person->pat_covariate[l++][r] = covariate[m];
         }
      }
      free(covariate);
   }
}

void add_c2p_contact_history_to_community(int h, int start_day, int stop_day, int contact_mode, double offset)
{
   int i, j, k, l, m, n, r, t;
   int found, match;
   CONTACT *ptr_contact;

   //printf("start:%d  stop:%d\n", start_day, stop_day);
   for(t = start_day; t <= stop_day; t++)
   {
      if(t >= community[h].day_epi_start && t <= community[h].day_epi_stop)
      {
      	 //printf("t=%d\n", t);
         r = t - community[h].day_epi_start;

         ptr_contact = (CONTACT *) malloc((size_t) sizeof(CONTACT));
         ptr_contact->next = NULL;
         ptr_contact->contact_mode = contact_mode;
         ptr_contact->offset = offset;

         if(community[h].contact_history[r].c2p_contact == NULL)
            community[h].contact_history[r].c2p_contact = ptr_contact;
         else
            community[h].contact_history[r].c2p_contact_rear->next = ptr_contact;
         community[h].contact_history[r].c2p_contact_rear = ptr_contact;
      }
   }   
}

void add_c2p_contact_history(int t, PEOPLE *person, int contact_mode, double offset)
{
   int h, i, j, k, l, m, n, r;
   int found, match;
   CONTACT *ptr_contact;

   h = person->community;
   r = t - community[h].day_epi_start;

   ptr_contact = (CONTACT *) malloc((size_t) sizeof(CONTACT));
   ptr_contact->next = NULL;
   ptr_contact->contact_mode = contact_mode;
   ptr_contact->offset = offset;

   if(person->contact_history[r].c2p_contact == NULL)
      person->contact_history[r].c2p_contact = ptr_contact;
   else
      person->contact_history[r].c2p_contact_rear->next = ptr_contact;
   person->contact_history[r].c2p_contact_rear = ptr_contact;
   
}

void add_c2p_risk_history(int t, PEOPLE *sus_person, int contact_mode, double offset)
{
   int h, i, j, k, l, m, n, r;
   int found, match;
   int n_c2p_covariate;
   double *c2p_covariate;
   RISK *ptr_risk;

   h = sus_person->community;
   r = t - community[h].day_epi_start;
   n_c2p_covariate = cfg_pars.n_c2p_covariate;

   c2p_covariate = NULL;
   n_c2p_covariate = cfg_pars.n_c2p_covariate;
   if(n_c2p_covariate > 0)  
   {
      c2p_covariate = (double *)malloc((size_t) (n_c2p_covariate * sizeof(double)));

      organize_c2p_covariate(t, sus_person, c2p_covariate);
   }
   
   found = 0;
   ptr_risk = sus_person->risk_history[r].c2p_risk;
   while(ptr_risk != NULL)
   {
      match = 1;
      for(l=0; l<n_c2p_covariate; l++)
      {
         if(ptr_risk->covariate[l] != c2p_covariate[l])
            match = 0;
      }

      if(ptr_risk->contact_mode != contact_mode)  match = 0;
      if(ptr_risk->offset != offset)  match = 0;
      if(match == 1)  { found = 1;  break;}
      ptr_risk = ptr_risk->next;
   }
   if(found == 1)
   {
      ptr_risk->size ++;
   }
   else
   {
      ptr_risk = (RISK *) malloc((size_t) sizeof(RISK));
      ptr_risk->next = NULL;
      ptr_risk->covariate = NULL;
      if(n_c2p_covariate > 0)
      {
         ptr_risk->covariate = (double *) malloc((size_t) n_c2p_covariate * sizeof(double));
         for(l=0; l<n_c2p_covariate; l++)
            ptr_risk->covariate[l] = c2p_covariate[l];
      }
      ptr_risk->contact_mode = contact_mode;
      ptr_risk->infective_prob = 1.0;
      ptr_risk->offset = offset;
      ptr_risk->size = 1;

      if(sus_person->risk_history[r].c2p_risk == NULL)
         sus_person->risk_history[r].c2p_risk = ptr_risk;
      else
         sus_person->risk_history[r].c2p_risk_rear->next = ptr_risk;
      sus_person->risk_history[r].c2p_risk_rear = ptr_risk;
   }
   free(c2p_covariate);
}

void add_c2p_risk_history_to_risk_class(int t, int h, int contact_mode, double offset)
{
   int i, j, k, l, m, n, r;
   int found, match;
   int n_c2p_covariate;
   double *c2p_covariate;
   RISK *ptr_risk;
   RISK_CLASS *ptr_class;
   PEOPLE *person;
   
   r = t - community[h].day_epi_start;
   c2p_covariate = NULL;
   n_c2p_covariate = cfg_pars.n_c2p_covariate;
   if(n_c2p_covariate > 0)  
      c2p_covariate = (double *)malloc((size_t) (n_c2p_covariate * sizeof(double)));

   ptr_class = community[h].risk_class;
   while(ptr_class != NULL)
   {   
      person = people + ptr_class->member[0].id;
         organize_c2p_covariate(t, person, c2p_covariate);
      found = 0;
      ptr_risk = ptr_class->risk_history[r].c2p_risk;
      while(ptr_risk != NULL)
      {
         match = 1;
         for(l=0; l<n_c2p_covariate; l++)
         {
            if(ptr_risk->covariate[l] != c2p_covariate[l])
               match = 0;
         }

         if(ptr_risk->contact_mode != contact_mode)  match = 0;
         if(ptr_risk->offset != offset)  match = 0;
         if(match == 1)  { found = 1;  break;}
         ptr_risk = ptr_risk->next;
      }
      if(found == 1)
         ptr_risk->size ++;
      else
      {
         ptr_risk = (RISK *) malloc((size_t) sizeof(RISK));
         ptr_risk->next = NULL;
         ptr_risk->covariate = NULL;
         if(n_c2p_covariate > 0)
         {
            ptr_risk->covariate = (double *) malloc((size_t) n_c2p_covariate * sizeof(double));
            for(l=0; l<n_c2p_covariate; l++)
               ptr_risk->covariate[l] = c2p_covariate[l];
         }
         ptr_risk->contact_mode = contact_mode;
         ptr_risk->infective_prob = 1.0;
         ptr_risk->offset = offset;
         ptr_risk->size = 1;

         if(ptr_class->risk_history[r].c2p_risk == NULL)
            ptr_class->risk_history[r].c2p_risk = ptr_risk;
         else
            ptr_class->risk_history[r].c2p_risk_rear->next = ptr_risk;
         ptr_class->risk_history[r].c2p_risk_rear = ptr_risk;
      }
      ptr_class = ptr_class->next;
   }  
   free(c2p_covariate);
}

// IMPORTANT: this function adds a p2p contact node to the contact history of
// inf_person, not sus_person, which means that the covariates in each contact
// node are ordered as if the host (to whom the node is linked) is an infective.
// This is different from the risk nodes inwhich the covariates are ordered
// as if the host of the node is a susceptible.
// The reason for this special order is that, when build risk history,
// we scan the infectious period of each case for his contact hisotry,
// and then use the contact node of the case to build the risk node of the susceptible
// indicated by that contact node.
void add_p2p_contact_history_to_community(int h, int start_day, int stop_day,
                             int contact_mode, double offset)
{
   int i, j, k, l, m, n, r, t;
   int found, match;
   int condition;
   CONTACT *ptr_contact;

   for(t = start_day; t <= stop_day; t++)
   {
      if(t >= community[h].day_epi_start && t <= community[h].day_epi_stop)
      {
         r = t - community[h].day_epi_start;
         for(k=0; k<community[h].size; k++)
         {
            i = community[h].member[k];
            // if no imputation, record contacts only for infected people; otherwise, recorde all of them.
            if(cfg_pars.EM == 1 || cfg_pars.simulation == 1)  condition = (people[i].ignore != 1 && t <= people[i].day_exit);
            else condition = (people[i].ignore != 1 && people[i].infection == 1 &&
                              people[i].day_infective_lower != MISSING && people[i].day_infective_upper != MISSING &&
                              t>=people[i].day_infective_lower && t<=people[i].day_infective_upper && t <= people[i].day_exit); 
            
            if(condition)
            {
               ptr_contact = (CONTACT *) malloc((size_t) sizeof(CONTACT));
               ptr_contact->next = NULL;
               ptr_contact->pair = NULL;
               ptr_contact->contact_mode = contact_mode;
               ptr_contact->contact_id = people[i].id;
               ptr_contact->offset = offset;

               if(community[h].contact_history[r].p2p_contact == NULL)
                  community[h].contact_history[r].p2p_contact = ptr_contact;
               else
                  community[h].contact_history[r].p2p_contact_rear->next = ptr_contact;
               community[h].contact_history[r].p2p_contact_rear = ptr_contact;
            }
         }
      }
   }         
}

CONTACT * add_p2p_contact_history(int t, PEOPLE *sus_person, PEOPLE *inf_person,
                             int contact_mode, double offset)
{
   int h, i, j, k, l, m, n, r;
   int found, match;
   CONTACT *ptr_contact;

   h = sus_person->community;
   r = t - community[h].day_epi_start;

   ptr_contact = (CONTACT *) malloc((size_t) sizeof(CONTACT));
   ptr_contact->next = NULL;
   ptr_contact->pair = NULL;
   ptr_contact->contact_mode = contact_mode;
   ptr_contact->contact_id = inf_person->id;
   ptr_contact->offset = offset;

   if(sus_person->contact_history[r].p2p_contact == NULL)
      sus_person->contact_history[r].p2p_contact = ptr_contact;
   else
      sus_person->contact_history[r].p2p_contact_rear->next = ptr_contact;
   sus_person->contact_history[r].p2p_contact_rear = ptr_contact;
   return(ptr_contact);
}

// NOTE: the pointer ptr_contact points to a contact node associated with an infective to whom sus_person
// was exposed to. 
void add_p2p_risk_history(int t, PEOPLE *sus_person, PEOPLE *inf_person, int contact_mode, double offset, double infective_prob, int symptom)
{
   int h, i, j, k, l, m, n, r;
   int found, match;
   int n_p2p_covariate;
   double *p2p_covariate;
   RISK *ptr_risk;

   h = sus_person->community;
   r = t - community[h].day_epi_start;
   
   p2p_covariate = NULL;
   n_p2p_covariate = cfg_pars.n_p2p_covariate;
   if(n_p2p_covariate > 0)  
   {
      p2p_covariate = (double *)malloc((size_t) (n_p2p_covariate * sizeof(double)));

      organize_p2p_covariate(t, sus_person, inf_person, p2p_covariate);
   }   
   // we are scanning the contact history of the infective person, but we are building the risk history
   // for the susceptible person; therefore, we need to
   //find the contact node of the susceptible person corresonding to the contact node of the infective person  

   found = 0;
   ptr_risk = sus_person->risk_history[r].p2p_risk;
   while(ptr_risk != NULL)
   {
      match = 1;
      for(l=0; l<n_p2p_covariate; l++)
      {
         if(ptr_risk->covariate[l] != p2p_covariate[l])
            match = 0;
      }

      if(ptr_risk->contact_mode != contact_mode)  match = 0;
      if(ptr_risk->infective_prob != infective_prob)  match = 0;
      if(ptr_risk->symptom != symptom)  match = 0;
      if(ptr_risk->offset != offset)  match = 0;
      if(match == 1)  { found = 1;  break;}
      ptr_risk = ptr_risk->next;
   }
   if(found == 1)
   {
      ptr_risk->size += inf_person->weight;
   }
   else
   {
      ptr_risk = (RISK *) malloc((size_t) sizeof(RISK));
      ptr_risk->next = NULL;
      ptr_risk->covariate = NULL;
      if(n_p2p_covariate > 0)
      {
         ptr_risk->covariate = (double *) malloc((size_t) n_p2p_covariate * sizeof(double));
         for(l=0; l<n_p2p_covariate; l++)
            ptr_risk->covariate[l] = p2p_covariate[l];
      }
      ptr_risk->contact_mode = contact_mode;
      ptr_risk->infective_prob = infective_prob;
      ptr_risk->symptom = symptom;
      ptr_risk->offset = offset;
      ptr_risk->size = inf_person->weight;

      if(sus_person->risk_history[r].p2p_risk == NULL)
         sus_person->risk_history[r].p2p_risk = ptr_risk;
      else
         sus_person->risk_history[r].p2p_risk_rear->next = ptr_risk;
      sus_person->risk_history[r].p2p_risk_rear = ptr_risk;
   }
   free(p2p_covariate);
}

void add_p2p_risk_history_to_risk_class(int t, PEOPLE *inf_person, int contact_mode, double offset, double infective_prob, int symptom)
{
   int h, i, j, k, l, m, n, r;
   int found, match;
   int n_p2p_covariate;
   double *p2p_covariate;
   RISK *ptr_risk;
   RISK_CLASS *ptr_class;
   PEOPLE *sus_person;

   h = inf_person->community;
   r = t - community[h].day_epi_start;
   p2p_covariate = NULL;
   n_p2p_covariate = cfg_pars.n_p2p_covariate;
   if(n_p2p_covariate > 0)  
      p2p_covariate = (double *)malloc((size_t) (n_p2p_covariate * sizeof(double)));

   // if inf_person is a contact for the community, add him or her to the risk history for each risk class in the community
   ptr_class = community[h].risk_class;
   while(ptr_class != NULL)
   {
      sus_person = people + ptr_class->member[0].id;
      organize_p2p_covariate(t, sus_person, inf_person, p2p_covariate);
   
      // we are scanning the contact history of the infective person, but we are building the risk history
      // for the susceptible person; therefore, we need to
      //find the contact node of the susceptible person corresonding to the contact node of the infective person  

      found = 0;
      ptr_risk = ptr_class->risk_history[r].p2p_risk;
      while(ptr_risk != NULL)
      {
         match = 1;
         for(l=0; l<n_p2p_covariate; l++)
         {
            if(ptr_risk->covariate[l] != p2p_covariate[l])
               match = 0;
         }

         if(ptr_risk->contact_mode != contact_mode)  match = 0;
         if(ptr_risk->infective_prob != infective_prob)  match = 0;
         if(ptr_risk->symptom != symptom)  match = 0;
         if(ptr_risk->offset != offset)  match = 0;
         if(match == 1)  { found = 1;  break;}
         ptr_risk = ptr_risk->next;
      }
      if(found == 1)
         ptr_risk->size += inf_person->weight;
      else
      {
         ptr_risk = (RISK *) malloc((size_t) sizeof(RISK));
         ptr_risk->next = NULL;
         ptr_risk->covariate = NULL;
         if(n_p2p_covariate > 0)
         {
            ptr_risk->covariate = (double *) malloc((size_t) n_p2p_covariate * sizeof(double));
            for(l=0; l<n_p2p_covariate; l++)
               ptr_risk->covariate[l] = p2p_covariate[l];
         }
         ptr_risk->contact_mode = contact_mode;
         ptr_risk->infective_prob = infective_prob;
         ptr_risk->symptom = symptom;
         ptr_risk->offset = offset;
         ptr_risk->size = inf_person->weight;

         if(ptr_class->risk_history[r].p2p_risk == NULL)
            ptr_class->risk_history[r].p2p_risk = ptr_risk;
         else
            ptr_class->risk_history[r].p2p_risk_rear->next = ptr_risk;
         ptr_class->risk_history[r].p2p_risk_rear = ptr_risk;
      }
      ptr_class = ptr_class->next;
   }   
   free(p2p_covariate);
}

// This function creates risk classes within each community. A risk class
// is composed of members whose covariates are identical so that they share the same risk
// given that they also share the same contact history. If people in the same community
// do not share contact history, or none of them share covariates, each individual is considered a class. 
// Risk classes are built in communities regardless of whether people in the same community share contact history or not.
// If they do share common contact history, then a risk history array is created at the community level.
// If not, then an individual level risk history is created in the core.h file, not here
void create_risk_class(void)
{
   int h, i, j, k, l, m, n, r, t;
   int *match, identical, n_match, n_idx, n_valid;
   int n_c2p_covariate, n_sus_p2p_covariate, n_int_p2p_covariate;
   int n_time_ind_covariate, n_time_dep_covariate, n_covariate;
   double *i_covariate, *j_covariate;
   RISK_CLASS *ptr_class;
   INTEGER_CHAIN *ptr_integer;

   n_c2p_covariate = cfg_pars.n_c2p_covariate;
   n_sus_p2p_covariate = cfg_pars.n_sus_p2p_covariate;
   n_int_p2p_covariate = cfg_pars.n_int_p2p_covariate; 
   n_time_ind_covariate = cfg_pars.n_time_ind_covariate;
   n_time_dep_covariate = cfg_pars.n_time_dep_covariate;
   n_covariate = n_time_ind_covariate + n_time_dep_covariate;

   i_covariate = (double *)malloc((size_t) (n_covariate * sizeof(double)));
   j_covariate = (double *)malloc((size_t) (n_covariate * sizeof(double)));

   for(h=0; h<n_community; h++)
   if(community[h].ignore == 0)
   {
      make_1d_array_int(&match, community[h].size, 0);
      n_match = 0;
      
      n_valid = community[h].size;

      for(m=0; m<community[h].size; m++)
      if(match[m] == 0)
      {
         i = community[h].member[m];
         if(people[i].ignore == 0)
         {
            match[m] = 1;
            n_match++;
            if(community[h].risk_class == NULL)
            {
               ptr_class = (RISK_CLASS *) malloc((size_t) sizeof(RISK_CLASS));
               community[h].risk_class = ptr_class; 
            }
            else
            {
               ptr_class = community[h].risk_class_rear;
               ptr_class->next = (RISK_CLASS *) malloc((size_t) sizeof(RISK_CLASS));
               ptr_class = ptr_class->next;
            }
            community[h].risk_class_rear = ptr_class;
            ptr_class->next = NULL;
            ptr_integer = (INTEGER_CHAIN *) malloc((size_t) sizeof(INTEGER_CHAIN));
            ptr_class->member = ptr_integer;
            ptr_integer->id = i;
            ptr_integer->next = NULL;
            ptr_class->member_rear = ptr_integer;
            ptr_class->risk_history = NULL;
            ptr_class->size = 1;
            people[i].risk_class = ptr_class;
            if(cfg_pars.common_contact_history_within_community == 1)
            {
               // risk hitroy is stored in risk class if community share common contact histry
               ptr_class->risk_history = (RISK_HISTORY *)malloc((size_t) (community[h].epi_duration * sizeof(RISK_HISTORY)));
               for(t=community[h].day_epi_start; t<=community[h].day_epi_stop; t++)
               {
                  r = t - community[h].day_epi_start;
                  ptr_class->risk_history[r].c2p_risk = NULL;
                  ptr_class->risk_history[r].c2p_risk_rear = NULL;
                  ptr_class->risk_history[r].p2p_risk = NULL;
                  ptr_class->risk_history[r].p2p_risk_rear = NULL;
               }
               // find members belonging to the same risk class in the community
               for(n=m+1; n<community[h].size; n++)
               if(match[n] == 0)
               {
                  j = community[h].member[n];
                  if(people[j].ignore == 0)
                  {
                     identical = 1;
                     for(k=0; k<n_time_ind_covariate; k++)
                     {
                        i_covariate[k] = people[i].time_ind_covariate[k];
                        j_covariate[k] = people[j].time_ind_covariate[k];
                     }
                     if(cfg_pars.n_time_dep_covariate == 0)
                     {
                        // compare covariates affecting c2p risk
                        for(k=0; k<n_c2p_covariate; k++)
                        {
                           l = cfg_pars.c2p_covariate[k] - 1;
                           if(i_covariate[l] != j_covariate[l])  
                           {  
                              identical = 0;
                              break;
                           }
                        }
                     
                        // compare covariates affecting p2p risk on the susceptible side
                        for(k=0; k<n_sus_p2p_covariate; k++)
                        {
                           l = cfg_pars.sus_p2p_covariate[k] - 1;
                           if(i_covariate[l] != j_covariate[l])  
                           {  
                              identical = 0;
                              break;
                           }
                        }
                        // compare covariates affecting p2p risk via interaction on the susceptible side
                        for(k=0; k<n_int_p2p_covariate; k++)
                        {
                           l = cfg_pars.interaction[k][0] - 1;
                           if(i_covariate[l] != j_covariate[l])  
                           {  
                              identical = 0;
                              break;
                           }
                        }
                     }
                     else
                     {
                        for(t=community[h].day_epi_start; t<= community[h].day_epi_stop; t++)
                        {
                           r = t - community[h].day_epi_start;
                           for(k=0; k<n_time_dep_covariate; k++)
                           {
                              i_covariate[n_time_ind_covariate + k] = people[i].time_dep_covariate[r][k];
                              j_covariate[n_time_ind_covariate + k] = people[j].time_dep_covariate[r][k];
                           }
                           // compare covariates affecting c2p risk
                           for(k=0; k<n_c2p_covariate; k++)
                           {
                              l = cfg_pars.c2p_covariate[k] - 1;
                              if(i_covariate[l] != j_covariate[l])  
                              {  
                                 identical = 0;
                                 break;
                              }
                           }
                        
                           // compare covariates affecting p2p risk on the susceptible side
                           for(k=0; k<n_sus_p2p_covariate; k++)
                           {
                              l = cfg_pars.sus_p2p_covariate[k] - 1;
                              if(i_covariate[l] != j_covariate[l])  
                              {  
                                 identical = 0;
                                 break;
                              }
                           }
                           // compare covariates affecting p2p risk via interaction on the susceptible side
                           for(k=0; k<n_int_p2p_covariate; k++)
                           {
                              l = cfg_pars.interaction[k][0] - 1;
                              if(i_covariate[l] != j_covariate[l])  
                              {  
                                 identical = 0;
                                 break;
                              }
                           }
                           if(identical == 0)  break;
                        }
                     }
                     if(identical == 1)
                     {
                        match[n] = 1;
                        n_match++;
                        ptr_integer = ptr_class->member_rear;
                        ptr_integer->next = (INTEGER_CHAIN *) malloc((size_t) sizeof(INTEGER_CHAIN));
                        ptr_integer = ptr_integer->next;
                        ptr_integer->id = j;
                        ptr_integer->next = NULL;
                        ptr_class->member_rear = ptr_integer;
                        ptr_class->size ++;
                        people[j].risk_class = ptr_class;
                     }
                  }
               }
            }
         }
         if(n_match >= n_valid)  break;
      }
      free(match);
   }
   free(i_covariate);
   free(j_covariate);
}
           
            
int delete_p2p_risk_history(int t, PEOPLE *sus_person, CONTACT *sus_ptr_contact, double infective_prob, int symptom)
{
   int h, i, j, k, l, m, n, r;
   int found, match;
   int n_p2p_covariate;
   double *p2p_covariate;
   PEOPLE *inf_person;
   RISK *ptr_risk, *ptr_prior, *ptr_next;

   h = sus_person->community;
   r = t - community[h].day_epi_start;
   n_p2p_covariate = cfg_pars.n_p2p_covariate;

   inf_person = people + sus_ptr_contact->contact_id;
   p2p_covariate = NULL;
   if(n_p2p_covariate > 0)  
   {
      p2p_covariate = (double *)malloc((size_t) (n_p2p_covariate * sizeof(double)));

      organize_p2p_covariate(t, sus_person, inf_person, p2p_covariate);
   }   

   found = 0;
   ptr_risk = sus_person->risk_history[r].p2p_risk;
   if(ptr_risk != NULL)
   {
      ptr_prior = NULL;
      ptr_next = ptr_risk->next;
      while(ptr_risk != NULL)
      {
         match = 1;
         for(l=0; l<n_p2p_covariate; l++)
         {
            if(ptr_risk->covariate[l] != p2p_covariate[l])
               match = 0;
         }

         if(ptr_risk->contact_mode != sus_ptr_contact->contact_mode)  match = 0;
         if(ptr_risk->infective_prob != infective_prob)  match = 0;
         if(ptr_risk->symptom != symptom)  match = 0;
         if(ptr_risk->offset != sus_ptr_contact->offset)  match = 0;
         if(match == 1)  { found = 1;  break;}
         ptr_prior = ptr_risk;
         ptr_risk = ptr_risk->next;
         ptr_next = ptr_risk->next;
      }
      if(found == 1)
      {
         ptr_risk->size -= inf_person->weight;
         //this node is now empty and should be deleted
         if(ptr_risk->size == 0)
         {
            //if the head node is to be deleted
            if(ptr_risk == sus_person->risk_history[r].p2p_risk)
               sus_person->risk_history[r].p2p_risk = ptr_next;
            else
               ptr_prior->next = ptr_next;
            //if the rear node is to be deleted
            if(ptr_risk == sus_person->risk_history[r].p2p_risk_rear)
               sus_person->risk_history[r].p2p_risk_rear = ptr_prior;
           
            if(ptr_risk->covariate != NULL)  free(ptr_risk->covariate);
            free(ptr_risk);
         }
      }
   }
   free(p2p_covariate);
   return(0);
}

int delete_p2p_risk_history_in_risk_class(int t, PEOPLE *inf_person, double infective_prob, int symptom)
{
   int h, i, j, k, l, m, n, r;
   int found, match;
   int n_p2p_covariate;
   double *p2p_covariate;
   PEOPLE *sus_person;
   RISK *ptr_risk, *ptr_risk_prior, *ptr_risk_next;
   RISK_CLASS *ptr_class;
   CONTACT *ptr_contact; 

   p2p_covariate = NULL;
   n_p2p_covariate = cfg_pars.n_p2p_covariate;
   if(n_p2p_covariate > 0)  
      p2p_covariate = (double *)malloc((size_t) (n_p2p_covariate * sizeof(double)));

   h = inf_person->community;
   r = t - community[h].day_epi_start;
   // loop over contact first. There could be multiple types of contact mode for the same person
   ptr_contact = community[h].contact_history[r].p2p_contact;
   while(ptr_contact != NULL)
   {
      if(ptr_contact->contact_id == inf_person->id)  
      {
         ptr_class = community[h].risk_class;
         while(ptr_class != NULL)
         {
            sus_person = people + ptr_class->member[0].id;
            if(n_p2p_covariate > 0)  organize_p2p_covariate(t, sus_person, inf_person, p2p_covariate); 
            found = 0;
            ptr_risk = ptr_class->risk_history[r].p2p_risk;
            if(ptr_risk != NULL)
            {
               ptr_risk_prior = NULL;
               ptr_risk_next = ptr_risk->next;
               while(ptr_risk != NULL)
               {
                  match = 1;
                  for(l=0; l<n_p2p_covariate; l++)
                  {
                     if(ptr_risk->covariate[l] != p2p_covariate[l])
                        match = 0;
                  }
                  if(ptr_risk->contact_mode != ptr_contact->contact_mode)  match = 0;
                  if(ptr_risk->infective_prob != infective_prob)  match = 0;
                  if(ptr_risk->symptom != symptom)  match = 0;
                  if(ptr_risk->offset != ptr_contact->offset)  match = 0;
                  if(match == 1)  { found = 1;  break;}
                  ptr_risk_prior = ptr_risk;
                  ptr_risk = ptr_risk->next;
                  ptr_risk_next = ptr_risk->next;
               }
            }   
            if(found == 1)
            {
               ptr_risk->size -= inf_person->weight;
               //this node is now empty and should be deleted
               if(ptr_risk->size == 0)
               {
                  //if the head node is to be deleted
                  if(ptr_risk == ptr_class->risk_history[r].p2p_risk)
                     ptr_class->risk_history[r].p2p_risk = ptr_risk_next;
                  else
                     ptr_risk_prior->next = ptr_risk_next;
                  //if the rear node is to be deleted
                  if(ptr_risk == ptr_class->risk_history[r].p2p_risk_rear)
                     ptr_class->risk_history[r].p2p_risk_rear = ptr_risk_prior;
                  if(ptr_risk->covariate != NULL)  free(ptr_risk->covariate);
                  free(ptr_risk);
               }
            }
            ptr_class = ptr_class->next;
         }
      }
      ptr_contact = ptr_contact->next;
   }
   free(p2p_covariate);
}

int centered_discrete(int length, double *log_prob)
{
   int j, k;
   double min_log_prob, max_log_prob, *prob;

   max_log_prob = log_prob[0];
   for(j = 0; j < length; j++)
   {
      max_log_prob = max(max_log_prob, log_prob[j]);
   }
   prob = (double *) malloc((size_t) length * sizeof(double));
   for(j = 0; j < length; j++)
   {
      prob[j] = exp(log_prob[j] - max_log_prob);
   }
   k = discrete(length, prob, &seed);
   free(prob);
   return(k);
}

// check individual-level c2p risk history
void show_c2p_risk(PEOPLE *person, int start, int stop)
{
   int h, l, r, t;
   RISK *ptr_risk;

   printf("c2p risk for person %d during days %d - %d \n", person->id, start, stop);
   h = person->community;
   for(t=start; t<=stop; t++)
   {
      r = t - community[h].day_epi_start;
      {
         ptr_risk = person->risk_history[r].c2p_risk;
         if(ptr_risk == NULL)  printf("No c2p risk on day %d for person %d\n", t, person->id);
         while(ptr_risk != NULL)
         {
            printf("i=%d  day=%d  contact_mode=%d  offset=%f  size=%d  covariate:",
                  person->id, t, ptr_risk->contact_mode, 
                  ptr_risk->offset, ptr_risk->size);
            for(l=0; l<cfg_pars.n_c2p_covariate; l++)
              printf("%f  ", ptr_risk->covariate[l]);
            printf("\n");
            ptr_risk = ptr_risk->next;
         }
      }
   }
}
// check individual-level p2p risk history
void show_p2p_risk(PEOPLE *person, int start, int stop)
{
   int h, l, r, t;
   RISK *ptr_risk;

   printf("p2p risk for person %d during days %d - %d \n", person->id, start, stop);
   h = person->community;
   for(t=start; t<=stop; t++)
   {
      r = t - community[h].day_epi_start;
      {
         ptr_risk = person->risk_history[r].p2p_risk;
         if(ptr_risk == NULL)  printf("No p2p risk on day %d for person %d\n", t, person->id);
         while(ptr_risk != NULL)
         {
            printf("i=%d  day=%d  contact_mode=%d  symptom=%d  infective_prob=%e offset=%f  size=%d  covariate:",
                  person->id, t, ptr_risk->contact_mode, ptr_risk->symptom, ptr_risk->infective_prob,
                  ptr_risk->offset, ptr_risk->size);
            for(l=0; l<cfg_pars.n_p2p_covariate; l++)
              printf("%f  ", ptr_risk->covariate[l]);
            printf("\n");
            ptr_risk = ptr_risk->next;
         }
      }
   }
}

// check risk-class-level c2p risk history
void show_c2p_risk_in_risk_class(int h, int start, int stop)
{
   int l, r, t;
   RISK *ptr_risk;
   RISK_CLASS *ptr_class;
   
   printf("c2p risk for commuity %d during days %d - %d \n", h, start, stop);
   ptr_class = community[h].risk_class;
   while(ptr_class != NULL)
   {
      for(t=start; t<=stop; t++)
      {
         r = t - community[h].day_epi_start;
         ptr_risk = ptr_class->risk_history[r].c2p_risk;
         if(ptr_risk == NULL)  printf("No c2p risk on day %d for community %d\n", t, h);
         while(ptr_risk != NULL)
         {
            printf("day=%d  contact_mode=%d  offset=%f  size=%d  covariate:",
                   t, ptr_risk->contact_mode, ptr_risk->offset, ptr_risk->size);
            for(l=0; l<cfg_pars.n_c2p_covariate; l++)
               printf("%f  ", ptr_risk->covariate[l]);
            printf("\n");
            ptr_risk = ptr_risk->next;
         }
      }
      ptr_class = ptr_class->next;
   }
}

// check risk-class-level p2p risk history
void show_p2p_risk_in_risk_class(int h, int start, int stop)
{
   int l, r, t;
   RISK *ptr_risk;
   RISK_CLASS *ptr_class;
   
   printf("p2p risk for commuity %d during days %d - %d \n", h, start, stop);
   ptr_class = community[h].risk_class;
   while(ptr_class != NULL)
   {
      printf("h=%d  size=%d\n", h, ptr_class->size);
      printf("start=%d  stop=%d\n", start, stop);
      for(t=start; t<=stop; t++)
      {
         r = t - community[h].day_epi_start;
         ptr_risk = ptr_class->risk_history[r].p2p_risk;
         if(ptr_risk == NULL)  printf("No p2p risk on day %d for community %d\n", t, h);
         while(ptr_risk != NULL)
         {
            printf("day=%d  contact_mode=%d  symptom=%d  infective_prob=%e offset=%f  size=%d  covariate:",
                   t, ptr_risk->contact_mode, ptr_risk->symptom, ptr_risk->infective_prob,
                   ptr_risk->offset, ptr_risk->size);
            for(l=0; l<cfg_pars.n_p2p_covariate; l++)
               printf("%f  ", ptr_risk->covariate[l]);
            printf("\n");
            ptr_risk = ptr_risk->next;
         }
      }
      ptr_class = ptr_class->next;
   }
}


void show_community_states(COMMUNITY *community)
{
   int i, j;

   for(j=0; j < community->size; j++)
   {
      i = community->member[j];
      printf("%8d  %8d  %1d  %1d  %1d  %8d  %1d\n",
            people[i].id, people[i].community,
            people[i].pre_immune, people[i].infection,
            people[i].symptom, people[i].day_ill, people[i].idx);
   }
}

void show_par_effective(double *par_effective, FILE *file)
{
   int i, j, m;
   if(file == NULL)
   {
      for(j=0; j<cfg_pars.n_par_equiclass; j++)
      {
         m = cfg_pars.par_equiclass[j].member[0] - 1;
         if(m < cfg_pars.n_b_mode + cfg_pars.n_p_mode + cfg_pars.n_u_mode + cfg_pars.n_q_mode)
            printf("%e  %e\n", par_effective[j], inv_logit(par_effective[j]));
         else printf("%e  %e\n", par_effective[j], exp(par_effective[j]));
      }
   }
   else
   {
      fprintf(file, "\n");
      for(j=0; j<cfg_pars.n_par_equiclass; j++)
      {
         m = cfg_pars.par_equiclass[j].member[0] - 1;
         if(m < cfg_pars.n_b_mode + cfg_pars.n_p_mode + cfg_pars.n_u_mode + cfg_pars.n_q_mode)
            fprintf(file, "%e  %e\n", par_effective[j], inv_logit(par_effective[j]));
         else fprintf(file, "%e  %e\n", par_effective[j], exp(par_effective[j]));
      }
      fprintf(file, "\n");
   }
}

/* this function check the admissibility of the model. By admissibility we mean that there should be some information about
each parameter. If there is no information for any parameter, the information matrix will not be inversable and will
cause error messages, e.g, things like log(0) or x/0 or likelihood=0 will occur. Pass of this admissibility checking does not 
guaranttee estimability of the parameters. Please note that without any information is not equivalent to without sufficient 
information. Without sufficient information, the transmission probabilities may not be estimable (e.g., go to 0), but we do 
not expect error messages from this situation, because we enforce a lower bound, close_to_0, for any estimate. 
Plus, there is no way to pre-check the sufficiency of information. 
If you considers asymptomatic infection, checking admissibility may not be informative.*/

int admissibility()
{
   int h, i, j, k, l, m, n, t, r;
   int error;
   int b_mode, p_mode;
   int n_b_mode, n_p_mode;
   int *b_admissible, *p_admissible, n_impossible_case;
   int *b_fixed, *p_fixed;
   int *n_case, **cases;

   CONTACT *ptr_contact;
   FILE *file;

   n_case = (int *) malloc((size_t) (n_community * sizeof(int)));

   for(h=0; h<n_community; h++)
      n_case[h] = 0;

   for(h=0; h<n_community; h++)
   {
      for(j=0; j<community[h].size; j++)
      {
         i = community[h].member[j];
         if(people[i].ignore == 0 && people[i].infection == 1 && people[i].symptom == 1)
            n_case[h]++;
      }
   }

   
   cases = (int **)malloc((size_t) (n_community * sizeof(int *)));
   for(h=0; h<n_community; h++)
   {
      cases[h] = NULL;
      if(n_case[h] > 0)  cases[h] = (int *) malloc((size_t) (n_case[h] * sizeof(int)));
   }

   for(h=0; h<n_community; h++)
   {
      l = 0;
      for(j=0; j<community[h].size; j++)
      {
         i = community[h].member[j];
         if(people[i].ignore == 0 && people[i].infection == 1 && people[i].symptom == 1)
            cases[h][l++] = i;
      }
   }

   n_b_mode = cfg_pars.n_b_mode;
   n_p_mode = cfg_pars.n_p_mode;

   if(n_b_mode > 0)
   {
      b_admissible = (int *)malloc((size_t) (n_b_mode * sizeof(int)));
      b_fixed = (int *)malloc((size_t) (n_b_mode * sizeof(int)));
   }
   if(n_p_mode > 0)
   {
      p_admissible = (int *)malloc((size_t) (n_p_mode * sizeof(int)));
      p_fixed = (int *)malloc((size_t) (n_p_mode * sizeof(int)));
   }

   for(i=0; i<n_b_mode; i++)
   {
       b_admissible[i] = 0;
       b_fixed[i] = 0;
   }
   for(i=0; i<n_p_mode; i++)
   {
      p_admissible[i] = 0;
      p_fixed[i] = 0;
   }

   //find out which b or p is fixed, i.e., not to be estimated
   for(i=0; i<cfg_pars.n_par_fixed; i++)
   {
      j = cfg_pars.par_fixed_id[i] - 1;
      if(j < n_b_mode)  b_fixed[j] = 1;
      else if(j < n_b_mode + n_p_mode)  p_fixed[j-n_b_mode] = 1;
   }

   n_impossible_case = 0;
   for(h=0; h<n_community; h++)
   {
      for(k=0; k<n_case[h]; k++)
      {
         i = cases[h][k];
         people[i].any_exposure = -1;

         if(!(cfg_pars.adjust_for_left_truncation == 1 && people[i].idx == 1))
         {
            //printf("id=%d  infeciton=%d  day_ill=%d  day_infection_lower=%d  day_infection_upper=%d\n", i, people[i].infection, people[i].day_ill, people[i].day_infection_lower, people[i].day_infection_upper);
            people[i].any_exposure = 0;

            for(t=people[i].day_infection_lower; t<=people[i].day_infection_upper; t++)
            {
               r = t - community[h].day_epi_start;
               if(cfg_pars.common_contact_history_within_community == 0)
                  ptr_contact = people[i].contact_history[r].c2p_contact;
               else   
                  ptr_contact = community[h].contact_history[r].c2p_contact;
               /*
               if(i == 1)
               {
                  printf("t=%d null contact=%d\n", t, (ptr_contact == NULL));
                  printf("j=%d contact_mode=%d\n", j, ptr_contact->contact_mode);
                  exit(0);
               }*/
               while(ptr_contact != NULL)
               {
                  b_mode = ptr_contact->contact_mode;
                  b_admissible[b_mode] = 1;
                  people[i].any_exposure = 1;
                  ptr_contact = ptr_contact->next;
               }
               if(cfg_pars.common_contact_history_within_community == 0)
                  ptr_contact = people[i].contact_history[r].p2p_contact;
               else   
                  ptr_contact = community[h].contact_history[r].p2p_contact;
               while(ptr_contact != NULL)
               {
                  p_mode = ptr_contact->contact_mode;
                  p_admissible[p_mode] = 1;
                  people[i].any_exposure = 1;
                  ptr_contact = ptr_contact->next;
               }
            }
            if(people[i].any_exposure != 1)  n_impossible_case ++;
         }
      }
   }

   error = 0;
   if(n_impossible_case > 0) error = 1;
   for(i=0; i<n_b_mode; i++)
      if(b_admissible[i] == 0 && b_fixed[i] == 0)  error = 1;

   for(i=0; i<n_p_mode; i++)
      if(p_admissible[i] == 0 && p_fixed[i] == 0)  error = 1;

   if(error == 1)
   {
      if(cfg_pars.silent_run == 0)  printf("The model is not admissible for the given data\n");
      if(cfg_pars.write_error_log == 1)
      {
         if((file = fopen("errors.txt", "a")) == NULL)
            file = fopen("errors.txt", "w");

         fprintf(file, "\n#1.1 Admissibility error: %d impossible cases\n", n_impossible_case);
         if(n_impossible_case > 0)
         {
            for(h=0; h<n_community; h++)
            {
               for(i=0; i<n_case[h]; i++)
               {
                  j = cases[h][i];
                  if(people[j].any_exposure == 0)
                  {
                     fprintf(file, "community:%d  person:%d\n", community[h].id, people[j].id);
                  }
               }
            }
         }
         fprintf(file, "\n#1.2 Admissibility error: Is there information about c2p transmission?\n");
         for(k=0; k<n_b_mode; k++)
         {
            if(b_admissible[k] == 1) fprintf(file, "%d  YES\n", k);
            else fprintf(file, "%d  NO\n", k);
         }
         fprintf(file, "\n#1.3 Admissibility error: Is  there information about p2p transmission?\n");
         for(k=0; k<n_p_mode; k++)
         {
            if(p_admissible[k] == 1) fprintf(file, "%d  YES\n", k);
            else fprintf(file, "%d  NO\n", k);
         }
         fclose(file);
      }   
   }
   else  if(cfg_pars.silent_run == 0)  printf("The model is admissible for the given data\n");


   if(n_b_mode > 0)  {free(b_admissible); free(b_fixed);}
   if(n_p_mode > 0)  {free(p_admissible); free(p_fixed);}
   free(n_case);
   if(cases!= NULL)
   {
      for(h=0; h<n_community; h++)
      {
         free(cases[h]);
      }
      free(cases);
   }

   return(error);
}

// create a chain of size n, each node has an array of size m to record the possible state combination of 
// the m members in the family who have multiple possible states. 
void create_state_chain(int n,  int m, STATE **head, STATE **rear)
{
   int i, j, k, l;
   int n_sample, size_impute;
   STATE *ptr_stat;

   ptr_stat = NULL;
   for(i=0; i<n; i++)
   {
      if(ptr_stat == NULL)
      {
         ptr_stat = (STATE *) malloc((size_t) sizeof(STATE));
         (* head) = ptr_stat;
      }
      else
      {
         ptr_stat->next = (STATE *) malloc((size_t) sizeof(STATE));
         ptr_stat = ptr_stat->next;
      }
      ptr_stat->states = (int *) malloc((size_t) m * sizeof(int));
      for(j=0; j<m; j++) ptr_stat->states[j] = 0;

      ptr_stat->log_L_ini = 0;
      ptr_stat->log_L_cur = 0;
      ptr_stat->log_L_pre = 0;
      ptr_stat->ratio = 1.0;
      ptr_stat->next = NULL;

   }
   (* rear) = ptr_stat;
}

void free_state_chain(STATE *head, STATE *rear)
{
   int i=0;
   STATE *ptr_stat, *ptr2_stat;

   ptr_stat = head;
   while(ptr_stat != rear)
   {
      ptr2_stat = ptr_stat->next;
      free(ptr_stat->states);
      free(ptr_stat);
      ptr_stat = ptr2_stat;
   }
   if(rear != NULL)
   {
      free(rear->states);
      free(rear);
   }
}

void restore_susceptibility()
{
   int h, i, j, k, r, s, t;
   RISK *ptr_risk, *ptr2_risk;
   STATE *ptr_stat, *ptr2_stat;
   RISK_CLASS *ptr_class;
   
   if(people != NULL)
   {
      for(i=0; i<p_size; i++)
      {
         people[i].ignore = 0;
         if(cfg_pars.preset_index == 0)
         {
            people[i].idx = 0;
            people[i].infection = 0;
            people[i].day_infection = MISSING;
            people[i].day_infection_lower = people[i].day_infection_upper = MISSING;
            people[i].symptom = 0;
            people[i].day_ill = MISSING;
            people[i].day_infective_lower = people[i].day_infective_upper = MISSING;
         }
         else
         {
            if(people[i].idx == 0)
            {
               people[i].infection = 0;
               people[i].day_infection = MISSING;
               people[i].day_infection_lower = people[i].day_infection_upper = MISSING;
               people[i].symptom = 0;
               people[i].day_ill = MISSING;
               people[i].day_infective_lower = people[i].day_infective_upper = MISSING;
            }
         }
               
         if(people[i].possible_states != NULL)
         {
            free_2d_array_int(people[i].possible_states);
            people[i].possible_states = NULL;
            people[i].size_possible_states = 0;
         }

         h = people[i].community;
         if(people[i].risk_history != NULL)
         {
            for(r=0; r<community[h].epi_duration; r++)
            {
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
               people[i].risk_history[r].c2p_risk = NULL;
               people[i].risk_history[r].p2p_risk = NULL;
               people[i].final_risk_day = MISSING;
            }
         }
      }
   }
   if(community != NULL)
   {
      for(h=0; h < n_community; h++)
      {
         community[h].ignore = 0;
         if(cfg_pars.preset_index == 0)
         {
            if(community[h].idx != NULL)
            {
               free(community[h].idx);
               community[h].idx = NULL;
            }      
            community[h].size_idx = 0;
            community[h].earliest_idx_day_ill = MISSING;
            community[h].latest_idx_day_ill = MISSING;
         }
         if(community[h].member_impute != NULL)
         {
            free(community[h].member_impute);
            community[h].member_impute = NULL;
            community[h].size_impute = 0;
            community[h].size_possible_states = 1;
         }
         if(community[h].sample_states != NULL) free_state_chain(community[h].sample_states, community[h].sample_states_rear);
         community[h].sample_states = community[h].sample_states_rear = NULL;
         if(community[h].list_states != NULL) free_state_chain(community[h].list_states, community[h].list_states_rear);
         community[h].list_states = community[h].list_states_rear = NULL;
         //reset risk hisotry built at risk class level 
         ptr_class = community[h].risk_class;
         while(ptr_class != NULL)
         {
            if(ptr_class->risk_history != NULL)
            {
               for(r=0; r<community[h].epi_duration; r++)
               {
                  ptr_risk = ptr_class->risk_history[r].c2p_risk;
                  while(ptr_risk != NULL)
                  {
                     ptr2_risk = ptr_risk->next;
                     if(ptr_risk->covariate != NULL)  free(ptr_risk->covariate);
                     free(ptr_risk);
                     ptr_risk = ptr2_risk;
                  }
                  ptr_risk = ptr_class->risk_history[r].p2p_risk;
                  while(ptr_risk != NULL)
                  {
                     ptr2_risk = ptr_risk->next;
                     if(ptr_risk->covariate != NULL)  free(ptr_risk->covariate);
                     free(ptr_risk);
                     ptr_risk = ptr2_risk;
                  }
                  ptr_class->risk_history[r].c2p_risk = NULL;
                  ptr_class->risk_history[r].p2p_risk = NULL;
               }
            }   
            ptr_class = ptr_class->next;
         }   
      }
   }
}


void set_index_cases(COMMUNITY *community)
{
   int i, j, k, earliest_day_ill;
      // delete old index cases and assign new index cases according to the new smaples of onset dates
   
   if(cfg_pars.preset_index == 1) // if index cases are given in pop.dat. This is true in most situations
   {
      community->size_idx = 0;
      for(j=0; j<community->size; j++)
      {
         i = community->member[j];
         if(people[i].idx == 1)  community->size_idx ++;
      }
      if(community->size_idx > 0)
      {
         community->idx = (int *) malloc((size_t) community->size_idx * sizeof(int));
         community->counter = 0;
         for(j=0; j<community->size; j++)
         {
            i = community->member[j];
            if(people[i].idx == 1)  
            {
               community->idx[community->counter] = i;
               community->counter++;
            }
         }
      }
      i = community->idx[0];
      community->earliest_idx_day_ill = community->latest_idx_day_ill = people[i].day_ill;
      for(j=1; j<community->size_idx; j++)
      {
         i = community->idx[j];
         if(people[i].day_ill < community->earliest_idx_day_ill)
            community->earliest_idx_day_ill = people[i].day_ill;
         if(people[i].day_ill > community->latest_idx_day_ill)
            community->latest_idx_day_ill = people[i].day_ill;
      }
      //printf("earliest_idx_day_ill=%d  latest_idx_day_ill=%d\n", community->earliest_idx_day_ill, community->latest_idx_day_ill);
   }
   // If user does not prespecify index cases or we are simulating data.
   // We assume index cases have the earliest day_ill. If more than one community members have the earliest day_ill,
   // all of them are index cases (an alternative wy is to randomly pick one case as index case, which is not programmed here). 
   // This assignment may not be ideal, depending on how the study enrolled communities or households.
   // However, it sould be rare to have multiple people with illness on the same day for a typical household study. 
   // Note that asymptomatic infection may also have earliest day_ill and classified as index cases.
   else
   {
      if(community->idx != NULL)
      {
         free(community->idx);
         community->idx = NULL;
      }      
      community->size_idx = 0;
      for(j=0; j<community->size; j++)
         people[community->member[j]].idx = 0;
   
      community->ignore = 0;
      for(j=0; j<community->size; j++)  people[community->member[j]].ignore = 0; 
      
      community->earliest_idx_day_ill = MISSING;
      community->latest_idx_day_ill = MISSING;
      for(j=0; j<community->size; j++)
      {
         i = community->member[j];
         if(people[i].infection == 1 && people[i].day_ill != MISSING)
         {
            if(community->earliest_idx_day_ill == MISSING) 
               community->earliest_idx_day_ill = people[i].day_ill;
            else
               community->earliest_idx_day_ill = min(community->earliest_idx_day_ill, people[i].day_ill);    
            if(community->latest_idx_day_ill == MISSING) 
               community->latest_idx_day_ill = people[i].day_ill;
            else
               community->latest_idx_day_ill = max(community->latest_idx_day_ill, people[i].day_ill);    
         }
      }
      for(j=0; j<community->size; j++)
      {
         i = community->member[j];
         if(people[i].infection == 1 && people[i].day_ill != MISSING)
         {
            //if(fabs(people[i].day_ill - earliest_day_ill) < cfg_pars.min_incubation)
            //if(people[i].day_ill == community->earliest_idx_day_ill && community->size_idx == 0) //pick only one person as index case
            if(people[i].day_ill == community->earliest_idx_day_ill) //all cases with illness onset on day earliest_idx_day_ill are aset as index case
            { 
               people[i].idx = 1;
               community->size_idx ++;
            }   
         }
      }
      // Families wihtout any infection are excluded from further analysis
      if(community->size_idx > 0)
      {
         community->idx = (int *) malloc((size_t) community->size_idx * sizeof(int));
         community->counter = 0;
         for(j=0; j<community->size; j++)
         {
            i = community->member[j];
            if(people[i].idx == 1)  
            {
               //printf("person %d: day ill=%d\n", i, people[i].day_ill);
               community->idx[community->counter] = i;
               community->counter++;
            }   
         }
      }
   }
}
