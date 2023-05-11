////// Author: Yang Yang
////// Department of Biostatistics and Emerging Pathogens Institute
////// University of Florida

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>


#define logit(x) ( log( (x) / ( 1.0 - (x) ) ) )
#define inv_logit(x) ( 1.0 / ( 1.0 + exp(-(x)) ) )

#define bipow(x,y) (((y)==0)? 1:(x))
#define max(x, y) ((x)>(y)?(x):(y))
#define min(x, y) ((x)<(y)?(x):(y))
#define logit(x) ( log( (x) / ( 1.0 - (x) ) ) )
#define inv_logit(x) ( 1.0 / ( 1.0 + exp(-(x)) ) )

#define MISSING 99999
#define INFINITY_INTEGER 100000
#define close_to_0 1e-20
#define close_to_1 0.99999
#define close_to_05 0.49999
#define num_ini_estimation 50
#define num_ini_test 10

#include "matrix.h"
#include "mathfunc.h"
#include "distribution.h"
#include "statistics.h"
#include "optimization.h"

#include "datastruct.h" /*specifying data structure*/
CFG_PARS cfg_pars;
PEOPLE *people; /* pointer to people */
int p_size;

COMMUNITY *community;
int n_community;

int DEBUG = 0;

int max_epi_duration;
long seed;
int serial_number;
double *pdf_incubation, *sdf_incubation;

int size_sample_states;
double *importance_weight;

int ITER;
int n_need_MCEM, n_need_OEM;

long TIME_SAMPLING_SEC, TIME_SAMPLING_CPU;
long TIME_ITERATION_SEC, TIME_ITERATION_CPU;
long TIME_VARIANCE_SEC, TIME_VARIANCE_CPU;

int get_size(FILE *in)
{
  int c, i=0;
  rewind(in);
  while((c=getc(in))!=EOF)
    if(c == '\n') i++;
  return i;
 }

#include "config.h" /*read in configuration parameters for flu and trial*/
#include "utility.h"
#include "simulate.h"
#include "preamble.h"
#include "derivatives.h"
#include "estimation.h"
#include "test_cs_hh_improve.h"
#include "core.h"


int main(int argc, char *argv[])
{
  int h, i, j, k, l, m, n, t, len, len_inc, len_inf;
  int *day_epi_start, *day_epi_stop;
  int R0_window_bound;
  FILE * file;
  char command[100], file_name[100];

  serial_number = 1;
  seed = 12345678; 
  if(argc > 1)
     serial_number = atoi(argv[1]);
  if(argc > 2)
     seed = atoi(argv[2]);

  mt_init(seed);
  
  printf("serial_number=%d\n", serial_number);
  if((file = fopen("config.file", "r")) == NULL)
  {
    printf("couldn't open config.file\n");
    exit(0);
  }
  cfg_pars = get_cfg_pars(file);
  fclose(file);
  
  sprintf(file_name, "%soutput_%d.txt", cfg_pars.path_out, serial_number);
  //if( fopen(file_name, "r") != NULL )  
  //{
  //   fclose(file); 
  //   if(cfg_pars.windows_os == 1)  sprintf(command, "del -f %soutput_%d.txt", cfg_pars.path_out, serial_number);
  //   else  sprintf(command, "rm -f %soutput_%d.txt", cfg_pars.path_out, serial_number);
  //   system(command);
  //}
  for(i=0; i<cfg_pars.n_inc; i++)
  {
     cfg_pars.max_incubation = cfg_pars.max_inc[i];
     cfg_pars.min_incubation =cfg_pars.min_inc[i];
     len_inc = cfg_pars.max_inc[i] - cfg_pars.min_inc[i] + 1;
     cfg_pars.prob_incubation = (double *)malloc((size_t) (len_inc * sizeof(double)));
     for(k = 0 ; k < len_inc ; k++)
        cfg_pars.prob_incubation[k] = cfg_pars.prob_inc[i][k];
     for(j=0; j<cfg_pars.n_inf; j++)
     {
        cfg_pars.upper_infectious = cfg_pars.upper_inf[j];
        cfg_pars.lower_infectious = cfg_pars.lower_inf[j];
        if(cfg_pars.lower_infectious < 0 && cfg_pars.common_contact_history_within_community == 1)
        {
           printf("Error: infectious period should not start before symptom onset when common contact history within community is assumed.\n");
           printf("This is because, the assumption of common contact history will lead to the construction of a common risk history\n");
           printf("for the whole community based on the infectious periods of all cases.\n");
           printf("The infectiousness before symptom onset of an infected person will then contribute to his own likelihood of infection,\n");
           printf("which is wrong as a person cannot be infected by himself.\n");
           free(cfg_pars.prob_incubation);
           goto main_end;
        }   
           
        printf("min_inc=%d  max_inc=%d  lower_inf=%d  upper_inf=%d\n", cfg_pars.min_incubation, cfg_pars.max_incubation, cfg_pars.lower_infectious, cfg_pars.upper_infectious);
        len_inf = cfg_pars.upper_inf[j] - cfg_pars.lower_inf[j] + 1;
        cfg_pars.prob_infectious = (double *)malloc((size_t) (len_inf * sizeof(double)));
        for(l = 0 ; l < len_inf ; l++)
           cfg_pars.prob_infectious[l] = cfg_pars.prob_inf[j][l];
        if(cfg_pars.effective_bounds_provided == 0)
        {
           for(l=0; l<cfg_pars.n_p_mode; l++)
           {
              cfg_pars.effective_lower_infectious[l] = cfg_pars.lower_infectious; 
              cfg_pars.effective_upper_infectious[l] = cfg_pars.upper_infectious; 
           }           
        }
        if(cfg_pars.simplify_output == 0)
        {
           sprintf(file_name, "%soutput.txt", cfg_pars.path_out);
           if((file = fopen(file_name, "a")) == NULL)
              file = fopen(file_name, "w");
           fprintf(file, "\n\n Serial number=%d\n", serial_number);
           fprintf(file, "=======================================================================================\n");
           fprintf(file, "%d: min_inc=%d  max_inc=%d  prob_inc=(", i+1, cfg_pars.min_incubation, cfg_pars.max_incubation);
           for(k = 0 ; k < len_inc ; k++)  fprintf(file, "%5.3lf, ", cfg_pars.prob_incubation[k]);
           fprintf(file, ")\n");
           fprintf(file, "\n%d: lower_inf=%d  upper_inf=%d  prob_inf=(", j+1, cfg_pars.lower_infectious, cfg_pars.upper_infectious);
           for(l = 0 ; l < len_inf ; l++)  fprintf(file, "%5.3lf, ", cfg_pars.prob_infectious[l]);
           fprintf(file, ")\n");
           fprintf(file, "=======================================================================================\n");
           fclose(file);
        }  
        // to get R0 estimates for a series of cut-point days. We assume p is constant before and after the cut-point but with different values.
        // 
        printf("R0_divide_by_time=%d\n", cfg_pars.R0_divide_by_time);
        if(cfg_pars.R0_divide_by_time == 1) 
        {
           if(!(cfg_pars.n_p_mode == 2 | cfg_pars.n_p_mode == 3))
           {
              printf("For time-varying reproductive numbers, the number of p2p transmission probabilities should be either 2 or 3\n");
              goto main_end;
           }   
           sprintf(file_name, "%scommunity.dat", cfg_pars.path_in);
           if((file = fopen(file_name, "r")) == NULL)
           {
              printf("No valid community file\n");
              goto main_end;
           }   
           n_community = get_size(file);
           day_epi_start = (int *)malloc((size_t) n_community * sizeof(int));
           day_epi_stop = (int *)malloc((size_t) n_community * sizeof(int));
           rewind(file); 
           for(h=0; h < n_community; h++)
           {
              fscanf(file, "%d  %d  %d", &m, &day_epi_start[h], &day_epi_stop[h]);
              if(cfg_pars.R0_divide_start <= day_epi_start[h] || cfg_pars.R0_divide_stop >= day_epi_stop[h])
              {
                  printf("Dividion days are out of epidemic day range of community %d\n", h);
                  free(day_epi_start); free(day_epi_stop);
                  fclose(file);
                  goto  main_end;
              }
           }      
           fclose(file);
              
           for(t=cfg_pars.R0_divide_start; t<=cfg_pars.R0_divide_stop; t++)
           {
              sprintf(file_name, "%sp2p_contact.dat", cfg_pars.path_in);
              file = fopen(file_name, "w");
              for(h=0; h < n_community; h++)
              { 
                 if(cfg_pars.n_p_mode == 2)
                 {
                    fprintf(file, "%8d  %8d  %8d  %2d  %6.2f  %2d\n", h, day_epi_start[h], t, 0, 0.0, 0);
                    fprintf(file, "%8d  %8d  %8d  %2d  %6.2f  %2d\n", h, t+1, day_epi_stop[h], 1, 0.0, 0);
                 }
                 else
                 {   
                    if(cfg_pars.R0_window_size == 0)
                    {
                       printf("Window size should be an integer greater than 0\n");
                       goto main_end;
                    }   
                    R0_window_bound = t - cfg_pars.R0_window_size + 1;
                 
                    if(R0_window_bound <= day_epi_start[h])  continue;
                    fprintf(file, "%8d  %8d  %8d  %2d  %6.2f  %2d\n", h, day_epi_start[h], R0_window_bound-1, 0, 0.0, 0);
                    fprintf(file, "%8d  %8d  %8d  %2d  %6.2f  %2d\n", h, R0_window_bound, t, 1, 0.0, 0);
                    fprintf(file, "%8d  %8d  %8d  %2d  %6.2f  %2d\n", h, t+1, day_epi_stop[h], 2, 0.0, 0);
                 }   
              }
              fclose(file);
              core(i, j, t);
           }
           free(day_epi_start); free(day_epi_stop);       
        }
        else    core(i, j, 0);
        free(cfg_pars.prob_infectious);
     }  
     free(cfg_pars.prob_incubation);
  }  

main_end:   
  free(cfg_pars.min_inc);
  free(cfg_pars.max_inc);
  free(cfg_pars.lower_inf);
  free(cfg_pars.upper_inf);
  for(i=0; i<cfg_pars.n_inc; i++)  free(cfg_pars.prob_inc[i]);
  free(cfg_pars.prob_inc);
  for(i=0; i<cfg_pars.n_inf; i++)  free(cfg_pars.prob_inf[i]);
  free(cfg_pars.prob_inf);

  free(cfg_pars.c2p_covariate);
  free(cfg_pars.sus_p2p_covariate);
  free(cfg_pars.inf_p2p_covariate);
  free(cfg_pars.pat_covariate);
  free(cfg_pars.imm_covariate);
  free_2d_array_int(cfg_pars.interaction);
  free_2d_array_double(cfg_pars.ini_par_effective);
  free(cfg_pars.lower_search_bound);
  free(cfg_pars.upper_search_bound);
  free(cfg_pars.converge_criteria);

  free_2d_array_double(cfg_pars.SAR_sus_time_ind_covariate);
  free_2d_array_double(cfg_pars.SAR_inf_time_ind_covariate);
  free_3d_array_double(cfg_pars.SAR_sus_time_dep_covariate);
  free_3d_array_double(cfg_pars.SAR_inf_time_dep_covariate);
  free(cfg_pars.R0_multiplier);
  free(cfg_pars.R0_multiplier_var);
  free(cfg_pars.par_fixed_id);
  free(cfg_pars.par_fixed_value);
  free(cfg_pars.sim_par_effective);
  free(cfg_pars.effective_lower_infectious);
  free(cfg_pars.effective_upper_infectious);
  if(cfg_pars.par_equiclass != NULL)
  {
     for(i=0; i<cfg_pars.n_par_equiclass; i++)
     {
        free(cfg_pars.par_equiclass[i].member);
     }
     free(cfg_pars.par_equiclass);
  }
  return(0);
}
