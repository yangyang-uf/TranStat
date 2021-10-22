#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int get_size(FILE *in)
{
  int c, i=0;
  rewind(in);
  while((c=getc(in))!=EOF)
    if(c == '\n') i++;
  return i;
 }

int main(int argc, char *argv[])
{
  FILE *file;
  int h, i, j, k, l, n, m, r, t;
  int epi_duration;
  int *n_case, *n_inf, *day_ill;
  int n_phase = 2;
  int shift = 50;
  int pop_multiplier = 50;
  int n_total;
  int n_total_inf, n_total_esc;
  int start_day, stop_day; 
  int infectious_start_day, infectious_stop_day; 
  int susceptible_start_day, susceptible_stop_day; 
  int seed, idx, idx_cut_day;
  int *phase_start_day, *phase_stop_day;
  int value_int;
  double value_double, p;

  idx_cut_day = 10;

  if((file = fopen("salvado.dat", "r")) == NULL)
  {
     printf("No valid data\n");
     exit(0);
  }   
  epi_duration = get_size(file);
  rewind(file);
  n_case = (int *) malloc((size_t) epi_duration * sizeof(int));
  for(t=0; t<epi_duration; t++)
     fscanf(file, "%d", &n_case[t]);
  fclose(file);
  printf("Epidemic duration=%d\n", epi_duration);
 
  phase_start_day = (int *) malloc((size_t) n_phase * sizeof(int));
  phase_stop_day = (int *) malloc((size_t) n_phase * sizeof(int));
  // for giaradot
  //phase_start_day[0] = 0; phase_stop_day[0] = 40;
  //phase_start_day[1] = 41; phase_stop_day[1] = epi_duration - 1;
  //for salvado
  phase_start_day[0] = 0; phase_stop_day[0] = 100;
  phase_start_day[1] = 101; phase_stop_day[1] = epi_duration - 1;

  n_inf = (int *) malloc((size_t) epi_duration * sizeof(int));
  
  n_total_inf = 0;
  for(t=0; t<epi_duration; t++)
  {
     if(t >=0 && t <= 2)  p = 1.0;
     if(t >=3 && t <= 5)  p = 1.0;
     if(t >=6 && t <= 8)  p = 1.0;
     if(t >=9 && t <= 11)  p = 1.0;
     if(t >=12)  p = 1.0;
     value_double = n_case[t] / p;
     value_int = (int) floor(n_case[t] / p);
     n_inf[t] = (value_double - value_int < 0.5)? value_int : (value_int + 1);
     n_total_inf += n_inf[t];
     printf("%d: %d\n", t, n_inf[t]);
  }
  n_total = n_total_inf * pop_multiplier;
  printf("n_total_inf=%d  n_total=%d\n", n_total_inf, n_total);

  n_total_esc = n_total - n_total_inf;
  i = 0;
  file = fopen("pop_1.dat", "w");
  {
     for(t=0; t<epi_duration; t++)
     if(n_inf[t] > 0)
     {
        idx = 0;
        if(t <= idx_cut_day)  idx = 1;
        fprintf(file, "%6d   %5d   %1d   %1d   %1d   %5d   %1d   %5d   %1d   %3d   %3d   %12d  %1d\n",
                   i, 0, 0, 1, 1, shift + t, 0, -1, idx, 0, 0, n_inf[t], 0);  
        i++;        
     }
     if(n_total_esc > 0)
        fprintf(file, "%6d   %5d   %1d   %1d   %1d   %5d   %1d   %5d   %1d   %3d   %3d   %12d  %1d\n",
                i, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, n_total_esc, 0);  
     i++;        
  }
  fclose(file);

  
  file = fopen("community_1.dat", "w");
  fprintf(file, "%d   %d   %d\n", 0, 1, shift + epi_duration - 1);
  fclose(file);


  file = fopen("p2p_contact_1.dat", "w");
  for(j=0; j<n_phase; j++)
  { 
     start_day = shift + phase_start_day[j];
     stop_day = shift + phase_stop_day[j];
     fprintf(file, "0  %8d  %8d  %2d  %3.1f  %2d\n", start_day, stop_day, j, 0.0, 0);
   }
   fclose(file);

  free(phase_start_day);
  free(phase_stop_day);
  free(n_case);
  free(n_inf);
}
  
  
