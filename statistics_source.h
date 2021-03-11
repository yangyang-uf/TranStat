/*================================
  Basic statistic summary functions
  ================================*/

double function(sum,BASE)(DATATYPE *array, unsigned int size)
{
  unsigned int i;
  double sum=0;
  for(i =0; i<size; i++)
    sum += (double) array[i];
  return sum;
}

double function(sumsquare,BASE)(DATATYPE *array, unsigned int size)
{
  unsigned int i;
  double sum=0;
  for(i =0; i<size; i++)
    sum += ((double) array[i]) * ((double) array[i]);
  return sum;
}

double function(cross_product,BASE)(DATATYPE *array1, DATATYPE *array2, unsigned int size)
{
  unsigned int i;
  double sum=0;
  for(i =0; i<size; i++)
    sum += ((double) array1[i]) * ((double) array2[i]);
  return sum;
}

double function(mean,BASE)(DATATYPE *array, unsigned int size)
{
  unsigned int i;
  double sum=0;
  for(i =0; i<size; i++)
    {
      sum += array[i];
    }  
  return sum/(double)size;
}

double function(variance,BASE)(DATATYPE *array, unsigned int size)
{
  unsigned int i;
  double sum, sum2, v;
  sum = sum2 = 0;
  for(i =0; i<size; i++)
  {
    sum += array[i];
    sum2 += array[i]*array[i];
  }
  v = (sum2-sum*sum/size)/(size-1);
  if(v < 0)  v = 0;
  return v;
}

double function(mse,BASE)(DATATYPE *array, double mu, unsigned int size)
{
  unsigned int i;
  double sum;
  sum = 0;
  for(i =0; i<size; i++)
  {
    sum += (array[i] - mu) * (array[i] - mu);
  }

  return sum/size;
}

double function(covariance,BASE)(DATATYPE *array1, DATATYPE *array2, unsigned int size)
{
  unsigned int i;
  double sum_array1, sum_array2, sum_product;
  sum_array1 = sum_array2 = sum_product = 0;
  for(i =0; i<size; i++)
  {
    sum_array1 += array1[i];
    sum_array2 += array2[i];
    sum_product += array1[i]*array2[i];
  }

  return (sum_product-sum_array1*sum_array2/size)/(size-1);
}


double function(sd,BASE)(DATATYPE *array, unsigned int size)
{
  unsigned int i;
  double sum, sum2, variance;
  sum = sum2 = 0;
  for(i =0; i<size; i++)
  {
    sum += array[i];
    sum2 += array[i]*array[i];
  }
    
  variance = (sum2-sum*sum/size)/(size-1);
  return sqrt(variance);
}
       
double function(skewness,BASE)(DATATYPE *array, unsigned int size)
{
  unsigned int i;
  double sum1, sum2, sum3;
  double std, sum;
  sum1 = sum2 = sum3 = 0;
  for(i =0; i<size; i++)
  {
    sum1 += array[i];
    sum2 += array[i]*array[i];
    sum3 += array[i]*array[i]*array[i];
  }
  std = sqrt((sum2-sum1*sum1/size)/(size-1));
  sum = sum3-3*sum1*sum2/size+2*sum1*sum1*sum1/(size*size);
  return ((sum*size)/(std*std*std*(size-1)*(size-2)));
}

DATATYPE function(max, BASE) (DATATYPE *array, unsigned int size)
{
  /* finds the largest member of a dataset */

  DATATYPE max = array[0];
  size_t i;
  for (i = 0; i < size; i++)
    {
      if (array[i] > max)
        max = array[i];
    }

  return max;
}

DATATYPE function(min, BASE) (DATATYPE *array, unsigned int size)
{
  /* finds the largest member of a dataset */

  DATATYPE min = array[0];
  size_t i;
  for (i = 0; i < size; i++)
    {
      if (array[i] < min)
        min = array[i];
    }

  return min;
}

DATATYPE function(quantile, BASE) (DATATYPE *array, unsigned int size, double percent)
{
  int i, *order;
  double *sorted_data;
  const double index = percent * (size - 1);
  const int lhs = (int) index;
  const double delta = index - lhs;
  double result;

  if (size == 0)
    return 0.0 ;

  sorted_data = (double *) malloc((size_t) size * sizeof(double));
  order = (int *) malloc((size_t) size * sizeof(int));
  for(i=0; i<size; i++)  
  {
     sorted_data[i] = (double) array[i];
     order[i] = i;
  }   
  quick_sort(1, size, sorted_data, order);
  
  if (lhs == size - 1)
  {
     result = sorted_data[lhs];
  }
  else 
  {
     result = (1.0 - delta) * sorted_data[lhs] + delta * sorted_data[lhs + 1];
  }

  free(sorted_data);
  free(order);
  return result ;
}
