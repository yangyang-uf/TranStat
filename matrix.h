 typedef struct{
  unsigned int row;
  unsigned int col;
  double ** data;
  } MATRIX;

void initialize_matrix(MATRIX *mat)
{
   mat->row = mat->col = 0;
   mat->data = NULL;
}

void deflate_matrix(MATRIX *mat)
{
  if( mat->row > 0 && mat->col > 0  && mat->data != NULL)
  {
     free(mat->data[0]);
     free(mat->data);
  }
  mat->row=0; 
  mat->col=0;
  mat->data = NULL;
}

void inflate_matrix(MATRIX *mat, unsigned int nrow, unsigned int ncol, double ini)
{
  unsigned int i, j;
  if(nrow == 0 || ncol == 0)
  {
    printf("\n colume and row size need to be positive integers!\n"); 
    exit(0);
  }
  if( mat->row > 0 && mat->col > 0  && mat->data != NULL)
  {
     free(mat->data[0]);
     free(mat->data);
  }
  mat->row = nrow;
  mat->col = ncol;
  mat->data = (double **)malloc((size_t) (nrow*sizeof(double *))); 
  mat->data[0] = (double *)malloc((size_t) ((nrow*ncol)*sizeof(double))); 
  for(i=1;i<nrow;i++) 
    mat->data[i] = mat->data[i-1] + ncol; 
  for(i=0;i<nrow;i++) 
    for(j=0;j<ncol;j++) 
      mat->data[i][j] = ini;
}
      
# define swap(a,b) {temp=(a); (a)=(b); (b)=temp;}
int inverse(MATRIX *mat1, MATRIX *mat2)
{
  double **a;
  int *indxc, *indxr, *ipiv;
  int i, icol, irow, j, k, l, ll, n;
  double big, dum, pivinv, temp;
  
  if(mat1->row != mat1->col || mat2->row != mat2->col)
  {
    printf("\n Unequal row and column size!\n");
    exit(0);
  }
    
  n = mat1->row;
  
  if(n == 1)
  {
     mat2->data[0][0] = 1.0 / mat1->data[0][0];
     return(0);
  } 

  if(n == 2)
  {
     temp = mat1->data[0][0] * mat1->data[1][1] - mat1->data[0][1] * mat1->data[1][0];
     if(fabs(temp) < 1e-40)
     {
        printf("Singular 2X2 matrix\n");  
        return(1);
     }
     mat2->data[0][0] = mat1->data[1][1] / temp;
     mat2->data[0][1] = -mat1->data[0][1] / temp;
     mat2->data[1][0] = -mat1->data[1][0] / temp;
     mat2->data[1][1] = mat1->data[0][0] / temp;
     return(0);
  }
 
    
  a = (double **)malloc((size_t) (n*sizeof(double *))); 
  a[0] = (double *)malloc((size_t) ((n*n)*sizeof(double))); 
  for(i=1;i<n;i++) 
    a[i] = a[i-1] + n; 
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      a[i][j] = mat1->data[i][j];
  indxc = (int *)malloc((size_t) ((n)*sizeof(int)));
  indxr = (int *)malloc((size_t) ((n)*sizeof(int)));
  ipiv = (int *)malloc((size_t) ((n)*sizeof(int)));
  for(j=0;j<n;j++) ipiv[j] = 0;
  for(i=0;i<n;i++)
  {
    big=0.0;
    for(j=0;j<n;j++)
      if(ipiv[j]!=1)
        for(k=0;k<n;k++)
        {
          if(ipiv[k] == 0)
          {
            if(fabs(a[j][k]) >= big)
            {
              big = fabs(a[j][k]);
              irow=j;
              icol=k;
            }
          }
          else if(ipiv[k]>1) 
          {
            printf("\n Singular Matrix\n");
	    free(a[0]);
	    free(a);
	    free(indxc);
	    free(indxr);
	    free(ipiv);
            return(1);
          }
        }
    ++(ipiv[icol]);
    if(irow != icol)
    {
      for(l=0;l<n;l++)  swap(a[irow][l], a[icol][l])
    }
    indxr[i] = irow;
    indxc[i] = icol;
    if(a[icol][icol] == 0.0)   
    {
        printf("\n Singular Matrix\n");
	free(a[0]);
	free(a);
	free(indxc);
	free(indxr);
	free(ipiv);
        return(1);
    }
    pivinv = 1.0/a[icol][icol];
    a[icol][icol] = 1.0;
    for(l=0;l<n;l++)  a[icol][l] *= pivinv;
    for(ll=0;ll<n;ll++)
      if(ll != icol)
      {
        dum=a[ll][icol];
        a[ll][icol]=0.0;
        for(l=0;l<n;l++)  a[ll][l] -= a[icol][l]*dum;
      }
  }
  for(l=n-1;l>=0;l--)
  {
    if(indxr[l] != indxc[l])
      for(k=0;k<n;k++)
        swap(a[k][indxr[l]], a[k][indxc[l]])
  }
  
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      mat2->data[i][j] = a[i][j];
  free(a[0]);
  free(a);
  free(ipiv);
  free(indxr);
  free(indxc);
  return(0);
}              

/*** Inversion using Cholesky decomposition. The input matrix must be
     a symmetric positive matrix ***/
int cholesky_inverse(MATRIX *mat1, MATRIX *mat2)
{
  double **a, **b, *p;
  int i, j, k, n, error;
  double sum;
  
  if(mat1->row != mat1->col || mat2->row != mat2->col)
  {
    printf("\n Unequal row and column size!\n");
    return(2);
  }
    
  n = mat1->row;
  
  if(n == 1)
  {
     mat2->data[0][0] = 1.0 / mat1->data[0][0];
     return(0);
  }
     
  a = (double **)malloc((size_t) (n*sizeof(double *))); 
  a[0] = (double *)malloc((size_t) ((n*n)*sizeof(double))); 
  for(i=1;i<n;i++) 
    a[i] = a[i-1] + n; 
  for(i=0;i<n;i++)
     for(j=i;j<n;j++)
     {
        a[i][j] = mat1->data[i][j];
        if(i != j) a[j][i] = 0;
     }   

  b = (double **)malloc((size_t) (n*sizeof(double *))); 
  b[0] = (double *)malloc((size_t) ((n*n)*sizeof(double))); 
  for(i=1;i<n;i++) 
    b[i] = b[i-1] + n; 
      
  p = (double *)malloc((size_t) (n*sizeof(double))); 
  
  for(i=0; i<n; i++)
  {
     for(j=i; j<n; j++)
     {
        for(sum=a[i][j],k=i-1; k>=0; k--)
           sum -= a[i][k] * a[j][k];
        if(i == j)
        {
           if(sum <= 0.0)
           {
              printf("cholesky decomposition failed\n");
              error = 1;
              goto end;
           }
           p[i] = sqrt(sum);
        }
        else
        {
           a[j][i] = sum / p[i];
        }
     }
  }
  
  
  for(i=0; i<n; i++)
  {
     a[i][i] = 1.0 / p[i];
     for(j=i+1; j<n; j++)
     {
        sum = 0.0;
        for(k=i; k<j; k++)
        {
           sum -= a[j][k] * a[k][i];
        }
        a[j][i] = sum / p[j];
     }
  }
  
  
  /*** a now is the lower triangle L^(-1) ***/
  
  for(i=0;i<n;i++)
     for(j=i;j<n;j++)
     {
        b[i][j] = a[j][i];
        if(i != j)
        {
           b[j][i] = 0;
           a[i][j] = 0;
        }
     }      
  
  for(i=0;i<n;i++)
     for(j=i;j<n;j++)
     {
        mat2->data[i][j] = 0;
        for(k=j; k<n; k++)
           mat2->data[i][j] += b[i][k] * a[k][j];
           
        mat2->data[j][i] = mat2->data[i][j];   
     }      

  error = 0;

end:
  free(a[0]);
  free(a);
  free(b[0]);
  free(b);
  free(p);
  return(error);
}  
                 

void product(MATRIX *mat1, MATRIX *mat2, MATRIX *mat3)
{
  unsigned int i, j, k;
  double sum;
  MATRIX temp;


  if(mat1->col != mat2->row || mat3->row != mat1->row || mat3->col != mat2->col) 
  {
     printf("\n colume and row size do not match for multiplication!\n");
     exit(0);
  }

  initialize_matrix(&temp);
  inflate_matrix(&temp, mat1->row, mat2->col, 0.0);

  for(i=0;i<mat1->row;i++)
     for(j=0;j<mat2->col;j++)
     {
        sum = 0;
        for(k=0;k<mat1->col;k++)
           sum += mat1->data[i][k]*mat2->data[k][j];
        temp.data[i][j] = sum;
     }

  for(i=0;i<mat1->row;i++)
     for(j=0;j<mat2->col;j++)
     {
        mat3->data[i][j] = temp.data[i][j];
     }
  deflate_matrix(&temp);
}

void triple_product(MATRIX *mat1, MATRIX *mat2, MATRIX *mat3, MATRIX *mat4)
{
  unsigned int i, j, k;
  double sum;
  MATRIX temp1, temp2;
  
  if(mat1->col != mat2->row || mat2->col != mat3->row || mat4->row != mat1->row || mat4->col != mat3->col) 
  {
     printf("\n colume and row size do not match for multiplication!\n");
     exit(0);
  }

  initialize_matrix(&temp1);
  inflate_matrix(&temp1, mat1->row, mat2->col, 0.0);
  
  for(i=0;i<mat1->row;i++)
     for(j=0;j<mat2->col;j++)
     {
        sum = 0;
        for(k=0;k<mat1->col;k++)
           sum += mat1->data[i][k] * mat2->data[k][j];
        temp1.data[i][j] = sum;
     }
     
  initialize_matrix(&temp2);
  inflate_matrix(&temp2, mat1->row, mat3->col, 0.0);

  for(i=0;i<temp1.row;i++)
     for(j=0;j<mat3->col;j++)
     {
        sum = 0;
        for(k=0;k<temp1.col;k++)
           sum += temp1.data[i][k] * mat3->data[k][j];
        temp2.data[i][j] = sum;
     }

  for(i=0;i<temp2.row;i++)
     for(j=0;j<temp2.col;j++)
     {
        mat4->data[i][j] = temp2.data[i][j];
     }

  deflate_matrix(&temp1);   
  deflate_matrix(&temp2);   

}
      
      
      
double sandwich(double *x, MATRIX *A, double *y)
{
  unsigned int i, j, n;
  double sum;
  
  if(A->col != A->row) 
  {
     printf("\n colume and row size do not match for sandwich()!\n");
     exit(0);
  }

  n = A->row;
  sum = 0.0;
  for(i=0;i<n;i++)
  {
     for(j=0;j<n;j++)
     {
        sum += x[i] * A->data[i][j] * y[j];
     }
  }   
  return(sum);   
}

double trace(MATRIX *A)
{
  unsigned int i, j, n;
  double sum;

  if(A->col != A->row)
  {
     printf("\n colume and row size do not equal for trace()!\n");
     exit(0);
  }

  n = A->row;
  sum = 0.0;
  for(i=0;i<n;i++)
  {
     sum += A->data[i][i];
  }
  return(sum);
}

/*this is only for self-product of symmetric matrices*/      
void selfproduct(MATRIX *mat1, MATRIX *mat2)
{
  unsigned int i, j, k;
  double sum;
  MATRIX mat;

  initialize_matrix(&mat);
  inflate_matrix(&mat, mat1->row, mat1->col, 0.0);
  for(i=0;i<mat1->row;i++)
     for(j=i;j<mat1->row;j++)
     {
        sum = 0;
        for(k=0;k<mat1->col;k++)
           sum += mat1->data[i][k]*mat1->data[j][k];
        mat.data[i][j] = sum;
        mat.data[j][i] = sum;  
     }
  for(i=0;i<mat.row;i++)
     for(j=i;j<mat.col;j++)
     {
        mat2->data[i][j] = mat.data[i][j];
        mat2->data[j][i] = mat.data[i][j];
     }
  deflate_matrix(&mat);
}

void addition(MATRIX *mat1, double coeff1, MATRIX *mat2, double coeff2,  MATRIX *mat3)
{
   unsigned int i, j;

   if(mat1->row != mat2->row || mat3->row != mat1->row ||
      mat1->col != mat2->col || mat3->col != mat1->col)
   {
      printf("\n colume and row size do not match for addition!\n");
      exit(0);
   }
   for(i=0;i<mat1->row;i++)
      for(j=0;j<mat1->col;j++)
      {
          mat3->data[i][j] = coeff1 * mat1->data[i][j] + coeff2 * mat2->data[i][j];
      }

}

void scalar_addition(MATRIX *mat1, double coeff1, double coeff2,  MATRIX *mat2)
{
   unsigned int i, j;

   if(mat2->row != mat1->row || mat2->col != mat1->col)
   {
      printf("\n colume and row size do not match for scalor addition!\n");
      exit(0);
   }
   for(i=0;i<mat1->row;i++)
      for(j=0;j<mat1->col;j++)
      {
          mat2->data[i][j] = coeff1 + coeff2 *  mat1->data[i][j];
      }
}



void substract(MATRIX *mat1, MATRIX *mat2, MATRIX *mat3)
{
   unsigned int i, j;

   if(mat1->row != mat2->row || mat3->row != mat1->row ||
      mat1->col != mat2->col || mat3->col != mat1->col)
   {
      printf("\n colume and row size do not match for substrcttion!\n");
      exit(0);
   }
   for(i=0;i<mat1->row;i++)
      for(j=0;j<mat1->col;j++)
      {
          mat3->data[i][j] = mat1->data[i][j] - mat2->data[i][j];
      }

}

int transpose(MATRIX *mat1, MATRIX *mat2)
{
   unsigned int i, j;

   if(mat1->row != mat2->col || mat1->col != mat2->row)
   { 
      printf("\n colume and row size do not match for transpose!\n"); 
      return(1);
   }

   for(i=0;i<mat1->row;i++)
      for(j=0;j<mat1->col;j++)
      {
         mat2->data[j][i] = mat1->data[i][j];
      }
   return(0);
}

int duplicate(MATRIX *mat1, MATRIX *mat2)
{
   unsigned int i, j;

   if(mat1->row != mat2->row || mat1->col != mat2->col)
   { 
      printf("\n colume and row size do not match for duplication!\n"); 
      return(1);
   }

   for(i=0;i<mat1->row;i++)
      for(j=0;j<mat1->col;j++)
      {
         mat2->data[i][j] = mat1->data[i][j];
      }
   return(0);
}

void diagonal(MATRIX * array, MATRIX * mat)
{
  unsigned int i, j;
  if(array->col != 1) {printf("\n must be colume vector!\n"); exit(0);}
  if(array->row == 1) {printf("\n vector must have a length >= 2!\n"); exit(0);}
  for(i=0;i<mat->row;i++)
    for(j=0;j<mat->col;j++)
    {
        mat->data[i][j] = (i==j)?array->data[i][0]:0;
    }

}

void reset(MATRIX * mat, double ini)
{
  unsigned int i, j;
  for(i=0;i<mat->row;i++)
    for(j=0;j<mat->col;j++)
    {
        mat->data[i][j] = ini;
    }

}

/*gives lower triangle of cholesky decomposition*/
int cholesky_decomp(MATRIX *A, MATRIX *B)
{
  int n = A->row;
  int m = A->col;
  int i,j,k;
  int error = 0;
  double *diag, sum;
  
  if (n != m)
  {
     printf("cholesky decomposition requires square matrix\n");
     exit(1);
  }
  else
  {
     if(n == 1) 
     {
        B->data[0][0] = sqrt(A->data[0][0]);
        return(0);
     }  
      
     duplicate(A, B);
     diag = (double *) malloc((size_t) n * sizeof(double));
     for(i=0; i < n; i++)
     {
        for(j=i; j<n; j++)
        {
           sum = B->data[i][j];
           for(k=i-1; k>= 0; k--)
              sum -= B->data[i][k] * B->data[j][k];
           if(i == j) 
           {
              if(sum <= 0.0)
              {
                 error = 1;
                 goto end;
              }
              diag[i] = sqrt(sum);
           }
           else
              B->data[j][i] = sum / diag[i];
        }
     }
     for(i=0; i < n; i++)
     {
        for(j=i; j<n; j++)
        {
           if(i == j)
              B->data[i][j] = diag[i];
           else   
              B->data[i][j] = 0.0;
        }
     }                  
  }

  end:
     free(diag);
     return(error);
}

double det_sympos(MATRIX *A)
{
   int i, n;
   double det;
   MATRIX B;
   
   n = A->row;
   
   if(n == 1)
   {
      return(A->data[0][0]);
   }

   if(n == 2)
   {
      return(A->data[0][0] * A->data[1][1] - A->data[0][1] * A->data[1][0]);
   }

   initialize_matrix(&B);
   inflate_matrix(&B, n, n, 0.0);
   
   if(cholesky_decomp(A, &B) == 1)
      det = 0.0;
   else
   {
      det = 1.0;
      for(i=0; i<n; i++)
         det *= B.data[i][i] * B.data[i][i];
   }         
   deflate_matrix(&B);
   return(det);
}
   
void fmprintf(MATRIX *mat)
{
  unsigned int i, j;
  printf("\n");
  for(i=0;i<mat->row;i++)
  {
    for(j=0;j<mat->col;j++)
    {
      printf("%e  ",mat->data[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

void fmfprintf(FILE *file, MATRIX *mat, int by_row)
{
  unsigned int i, j;
  if(by_row == 1)
  {
     for(i=0;i<mat->row;i++)
     {
        for(j=0;j<mat->col;j++)
        {
           fprintf(file, "%e  ",mat->data[i][j]);
        }
        fprintf(file, "\n");
     }
     fprintf(file, "\n");
  }
  else
  {
     for(i=0;i<mat->col;i++)
     {
        for(j=0;j<mat->row;j++)
        {
           fprintf(file, "%e  ",mat->data[j][i]);
        }
        fprintf(file, "\n");
     }
     fprintf(file, "\n");
  }
}

void dmprintf(MATRIX *mat)
{
  unsigned int i, j;
  for(i=0;i<mat->row;i++)
  {
    for(j=0;j<mat->col;j++)
    {
      printf("%5.0f  ",mat->data[i][j]);
    }
    printf("\n");
  }
}


void dmfprintf(FILE *file, MATRIX *mat)
{
  unsigned int i, j;
  for(i=0;i<mat->row;i++)
  {
    for(j=0;j<mat->col;j++)
    {
      fprintf(file, "%5.0f  ",mat->data[i][j]);
    }
    fprintf(file, "\n");
  }
}


void make_1d_array_double(double **x, int dim, double ini)
{
  int j;
  if(dim > 0)
  {
     (*x) = (double *) malloc((size_t) (dim * sizeof(double)));
     for(j=0; j<dim; j++)  (*x)[j] = ini;
  }
  else  (*x) = NULL;  
}

void make_1d_array_int(int **x, int dim, int ini)
{
  int j;
  if(dim > 0)
  {
     (*x) = (int *) malloc((size_t) (dim * sizeof(int)));
     for(j=0; j<dim; j++)  (*x)[j] = ini;
  }
  else  (*x) = NULL;  
}

void reset_1d_array_double(double *x, int dim, double ini)
{
  int j;
  if(x != NULL)
  {
     for(j=0; j<dim; j++)  x[j] = ini;
  } 
}

void reset_1d_array_int(int *x, int dim, int ini)
{
  int j;
  if(x != NULL)
  {
     for(j=0; j<dim; j++)  x[j] = ini;
  } 
}


void make_2d_array_double(double ***x, int d1, int d2, double ini)
{
  int j, m;
  if(d1 > 0 && d2 > 0)
  {
     (*x) = (double **) malloc((size_t) (d1 * sizeof(double *)));
     (*x)[0] = (double *) malloc((size_t) ((d1 * d2) * sizeof(double)));
     for(j=1; j<d1; j++)  (*x)[j] = (*x)[j-1] + d2;
     m = d1 * d2;
     for(j=0; j<m; j++)  (*x)[0][j] = ini;
  } 
  else  (*x) = NULL;  
}

void make_2d_array_int(int ***x, int d1, int d2, int ini)
{
  int j, m;
  if(d1 > 0 && d2 > 0)
  {
     (*x) = (int **) malloc((size_t) (d1 * sizeof(int *)));
     (*x)[0] = (int *) malloc((size_t) ((d1 * d2) * sizeof(int)));
     for(j=1; j<d1; j++)  (*x)[j] = (*x)[j-1] + d2;
     m = d1 * d2;
     for(j=0; j<m; j++)  (*x)[0][j] = ini;
  } 
  else  (*x) = NULL;  
}

void reset_2d_array_double(double **x, int d1, int d2, double ini)
{
  int j, m;
  if(x != NULL)
  {
     m = d1 * d2;
     for(j=0; j<m; j++)  x[0][j] = ini;
  } 
}

void reset_2d_array_int(int **x, int d1, int d2, int ini)
{
  int j, m;
  if(x != NULL)
  {
     m = d1 * d2;
     for(j=0; j<m; j++)  x[0][j] = ini;
  }
}

void free_2d_array_double(double **x)
{
   if(x != NULL)
   {
      if(x[0] != NULL) free(x[0]);
      free(x);
   }
}

void free_2d_array_int(int **x)
{
   if(x != NULL)
   {
      if(x[0] != NULL) free(x[0]);
      free(x);
   }
}

void print_2d_int(int **x, int d1, int d2)
{
  unsigned int i, j;
  for(i=0;i<d1;i++)
  {
    for(j=0;j<d2;j++)
    {
      printf("%12d  ", x[i][j]);
    }
    printf("\n");
  }
}

void print_2d_double(double **x, int d1, int d2)
{
  unsigned int i, j;
  for(i=0;i<d1;i++)
  {
    for(j=0;j<d2;j++)
    {
      printf("%e  ", x[i][j]);
    }
    printf("\n");
  }
}

void make_3d_array_double(double ****x, int d1, int d2, int d3, double ini)
{
   int j, k, l, m;
   if(d1 > 0 && d2 > 0 && d3 > 0)
   {
      (*x) = (double ***) malloc((size_t) (d1 * sizeof(double **)));
      (*x)[0] = (double **) malloc((size_t) ((d1 * d2) * sizeof(double *)));
      for(j=1; j<d1; j++)
         (*x)[j] = (*x)[j-1] + d2;
      (*x)[0][0] = (double *) malloc((size_t) ((d1 * d2 * d3) * sizeof(double)));
      m = d1 * d2 * d3;
      for(j=0; j<m; j++)  (*x)[0][0][j] = ini;
      for(j=0; j<d1; j++)
      {
         for(k=0; k<d2; k++)
         {
            l = j * d2 + k;
            if(l > 0)  (* ((*x)[0] + l)) = (* ((*x)[0] + l - 1)) + d3;
         }
      }
   }
   else  (*x) = NULL;
}

void make_3d_array_int(int ****x, int d1, int d2, int d3, int ini)
{
   int j, k, l, m;
   if(d1 > 0 && d2 > 0 && d3 > 0)
   {
      (*x) = (int ***) malloc((size_t) (d1 * sizeof(int **)));
      (*x)[0] = (int **) malloc((size_t) ((d1 * d2) * sizeof(int *)));
      for(j=1; j<d1; j++)
         (*x)[j] = (*x)[j-1] + d2;
      (*x)[0][0] = (int *) malloc((size_t) ((d1 * d2 * d3) * sizeof(int)));
      m = d1 * d2 * d3;
      for(j=0; j<m; j++)  (*x)[0][0][j] = ini;
      for(j=0; j<d1; j++)
      {
         for(k=0; k<d2; k++)
         {
            l = j * d2 + k;
            if(l > 0)  (* ((*x)[0] + l)) = (* ((*x)[0] + l - 1)) + d3;
         }
      }
   }
   else  (*x) = NULL;
}

void reset_3d_array_double(double ***x, int d1, int d2, int d3, double ini)
{
  int j, m;
  if(x != NULL)
  {
     m = d1 * d2 * d3;
     for(j=0; j<m; j++)  x[0][0][j] = ini;
  }
}

void reset_3d_array_int(int ***x, int d1, int d2, int d3, int ini)
{
  int j, m;
  if(x != NULL)
  {
     m = d1 * d2 * d3;
     for(j=0; j<m; j++)  x[0][0][j] = ini;
  } 
}

void free_3d_array_double(double ***x)
{
   if(x != NULL)
   {
      if(x[0] != NULL)
      {
         if(x[0] != NULL)  free(x[0][0]);
         free(x[0]);
      }
      free(x);
   }
}

void free_3d_array_int(int ***x)
{
   if(x != NULL)
   {
      if(x[0] != NULL)
      {
         if(x[0] != NULL)  free(x[0][0]);
         free(x[0]);
      }
      free(x);
   }
}

