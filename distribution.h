#include <math.h>
#define runiform mt32real
#ifndef M_PI
#define M_PI 3.1415926535897932
#endif
#ifndef M_E
#define M_E 2.718281828
#endif

#ifndef max
#define max(x, y) ((x)>(y)?(x):(y))
#endif
#ifndef min
#define min(x, y) ((x)<(y)?(x):(y))
#endif

/*\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  Continuous distribution                          \\\\
  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/

/*================================================
  Uniform(0, 1) 
  ================================================*/

/* Uniform random number generator*/
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum)
{

int j;
long k;
static long idum2=123456789;
static long iy=0;
static long iv[NTAB];
double temp;

if(*idum <= 0)
{
    if(-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for(j=NTAB+7;j>=0;j--)
    {
        k=(*idum)/IQ1;
        *idum=IA1*(*idum-k*IQ1)-k*IR1;
        if(*idum < 0) *idum += IM1;
        if(j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
}
k=(*idum)/IQ1;
*idum=IA1*(*idum-k*IQ1)-k*IR1;
if(*idum < 0) *idum += IM1;
k=idum2/IQ2;
idum2=IA2*(idum2-k*IQ2)-k*IR2;
if(idum2 < 0) idum2 += IM2;
j=iy/NDIV;
iy=iv[j]-idum2;
iv[j] = *idum;
if(iy < 1) iy+= IMM1;
if ((temp=AM*iy) > RNMX) return RNMX;
else return temp;
}


/* Uniform random number generator, simple version*/
#define IA 16807
#define IM 2147483647
#ifdef AM
#undef AM
#endif
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

double ran1(long *idum)
{
  long k;
  double ans;

  *idum ^= MASK;
  k=(*idum)/IQ;
  *idum = IA*(*idum-k*IQ)-IR*k;
  if(*idum < 0) *idum += IM;
  ans = AM*(*idum);
  *idum ^= MASK;
  return(ans);
}

/**********************************************************************
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using mt_init(seed)
   or init_by_array(init_key, key_length).
************************************************************************/

#include <sys/types.h>
#include <inttypes.h>

/* Period parameters */
#define N_CON 624
#define M_CON 397
#define MATRIX_A   0x9908b0df /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

static unsigned int mt[N_CON]; /* the array for the state vector  */
static int mti; 

/* initializes mt[N] with a seed */
void mt_init(long seed)
{
    mt[0]= seed>0? seed:(- seed); /* & 0xffffffff; see comment [XU] below */
    for (mti=1; mti<N_CON; mti++) {
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] = (1812433253 * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);

        /* for >32 bit machines */
                /* not necessary when using uint32_t, and 32bit integer
                 * operations are faster than 64bit ones, so comment out
                 * the next statement [XU] */
        /* mt[mti] &= 0xffffffff;*/
    }
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned int mt32uint(long *seed)
{
    unsigned int y;
    unsigned int mag01[2]={0x0, MATRIX_A};
    int kk; 

    if (mti >= N_CON) 
    { /* generate N words at one time */

        for (kk=0;kk<N_CON - M_CON; kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M_CON] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<N_CON - 1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk + (M_CON - N_CON)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[N_CON-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N_CON-1] = mt[M_CON-1] ^ (y >> 1) ^ mag01[y & 0x1];

        mti = 0;
    }
    *seed = mt[mti];
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y <<  7) & 0x9d2c5680;
    y ^= (y << 15) & 0xefc60000;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on (0,1) with a 32-bit resolution */
double mt32real(long *seed)
{
    unsigned int y;
    unsigned int mag01[2]={0x0, MATRIX_A};
    int kk; 
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N_CON) 
    {

                /* generate N words at one time */
        for (kk=0; kk<N_CON - M_CON; kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk + M_CON] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<N_CON - 1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk + (M_CON - N_CON)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[N_CON-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N_CON-1] = mt[M_CON-1] ^ (y >> 1) ^ mag01[y & 0x1];

        mti = 0;
    }
    *seed = mt[mti];
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y <<  7) & 0x9d2c5680;
    y ^= (y << 15) & 0xefc60000;
    y ^= (y >> 18);

    return (y+0.5)*(1.0/4294967296.0);
    /* divided by 2^32 */
}


/*================================================
  Normal(mean, sigma^2) 
  ================================================*/

/* sampling from normal*/
double rnorm(double mean, double sigma, long *seed)
{
  double x, y, rsq, fac;
  static int iset = 0;
  static double extra;
  if(*seed < 0)
     iset = 0;
  if(iset == 0)
  {
    do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */

      x = -1.0 + 2.0 * runiform (seed);
      y = -1.0 + 2.0 * runiform (seed);

      /* see if it is in the unit circle */
      rsq = x * x + y * y;
    }
    while (rsq > 1.0 || rsq == 0);
    /* Box-Muller transform */
    fac =  sqrt (-2.0 * log (rsq) / rsq);
    extra = x*fac;
    iset = 1;
    return  mean + sigma*y*fac;
  }
  else
  {
    iset = 0;
    return mean + sigma*extra;
  }
}

/* normal density function*/
double dnorm(double t, double mean, double sigma) 
{ return (0.398942280401433/sigma)*exp(-(t-mean)*(t-mean)/(2.0*sigma*sigma)); }

/* normal cumulative distribution function*/
double pnorm(double x, double mean, double sigma) 
{ 
  double result; 
  static double a[5] = {0.31938153,-0.356563782,1.781477937,-1.821255978,1.330274429};
  x = (x-mean)/sigma;
  if (x<-7.0)
    result = dnorm(x, 0, 1)/sqrt(1.0+x*x);
  else if (x>7.0)
    result = 1.0 - pnorm(-x,0,1);
  else
  {
    result = 0.2316419;
    result=1.0/(1+result*fabs(x));
    result=1.0 - 0.398942280401433*exp(-x*x/2.0)
           *(result*(a[0]+result*(a[1]+result*(a[2]+result*(a[3]+result*a[4])))));
    if (x<=0.0)
      result=1.0-result;
  } 
  return result;
}

/* inverse of normal cumulative distribution function*/
double qnorm(double u, double mean, double sigma)
{
 /* returns the inverse of cumulative normal distribution function
  Reference> The Full Monte, by Boris Moro, Union Bank of Switzerland
                        RISK 1995(2)*/

  static double a[4]={2.50662823884, -18.61500062529, 41.39119773534, -25.44106049637};
  static double b[4]={-8.47351093090, 23.08336743743,-21.06224101826,   3.13082909833};
  static double c[9]={0.3374754822726147,
                                        0.9761690190917186,
                                        0.1607979714918209,
                                        0.0276438810333863,
                                        0.0038405729373609,
                                        0.0003951896511919,
                                        0.0000321767881768,
                                        0.0000002888167364,
                                        0.0000003960315187};
  double x,r;
  x=u-0.5;
  if (fabs(x)<0.42)
  {
    r=x*x;
    r=x*(((a[3]*r+a[2])*r+a[1])*r+a[0])/
       ((((b[3]*r+b[2])*r+b[1])*r+b[0])*r+1.0);
    return(r*sigma+mean);
  }
  r=u;
  if(x>0.0) 
    r=1.0-u;
  r=log(-log(r));
  r=c[0]+r*(c[1]+r*(c[2]+r*(c[3]+r*(c[4]+r*(c[5]+r*(c[6]+
                r*(c[7]+r*c[8])))))));
  if(x<0.0) r=-r;
  return(r*sigma+mean);
}

/*================================================
  Generate a bi-variate normal random value
  ================================================*/
/* sampling from binormal*/
void rbinorm (double mean_x, double sigma_x, double mean_y, double sigma_y, double rho,
                            double *x, double *y, long *seed)
{
  double u, v, r2, scale;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */

      u = -1 + 2 * runiform (seed);
      v = -1 + 2 * runiform (seed);

      /* see if it is in the unit circle */
      r2 = u * u + v * v;
    }
  while (r2 > 1.0 || r2 == 0);

  scale = sqrt (-2.0 * log (r2) / r2);

  *x = mean_x + sigma_x * u * scale;
  *y = mean_y + sigma_y * (rho * u + sqrt(1 - rho*rho) * v) * scale;
}

/*density of binormal*/
double dbinorm (double x, double y, double mean_x, double mean_y,  
               double sigma_x, double sigma_y, double rho)
{
  double u, v, c, p;
  u = (x-mean_x) / sigma_x ;
  v = (y-mean_y) / sigma_y ;
  c = 1 - rho*rho ;
  p = (1 / (2 * M_PI * sigma_x * sigma_y * sqrt(c))) 
    * exp (-(u * u - 2 * rho * u * v + v * v) / (2 * c));
  return p;
}

/*================================================
  Generate a normal random vector
  ================================================*/
/* sampling from multivariate normal*/
void rmultinorm (double *y, int n, double *mean, MATRIX *sigma, long *seed)
{
  int i, j, k;
  double *u;
  MATRIX L;
  
  initialize_matrix(&L);
  inflate_matrix(&L, n, n, 0.0);
  
  u = (double *) malloc((size_t) n * sizeof(double));
  /* 
  for(i=0; i<n; i++)
     u[i] = rnorm(0.0, 1.0, seed);
  */
  
  if((n%2) == 0)  
  {
     k = n/2;
     for(i=0; i<k; i++)
     {
        rbinorm (0.0, 1.0, 0.0, 1.0, 0.0, &u[2*i], &u[2*i+1], seed);
     }
  }
  else  
  {
     k = (n-1)/2;
     for(i=0; i<k; i++)
     {
        rbinorm (0.0, 1.0, 0.0, 1.0, 0.0, &u[2*i], &u[2*i+1], seed);
     }
     u[n-1] = rnorm(0.0, 1.0, seed);
  } 
    
  cholesky_decomp(sigma, &L);
  for(i=0; i<n; i++)
  {
     y[i] = 0.0;
     for(j=0; j<n; j++)
        y[i] += u[j] * L.data[i][j];
  }        
  for(i=0; i<n; i++)  y[i] += mean[i];
  
  free(u);
  deflate_matrix(&L);

}

/*density of binormal*/
double dmultinorm (double *y, int n, double *mean, MATRIX *sigma)
{
  int i, j;
  MATRIX inv;
  double *c_mean, d, log_f;
  
  initialize_matrix(&inv);
  inflate_matrix(&inv, n, n, 0.0);
  
  c_mean = (double *) malloc((size_t) n * sizeof(double));
  for(i=0; i<n; i++)  
     c_mean[i] = y[i] - mean[i];
  
  inverse(sigma, &inv);
  
  d = 0.0;
  for(i=0; i<n; i++)  
     for(j=0; j<n; j++)
        d += c_mean[i] * inv.data[i][j] * c_mean[j];
  
  log_f = -0.5 * (n * log(2.0 * 3.14159265359) + log(det_sympos(sigma)) + d);
  
  free(c_mean);
  deflate_matrix(&inv);
  return (exp(log_f));
}
/*================================================
  Lognormal(zeta, sigma)
  p(x) = 1/(x * sqrt(2 pi sigma^2)) exp(-(ln(x) - zeta)^2/2 sigma^2)
  ================================================*/
/* sampling lognormal*/
double rlognorm (double zeta, double sigma, long *seed)
{
  double u, v, r2, normal, z;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */

      u = -1 + 2 * runiform (seed);
      v = -1 + 2 * runiform (seed);

      /* see if it is in the unit circle */
      r2 = u * u + v * v;
    }
  while (r2 > 1.0 || r2 == 0);

  normal = u * sqrt (-2.0 * log (r2) / r2);

  z =  exp (sigma * normal + zeta);

  return z;
}

/* lognormal density*/
double dlognormal (double x, double zeta, double sigma)
{
  double u, p;
  if (x <= 0)
    {
      return 0 ;
    }
  else
    {
      u = (log (x) - zeta)/sigma;
      p = 1 / (x * fabs(sigma) * sqrt (2 * M_PI)) * exp (-(u * u) /2);
      return p;
    }
}

double plognorm(double x, double zeta, double sigma)
{ 
  double result; 
  static double a[5] = {0.31938153,-0.356563782,1.781477937,-1.821255978,1.330274429};
  x = (log(x) - zeta) / sigma;
  if (x<-7.0)
    result = dnorm(x, 0, 1)/sqrt(1.0+x*x);
  else if (x>7.0)
    result = 1.0 - pnorm(-x,0,1);
  else
  {
    result = 0.2316419;
    result=1.0/(1+result*fabs(x));
    result=1.0 - 0.398942280401433*exp(-x*x/2.0)
           *(result*(a[0]+result*(a[1]+result*(a[2]+result*(a[3]+result*a[4])))));
    if (x<=0.0)
      result=1.0-result;
  } 
  return result;
}

/*================================================
  Generate a log-normal random vector
  ================================================*/
/* sampling from multivariate normal*/
void rmultilognorm (double *y, int n, double *mean, MATRIX *sigma, long *seed)
{
  int i, j, k;
  double *u;
  MATRIX L;
  
  initialize_matrix(&L);
  inflate_matrix(&L, n, n, 0.0);
  
  u = (double *) malloc((size_t) n * sizeof(double));
  
  for(i=0; i<n; i++)
     u[i] = rnorm(0.0, 1.0, seed);
  cholesky_decomp(sigma, &L);
  for(i=0; i<n; i++)
  {
     y[i] = 0.0;
     for(j=0; j<n; j++)
        y[i] += u[j] * L.data[i][j];
  }        
  for(i=0; i<n; i++)  y[i] += mean[i];
  
  for(i=0; i<n; i++)
     y[i] = exp(y[i]);
  free(u);
  deflate_matrix(&L);

}

/*density of binormal*/
double dmultilognorm (double *y, int n, double *mean, MATRIX *sigma)
{
  int i, j;
  MATRIX inv;
  double *c_mean, d, log_f;
  
  initialize_matrix(&inv);
  inflate_matrix(&inv, n, n, 0.0);
  
  c_mean = (double *) malloc((size_t) n * sizeof(double));
  for(i=0; i<n; i++)  
     c_mean[i] = log(y[i]) - mean[i];
  
  inverse(sigma, &inv);
  
  d = 0.0;
  for(i=0; i<n; i++)  
     for(j=0; j<n; j++)
        d += c_mean[i] * inv.data[i][j] * c_mean[j];
  
  log_f = -0.5 * (n * log(2.0 * 3.14159265359) + log(det_sympos(sigma)) + d);
  for(i=0; i<n; i++)  
     log_f += - log(y[i]); 
  free(c_mean);
  deflate_matrix(&inv);
  return (exp(log_f));
}
/*================================================
  Exponential (lambda)
  p(x) = lambda * exp(-x*lambda)

  ================================================*/
/* sampling exponential*/
double rexponential(double lambda, long *seed)
{
  double x;
  if(lambda <= 0)
  {
    printf("parameter for exponential random variable must be positive\n");
    exit(0);
  }   
  do
  {
      x = runiform(seed);
  }
  while (x == 0);

  /* Box-Muller transform */
  return -log(x)/lambda;
}

/* density of exponential */
double dexponential (double x, double lambda)
{
  double p;
  if (x < 0)
    {
      return 0 ;
    }
  else
    {
      p = lambda * exp (-x*lambda);

      return p;
    }
}

/*CDF of exponential*/
double pexponential (double x, double lambda)
{
   if(x <= 0.0)  return (0.0);
   else   
     return InGamma(1, lambda * x);
}

/*================================================
  Gamma (a, b)
  p(x)  = {1 / [Gamma(a) b^a] } x^{a-1} e^{-x/b}
  ================================================*/
/* sampling gamma*/
double rgamma (double a, double b, long *seed)
{
  unsigned int i;
  unsigned int a_int;
  double a_frac, a_int_f;
  double pos_uniform, prod = 1.0;
  double am, p, q, s, x, y, u, v;
  double z_int, z_frac;

  if(a <= 0 || b <= 0)
  {
    printf("parameters for gamma random variable must be positive\n");
    exit(0);
  }
  a_int = floor (a); /* a_int is the integer part of a*/
  /* a_int_f is the floating version of a_int*/
  a_int_f = (double) a_int;
  if(a_int > 0)
  {
     if (a_int < 12)
     {
        for (i = 0; i < a_int; i++)
        {
           do
           {
              pos_uniform = runiform(seed);
           } while(pos_uniform == 0);
           prod *= pos_uniform;
        }

        /* Note: for 12 iterations we are safe against underflow, since
         the smallest positive random number is O(2^-32). This means
         the smallest possible product is 2^(-12*32) = 10^-116 which
         is within the range of double precision. */

        z_int = -log(prod);
     }
     else
     {
        am = a_int_f - 1.0;
        s = sqrt(2.0*am + 1.0);
        do
        {
      	    do
      	    {
               y = tan (M_PI * runiform(seed));
               x = s*y + am;
           } while(x <= 0.0);
           q = (1.0+y*y) * exp(am * log(x/am) - s*y);
        } while(runiform(seed) > q);
        z_int = x;
     }
  }
  /*if a is an integer*/
  if(a == a_int_f) 
    return b*z_int;

  a_frac = a - a_int_f; /*a_frac is the fraction part of a*/
  p = M_E / (a_frac + M_E);
  am = 1.0/a_frac;
  do
  {
      u = runiform (seed);
      do
      {
        v = runiform(seed);
      } while(v == 0);
      if (u < p)
        {
          x = exp (am * log (v));
          q = exp (-x);
        }
      else
        {
          x = 1.0 - log (v);
          q = exp ((a_frac - 1.0) * log (x));
        }
  }
  while (runiform (seed) >= q);

  z_frac = x;

  if(a_int == 0)
      return b*z_frac;
  else
      return b*(z_int + z_frac) ;
}

/*density of gamma*/
double dgamma (double x, double a, double b)
{
  double p, ln_gamma;
  if (x < 0)
    {
      return 0 ;
    }
  else if (x == 0)
    {
      if (a == 1.0)
        return 1/b ;
      else
        return 0 ;
    }
  else if (a == 1.0)
    {
      return exp(-x/b)/b ;
    }
  else
    {
      ln_gamma = lgamma(a);
      p = exp ((a - 1.0) * log (x/b) - x/b - ln_gamma)/b;
      return p;
    }
}

/*CDF of gamma*/
double pgamma (double x, double a, double b)
{
   if(x <= 0.0)  return (0.0);
   else   
     return InGamma(a, x/b);
}

/*================================================
   Beta (a, b)
   p(x) = (Gamma(a + b)/(Gamma(a) Gamma(b))) x^(a-1) (1-x)^(b-1)

  ================================================*/
/* sampling beta. when a is large and b is small, say, a=4 and b=0.07,
   it is possible that x1/(x1+x2) will be rounded to 1 by the computer.*/
double rbeta (double a, double b, long *seed)
{
  double x1;
  double x2;
  if(a <= 0 || b <= 0)
  {
    printf("parameters for beta random variable must be positive\n");
    exit(0);
  }
  x1 = rgamma(a, 1.0, seed);
  x2 = rgamma(b, 1.0, seed);
  return x1/(x1+x2);
}

/* density of beta*/
double dbeta (double x, double a, double b)
{
  double p, gab, ga, gb;
  if (x < 0 || x > 1)
    {
      return 0 ;
    }
  else
    {
      gab = lgamma (a + b);
      ga = lgamma (a);
      gb = lgamma (b);

      p = exp (gab - ga - gb) * pow (x, a - 1.0) * pow (1.0 - x, b - 1.0);

      return p;
    }
}

/* CDF of beta*/
double pbeta (double x, double a, double b)
{
  if (x <= 0.0) return(0.0);
  else if(x >= 1.0)   return(1.0);
  else
      return(InBeta (a, b, x));

}

/* Inverse CDF of beta*/
double qbeta_processor(double P, double a, double b)
{
   double x, lx, lg_a, lg_b, lg_ab, mean;
   double lambda, dP, phi;
   double step, step0, step1;
   unsigned int n = 0;

   if (P < 0.0 || P > 1.0)
   {
      printf("Error in inversing Beta CDF: P must be in range 0 < P < 1\n");
      exit(0);
   }

   if (a < 0.0 || b < 0.0)
   {
      printf("Error in inversing Beta CDF: wrong parameters a = %f  b=%f\n", a, b);
      exit(0);
   }

   if (P == 0.0)
   {
      return 0.0;
   } 

   if (P == 1.0)
   {
       return 1.0;
   }

   mean = a / (a + b);

   if (P < 0.1)
   {
      /* small x */

      lg_ab = lgamma (a + b);
      lg_a = lgamma (a);
      lg_b = lgamma (b);

      lx = (log (a) + lg_a + lg_b - lg_ab + log (P)) / a;
      x = exp (lx);             /* first approximation */
      x *= pow (1 - x, -(b - 1) / a);   /* second approximation */

      if (x > mean)
         x = mean;
   }
   else
   {
       /* Use expected value as first guess */
       x = mean;
   }

start:
   dP = P - pbeta(x, a, b);
   phi = dbeta(x, a, b);

   if (dP == 0.0 || n++ > 64)
      goto end;

   lambda = dP / max(2 * fabs (dP / x), phi);

   step0 = lambda;
   step1 = -((a - 1) / x - (b - 1) / (1 - x)) * lambda * lambda / 2;

   step = step0;

   if (fabs (step1) < fabs (step0))
   {
      step += step1;
   }
   else
   {
       /* scale back step to a reasonable size when too large */
       step *= 2 * fabs (step0 / step1);
   }

   if (x + step > 0 && x + step < 1)
   {
      x += step;
   }
   else
   {
       x = sqrt (x) * sqrt (mean);   /* try a new starting point */
   }

   if (fabs (step0) > 1e-10 * x)
      goto start;

end:
    return x;

}

double qbeta(double P, double a, double b)
{
   if (P <= 0.5)
   {
      return(qbeta_processor(P, a, b));
   }
   else
   {
      return(1.0 - qbeta_processor(1.0 - P, b, a));
   }   
}
/*================================================
  Generate a caucy random value.
  p(x) = (1/(pi a)) (1 + (x/a)^2)^(-1)
  ================================================*/

double rcauchy (double a, long *seed)
{
  double u;
  do
    {
      u = runiform (seed);
    }
  while (u == 0.5);

  return a * tan (M_PI * u);
}

double dcauchy (double x, double a)
{
  return (1 / (M_PI * a)) / (1.0 + (x * x)/(a * a));
}

/*================================================
  Chisquare (df=nu).
  p(x) = (1/(2*Gamma(nu/2))) (x/2)^(nu/2 - 1) exp(-x/2)
  ================================================*/
/* sampling chisquare*/
double rchisq (double nu, long *seed)
{
  double chisq = 2.0 * rgamma (nu / 2.0, 1.0, seed);
  return chisq;
}

/* density of chisquare*/
double dchisq (double x, double nu)
{
  double p, lngamma;
  if (x <= 0)
    {
      return 0 ;
    }
  else
    {
      lngamma = lgamma (nu / 2.0);

      p = exp ((nu / 2.0 - 1.0) * log (x/2.0) - x/2.0 - lngamma) / 2.0;
      return p;
    }
}

/* cdf of chisquare*/
double pchisq(double nu, double x)
{
   return (InGamma(nu / 2.0, x / 2.0));
}

/*================================================
  F(nu1, nu2)
  p(x) = (nu1^(nu1/2) nu2^(nu2/2) Gamma((nu1 + nu2)/2) /
   (Gamma(nu1/2) Gamma(nu2/2))) *
   x^(nu1/2 - 1) (nu2 + nu1 * x)^(-nu1/2 -nu2/2)
  ================================================*/
/* sampling F*/
double fdist (double nu1, double nu2, long *seed)
{

  double Y1 ;
  double Y2 ;

  if(nu1 <= 0 || nu2 <= 0)
  {
    printf("degree of freedom for F random variable must be positive\n");
    exit(0);
  }
  Y1 =  rgamma (nu1 / 2.0, 2.0, seed);
  Y2 =  rgamma (nu2 / 2.0, 2.0, seed);

  return (Y1 * nu2) / (Y2 * nu1);
}

/* density of F*/
double dfdist (double x, double nu1, double nu2)
{
  double p, lglg, lg12, lg1, lg2;
  if (x < 0)
    {
      return 0 ;
    }
  else
    {
      lglg = (nu1 / 2.0) * log (nu1) + (nu2 / 2.0) * log (nu2) ;

      lg12 = lgamma ((nu1 + nu2) / 2.0);
      lg1 = lgamma (nu1 / 2.0);
      lg2 = lgamma (nu2 / 2.0);

      p = exp (lglg + lg12 - lg1 - lg2)
        * pow (x, nu1 / 2.0 - 1.0) * pow (nu2 + nu1 * x, -nu1 / 2.0 - nu2 / 2.0);

      return p;
    }
}

/*================================================
  Generate a T(nu) random value.
  p(x) = (Gamma((nu + 1)/2)/(sqrt(pi nu) Gamma(nu/2))
   * (1 + (x^2)/nu)^-((nu + 1)/2)
  ================================================*/
/* sampling T*/
double rtdist (double nu, long *seed)
{
  double Y1, Y2, t, Z;
  if (nu <= 2.0)
    {
      Y1 = rnorm (0.0,1.0,seed);
      Y2 = rchisq (nu, seed);

      t = Y1 / sqrt (Y2 / nu);

      return t;
    }
  else
    {
      do
      {
          Y1 = rnorm (0.0,1.0,seed);
          Y2 = rexponential (1.0 / (nu/2.0 - 1.0), seed);

          Z = Y1 * Y1 / (nu - 2.0);
      }
      while (1.0 - Z < 0.0 || exp (-Y2 - Z) > (1.0 - Z));

      /* Note that there is a typo in Knuth's formula, the line below
         is taken from the original paper of Marsaglia, Mathematics of
         Computation, 34 (1980), p 234-256 */

      t = Y1 / sqrt ((1.0 - 2.0 / nu) * (1.0 - Z));
      return t;
    }
}

/* density of T*/
double dtdist (double x, double nu)
{
  double p, lg1, lg2;

  lg1 = lgamma (nu / 2.0);
  lg2 = lgamma ((nu + 1.0) / 2.0);

  p = ((exp (lg2 - lg1) / sqrt (M_PI * nu))
       * pow ((1.0 + x * x / nu), -(nu + 1.0) / 2.0));
  return p;
}

   double ptdist(long n, double x)
/* =================================== 
 * NOTE: use n >= 1 and x > 0.0 
 * ===================================
 */
{ 
   double s, t;

   t = (x * x) / (n + x * x);
   s = InBeta(0.5, n / 2.0, t);
   if (x >= 0.0)
     return (0.5 * (1.0 + s));
   else
     return (0.5 * (1.0 - s));
}

/*================================================
  Generate multivariate T distribution random vector.
  p(x) = Gamma((n + p)/2
        --------------------------------------------------------------------------
        Gamma(n/2) (n PI)^{p/2} |Sigma|^{1/2} [1+(1/n) (x-u)'Sigma^(-1) (x-u)]^((n+p)/2)
  ================================================*/
void rmultitdist (double *y, int n, double df, double *mean, MATRIX *sigma, long *seed)
{
   int i;
   double *zero, x;
   zero = (double *) malloc((size_t) n * sizeof(double));
   for(i=0; i<n; i++)  zero[i] = 0.0;
   rmultinorm (y, n, zero, sigma, seed);
   x = rchisq (df, seed);
   for(i=0; i<n; i++)
      y[i] = y[i] / sqrt(x/df) + mean[i];

   free(zero);
}
   
/*================================================
  Generate a Weibull random value. 
  p(x) = (b/a) (x/a)^(b-1) exp(-(x/a)^b) 
  ================================================*/
/* sampling weibull*/
double rweibull (double a, double b, long *seed)
{
  double x, z;
  do
  {
    x = runiform (seed);
  }  
  while( x == 0.0);
  z = pow (-log (x), 1.0 / b);

  return a * z;
}

/* density of weibull*/
double dweibull (double x, double a, double b)
{
  double p;
  if (x < 0.0)
    {
      return 0.0 ;
    }
  else if (x == 0.0)
    {
      if (b == 1.0)
        return 1.0/a ;
      else
        return 0.0 ;
    }
  else if (b == 1.0)
    {
      return exp(-x/a)/a ;
    }
  else
    {
      p = (b/a) * exp (-pow (x/a, b) + (b - 1.0) * log (x/a));
      return p;
    }
}

/*================================================
  Generate a bi-variate student-T random value
  ================================================*/
/* sampling from binormal*/
void rbitdist (double mean_x, double mean_y, double sigma_x, double sigma_y, double rho, double df,
                            double *x, double *y, long *seed)
{
  double u, v, z, r2, scale;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */

      u = -1 + 2 * runiform (seed);
      v = -1 + 2 * runiform (seed);

      /* see if it is in the unit circle */
      r2 = u * u + v * v;
    }
  while (r2 > 1.0 || r2 == 0);

  scale = sqrt (-2.0 * log (r2) / r2);

  *x = sigma_x * u * scale;
  *y = sigma_y * (rho * u + sqrt(1 - rho*rho) * v) * scale;

  z = rchisq(df, seed);

  *x = (*x) / sqrt(z/df) + mean_x;
  *y = (*y) / sqrt(z/df) + mean_y;
}

/*density of bi-variate */
double dbitdist(double x, double y, double mean_x, double mean_y,
               double sigma_x, double sigma_y, double rho, double df)
{
  double u, v, c, p;

  u = (x-mean_x) / sigma_x ;
  v = (y-mean_y) / sigma_y ;
  c = 1 - rho*rho ;
  
  p = lgamma(0.5 * df + 1.0) - lgamma(0.5 * df) - log(df * M_PI)
      - 0.5 * log(sigma_x * sigma_x * sigma_y * sigma_y * c)
      - (0.5 * df + 1.0) * log(1.0 + (u * u - 2 * rho * u * v + v * v) / (df * c) );
  return(exp(p));
}

/*\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  Discrete random sampling                          \\\\
  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/


/*================================================
  Binomial (n, p) 
  ================================================*/
/* sampling binomial*/
/*this is numerical recipe algorithm*/
unsigned int rbinomial (unsigned int n, double pp, long *seed)
{
  unsigned int i, j;
  unsigned int bnl;
  double am, em, en, g, oldg, p, q, s, y, t;
  double log_p, log_q;

  if(pp == 0.0)  return 0;
  if(pp == 1.0)  return n;
  if(pp < 0 || pp > 1)
  {
    printf("parameter p for binomial random variable must be in (0, 1)\n");
    printf("p=%f\n", p);
    exit(0);
  }
  if( n == 0)  return 0;
  p = (pp <= 0.5)? pp:1.0 - pp;
  en = (double) n;
  am = en*p;
  if(n < 25)        /* This parameter is tunable */
  {
      bnl = 0;
      for(j = 1; j <= n; j++)
         if( runiform(seed) < p)  ++bnl;
  }
  else if(am < 1.0)
  {
      g = exp(-am);
      t = 1.0;
      for(j = 0; j <=n; j++)
      {
      	 t *= runiform(seed);
      	 if(t < g) break;
      }
      bnl = (j <= n)? j:n;
  }
  else
  {
      oldg = lgamma(en+1.0);
      q = 1.0-p;
      log_p = log(p);
      log_q = log(q);
      s = sqrt(2.0*am*q);
      do
      {
      	 do
      	 {
      	    y = tan(M_PI * runiform(seed));
      	    em = s * y + am;
      	 } while(em < 0.0 || em >= (en + 1.0));
         em = floor(em);
         t = 1.2 * s * (1.0+y*y) * exp(oldg - lgamma(em+1.0)
            - lgamma(en-em+1.0) + em * log_p + (en-em) * log_q);
      } while( runiform(seed) > t);
      bnl = (unsigned int) em;
  }
  if(p != pp)  bnl = n - bnl;

  return bnl;
}
/* mass probability function of binomial*/
double dbinomial(unsigned int n, double p, unsigned int k)
{
  double P, ln_Cnk;
  if (k > n)
    {
      return 0;
    }
  else
    {
      ln_Cnk = lchoose (n, k);
      P = ln_Cnk + k * log (p) + (n - k) * log (1 - p);
      P = exp (P);

      return P;
    }
}

/* cumulative probability function of binomial*/
double pbinomial(unsigned int n, double p, unsigned int x)
/* ============================================ 
 * NOTE: use 0 <= x <= n and 0.0 < p < 1.0 
 * ============================================
 */
{
   if (x < n)
     return 1.0 - InBeta(x + 1.0, (double)(n - x), p);
   else
     return 1.0;
}


/*================================================
  Generate a bernoulli random value. 
  ================================================*/

unsigned int rbernoulli (double p, long *seed)
{
  double u = runiform (seed) ;

  if(p <= 0 || p >= 1)
  {
    printf("parameter p for bernoulli random variable must be in (0, 1)\n");
    exit(0);
  }   
  if (u < p)
    {
      return 1 ;
    }
  else
    {
      return 0 ;
    }
}


/*================================================
  Generate a geometric random value. 
  prob(k) =  p (1 - p)^(k-1) for n = 1, 2, 3, ...
  ================================================*/
/* sampling geometric*/
unsigned int rgeometric (double p, long *seed)
{
  double u;

  unsigned int k;

  do
  {
    u = runiform(seed);
  }
  while(u == 0);
  if (p == 1)
    {
      k = 1;
    }
  else
    {
      k = log (u) / log (1 - p) + 1;
    }

  return k;
}

/* mass prob of geometric*/
double dgeometric (unsigned int k, double p)
{
  double P;
  if (k == 0)
    {
      return 0 ;
    }
  else if (k == 1)
    {
      return p ;
    }
  else
    {
      P = p * ipow (1 - p, k - 1);
      return P;
    }
}
/*================================================
  Generate a hypergeometric random value. 
  prob(k) =  choose(n1,t) choose(n2, t-k) / choose(n1+n2,t)

   where choose(a,b) = a!/(b!(a-b)!) 

   n1 + n2 is the total population (tagged plus untagged)
   n1 is the tagged population
   t is the number of samples taken (without replacement)
   k is the number of tagged samples found
  ================================================*/
  /* sampling hyper geometric*/
unsigned int rhyperg (unsigned int n1, unsigned int n2, 
                        unsigned int t, long *seed)
{
  const unsigned int n = n1 + n2;

  unsigned int i;
  unsigned int a = n1;
  unsigned int b = n1 + n2;
  unsigned int k = 0;
  double u;

  if (t > n)
    {
      t = n ;
    }

  if (t < n / 2) 
    {
      for (i = 0 ; i < t ; i++)
        {
          u = runiform(seed) ;
          
          if (b * u < a)
            {
              k++ ;
              if (k == n1)
                return k ;
              a-- ;
            }
          b-- ;
        }
      return k;
    }
  else
    {
      for (i = 0 ; i < n - t ; i++)
        {
          u = runiform(seed) ;
          
          if (b * u < a)
            {
              k++ ;
              if (k == n1)
                return n1 - k ;
              a-- ;
            }
          b-- ;
        }
      return n1 - k;
    }
}

/* mass prob of hypergeometric*/
double dhyperg (unsigned int k, unsigned int n1, 
                unsigned int n2, unsigned int t)
{
  double p, c1, c2, c3;
  if (t > n1 + n2)
    {
      t = n1 + n2 ;
    }

  if (k > n1 || k > t)
    {
      return 0 ;
    }
  else if (t > n2 && k + n2 < t )
    {
      return 0 ;
    }
  else 
    {
      c1 = lchoose(n1,k);
      c2 = lchoose(n2,t-k);
      c3 = lchoose(n1+n2,t);

      p = exp(c1 + c2 - c3) ;

      return p;
    }
}

/*================================================
  Generate a multinomoial random value. 
                                      N!           n_1  n_2      n_K
   prob(n_1, n_2, ... n_K) = -------------------- p_1  p_2  ... p_K
                             (n_1! n_2! ... n_K!) 
  ================================================*/
/* sampling multinomial*/
void rmultinomial (unsigned int K, unsigned int N, double p[],
                   unsigned int n[], long *seed)
{
  unsigned int k;
  double norm = 0.0;
  double sum_p = 0.0;
  double temp;

  unsigned int sum_n = 0;

  /* p[k] may contain non-negative weights that do not sum to 1.0.
   * Even a probability distribution will not exactly sum to 1.0
   * due to rounding errors. 
   */

  for (k = 0; k < K; k++)
    {
      norm += p[k];
      n[k] = 0;
    }

  for (k = 0; k < K; k++)
    {
      if (p[k] > 0.0)
        {
          temp = p[k] / (norm - sum_p);
          if(temp > 1.0)  temp = 1.0;
          n[k] = rbinomial (N - sum_n, temp, seed);
        }
      else
        {
          n[k] = 0;
        }

      sum_p += p[k];
      sum_n += n[k];
    }

}

/*mass prob of multinomial*/
double dmultinomial (unsigned int K, double p[], unsigned int n[])
{
  unsigned int k;
  unsigned int N = 0;
  double log_pdf;
  double norm = 0.0;

  for (k = 0; k < K; k++)
    {
      N += n[k];
    }

  for (k = 0; k < K; k++)
    {
      norm += p[k];
    }

  /* Note: n! == gamma(n+1) */
  log_pdf = lgamma (N + 1);

  for (k = 0; k < K; k++)
    {
      log_pdf -= lgamma (n[k] + 1);
    }

  for (k = 0; k < K; k++)
    {
      log_pdf += log (p[k] / norm) * n[k];
    }

  return exp(log_pdf);
}

/*================================================
  Generate a Poisson random value.
  p(n) = (mu^n / n!) exp(-mu)
  ================================================*/
/*this is numerical recipe algorithm*/
unsigned int rpoisson (double mu, long *seed)
{
   double s, alxm, g;
   double em, t, y;
   if(mu < 12.0)
   {
        g = exp(-mu);
   	em = -1.0;
   	t = 1.0;
   	do
   	{
   	   ++em;
   	   t *= runiform(seed);
   	} while(t > g);
    }
    else
    {
    	s = sqrt(2.0*mu);
    	alxm = log(mu);
    	g = mu * alxm - lgamma(mu+1.0);
        do
        {
           do
           {
              y = tan(M_PI * runiform(seed));
              em = s * y + mu;
           } while(em < 0.0);
           em = floor(em);
           t = 0.9 * (1.0 + y*y) * exp(em * alxm - lgamma(em+1.0) - g);
        } while(runiform(seed) > t);
    }
    return (unsigned int) em;
}

/*This is gsl algorithm*/
/*
unsigned int rpoisson (double mu, long *seed)
{
  double lambda, emu;
  double prod = 1.0;
  unsigned int m, k = 0;
  double X;
  lambda = mu;
  while (lambda > 10)
    {
      m = (unsigned int) (lambda * (7.0 / 8.0));

      X = rgamma(m,1.0,seed);

      if (X >= lambda)
        {
          return k + rbinomial (m-1, lambda / X, seed);
        }
      else
        {
          k += m;
          lambda -= X;
        }
    }

  emu = exp (-lambda);

  do
    {
      prod *= runiform (seed);
      k++;
    }
  while (prod > emu);

  return k - 1;
}
*/
double dpoisson (unsigned int k, double mu)
{
  double p;
  double lf = lgamma((double)k+1.0);

  p = exp (log (mu) * k - lf - mu);
  return p;
}

double ppoisson(double lambda, unsigned int x)
/* ===================================
 * NOTE: use m > 0 and x >= 0
 * ===================================
 */
{
   return (1.0 - InGamma(x + 1.0, lambda));
}

/*================================================
  Generate a negative binomial random value.
   prob(k) =  Gamma(n + k)/(Gamma(n) Gamma(k + 1))  p^n (1-p)^k
   for k = 0, 1, ... . Note that n does not have to be an integer.
  ================================================*/

unsigned int rnbinomial (double n, double p, long *seed)
{
  double X = rgamma (n, 1.0, seed) ;
  unsigned int k = rpoisson (X*(1-p)/p, seed) ;
  return k ;
}

double dnbinomial (unsigned int k, double n, double p)
{
  double P, f, a, b;

  f = lgamma ((double) (k + n)) ;
  a = lgamma ((double) n) ;
  b = lgamma ((double)k + 1.0) ;

  // n could be fractional, so we use pow instead of ipow.
  P = exp(f-a-b) * pow (p, n) * pow (1 - p, k);

  return P;
}

/*==========================================================================
  Generate a dirichlet random value. 
  The Dirichlet probability distribution of order K-1 is 

     p(\theta_1,...,\theta_K) d\theta_1 ... d\theta_K = 
        (1/Z) \prod_i=1,K \theta_i^{alpha_i - 1} \delta(1 -\sum_i=1,K \theta_i)

   The normalization factor Z can be expressed in terms of gamma functions:

      Z = {\prod_i=1,K \Gamma(\alpha_i)} / {\Gamma( \sum_i=1,K \alpha_i)}  

   The K constants, \alpha_1,...,\alpha_K, must be positive. The K parameters, 
   \theta_1,...,\theta_K are nonnegative and sum to 1.

   The random variates are generated by sampling K values from gamma
   distributions with parameters a=\alpha_i, b=1, and renormalizing. 
   See A.M. Law, W.D. Kelton, Simulation Modeling and Analysis (1991)  
   ========================================================================*/

void rdirichlet (unsigned int K,
                   double alpha[], double theta[], long *seed)
{
  unsigned int i;
  double norm = 0.0;

  for (i = 0; i < K; i++)
    {
      theta[i] = rgamma (alpha[i], 1.0, seed);
    }

  for (i = 0; i < K; i++)
    {
      norm += theta[i];
    }

  for (i = 0; i < K; i++)
    {
      theta[i] /= norm;
    }
}

/* dirichlet density*/
double ddirichlet (int K, double alpha[], double theta[])
{
  int i;
  double sum_alpha, sum_log, sum_gamma;

  sum_alpha = sum_log = sum_gamma = 0.0;
  for (i = 0; i < K; i++)
  {
      sum_alpha += alpha[i];
      sum_log += (alpha[i] - 1.0) * log(theta[i]);
      sum_gamma += lgamma(alpha[i]);
  }

  return( exp(sum_log + lgamma(sum_alpha) - sum_gamma) );

}

unsigned int discrete(unsigned int K, double p[], long *seed)
{
   unsigned int i, k;
   double norm, u, sum=0.0;
   double *q;
   q = (double *) malloc((size_t) K * sizeof(double));
   if(p == NULL) for(k=0; k<K; k++) q[k] = 1;
   else  for(k=0; k<K; k++) q[k] = p[k];
   norm = 0.0;
   for (k = 0; k < K; k++)
   {
      norm += q[k];
   }
   
   u = norm * runiform(seed);
   for(i=0; i<K; i++)
   {
     sum += q[i];
     if(u < sum)
     {
        free(q);
        return(i);
     }
   }
   free(q);
   return(K-1);
}        

/* very slow version, but provide some hint for sorting
void permute(unsigned int n, unsigned int *array, long *seed)
{
  int i,j;
  int n_small, relocate;
  float current, *ran;
  int *rank, *small;

  ran = (float *)malloc((size_t) (n*sizeof(float)));
  rank = (int *)malloc((size_t) (n*sizeof(int)));
  small = (int *)malloc((size_t) (n*sizeof(int)));
  
  for(i=0; i<n; i++)
  {
     ran[i] = runiform(seed);
     small[i] = 0;
  } 
  for(i=0; i<n; i++)
  {
     current = ran[i];
     n_small = small[i];
     for(j=i+1; j<n; j++)
     {
	if(current >= ran[j]) 
	   n_small++;
	else   
	   small[j]++;
     }
     array[i] = n_small;
  }
  free(ran);
  free(rank);
  free(small);
}  

struct ring {
   struct ring *next;
   struct ring *prev;
   struct ring *next2;
   struct ring *prev2;
   unsigned int index;
};


void permute(unsigned int n, unsigned int *array, long *seed)
{
  unsigned int i, j, odd;
  unsigned int ringsize, locate, relocate, fastmove;
  struct ring *head, *p0, *p1, *p2, *p3, *p4;

  head = (struct ring *) malloc(sizeof(struct ring));
  head->index=0;
  p1 = head;
  p0 = NULL;
  for(i=1; i<n; i++)
  {
     p1->next  = (struct ring *) malloc(sizeof(struct ring));
     
     p2 = p1->next;
     p2->index = i;
     p2->prev = p1;
     p2->prev2 = p0;
     
     if(p0 != NULL) p0->next2 = p2;
     
     p0 = p1;
     p1 = p2;
  } 
  p1->next = head;
  p1->next2 = head->next;
  p0->next2 = head;
  head->prev = p1;
  head->prev2 = p0;
  head->next->prev2 = p1;
  
  ringsize = n;
  i = 0;
  while(ringsize > 1)
  {
     if(ringsize > 100)
     {
        locate = floor(ringsize * runiform(seed));
        p2 = head;
        if(2 * locate <= ringsize)
        {
           odd = locate%2;
           fastmove = (locate - odd) / 2;
           for(j=0; j<fastmove; j++) p2 = p2->next2;
           if(odd)   p2 = p2->next;
        }
        else   
        {
           relocate = ringsize - locate;
           odd = relocate%2;
           fastmove = (relocate - odd) / 2 ;
           for(j=0; j<fastmove; j++) p2 = p2->prev2;
           if(odd)   p2 = p2->next;
        }
     
        p0 = p2->prev2;
        p1 = p2->prev;
        p3 = p2->next;
        p4 = p2->next2;
          
        array[i++] = p2->index;
        p0->next2 = p3;
        p1->next = p3;
        p1->next2 = p4;
        p3->prev = p1;
        p3->prev2 = p0;
        p4->prev2 = p1;
     
        head = p3;
        free(p2);
        ringsize--;
     }   
     else
     {
        locate = floor(ringsize * runiform(seed));
     
        p2 = head;
        if(2 * locate <= ringsize)
        {
           for(j=0; j<locate; j++) p2 = p2->next;
           
        }
        else   
        {
           relocate = ringsize - locate;
           for(j=0; j<relocate; j++) p2 = p2->prev;
        }
        p1 = p2->prev;
        p3 = p2->next;
     
        array[i++] = p2->index;
        p1->next = p3;
        p3->prev = p1;
     
        head = p3;
        free(p2);
        ringsize--;
     }
  }
  array[n-1] = head->index; 
  free(head);
}  
*/
/* very fast version*/
void permute(unsigned int pop_size, unsigned int sample_size, unsigned int *array, long *seed)
{
  unsigned int i, j, k, temp;
  unsigned int *pop;
 
  pop = (unsigned int *) malloc((size_t) (pop_size*sizeof(unsigned int)));
  for(i=0; i<pop_size; i++)  pop[i] = i;
  i = pop_size;
  k = 0;
  while(i > 0 && k < sample_size)
  {
     j = floor(i * runiform(seed));
     temp = pop[i-1];
     pop[i-1] = pop[j];
     pop[j] = temp;
     array[k] = pop[i-1];
     i--;
     k++;
  }
  free(pop);  
}


void bootstrap(unsigned int pop_size, unsigned int sample_size, unsigned int *array, long *seed)
{
  int i;
                                                                                                                                                             
  for(i=0; i<sample_size; i++)
     array[i] = floor(pop_size * runiform(seed));
}

void combination(int n_factor, int *dim, int **input, int **output)
{
   int i, j, k, l;
   int prior_block_size, after_block_size, full_size;
   full_size = 1;
   for(i=0; i<n_factor; i++)  full_size *= dim[i];
   for(i=0; i<n_factor; i++)
   {
      prior_block_size = after_block_size = 1;
      for(j=i; j<n_factor; j++)  
      {
         prior_block_size *= dim[j];
         if(j > i)  after_block_size *= dim[j];
      }   
      for(j=0; j<full_size; j++) 
      {
         k = j % prior_block_size;
         l = floor(k/after_block_size);
         output[i][j] = input[i][l]; 
      }   
   }
}
   
/*****************************************************************************************
  output is a vector of sorted values. direction denotes ascending (1) or descending (0). 
  old_position is a vector of the original subscript of the sorted values.
  This method is not very efficient, about N*(N+1)/2.
  quick_sort is a faster sorting algorithm from numerical recipe.
 *****************************************************************************************/
void sort_vector(unsigned int n, double *input, unsigned int direction, double *output, unsigned int *position)
{
  unsigned int i,j;
  unsigned int n_small, relocate;
  double current;
  unsigned int *small;
  
  small = (unsigned int *)malloc((size_t) (n*sizeof(unsigned int)));
  
  for(i=0; i<n; i++)
  {
     small[i] = 0;
  } 
  for(i=0; i<n; i++)
  {
     current = input[i];
     n_small = small[i];
     for(j=i+1; j<n; j++)
     {
	if(current >= input[j]) 
	   n_small++;
	else   
	   small[j]++;
     }
     if(direction == 0) /*descending*/
     {
        position[n-1-n_small] = i;
        output[n-1-n_small] = current;
     }
     else /*ascending*/
     {
        position[n_small] = i;
        output[n_small] = current;
     }
        
  }
  free(small);
}  

/*****************************************************************************************
  quick_sort is a faster sorting algorithm from numerical recipe. direction=0 means
  descending, and 1 means ascending. "arr" is the input array, and will be sorted
  itself. The original positions are stored in "brr". "brr" should store the original
  cardinal order first, i.e., from 1 to n.
 *****************************************************************************************/
#define fswab(a, b)  ftemp=(a); (a) = (b); (b) = ftemp;
#define iswab(a, b)  itemp=(a); (a) = (b); (b) = itemp;
#define SMALL_SIZE  7
#define STACK_SIZE  500


void quick_sort(unsigned int direction, unsigned int n, double *arr, int *brr)
{
   long h, i, ir, j, k, l, *istack, jstack=-1;
   double a, ftemp;
   unsigned int b, itemp;
   FILE *file;
   
   ir = n-1;
   l = 0;
   istack = (long *) malloc((size_t) STACK_SIZE*sizeof(long));
   
   if(direction == 1) /*if ascending */
   for(; ;)
   {
      if(ir-l < SMALL_SIZE)
      {
         /*printf("start  l:%d, ir:%d\n", l, ir);*/
         for(j=l+1; j<=ir; j++)
         {
            a = arr[j];
            b = brr[j];
            for(h=j-1; h>=l; h--)
            {
               if(arr[h] <= a)  break;
               arr[h+1] = arr[h];
               brr[h+1] = brr[h];
            }
            arr[h+1] = a;
            brr[h+1] = b;
         }
         if(jstack <= 0)  break;
         ir = istack[jstack];
         l = istack[jstack-1];
         jstack -= 2;
         /*printf("get from stack  l:%d, ir:%d\n", l, ir);*/
      }
      else
      {
         /*printf("start  l:%d, i:%d, j:%d, ir:%d\n", l, i, j, ir);*/
         k = (l+ir) >> 1;
         
         fswab(arr[k], arr[l+1])
         iswab(brr[k], brr[l+1])
         
         if(arr[l] > arr[ir])
         {
            fswab(arr[l], arr[ir])
            iswab(brr[l], brr[ir])
         }
         if(arr[l+1] > arr[ir])
         {
            fswab(arr[l+1], arr[ir])
            iswab(brr[l+1], brr[ir])
         }
         if(arr[l] > arr[l+1])
         {
            fswab(arr[l], arr[l+1])
            iswab(brr[l], brr[l+1])
         } 
           
         i = l+1;
         j = ir;
         a = arr[l+1];
         b = brr[l+1];
         for(; ;)
         {
            do  i++; while(arr[i] < a);
            do  j--; while(arr[j] > a);
            if(j < i) break;
            fswab(arr[i], arr[j])
            iswab(brr[i], brr[j])
         }
         /*printf("partition at  l=%d  i=%d  j=%d  ir=%d\n", l, i, j, ir);*/
         arr[l+1] = arr[j];
         arr[j] = a;
         brr[l+1] = brr[j];
         brr[j] = b;
         jstack += 2;
         if(jstack >= STACK_SIZE)
         {
            printf("stacj size is too small for sorting\n");
            exit(0);
         }
         if(ir-i+1 >= j-l)
         {
            istack[jstack] = ir;
            istack[jstack-1] = i;
            /*printf("put on stack  %d, %d  ", i, ir);*/
            ir = j-1;
            /*printf("leave for partition  %d, %d\n", l, j-1);*/
         }
         else
         {
            istack[jstack] = j-1;
            istack[jstack-1] = l;
            /*printf("put on stack  %d, %d  ", l, j-1);*/
            l = i;
            /*printf("leave for partition  %d, %d\n", i, ir);*/
         }
         
      }
   }
   else  /* if descending*/
   for(; ;)
   {
      if(ir-l < SMALL_SIZE)
      {
         for(j=l+1; j<=ir; j++)
         {
            a = arr[j];
            b = brr[j];
            for(h=j-1; h>=l; h--)
            {
               if(arr[h] >= a)  break;
               arr[h+1] = arr[h];
               brr[h+1] = brr[h];
            }
            arr[h+1] = a;
            brr[h+1] = b;
         }
         if(jstack == -1)  break;
         ir = istack[jstack];
         l = istack[jstack-1];
         jstack -= 2;
      }
      else
      {
         k = (l+ir) >> 1;
         fswab(arr[k], arr[l+1])
         iswab(brr[k], brr[l+1])
         if(arr[l] < arr[ir])
         {
            fswab(arr[l], arr[ir])
            iswab(brr[l], brr[ir])
         }
         if(arr[l+1] < arr[ir])
         {
            fswab(arr[l+1], arr[ir])
            iswab(brr[l+1], brr[ir])
         }
         if(arr[l] < arr[l+1])
         {
            fswab(arr[l], arr[l+1])
            iswab(brr[l], brr[l+1])
         }   
         i = l+1;
         j = ir;
         a = arr[l+1];
         b = brr[l+1];
         for(; ;)
         {
            do  i++; while(arr[i] > a);
            do  j--; while(arr[j] < a);
            if(j < i)  break;
            fswab(arr[i], arr[j])
            iswab(brr[i], brr[j])
         }
         arr[l+1] = arr[j];
         arr[j] = a;
         brr[l+1] = brr[j];
         brr[j] = b;
         jstack += 2;
         if(jstack >= STACK_SIZE)
         {
            printf("stacj size is too small for sorting\n");
            exit(0);
         }
         if(ir-i+1 >= j-l)
         {
            istack[jstack] = ir;
            istack[jstack-1] = i;
            ir = j-1;
         }
         else
         {
            istack[jstack] = j-1;
            istack[jstack-1] = l;
            l = i;
         }
      }
   }
   free(istack);
}  

/***************************************************
 Binary searching in a sorted array. It actually
 returns the number of values in array a smaller than
 x. So the number of values at least as extreme (large)
 as x is n subtracting the returned number. This is useful
 for calculating p-value. 
***************************************************/
int binary_search(double x,  double *a, int n)
{
   int i, j;
   int head, rear, middle;
   
   head = 0;
   rear = n - 1;
   if(x <= a[head])  return 0;
   if(x > a[rear]) return n;

   while(1)
   {
      if(head == rear - 1)
      {
         if(x < a[head] || x > a[rear])
         {
            printf("error in binary searching for %e\n", x);
            printf("head:%d  value:%e\n", head, a[head]);
            printf("rear:%d  value:%e\n", rear, a[rear]);
            exit(0);
         }
         if(x == a[head])
            return head;
         else
            return rear;
      } 
      middle = floor((head + rear) / 2.0);
      if(x <= a[middle])
      {
         rear = middle;
      }
      else
      {
         head = middle;
      }
   }
}

/***************************************************
 A different context of binary searching. Here we have
 a monotonically increasing function F(such as CDF) with
 all natural integers as its domain.
 The returned value is the smallest non-negative integer i that
 satisfies F(i)>=f for a given f. It is useful for
 calculating quantiles for discrete random variables
 related to binomial or poisson.
***************************************************/
int binary_search_dyneval(double x,  double (*func)(int m, double *extra_par, int n_extra_par), 
                          double *extra_par, int n_extra_par, int lower, int upper, int *location)
{
   int i, j, loop;
   int num_search = 100;
   int lower_bound, upper_bound;
   int head, rear, middle;
   double func_lower, func_mid, func_upper;
   
   if(lower > upper)
   {
      printf("Upper bound should be larger than lower bound\n");
      return(1);
   }    
   lower_bound = lower;
   if(lower_bound < 0)
   {
       printf("The lower bound can not be negative\n");
       return(1);
   }
   if((func_lower = (*func)(lower_bound, extra_par, n_extra_par)) >= x)
   {
      (* location) = lower_bound;
      return(0);
   }          
   upper_bound = upper;
   loop = 0;
   while((func_upper = (*func)(upper_bound, extra_par, n_extra_par)) < x)
   {
      upper_bound *= 2;
      if((++loop) > num_search)
      {
         printf("The given value %f can not be reached by the given function\n", x);
         printf("upper bound=%d  func_upper=%e\n", upper_bound, func_upper);
         return(1);
      }
   }   
   
   if(func_lower > func_upper)
   {
      printf("The function is not monotonically increasing\n");
      return(1);
   }   
   
      
   head = lower_bound;
   rear = upper_bound;
   while(1)
   {
      
      //if(DEBUG == 1)
      //   printf("lower:%d(%25.20f)  middle:%d(%25.20f)  upper:%d(%25.20f)\n", 
      //          head, func_lower, middle, func_mid, rear, func_upper);
      
      if(rear <= head)
      {
         printf("Error in dynamic binary searching\n");
         printf("lower:%d(%25.20f)  middle:%d(%25.20f)  upper:%d(%25.20f)\n", 
                head, func_lower, middle, func_mid, rear, func_upper);
         return(1);
      }
      if(head == rear - 1)
      {
         if(x - func_lower <= 0.0)
            (*location) = head;
         else
            (*location) = rear;
         return(0);
      }
      middle = (int) (head + rear) / 2;
      func_mid = (*func)(middle, extra_par, n_extra_par);
      if(x - func_mid > 0.0)
      {
         head = middle;
         func_lower = func_mid;
      }
      else
      {
         rear = middle;
         func_upper = func_mid;
      }
   }
}

/************************************************************************
  This is Kolmogorov-Smirnov test from numerical recipe.
  IT is testing deviation from standard normal. It could be
  revised for other distributions.
 ************************************************************************/           
 
double probks(double alam)
{
   int j;
   double EPS1 = 0.001;
   double EPS2 = 1.0e-8;
   double a2, fac=2.0, sum=0.0, term, termbf=0.0, term_pos;
   
   a2 = -2.0 * alam * alam;
   for(j=1; j<=100; j++)
   {
      term = fac * exp(a2 * j * j);
      sum += term;
      term_pos = fabs(term);
      if(term_pos <= EPS1 * termbf || term_pos <= EPS2 * sum)
         return sum;
      fac = -fac;
      termbf = term_pos;
   }
   return 1.0;
}          


void ks_test(double *data, unsigned int n, double *d, double *prob)              
{
   double probks(double alam);
   /*void quick_sort(unsigned int direction, unsigned int n, double *arr, unsigned int *brr);*/
   void sort_vector(unsigned int n, double *input, unsigned int direction, double *output, unsigned int *position);
   unsigned int j;
   double dt, en, ff, fn, fo=0.0;
   double temp1, temp2;
   double *data_sorted;
   int *index;
   
   data_sorted = (double *) malloc((size_t) n*sizeof(double));
   index = (int *) malloc((size_t) n*sizeof(int));
   
   for(j=0; j<n; j++)
   {
      data_sorted[j] = data[j];
      index[j] = j;
   }   
   quick_sort(1, n, data_sorted, index);
   /*sort_vector(n, data, 1, data_sorted, index);*/
   en = n;
   *d = 0;
   for(j=1; j<=n; j++)
   {
      fn = (double) j/en;
      ff = pnorm(data_sorted[j-1], 0, 1);
      /*printf("%e   %e\n", data_sorted[j-1], ff);*/
      temp1 = fabs(fo-ff);
      temp2 = fabs(fn-ff);
      dt = (temp1 > temp2)? temp1:temp2;
      if(dt > *d) *d = dt;
      fo = fn;
   }
   en = sqrt(en);
   *prob = probks((en+0.12+0.11/en)*(*d));
   
   free(data_sorted);
   free(index);
}      
