/*  gamma.cpp -- computation of gamma function.
      Algorithms and coefficient values from "Computation of Special
      Functions", Zhang and Jin, John Wiley and Sons, 1996.

  (C) 2003, C. Bond. All rights reserved.

 Returns gamma function of argument 'x'.

 NOTE: Returns 1e308 if argument is a negative integer or 0,
      or if argument exceeds 171.
*/
#include <math.h>


#define TINY    1.0e-30
#define SQRT2PI 2.506628274631               /* sqrt(2 * pi) */

   
/* ============
   integer power
   ============*/
double ipow(double x, long n)
{
  long nplus;
  int sign;
  double xplus, logx, product;
  if(n == 0) return(1.0);
  if(x == 1.0) return(1.0);
  if(x == 0.0)
  {
    /* if negative n, indicate overflow*/
    if(n < 0) 
    { 
      printf("zero denominator in function ipow()\n");
      exit(0);
    }  
    else return(0.0);
  }    
  xplus = fabs(x);
  logx = log(xplus);
  nplus = fabs(n);
  product = exp(nplus*logx);  
  if(x < 0.0)
  { 
    sign = ((nplus%2)==0)?1:-1;
    product *= sign;  
  }   
  return (n<0)?(1.0/product):product;    
}

/*===================
  gamma function
  ===================*/
double gamma(double x)
{
    int i,k,m;
    double ga,gr,r,z;

    static double g[] = {
        1.0,
        0.5772156649015329,
       -0.6558780715202538,
       -0.420026350340952e-1,
        0.1665386113822915,
       -0.421977345555443e-1,
       -0.9621971527877e-2,
        0.7218943246663e-2,
       -0.11651675918591e-2,
       -0.2152416741149e-3,
        0.1280502823882e-3,
       -0.201348547807e-4,
       -0.12504934821e-5,
        0.1133027232e-5,
       -0.2056338417e-6,
        0.6116095e-8,
        0.50020075e-8,
       -0.11812746e-8,
        0.1043427e-9,
        0.77823e-11,
       -0.36968e-11,
        0.51e-12,
       -0.206e-13,
       -0.54e-14,
        0.14e-14};

    if (x > 171.0) return 1e308;    /* This value is an overflow flag.*/
    if (x == (int)x) {
        if (x > 0.0) {
            ga = 1.0;               /* use factorial*/
            for (i=2;i<x;i++) {
               ga *= i;
            }
         }
         else
            ga = 1e308;
     }
     else {
        if (fabs(x) > 1.0) {
            z = fabs(x);
            m = (int)z;
            r = 1.0;
            for (k=1;k<=m;k++) {
                r *= (z-k);
            }
            z -= m;
        }
        else
            z = x;
        gr = g[24];
        for (k=23;k>=0;k--) {
            gr = gr*z+g[k];
        }
        ga = 1.0/(gr*z);
        if (fabs(x) > 1.0) {
            ga *= r;
            if (x < 0.0) {
                ga = -M_PI/(x*ga*sin(M_PI*x));
            }
        }
    }
    return ga;
}

/*=====================
  log gamma function
  =====================*/
double lgamma(double x)
{
    double x0,x2,xp,gl,gl0;
    int n,k;
    static double a[] = {
        8.333333333333333e-02,
       -2.777777777777778e-03,
        7.936507936507937e-04,
       -5.952380952380952e-04,
        8.417508417508418e-04,
       -1.917526917526918e-03,
        6.410256410256410e-03,
       -2.955065359477124e-02,
        1.796443723688307e-01,
       -1.39243221690590};

    x0 = x;
    if (x <= 0.0) return 1e308;
    else if ((x == 1.0) || (x == 2.0)) return 0.0;
    else if (x <= 7.0) {
        n = (int)(7-x);
        x0 = x+n;
    }
    x2 = 1.0/(x0*x0);
    xp = 2.0*M_PI;
    gl0 = a[9];
    for (k=8;k>=0;k--) {
        gl0 = gl0*x2 + a[k];
    }
    gl = gl0/x0+0.5*log(xp)+(x0-0.5)*log(x0)-x0;
    if (x <= 7.0) {
        for (k=1;k<=n;k++) {
            gl -= log(x0-1.0);
            x0 -= 1.0;
        }
    }
    return gl;
}


/*=====================
  digamma function
  this algorithm comes from
  http://lib.stat.cmu.edu/apstat/
  =====================*/
double digamma(double x)
{
    double Y, R, DIGAMA;
    double S, C, S3, S4, S5, D1;

    S = 1E-06;
    C = 8.5;
    S3 = 8.333333333E-02;
    S4 = 8.3333333333E-03;
    S5 = 3.968253968E-03;
    D1 = -0.5772156649;
    DIGAMA = 0;
    Y = x;
    if(Y <= 0)
    {
      printf("parameter of digamma function must be positive!\n");
      exit(0);
    }
    if(Y <= S)
    {
	return D1 - 1.0 / Y;
    }
    while(Y < C)
    {
      DIGAMA = DIGAMA - 1.0/Y;
      Y = Y + 1.0;
    }
    R = 1.0 / Y;
    DIGAMA = DIGAMA + log(Y) - 0.5*R;
    R = R * R;
    DIGAMA = DIGAMA - R*(S3 - R*(S4 - R*S5));
    return DIGAMA;
}


/*=====================
  trigamma function
  calculates trigamma(x) = d**2(log(gamma(x))) / dx**2
  this algorithm comes from
  http://lib.stat.cmu.edu/apstat/
  algorithm as121   Appl. Statist. (1978) vol 27, no. 1
  =====================*/

 double trigamma(double x)
 {
    double a, b, b2, b4, b6,b8, y, z, trigam;

    a = 1.0e-6;
    b = 5.0;
    /*b2, b4, b6 and b8 are Bernoulli numbers*/

    /*b2 = 0.1666666667;
    b4 = -0.03333333333;
    b6 = 0.02380952381;
    b8 =  -0.03333333333; */
    b2 = 1.0/6.0;
    b4 = -1.0/30.0;
    b6 = 0.02380952381;
    b8 = -1.0/30.0;

    /*check for positive value of x */

    trigam = 0;
    if(x <= 0)
    {
      printf("parameter of trigamma function must be positive!\n");
      exit(0);
    }
    z = x;
    /*use small value approximation if x <= a */
    if(z <= a)
      return 1.0 / (z * z);
    /* increase argument to (x+i) >= b*/
    while(z < b)
    {
       trigam = trigam + 1.0 / (z * z);
       z = z + 1.0;
    }
    /*apply asymptotic formula if argument >= b*/

    y = 1.0 / (z * z);
    trigam = trigam + 0.5 * y +
           (1.0 + y * (b2 + y * (b4 + y * (b6 + y * b8)))) / z;
    return trigam;
 }


   double InGamma(double a, double x)
/* ========================================================================
 * Evaluates the incomplete gamma function.
 * NOTE: use a > 0.0 and x >= 0.0
 *
 * The algorithm used to evaluate the incomplete gamma function is based on
 * Algorithm AS 32, J. Applied Statistics, 1970, by G. P. Bhattacharjee.
 * See also equations 6.5.29 and 6.5.31 in the Handbook of Mathematical
 * Functions, Abramowitz and Stegum (editors).  The absolute error is less 
 * than 1e-10 for all non-negative values of x.
 * ========================================================================
 */
{ 
   double t, sum, term, factor, f, g, c[2], p[3], q[3];
   long   n;

   if (x > 0.0)
     factor = exp(-x + a * log(x) - lgamma(a));
   else
     factor = 0.0;
   if (x < a + 1.0) {                 /* evaluate as an infinite series - */
     t    = a;                        /* A & S equation 6.5.29            */
     term = 1.0 / a;
     sum  = term;
     while (term >= TINY * sum) {     /* sum until 'term' is small */
       t++;
       term *= x / t;
       sum  += term;
     } 
     return (factor * sum);
   }
   else {                             /* evaluate as a continued fraction - */
     p[0]  = 0.0;                     /* A & S eqn 6.5.31 with the extended */
     q[0]  = 1.0;                     /* pattern 2-a, 2, 3-a, 3, 4-a, 4,... */
     p[1]  = 1.0;                     /* - see also A & S sec 3.10, eqn (3) */
     q[1]  = x;
     f     = p[1] / q[1];
     n     = 0;
     do {                             /* recursively generate the continued */
       g  = f;                        /* fraction 'f' until two consecutive */
       n++;                           /* values are small                   */
       if ((n % 2) > 0) {
         c[0] = ((double) (n + 1) / 2) - a;
         c[1] = 1.0;
       }
       else {
         c[0] = (double) n / 2;
         c[1] = x;
       }
       p[2] = c[1] * p[1] + c[0] * p[0];
       q[2] = c[1] * q[1] + c[0] * q[0];
       if (q[2] != 0.0) {             /* rescale to avoid overflow */
         p[0] = p[1] / q[2];
         q[0] = q[1] / q[2];
         p[1] = p[2] / q[2];
         q[1] = 1.0;
         f    = p[1];
       }
       if(n > 1000)  break;
     } while ((fabs(f - g) >= TINY) || (q[1] != 1.0));
     return (1.0 - factor * f);
   }
}

   double Beta(double a, double b)
/* ======================================================================
 * LogBeta returns the natural log of the beta function.
 * NOTE: use a > 0.0 and b > 0.0
 *
 * The algorithm used to evaluate the natural log of the beta function is
 * based on a simple equation which relates the gamma and beta functions.
 *
 */
{
   return exp(lgamma(a) + lgamma(b) - lgamma(a + b));
}


   double lbeta(double a, double b)
/* ======================================================================
 * LogBeta returns the natural log of the beta function.
 * NOTE: use a > 0.0 and b > 0.0
 *
 * The algorithm used to evaluate the natural log of the beta function is
 * based on a simple equation which relates the gamma and beta functions.
 *
 */
{
   return (lgamma(a) + lgamma(b) - lgamma(a + b));
}

/* =======================================================================
 * Evaluates the incomplete beta function.
 * NOTE: use a > 0.0, b > 0.0 and 0.0 <= x <= 1.0
 *
 * The algorithm used to evaluate the incomplete beta function is based on
 * equation 26.5.8 in the Handbook of Mathematical Functions, Abramowitz
 * and Stegum (editors).  The absolute error is less than 1e-10 for all x
 * between 0 and 1.
 * =======================================================================
 */

double InBeta(double a, double b, double x)
{
   double t, factor, f, g, c, p[3], q[3];
   int    swap;
   long   n;

   if (x > (a + 1.0) / (a + b + 1.0)) { // to accelerate convergence   
     swap = 1;                          // complement x and swap a & b 
     x    = 1.0 - x;
     t    = a;
     a    = b;
     b    = t;
   }
   else                                 // do nothing
     swap = 0;
   if (x > 0)
     factor = exp(a * log(x) + b * log(1.0 - x) - lbeta(a,b)) / a;
   else
     factor = 0.0;
   p[0] = 0.0;
   q[0] = 1.0;
   p[1] = 1.0;
   q[1] = 1.0;
   f    = p[1] / q[1];
   n    = 0;
   do {                               // recursively generate the continued 
     g = f;                           // fraction 'f' until two consecutive 
     n++;                             // values are small                   
     if ((n % 2) > 0) {
       t = (double) (n - 1) / 2;
       c = -(a + t) * (a + b + t) * x / ((a + n - 1.0) * (a + n));
     }
     else {
       t = (double) n / 2;
       c = t * (b - t) * x / ((a + n - 1.0) * (a + n));
     }
     p[2] = p[1] + c * p[0];
     q[2] = q[1] + c * q[0];
     if (q[2] != 0.0) {                 // rescale to avoid overflow 
       p[0] = p[1] / q[2];
       q[0] = q[1] / q[2];
       p[1] = p[2] / q[2];
       q[1] = 1.0;
       f    = p[1];
     }
     if(n > 1000)  break;
   } while ((fabs(f - g) >= TINY) || (q[1] != 1.0));
   if (swap) 
     return (1.0 - factor * f);
   else
     return (factor * f);
}


/* =======================================================================
      double precision function betain(x, p, q, beta, ifault)
      implicit double precision (a-h, o-z)
      algorithm as 63  appl. statist. (1973), vol.22, no.3
      computes incomplete beta function ratio for arguments
      x between zero and one, p and q positive.
      log of complete beta function, beta, is assumed to be known
   ======================================================================= */
/*
double InBeta2(double p, double q, double x)
{
      int indx, ns, n;
      double psq, cx, xx, rx, pp, qq;
      double term, temp, ai, betain;
      
      if(x <= 0)  return(0.0);
      if(x >= 1)  return(1.0);
      betain=x;

//     change tail if necessary and determine s
      psq = p + q;
      cx = 1.0 - x;
      if(p < psq * x) 
      {
         xx = cx;
         cx = x;
         pp = q;
         qq = p;
         indx = 1;
      }
      else
      {
         xx = x;
         pp = p;
         qq = q;
         indx= 0;
      }   
      
      term = 1.0;
      ai = 1.0;
      betain = 1.0;
      ns = (int) (qq + cx * psq);

//     user soper's reduction formulae.

      rx = xx / cx;
      temp = qq - ai;
      if(ns == 0) rx = xx;
      n = 0;
      while(n < 1000)
      {
         term = term * temp * rx / (pp + ai);
         betain += term;
         temp = (term < 0)? (-term) : term;
         if(temp <= TINY && temp <= TINY * betain) break;
         ai += 1.0;
         ns --;
         if(ns >= 0)
         {
            temp = qq - ai;
            if(ns == 0) rx = xx;
         }
         else
         {   
            temp = psq;
            psq += 1.0;
         }  
         n ++;
      }

//     calculate result

      betain = betain * exp(pp * log(xx) + (qq - 1.0) * log(cx) - lbeta(p,q)) / pp;
      if(indx) betain = 1.0 - betain;
      return(betain);
}
*/

/*================================
  choose m items from n items
  ================================*/
double choose(double n, double m)
{
  double y;
  if(n<m)
  {
    printf("invalid input for choose function\n");
    exit(0);
  }
  if(n == m || m == 0.0) return(1.0);
  if(m == 1.0 || m == n - 1.0) return(n);
  if(n == 0.0) return(1.0);
  y = lgamma(n + 1.0) - lgamma(m + 1.0) - lgamma(n - m + 1.0);
  y = exp(y);
  return y ;
}  

/*================================
  log of choose(n, m)
  ================================*/
double lchoose(double n, double m)
{
  double y;
  if(n<m)
  {
    printf("invalid input for choose function\n");
    exit(0);
  }  
  if(n == m || m == 0.0) return(0.0);
  if(m == 1.0 || m == n - 1.0) return(log(n));
  if(n == 0.0) return(0.0);
  y = lgamma(n + 1.0) - lgamma(m + 1.0) - lgamma(n - m + 1.0);
  return y;
}

/*======================================
  finite harmonic sum, 1+1/2+1/3+...+1/n
  ======================================*/
double harmonic_sum(int n)
{
  int i;
  double s;
  if(n<=100)
  {
    for(i=1; i<=n; i++)  s += 1.0 / i; 
  }
  else if (n <= 500)
  {
     s = 2.45 + log((n + 1)/ 7.0) + 0.5 / 7; // this approximation works better for n <= 500. this is based on 1/n approximated by log(1+1/n) - 0.5 * (1/n) * (1/n)
  }
  else
     s = log(n) - 0.5772156649; // 0.5772 is the Euler constant
  return s;
}
