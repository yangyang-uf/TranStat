#include <math.h>

#define SIGN(a, b) ( (b) >= 0.0? fabs(a): -fabs(a) )
#ifndef max
#define max(a,b) ((a) > (b)? (a):(b) )
#endif
#ifndef min
#define min(a,b) ((a) < (b)? (a):(b) )
#endif
#define SHIFT(a,b,c,d) (a)=(b); (b)=(c); (c)=(d);


void frprmn(double p[], int n, double ftol, int *iter, double *fret,
           double (*func)(double []), void (*dfunc)(double [], double []))
{
   void dlinmin(double p[], double xi[], int n, double *fret, double (*func)(double []),
               void (*dfunc)(double [], double []));
   int j, its;
   double gg, gam, fp, dgg;
   double *g, *h, *xi;
   int ITMAX = 500;
   double eps = 1.0e-10;

   g = (double *) malloc((size_t) n * sizeof(double));
   h = (double *) malloc((size_t) n * sizeof(double));
   xi = (double *) malloc((size_t) n * sizeof(double));

   fp = (*func)(p);
   (*dfunc)(p, xi); 
   for(j=0; j<n; j++)
   {
      g[j] = -xi[j];
      xi[j] = h[j] = g[j];
   }
   for(its=0; its<ITMAX; its++)
   {
      *iter = its + 1;
      dlinmin(p, xi, n, fret, func, dfunc);
      if(2.0*fabs(*fret-fp) <= ftol * (fabs(*fret) + fabs(fp) + eps))
      {
         goto end;
      }
      fp = *fret;
      (*dfunc)(p, xi);
      dgg = gg = 0.0;
      for(j=0; j<n; j++)
      {
         gg += g[j] * g[j];
         /*dgg += xi[j] * xi[j];*/
         dgg += (xi[j] + g[j]) * xi[j];
      }
      if(gg == 0.0)
      {
         goto end;
      }
      gam = dgg / gg;
      for(j=0; j<n; j++)
      {
         g[j] = -xi[j];
         xi[j] = h[j] = g[j] + gam * h[j];
      }
   }
   printf("too many iterations in frprmn\n");

end:         
   free(g);
   free(h);
   free(xi);
   return;
}

#define TOL 2.0e-4

int ncom;
double *pcom, *xicom;
double (*nrfunc)(double []);
void (*nrdfunc)(double [], double []);

void dlinmin(double p[], double xi[], int n, double *fret, double (*func)(double []),
            void (*dfunc)(double [], double []))
{
   double brent(double ax, double bx, double cx, double (*func)(double),
                double tol, double *xmin);
   double f1dim(double x);
   double df1dim(double x);
   void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
               double (*func)(double));
   int j, error = 0;
   double xx, xmin, fx, fb, fa, bx, ax;

   ncom = n;
   pcom = (double *) malloc((size_t) n * sizeof(double));
   xicom = (double *) malloc((size_t) n * sizeof(double));
   nrfunc = func;
   nrdfunc = dfunc;

   for(j=0; j<n; j++)
   {
      pcom[j] = p[j];
      xicom[j] = xi[j];
   }
   ax = 0.0;
   xx = 1.0;
   mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, f1dim);  
   printf("a=%f  x=%f  b=%f  fa=%e  fx=%e  fb=%e\n", ax, xx, bx, fa, fx, fb); 
      
   (*fret) = brent(ax, xx, bx, f1dim, TOL, &xmin);
   
   for(j=0; j<n; j++)
   {
      xi[j] *= xmin;
      p[j] += xi[j];
   }
   free(xicom);
   free(pcom);
   
}

double df1dim(double x)
{
   int j;
   double df1=0.0;
   double *xt, *df;

   xt = (double *) malloc((size_t) ncom * sizeof(double));
   df = (double *) malloc((size_t) ncom * sizeof(double));
   for(j=0; j<ncom; j++)
   {
      xt[j] = pcom[j] + x *xicom[j];
   }
   (*nrdfunc)(xt, df);
   for(j=0; j<ncom; j++)
   {
      df1 += df[j]*xicom[j];
   }
   free(df);
   free(xt);
   return(df1);
}


double f1dim(double x)
{
   int j;
   double f, *xt;

   xt = (double *) malloc((size_t) ncom * sizeof(double));
   for(j=0; j<ncom; j++)
   {
      xt[j] = pcom[j] + x *xicom[j];
   }
   f = (*nrfunc)(xt);
   free(xt);
   return(f);
}

double brent(double ax, double bx, double cx, double (*func)(double),
             double tol, double *xmin)
{
   int ITMAX = 500;
   double ZEPS = 1.0e-10;
   double CGOLD = 0.381966;
   int iter;
   double a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
   double e = 0.0;
   
   a = (ax < cx? ax:cx);
   b = (ax > cx? ax:cx);
   x = w = v = bx;
   fw = fv= fx = (*func)(x);
   
   for(iter=0; iter<ITMAX; iter++)
   {
      xm = 0.5 * (a + b);
      tol1 = tol * fabs(x) + ZEPS;
      tol2 = 2.0 * tol1;
      if(fabs(x-xm) <= (tol2-0.5*(b-a)))
      {
         *xmin = x;
         return(fx);
      }
      if(fabs(e) > tol1)
      {
         r = (x-w) * (fx - fv);
         q = (x-v) * (fx - fw);
         p = 2.0 * (q - r);
         if(q > 0.0)  p = -p;
         q = fabs(q);
         etemp = e;
         e = d;
         if(fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
         {
            d = CGOLD *(e = (x >= xm? a-x : b-x));
         }
         else
         {
            d = p/q;
            u = x+d;
            if(u-a < tol2 || b-u < tol2)
               d = SIGN(tol1, xm-x);
         }
      }
      else
      {
         d = CGOLD * (e = (x >= xm)? a-x : b-x);
      }
      if(fabs(d) >= tol1)
      {
         u = x + d;
         fu = (*func)(u);
      }
      else
      {
         u = x + SIGN(tol1, d);
         fu = (*func)(u);
      }   
      if(fu <= fx)
      {
         if(u >= x)  a = x;
         else  b = x;
         SHIFT(v, w, x, u);
         SHIFT(fv, fw, fx, fu);
      }
      else
      {
         if(u < x)  a = u;
         else  b = u;
         if(fu <= fw || w == x)
         {
            v = w;
            w = u;
            fv = fw;
            fw = fu;
         }   
         else if (fu < fv || v == x || v == w)
         {
            v = u;
            fv = fu;
         }   
      }
   }
   printf("too many iterations in dbrent\n");
   *xmin = x;
   return(fx);
}


void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
            double (*func)(double))
{
   double ulim, u, r, q, fu, dum;
   double GOLD = 1.618034;
   double GLIMIT=100.0;
   double tiny=1e-20;

   *fa = (*func)(*ax);
   *fb = (*func)(*bx);
   
   
   if(*fb > *fa)
   {
      SHIFT(dum, *ax, *bx, dum);
      SHIFT(dum, *fb, *fa, dum);
   }
   *cx = (*bx) + GOLD * (*bx - *ax);
   *fc = (*func)(*cx);
   while(*fb > *fc)
   {
      r = (*bx - *ax) * (*fb - *fc);
      q = (*bx - *cx) * (*fb - *fa);
      u = (*bx) - ((*bx - *cx) * q - (*bx - *ax) * r) / (2.0 * SIGN(max(fabs(q-r), tiny), q-r));
      ulim = (*bx) + GLIMIT * (*cx - *bx);
      if( (*bx - u) * (u - *cx) > 0.0)
      {
         fu = (*func)(u);
         if(fu < *fc)
         {
            *ax = (*bx);
            *bx = u;
            *fa = (*fb);
            *fb = fu;
            return;
         }
         else if(fu > *fb)
         {
            *cx = u;
            *fc = fu;
            return;
         }
         u = (*cx) + GOLD * (*cx - *bx);
         fu = (*func)(u);
      }
      else if ((*cx - u) * (u - ulim) > 0.0)
      {
         fu = (*func)(u);
         if(fu < *fc)
         {
            SHIFT(*bx, *cx, u, (*cx) + GOLD * (*cx - *bx));
            SHIFT(*fb, *fc, fu, (*func)(u));
         }
      }
      else if((u - ulim) * (ulim - *cx) >= 0.0)
      {
         u = ulim;
         fu = (*func)(u);
      }
      else
      {
         u = (*cx) + GOLD * (*cx - *bx);
         fu = (*func)(u);
      }
      SHIFT(*ax, *bx, *cx, u);
      SHIFT(*fa, *fb, *fc, fu);
   }
   
}


#define MOV3(a,b,c,d,e,f)  (a)=(d); (b)=(e); (c)=(f);
double dbrent(double ax, double bx, double cx, double (*func)(double), double (*dfunc)(double),
             double tol, double *xmin)
{
   int ITMAX = 500;
   double ZEPS = 1.0e-10;
   int iter, ok1, ok2;
   double a, b, d, d1, d2, du, dv, dw, dx, e=0.0;
   double fu, fv, fw, fx, olde, tol1, tol2, u, u1, u2, v, w, x, xm;

   a = (ax < cx? ax:cx);
   b = (ax > cx? ax:cx);
   x = w = v = bx;
   fw = fv= fx = (*func)(x);
   dw = dv= dx = (*dfunc)(x);
   
   for(iter=0; iter<ITMAX; iter++)
   {
      xm = 0.5 * (a + b);
      tol1 = tol * fabs(x) + ZEPS;
      tol2 = 2.0 * tol1;
      if(fabs(x-xm) <= (tol2-0.5*(b-a)))
      {
         *xmin = x;
         return(fx);
      }
      if(fabs(e) > tol1)
      {
         d1 = 2.0 * (b-a);
         d2 = d1;
         if(dw != dx) d1 = (w-x)*dx/(dx-dw);
         if(dv != dx) d2 = (v-x)*dx/(dx-dv);
         u1 = x + d1;
         u2 = x + d2;
         ok1 = (a-u1) * (u1-b) > 0.0 && dx*d1 <= 0.0;
         ok2 = (a-u2) * (u2-b) > 0.0 && dx*d2 <= 0.0;
         olde = e;
         e = d;
         if(ok1 || ok2)
         {
            if(ok1 && ok2)
               d = (fabs(d1) < fabs(d2) ? d1 : d2);
            else if (ok1)
               d = d1;
            else
               d = d2;
            if(fabs(d) <= fabs(0.5*olde))
            {
               u = x+d;
               if(u - a < tol2 || b - u < tol2)
                 d = SIGN(tol1, xm-x);
            }
            else
            {
               d = 0.5 * (e = (dx >= 0.0)? a-x : b-x);
            }
         }
         else
         {
            d = 0.5 * (e = (dx >= 0.0)? a-x : b-x);
         }
      }
      else
      {
         d = 0.5 * (e = (dx >= 0.0)? a-x : b-x);
      }
      if(fabs(d) >= tol1)
      {
         u = x + d;
         fu = (*func)(u);
      }
      else
      {
         u = x + SIGN(tol1, d);
         fu = (*func)(u);
         
         if(fu > fx)
         {
            *xmin = x;
            return(fx);
         }
      }
      du = (*dfunc)(u);
      if(fu <= fx)
      {
         if(u >= x)  a = x;
         else  b = x;
         MOV3(v, fv, dv, w, fw, dw);
         MOV3(w, fw, dw, x, fx, dx);
         MOV3(x, fx, dx, u, fu, du);
      }
      else
      {
         if(u < x)  a = u;
         else  b = u;
         if(fu <= fw || w == x)
         {
            MOV3(v, fv, dv, w, fw, dw);
            MOV3(w, fw, dw, u, fu, du);
         }
         else if (fu < fv || v == x || v == w)
            MOV3(v, fv, dv, u, fu, du);
      }
   }
   printf("too many iterations in dbrent\n");
   return(0);
}

#define SIMPLEX_MAX 5000
#define GET_PSUM \
   for(j=0; j<ndim; j++)\
   {\
      for(sum=0.0, i=0; i<mpts; i++)  sum += par.data[i][j];\
      psum[j] = sum;\
   }
#ifndef swap
#define swap(a, b) {swap = (a); (a) = (b); (b) = swap;}
#endif
int down_hill_simplex(double *inipar, double *lower, double *upper, double *step, double *value, int ndim, int minimization,
                       double ftol, double (*func)(double *, double *), double *extra_par, int *nfunc)
{
   double amotry(MATRIX *par, double *lower, double *upper, double *y, double *psum, int ndim, 
                 double (*func)(double *, double*), double *extra_par, double sign, int ihi, double fac);
   int i, ihi, ilo, inhi, j, mpts = ndim+1;
   int error= 0;
   double rtol, sum, swap, ysave, ytry, *psum, *y, sign;
   double tiny = 1e-20;
   MATRIX par;
   
   sign = (minimization == 1)? 1.0 : -1.0;
   initialize_matrix(&par);
   inflate_matrix(&par, mpts, ndim, 0);
   y = (double *) malloc((size_t) mpts * sizeof(double));
   
   for(j=0; j<ndim; j++)
   {
      par.data[0][j] = inipar[j];
      /* if initial value is on boundary, move it off the boundary, so that a correct simplex can be constructed*/
      if(lower != NULL)
      {
         if(par.data[0][j] <= lower[j])  par.data[0][j] = lower[j] + (upper[j] - lower[j]) / 1000.0;
      }   
      if(upper != NULL)
      {
         if(par.data[0][j] >= upper[j])  par.data[0][j] = upper[j] - (upper[j] - lower[j]) / 1000.0;
      }   
   }
   
   for(i=1; i<mpts; i++)
   {
      for(j=0; j<ndim; j++)  par.data[i][j] = par.data[0][j];
      par.data[i][i-1] += step[i-1];
      if(lower != NULL)
      {
         if(par.data[i][i-1] < lower[i-1])  par.data[i][i-1] = lower[i-1];
      }   
      if(upper != NULL)
      {
         if(par.data[i][i-1] > upper[i-1])  par.data[i][i-1] = upper[i-1];
      }   
   }
   for(i=0; i<mpts; i++)
      y[i] = sign * (*func)(par.data[i], extra_par);
   
   psum = (double *) malloc((size_t) ndim * sizeof(double));
   *nfunc = 0;
   GET_PSUM
   
   for(;;)
   {
      /* ilo will be the lowest point (best), ihi the highest (worst) and inhi the next highest*/
      ilo = 0;
      ihi = y[0] > y[1] ? (inhi=1, 0) : (inhi=0, 1);
      for(i=0; i<mpts; i++)
      {
         if(y[i] <= y[ilo])  ilo = i;
         if(y[i] > y[ihi])
         {
            inhi = ihi;
            ihi = i;
         }
         else if(y[i] > y[inhi] && i != ihi)  inhi = i;
      }
      rtol = 2.0 * fabs(y[ihi] - y[ilo])/(fabs(y[ihi]) + fabs(y[ilo]) + tiny);
      if(rtol < ftol)
      {
         //printf("rtol=%e  ftol=%e\n", rtol, ftol);
         //for(i=0; i<mpts; i++) printf("%e  ", y[i]);
         //printf("\n");
         //for(i=0; i<mpts; i++)
         //{
         //   for(j=0; j<ndim; j++) 
         //      printf("%e  ", par.data[i][j]);
         //   printf("\n");
         //}
         for(i=0; i<ndim; i++)  inipar[i] = par.data[ilo][i];
         (* value) = sign * y[ilo];
         goto end;
      }
      if(*nfunc > SIMPLEX_MAX)
      {
         /*printf("Max number of iterations exceeded\n");*/
         for(i=0; i<ndim; i++)  inipar[i] = par.data[ilo][i];
         (* value) = sign * y[ilo];
         error = 1;
         goto end;
      }
      *nfunc += 2;
      //Reflex
      ytry = amotry(&par, lower, upper, y, psum, ndim, func, extra_par, sign, ihi, -1.0);
      //extenstion
      if(ytry <= y[ilo])
         ytry = amotry(&par, lower, upper, y, psum, ndim, func, extra_par, sign, ihi, 2.0);
      else if (ytry >= y[inhi])
      {
         ysave = y[ihi];
         //contraction
         ytry = amotry(&par, lower, upper, y, psum, ndim, func, extra_par, sign, ihi, 0.5);
         if(ytry >= ysave)
         {
            for(i=0; i<mpts; i++)
            {
               if(i != ilo)
               {
                  for(j=0; j<ndim; j++)
                     par.data[i][j] = psum[j] = 0.5 * (par.data[i][j] + par.data[ilo][j]);
                  y[i] = sign * (*func)(psum, extra_par);
               }
            }
            *nfunc += ndim;
            GET_PSUM
         }
      }
      else
         --(*nfunc);
   }
end:  
   free(psum);
   free(y);
   deflate_matrix(&par);
   return(error);
}
                        
double amotry(MATRIX *par, double *lower, double *upper, double *y, double *psum, int ndim, 
              double (*func)(double *, double *), double *extra_par, double sign, int ihi, double fac)
{
   int j;
   double fac1, fac2, ytry, *ptry;
   ptry = (double *) malloc((size_t) ndim * sizeof(double));
   fac1 = (1 - fac) / ndim;
   fac2 = fac1 - fac;
   for(j=0; j<ndim; j++)
   {
      ptry[j] = psum[j] * fac1 - par->data[ihi][j] * fac2;
      if(lower != NULL)
      {
         if(ptry[j] < lower[j])  ptry[j] = lower[j];
      }   
      if(upper != NULL)
      {
         if(ptry[j] > upper[j])  ptry[j] = upper[j];
      }   
   }   
   ytry = sign * (*func)(ptry, extra_par);
   /* 
   printf("response=%f  par=", ytry);
   for(j=0; j<ndim; j++)
      printf("  %e", ptry[j]);
   printf("\n\n");
   */ 
   if(ytry < y[ihi])
   {
      y[ihi] = ytry;
      for(j=0; j<ndim; j++)
      {
         psum[j] += ptry[j] - par->data[ihi][j];
         par->data[ihi][j] = ptry[j];
      }
      
   }
   free(ptry);
   return(ytry);
}   
   
            
         
