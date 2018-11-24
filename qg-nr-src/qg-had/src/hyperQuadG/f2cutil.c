/*-------------------------------------------------------------------------
 *
 *  $Id: f2cutil.c,v 1.1 2010-07-17 22:09:42 carlos Exp $
 *
 *-------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define LTRACE 0

/*-------------------------------------------------------------------------
 *
 *  MAP_FUNC_PTR_
 *
 *  This routine simply converts an F90 assumed shape pointer to a C 
 *  pointer, which will be stored in a Fortran integer array.  

 *  This solution is not elegant, calling this routine for every function,
 *  but was the first thing that I could get working.  I tried reading the
 *  gridfunction module directly, but apparently had problems because
 *  the pointers are for assumed shape arrays, which are not equivalent
 *  to C pointers. 
 *
 *-------------------------------------------------------------------------*/
#ifdef AIX
void map_func_ptr(double **p, double *f) {
#else
void map_func_ptr_(double **p, double *f) {
#endif
  *p = f;
/*
  if (LTRACE) printf(" map_func_ptr: f = %p, %d\n",f,f);
*/
}

/*-------------------------------------------------------------------------
 *
 *  
 *
 *-------------------------------------------------------------------------*/
#ifdef AIX
void modify_func(double **ff, int *nx, int *ny, int *nz) 
#else
void modify_func_(double **ff, int *nx, int *ny, int *nz) 
#endif
{
  int i,j,n;
  double *f;

  f = *ff;
  n = (*nx)*(*ny)*(*nz);
  printf("...modify_func_: f=%p, ff=%p\n",f,ff);
  for (i = 0; i < n; i++) {
    f[i] = -9.0;
  }
}

/*-------------------------------------------------------------------------
 *
 *  
 *
 *-------------------------------------------------------------------------*/
#ifdef AIX
void set_int(int *p) 
#else
void set_int_(int *p) 
#endif
{
  *p = 201;
}

/*-------------------------------------------------------------------------
 *
 *  imin
 *
 *-------------------------------------------------------------------------*/
int imin(int a, int b) {
  return (a < b)? a : b;
} 

/*-------------------------------------------------------------------------
 *
 *  imax
 *
 *-------------------------------------------------------------------------*/
int imax(int a, int b) {
  return (a > b)? a : b;
}

/*-------------------------------------------------------------------------
 *
 *  min -- min is included in newer C++ compilers.  For older compilers
 *         we define it here.
 *
 *-------------------------------------------------------------------------*/
double min( double a,  double b) {
  return (a < b)? a : b;
}

/*-------------------------------------------------------------------------
 *
 *  max -- max is included in newer C++ compilers.  For older compilers
 *         we define it here.
 *
 *-------------------------------------------------------------------------*/
double max( double a,  double b) {
  return (a > b)? a : b;
}

/*-------------------------------------------------------------------------
 *
 *  L2norm
 *
 *-------------------------------------------------------------------------*/
double L2Norm(const double * const f, const int n) {

  int    i;
  double sum = 0.0;

  for (i = 0; i < n; ++i) {
    sum += f[i]*f[i];
  }

  return sqrt(sum/n);

}

/*-------------------------------------------------------------------------
 *
 *  L2normInt
 *
 *-------------------------------------------------------------------------*/
double L2NormInt(const double * const f, const double dx, const int n) 
{

  int    i;
  double sum = 0.0;

  for (i = 0; i < n; ++i) {
    sum += f[i]*f[i]*dx*dx*dx;
  }

  return sqrt(sum);

}

/*-------------------------------------------------------------------------
 *
 *  AveArray
 *
 *-------------------------------------------------------------------------*/
double AveArray(const double * const f, const int n) {

  int    i;
  double sum = 0.0;

  for (i = 0; i < n; ++i) {
    sum += f[i];
  }

  return sum/n;

}


/*-------------------------------------------------------------------------
 *
 *  SetToZero
 *
 *-------------------------------------------------------------------------*/
void SetToZero(double *f, int n) 
{
  int i;
  for(i = 0; i < n; i++) f[i] = 0.0;
} 

/*-------------------------------------------------------------------------
 *
 *  for_isfinite
 *
 *  To be called from fortran, which does not implement isfinite in
 *  a general way, yet.  Fortran 2003 will include this.
 *
 *-------------------------------------------------------------------------*/
#ifdef AIX
void    for_isfinite( double *x, int *rc)
#else
void    for_isfinite_(double *x, int *rc)
#endif  
{

  *rc = isfinite(*x);

}


/*-------------------------------------------------------------------------
 *
 *
 *
 *-------------------------------------------------------------------------*/
#ifdef AIX
void check_isfinite_1d(double *f, int *rc, int *n) {
#else
void check_isfinite_1d_(double *f, int *rc, int *n) {
#endif

  int  i, ret;
  ret = 1;

  for (i = 0; i < *n; i++) {
    ret *= isfinite(f[i]);
  }

  *rc = ret;
}
