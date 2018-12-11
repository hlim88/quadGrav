/*-------------------------------------------------------------------------
 *
 *  $Id: cook_bbh.cc,v 1.2 2008-07-08 14:32:32 dneilsen Exp $
 *
 *-------------------------------------------------------------------------*/

#include <iostream>
#include <cmath>
#include <vector>

// Definition of the Bin_BH class
#include "./include/PublicID.hpp"

#ifdef AIX
#define COOK_BBH cook_bbh
#else
#define COOK_BBH cook_bbh_
#endif

using std::vector;  
using std::cout;
using std::endl;

extern "C" void COOK_BBH (
                double *gxx,double *gxy,double *gxz,
                double *gyy,double *gyz,double *gzz,
                double *Kxx,double *Kxy,double *Kxz,
                double *Kyy,double *Kyz,double *Kzz,
                double *Alpha, 
                double *Betaux,double *Betauy, double *Betauz,
                double *x3d, double *y3d, double *z3d,
                int *nnx,int *nny,int *nnz);

/*--------------------------------------------------------------------------
 *
 *
 *
 *-------------------------------------------------------------------------*/
void COOK_BBH (double *gxx3d,double *gxy3d,double *gxz3d,
                double *gyy3d,double *gyz3d,double *gzz3d,
                double *Kxx3d,double *Kxy3d,double *Kxz3d,
                double *Kyy3d,double *Kyz3d,double *Kzz3d,
                double *Alpha3d, 
                double *Betaux3d,double *Betauy3d, double *Betauz3d,
                double *x3d, double *y3d, double *z3d,
                int *nnx,int *nny,int *nnz) {


  
  int    i;
  char * dfile  = "resu.d";

  int nx = *nnx;
  int ny = *nny;
  int nz = *nnz;
  int nd = nx*ny*nz;

  vector<double> x(nd),y(nd),z(nd);

  const double CUtoKM = 1.0;
  const double InvCUtoKM = 1.0/CUtoKM;

  // rescale the coordinates
  for (i = 0; i < nd; i++) {
    x[i] = x3d[i]*CUtoKM;
    y[i] = y3d[i]*CUtoKM;
    z[i] = z3d[i]*CUtoKM;
  }

  // variables returning the interpolated data 
  vector<double> gxx, gxy, gxz, gyy, gyz, gzz;
  vector<double> Kxx, Kxy, Kxz, Kyy, Kyz, Kzz;
  vector<double> Shiftx, Shifty, Shiftz, Lapse;

#if 0

  // Interpolate!
  InterpolateData(x,y,z,
		  gxx, gxy, gxz, gyy, gyz, gzz,
		  Kxx, Kxy, Kxz, Kyy, Kyz, Kzz,
		  Shiftx, Shifty, Shiftz, Lapse);
  
#endif


  for (i=0;i<nd;i++) {
    gxx3d[i] = gxx[i];
    gxy3d[i] = gxy[i];
    gxz3d[i] = gxz[i];
    gyy3d[i] = gyy[i];
    gyz3d[i] = gyz[i];
    gzz3d[i] = gzz[i];

    Kxx3d[i] = Kxx[i] * CUtoKM;
    Kxy3d[i] = Kxy[i] * CUtoKM;
    Kxz3d[i] = Kxz[i] * CUtoKM;
    Kyy3d[i] = Kyy[i] * CUtoKM;
    Kyz3d[i] = Kyz[i] * CUtoKM;
    Kzz3d[i] = Kzz[i] * CUtoKM;

    Alpha3d[i]  = Lapse[i];
    Betaux3d[i] = Shiftx[i];
    Betauy3d[i] = Shifty[i];
    Betauz3d[i] = Shiftz[i];
  }

  // unscale the coordinates
  for (i = 0; i < nd; i++) {
    x3d[i] *= InvCUtoKM;
    y3d[i] *= InvCUtoKM;
    z3d[i] *= InvCUtoKM;
  }
}


