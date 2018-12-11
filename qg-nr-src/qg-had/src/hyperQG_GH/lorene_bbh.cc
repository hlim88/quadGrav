/*-------------------------------------------------------------------------
 *
 *  $Id: lorene_bbh.cc,v 1.6 2009-10-05 21:54:07 matt Exp $
 *
 *-------------------------------------------------------------------------*/

#include <iostream>
#include <cmath>

// Definition of the Bin_BH class
#include "bin_bh.h"

#ifdef AIX
#define LORENE_BBH lorene_bbh
#else
#define LORENE_BBH lorene_bbh_
#endif

using namespace std;

extern "C" void LORENE_BBH (
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
void LORENE_BBH (double *gxx,double *gxy,double *gxz,
                double *gyy,double *gyz,double *gzz,
                double *Kxx,double *Kxy,double *Kxz,
                double *Kyy,double *Kyz,double *Kzz,
                double *Alpha, 
                double *Betaux,double *Betauy, double *Betauz,
                double *x3d, double *y3d, double *z3d,
                int *nnx,int *nny,int *nnz) {

  int    i;
  char * dfile  = "resu.d";

  int nx = *nnx;
  int ny = *nny;
  int nz = *nnz;
  int nd = nx*ny*nz;

  int fill = 1;   // fill = 0 sets the fields to zero inside the mask
                  // fill = 1 extrapolates the fields into the mask.

  const double CUtoKM = 1.0;
  const double InvCUtoKM = 1.0/CUtoKM;

  // rescale the coordinates
  for (i = 0; i < nd; i++) {
    x3d[i] *= CUtoKM;
    y3d[i] *= CUtoKM;
    z3d[i] *= CUtoKM;
  }

  Bin_BH binary(nd, x3d, y3d, z3d, fill, dfile);

  for (i=0;i<nd;i++) {
    gxx[i] = binary.g_xx[i];
    gxy[i] = binary.g_xy[i];
    gxz[i] = binary.g_xz[i];
    gyy[i] = binary.g_yy[i];
    gyz[i] = binary.g_yz[i];
    gzz[i] = binary.g_zz[i];

    Kxx[i] = binary.k_xx[i] * CUtoKM;
    Kxy[i] = binary.k_xy[i] * CUtoKM;
    Kxz[i] = binary.k_xz[i] * CUtoKM;
    Kyy[i] = binary.k_yy[i] * CUtoKM;
    Kyz[i] = binary.k_yz[i] * CUtoKM;
    Kzz[i] = binary.k_zz[i] * CUtoKM;

    Alpha[i]  =  binary.nnn[i];
    Betaux[i] = -binary.beta_x[i];
    Betauy[i] = -binary.beta_y[i];
    Betauz[i] = -binary.beta_z[i];
  }

  // unscale the coordinates
  for (i = 0; i < nd; i++) {
    x3d[i] *= InvCUtoKM;
    y3d[i] *= InvCUtoKM;
    z3d[i] *= InvCUtoKM;
  }
}


