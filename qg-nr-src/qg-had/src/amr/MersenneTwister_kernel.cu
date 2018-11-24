/*
 * Copyright 1993-2007 NVIDIA Corporation.  All rights reserved.
 *
 * NOTICE TO USER:
 *
 * This source code is subject to NVIDIA ownership rights under U.S. and
 * international Copyright laws.  Users and possessors of this source code
 * are hereby granted a nonexclusive, royalty-free license to use this code
 * in individual and commercial software.
 *
 * NVIDIA MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE
 * CODE FOR ANY PURPOSE.  IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR
 * IMPLIED WARRANTY OF ANY KIND.  NVIDIA DISCLAIMS ALL WARRANTIES WITH
 * REGARD TO THIS SOURCE CODE, INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.
 * IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL,
 * OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS,  WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION,  ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOURCE CODE.
 *
 * U.S. Government End Users.   This source code is a "commercial item" as
 * that term is defined at  48 C.F.R. 2.101 (OCT 1995), consisting  of
 * "commercial computer  software"  and "commercial computer software
 * documentation" as such terms are  used in 48 C.F.R. 12.212 (SEPT 1995)
 * and is provided to the U.S. Government only as a commercial end item.
 * Consistent with 48 C.F.R.12.212 and 48 C.F.R. 227.7202-1 through
 * 227.7202-4 (JUNE 1995), all U.S. Government End Users acquire the
 * source code with only those rights set forth herein.
 *
 * Any use of this source code in individual and commercial software must
 * include, in the user documentation and internal comments to the code,
 * the above Disclaimer and U.S. Government End Users Notice.
 */

/* 15 Apr 2009 
 * Matt Anderson
 * Code to do neutron scattering for tests
 */

#include "MersenneTwister.h"



__device__ static mt_struct_stripped ds_MT[MT_RNG_COUNT];
static mt_struct_stripped h_MT[MT_RNG_COUNT];



void initMTGPU(const char *fname){
    FILE *fd = fopen(fname, "rb");
    if(!fd){
        printf("initMTGPU(): failed to open %s\n", fname);
        printf("TEST FAILED\n");
        exit(0);
    }
    if( !fread(h_MT, sizeof(h_MT), 1, fd) ){
        printf("initMTGPU(): failed to load %s\n", fname);
        printf("TEST FAILED\n");
        exit(0);
    }
    fclose(fd);

    CUDA_SAFE_CALL( cudaMemcpyToSymbol(ds_MT, h_MT, sizeof(h_MT)) );
}



////////////////////////////////////////////////////////////////////////////////
// Write MT_RNG_COUNT vertical lanes of NPerRng random numbers to *d_Random.
// For coalesced global writes MT_RNG_COUNT should be a multiple of warp size.
// Initial states for each generator are the same, since the states are
// initialized from the global seed. In order to improve distribution properties
// on small NPerRng supply dedicated (local) seed to each twister.
// The local seeds, in their turn, can be extracted from global seed
// by means of any simple random number generator, like LCG.
////////////////////////////////////////////////////////////////////////////////
__global__ void RandomGPU(
    float *d_Random,
    int NPerRng,
    unsigned int seed
){
    const int      tid = blockDim.x * blockIdx.x + threadIdx.x;
    const int THREAD_N = blockDim.x * gridDim.x;

    int iState, iState1, iStateM, iOut;
    unsigned int mti, mti1, mtiM, x;
    unsigned int mt[MT_NN];

    for(int iRng = tid; iRng < MT_RNG_COUNT; iRng += THREAD_N){
        //Load bit-vector Mersenne Twister parameters
        mt_struct_stripped config = ds_MT[iRng];

        //Initialize current state
        mt[0] = seed;
        for(iState = 1; iState < MT_NN; iState++)
            mt[iState] = (1812433253U * (mt[iState - 1] ^ (mt[iState - 1] >> 30)) + iState) & MT_WMASK;

        iState = 0;
        mti1 = mt[0];
        for(iOut = 0; iOut < NPerRng; iOut++){
            //iState1 = (iState +     1) % MT_NN
            //iStateM = (iState + MT_MM) % MT_NN
            iState1 = iState + 1;
            iStateM = iState + MT_MM;
            if(iState1 >= MT_NN) iState1 -= MT_NN;
            if(iStateM >= MT_NN) iStateM -= MT_NN;
            mti  = mti1;
            mti1 = mt[iState1];
            mtiM = mt[iStateM];

            x    = (mti & MT_UMASK) | (mti1 & MT_LMASK);
            x    =  mtiM ^ (x >> 1) ^ ((x & 1) ? config.matrix_a : 0);
            mt[iState] = x;
            iState = iState1;

            //Tempering transformation
            x ^= (x >> MT_SHIFT0);
            x ^= (x << MT_SHIFTB) & config.mask_b;
            x ^= (x << MT_SHIFTC) & config.mask_c;
            x ^= (x >> MT_SHIFT1);

            //Convert to (0, 1] float and write to global memory
            d_Random[iRng + iOut * MT_RNG_COUNT] = ((float)x + 1.0f) / 4294967296.0f;
        }
    }
}



////////////////////////////////////////////////////////////////////////////////
// Transform each of MT_RNG_COUNT lanes of NPerRng uniformly distributed 
// random samples, produced by RandomGPU(), to normally distributed lanes
// using Cartesian form of Box-Muller transformation.
// NPerRng must be even.
////////////////////////////////////////////////////////////////////////////////
#define PI 3.14159265358979323846264338327950288f

__device__ void BoxMuller(float& u1, float& u2,
                          float& u3, float& u4,
                          float& u5, float& u6,
                          float& u7, float& u8,
                          float& u9, float& u10,
                          float& u11, float& u12,
                          float E0, float Ef,
                          float hdens, float odens,
                          float radius){
    float u,v,w,x,y,z;
    float xtest,ytest,E;
    float hxs,oxs,xs,d,r2,r;
    float costh,sinth,phi;
    float ul,vl,wl,a;
    float sr,ux,vx,wx,vec;
    int i,ie,j,step;
    float myrand[12];

    myrand[0] = u1;
    myrand[1] = u2;
    myrand[2] = u3;
    myrand[3] = u4;
    myrand[4] = u5;
    myrand[5] = u6;
    myrand[6] = u7;
    myrand[7] = u8;
    myrand[8] = u9;
    myrand[9] = u10;
    myrand[10] = u11;
    myrand[11] = u12;

    // hydrogen and oxygen data {{{
    int found = 0;
    float h_xspoint[11];
    float h_epoint[11];
    float h_slope;

    float o_xspoint[84];
    float o_epoint[84];
    float o_slope;

    h_epoint[0] = 0.01;
    h_epoint[1] = 0.1;
    h_epoint[2] = 1.0;
    h_epoint[3] = 1000.0;
    h_epoint[4] = 10000.0;
    h_epoint[5] = 50000.0;
    h_epoint[6] = 1.e5;
    h_epoint[7] = 1.e6;
    h_epoint[8] = 5.e6;
    h_epoint[9] = 1.e7;
    h_epoint[10] = 3.e7;

    h_xspoint[0] = 68.0;
    h_xspoint[1] = 26.0;
    h_xspoint[2] = 20.0;
    h_xspoint[3] = 20.0;
    h_xspoint[4] = 19.0;
    h_xspoint[5] = 15.0;
    h_xspoint[6] = 7.0;
    h_xspoint[7] = 4.0;
    h_xspoint[8] = 1.5;
    h_xspoint[9] = 0.9;
    h_xspoint[10] = 0.3;

    o_xspoint[0] = 3.7;
    o_xspoint[1] = 3.7;
    o_xspoint[2] = 3.7;
    o_xspoint[3] = 3.5;
    o_xspoint[4] = 3.8;
    o_xspoint[5] = 16.0;
    o_xspoint[6] = 4.0;
    o_xspoint[7] = 3.1;
    o_xspoint[8] = 2.7;

    o_xspoint[9] = 8.0;
    o_xspoint[10] = 2.7;
    o_xspoint[11] = 7.0;
    o_xspoint[12] = 2.5;
    o_xspoint[13] = 1.9;
    o_xspoint[14] = 9.0;
    o_xspoint[15] = 1.9;
    o_xspoint[16] = 1.9;
    o_xspoint[17] = 5.0;
    o_xspoint[18] = 1.8;
    o_xspoint[19] = 1.6;
    o_xspoint[20] = 5.0;
    o_xspoint[21] = 2.0;
    o_xspoint[22] = 3.0;
    o_xspoint[23] = 1.5;

    o_xspoint[24] = 1.0;
    o_xspoint[25] = 0.09;
    o_xspoint[26] = 0.7;
    o_xspoint[27] = 1.1;
    o_xspoint[28] = 1.1;
    o_xspoint[29] = 1.4;
    o_xspoint[30] = 5.0;
    o_xspoint[31] = 3.0;
    o_xspoint[32] = 3.0;
    o_xspoint[33] = 5.5;
    o_xspoint[34] = 3.0;
    o_xspoint[35] = 3.0;
    o_xspoint[36] = 5.6;
    o_xspoint[37] = 3.0;
    o_xspoint[38] = 2.0;
    o_xspoint[39] = 1.1;
    o_xspoint[40] = 2.1;

    o_xspoint[41] = 1.4;
    o_xspoint[42] = 2.0;
    o_xspoint[43] = 3.0;
    o_xspoint[44] = 1.3;
    o_xspoint[45] = 2.0;
    o_xspoint[46] = 0.9;
    o_xspoint[47] = 2.0;
    o_xspoint[48] = 0.9;
    o_xspoint[49] = 3.95;
    o_xspoint[50] = 1.4;
    o_xspoint[51] = 1.2;
    o_xspoint[52] = 2.1;
    o_xspoint[53] = 0.8;
    o_xspoint[54] = 1.05;
    o_xspoint[55] = 1.2;
    o_xspoint[56] = 2.6;
    o_xspoint[57] = 1.7;

    o_xspoint[58] = 1.4;
    o_xspoint[59] = 0.8;
    o_xspoint[60] = 1.6;
    o_xspoint[61] = 1.0;
    o_xspoint[62] = 0.8;
    o_xspoint[63] = 1.7;
    o_xspoint[64] = 0.6;
    o_xspoint[65] = 1.1;
    o_xspoint[66] = 1.1;
    o_xspoint[67] = 1.5;
    o_xspoint[68] = 0.85;
    o_xspoint[69] = 0.85;
    o_xspoint[70] = 1.9;
    o_xspoint[71] = 1.4;
    o_xspoint[72] = 1.1;

    o_xspoint[73] = 1.5;
    o_xspoint[74] = 0.8;
    o_xspoint[75] = 1.4;
    o_xspoint[76] = 1.1;
    o_xspoint[77] = 1.5;
    o_xspoint[78] = 1.1;
    o_xspoint[79] = 1.35;
    o_xspoint[80] = 1.08;
    o_xspoint[81] = 1.4;
    o_xspoint[82] = 1.2;
    o_xspoint[83] = 1.6;

    o_epoint[0] = 0.01;
    o_epoint[1] = 0.1;
    o_epoint[2] = 1.e4;
    o_epoint[3] = 3.e5;
    o_epoint[4] = 3.5e5;
    o_epoint[5] = 4.4e5;
    o_epoint[6] = 5.2e5;
    o_epoint[7] = 6.e5;

    o_epoint[8] = 8.5e5;
    o_epoint[9] = 1.e6;
    o_epoint[10] = 1.2e6;
    o_epoint[11] = 1.31e6;
    o_epoint[12] = 1.4e6;
    o_epoint[13] = 1.64e6;
    o_epoint[14] = 1.65e6;

    o_epoint[15] = 1.66e6;
    o_epoint[16] = 1.67e6;
    o_epoint[17] = 1.68e6;
    o_epoint[18] = 1.69e6;

    o_epoint[19] = 1.8e6;
    o_epoint[20] = 1.82e6;
    o_epoint[21] = 1.83e6;
    o_epoint[22] = 1.9e6;
    o_epoint[23] = 2.0e6;
    o_epoint[24] = 2.25e6;
    o_epoint[25] = 2.35e6;

    o_epoint[26] = 2.4e6;
    o_epoint[27] = 2.5e6;
    o_epoint[28] = 2.9e6;
    o_epoint[29] = 3.1e6;

    o_epoint[30] = 3.2e6;
    o_epoint[31] = 3.21e6;
    o_epoint[32] = 3.4e6;
    o_epoint[33] = 3.41e6;
    o_epoint[34] = 3.42e6;
    o_epoint[35] = 3.78e6;
    o_epoint[36] = 3.79e6;
    o_epoint[37] = 3.8e6;

    o_epoint[38] = 4.e6;
    o_epoint[39] = 4.1e6;
    o_epoint[40] = 4.4e6;
    o_epoint[41] = 4.45e6;
    o_epoint[42] = 4.5e6;
    o_epoint[43] = 4.51e6;

    o_epoint[44] = 4.52e6;
    o_epoint[45] = 4.6e6;
    o_epoint[46] = 4.7e6;
    o_epoint[47] = 4.81e6;
    o_epoint[48] = 5.e6;
    o_epoint[49] = 5.15e6;
    o_epoint[50] = 5.2e6;

    o_epoint[51] = 5.35e6;
    o_epoint[52] = 5.36e6;
    o_epoint[53] = 5.37e6;
    o_epoint[54] = 5.4e6;
    o_epoint[55] = 5.65e6;
    o_epoint[56] = 5.66e6;
    o_epoint[57] = 5.67e6;

    o_epoint[58] = 5.98e6;
    o_epoint[59] = 5.99e6;

    o_epoint[60] = 6.05e6;
    o_epoint[61] = 6.06e6;
    o_epoint[62] = 6.3e6;
    o_epoint[63] = 6.4e6;
    o_epoint[64] = 6.5e6;
    o_epoint[65] = 6.8e6;
    o_epoint[66] = 6.85e6;

    o_epoint[67] = 6.86e6;
    o_epoint[68] = 6.88e6;
    o_epoint[69] = 7.e6;
    o_epoint[70] = 7.2e6;
    o_epoint[71] = 7.4e6;
    o_epoint[72] = 7.5e6;
    o_epoint[73] = 7.8e6;

    o_epoint[74] = 8.05e6;
    o_epoint[75] = 8.3e6;
    o_epoint[76] = 8.6e6;
    o_epoint[77] = 8.7e6;
    o_epoint[78] = 9.e6;
    o_epoint[79] = 9.1e6;
    o_epoint[80] = 9.3e6;

    o_epoint[81] = 11.e6;
    o_epoint[82] = 12.e6;
    o_epoint[83] = 18.e6;
    // }}}

    u = 0.0;
    v = 0.0;
    w = 1.0;

    // source location
    x = 0.0;
    y = 0.0;
    z = 0.0;

    step = 0;
    for (j=0;j<11;j++) {
      switch(step)
      {
        case 0: 
          xtest = -logf(myrand[j]);  
          step = 1;
          break;
        case 1: 
          ytest = -logf(myrand[j]);  

          if ( powf(ytest-xtest-1.0,2) <= 4.*xtest ) {
            E = 2.*xtest*1.e6; // we are working in eV
            step = 2;
          } else {
            step = 0;
          } 
          break;
        case 2: 
          // hydrogen {{{
          found = 0;
          for (i=1;i<11;i++) {
            ie = i;
            if ( E < h_epoint[ie] ) {
              found = 1;
              break;
            }
          }
          if ( found == 0 ) {
            // E >= 3E7
            hxs = hdens*h_xspoint[ie];
          } else {
            h_slope = (h_xspoint[ie] - h_xspoint[ie-1])/logf(h_epoint[ie]/h_epoint[ie-1]);
            hxs = h_slope*logf(E/h_epoint[ie-1])+h_xspoint[ie-1]; // extrapolate between ie-1 and ie
            hxs = hxs*hdens;
          }
          // }}}

          // oxygen {{{
          found = 0;
          for (i=1;i<84;i++) {
            ie = i;
            if ( E < o_epoint[ie] ) {
              found = 1;
              break;
            }
          }
          if ( found == 0 ) {
            // E >= 1.8E7
            oxs = odens*o_xspoint[ie];
          } else {
            o_slope = (o_xspoint[ie] - o_xspoint[ie-1])/logf(o_epoint[ie]/o_epoint[ie-1]);
            oxs = o_slope*logf(E/o_epoint[ie-1])+o_xspoint[ie-1]; // extrapolate between ie-1 and ie
            oxs = oxs*odens;
          }
          // }}}

          // macro cross section xs in 1/cm
          xs = hxs + oxs;

          // get a flight path d in cm
          d = -logf(myrand[j])/xs;

          // collision site
          x = x + d*u;
          y = y + d*v;
          z = z + d*w;

          // radius of collision site
          r2 = x*x + y*y + z*z;
          r = sqrtf(r2);

          if ( r >= radius ) {
            // the particle escaped
            step = 0;

            u = 0.0;
            v = 0.0;
            w = 1.0;

            // source location
            x = 0.0;
            y = 0.0;
            z = 0.0;
          } else {
            step = 3;
          } 
          break;
        case 3: 
          // random number to select type of atom
          if (myrand[j] < hxs/xs) {
            // if collision in hydrogen
            a = 1.0;
          } else {
            // if collision in oxygen
            a = 16.0;
          }
          step = 4;
          break;
        case 4: 
          costh = 2.*myrand[j] - 1.0;
          step = 5;
          break;
        case 5: 
          phi = 2.*myrand[j]*PI;
          step = 6;

          wl = (a*costh+1.0)/sqrtf(a*a+2.*a*costh+1.0);
          sinth = sqrtf(1.-wl*wl);

          ul = sinth*cosf(phi);
          vl = sinth*sinf(phi);
 
          // post-collision energy
          E = E*(a*a + 2.*a*costh+1.0)/powf(a+1.0,2);

          if ( fabs(u) < 0.9 ) {
            // x-axis transformation Eqn. 4.43
            sr = sqrtf(1.0-u*u);
            ux = sr*ul + u*wl;
            vx = -u*v*ul/sr + w*vl/sr + v*wl;
            wx = -w*u*ul/sr - v*vl/sr + w*wl;
          } else {
            // y-axis transformation Eqn. 4.44
            sr = sqrt(1.0-v*v);
            ux = w*ul/sr - u*v*vl/sr + u*wl;
            vx = vl*sr + v*wl;
            wx = -u*ul/sr - v*w*vl/sr + w*wl;
          }
          vec = sqrtf(ux*ux + vx*vx + wx*wx);
               
          // normalized lab direction cosines
          u = ux/vec;
          v = vx/vec;
          w = wx/vec;

          if ( E > Ef ) {
            step = 2;
          } else {
            step = 6;
          }
          break;
      }
      if ( step == 6 ) break;
    }

    if ( step == 6 ) {
      u1 = r2;
    } else {
      u1 = -1.0;
    }
    u2 = 15.0;
    u3 = -1.0;
    u4 = -1.0;
    u5 = -1.0;
    u6 = -1.0;
    u7 = -1.0;
    u8 = -1.0;
    u9 = -1.0;
    u10 = -1.0;
    u11 = -1.0;
    u12 = -1.0;
}

__global__ void BoxMullerGPU(float *d_Random, int NPerRng,
                             float E0,float Ef,float hdens,
                             float odens, float radius){
    const int      tid = blockDim.x * blockIdx.x + threadIdx.x;
    const int THREAD_N = blockDim.x * gridDim.x;

    for(int iRng = tid; iRng < MT_RNG_COUNT; iRng += THREAD_N)
        for(int iOut = 0; iOut < NPerRng; iOut += 12)
            BoxMuller(
                d_Random[iRng + (iOut + 0) * MT_RNG_COUNT],
                d_Random[iRng + (iOut + 1) * MT_RNG_COUNT],
                d_Random[iRng + (iOut + 2) * MT_RNG_COUNT],
                d_Random[iRng + (iOut + 3) * MT_RNG_COUNT],
                d_Random[iRng + (iOut + 4) * MT_RNG_COUNT],
                d_Random[iRng + (iOut + 5) * MT_RNG_COUNT],
                d_Random[iRng + (iOut + 6) * MT_RNG_COUNT],
                d_Random[iRng + (iOut + 7) * MT_RNG_COUNT],
                d_Random[iRng + (iOut + 8) * MT_RNG_COUNT],
                d_Random[iRng + (iOut + 9) * MT_RNG_COUNT],
                d_Random[iRng + (iOut + 10) * MT_RNG_COUNT],
                d_Random[iRng + (iOut + 11) * MT_RNG_COUNT],
                E0,Ef,hdens,odens,radius
            );
}
