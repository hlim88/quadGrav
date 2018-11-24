/* Matt Anderson
 * 21 Mar 2009
 * Neutrino transport via Monte Carlo on the GPU
 */



#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <cutil.h>
#include "MersenneTwister.h"

////////////////////////////////////////////////////////////////////////////////
// Common host and device functions
////////////////////////////////////////////////////////////////////////////////
//Round a / b to nearest higher integer value
int iDivUp(int a, int b){
    return (a % b != 0) ? (a / b + 1) : (a / b);
}

//Round a / b to nearest lower integer value
int iDivDown(int a, int b){
    return a / b;
}

//Align a to nearest higher multiple of b
int iAlignUp(int a, int b){
    return (a % b != 0) ?  (a - a % b + b) : a;
}

//Align a to nearest lower multiple of b
int iAlignDown(int a, int b){
    return a - a % b;
}

float RandFloat(float low, float high){
    float t = (float)rand() / (float)RAND_MAX;
    return (1.0f - t) * low + t * high;
}


///////////////////////////////////////////////////////////////////////////////
// GPU code
///////////////////////////////////////////////////////////////////////////////
#include "MersenneTwister_kernel.cu"
#include "Neutrino_kernel.cu"



///////////////////////////////////////////////////////////////////////////////
// Data configuration
///////////////////////////////////////////////////////////////////////////////
//Simulation paths (random samples) count 
//const int PATH_N = 24000000;
//Number of outputs per generator; align to even for Box-Muller transform
//const int N_PER_RNG = iAlignUp(iDivUp(PATH_N, MT_RNG_COUNT), 2);
//Total numbers of sample to generate
//const int    RAND_N = MT_RNG_COUNT * N_PER_RNG;

//Reduce problem size to have reasonable emulation time
const int  OPT_N = 1;
//#ifndef __DEVICE_EMULATION__
//const int  OPT_N = 128;
//#else
//const int  OPT_N = 4;
//#endif


///////////////////////////////////////////////////////////////////////////////
// Start off point for neutrino transport calculation
///////////////////////////////////////////////////////////////////////////////
extern "C" void neutrino_(double *etemp,double *edens,double *xx,double *yy, double *zz)
{
    float
        *d_Random;

    float
        h_particlesGPU[OPT_N],
        h_r2barsqGPU[OPT_N],
        x[OPT_N];

    double gpuTime;

    //Simulation paths (random samples) count 
    int PATH_N = 10000;
    //Number of outputs per generator; align to even for Box-Muller transform
    int N_PER_RNG = iAlignUp(iDivUp(PATH_N, MT_RNG_COUNT), 12);
    //Total numbers of sample to generate
    int RAND_N = MT_RNG_COUNT * N_PER_RNG;
    
    // the start particle energy -- eV
    double E0 = 0.0;

    // the cutoff particle energy -- eV
    double Ef = 1.4;

    // atomic density of hydrogen and oxygen: units g/cm^3
    double hdens, odens;
    hdens = 0.06692;
    odens = 0.03346;

    // radius of fluid
    double radius = 300.0;

    // random seed
    int iseed = 5;

    unsigned int hTimer;

    CUT_DEVICE_INIT();
    CUT_SAFE_CALL( cutCreateTimer(&hTimer) );

    printf("Loading GPU twisters configurations...\n");
    
    //const char *dat_path = cutFindFilePath("MersenneTwister.dat", argv[0]);
    //  initMTGPU(dat_path);
      initMTGPU("/home/matt/had_cvs/had/src/hyperGHMHD/MersenneTwister.dat");

    printf("Generating random options...\n");
    
    CUDA_SAFE_CALL( cudaMalloc((void **)&d_Random, 12*RAND_N  * sizeof(float)) );

    srand(iseed);

    printf("Data init done.\n");

    printf("RandomGPU()...\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
    CUT_SAFE_CALL( cutResetTimer(hTimer) );
    CUT_SAFE_CALL( cutStartTimer(hTimer) );
    RandomGPU<<<32, 128>>>(d_Random, N_PER_RNG, 777);
    CUT_CHECK_ERROR("RandomGPU() execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
    CUT_SAFE_CALL( cutStopTimer(hTimer) );
    gpuTime = cutGetTimerValue(hTimer);
    printf("Generated samples : %i \n", RAND_N);
    printf("RandomGPU() time  : %f ms\n", gpuTime);
    printf("Samples per second: %E \n", RAND_N / (gpuTime * 0.001));

    printf("Scattering...\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
    CUT_SAFE_CALL( cutResetTimer(hTimer) );
    CUT_SAFE_CALL( cutStartTimer(hTimer) );
    BoxMullerGPU<<<32, 128>>>(d_Random, N_PER_RNG,
                              E0,Ef,hdens,odens,radius);
    CUT_CHECK_ERROR("ScatteringGPU() execution failed\n");
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
    CUT_SAFE_CALL( cutStopTimer(hTimer) );
    gpuTime = cutGetTimerValue(hTimer);
    printf("Scattering time : %f ms\n", gpuTime);
    printf("Samples per second  : %E \n", RAND_N / (gpuTime * 0.001));


    printf("GPU Monte-Carlo simulation...\n");
    initMonteCarloGPU(OPT_N, PATH_N,x);
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
    CUT_SAFE_CALL( cutResetTimer(hTimer) );
    CUT_SAFE_CALL( cutStartTimer(hTimer) );
    MonteCarloGPU(
             h_r2barsqGPU, h_particlesGPU,
             OPT_N, d_Random, PATH_N);

    CUDA_SAFE_CALL( cudaThreadSynchronize() );
    CUT_SAFE_CALL( cutStopTimer(hTimer) );
    closeMonteCarloGPU();
    gpuTime = cutGetTimerValue(hTimer);
    printf("Total GPU time     : %f ms\n", gpuTime);
    printf("Neutron lifetime: : %f %f \n",h_r2barsqGPU[0]/6.0,h_particlesGPU[0]);
    
    printf("Shutting down...\n");
    CUDA_SAFE_CALL( cudaFree(d_Random) );

    CUT_SAFE_CALL( cutDeleteTimer( hTimer) );

}
