/* 21 Mar 2009
 * Matt Anderson
 * Code for neutrino transport
 */
static const int MAX_OPTIONS = 256;
static const int THREAD_N = 256;

// Determined empirically for G80 GPUs
// Experiment with this value to obtain higher performance on 
// GPUs with more or fewer multiprocessors
static const int MULTIBLOCK_THRESHOLD = 8192;

static float *d_Sum;
static float *h_Sum;

__constant__ float d_X[MAX_OPTIONS];

////////////////////////////////////////////////////////////////////////////////////
// Compute an efficient number of CTAs to use per option for the multiblock
// version (MonteCarloKernel()).  These numbers were determined via experimentation
// on G80 GPUs.  Optimal values may be different on other GPUs.
////////////////////////////////////////////////////////////////////////////////////
unsigned int computeNumCTAs(unsigned int optN)
{
    return (optN < 16) ? 64 : 16;
}

////////////////////////////////////////////////////////////////////////////////////
// Allocate intermediate strage for the Monte Carlo integration
////////////////////////////////////////////////////////////////////////////////////
void initMonteCarloGPU(unsigned int optN, unsigned int pathN, float *h_X){
    
    unsigned int ratio = pathN / optN;
    unsigned int accumSz = 2 * sizeof(float);
    if (ratio >= MULTIBLOCK_THRESHOLD) 
    {
        // in this case we need to store a number of partial sums per thread block
        unsigned int accumN = computeNumCTAs(optN) * THREAD_N;
        accumSz *= accumN;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Use OS-pinned memory on host side. Allocation takes slightly more time,
    // But OS-pinned<==>device memory transfers are faster depending on 
    // the system configuration. Refer to the programming guide and 
    // bandwidthTest CUDA SDK sample for performance comparisons on the
    // particular system.
    ////////////////////////////////////////////////////////////////////////////
    CUDA_SAFE_CALL( cudaMallocHost((void **)&h_Sum,  optN * 2 * sizeof(float)) );
    CUDA_SAFE_CALL( cudaMalloc((void **)&d_Sum,  accumSz * optN) );

    // Initialize the per-option data in constant arrays accessible by MonteCarloKernel()
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(d_X, h_X, MAX_OPTIONS * sizeof(float)) );
}

void closeMonteCarloGPU(void){
    CUDA_SAFE_CALL( cudaFree(d_Sum)      );
    CUDA_SAFE_CALL( cudaFreeHost(h_Sum)  );
}

// Needed by the optimized sum reduction for correct execution in device emulation
#ifdef __DEVICE_EMULATION__
#define SYNC __syncthreads()
#else
#define SYNC
#endif

////////////////////////////////////////////////////////////////////////////////////
// Given shared memory with blockSize valus and blockSize squared values,
// This function computes the sum of each array.  The result for each array
// is stored in element 0 of tha array.
////////////////////////////////////////////////////////////////////////////////////
template <unsigned int blockSize>
__device__ void 
sumReduceSharedMem(float *sum, float *sum2)
{
    unsigned int tid = threadIdx.x;

    // do reduction in shared mem
    if (blockSize >= 512) { if (tid < 256) { sum[tid] += sum[tid + 256]; sum2[tid] += sum2[tid + 256]; } __syncthreads(); }
    if (blockSize >= 256) { if (tid < 128) { sum[tid] += sum[tid + 128]; sum2[tid] += sum2[tid + 128]; } __syncthreads(); }
    if (blockSize >= 128) { if (tid <  64) { sum[tid] += sum[tid +  64]; sum2[tid] += sum2[tid +  64]; } __syncthreads(); }
    
#ifndef __DEVICE_EMULATION__
    if (tid < 32)
#endif
    {
        if (blockSize >=  64) { sum[tid] += sum[tid + 32]; sum2[tid] += sum2[tid + 32]; SYNC; }
        if (blockSize >=  32) { sum[tid] += sum[tid + 16]; sum2[tid] += sum2[tid + 16]; SYNC; }
        if (blockSize >=  16) { sum[tid] += sum[tid +  8]; sum2[tid] += sum2[tid +  8]; SYNC; }
        if (blockSize >=   8) { sum[tid] += sum[tid +  4]; sum2[tid] += sum2[tid +  4]; SYNC; }
        if (blockSize >=   4) { sum[tid] += sum[tid +  2]; sum2[tid] += sum2[tid +  2]; SYNC; }
        if (blockSize >=   2) { sum[tid] += sum[tid +  1]; sum2[tid] += sum2[tid +  1]; SYNC; }
    }
}

////////////////////////////////////////////////////////////////////////////////////
// Compute the final sum and sum-of-squares of ACCUM_N values for each option using 
// an optimized parallel tree reduction.  Calls sumReduceSharedMem
////////////////////////////////////////////////////////////////////////////////////
template <unsigned int blockSize>
__global__ void
sumReduction(float *g_odata, float *g_idata, unsigned int blockDataSize)
{
    __shared__ float sum[blockSize];
    __shared__ float sum2[blockSize]; // sum of squares

    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*2*blockDataSize + threadIdx.x;
    sum[tid] = 0;
    sum2[tid] = 0;

    // we reduce multiple elements per thread.  The number is determined by the 
    // number of active thread blocks (via gridSize).  More blocks will result
    // in a larger gridSize and therefore fewer elements per thread
    for (int count = 0; count < blockDataSize/blockSize; count++)
    {
        sum[tid]  += g_idata[i];  
        sum2[tid] += g_idata[i + blockDataSize];
        i += blockSize;        
    } 
    __syncthreads();

    // do reduction in shared mem
    sumReduceSharedMem<blockSize>(sum, sum2);
    
    // write result for this block to global mem 
    if (tid == 0) 
    {
        g_odata[2 * blockIdx.x]     =  sum[0];
        g_odata[2 * blockIdx.x + 1] = sum2[0];

    }
}

////////////////////////////////////////////////////////////////////////////////////
// This kernel computes partial integrals over all paths using a multiple thread 
// blocks per option.  It is used when a single thread block per option would not
// be enough to keep the GPU busy.  Execution of this kernel is followed by
// a sumReduction() to get the complete integral for each option.
////////////////////////////////////////////////////////////////////////////////////
__global__ void MonteCarloKernel(
    float *d_Sum,    //Partial sums (+sum of squares) destination
    int   accumN,    //Partial sums (sum of squares) count
    float *d_Random, //N(0, 1) random samples array
    int   pathN
){
    const int tid      = blockDim.x * blockIdx.x + threadIdx.x;
    const int optIndex = blockIdx.y;
    const int threadN  = blockDim.x * gridDim.x;
    float r2;

    for(int iAccum = tid; iAccum < accumN; iAccum += threadN) {
      float sum = 0, sum2 = 0;
      for (int iPath = iAccum; iPath < pathN; iPath += accumN) {
        r2 = d_Random[iPath];
        if ( r2 > 0.0 ) {
          sum += r2;
          sum2 += 1.0;
        }
      }
      d_Sum[optIndex * 2 * accumN + iAccum +      0] = sum;
      d_Sum[optIndex * 2 * accumN + iAccum + accumN] = sum2;
    }


  //  for(int iAccum = tid; iAccum < accumN; iAccum += threadN) {
  //    float sum = 0, sum2 = 0;
  //    for (int iPath = iAccum; iPath < pathN; iPath += accumN) {
  //      myrand = d_Random[iPath];
  //
  //       sum += myrand;
  //      sum2 += 1.0;
  //    }
  //    d_Sum[optIndex * 2 * accumN + iAccum +      0] = sum;
  //    d_Sum[optIndex * 2 * accumN + iAccum + accumN] = sum2;
  //  }
}

////////////////////////////////////////////////////////////////////////////////////
// Here we choose between two different methods for performing Monte Carlo 
// integration on the GPU.  When the ratio of paths to options is lower than a 
// threshold (8192 determined empirically for G80 GPUs -- a different threshold is 
// likely to be better on other GPUs!), we run a single kernel that runs one thread 
// block per option and integrates all samples for that option.  This is 
// MonteCarloKernelOneBlockPerOption().  When the ratio is high, then we need more
// threads to better hide memory latency.  In this case we run multiple thread
// blocks per option and compute partial sums stored in the d_Sum array.  This is
// MonteCarloKernel().  These partial sums are then reduced to a final sum using  
// a parallel reduction (sumReduction()).  In both cases, the sum and sum of 
// squares for each option is read back to the host where the final callResult and
// confidenceWidth are computed.  These are computed on the CPU because doing so on
// the GPU would leave most threads idle.
////////////////////////////////////////////////////////////////////////////////////
void MonteCarloGPU(
    float *r2barsq,      // r^2 bar squared
    float *particles,    // number of particles
    int optN,               //Input options count
    float *d_Random,        //(0, 1) random samples array
    int pathN
){

    int ctaN = computeNumCTAs(optN);
    int accumN = ctaN * THREAD_N;
    dim3 gridDim(ctaN, optN, 1);
    MonteCarloKernel<<<gridDim, THREAD_N>>>(
        d_Sum, accumN, d_Random, pathN);
    CUT_CHECK_ERROR("MonteCarloKernel() execution failed\n");
       
    // Perform a parallel sum reduction on the device to reduce the ACCUM_N values
    // generated per option to a single value (actually two values -- sum and sum of 
    // squares).  This reduction is very efficient on the GPU.
    sumReduction<128><<<optN, 128>>>(d_Sum, d_Sum, accumN);
    CUT_CHECK_ERROR("sumReduction() execution failed\n");
    
    // Read back only the sum and sum of squares for each option to the CPU.
    CUDA_SAFE_CALL( cudaMemcpy(h_Sum, d_Sum, 2 * optN * sizeof(float), cudaMemcpyDeviceToHost) );

    // Compute final statistics
    for(int opt = 0; opt < optN; opt++){
        float sum  = h_Sum[2*opt];
        float sum2 = h_Sum[2*opt+1];

        //Derive average from the total sum and discount by riskfree rate 
        r2barsq[opt] = (float)( sum / sum2);
        //r2barsq[opt] = (float)( sum );
        //waveResult[opt] = (float)( sum );

        //Standard deviation
        //double stdDev = sqrt(((double)pathN * sum2 - sum * sum)/ ((double)pathN * (double)(pathN - 1)));

        //Confidence width; in 95% of all cases theoretical value lies within these borders
        //confidenceWidth[opt] = (float)(1.96 * stdDev / sqrt((double)pathN));
        particles[opt] = (float)(sum2);
    }
}
