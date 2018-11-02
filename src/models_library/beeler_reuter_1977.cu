#include "beeler_reuter_1977.h"
#include <stddef.h>
#include <stdint.h>
#include "model_gpu_utils.h"

extern "C" SET_ODE_INITIAL_CONDITIONS_GPU(set_model_initial_conditions_gpu) 
{

    print_to_stdout_and_file("Using beller_reuter_1977 GPU model\n");

    // execution configuration
    const int GRID  = (num_volumes + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t size = num_volumes*sizeof(real);

    check_cuda_error(cudaMallocPitch((void **) &(*sv), &pitch_h, size, (size_t )NEQ));
    check_cuda_error(cudaMemcpyToSymbol(pitch, &pitch_h, sizeof(size_t)));


    kernel_set_model_inital_conditions <<<GRID, BLOCK_SIZE>>>(*sv, num_volumes);

    check_cuda_error( cudaPeekAtLastError() );
    cudaDeviceSynchronize();
    return pitch_h;

}

extern "C" SOLVE_MODEL_ODES_GPU(solve_model_odes_gpu) {

    // execution configuration
    const int GRID  = ((int)num_cells_to_solve + BLOCK_SIZE - 1)/BLOCK_SIZE;


    size_t stim_currents_size = sizeof(real)*num_cells_to_solve;
    size_t cells_to_solve_size = sizeof(uint32_t)*num_cells_to_solve;

    real *stims_currents_device;
    check_cuda_error(cudaMalloc((void **) &stims_currents_device, stim_currents_size));
    check_cuda_error(cudaMemcpy(stims_currents_device, stim_currents, stim_currents_size, cudaMemcpyHostToDevice));


    //the array cells to solve is passed when we are using and adapative mesh
    uint32_t *cells_to_solve_device = NULL;
    if(cells_to_solve != NULL) {
        check_cuda_error(cudaMalloc((void **) &cells_to_solve_device, cells_to_solve_size));
        check_cuda_error(cudaMemcpy(cells_to_solve_device, cells_to_solve, cells_to_solve_size, cudaMemcpyHostToDevice));
    }
    solve_gpu <<<GRID, BLOCK_SIZE>>>(dt, sv, stims_currents_device, cells_to_solve_device, num_cells_to_solve, num_steps);

    check_cuda_error( cudaPeekAtLastError() );

    check_cuda_error(cudaFree(stims_currents_device));
    if(cells_to_solve_device) check_cuda_error(cudaFree(cells_to_solve_device));

}

__global__ void kernel_set_model_inital_conditions(real *sv, int num_volumes) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadID < num_volumes) {

        *((real * )((char *) sv + pitch * 0) + threadID) = -84.624;        // V
        *((real * )((char *) sv + pitch * 1) + threadID) = 0.011;          // m
        *((real * )((char *) sv + pitch * 2) + threadID) = 0.988;          // h
        *((real * )((char *) sv + pitch * 3) + threadID) = 0.975;          // j
        *((real * )((char *) sv + pitch * 4) + threadID) = 1e-4;           // Cai
        *((real * )((char *) sv + pitch * 5) + threadID) = 0.003;          // d
        *((real * )((char *) sv + pitch * 6) + threadID) = 0.994;          // f
        *((real * )((char *) sv + pitch * 7) + threadID) = 0.0001;         // x1
    }
}

// Solving the model for each cell in the tissue matrix ni x nj
__global__ void solve_gpu(real dt, real *sv, real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps)
{
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;
    int sv_id;

    // Each thread solves one cell model
    if(threadID < num_cells_to_solve) {
        if(cells_to_solve)
            sv_id = cells_to_solve[threadID];
        else
            sv_id = threadID;

        real rDY[NEQ];

        for (int n = 0; n < num_steps; ++n) {

            RHS_gpu(sv, rDY, stim_currents[threadID], sv_id);

            for(int i = 0; i < NEQ; i++) {
                *((real *) ((char *) sv + pitch * i) + sv_id) = dt * rDY[i] + *((real *) ((char *) sv + pitch * i) + sv_id);
            }            

        }

    }
}

inline __device__ void RHS_gpu(real *sv_, real *rDY_, real stim_current, int threadID_) {

    // State variables
    const real V_old_ = *((real*)((char*)sv_ + pitch * 0) + threadID_);
    const real m_old_ = *((real*)((char*)sv_ + pitch * 1) + threadID_);
    const real h_old_ = *((real*)((char*)sv_ + pitch * 2) + threadID_);
    const real j_old_ = *((real*)((char*)sv_ + pitch * 3) + threadID_);
    const real Cai_old_ = *((real*)((char*)sv_ + pitch * 4) + threadID_);
    const real d_old_ = *((real*)((char*)sv_ + pitch * 5) + threadID_);
    const real f_old_ = *((real*)((char*)sv_ + pitch * 6) + threadID_);
    const real x1_old_ = *((real*)((char*)sv_ + pitch * 7) + threadID_);

    // Constants
    const real C = 0.01;
    const real g_na = 4e-2;
    const real E_na = 50;
    const real g_nac = 3e-5;
    const real g_s = 9e-4;

    // Algebraics
    real alpha_m = ( - 1.00000*(V_old_+47.0000))/(exp( - 0.100000*(V_old_+47.0000)) - 1.00000);
    real beta_m =  40.0000*exp( - 0.0560000*(V_old_+72.0000));
    real alpha_h =  0.126000*exp( - 0.250000*(V_old_+77.0000));
    real beta_h = 1.70000/(exp( - 0.0820000*(V_old_+22.5000))+1.00000);
    real alpha_j = ( 0.0550000*exp( - 0.250000*(V_old_+78.0000)))/(exp( - 0.200000*(V_old_+78.0000))+1.00000);
    real beta_j = 0.300000/(exp( - 0.100000*(V_old_+32.0000))+1.00000);
    real alpha_d = ( 0.0950000*exp(- (V_old_ - 5.00000)/100.000))/(1.00000+exp(- (V_old_ - 5.00000)/13.8900));
    real beta_d = ( 0.0700000*exp(- (V_old_+44.0000)/59.0000))/(1.00000+exp((V_old_+44.0000)/20.0000));
    real alpha_f = ( 0.0120000*exp(- (V_old_+28.0000)/125.000))/(1.00000+exp((V_old_+28.0000)/6.67000));
    real beta_f = ( 0.00650000*exp(- (V_old_+30.0000)/50.0000))/(1.00000+exp(- (V_old_+30.0000)/5.00000));
    real alpha_x1 = ( 0.000500000*exp((V_old_+50.0000)/12.1000))/(1.00000+exp((V_old_+50.0000)/17.5000));
    real beta_x1 = ( 0.00130000*exp(- (V_old_+20.0000)/16.6700))/(1.00000+exp(- (V_old_+20.0000)/25.0000));
    real E_s = - 82.3000 -  13.0287*log( Cai_old_*0.00100000);
    real i_s =  g_s*d_old_*f_old_*(V_old_ - E_s);
    real i_na =  ( g_na*powf(m_old_, 3.00000)*h_old_*j_old_+g_nac)*(V_old_ - E_na);
    real i_x1 = ( x1_old_*0.00800000*(exp( 0.0400000*(V_old_+77.0000)) - 1.00000))/exp( 0.0400000*(V_old_+35.0000));
    real i_k1 =  0.00350000*(( 4.00000*(exp( 0.0400000*(V_old_+85.0000)) - 1.00000))/(exp( 0.0800000*(V_old_+53.0000))+exp( 0.0400000*(V_old_+53.0000)))+( 0.200000*(V_old_+23.0000))/(1.00000 - exp( - 0.0400000*(V_old_+23.0000))));
    real i_stim = stim_current;

    // Rates
    rDY_[0] = (i_stim - (i_na+i_s+i_x1+i_k1))/C;
    rDY_[1] = alpha_m*(1.00000 - m_old_) -  beta_m*m_old_;
    rDY_[2] = alpha_h*(1.00000 - h_old_) -  beta_h*h_old_;
    rDY_[3] = alpha_j*(1.00000 - j_old_) -  beta_j*j_old_;
    rDY_[4] = ( - 0.0100000*i_s)/1.00000+ 0.0700000*(0.000100000 - Cai_old_);
    rDY_[5] = alpha_d*(1.00000 - d_old_) -  beta_d*d_old_;
    rDY_[6] = alpha_f*(1.00000 - f_old_) -  beta_f*f_old_;
    rDY_[7] = alpha_x1*(1.00000 - x1_old_) -  beta_x1*x1_old_;

}
