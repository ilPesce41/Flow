/*
Cole Hill
University of South Florida
Computer Science and Engineering
Programming Massively Parallel Systems
Fall 2021
*/
#include <iostream>
#include <math.h>
#include "flow_k.cuh"
#include "spatial_distance.cuh"
#include "flow.hpp"
#include <algorithm>
#include <curand.h>
#include <curand_kernel.h>
#include <stdio.h>
#include "utils.hpp"

#define PI 3.1416
#define BLOCK_DIM 128


/*
Implements ipleys K-function adapted to flows
Calculation taken from http://ceadserv1.nku.edu/longa//geomed/ppa/doc/localK/localK.htm

distance_matrix is a precomputed distance matrix between all flows
*/
__device__ float local_k_function(float * distance_matrix,int index, float r,float A,int length)
{   
    //Count number of flows located within spatial radius threshold r
    int count = 0;
    for(int i=0;i<length;i++)
    {
        if(distance_matrix[index*length + i]<r)
        {
            count++;
        }
    }
    count--; // ignore self
    //return k value
    return sqrt(A*count/(PI*length-1));
}

/*
Implements cross K-function adapted to flows
Flow Cross K-function: a bivariate flow analytical method
https://onlinelibrary.wiley.com/doi/epdf/10.1111/gean.12100

distance_matrix is a precomputed distance matrix between all flows
piv is index where classes switch
*/
__device__ float local_cross_k_function(float * distance_matrix,int index, float r,float A,int length,int piv)
{
    int count = 0;
    int start,end;
    //check first half of matrix
    if(index>piv)
    {
        start = 0;
        end = piv;
    }
    //check second half of matrix
    else
    {
        start = piv;
        end = length;
    }

    //Calculate ripley's k function
    for(int i=start;i<end;i++)
    {
        if(distance_matrix[index*length + i]<r)
        {
            count++;
        }
    }
    return sqrt(A*count/(PI*length-1));
}

/*
GPU kernel for calculating K funtion for all flows
*/
__global__ void calculate_k_function(float * dist_matrix, float * sig_vec, float r,float A, int length){
	int i = blockDim.x*blockIdx.x + threadIdx.x;
    if(i>=length)
        return;
    float sig_value = local_k_function(dist_matrix,i,r,A,length);
    sig_vec[i] = sig_value;
}

/*
GPU kernel for calculating cross-K funtion for all flows
*/
__global__ void calculate_cross_k_function(float * dist_matrix, float * sig_vec, float r,float A, int length, int piv){
	int i = blockDim.x*blockIdx.x + threadIdx.x;
    if(i>=length)
        return;
    float sig_value = local_cross_k_function(dist_matrix,i,r,A,length,piv);
    sig_vec[i] = sig_value;
}

/*
GPU kernel for shuffling vectors
*/
__global__ void vec_shuffle(float * vec, int length){
	int i = blockDim.x*blockIdx.x + threadIdx.x;
    if(i>=length){
        return;
    }
    // Init random number generator
    curandState state;
    curand_init((unsigned long long)clock() + i, 0, 0, &state);
    
    // Pick alternate value to take
    int j = (int)(curand_uniform(&state)*(length-1));
    float eps = curand_uniform(&state)*0.1 - 0.05;
    float value = vec[j] + eps;
    __syncthreads();
    //Update value
    vec[i] = value;
    __syncthreads();
}

/*
GPU kernel for shuffling class labels of flows
*/
__global__ void shuffle_label(float * d_sx, float * d_sy, float * d_dx, float * d_dy, float * d_L, int length){
	int i = blockDim.x*blockIdx.x + threadIdx.x;
    if(i>=length){
        return;
    }
    // Init random number generator
    curandState state;
    curand_init((unsigned long long)clock() + i, 0, 0, &state);
    
    // Pick alternate value to take
    int j = (int)(curand_uniform(&state)*(length-1));
    float eps = curand_uniform(&state)*0.1 - 0.05;
    float sx,sy,dx,dy,L;
    sx = d_sx[j];
    sy = d_sy[j];
    dx = d_dx[j];
    dy = d_dy[j];
    L = d_L[j];
    //Update value
    __syncthreads();
    d_sx[i] = sx;
    d_sy[i] = sy;
    d_dx[i] = dx;
    d_dy[i] = dy;
    d_L[i] = L;
    __syncthreads();
}

/*
Calculate Ripleys k value for each flow and run monte carlo simulation to 
test significance. Returns significant flows
*/
std::vector<int> process_flow_k(FlowData flow, float *d_sx, float *d_sy, float *d_dx, float *d_dy, float *d_L, float * dist_matrix_gpu, int num_iter, int func_type,float alpha,size_t shared_mem_size,float radius)
{ 
    // K-value vectors
    float * sig_vec_real;
    float * sig_vec_cpu;
    sig_vec_real = (float *)malloc(sizeof(float)*flow.length);
    sig_vec_cpu = (float *)malloc(sizeof(float)*flow.length);

    // Copy flow data to GPU
    cudaMemcpy(d_sx, flow.sx, flow.length*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sy, flow.sy, flow.length*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_dx, flow.dx, flow.length*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_dy, flow.dy, flow.length*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_L, flow.L, flow.length*sizeof(float), cudaMemcpyHostToDevice);
    // Calculate spatial distance matrix
    calculate_spatial_distance_matrix <<< ceil((float)flow.length/BLOCK_DIM) , BLOCK_DIM, shared_mem_size >>> (d_sx,d_sy,d_dx,d_dy,d_L,dist_matrix_gpu,flow.length,func_type,alpha);

    // Calculate k values for real data
    float * sig_vec;
    cudaMalloc((void**)&sig_vec, flow.length*sizeof(float));
    calculate_k_function<<< ceil((float)flow.length/128) , 128 >>>(dist_matrix_gpu,sig_vec,radius,flow.area,flow.length);
	cudaMemcpy(sig_vec_real, sig_vec, flow.length*sizeof(float), cudaMemcpyDeviceToHost);

    //Generate synthetic data (monte carlo step)
    float upper_envelope, lower_envelope;
    for (int ii=0;ii<num_iter;ii++)
    {
        // Copy flow data to GPU
        cudaMemcpy(d_sx, flow.sx, flow.length*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_sy, flow.sy, flow.length*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_dx, flow.dx, flow.length*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_dy, flow.dy, flow.length*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_L, flow.L, flow.length*sizeof(float), cudaMemcpyHostToDevice);
        cudaDeviceSynchronize();
        // Shuffle endpoints
        vec_shuffle<<< ceil((float)flow.length/64) , 64 >>>(d_sx, flow.length);
        cudaDeviceSynchronize();
        vec_shuffle<<< ceil((float)flow.length/64) , 64 >>>(d_sy, flow.length);
        vec_shuffle<<< ceil((float)flow.length/64) , 64 >>>(d_dx, flow.length);
        vec_shuffle<<< ceil((float)flow.length/64) , 64 >>>(d_dy, flow.length);
        // Calculate spatial distance matrix
        calculate_spatial_distance_matrix <<< ceil((float)flow.length/128) , 128, shared_mem_size >>> (d_sx,d_sy,d_dx,d_dy,d_L,dist_matrix_gpu,flow.length,alpha,func_type);
        // Calculate k values for synthetic data
        calculate_k_function<<< ceil((float)flow.length/128) , 128 >>>(dist_matrix_gpu,sig_vec,radius,flow.area,flow.length);
        cudaMemcpy(sig_vec_cpu, sig_vec, flow.length*sizeof(float), cudaMemcpyDeviceToHost);
        // Find upper and lower significance values
        if(ii>0)
        {
            upper_envelope = max(upper_envelope,get_max(sig_vec_cpu,flow.length));
            lower_envelope = min(lower_envelope,get_min(sig_vec_cpu,flow.length));
        }
        else
        {
            upper_envelope = get_max(sig_vec_cpu,flow.length);
            lower_envelope = get_min(sig_vec_cpu,flow.length);
        }
    }
    //Return flows with k-value above significance threshold
    std::vector<int> output;
    for(int i=0;i<flow.length;i++)
    {
        if(sig_vec_real[i]>upper_envelope)
        {
            output.push_back(i);
        }
    }
    return output;
}

/*
Calculate cross flow Ripleys k value for each flow and run monte carlo simulation to 
test significance. Returns significant flows
*/
std::vector<int> process_cross_flow_k(FlowData flow_1, FlowData flow_2, int num_iter, int func_type,float alpha,size_t shared_mem_size,float radius)
{ 
    
    cudaDeviceReset();
    
    float * sig_vec_real;
    float * sig_vec_cpu;
    int length = flow_1.length + flow_2.length;

    // copy data to GPU 
	float *d_sx, *d_sy, *d_dx, *d_dy, *d_L;
	cudaMalloc((void**)&d_sx, length*sizeof(float));
	cudaMalloc((void**)&d_sy, length*sizeof(float));
	cudaMalloc((void**)&d_dx, length*sizeof(float));
	cudaMalloc((void**)&d_dy, length*sizeof(float));
	cudaMalloc((void**)&d_L, length*sizeof(float));
    //Copy first flow
	cudaMemcpy(d_sx, flow_1.sx, flow_1.length*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_sy, flow_1.sy, flow_1.length*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_dx, flow_1.dx, flow_1.length*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_dy, flow_1.dy, flow_1.length*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_L, flow_1.L, flow_1.length*sizeof(float), cudaMemcpyHostToDevice);
    //Copy second flow
    cudaMemcpy(d_sx + flow_1.length, flow_2.sx, flow_2.length*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_sy + flow_1.length, flow_2.sy, flow_2.length*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_dx + flow_1.length, flow_2.dx, flow_2.length*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_dy + flow_1.length, flow_2.dy, flow_2.length*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_L + flow_1.length, flow_2.L, flow_2.length*sizeof(float), cudaMemcpyHostToDevice);
    
    // K-value vectors
    sig_vec_real = (float *)malloc(sizeof(float)*length);
    sig_vec_cpu = (float *)malloc(sizeof(float)*length);
    
    // Calculate study area
    float max_x = max(flow_1.max_x,flow_2.max_x);
    float max_y = max(flow_1.max_y,flow_2.max_y);
    float min_x = min(flow_1.min_x,flow_2.min_x);
    float min_y = min(flow_1.min_y,flow_2.min_y);
    float A = (max_x-min_x)*(max_y-min_y);

    // Calculate spatial distance matrix
    float * dist_matrix_gpu;
	cudaMalloc((void**)&dist_matrix_gpu, length*length*sizeof(float));
	cudaMemset(dist_matrix_gpu, 0, length*length*sizeof(float));
    calculate_spatial_distance_matrix <<< ceil((float)length/128) , 128, shared_mem_size >>> (d_sx,d_sy,d_dx,d_dy,d_L,dist_matrix_gpu,length,func_type,alpha);

    // Calculate k values for real data
    float * sig_vec;
    cudaMalloc((void**)&sig_vec, length*sizeof(float));
    calculate_cross_k_function<<< ceil((float)length/128) , 128 >>>(dist_matrix_gpu,sig_vec,radius,A,length,flow_1.length);
	cudaMemcpy(sig_vec_real, sig_vec, length*sizeof(float), cudaMemcpyDeviceToHost);

    //Generate synthetic data
    float upper_envelope_1, lower_envelope_1;
    float upper_envelope_2, lower_envelope_2;
    for (int ii=0;ii<num_iter;ii++)
    {
        //Copy first flow
        cudaMemcpy(d_sx, flow_1.sx, flow_1.length*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_sy, flow_1.sy, flow_1.length*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_dx, flow_1.dx, flow_1.length*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_dy, flow_1.dy, flow_1.length*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_L, flow_1.L, flow_1.length*sizeof(float), cudaMemcpyHostToDevice);
        //Copy second flow
        cudaMemcpy(d_sx + flow_1.length, flow_2.sx, flow_2.length*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_sy + flow_1.length, flow_2.sy, flow_2.length*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_dx + flow_1.length, flow_2.dx, flow_2.length*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_dy + flow_1.length, flow_2.dy, flow_2.length*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_L + flow_1.length, flow_2.L, flow_2.length*sizeof(float), cudaMemcpyHostToDevice);
        cudaDeviceSynchronize();
        // Shuffle flow labels
        shuffle_label<<< ceil((float)length/64) , 64 >>>(d_sx, d_sy, d_dx, d_dy, d_L, length);
        cudaDeviceSynchronize();
        // Recalate spatial distance matrix and cross k functions for synth
        // data
        calculate_spatial_distance_matrix <<< ceil((float)length/128) , 128, shared_mem_size >>> (d_sx,d_sy,d_dx,d_dy,d_L,dist_matrix_gpu,length,alpha,func_type);
        calculate_cross_k_function<<< ceil((float)length/128) , 128 >>>(dist_matrix_gpu,sig_vec,radius,A,length,flow_1.length);
        cudaMemcpy(sig_vec_cpu, sig_vec, length*sizeof(float), cudaMemcpyDeviceToHost);
        // Find upper and lower significance values
        if(ii>0)
        {
            upper_envelope_1 = max(upper_envelope_1,get_max(sig_vec_cpu,0,flow_1.length));
            lower_envelope_1 = min(lower_envelope_1,get_min(sig_vec_cpu,0,flow_1.length));
            upper_envelope_2 = max(upper_envelope_2,get_max(sig_vec_cpu,flow_1.length,length));
            lower_envelope_2 = min(lower_envelope_2,get_min(sig_vec_cpu,flow_1.length,length));
        }
        else
        {
            upper_envelope_1 = get_max(sig_vec_cpu,0,flow_1.length);
            lower_envelope_1 = get_min(sig_vec_cpu,0,flow_1.length);
            upper_envelope_2 = get_max(sig_vec_cpu,flow_1.length,length);
            lower_envelope_2 = get_min(sig_vec_cpu,flow_1.length,length);
        }
    }
    
    //Return flows with k-value above significance threshold
    std::vector<int> output;
    for(int i=0;i<flow_1.length;i++)
    {
        if(sig_vec_real[i]>upper_envelope_1)
        {
            output.push_back(i);
        }
    }
    output.push_back(-1);
    for(int i=flow_1.length;i<length;i++)
    {
        if(sig_vec_real[i]>upper_envelope_2)
        {
            output.push_back(i-flow_1.length);
        }
    }
    return output;
}
