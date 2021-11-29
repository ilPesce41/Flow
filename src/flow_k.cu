#include <iostream>
#include <math.h>
#include "flow_k.cuh"
#include "spatial_distance.cuh"
#include "flow.hpp"
#include <algorithm>
#include <curand.h>
#include <curand_kernel.h>
#include <stdio.h>


#define PI 3.1416


__device__ float local_k_function(float * distance_matrix,int index, float r,float A,int length)
{
    int count = 0;
    for(int i=0;i<length;i++)
    {
        if(distance_matrix[index*length + i]<r)
        {
            count++;
        }
    }
    count--; // ignore self
    return sqrt(A*count/(PI*length-1));
}

__global__ void calculate_k_function(float * dist_matrix, float * sig_vec, float r,float A, int length){
	int i = blockDim.x*blockIdx.x + threadIdx.x;
    if(i>=length)
        return;
    float sig_value = local_k_function(dist_matrix,i,r,A,length);
    sig_vec[i] = sig_value;
}

__global__ void vec_shuffle(float * vec, int length){
	int i = blockDim.x*blockIdx.x + threadIdx.x;
    if(i>=length){
        return;
    }
    curandState state;
    curand_init((unsigned long long)clock() + i, 0, 0, &state);
    int j = (int)(curand_uniform(&state)*(length-1));
    float eps = curand_uniform(&state)*0.1 - 0.05;
    float value = vec[j] + eps;
    __syncthreads();
    vec[i] = value;
    __syncthreads();

}

float get_min(float * arr,int length)
{
    float min_value = arr[0];
    for(int i=0;i<length;i++)
    {
        if(arr[i]<min_value && abs(arr[i])>1e-3)
            min_value = arr[i];
    }
    return min_value;
}

float get_max(float * arr,int length)
{
    float max_value = arr[0];
    for(int i=0;i<length;i++)
    {
        if(arr[i]>max_value && abs(arr[i])>1e-3)
            max_value = arr[i];
    }
    return max_value;
}


// void process_flow_k(FlowData flow, int num_iter, int func_type,float alpha,size_t shared_mem_size)
std::vector<int> process_flow_k(FlowData flow, float *d_sx, float *d_sy, float *d_dx, float *d_dy, float *d_L, float * dist_matrix_cpu, float * dist_matrix_gpu, int num_iter, int func_type,float alpha,size_t shared_mem_size,float radius)
{ 

    float * sig_vec_real;
    float * sig_vec_cpu;
    sig_vec_real = (float *)malloc(sizeof(float)*flow.length);
    sig_vec_cpu = (float *)malloc(sizeof(float)*flow.length);

    float min_sx = get_min(flow.sx, flow.length);
    float min_sy = get_min(flow.sy, flow.length);
    float min_dx = get_min(flow.dx, flow.length);
    float min_dy = get_min(flow.dy, flow.length);
    float max_sx = get_max(flow.sx, flow.length);
    float max_sy = get_max(flow.sy, flow.length);
    float max_dx = get_max(flow.dx, flow.length);
    float max_dy = get_max(flow.dy, flow.length);

    
    float min_x, min_y, max_x, max_y;
    min_x = std::min(min_sx,min_dx);
    min_y = std::min(min_sy,min_dy);
    max_x = std::max(max_sx,max_dx);
    max_y = std::max(max_sy,max_dy);
    std::cout << max_x<< " " << min_x << std::endl;
    std::cout << max_y<< " " << min_y << std::endl;
    float A = (max_x-min_x)*(max_y-min_y);
    std::cout << "Area: " << A << std::endl;

    float * sig_vec;
    cudaMalloc((void**)&sig_vec, flow.length*sizeof(float));
    calculate_k_function<<< ceil((float)flow.length/128) , 128 >>>(dist_matrix_gpu,sig_vec,radius,A,flow.length);
	cudaMemcpy(sig_vec_real, sig_vec, flow.length*sizeof(float), cudaMemcpyDeviceToHost);
    for(int j=0;j<10;j++)
        {
            std::cout << sig_vec_real[j] << " ";
        }
    std::cout << std::endl;
    std::cout << std::endl;

    //Generate synthetic data
    float upper_envelope, lower_envelope;
    for (int ii=0;ii<num_iter;ii++)
    {
        cudaMemcpy(d_sx, flow.sx, flow.length*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_sy, flow.sy, flow.length*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_dx, flow.dx, flow.length*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_dy, flow.dy, flow.length*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_L, flow.L, flow.length*sizeof(float), cudaMemcpyHostToDevice);
        cudaDeviceSynchronize();
        vec_shuffle<<< ceil((float)flow.length/64) , 64 >>>(d_sx, flow.length);
        float * tmp_arr = (float *)malloc(flow.length*sizeof(float));
        cudaMemcpy(tmp_arr, d_sx, flow.length*sizeof(float), cudaMemcpyDeviceToHost);
        // for(int i=0;i<10;i++)
        // {
        //     std::cout << tmp_arr[i] << " ";
        // }
        // cout << endl;
        cudaDeviceSynchronize();
        // printf("GPUassert: %s \n", cudaGetErrorString(cudaGetLastError()));
        vec_shuffle<<< ceil((float)flow.length/64) , 64 >>>(d_sy, flow.length);
        vec_shuffle<<< ceil((float)flow.length/64) , 64 >>>(d_dx, flow.length);
        vec_shuffle<<< ceil((float)flow.length/64) , 64 >>>(d_dy, flow.length);
        // std::cout<< shared_mem_size << std::endl;
        calculate_spatial_distance_matrix <<< ceil((float)flow.length/128) , 128, shared_mem_size >>> (d_sx,d_sy,d_dx,d_dy,d_L,dist_matrix_gpu,flow.length,alpha,func_type);
        // cudaMemcpy(dist_matrix_cpu, dist_matrix_gpu, flow.length*flow.length*sizeof(float), cudaMemcpyDeviceToHost);
        // for(int i=0;i<3;i++)
        // {
        //     for(int j=0;j<10;j++)
        //     {
        //         std::cout << dist_matrix_cpu[i*flow.length+j] << " ";
        //     }
        //     std::cout << std::endl;
        // }
        calculate_k_function<<< ceil((float)flow.length/128) , 128 >>>(dist_matrix_gpu,sig_vec,radius,A,flow.length);
        cudaMemcpy(sig_vec_cpu, sig_vec, flow.length*sizeof(float), cudaMemcpyDeviceToHost);
        // for(int j=0;j<10;j++)
        // {
        //     std::cout << sig_vec_cpu[j] << " ";
        // }
        // std::cout << std::endl;
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

    std::vector<int> output;
    for(int i=0;i<flow.length;i++)
    {
        if(sig_vec_real[i]<upper_envelope)
        {
            output.push_back(i);
        }
    }
    return output;
}
