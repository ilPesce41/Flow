#pragma once
#include "flow.hpp"

__device__ float local_k_function(float * distance_matrix,int index, float r,float A,int length);
__global__ void calculate_k_function(float * dist_matrix, float * sig_vec, float r,float A, int length);
__global__ void vec_shuffle(float * vec, int length);
float get_min(float * arr,int length);
float get_max(float * arr,int length);
std::vector<int>  process_flow_k(FlowData flow, float *d_sx, float *d_sy, float *d_dx, float *d_dy, float *d_L, float * dist_matrix_gpu, int num_iter, int func_type,float alpha,size_t shared_mem_size,float radius);
std::vector<int> process_cross_flow_k(FlowData flow_1, FlowData flow_2, int num_iter, int func_type,float alpha,size_t shared_mem_size,float radius);