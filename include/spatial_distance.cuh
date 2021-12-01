/*
Cole Hill
University of South Florida
Computer Science and Engineering
Programming Massively Parallel Systems
Fall 2021
*/
#pragma once

__device__ float flow_distance(float alpha, float xi, float yi, float ui, float vi,float xj, float yj, float uj, float vj);
__device__ float flow_disimilarity(float alpha, float xi, float yi, float ui, float vi,float Li, float xj, float yj, float uj, float vj,float Lj);
__global__ void calculate_spatial_distance_matrix(float * d_sx,float * d_sy,float * d_dx,float * d_dy,float * d_L,float * dist_matrix, int length, float alpha, int func);
