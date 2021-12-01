/*
Cole Hill
University of South Florida
Computer Science and Engineering
Programming Massively Parallel Systems
Fall 2021
*/
#include <stdio.h>

/*
  distance function on GPU
*/
__device__ float flow_distance(float alpha, float xi, float yi, float ui, float vi,float xj, float yj, float uj, float vj) {
	float x_diff = xi-xj;
	float y_diff = yi-yj;
	float u_diff = ui-uj;
	float v_diff = vi-vj;
    
    return sqrt(alpha*(x_diff*x_diff + y_diff*y_diff) + (2-alpha)*(u_diff*u_diff + v_diff*v_diff));
}

/*
  distance function on GPU
*/
__device__ float flow_disimilarity(float alpha, float xi, float yi, float ui, float vi,float Li, float xj, float yj, float uj, float vj,float Lj) {
	float x_diff = xi-xj;
	float y_diff = yi-yj;
	float u_diff = ui-uj;
	float v_diff = vi-vj;
    
    return sqrt((alpha*(x_diff*x_diff + y_diff*y_diff) + (2-alpha)*(u_diff*u_diff + v_diff*v_diff))/(Li*Lj));
}

/*
  distance function on GPU
*/
__device__ float flow_max_distance(float xi, float yi, float ui, float vi, float xj, float yj, float uj, float vj) {
	float x_diff = xi-xj;
	float y_diff = yi-yj;
	float u_diff = ui-uj;
	float v_diff = vi-vj;
    
    return max(sqrt(x_diff*x_diff + y_diff*y_diff),sqrt(u_diff*u_diff + v_diff*v_diff));
}

__global__ void calculate_spatial_distance_matrix(float * d_sx,float * d_sy,float * d_dx,float * d_dy,float * d_L,float * dist_matrix, int length, float alpha, int func){
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	float dist;
	if(i>=length)
    {
        return;
    }

	// Initialize the b-th block into memory (need j)
	extern __shared__ float s[];
	float *B_x = s;
	float *B_y = (float*)&s[blockDim.x*1];
	float *B_u = (float*)&s[blockDim.x*2];
	float *B_v = (float*)&s[blockDim.x*3];
	float *B_L = (float*)&s[blockDim.x*4];
	float *R_x = (float*)&s[blockDim.x*5];
	float *R_y = (float*)&s[blockDim.x*6];
	float *R_u = (float*)&s[blockDim.x*7];
	float *R_v = (float*)&s[blockDim.x*8];
	float *R_L = (float*)&s[blockDim.x*9];

	__syncthreads();

	// Load B-th block into memory
	float x1 = d_sx[i];
	float y1 = d_sy[i];
	float u1 = d_dx[i];
	float v1 = d_dy[i];
	float L1 = d_L[i];
	B_x[threadIdx.x] = x1;
	B_y[threadIdx.x] = y1;
	B_u[threadIdx.x] = u1;
	B_v[threadIdx.x] = v1;
	B_L[threadIdx.x] = L1;

	// inter-block differences
	//for i = b+1 to M 
	// load block R (next block along grid)
	// syncthreads
	// for 0 to block size
	//   calculate distance and update histogram
	for(int iter=blockIdx.x+1;iter<((int)(length/blockDim.x)+1);iter++)
	{
		// Load B-th block into memory
		int j = blockDim.x*iter + threadIdx.x;
		if(j<length)
		{
			R_x[threadIdx.x] = d_sx[j];
			R_y[threadIdx.x] = d_sy[j];
			R_u[threadIdx.x] = d_dx[j];
			R_v[threadIdx.x] = d_dy[j];
			R_L[threadIdx.x] = d_L[j];
		}
		__syncthreads();

		for(int t=0;t<blockDim.x;t++)
		{
			int j = blockDim.x*iter + t;
			if (i<length && j<length)
			{
                if(func==0)
                {
				    dist = flow_distance(alpha,x1,y1,u1,v1,R_x[t],R_y[t],R_u[t],R_v[t]);
                }
                else if(func==1)
                {
				    dist = flow_disimilarity(alpha,x1,y1,u1,v1,L1,R_x[t],R_y[t],R_u[t],R_v[t],R_L[t]);
                }
				else
				{
				    dist = flow_max_distance(x1,y1,u1,v1,R_x[t],R_y[t],R_u[t],R_v[t]);
				}
                dist_matrix[i*length + j] = dist;
                dist_matrix[j*length + i] = dist;
			}
		}
		__syncthreads();
	}

	// intra-block distances
	// for i = t+1 to B do
	//  calculate distance and update histogram
	for(int t=threadIdx.x+1;t<blockDim.x;t++)
		{
			int j = blockDim.x*blockIdx.x + t;
			if (i<length && j<length)
			{
                if(func==0)
                {
				    dist = flow_distance(alpha,x1,y1,u1,v1,B_x[t],B_y[t],B_u[t],B_v[t]);
                }
                else if(func==1)
                {
				    dist = flow_disimilarity(alpha,x1,y1,u1,v1,L1,B_x[t],B_y[t],B_u[t],B_v[t],B_L[t]);
                }
				else
				{
				    dist = flow_max_distance(x1,y1,u1,v1,B_x[t],B_y[t],B_u[t],B_v[t]);
				}
                dist_matrix[i*length + j] = dist;
                dist_matrix[j*length + i] = dist;
			}
		}
	__syncthreads();
	
}