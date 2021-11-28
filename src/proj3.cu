#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdint.h>
#include "csv.hpp"
#include <vector>
#include <string>
#include <iostream>
#include "flow.hpp"

using namespace std;



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
				    dist = flow_distance(1.0f,x1,y1,u1,v1,R_x[t],R_y[t],R_u[t],R_v[t]);
                }
                else
                {
				    dist = flow_disimilarity(1.0f,x1,y1,u1,v1,L1,R_x[t],R_y[t],R_u[t],R_v[t],R_L[t]);
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
                else
                {
				    dist = flow_disimilarity(alpha,x1,y1,u1,v1,L1,B_x[t],B_y[t],B_u[t],B_v[t],B_L[t]);
                }
                dist_matrix[i*length + j] = dist;
                dist_matrix[j*length + i] = dist;
			}
		}
	__syncthreads();
	
}

/*
Insertion sort, keeps track of indices
*/
void insertion_sort(float * dist, int * indices,int len)
{
    for(int i=1;i<len;i++)
    {
        float dist_val = dist[i];
        int index_val = indices[i];
        for(int j=1;j<=i;j++)
        {
            if(dist[i-j]>dist_val)
            {
                dist[i-j+1] = dist[i-j];
                indices[i-j+1] = indices[i-j];
                dist[i-j] = dist_val;
                indices[i-j] = index_val;
            }
            else
            {
                break;
            }
        }
    }
}

/*
Naiive KNN on CPU using distance matrix
*/
void get_k_nearest_neighbors(int k,int flow_idx,int length,float * distance_matrix)
{
    int neighbors[k+1];
    float distances[k+1];
    for(int i=0;i<k+1;i++)
    {
        neighbors[i] = i;
        distances[i] = distance_matrix[i];
    }
    insertion_sort(distances,neighbors,k+1);
    for(int i=k;i<length;i++)
    {
        if(distance_matrix[i]<distances[k-1])
        {
            distances[k] = distance_matrix[i];
            neighbors[k] = i;
            insertion_sort(distances,neighbors,k+1);
        }
    }
    for (int i=0;i<k;i++)
    {
        cout << i << " Distance: " << distances[i] << " Index: " << neighbors[i]<< endl;
    }
}

int main(int argc, char **argv)
{

    /*
    Parse User Arguments
    */
    // CSVParser parser("/home/c/coleh/COP6527_project_3/data/Taxi/green_tripdata_2015-06.csv");
    // vector<string> columns = parser.get_columns();
    // for (string col: columns){
    //     cout<<col<<" ";
    // }
    const char * filename = "/home/c/coleh/COP6527_project_3/data/Taxi/green_tripdata_2015-06_subset.csv";
    string sx_col("Pickup_latitude");
    string sy_col("Pickup_longitude");
    string dx_col("Dropoff_latitude");
    string dy_col("Dropoff_longitude");
    string length_col("Trip_distance");
    int block_dim = 128;
    float alpha = 1.0; //weighting between origin and destination 1.0=even
    int func_type = 0; // 0=flow_distance, 1=flow_dissimilarity


    /*
    Load Data
    */
    FlowData flow(filename, sx_col, sy_col, dx_col, dy_col, length_col);
    cout << "Loaded Data:" << flow.length <<endl;

    // copy data to GPU 
	float *d_sx, *d_sy, *d_dx, *d_dy, *d_L;
	cudaMalloc((void**)&d_sx, flow.length*sizeof(float));
	cudaMalloc((void**)&d_sy, flow.length*sizeof(float));
	cudaMalloc((void**)&d_dx, flow.length*sizeof(float));
	cudaMalloc((void**)&d_dy, flow.length*sizeof(float));
	cudaMalloc((void**)&d_L, flow.length*sizeof(float));
	cudaMemcpy(d_sx, flow.sx, flow.length*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_sy, flow.sy, flow.length*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_dx, flow.dx, flow.length*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_dy, flow.dy, flow.length*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_L, flow.L, flow.length*sizeof(float), cudaMemcpyHostToDevice);
    cout << "Data on GPU" <<endl;

    size_t shared_mem_size = 10*block_dim*sizeof(float);

	// create output space
	float * dist_matrix_gpu;
	cudaMalloc((void**)&dist_matrix_gpu, flow.length*flow.length*sizeof(float));
	cudaMemset(dist_matrix_gpu, 0, flow.length*flow.length*sizeof(float));

    cout << "Starting Kernel" <<endl;
    calculate_spatial_distance_matrix <<< ceil((float)flow.length/block_dim) , block_dim, shared_mem_size >>> (d_sx,d_sy,d_dx,d_dy,d_L,dist_matrix_gpu,flow.length,alpha,func_type);
    cout << "Kernel Done" <<endl;

    float * dist_matrix_cpu = (float*)malloc(flow.length*flow.length*sizeof(float));
	cudaMemcpy(dist_matrix_cpu, dist_matrix_gpu, flow.length*flow.length*sizeof(float), cudaMemcpyDeviceToHost);
    // for(int i=0;i<3;i++)
    // {
    //     for(int j=0;j<flow.length;j++)
    //     {
    //         cout << dist_matrix_cpu[i*flow.length+j] << " ";
    //     }
    //     cout << endl;
    // }
    get_k_nearest_neighbors(5,0,flow.length,dist_matrix_cpu);

    /*
    Process Data
    */


    /*
    Write Output
    */

}