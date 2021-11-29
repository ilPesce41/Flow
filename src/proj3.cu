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
#include "spatial_distance.cuh"
#include "knn.hpp"
#include "flow_k.cuh"

using namespace std;

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
    const char * filename_1 = "/home/c/coleh/COP6527_project_3/data/Taxi/green_tripdata_2015-06_subset.csv";
    const char * filename_2 = "/home/c/coleh/COP6527_project_3/data/Taxi/yellow_tripdata_2015-06_subset.csv";
    string sx_col("Pickup_latitude");
    string sy_col("Pickup_longitude");
    string dx_col("Dropoff_latitude");
    string dy_col("Dropoff_longitude");
    string length_col("Trip_distance");
    int block_dim = 128;
    float alpha = 1.0; //weighting between origin and destination 1.0=even
    int func_type = 0; // 0=flow_distance, 1=flow_dissimilarity
    size_t shared_mem_size = 10*block_dim*sizeof(float);

    FlowData flow_1(filename_1, sx_col, sy_col, dx_col, dy_col, length_col);
    FlowData flow_2(filename_2, sx_col, sy_col, dx_col, dy_col, length_col);


    /*
    Load Data
    */
    // FlowData flow(filename, sx_col, sy_col, dx_col, dy_col, length_col);
    // cout << "Loaded Data:" << flow.length <<endl;

    // // copy data to GPU 
	// float *d_sx, *d_sy, *d_dx, *d_dy, *d_L;
	// cudaMalloc((void**)&d_sx, flow.length*sizeof(float));
	// cudaMalloc((void**)&d_sy, flow.length*sizeof(float));
	// cudaMalloc((void**)&d_dx, flow.length*sizeof(float));
	// cudaMalloc((void**)&d_dy, flow.length*sizeof(float));
	// cudaMalloc((void**)&d_L, flow.length*sizeof(float));
	// cudaMemcpy(d_sx, flow.sx, flow.length*sizeof(float), cudaMemcpyHostToDevice);
	// cudaMemcpy(d_sy, flow.sy, flow.length*sizeof(float), cudaMemcpyHostToDevice);
	// cudaMemcpy(d_dx, flow.dx, flow.length*sizeof(float), cudaMemcpyHostToDevice);
	// cudaMemcpy(d_dy, flow.dy, flow.length*sizeof(float), cudaMemcpyHostToDevice);
	// cudaMemcpy(d_L, flow.L, flow.length*sizeof(float), cudaMemcpyHostToDevice);
    // cout << "Data on GPU" <<endl;

    // float * tmp_arr;
    // tmp_arr = (float *)malloc(flow.length*sizeof(float));
	// cudaMemcpy(tmp_arr, d_sx, flow.length*sizeof(float), cudaMemcpyDeviceToHost);
    // for(int i=0;i<10;i++)
    // {
    //     std::cout << tmp_arr[i] << " ";
    // }
    // cout << endl;
    // for(int i=0;i<10;i++)
    // {
    //     std::cout << flow.sx[i] << " ";
    // }


	// // create output space
	// float * dist_matrix_gpu;
	// cudaMalloc((void**)&dist_matrix_gpu, flow.length*flow.length*sizeof(float));
	// cudaMemset(dist_matrix_gpu, 0, flow.length*flow.length*sizeof(float));

    // cout << "Starting Kernel" <<endl;
    // std::cout<< shared_mem_size << std::endl;
    // calculate_spatial_distance_matrix <<< ceil((float)flow.length/block_dim) , block_dim, shared_mem_size >>> (d_sx,d_sy,d_dx,d_dy,d_L,dist_matrix_gpu,flow.length,func_type,alpha);
    // printf("GPUassert: %s \n", cudaGetErrorString(cudaGetLastError()));
    // cout << "Kernel Done" <<endl;

    // float * dist_matrix_cpu = (float*)malloc(flow.length*flow.length*sizeof(float));
	// cudaMemcpy(dist_matrix_cpu, dist_matrix_gpu, flow.length*flow.length*sizeof(float), cudaMemcpyDeviceToHost);
    // for(int i=0;i<3;i++)
    // {
    //     for(int j=0;j<flow.length;j++)
    //     {
    //         cout << dist_matrix_cpu[i*flow.length+j] << " ";
    //     }
    //     cout << endl;
    // }

    // NN
    // get_k_nearest_neighbors(5,0,flow.length,dist_matrix_cpu);
    
    //flow-k
    // vector<int> clusters = process_flow_k(flow, d_sx, d_sy, d_dx, d_dy, d_L, dist_matrix_cpu, dist_matrix_gpu, 100,func_type,alpha,shared_mem_size,0.17);
    // cout << clusters.size() << endl;

    vector<int> clusters = process_cross_flow_k(flow_1, flow_2, 1000,func_type,alpha,shared_mem_size,0.15);
    cout << clusters.size() << endl;
    for(int i=0;i<clusters.size();i++)
    {
        cout << clusters[i] << endl;
    }
}