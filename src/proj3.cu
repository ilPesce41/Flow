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


    /*
    Load Data
    */
    FlowData flow(filename, sx_col, sy_col, dx_col, dy_col, length_col);
    cout << flow.length <<endl;

    // copy data to GPU 
	float *d_sx, *d_sy, *d_dx, *d_dy, *d_L;
	cudaMalloc((void**)&d_sx, flow.length*sizeof(double));
	cudaMalloc((void**)&d_sy, flow.length*sizeof(double));
	cudaMalloc((void**)&d_dx, flow.length*sizeof(double));
	cudaMalloc((void**)&d_dy, flow.length*sizeof(double));
	cudaMalloc((void**)&d_L, flow.length*sizeof(double));
	cudaMemcpy(d_sx, flow.sx, flow.length*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_sy, flow.sy, flow.length*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_dx, flow.dx, flow.length*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_dy, flow.dy, flow.length*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_L, flow.L, flow.length*sizeof(double), cudaMemcpyHostToDevice);

	// create output space
	// BUCKET_TYPE *d_hist;
	// cudaMalloc((void**)&d_hist, num_buckets*sizeof(BUCKET_TYPE));
	// cudaMemset(d_hist, 0, num_buckets*sizeof(BUCKET_TYPE));

    
    /*
    Process Data
    */


    /*
    Write Output
    */

}