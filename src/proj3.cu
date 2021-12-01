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
#include "colocation.cuh"
#include <fstream>
#include<sstream>

#define BLOCK_DIM 128

using namespace std;

void knn(FlowData flow, int k, int flow_idx,int func_type, float alpha,string filename)
{
    cudaDeviceReset();
    float *d_sx, *d_sy, *d_dx, *d_dy, *d_L;
    size_t shared_mem_size = 10*BLOCK_DIM*sizeof(float);

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
    
    float * dist_matrix_gpu;
	cudaMalloc((void**)&dist_matrix_gpu, flow.length*flow.length*sizeof(float));
	cudaMemset(dist_matrix_gpu, 0, flow.length*flow.length*sizeof(float));

    calculate_spatial_distance_matrix <<< ceil((float)flow.length/BLOCK_DIM) , BLOCK_DIM, shared_mem_size >>> (d_sx,d_sy,d_dx,d_dy,d_L,dist_matrix_gpu,flow.length,func_type,alpha);

    float * dist_matrix_cpu = (float*)malloc(flow.length*flow.length*sizeof(float));
	cudaMemcpy(dist_matrix_cpu, dist_matrix_gpu, flow.length*flow.length*sizeof(float), cudaMemcpyDeviceToHost);
    get_k_nearest_neighbors(k, flow_idx, flow.length, dist_matrix_cpu);
    cudaDeviceReset();
}

void flow_k(FlowData flow, int num_iter, int func_type,float alpha,float radius,string filename)
{
    cudaDeviceReset();
    float *d_sx, *d_sy, *d_dx, *d_dy, *d_L;
    size_t shared_mem_size = 10*BLOCK_DIM*sizeof(float);
    cudaMalloc((void**)&d_sx, flow.length*sizeof(float));
	cudaMalloc((void**)&d_sy, flow.length*sizeof(float));
	cudaMalloc((void**)&d_dx, flow.length*sizeof(float));
	cudaMalloc((void**)&d_dy, flow.length*sizeof(float));
	cudaMalloc((void**)&d_L, flow.length*sizeof(float));

    float * dist_matrix_gpu;
	cudaMalloc((void**)&dist_matrix_gpu, flow.length*flow.length*sizeof(float));
	cudaMemset(dist_matrix_gpu, 0, flow.length*flow.length*sizeof(float));

    process_flow_k(flow, d_sx, d_sy, d_dx, d_dy, d_L, dist_matrix_gpu, num_iter, func_type, alpha, shared_mem_size, radius);
    cudaDeviceReset();
}

void cross_flow_k(FlowData flow_1, FlowData flow_2, int num_iter, int func_type,float alpha,float radius,string filename)
{
    cudaDeviceReset();
    size_t shared_mem_size = 10*BLOCK_DIM*sizeof(float);
    process_cross_flow_k(flow_1, flow_2, num_iter, func_type, alpha, shared_mem_size, radius);
    cudaDeviceReset();
}

void extract_colocation_patterns(vector<FlowData> flows,float frequency_threshold, float spatial_threshold,string filename)
{
    cudaDeviceReset();
    size_t shared_mem_size = 10*BLOCK_DIM*sizeof(float);
    colocate(flows,frequency_threshold, spatial_threshold,shared_mem_size);
    cudaDeviceReset();
}

void check_arguments(vector<string> args,int number,int line)
{
    if (args.size()<number)
    {
        cout << "Invalid number of arguments on line "<< line << endl;
        cout << "There are "<<args.size()<< " when there should be " << number << endl;
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char **argv)
{

    if(argc!=2)
    {
        cout << "Invalid Arguments. Exiting..." << endl;
        exit(1);
    }

    vector<FlowData> flows;
    ifstream config_file;
    config_file.open(argv[1]);
    if(config_file.fail())
    {
        cout << "File does not exist. Exiting..." << endl;
        exit(1);
    }

    string line;
    int line_number = 1;
    while(getline(config_file,line))
    {
        vector<string> tokens;
        stringstream linestream(line); 
        while(linestream.good())
        {
            string substr;
            getline(linestream, substr, ' ');
            tokens.push_back(substr);
        }
        if(tokens.size()==0){continue;line_number++;}
        if(tokens[0]=="load")
        {
            check_arguments(tokens,7,line_number);
            FlowData flow = FlowData(tokens[1],tokens[2],tokens[3],tokens[4],tokens[5],tokens[6]);
            flows.push_back(flow);
        }
        else if (tokens[0]=="filter")
        {
            check_arguments(tokens,5,line_number);
            int flow_index = stoi(tokens[1]);
            if(flow_index>tokens.size()){cout<<"Invalid index on line "<<line_number<<" exiting."<<endl;exit(EXIT_FAILURE);}
            flows[flow_index].filter(tokens[2],stod(tokens[3]),stod(tokens[4]));
        }
        else if (tokens[0]=="apply")
        {
            check_arguments(tokens,2,line_number);
            int flow_index = stoi(tokens[1]);
            if(flow_index>tokens.size()){cout<<"Invalid index on line "<<line_number<<" exiting."<<endl;exit(EXIT_FAILURE);}
            flows[flow_index].apply_mask();
        }
        else if (tokens[0]=="knn")
        {
            check_arguments(tokens,7,line_number);
            int flow_index = stoi(tokens[1]);
            if(flow_index>tokens.size()){cout<<"Invalid index on line "<<line_number<<" exiting."<<endl;exit(EXIT_FAILURE);}
            knn(flows[flow_index],stoi(tokens[2]),stoi(tokens[3]),stoi(tokens[4]),stod(tokens[5]),tokens[6]);
        }
        else if (tokens[0]=="flow_k")
        {
            check_arguments(tokens,7,line_number);
            int flow_index = stoi(tokens[1]);
            if(flow_index>tokens.size()){cout<<"Invalid index on line "<<line_number<<" exiting."<<endl;exit(EXIT_FAILURE);}
            flow_k(flows[flow_index],stoi(tokens[2]),stoi(tokens[3]),stod(tokens[4]),stod(tokens[5]),tokens[6]);
        }
        else if (tokens[0]=="cross_flow_k")
        {
            check_arguments(tokens,8,line_number);
            int flow_index_1 = stoi(tokens[1]);
            if(flow_index_1>tokens.size()){cout<<"Invalid index on line "<<line_number<<" exiting."<<endl;exit(EXIT_FAILURE);}
            int flow_index_2 = stoi(tokens[2]);
            if(flow_index_2>tokens.size()){cout<<"Invalid index on line "<<line_number<<" exiting."<<endl;exit(EXIT_FAILURE);}
            cross_flow_k(flows[flow_index_1],flows[flow_index_2],stoi(tokens[3]),stoi(tokens[4]),stod(tokens[5]),stod(tokens[6]),tokens[7]);            
        }
        else if (tokens[0]=="colocation")
        {
            check_arguments(tokens,4,line_number);
            if(flows.size()<2)
            {
                cout << "Insufficient number of flows loaded, need at least 2"<<endl;
                cout << "Exiting" << endl;
                exit(EXIT_FAILURE);
            }
            extract_colocation_patterns(flows, stod(tokens[1]),stod(tokens[2]),tokens[3]);
            
        }
        line_number++;
    }




    
    

}