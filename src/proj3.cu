/*
Cole Hill
University of South Florida
Computer Science and Engineering
Programming Massively Parallel Systems
Fall 2021
*/
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


/*
Find KNN of a flow in a dataset of flows
*/
void knn(FlowData flow, int k, int flow_idx,int func_type, float alpha,string filename)
{
    cudaDeviceReset();
    // Copy flow data to GPU
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
    
    // Calculate spatial distance matrix
    float * dist_matrix_gpu;
	cudaMalloc((void**)&dist_matrix_gpu, flow.length*flow.length*sizeof(float));
	cudaMemset(dist_matrix_gpu, 0, flow.length*flow.length*sizeof(float));
    calculate_spatial_distance_matrix <<< ceil((float)flow.length/BLOCK_DIM) , BLOCK_DIM, shared_mem_size >>> (d_sx,d_sy,d_dx,d_dy,d_L,dist_matrix_gpu,flow.length,func_type,alpha);

    // Copy distance matrix back to CPU
    float * dist_matrix_cpu = (float*)malloc(flow.length*flow.length*sizeof(float));
	cudaMemcpy(dist_matrix_cpu, dist_matrix_gpu, flow.length*flow.length*sizeof(float), cudaMemcpyDeviceToHost);
    // Get KNN
    KNNResult result = get_k_nearest_neighbors(k, flow_idx, flow.length, dist_matrix_cpu);

    // Output KNN to csv file
    ofstream output_file;
    output_file.open(filename);
    output_file << "sx,sy,dx,dy,L,distance\n";
    for(int i=0;i<k;i++)
    {
        output_file << flow.sx[result.neighbors[i]] << ",";
        output_file << flow.sy[result.neighbors[i]] << ",";
        output_file << flow.dx[result.neighbors[i]] << ",";
        output_file << flow.dy[result.neighbors[i]] << ",";
        output_file << flow.L[result.neighbors[i]] << ",";
        output_file << result.distances[i] << "\n";
    }
    cudaDeviceReset();
}

/*
Find flow clusters using ripleys k value
*/
void flow_k(FlowData flow, int num_iter, int func_type,float alpha,float radius,string filename)
{
    cudaDeviceReset();
    // Copy flow data to GPU
    float *d_sx, *d_sy, *d_dx, *d_dy, *d_L;
    size_t shared_mem_size = 10*BLOCK_DIM*sizeof(float);
    cudaMalloc((void**)&d_sx, flow.length*sizeof(float));
	cudaMalloc((void**)&d_sy, flow.length*sizeof(float));
	cudaMalloc((void**)&d_dx, flow.length*sizeof(float));
	cudaMalloc((void**)&d_dy, flow.length*sizeof(float));
	cudaMalloc((void**)&d_L, flow.length*sizeof(float));

    // Allocate space for spatial distance matrix on GPU
    float * dist_matrix_gpu;
	cudaMalloc((void**)&dist_matrix_gpu, flow.length*flow.length*sizeof(float));
	cudaMemset(dist_matrix_gpu, 0, flow.length*flow.length*sizeof(float));
    
    // Calcuate flow clusters
    vector<int> result = process_flow_k(flow, d_sx, d_sy, d_dx, d_dy, d_L, dist_matrix_gpu, num_iter, func_type, alpha, shared_mem_size, radius);
    
    // Output flow clusters to CSV file
    ofstream output_file;
    output_file.open(filename);
    output_file << "sx,sy,dx,dy,L\n";
    for(int i=0;i<result.size();i++)
    {
        output_file << flow.sx[result[i]] << ",";
        output_file << flow.sy[result[i]] << ",";
        output_file << flow.dx[result[i]] << ",";
        output_file << flow.dy[result[i]] << ",";
        output_file << flow.L[result[i]] << "\n";
    }
    cudaDeviceReset();
}

/*
Find cross flow clusters using ripleys k value
*/
void cross_flow_k(FlowData flow_1, FlowData flow_2, int num_iter, int func_type,float alpha,float radius,string filename)
{
    cudaDeviceReset();
    size_t shared_mem_size = 10*BLOCK_DIM*sizeof(float);
    // Calculate cross flow k values
    vector<int> result = process_cross_flow_k(flow_1, flow_2, num_iter, func_type, alpha, shared_mem_size, radius);
    
    // Output significant cross flows to csv file
    ofstream output_file;
    output_file.open(filename);
    output_file << "class,sx,sy,dx,dy,L\n";
    bool in_flow_1 = true;
    for(int i=0;i<result.size();i++)
    {
        int res = result[i];
        if(res==-1)
        {
            in_flow_1=false;
        }
        else
        {
            if(in_flow_1)
            {
                output_file << "0,";
                output_file << flow_1.sx[result[i]] << ",";
                output_file << flow_1.sy[result[i]] << ",";
                output_file << flow_1.dx[result[i]] << ",";
                output_file << flow_1.dy[result[i]] << ",";
                output_file << flow_1.L[result[i]] << "\n";
            }
            else
            {
                output_file << "1,";
                output_file << flow_2.sx[result[i]] << ",";
                output_file << flow_2.sy[result[i]] << ",";
                output_file << flow_2.dx[result[i]] << ",";
                output_file << flow_2.dy[result[i]] << ",";
                output_file << flow_2.L[result[i]] << "\n";
            }
        }
    }
    cudaDeviceReset();
}

/*
Extract multivariate flow colocation patterns
*/
void extract_colocation_patterns(vector<FlowData> flows,float frequency_threshold, float spatial_threshold,string filename)
{
    cudaDeviceReset();
    size_t shared_mem_size = 10*BLOCK_DIM*sizeof(float);
    ColocationResult result = colocate(flows,frequency_threshold, spatial_threshold,shared_mem_size);
    ofstream output_file;
    output_file.open(filename);
    output_file << "index,class,sx,sy,dx,dy,L,k\n";
    int length = sizeof(result.indices)/sizeof(int);
    for(int i=0;i<result.length;i++)
    {
        int index = result.indices[i];
        if(index<result.flow_length && index>=0)
        {
            int class_val = result.class_lookup[index];
            index = result.index_lookup[index];
            output_file << index << ",";
            output_file << class_val << ",";
            FlowData flow = flows[class_val];
            output_file << flow.sx[index] << ",";
            output_file << flow.sy[index] << ",";
            output_file << flow.dx[index] << ",";
            output_file << flow.dy[index] << ",";
            output_file << flow.L[index] << ",";
            output_file << result.k << "\n";
        }
    }
    cudaDeviceReset();
}

/*
Assert number of arguments in configuration file is correct
*/
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

    // Check configuraiton file is provided
    if(argc!=2)
    {
        cout << "Invalid Arguments. Exiting..." << endl;
        exit(1);
    }

    // Initialize flow vector
    vector<FlowData> flows;

    //Open configuration file
    ifstream config_file;
    config_file.open(argv[1]);
    if(config_file.fail())
    {
        cout << "File does not exist. Exiting..." << endl;
        exit(1);
    }

    // Iterate through configuration file line by line
    // Perform each command as supplied
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
        
        // load flow csv
        if(tokens[0]=="load")
        {
            check_arguments(tokens,7,line_number);
            FlowData flow = FlowData(tokens[1],tokens[2],tokens[3],tokens[4],tokens[5],tokens[6]);
            flows.push_back(flow);
        }
        // add filter to flow file
        else if (tokens[0]=="filter")
        {
            check_arguments(tokens,5,line_number);
            int flow_index = stoi(tokens[1]);
            if(flow_index>tokens.size()){cout<<"Invalid index on line "<<line_number<<" exiting."<<endl;exit(EXIT_FAILURE);}
            flows[flow_index].filter(tokens[2],stod(tokens[3]),stod(tokens[4]));
        }
        
        // activate filters on flow file
        else if (tokens[0]=="apply")
        {
            check_arguments(tokens,2,line_number);
            int flow_index = stoi(tokens[1]);
            if(flow_index>flows.size()){cout<<"Invalid index on line "<<line_number<<" exiting."<<endl;exit(EXIT_FAILURE);}
            flows[flow_index].apply_mask();
        }
        
        // find knn for a flow in flow file
        else if (tokens[0]=="knn")
        {
            check_arguments(tokens,7,line_number);
            int flow_index = stoi(tokens[1]);
            if(flow_index>flows.size()){cout<<"Invalid index on line "<<line_number<<" exiting."<<endl;exit(EXIT_FAILURE);}
            knn(flows[flow_index],stoi(tokens[2]),stoi(tokens[3]),stoi(tokens[4]),stod(tokens[5]),tokens[6]);
        }
        
        // find flow clusters in flow file
        else if (tokens[0]=="flow_k")
        {
            check_arguments(tokens,7,line_number);
            int flow_index = stoi(tokens[1]);
            if(flow_index>flows.size()){cout<<"Invalid index on line "<<line_number<<" exiting."<<endl;exit(EXIT_FAILURE);}
            flow_k(flows[flow_index],stoi(tokens[2]),stoi(tokens[3]),stod(tokens[4]),stod(tokens[5]),tokens[6]);
        }
        
        // find cross flow clusters in two flows
        else if (tokens[0]=="cross_flow_k")
        {
            check_arguments(tokens,8,line_number);
            int flow_index_1 = stoi(tokens[1]);
            if(flow_index_1>flows.size()){cout<<"Invalid index on line "<<line_number<<" exiting."<<endl;exit(EXIT_FAILURE);}
            int flow_index_2 = stoi(tokens[2]);
            if(flow_index_2>flows.size()){cout<<"Invalid index on line "<<line_number<<" exiting."<<endl;exit(EXIT_FAILURE);}
            cross_flow_k(flows[flow_index_1],flows[flow_index_2],stoi(tokens[3]),stoi(tokens[4]),stod(tokens[5]),stod(tokens[6]),tokens[7]);            
        }
        
        // find flow colocation patterns in all loaded flows
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