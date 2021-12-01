#include "colocation.cuh"
#include "spatial_distance.cuh"
#include "flow.hpp"
#include <vector>
#include <iostream>
#include "utils.hpp"
#include "feature_set.hpp"
#include <algorithm>

using namespace std;

#define SUM_BLK_SIZE 128

/*
Take intersection of two sets of neighbors
Excludes neighbors who are in the excluded class list
*/
void merge_neighbors(int * combined_list, int * list_1, int * list_2, int * class_lookup, int max_degree,int class_filter)
{
    int index = 0;
    for(int i=0;i<max_degree;i++)
    {
        int candidate = list_1[i];
        //End of list
        if(candidate==-1){break;}
        if(class_lookup[candidate]!=class_filter)
        {
            bool mutual = false;
            for(int j=0;j<max_degree;j++)
            {
                int query = list_2[j];
                if(query==-1){break;}
                if(query==candidate){mutual=true;break;}
            }
            if(mutual)
            {
                combined_list[index] = candidate;
                index++;
            }
        }
    }
}

int build_fclp(int * members,int * features, int * neighbors,int * members_,int * features_,int * neighbors_,int * adj_list, int k, int pair_count, int max_degree,int * class_lookup)
{
    int index = 0;
    for(int i=0;i<pair_count;i++)
    {
        for(int j=0;j<max_degree;j++)
        {
            int neighbor = neighbors[i*max_degree + j];
            //No neighbor case
            if(neighbor==-1)
            {
                break;
            }
            
            for(int w=0;w<k-1;w++)
            {
                members_[index*k+w] = members[i*(k-1)+w];
                features_[index*k+w] = features[i*(k-1)+w];
            }
            members_[index*k+k-1] = neighbor;
            features_[index*k+k-1] = class_lookup[neighbor];
            merge_neighbors(neighbors_+index*max_degree,neighbors+i*max_degree,adj_list+neighbor*max_degree,class_lookup,max_degree,class_lookup[neighbor]);
            index++;      
        }
    }
    return index;
}

/*
Function to remove patterns with a FCI below specified threshold
*/
int purge_fclp(int *members,int *features,int *neighbors,float frequency_threshold,int k,int max_degree,int table_length, int * class_frequency)
{
    vector<FeatureSet> feature_patterns;
    vector<int> feature_patterns_count;
    int feature_index[table_length];
    vector<int> bad_features;

    for(int i=0;i<table_length;i++)
    {
        FeatureSet set = FeatureSet(features+i*k,k);
        bool in_set = false;
        for(int j=0;j<feature_patterns.size();j++){
            if (feature_patterns[j].is_equivalent(set))
            {
                in_set = true;
                feature_patterns_count[j] = feature_patterns_count[j] + 1;
                feature_index[i] = j;
            }
        }
        if(!in_set)
        {
            feature_patterns.push_back(set);
            feature_patterns_count.push_back(0);
            feature_index[i] = feature_patterns.size()-1;
        }
    }

    for(int i=0;i<feature_patterns.size();i++)
    {
        FeatureSet set = feature_patterns[i];
        float fci = (float)feature_patterns_count[i]/class_frequency[set.features[0]];
        for(int j=1;j<k;j++)
        {
            fci = min((float)feature_patterns_count[i]/class_frequency[set.features[j]],fci);
        }
        cout << fci << endl;
        if(fci<frequency_threshold)
        {
            bad_features.push_back(i);
        }
    }

    int idx = 0;
    for(int i=0;i<table_length;i++)
    {
        if (find(bad_features.begin(), bad_features.end(), feature_index[i]) != bad_features.end()) {
            //do nothing
        }
        else {
            //copy row to first open spot and increment table length
            for(int w=0;w<k;w++)
            {
                members[idx*k+w] = members[i*k+w];
                features[idx*k+w] = features[i*k+w];
            }
            for(int w=0;w<max_degree;w++)
            {
                neighbors[idx*max_degree+w] = neighbors[i*max_degree+w];
            }
            idx++;
        }
    }
    return idx;
}


/*
Function to get the size of the k+1 table
*/
int get_max_pair_count(int * neighbors,int length,int max_degree)
{
    int pair_count = 0;
    for(int i=0;i<length;i++)
    {
        for(int j=0;j<max_degree;j++)
        {
            int neighbor = neighbors[max_degree*i + j];
            if(neighbor>-1)
            {
                pair_count++;
            }
        }
    }
    return pair_count;
}

/*
Given a threshold value converts every distance in matrix 
to
 0 - not neighbors
 1 - neighbors
Essentially creates incidence matrix
*/
__global__ void threshold_distance_matrix(float * dist_matrix,float threshold,int length){
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if(i>=length*length)
    {
        return;
    }
    if(dist_matrix[i]>threshold)
    {
        dist_matrix[i] = 0;
    }
    else
    {
        dist_matrix[i] = 1;
    }
}

__global__ void get_neighbors(float * dist_matrix,int * neighbor_table, int max_degree ,int length){
	
    // Index number
    int i = blockDim.x*blockIdx.x + threadIdx.x;
    if(i>=length)
    {
        return;
    }
    int row_index = i * length;
    int idx = 0;
    for(int j=0;j<length;j++)
    {
        if(dist_matrix[row_index + j]>0 && j!=i)
        {
            neighbor_table[i*max_degree + idx] = j;
            idx++;
        }
    }
    __syncthreads();
}

__global__ void sum_rows(float * matrix,int * sum_arr,int length)
{
    __shared__ float partial_sum[2*SUM_BLK_SIZE];

    int row = (int)((blockIdx.x*blockDim.x)/length);

    unsigned int t = threadIdx.x;
    unsigned int start = 2*blockIdx.x*blockDim.x;

    if(blockIdx.x*blockDim.x+t >=length*length)
    {
        return;
    }
    if(start + t >=(row+1)*length)
    {
        partial_sum[t] = 0;
    }
    else
    {
        partial_sum[t] = matrix[start + t];
    }

    if(blockDim.x+t >=(row+1)*length)
    {
        partial_sum[blockDim.x+t] = 0;
    }
    else
    {
        partial_sum[blockDim.x+t] = matrix[start + blockDim.x+t];
    }
    for (unsigned int stride = blockDim.x; 
        stride > 0;  stride /= 2) 
    {
        __syncthreads();
        if (t < stride)
        partial_sum[t] += partial_sum[t+stride];
    }
    __syncthreads();

    if(t==0)
    {
        atomicAdd(sum_arr+row,(int)partial_sum[0]);
    }
    __syncthreads();

}

ColocationResult colocate(vector<FlowData> flows,float frequency_threshold, float spatial_threshold, size_t shared_mem_size)
{

    //Ensure we have a fresh device
    cudaDeviceReset();
    
    //Determine total number of points
    int length = 0;
    for(int i=0;i<flows.size();i++)
    {
        length += flows[i].length;
    }

    // copy data to GPU 
	float *d_sx, *d_sy, *d_dx, *d_dy, *d_L;
	cudaMalloc((void**)&d_sx, length*sizeof(float));
	cudaMalloc((void**)&d_sy, length*sizeof(float));
	cudaMalloc((void**)&d_dx, length*sizeof(float));
	cudaMalloc((void**)&d_dy, length*sizeof(float));
	cudaMalloc((void**)&d_L, length*sizeof(float));
    int * class_lookup = (int *)malloc(length*sizeof(int));
    int * index_lookup = (int *)malloc(length*sizeof(int));
    int * class_frequency = (int * )malloc(flows.size()*sizeof(int));
    for(int i=0;i<flows.size();i++){class_lookup[i] = 0;}

    int start_idx = 0;
    for(int i=0;i<flows.size();i++)
    { 
        FlowData flow = flows[i];
        cudaMemcpy(d_sx + start_idx, flow.sx, flow.length*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_sy + start_idx, flow.sy, flow.length*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_dx + start_idx, flow.dx, flow.length*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_dy + start_idx, flow.dy, flow.length*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_L + start_idx, flow.L, flow.length*sizeof(float), cudaMemcpyHostToDevice);
        // printf("GPUassert: %s \n", cudaGetErrorString(cudaGetLastError()));

        for(int j=0;j<flow.length;j++)
        {
            class_lookup[start_idx+j] = i;
            index_lookup[start_idx+j] = j;
            class_frequency[i]++;
        }
        start_idx += flow.length;
    }

    float * dist_matrix_gpu;
	cudaMalloc((void**)&dist_matrix_gpu, length*length*sizeof(float));
	cudaMemset(dist_matrix_gpu, 0, length*length*sizeof(float));

    calculate_spatial_distance_matrix <<< ceil((float)length/128) , 128, shared_mem_size >>> (d_sx,d_sy,d_dx,d_dy,d_L,dist_matrix_gpu,length,2,1.0f);
    threshold_distance_matrix <<< ceil((float)length*length/128) , 128 >>>(dist_matrix_gpu,spatial_threshold,length);
    cudaFree((void**)&d_sx);
    cudaFree((void**)&d_sy);
    cudaFree((void**)&d_dx);
    cudaFree((void**)&d_sx);
    cudaFree((void**)&d_dy);
    
	
    // Get maximum degree of flow neighbor graph
    int * number_neighbors;
    cudaMalloc((void**)&number_neighbors, length*sizeof(int));
	cudaMemset(number_neighbors, 0, length*sizeof(int));

    // sum_rows <<< ceil((float)(length*length)/(2*128)) , 128 >>>(dist_matrix_gpu,number_neighbors,length);

    float * dist_matrix_cpu = (float*)malloc(length*length*sizeof(float));
	cudaMemcpy(dist_matrix_cpu, dist_matrix_gpu, length*length*sizeof(float), cudaMemcpyDeviceToHost);

    int max_degree = 0;
    long int total = 0;
    for(int i=0;i<length;i++)
    {
        int tmp = 0;
        for(int j=0;j<length;j++)
            tmp += dist_matrix_cpu[i*length + j];
            total++;
        if(tmp>max_degree)
            max_degree = tmp;
    }
    
    /*
    Initialize adjacency list structure
    */
    int *adj_list_cpu;
    // cudaMalloc((void**)&adj_list_gpu, length*max_degree*sizeof(int));
	// cudaMemset(number_neighbors, -1, length*max_degree*sizeof(int));
    adj_list_cpu = (int*)malloc(length*max_degree*sizeof(int));
    // get_neighbors<<< ceil((float)length/128) , 128 >>>(dist_matrix_gpu,adj_list_gpu, max_degree ,length);
    for(int i=0;i<length*max_degree;i++)
        {adj_list_cpu[i] = -1;}
    for(int i=0;i<length;i++)
    {
        int row = i*length;
        int index = 0;
        for(int j=0;j<max_degree;j++)
        {
            if(dist_matrix_cpu[row+j]>0)
            {
                if(class_lookup[j]!=class_lookup[i])
                {
                    adj_list_cpu[i*max_degree + index] = j;
                    index++;
                }
            }
        }
    }

    // cudaFree(dist_matrix_gpu);
	// cudaMemcpy(adj_list_cpu, adj_list_gpu, length*max_degree*sizeof(int), cudaMemcpyDeviceToHost);
    
    /* Initialize k=1 table
    */
    int * members_, * features_, *neighbors_;
    int *members, *features, *neighbors;

    members = (int *)malloc(length*sizeof(int));
    features = (int *)malloc(length*sizeof(int));
    neighbors = (int *)malloc(max_degree*length*sizeof(int));

    for(int i=0;i<max_degree*length;i++){neighbors[i] = -1;}
    int pair_count = 0;
    
    for(int i=0;i<length;i++)
    {
        int row = i*max_degree;
        int index = 0;
        members[i] = i;
        features[i] = class_lookup[i];
        for(int j=0;j<max_degree;j++)
        {
            if(adj_list_cpu[row + j]>i)
            {
                pair_count++;
                neighbors[i*max_degree+index] = adj_list_cpu[row+j];
                index++;
            }
        }
    }


    /*
    Iteratively build FCLP tables
    */
    int pair_count_old=length;
    int k_val=1;
    for(int k=2;k<flows.size()+1;k++)
    {
        cout << k << " " << flows.size() << endl;
        if(k>2)
        {
            free(members_);
            free(features_);
            free(neighbors_);
        }
        members_ = members;
        features_ = features;
        neighbors_ = neighbors;
        members = (int *)malloc(k*pair_count*sizeof(int));
        features = (int *)malloc(k*pair_count*sizeof(int));
        neighbors = (int *)malloc(max_degree*pair_count*sizeof(int));
        cout << max_degree*pair_count*sizeof(int) << endl;
        for(int i=0;i<max_degree*pair_count;i++){neighbors[i]=-1;}
        int table_length = build_fclp(members_,features_,neighbors_,members,features,neighbors,adj_list_cpu,k,pair_count_old,max_degree,class_lookup);
        table_length = purge_fclp(members,features,neighbors,frequency_threshold,k,max_degree,table_length,class_frequency);
        pair_count_old=pair_count;
        pair_count = get_max_pair_count(neighbors,table_length,max_degree);
        if(pair_count==0)
        {
            if(table_length>0)
            {
                pair_count_old = table_length;
                free(members_);
                free(features_);
                free(neighbors_);
                members_ = members;
                features_ = features;
                neighbors_ = neighbors;
            }
            break;
        }else{k_val++;}
    }
    ColocationResult result(k_val);
    result.indices = members_;
    result.class_lookup = class_lookup;
    result.index_lookup = index_lookup;
    result.length = k_val*pair_count_old;
    result.flow_length = length;
    return result;
}