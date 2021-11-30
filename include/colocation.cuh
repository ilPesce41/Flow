#pragma once
#include <vector>
#include "flow.hpp"
using namespace std;

void merge_neighbors(int * combined_list, int * list_1, int * list_2, int * class_lookup, int max_degree,int class_filter);
int build_fclp(int * members_,int * features_,int * neighbors_,int * members,int * features,int * neighbors,int * adj_list, int k, int pair_count, int max_degree,int * class_lookup);
int purge_fclp(int *members,int *features,int *neighbors,float frequency_threshold,int k,int max_degree,int table_length, int * class_frequency);
int get_max_pair_count(int * neighbors,int length,int max_degree);
__global__ void threshold_distance_matrix(float * dist_matrix,float threshold,int length);
__global__ void get_neighbors(float * dist_matrix,int * neighbor_table, int max_degree ,int length);
__global__ void sum_rows(float * matrix,int * sum_arr,int length);
void colocate(vector<FlowData> flows,float frequency_threshold, float spatial_threshold, size_t shared_mem_size);