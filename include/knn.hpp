#pragma once

void insertion_sort(float * dist, int * indices,int len);

/*
Class for passing k-nearest neighbors results
*/
class KNNResult{
    
    public:
    KNNResult(int k);
    int * neighbors;
    float * distances;
};

KNNResult get_k_nearest_neighbors(int k,int flow_idx,int length,float * distance_matrix);