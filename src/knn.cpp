#include "knn.hpp"
#include <iostream>

using namespace std;


KNNResult::KNNResult(int k)
{
    distances = (float *)malloc(k*sizeof(float));
    neighbors = (int *)malloc(k*sizeof(int));
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
Naiive KNN on CPU using distance matrix from GPU
*/
KNNResult get_k_nearest_neighbors(int k,int flow_idx,int length,float * distance_matrix)
{
    // Initialize neighbor list with big values
    int neighbors[k+1];
    float distances[k+1];
    for(int i=0;i<k+1;i++)
    {
        neighbors[i] = 1000;
        distances[i] = 1000;
    }
    
    // Iterate through distance matrix and find nearest neighbors
    for(int i=0;i<length;i++)
    {
        if(distance_matrix[flow_idx*length+i]<distances[k-1] && distance_matrix[flow_idx*length+i]>0) //ignore self
        {
            distances[k] = distance_matrix[flow_idx*length+i];
            neighbors[k] = i;
            insertion_sort(distances,neighbors,k+1);
        }
    }

    // Return KNN result
    KNNResult result(k);
    for (int i=0;i<k;i++)
    {
        result.distances[i] = distances[i];
        result.neighbors[i] = neighbors[i];
    }
    return result;
}

