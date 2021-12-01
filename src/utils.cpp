# include "utils.hpp"
# include <vector>
# include <string>
# include <iostream>

using namespace std;

/*
Convert vector of strings to a float array
*/
void convert_vec_float(vector<string> vec, float * output)
{
    for(int i=0;i<vec.size();i++)
    {
        output[i] = stod(vec[i].c_str());
    }
}

/*
Apply a mask to an array
*/
void mask_array(float * input_array, float * output_array, int * mask_array, int length)
{
    int idx = 0;
    for(int i=0;i<length;i++)
    {
        if(mask_array[i])
        {
            output_array[idx] = input_array[i];
            idx++;
        }
    }
}

/*
Get min of array
*/
float get_min(float * arr,int length)
{
    float min_value = arr[0];
    for(int i=0;i<length;i++)
    {
        if(arr[i]<min_value && abs(arr[i])>1e-3)
            min_value = arr[i];
    }
    return min_value;
}

/*
Get min of array
*/
float get_min(float * arr,int start, int stop)
{
    float min_value = arr[0];
    for(int i=start;i<stop;i++)
    {
        if(arr[i]<min_value && abs(arr[i])>1e-3)
            min_value = arr[i];
    }
    return min_value;
}

/*
Get max of array
*/
float get_max(float * arr,int length)
{
    float max_value = arr[0];
    for(int i=0;i<length;i++)
    {
        if(arr[i]>max_value && abs(arr[i])>1e-3)
            max_value = arr[i];
    }
    return max_value;
}

/*
Get max of array
*/
float get_max(float * arr,int start, int stop)
{
    float max_value = arr[0];
    for(int i=start;i<stop;i++)
    {
        if(arr[i]>max_value && abs(arr[i])>1e-3)
            max_value = arr[i];
    }
    return max_value;
}