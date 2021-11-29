# include "utils.hpp"
# include <vector>
# include <string>
# include <iostream>

using namespace std;

void convert_vec_float(vector<string> vec, float * output)
{
    for(int i=0;i<vec.size();i++)
    {
        output[i] = stod(vec[i].c_str());
    }
}

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