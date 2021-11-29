#pragma once
# include <vector>
# include <string>
using namespace std;

void convert_vec_float(vector<string> vec, float * output);
float get_min(float * arr,int length);
float get_max(float * arr,int length);
float get_min(float * arr,int start, int stop);
float get_max(float * arr,int start, int stop);