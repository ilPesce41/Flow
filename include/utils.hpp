/*
Cole Hill
University of South Florida
Computer Science and Engineering
Programming Massively Parallel Systems
Fall 2021
*/
#pragma once
# include <vector>
# include <string>
using namespace std;

void convert_vec_float(vector<string> vec, float * output);
float get_min(float * arr,int length);
float get_max(float * arr,int length);
float get_min(float * arr,int start, int stop);
float get_max(float * arr,int start, int stop);
void mask_array(float * input_array, float * output_array, int * mask_array, int length);
