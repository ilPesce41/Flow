/*
Cole Hill
University of South Florida
Computer Science and Engineering
Programming Massively Parallel Systems
Fall 2021
*/
#include "csv.hpp"
#include <string>
#pragma once

using namespace std;
/*
Class for loading and manipulating flow data from
a csv file
*/
class FlowData{
    private:
        CSVParser * parser;
        float * sx_;
        float * sy_;
        float * dx_;
        float * dy_;
        float * L_;
        int * mask_;
        long int length_;

    public:
        FlowData(string filename, string sx_col, string sy_col, string dx_col, string dy_col, string length_col);
        long int length;
        float * sx;
        float * sy;
        float * dx;
        float * dy;
        float * L;
        float min_x;
        float max_x;
        float min_y;
        float max_y;
        float area;
        void update_bounds();
        void apply_mask();
        void clear_mask();
        void filter(string column,float lb, float ub);
};