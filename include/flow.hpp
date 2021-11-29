#include "csv.hpp"
#include <string>
#pragma once

using namespace std;
class FlowData{
    private:
        CSVParser * parser;

    public:
        FlowData(const char* filename, string sx_col, string sy_col, string dx_col, string dy_col, string length_col);
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
};