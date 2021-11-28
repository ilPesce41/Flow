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
        double * sx;
        double * sy;
        double * dx;
        double * dy;
        double * L;
};