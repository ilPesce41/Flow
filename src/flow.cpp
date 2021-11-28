#include "flow.hpp"
#include <string>
#include "utils.hpp"
using namespace std;
FlowData::FlowData(const char* filename, string sx_col, string sy_col, string dx_col, string dy_col, string length_col)
{
    parser = new CSVParser(filename);
    length = parser->length;

    //Allocate memory for flow data
    sx = (double *)malloc(sizeof(double)*length);
    sy = (double *)malloc(sizeof(double)*length);
    dx = (double *)malloc(sizeof(double)*length);
    dy = (double *)malloc(sizeof(double)*length);
    L  = (double *)malloc(sizeof(double)*length);

    convert_vec_double(parser->get_column_data(sx_col),sx);
    convert_vec_double(parser->get_column_data(sy_col),sy);
    convert_vec_double(parser->get_column_data(dx_col),dx);
    convert_vec_double(parser->get_column_data(dy_col),dy);
    convert_vec_double(parser->get_column_data(length_col),L);


}
