#include <string>
#include "flow.hpp"
#include "utils.hpp"

using namespace std;
FlowData::FlowData(const char* filename, string sx_col, string sy_col, string dx_col, string dy_col, string length_col)
{
    parser = new CSVParser(filename);
    length = parser->length;

    //Allocate memory for flow data
    sx = (float *)malloc(sizeof(float)*length);
    sy = (float *)malloc(sizeof(float)*length);
    dx = (float *)malloc(sizeof(float)*length);
    dy = (float *)malloc(sizeof(float)*length);
    L  = (float *)malloc(sizeof(float)*length);

    convert_vec_float(parser->get_column_data(sx_col),sx);
    convert_vec_float(parser->get_column_data(sy_col),sy);
    convert_vec_float(parser->get_column_data(dx_col),dx);
    convert_vec_float(parser->get_column_data(dy_col),dy);
    convert_vec_float(parser->get_column_data(length_col),L);

    float min_sx = get_min(sx, length);
    float min_sy = get_min(sy, length);
    float min_dx = get_min(dx, length);
    float min_dy = get_min(dy, length);
    float max_sx = get_max(sx, length);
    float max_sy = get_max(sy, length);
    float max_dx = get_max(dx, length);
    float max_dy = get_max(dy, length);

    min_x = std::min(min_sx,min_dx);
    min_y = std::min(min_sy,min_dy);
    max_x = std::max(max_sx,max_dx);
    max_y = std::max(max_sy,max_dy);
    area = (max_x-min_x)*(max_y-min_y);
}
