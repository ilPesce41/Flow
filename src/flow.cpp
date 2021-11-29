#include <string>
#include "flow.hpp"
#include "utils.hpp"
#include <iostream>
using namespace std;
FlowData::FlowData(const char* filename, string sx_col, string sy_col, string dx_col, string dy_col, string length_col)
{
    parser = new CSVParser(filename);
    length_ = parser->length;
    //Allocate memory for flow data
    sx_ = (float *)malloc(sizeof(float)*length_);
    sy_ = (float *)malloc(sizeof(float)*length_);
    dx_ = (float *)malloc(sizeof(float)*length_);
    dy_ = (float *)malloc(sizeof(float)*length_);
    L_  = (float *)malloc(sizeof(float)*length_);
    sx = (float *)malloc(sizeof(float)*length_);
    sy = (float *)malloc(sizeof(float)*length_);
    dx = (float *)malloc(sizeof(float)*length_);
    dy = (float *)malloc(sizeof(float)*length_);
    L  = (float *)malloc(sizeof(float)*length_);
    mask_ = (int *)malloc(sizeof(int)*length_);
    clear_mask();

    convert_vec_float(parser->get_column_data(sx_col),sx_);
    convert_vec_float(parser->get_column_data(sy_col),sy_);
    convert_vec_float(parser->get_column_data(dx_col),dx_);
    convert_vec_float(parser->get_column_data(dy_col),dy_);
    convert_vec_float(parser->get_column_data(length_col),L_);
    
    apply_mask();
    
}

void FlowData::update_bounds()
{
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

void FlowData::apply_mask()
{
    length=0;
    for(int i=0;i<length_;i++)
    {   
        length += mask_[i];
    }
    //free visible arrays
    free(sx);
    free(sy);
    free(dx);
    free(dy);
    free(L);
    sx = (float *)malloc(sizeof(float)*length);
    sy = (float *)malloc(sizeof(float)*length);
    dx = (float *)malloc(sizeof(float)*length);
    dy = (float *)malloc(sizeof(float)*length);
    L  = (float *)malloc(sizeof(float)*length);

    mask_array(sx_, sx, mask_, length_);
    mask_array(sy_, sy, mask_, length_);
    mask_array(dx_, dx, mask_, length_);
    mask_array(dy_, dy, mask_, length_);
    mask_array(L_, L, mask_, length_);
    update_bounds();
}

void FlowData::clear_mask()
{
    for(int i=0;i<length_;i++)
    {
        mask_[i] = 1;
    }
}

void FlowData::filter(string column,float lb, float ub)
{

    float *tmp;
    tmp = (float *)malloc(sizeof(float)*length_);
    convert_vec_float(parser->get_column_data(column),tmp);
    for(int i=0;i<length_;i++)
    {
        if(tmp[i]<lb || tmp[i]>ub)
        {
            mask_[i] = 0;
        }
    }
}
