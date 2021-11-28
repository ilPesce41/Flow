# include "utils.hpp"
# include <vector>
# include <string>
# include <iostream>
using namespace std;
void convert_vec_double(vector<string> vec, double * output)
{
    for(int i=0;i<vec.size();i++)
    {
        output[i] = stod(vec[i].c_str());
    }
}
