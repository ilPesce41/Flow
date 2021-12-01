#include <vector>
#include <string>
#pragma once
using namespace std;
/*
Class for Parsing Data from CSV Files
*/
class CSVParser
{
    private:
        string filename;
        vector<string> parse_line(string line);
        vector<string> columns;

    public:
        CSVParser(string filename);
        vector<string> get_columns();
        vector<string> get_column_data(string column_name);
        long int length;
};