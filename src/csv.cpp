#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "csv.hpp"
#include <algorithm>

using namespace std;

/*
Initializer function
*/
CSVParser::CSVParser(string filename_) {
    filename = filename_;
    columns = get_columns();
    length = get_column_data(columns[0]).size();
}

/*
Split line into a vector sting based on ',' delimiter
*/
vector<string> CSVParser::parse_line(string line){
    
    vector<string> data;
    stringstream stream(line);
    while(stream.good())
    {
        string token;
        getline(stream,token,',');
        data.push_back(token);
    }
    return data;
}

/*
Test is data file is accessible and gets the columns of 
the .csv file
*/
vector<string> CSVParser::get_columns(){
    ifstream input_stream(filename,ifstream::in);
    if (input_stream.is_open())
    {
        string column_line;
        getline(input_stream,column_line,'\n');
        return CSVParser::parse_line(column_line);
    }
    else
    {
        cout<< "Failed to open file!"<<endl;
        cout<< "Exiting"<< endl;
        exit(EXIT_FAILURE);
    }
    input_stream.close();
}  

/*
Returns all of the values in a column of a csv file as a vector<string>
Exits program if column not in data
*/
vector<string> CSVParser::get_column_data(string column_name){

    //Output data structure
    vector<string> column_data;
    
    //Check if column in csv and get its index
    auto col_pointer = find(columns.begin(), columns.end(), column_name);
    
    //Column in csv
    if (col_pointer != columns.end()) {
        int index = col_pointer - columns.begin();

        //Grab all of the values by iterating through every line
        ifstream input_stream(filename,ifstream::in);
        if (input_stream.is_open())
        {
            string column_line;
            while(getline(input_stream,column_line,'\n'))
            {
                column_data.push_back(CSVParser::parse_line(column_line)[index]);
            }
        }
        //IO issue, can't open file
        else
        {
            cout<< "Failed to open file!"<<endl;
            cout<< "Exiting"<< endl;
            exit(EXIT_FAILURE);
        }

    }
    // Column not in csv
    else {
        cout << "Element not found" << endl;
        cout << "Exiting" << endl;
        exit(EXIT_FAILURE);
    }
    column_data.erase(column_data.begin());
    return column_data;
} 
