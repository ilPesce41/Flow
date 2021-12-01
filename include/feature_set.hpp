/*
Cole Hill
University of South Florida
Computer Science and Engineering
Programming Massively Parallel Systems
Fall 2021
*/
#pragma once

/*
Class for holding feature sets and comparing them
*/
class FeatureSet{
    public:
    int k;
    int features[100];
    FeatureSet(int * features_, int k_);
    bool is_equivalent(FeatureSet set);
};