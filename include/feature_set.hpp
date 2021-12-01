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