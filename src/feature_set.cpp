/*
Cole Hill
University of South Florida
Computer Science and Engineering
Programming Massively Parallel Systems
Fall 2021
*/
#include "feature_set.hpp"
#include <iostream>

/*
Initializer function
*/
FeatureSet::FeatureSet(int * features_, int k_)
{
    k = k_;
    for(int i=0;i<k;i++)
    {
        features[i] = features_[i];
    }
}

/*
Compares if two sets are equivalent
Innefficient in its current implementations (could use hashing)
*/
bool FeatureSet::is_equivalent(FeatureSet set)
{
    //Iterate through features in set
    for(int i=0;i<k;i++)
    {
        bool in_list=false;
        int feat = features[i];
        //query if feature in other set
        for(int j=0;j<k;j++)
        {
            if(set.features[j]==feat)
            {
                //feature was in list 
                //don't need to check rest
                in_list = true;
                break;
            }
        }
        //Feature wasn't in set, sets are different
        if(!in_list){return false;}
    }
    //All checks passed, sets are the same
    return true;
}