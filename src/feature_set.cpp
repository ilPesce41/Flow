#include "feature_set.hpp"

FeatureSet::FeatureSet(int * features_, int k_)
{
    k = k_;
    for(int i=0;i<k;i++)
    {
        features[i] = features_[i];
    }
}

bool FeatureSet::is_equivalent(FeatureSet set)
{
    for(int i=0;i<k;i++)
    {
        bool in_list=false;
        int feat = features[i];
        for(int j=0;j<k;j++)
        {
            if(set.features[j]==feat)
            {
                in_list = true;
                break;
            }
        if(!in_list){return false;}
        }
    }
    return true;
}