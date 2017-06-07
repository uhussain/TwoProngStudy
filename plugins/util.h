// -*- C++ -*-                
//                            
// Package:     -
// File:        util.h       
//                            
//                            

#include<cmath>
//template<typename Comp1, typename Comp2, class T1, class T2>
//bool pairCompare(const std::pair<Comp1, T1>&p1, const std::pair<Comp2, T2>& p2) 
//{return p1.first < p2.first;}

bool pairCompare(const std::pair<double, const pat::PackedCandidate*>& firstElem, const std::pair<double, const pat::PackedCandidate*>& secondElem) {
    return firstElem.first < secondElem.first;

}
