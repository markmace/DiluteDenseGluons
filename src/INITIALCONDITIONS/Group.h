// Group.h is adapted from the IP-Glasma solver.
// Copyright (C) 2018 Bjoern Schenke.
// Reproduced with permission by Mark Mace 2019 for Dilute-Dense gluon solver

#ifndef Group_h
#define Group_h

#include <complex>
#include <iostream>
#include <cstdlib>
#include "Matrix.h"

using namespace std;

class Group
{
    private:
    Matrix** t;  // generators of the group
    Matrix** tA; // adjoint representation of generators of the group
    int locNc; //number of colors
    
    public:
    
    //constructor(s)
    Group(int N);
    ~Group();
    
    Matrix& getT(int i) const { return *t[i]; };
    Matrix& getTA(int i) const { return *tA[i]; };
    
};
#endif

