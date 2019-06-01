// Setup.h is adapted from of the IP-Glasma solver.
// Copyright (C) 2018 Bjoern Schenke.
// Reproduced with permission by Mark Mace 2019 for Dilute-Dense gluon solver

#ifndef Setup_H
#define Setup_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

class Setup {
    
    public:
    
    // Constructor.
    Setup() {}
    
    string StringFind(string file_name, string st);
    int IFind(string file_name, string st);
    unsigned long long int ULLIFind(string file_name, string st);
    double DFind(string file_name, string st);
    int IsFile(string file_name);
    
};

#endif // Setup_H
