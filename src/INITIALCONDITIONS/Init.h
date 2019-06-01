// Init.h is adapted from the IP-Glasma solver.
// Copyright (C) 2018 Bjoern Schenke.
// Reproduced with permission by Mark Mace 2019 for Dilute-Dense gluon solver

#ifndef Init_H
#define Init_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <limits>
#include <ctime>

#include "GIPSLattice.h"
#include "Parameters.h"
#include "Matrix.h"
#include "Random.h"
#include "Group.h"
#include "FFT.h"
#include "Glauber.h"
//#include "GaugeFix.h"
#include "gsl/gsl_linalg.h"

using namespace std;

class Init {
    
    private:
    int const static iymaxNuc = 44; // for the Tp-y table
    
    int const static iTpmax = 100;
    
    double const deltaYNuc = 0.25; // for the new table
    
    FFT *fft;
    //  Matrix** A;
    //  Glauber *glauber;
    double Qs2Nuclear[iTpmax][iymaxNuc];
    double Tlist[iTpmax];
    
    double As[1];
    
    vector<ReturnValue> nucleusA;  // list of x and y coordinates of nucleons in nucleus A
    vector<ReturnValue> nucleusB;  // list of x and y coordinates of nucleons in nucleus B
    
    public:
    
    // Constructor.
    
    Init(const int nn[])
    {
        //  random = new Random;
        fft = new FFT(nn);
        //      glauber = new Glauber;
    };
    
    ~Init()
    {
        delete fft;
        //delete glauber;
        //delete random;
    };
    
    void init(GIPSLattice *lat, Group *group, Parameters *param, Random *random, Glauber* glauber, int READFROMFILE,string OutDirectory,int RNGSEED,int OUTPUT_FLAG);
    void sampleTA(Parameters *param, Random *random, Glauber* glauber);
    void sampleProjectileTA(Parameters *param, Random *random, Glauber* glauber);
    void readNuclearQs(Parameters *param);
    //vector <complex<double> > solveAxb(Parameters *param, complex<double>* A, complex<double>* b);
    double getNuclearQs2(Parameters *param, Random *random, double Qs2atZeroY, double y);
    void setColorChargeDensity(GIPSLattice *lat, Parameters *param, Random *random, Glauber *glauber,string OutDirectory,int RNGSEED,int OUTPUT_FLAG);
    void setColorChargeDensityMVTarget(GIPSLattice *lat, Parameters *param, Random *random, Glauber *glauber,string OutDirectory,int RNGSEED,int OUTPUT_FLAG);

//    void setV(Lattice *lat, Group* group, Parameters *param, Random* random, Glauber *glauber);
//    void readV(Lattice *lat, Group* group, Parameters *param);
    // void eccentricity(Lattice *lat, Group *group, Parameters *param, Random *random, Glauber *glauber);
//    void multiplicity(Lattice *lat, Group *group, Parameters *param, Random *random, Glauber *glauber);
    
};

#endif // Init_H
