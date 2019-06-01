// GIPSLattice.cpp is adapted from the IP-Glasma solver.
// Copyright (C) 2018 Bjoern Schenke.
// Reproduced with permission by Mark Mace 2019 for Dilute-Dense gluon solver

#include "GIPSLattice.h"

//constructor
GIPSLattice::GIPSLattice(Parameters *param, int N, int length){
    
    locNc = N;
    size = length*length;
    double a = param->getL()/static_cast<double>(length);
    
    std::cerr << "# ALLOCATING LATTICE FOR INITIAL FIELDS OF SIZE " << length << "x" << length << " WITH a=" << a << " fm" << std::endl;;
    
    // initialize the array of cells
    //  cells = new Cell*[size];
    for(int i=0; i<size; i++)
    {
        Cell* cell;
        cell = new Cell(locNc);
        cells.push_back(cell);
    }
    
    //std::cerr << "# FINISHED ALLOCATING INITAL LATTICE" << endl;
}

GIPSLattice::~GIPSLattice()
{
    for(int i=0; i<size; i++)
    delete cells[i];
    cells.clear();
    //delete[] cells;
}

