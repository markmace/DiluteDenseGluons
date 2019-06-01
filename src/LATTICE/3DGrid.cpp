// 3DGrid.cpp is part of the Dilute-Dense Gluon solver //
// Copyright (C) 2019 Mark Mace //

#ifndef __3DGRID_CPP__
#define __3DGRID_CPP__

namespace Lattice{

    //DIMENSION OF THE TRANSVERSE LATTICE
    static const int Dimension=2;
    
    INT N[2]={512,512};    DOUBLE a[2]={0.0625,0.0625};
    INT NRap=100; DOUBLE aeta=1;
    
    DOUBLE Area=N[0]*N[1];
    
    DOUBLE SizeX=N[0]*a[0];
    DOUBLE SizeY=N[1]*a[1];
    
}

#endif
