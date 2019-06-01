// Cell.cpp is adapted from the IP-Glasma solver.
// Copyright (C) 2018 Bjoern Schenke.
// Reproduced with permission by Mark Mace 2019 for Dilute-Dense gluon solver

#include "Cell.h"

Cell::Cell(int N)
{
    locNc = N;
    int Nc2m1 = locNc*locNc-1;
    g2mu2A = new double;
    g2mu2B = new double;
    TpA = new double;
    TpB = new double;
    U = new Matrix(locNc,1.);
    U2 = new Matrix(locNc,1.);
    Ux = new Matrix(locNc,1.);
    Uy = new Matrix(locNc,1.);
    Ux1 = new Matrix(locNc,1.);
    Uy1 = new Matrix(locNc,1.);
    Ux2 = new Matrix(locNc,1.);
    Uy2 = new Matrix(locNc,1.);
}

Cell::~Cell()
{
    delete g2mu2A;
    delete g2mu2B;
    delete TpA;
    delete TpB;
    delete U;
    delete U2;
    delete Ux;
    delete Uy;
    delete Ux1;
    delete Ux2;
    delete Uy1;
    delete Uy2;
}
