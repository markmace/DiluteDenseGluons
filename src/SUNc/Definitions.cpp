// Definitions.cpp is part of the Dilute-Dense Gluon solver //
// Copyright (C) 2019 Mark Mace //

#ifndef __SUNc__DEFINITIONS_CPP__
#define __SUNc__DEFINITIONS_CPP__

//COPYING GAUGE LINKS 
#define COPY_SUNcMatrix(Destination,Origin)    std::memcpy((Destination),(Origin),(sizeof(SU_Nc_FUNDAMENTAL_FORMAT)*SUNcGroup::MatrixSize))
#define COPY_SUNcMatrixAdj(Destination,Origin)    std::memcpy((Destination),(Origin),(sizeof(SU_Nc_ADJOINT_FORMAT)*SUNcGroup::AdjointSize))

//DEFINITIONS OF TYPES FOR SU(3) GAUGE GROUP
#if SU_Nc_FLAG==SU3_FLAG

#define Nc 3

typedef DOUBLE SU_Nc_ALGEBRA_FORMAT;
typedef COMPLEX SU_Nc_FUNDAMENTAL_FORMAT;
typedef DOUBLE SU_Nc_ADJOINT_FORMAT;

#include "SU3/util/Diagonalization.cpp"
#include "SU3/util/ComplexEigensystem.cpp"
#include "SU3/GroupOperations.cpp"
#include "SU3/AlgebraOperations.cpp"

#endif

// NO OTHER GAUGE GROUPS ARE INCLUDED AT THIS TIME //

//INCLUDE COMPOSITE MATRIX PRODUCTS
#include "AdvancedOperations.cpp"


#endif
