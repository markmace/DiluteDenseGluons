// SetMV.cpp is part of the Dilute-Dense Gluon solver //
// Copyright (C) 2019 Mark Mace //

#ifndef __SETMV__CPP__
#define __SETMV__CPP__

// NEEDS TO BE CHECKED !!! //

////////////////////////////////////////////
//   GAUSSIAN MV MODEL INITIAL CONDITIONS //
////////////////////////////////////////////

#include "../LATTICE/FourierSpace.cpp"

namespace InitialConditions{

    ///////////////////////////////
    // SET TARGET WILSON LINES U //
    ///////////////////////////////
    
    void SetTarget(WilsonLines *U,WilsonLines *Utemp){
        
        // RESET TARGET LINKS //
        U->SetIdentity();
        
        // LOOP OVER SLICES IN RAPIDITY //
        for(INT RapSlice=0;RapSlice<Lattice::NRap;RapSlice++){
            if(RapSlice%10==0){
                std::cerr << "# BEGINNING COMPUTATION FOR SLICE NRap=" << RapSlice << " OF " << Lattice::NRap << std::endl;
            }
            
            // SET RHO FIELDS /
            for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                
                for(INT y=0;y<U->N[1];y++){
                    for(INT x=0;x<U->N[0];x++){
                        // OPTION FOR INFINITE MV //
                        FourierSpace::RhoT->SetXc(x,y,a,RandomNumberGenerator::Gauss(g2muT/(Lattice::a[0]*sqrt(Lattice::NRap))));
                        // END OPTION //
                        
                        // OPTION FOR INFINITE MV //
                        //FourierSpace::RhoT->SetXc(x,y,a,3.24*x*x+y*y*10.2+48.*x*y+12.123*a*a+420.*a*x);
                        // END OPTION //

                    }
                }
                
            }
            
            //std::cerr << "# PERFORMING FOURIER TRANSFORM OF TARGET RHO VARIABLES" << std::endl;
            
            // FOURIER TRANSFORM RHO TO MOMENTUM SPACE //
            FourierSpace::RhoT->ExecuteXtoP();
            
            // SET RHO OVER DERIVATIVE //
            for(INT pXIndex=0;pXIndex<Lattice::N[0];pXIndex++){
                for(INT pYIndex=0;pYIndex<Lattice::N[1];pYIndex++){
                    
                    // DEFINE MOMENTUM FACTOR // PAGE 1054 NR //
                    DOUBLE DerivativeFactor=(cos(2.0*PI*pXIndex/Lattice::N[0])+cos(2.0*PI*pYIndex/Lattice::N[1])-2.0);
                    
                    if(RegMass==0){
                        if(DerivativeFactor!=0.0){
                            
                            for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                                // SAVE MOMENTUM SPACE RHO //
                                COMPLEX rhoTemp=FourierSpace::RhoT->GetP(pXIndex,pYIndex,a);
                                
                                // SAVE RHO/DERIVATIVE //
                                FourierSpace::RhoT->SetP(pXIndex,pYIndex,a,rhoTemp*(-0.5*Lattice::a[0]*Lattice::a[1])/DerivativeFactor);
                                
                            }
                        
                        }
                        else{
                            
                            // SET ZERO POINT //
                            for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                                FourierSpace::RhoT->SetP(0,0,a,COMPLEX(0.0));
                            }
                        
                        }
                    }
                    else{
                        
                        for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                            // SAVE MOMENTUM SPACE RHO //
                            COMPLEX rhoTemp=FourierSpace::RhoT->GetP(pXIndex,pYIndex,a);
                            
                            // SAVE RHO/DERIVATIVE //
                            FourierSpace::RhoT->SetP(pXIndex,pYIndex,a,rhoTemp*(-0.5*Lattice::a[0]*Lattice::a[1])/(DerivativeFactor+RegMass*RegMass*Lattice::a[0]*Lattice::a[1]));
                            
                        }
                    
                    }
                }
            }

            // FOURIER TRANSFORM BACK //
            FourierSpace::RhoT->ExecutePtoX();
            
            // RESET TEMP TARGET FIELDS //
            Utemp->SetIdentity();
            
            // EXPONENTIATE SOLUTION AND SAVE TO TEMP TARGET FIELD //
            for(INT y=0;y<U->N[1];y++){
                for(INT x=0;x<U->N[0];x++){
                    
                    // CREATE ALGEBRA alpha ARRAY FOR EXPONENTIATION//
                    SU_Nc_ALGEBRA_FORMAT alphaTemp[SUNcAlgebra::VectorSize];
                    for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                        
                        // SET TARGET FIELDS //
                        alphaTemp[a]=real(FourierSpace::RhoT->GetX(x,y,a))/(U->N[0]*U->N[1]);
                        
                    }
                    // EXPONENTIATE TO GAUGE LINK //
                    SUNcAlgebra::Operations::MatrixIExp(1.0,alphaTemp,Utemp->Get(x,y,0));
                    
                    //MATRIX BUFFERS
                    SU_Nc_FUNDAMENTAL_FORMAT Buffer[SUNcGroup::MatrixSize];
                    
                    //COMPUTE UPDATED LINK AND SAVE TO PATH ORDER //
                    SUNcGroup::Operations::UU(U->Get(x,y,0),Utemp->Get(x,y,0),Buffer);
                    
                    // COPY BUFFER TO UPDATE PATH ORDERED LINK V //
                    COPY_SUNcMatrix(U->Get(x,y,0),Buffer);
                }
            }
            
            // OPTION TO MEASURE Qs AT EACH ADDITIONAL SLICE IN Qs //
            //GetQsT(TargetFields::U); // OBSOLETE //
            // END OPTION //
            
        }// END RAPIDITY LOOP

    }
    
    
    /////////////////////////////
    // SET PROJECTILE FIELDS A //
    /////////////////////////////
    
    void SetProjectile(){
        
        // SET PROJECTILE PARAMETERS //
        // ECCENTRICITY //
        DOUBLE aP=1.0;
        DOUBLE bP=1.0;
        // PROTON SIZE PARAMETER -- R_p~\Lambda~1/(g2muP)//
        DOUBLE Rp=DOUBLE(1.0/g2muP);
        
        // SET RHO FIELDS /
        for(INT y=0;y<ProjSolution::A->N[1];y++){
            for(INT x=0;x<ProjSolution::A->N[0];x++){
                
                // SET MIDPOINT FOR LOCALIZED COLOR CHARGE DISTRIBUTION//
                DOUBLE xMid=Lattice::a[0]*Lattice::N[0]/2;
                DOUBLE yMid=Lattice::a[1]*Lattice::N[1]/2;
                DOUBLE xPhys=x*Lattice::a[0];
                DOUBLE yPhys=y*Lattice::a[1];

                for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                    
                    // OPTION FOR LOCALIZED COLOR CHARGE DISTRIBUTION  //
                    DOUBLE arg=std::pow((xPhys-xMid)/(aP*Rp),2.0)+std::pow((yPhys-yMid)/(bP*Rp),2.0);
                    FourierSpace::RhoP->SetX(x,y,a,RandomNumberGenerator::Gauss(g2muP/Lattice::a[0])*exp(-arg/2.0));
                    // END OPTION //

                    // OPTION FOR UNIFORM COLOR CHARGE DISTRIBUTION  //
                    //FourierSpace::RhoP->SetX(x,y,a,RandomNumberGenerator::Gauss(g2muP/Lattice::a[0]));
                    // END OPTION //
                    
                    // OPTION FOR DEBUGGING //
                    //FourierSpace::RhoP->SetX(x,y,a,5.3*x*x+y*y*9.4+3.431*a*a+69.*y*a);
                    // END OPTION //
                }
                
            }
        }
        
        //std::cerr << "#PERFORMING FOURIER TRANSFORM OF PROJECTILE RHO VARIABLES" << std::endl;
        
        // FOURIER TRANSFORM RHOS TO MOMENTUM SPACE //
        FourierSpace::RhoP->ExecuteXtoP();
        
        // SET RHO OVER DERIVATIVE //
        for(INT pXIndex=0;pXIndex<Lattice::N[0];pXIndex++){
            for(INT pYIndex=0;pYIndex<Lattice::N[1];pYIndex++){
                
                // DEFINE MOMENTUM FACTOR // PAGE 1054 NR //
                DOUBLE DerivativeFactor=(cos(2.0*PI*pXIndex/Lattice::N[0])+cos(2.0*PI*pYIndex/Lattice::N[1])-2.0);
                
                if(RegMass==0){
                    if(DerivativeFactor!=0.0){
                        
                        for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                            
                            // SAVE MOMENTUM SPACE RHO //
                            COMPLEX rhoTemp=FourierSpace::RhoP->GetP(pXIndex,pYIndex,a);
                            
                            // SAVE RHO/DERIVATIVE //
                            FourierSpace::RhoP->SetP(pXIndex,pYIndex,a,rhoTemp*(-0.5*Lattice::a[0]*Lattice::a[1])/DerivativeFactor);
                            
                        }
                        
                    }
                    else{
                        
                        // SET ZERO POINT //
                        for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                            FourierSpace::RhoP->SetP(0,0,a,COMPLEX(0.0));
                        }
                        
                    }
                }
                else{
                    
                    for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                        // SAVE MOMENTUM SPACE RHO //
                        COMPLEX rhoTemp=FourierSpace::RhoP->GetP(pXIndex,pYIndex,a);
                        
                        // SAVE RHO/DERIVATIVE //
                        FourierSpace::RhoP->SetP(pXIndex,pYIndex,a,rhoTemp*(-0.5*Lattice::a[0]*Lattice::a[1])/(DerivativeFactor+RegMass*RegMass*Lattice::a[0]*Lattice::a[1]));
                        
                    }
                    
                }
            }
        }
        
        // FOURIER TRANSFORM BACK //
        FourierSpace::RhoP->ExecutePtoX();
        
        // SAVE TO SUNc VECTOR //
        for(INT y=0;y<ProjSolution::A->N[1];y++){
            for(INT x=0;x<ProjSolution::A->N[0];x++){
                for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                    ProjSolution::A->Get(x,y,0,a)[0]=real(FourierSpace::RhoP->GetX(x,y,a))/(ProjSolution::A->N[0]*ProjSolution::A->N[1]);
                }
            }
        }
  
    }
    
    void SetMVTarget(){
        
        ///////////////////////////////
        //   SET INITIAL CONDITIONS  //
        ///////////////////////////////
        
        std::cerr << "## SETTING MV TARGET" << std::endl;
        
        SetTarget(TargetFields::U,TempTargetFields::U);
        
        std::cerr << "## FINISHED SETTING MV TARGET" << std::endl;
        
        //GetQsT(TargetFields::U);
        
    }
    
    void SetMVProjectile(){
        
        ///////////////////////////////
        //   SET INITIAL CONDITIONS  //
        ///////////////////////////////
        
        std::cerr << "## SETTING MV PROJECTILE" << std::endl;
        
        SetProjectile();
        
        std::cerr << "## FINISHED SETTING MV PROJECTILE" << std::endl;
        
    }
    
}


#endif
