// SetGlauberIPSat.cpp is part of the Dilute-Dense Gluon solver //
// Copyright (C) 2019 Mark Mace //

#ifndef __SETGLAUBERIPSAT__CPP__
#define __SETGLAUBERIPSAT__CPP__

////////////////////////////////////////////
//    GLAUBER + IP-SAT INITIAL CONDITIONS //
////////////////////////////////////////////

#include "../LATTICE/FourierSpace.cpp"

double Initialize(std::string InputFile,std::string OutDirectory,long int RNG_SEED, double *PROJ_g2mu2,double *TARG_g2mu2,int OUTPUT_FLAG);


namespace InitialConditions{

    ///////////////////////////////
    // SET TARGET WILSON LINES U //
    ///////////////////////////////

    void SetTarget(WilsonLines *U,WilsonLines *Utemp,VectorFields *gSQRmu){
        
        // RESET TARGET LINKS //
        U->SetIdentity();
        
        // LOOP OVER SLICES IN RAPIDITY //
        for(INT RapSlice=0;RapSlice<Lattice::NRap;RapSlice++){
            if(RapSlice%20==0){
                std::cerr << "# BEGINNING COMPUTATION FOR SLICE NRap=" << RapSlice << " OF " << Lattice::NRap << std::endl;
            }
            
            // SET TARGET RHO FIELDS /
            for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                for(INT y=0;y<U->N[1];y++){
                    for(INT x=0;x<U->N[0];x++){
                        
                        // OPTION FOR UNIFORM COLOR CHARGE DISTRIBUTION  //
                        //FourierSpace::RhoT->SetXc(x,y,a,RandomNumberGenerator::Gauss(g2muT/(Lattice::a[0]*sqrt(Lattice::NRap))));
                        // END OPTION //
                        
                        // OPTION FOR LOCALIZED COLOR CHARGE DISTRIBUTION  //
                        // [g \mu]=fm^-1 -- NATIVELY [g \mu]= DIMLESS //
                        DOUBLE Rho0=gSQRmu->Get(x,y,0,0)[0]/(Lattice::a[0]*sqrt(Lattice::NRap)); // fm^-2 //
                        FourierSpace::RhoT->SetXc(x,y,a,Rho0*RandomNumberGenerator::Gauss(1.0)); // fm^-2 //
                        // END OPTION //

                    }
                }
                
            }
            
            // FOURIER TRANSFORM RHO TO MOMENTUM SPACE //
            FourierSpace::RhoT->ExecuteXtoP();
            
            // SET TARGET RHO OVER DERIVATIVE //
            #pragma omp parallel for collapse(2)
            for(INT pXIndex=0;pXIndex<Lattice::N[0];pXIndex++){
                for(INT pYIndex=0;pYIndex<Lattice::N[1];pYIndex++){
                    
                    // DEFINE MOMENTUM FACTOR // PAGE 1054 NR //
                    DOUBLE DerivativeFactorSqr=(cos(2.0*PI*pXIndex/Lattice::N[0])+cos(2.0*PI*pYIndex/Lattice::N[1])-2.0)/(-0.5*Lattice::a[0]*Lattice::a[1]);  // fm^-2 //
                    
                    if(RegMass==0){
                        if(DerivativeFactorSqr!=0.0){
                            
                            for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                                // SAVE MOMENTUM SPACE RHO //
                                COMPLEX rhoTemp=FourierSpace::RhoT->GetP(pXIndex,pYIndex,a); // fm^-2 //
                                
                                // SAVE RHO/DERIVATIVE^2 //
                                FourierSpace::RhoT->SetP(pXIndex,pYIndex,a,rhoTemp/DerivativeFactorSqr); // DIMLESS //
                                
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
                            COMPLEX rhoTemp=FourierSpace::RhoT->GetP(pXIndex,pYIndex,a);  // fm^-2 //
                            
                            // SAVE RHO/DERIVATIVE //
                            DOUBLE mLatSqr=RegMass*RegMass*GeVtoinvfm*GeVtoinvfm; // fm^-2 //
                            FourierSpace::RhoT->SetP(pXIndex,pYIndex,a,rhoTemp/(DerivativeFactorSqr+mLatSqr)); // DIMLESS //
                            
                        }
                        
                    }
                }
            }
            
            // FOURIER TRANSFORM BACK //
            FourierSpace::RhoT->ExecutePtoX();
            
            // RESET TEMP TARGET FIELDS //
            Utemp->SetIdentity();
            
            // EXPONENTIATE SOLUTION AND SAVE TO TEMP TARGET FIELD //
            #pragma omp parallel for collapse(2)
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
            //GetQsT(TargetFields::U);
            // END OPTION //
            
        }
        
    }

    /////////////////////////////
    // SET PROJECTILE FIELDS A //
    /////////////////////////////

    void SetProjectile(VectorFields *gSQRmu){
        
        // SET PROJECTILE RHO FIELDS /
        for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
            for(INT y=0;y<ProjSolution::A->N[1];y++){
                for(INT x=0;x<ProjSolution::A->N[0];x++){
                
                
                    // OPTION FOR LOCALIZED COLOR CHARGE DISTRIBUTION  //
                    // [g \mu]=fm^-1 -- NATIVELY [g \mu]= DIMLESS //
                    DOUBLE Rho0=gSQRmu->Get(x,y,0,0)[0]/Lattice::a[0];  // fm^-2 //
                    FourierSpace::RhoP->SetX(x,y,a,Rho0*RandomNumberGenerator::Gauss(1.0));  // fm^-2 //

                    // END OPTION //
                    
                    // OPTION FOR UNIFORM COLOR CHARGE DISTRIBUTION  //
                    //FourierSpace::RhoP->SetX(x,y,a,RandomNumberGenerator::Gauss(g2muP/Lattice::a[0]));
                    // END OPTION //

                }
                
            }
        }
        
        //std::cerr << "# PERFORMING FOURIER TRANSFORM OF PROJECTILE RHO VARIABLES" << std::endl;
        
        // FOURIER TRANSFORM RHOS TO MOMENTUM SPACE //
        FourierSpace::RhoP->ExecuteXtoP();
        
        // SET PROJECTILE RHO OVER DERIVATIVE //
        #pragma omp parallel for collapse(2)
        for(INT pXIndex=0;pXIndex<Lattice::N[0];pXIndex++){
            for(INT pYIndex=0;pYIndex<Lattice::N[1];pYIndex++){
                
                // DEFINE MOMENTUM FACTOR // PAGE 1054 NR //
                DOUBLE DerivativeFactorSqr=(cos(2.0*PI*pXIndex/Lattice::N[0])+cos(2.0*PI*pYIndex/Lattice::N[1])-2.0)/(-0.5*Lattice::a[0]*Lattice::a[1]);  // fm^-2 //
                
                if(RegMass==0){
                    if(DerivativeFactorSqr!=0.0){
                        
                        for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                            // SAVE MOMENTUM SPACE RHO //
                            COMPLEX rhoTemp=FourierSpace::RhoP->GetP(pXIndex,pYIndex,a); // fm^-2 //
                            
                            // SAVE RHO/DERIVATIVE^2 //
                            FourierSpace::RhoP->SetP(pXIndex,pYIndex,a,rhoTemp/DerivativeFactorSqr); // DIMLESS //
                            
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
                        COMPLEX rhoTemp=FourierSpace::RhoP->GetP(pXIndex,pYIndex,a);  // fm^-2 //
                        
                        // SAVE RHO/DERIVATIVE //
//                        FourierSpace::RhoP->SetP(pXIndex,pYIndex,a,rhoTemp*(-0.5*Lattice::a[0]*Lattice::a[1])/(DerivativeFactor+RegMass*RegMass*Lattice::a[0]*Lattice::a[1]*GeVtoinvfm*GeVtoinvfm));
                        DOUBLE mLatSqr=RegMass*RegMass*GeVtoinvfm*GeVtoinvfm;  // fm^-2 //
                        FourierSpace::RhoP->SetP(pXIndex,pYIndex,a,rhoTemp/(DerivativeFactorSqr+mLatSqr)); // DIMLESS //

                        
                    }
                    
                }
            }
        }
        
        // FOURIER TRANSFORM BACK //
        FourierSpace::RhoP->ExecutePtoX();
        
        // SAVE TO SUNc VECTOR //
        #pragma omp parallel for collapse(3)
        for(INT y=0;y<ProjSolution::A->N[1];y++){
            for(INT x=0;x<ProjSolution::A->N[0];x++){
                for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                    ProjSolution::A->Get(x,y,0,a)[0]=real(FourierSpace::RhoP->GetX(x,y,a))/(ProjSolution::A->N[0]*ProjSolution::A->N[1]);
                }
            }
        }
        
    }

    // GENERAL ROUTINE FOR SETTING GLAUBER + IP-Sat FIELDS //
    void SetbGlauberIPSat(long int RNG_SEED){
        
        ///////////////////////////////
        //   SET INITIAL CONDITIONS  //
        ///////////////////////////////
        
        //std::cerr << "## SETTING GLAUBER+IP-Sat CONDITIONS " << std::endl;
                
        // SAMPLE NUCLEON POSITIONS AND DETMERINE RHOS //
        //std::cerr << "# SAMPLING NUCLEI " << std::endl;
        
        DOUBLE *Projectile_g2mu2=new DOUBLE[Lattice::N[0]*Lattice::N[1]];
        DOUBLE *Target_g2mu2=new DOUBLE[Lattice::N[0]*Lattice::N[1]];
        
        bImpact=Initialize(IO::InputFile,IO::OutputDirectory,RNG_SEED,Projectile_g2mu2,Target_g2mu2,OUTPUT_FLAG);
                
        // COPY TO LOCAL FIELDS FROM GLAUBER IPSAT //
        for(INT yLoc=0;yLoc<Lattice::N[1];yLoc++){
            for(INT xLoc=0;xLoc<Lattice::N[0];xLoc++){
                
                INT Index=xLoc*Lattice::N[0]+yLoc;
                
                // g_s==1 //
                // [Projectile_g2mu2] = DIMLESS //
                // [g2mu]=fm^-2 //
                g2mu::Proj->Get(xLoc,yLoc,0,0)[0]=sqrt(Projectile_g2mu2[Index]/(Lattice::a[0]*Lattice::a[1]));
                g2mu::Targ->Get(xLoc,yLoc,0,0)[0]=sqrt(Target_g2mu2[Index]/(Lattice::a[0]*Lattice::a[1]));
                
            }
        }
        
        // CLEAN-UP //
        delete[] Projectile_g2mu2;
        delete[] Target_g2mu2;
        
        //std::cerr << "# FINSIHED SAMPLING NUCLEI " << std::endl;
        
        // OPTION FOR FILE OUTPUT OF RHOS //
        if(OUTPUT_FLAG==2){
            std::cerr << "# OUTPUTTING RHOS " << std::endl;
            IO::SaveConfigurationE("ProjRho",g2mu::Proj);
            IO::SaveConfigurationE("TargetRho",g2mu::Targ);
        }
        // END OPTION //
        
        // SET FIELDS FROM g2mu RHO FIELDS //
        //std::cerr << "# SETTING TARGET FIELDS " << std::endl;
        
        SetTarget(TargetFields::U,TempTargetFields::U,g2mu::Targ);
        
        //std::cerr << "# FINSIHED SETTING TARGET FIELDS " << std::endl;
        
        //std::cerr << "# SETTING TARGET FIELDS " << std::endl;
        
        SetProjectile(g2mu::Proj);
        
        //std::cerr << "# FINISHED SETTING TARGET FIELDS " << std::endl;
        
        std::cerr << "## FINISHED SETTING GLAUBER+IP-Sat INITIAL CONDITIONS" << std::endl;
        
        
    }


}

#endif