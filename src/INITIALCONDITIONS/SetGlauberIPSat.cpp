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
                std::cout << "# BEGINNING COMPUTATION FOR SLICE NRap=" << RapSlice << " OF " << Lattice::NRap << std::endl;
            }
            
            // SET TARGET RHO FIELDS /
            for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                for(INT y=0;y<U->N[1];y++){
                    for(INT x=0;x<U->N[0];x++){
                        
                        // OPTION FOR UNIFORM COLOR CHARGE DISTRIBUTION  //
                        //FourierSpace::RhoT->SetXc(x,y,a,RandomNumberGenerator::Gauss(g2muT/(Lattice::a[0]*sqrt(Lattice::NRap))));
                        // END OPTION //
                        
                        // OPTION FOR LOCALIZED COLOR CHARGE DISTRIBUTION  //
                        // [g \mu]=fm^-1 //
                        DOUBLE Rho0=gSQRmu->Get(x,y,0,0)[0]/(Lattice::a[0]*sqrt(Lattice::NRap)); // fm^-2 //
                        FourierSpace::RhoT->SetXc(x,y,a,Rho0*RandomNumberGenerator::Gauss(1.0)); // fm^-2 //
                        // END OPTION //
                        
                    }
                }
                
            }
            
            // FOURIER TRANSFORM RHO TO MOMENTUM SPACE //
            FourierSpace::RhoT->ExecuteXtoP();
            
            // RENORMALIZE AFTER X->P -- [d^2x] = fm^2 //
            for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                for(INT y=0;y<U->N[1];y++){
                    for(INT x=0;x<U->N[0];x++){
                        
                        COMPLEX RhoTTEMP=FourierSpace::RhoT->GetP(x,y,a); // fm^-2 //
                        FourierSpace::RhoT->SetP(x,y,a,RhoTTEMP*COMPLEX(Lattice::a[0]*Lattice::a[1])); // DIMLESS //
                    }
                }
            }
            
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
                                COMPLEX rhoTemp=FourierSpace::RhoT->GetP(pXIndex,pYIndex,a); // DIMLESS //
                                
                                // SAVE RHO/DERIVATIVE^2 //
                                FourierSpace::RhoT->SetP(pXIndex,pYIndex,a,rhoTemp/DerivativeFactorSqr); // fm^2 //
                                
                            }
                            
                        }
                        else{
                            
                            // SET ZERO POINT //
                            for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                                FourierSpace::RhoT->SetP(0,0,a,COMPLEX(0.0));
                            }
                            
                        }
                    }
                    else{ // REGULATOR MASS OPTION //
                        
                        for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                            // SAVE MOMENTUM SPACE RHO //
                            COMPLEX rhoTemp=FourierSpace::RhoT->GetP(pXIndex,pYIndex,a);  // DIMLESS //
                            
                            // SAVE RHO/DERIVATIVE //
                            DOUBLE mLatSqr=RegMass*RegMass*GeVtoinvfm*GeVtoinvfm; // fm^-2 //
                            FourierSpace::RhoT->SetP(pXIndex,pYIndex,a,rhoTemp/(DerivativeFactorSqr+mLatSqr)); // fm^2 //
                            
                        }
                        
                    }
                }
            }
            
            // FOURIER TRANSFORM BACK //
            FourierSpace::RhoT->ExecutePtoX();
            
            // RENORMALIZE AFTER P->X -- [d^2p] = [1/(Na)^2] = fm^-2 //
            for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                for(INT y=0;y<U->N[1];y++){
                    for(INT x=0;x<U->N[0];x++){
                        
                        COMPLEX RhoTTEMP=FourierSpace::RhoT->GetX(x,y,a); // fm^-2
                        FourierSpace::RhoT->SetXc(x,y,a,RhoTTEMP/COMPLEX(Lattice::N[0]*Lattice::N[1]*Lattice::a[0]*Lattice::a[1])); // DIMLESS //
                        
                    }
                }
            }
            
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
                        alphaTemp[a]=real(FourierSpace::RhoT->GetX(x,y,a)); // DIMLESS
                        
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
                    // [g \mu]=fm^-1 //
                    DOUBLE Rho0=gSQRmu->Get(x,y,0,0)[0]/Lattice::a[0];  // fm^-2 //
                    FourierSpace::RhoP->SetX(x,y,a,Rho0*RandomNumberGenerator::Gauss(1.0));  // fm^-2 //
                    
                    // END OPTION //
                    
                    // OPTION FOR UNIFORM COLOR CHARGE DISTRIBUTION  //
                    //FourierSpace::RhoP->SetX(x,y,a,RandomNumberGenerator::Gauss(g2muP/Lattice::a[0]));
                    // END OPTION //
                    
                }
                
            }
        }
        
        //std::cout << "# PERFORMING FOURIER TRANSFORM OF PROJECTILE RHO VARIABLES" << std::endl;
        
        // FOURIER TRANSFORM RHOS TO MOMENTUM SPACE //
        FourierSpace::RhoP->ExecuteXtoP();
        
        // RENORMALIZE AFTER X->P -- [d^2x] = [a^2] = fm^2 //
        for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
            for(INT y=0;y<ProjSolution::A->N[1];y++){
                for(INT x=0;x<ProjSolution::A->N[0];x++){
                    
                    COMPLEX RhoPTEMP=FourierSpace::RhoP->GetP(x,y,a); // fm^-2 //
                    FourierSpace::RhoP->SetP(x,y,a,RhoPTEMP*COMPLEX(Lattice::a[0]*Lattice::a[1])); // DIMLESS //
                }
            }
        }
        
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
                            COMPLEX rhoTemp=FourierSpace::RhoP->GetP(pXIndex,pYIndex,a); // DIMLESS //
                            
                            // SAVE RHO/DERIVATIVE^2 //
                            FourierSpace::RhoP->SetP(pXIndex,pYIndex,a,rhoTemp/DerivativeFactorSqr); // fm^2 //
                            
                        }
                        
                    }
                    else{
                        
                        // SET ZERO POINT //
                        for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                            FourierSpace::RhoP->SetP(0,0,a,COMPLEX(0.0));
                        }
                        
                    }
                }
                else{ // REGULATOR MASS OPTION //
                    
                    for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                        // SAVE MOMENTUM SPACE RHO //
                        COMPLEX rhoTemp=FourierSpace::RhoP->GetP(pXIndex,pYIndex,a);  // DIMLESS //
                        
                        // SAVE RHO/DERIVATIVE //
                        DOUBLE mLatSqr=RegMass*RegMass*GeVtoinvfm*GeVtoinvfm;  // fm^-2 //
                        FourierSpace::RhoP->SetP(pXIndex,pYIndex,a,rhoTemp/(DerivativeFactorSqr+mLatSqr)); // fm^2 //
                        
                    }
                    
                }
            }
        }
        
        // FOURIER TRANSFORM BACK //
        FourierSpace::RhoP->ExecutePtoX();
        
        // RENORMALIZE AFTER P->X -- [d^2p] = [1/(Na)^2] = fm^-2 //
        for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
            for(INT y=0;y<ProjSolution::A->N[1];y++){
                for(INT x=0;x<ProjSolution::A->N[0];x++){
                    
                    COMPLEX RhoPTEMP=FourierSpace::RhoP->GetX(x,y,a); // fm^2
                    FourierSpace::RhoP->SetXc(x,y,a,RhoPTEMP/COMPLEX(Lattice::N[0]*Lattice::N[1]*Lattice::a[0]*Lattice::a[1])); // DIMLESS //
                    
                }
            }
        }
        
        // SAVE TO SUNc VECTOR //
        #pragma omp parallel for collapse(3)
        for(INT y=0;y<ProjSolution::A->N[1];y++){
            for(INT x=0;x<ProjSolution::A->N[0];x++){
                for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                    ProjSolution::A->Get(x,y,0,a)[0]=real(FourierSpace::RhoP->GetX(x,y,a)); // DIMLESS //
                }
            }
        }
        
    }
    
    // GENERAL ROUTINE FOR SETTING GLAUBER + IP-Sat FIELDS //
    void SetGlauberIPSat(long int RNG_SEED){
        
        ///////////////////////////////
        //   SET INITIAL CONDITIONS  //
        ///////////////////////////////
        
        //std::cout << "# SETTING GLAUBER+IP-Sat CONDITIONS " << std::endl;
        
        // SAMPLE NUCLEON POSITIONS AND DETMERINE RHOS //
        //std::cout << "# SAMPLING NUCLEI " << std::endl;
        
        DOUBLE *Projectile_g2mu2=new DOUBLE[Lattice::N[0]*Lattice::N[1]];
        DOUBLE *Target_g2mu2=new DOUBLE[Lattice::N[0]*Lattice::N[1]];
        
        bImpact=Initialize(IO::InputFile,IO::OutputDirectory,RNG_SEED,Projectile_g2mu2,Target_g2mu2,OUTPUT_FLAG);
        
        // COPY TO LOCAL FIELDS FROM GLAUBER IPSAT //
        for(INT yLoc=0;yLoc<Lattice::N[1];yLoc++){
            for(INT xLoc=0;xLoc<Lattice::N[0];xLoc++){
                
                INT Index=xLoc*Lattice::N[0]+yLoc;
                
                // g_s==1 //
                // [Projectile_g2mu2] = DIMLESS //
                // [g2mu]=fm^-1 //
                g2mu::Proj->Get(xLoc,yLoc,0,0)[0]=sqrt(Projectile_g2mu2[Index]/(Lattice::a[0]*Lattice::a[1]));
                g2mu::Targ->Get(xLoc,yLoc,0,0)[0]=sqrt(Target_g2mu2[Index]/(Lattice::a[0]*Lattice::a[1]));
                
            }
        }
        
        // CLEAN-UP //
        delete[] Projectile_g2mu2;
        delete[] Target_g2mu2;
        
        //std::cout << "# FINSIHED SAMPLING NUCLEI " << std::endl;
        
        // OPTION FOR FILE OUTPUT OF RHOS //
        if(OUTPUT_FLAG==2){
            std::cout << "# OUTPUTTING RHOS " << std::endl;
            IO::SaveConfigurationE("ProjRho",g2mu::Proj);
            IO::SaveConfigurationE("TargetRho",g2mu::Targ);
        }
        // END OPTION //
        
        // SET FIELDS FROM g2mu RHO FIELDS //
        //std::cout << "# SETTING TARGET FIELDS " << std::endl;
        
        SetTarget(TargetFields::U,TempTargetFields::U,g2mu::Targ);
        
        //std::cout << "# FINSIHED SETTING TARGET FIELDS " << std::endl;
        
        //std::cout << "# SETTING TARGET FIELDS " << std::endl;
        
        SetProjectile(g2mu::Proj);
        
        //std::cout << "# FINISHED SETTING TARGET FIELDS " << std::endl;
        
        std::cout << "# FINISHED SETTING GLAUBER+IP-Sat INITIAL CONDITIONS" << std::endl;
        
        
    }
    
    
}

#endif
