// SetMV.cpp is part of the Dilute-Dense Gluon solver //
// Copyright (C) 2019 Mark Mace //

#ifndef __SETMV__CPP__
#define __SETMV__CPP__

////////////////////////////////////////////
//   GAUSSIAN MV MODEL INITIAL CONDITIONS //
////////////////////////////////////////////

#include "../LATTICE/FourierSpace.cpp"

// HULTHEN WAVEFUNCTION FOR DEUTERON -- r IS N-P DISTANCE //
DOUBLE Hulthen(DOUBLE r){
    
    DOUBLE a=0.228; // fm^-1 -- arXiv 1408.2549 //
    DOUBLE b=1.177; // fm^-1 -- arXiv 1408.2549 //
    DOUBLE diff=exp(-a*r)-exp(-b*r);
    DOUBLE norm=1.0/sqrt(2.0*PI)*sqrt(a*b*(a+b))/(b-a);
    DOUBLE density=pow(norm*diff/(r+1e-8),2);
    return density;
}

// SAMPLE DEUTERON Rnp DISTANCE -- RETURN POSITION FROM ORIGIN //
void SampleDeuteronNucleon(DOUBLE &x,DOUBLE &y,DOUBLE &z){
    // REJECTION SAMPLING //
    bool cont=true;
    do{
        DOUBLE trial=RandomNumberGenerator::rng();
        x=RandomNumberGenerator::rng()*20.0;
        y=RandomNumberGenerator::rng()*20.0;
        z=RandomNumberGenerator::rng()*20.0;
        DOUBLE r=sqrt(x*x+y*y+z*z);
        DOUBLE density=Hulthen(r);
        // SUCCESS //
        if(trial<=density){
            cont=false;
        }
    }while(cont);
}

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
                        FourierSpace::RhoT->SetXc(x,y,a,RandomNumberGenerator::Gauss(g2muT*GeVtoinvfm/(Lattice::a[0]*sqrt(Lattice::NRap)))); // fm^-2 //
                        // END OPTION //

                    }
                }
                
            }
            
            //std::cerr << "# PERFORMING FOURIER TRANSFORM OF TARGET RHO VARIABLES" << std::endl;
            
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
            
            // SET RHO OVER DERIVATIVE //
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
                                
                                // SAVE RHO/DERIVATIVE //
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
                    else{
                        
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
                        alphaTemp[a]=real(FourierSpace::RhoT->GetX(x,y,a));
                        
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
            
        }// END RAPIDITY LOOP

    }
    
    
    /////////////////////////////
    // SET PROJECTILE FIELDS A //
    /////////////////////////////
    
    void SetProjectile(){

        // OPTION FOR 'PROTON' //
        if(PROJ_FLAG==0){
            
            std::cerr << "# SETTING MV PROTON PROJECTILE " << std::endl;
            
            // SET PROJECTILE PARAMETERS -- SINGLE BOUNDED NUCLEUS //
            // ECCENTRICITY //
            DOUBLE aP=1.0; // DIMLESS //
            DOUBLE bP=1.0; // DIMLESS //
            // 'PROTON' SIZE PARAMETER -- R_p~\Lambda~1/(g2muP)//
            //DOUBLE Rp=DOUBLE(1.0/g2muP)*invGeVtofm; // fm //
            DOUBLE Rp=sqrt(4.0/PI)*invGeVtofm; // ASSUMING Bp=4.0GeV^-2 = PI Rp^2 -- fm //

            // SET RHO FIELDS /
            for(INT y=0;y<ProjSolution::A->N[1];y++){
                for(INT x=0;x<ProjSolution::A->N[0];x++){
                    
                    // SET MIDPOINT FOR LOCALIZED COLOR CHARGE DISTRIBUTION//
                    DOUBLE xMid=Lattice::a[0]*Lattice::N[0]/2; // fm //
                    DOUBLE yMid=Lattice::a[1]*Lattice::N[1]/2; // fm //
                    DOUBLE xPhys=x*Lattice::a[0]; // fm //
                    DOUBLE yPhys=y*Lattice::a[1]; // fm //

                    for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                        
                        DOUBLE arg=std::pow((xPhys-xMid)/(aP*Rp),2.0)+std::pow((yPhys-yMid)/(bP*Rp),2.0); // DIMLESS //
                        FourierSpace::RhoP->SetX(x,y,a,RandomNumberGenerator::Gauss(g2muP*GeVtoinvfm/Lattice::a[0])*exp(-arg/2.0)); // fm^-2 //


                    }
                    
                }
             }
        }
        // END 'PROTON' OPTION //
        
        // OPTION FOR 'DEUTERON ' //
        else if(PROJ_FLAG==1){
            
            std::cerr << "# SETTING MV DEUTERON PROJECTILE " << std::endl;

            // POSITION BUFFERS //
            DOUBLE Nucleon1x,Nucleon1y,Nucleon1z,Nucleon2x,Nucleon2y,Nucleon2z;
            // SAMPLE POSITIONS OF NUCLEONS FROM HULTHEN //
            Nucleon1x=0.0; // SET FIRST NUCLEON AT ORIGIN //
            Nucleon1y=0.0; // SET FIRST NUCLEON AT ORIGIN //
            Nucleon1z=0.0; // SET FIRST NUCLEON AT ORIGIN //

            SampleDeuteronNucleon(Nucleon2x,Nucleon2y,Nucleon2z); // SAMPLE SECOND NUCLEON //
            
            // NUCLEON TRANSVERSE SEPARATION //
            DOUBLE RnpTrans=sqrt(pow(Nucleon1x-Nucleon2x,2)+pow(Nucleon1y-Nucleon2y,2));
            DOUBLE Rnp3Vol=sqrt(pow(Nucleon1x-Nucleon2x,2)+pow(Nucleon1y-Nucleon2y,2)+pow(Nucleon1z-Nucleon2z,2));

            // CREATE OUTPUT STREAM FOR DEUTERON NUCLEONS //
            if(OUTPUT_FLAG==1 || OUTPUT_FLAG==2){
                
                std::ofstream DOutStream; DOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"DeuteronPositionID",MY_MPI_RNG_SEED,".txt").c_str());

                DOutStream << "# Nucleon1x Nucleon1y Nucleon1z Nucleon2x Nucleon2y Nucleon2z RnpTrans Rnp3Vol" << std::endl;

                DOutStream << Nucleon1x << " " << Nucleon1y  << " " << Nucleon1z  << " " << Nucleon2x  << " " << Nucleon2y << " " << Nucleon2z << " " << RnpTrans << " " << Rnp3Vol << std::endl;
                
                DOutStream.close();
                
            }
            
            // NUCLEON ECCENTRICITY //
            DOUBLE aP=1.0; // DIMLESS //
            DOUBLE bP=1.0; // DIMLESS //
            
            // 'NUCLEON' SIZE PARAMETER -- R_p~\Lambda~1/(g2muP)//
            //DOUBLE Rp=DOUBLE(1.0/g2muP)*invGeVtofm; // fm //
            DOUBLE Rp=sqrt(4.0/PI)*invGeVtofm; // ASSUMING Bp=4.0GeV^-1 = PI Rp^2 -- fm //

            // SET RHO FIELDS /
            for(INT y=0;y<ProjSolution::A->N[1];y++){
                for(INT x=0;x<ProjSolution::A->N[0];x++){
                    
                    // SET MIDPOINT FOR LOCALIZED COLOR CHARGE DISTRIBUTION//
                    DOUBLE xMid=Lattice::a[0]*Lattice::N[0]/2; // fm //
                    DOUBLE yMid=Lattice::a[1]*Lattice::N[1]/2; // fm //
                    DOUBLE xPhys=x*Lattice::a[0]; // fm //
                    DOUBLE yPhys=y*Lattice::a[1]; // fm //
                    
                    for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                        
                        DOUBLE argNuc1=std::pow(((xPhys-xMid)-Nucleon1x)/(aP*Rp),2.0)+std::pow(((yPhys-yMid)-Nucleon1y)/(bP*Rp),2.0); // DIMLESS //
                        DOUBLE argNuc2=std::pow(((xPhys-xMid)-Nucleon2x)/(aP*Rp),2.0)+std::pow(((yPhys-yMid)-Nucleon2y)/(bP*Rp),2.0); // DIMLESS //
                        // NEED TO CHECK //
                        DOUBLE deuteronProfile=sqrt(exp(-argNuc1)+exp(-argNuc2));

                        FourierSpace::RhoP->SetX(x,y,a,RandomNumberGenerator::Gauss(g2muP*GeVtoinvfm/Lattice::a[0])*deuteronProfile); // fm^-2 //
          
                        
                    }
                    
                }
            }
        }
        // END 'DEUTERON' OPTION //
        
        // OPTION FOR INFINITE UNIFORM COLOR CHARGE DISTRIBUTION  //
        else if(PROJ_FLAG==2){
            
            std::cerr << "# SETTING MV PROTON PROJECTILE " << std::endl;

            for(INT y=0;y<ProjSolution::A->N[1];y++){
                for(INT x=0;x<ProjSolution::A->N[0];x++){
                    
                    for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
     
                        FourierSpace::RhoP->SetX(x,y,a,RandomNumberGenerator::Gauss(g2muP*GeVtoinvfm/Lattice::a[0])); // fm^-2 //
             
                    }
                    
                }
             }
        }
        // END INFINITE NUCLEUS OPTION //
        else{
            
            std::cerr << "# NO MV PROJECTILE SET -- FATAL ERROR " << std::endl;
            exit(0);
        }
        //std::cerr << "#PERFORMING FOURIER TRANSFORM OF PROJECTILE RHO VARIABLES" << std::endl;
        
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
                    
                    COMPLEX RhoPTEMP=FourierSpace::RhoP->GetX(x,y,a); // fm^-2
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
    
    void SetMVTarget(){
        
        ///////////////////////////////
        //   SET INITIAL CONDITIONS  //
        ///////////////////////////////
        
        std::cerr << "# SETTING MV TARGET" << std::endl;
        
        SetTarget(TargetFields::U,TempTargetFields::U);
        
        std::cerr << "# FINISHED SETTING MV TARGET" << std::endl;
        
    }
    
    void SetMVProjectile(){
        
        ///////////////////////////////
        //   SET INITIAL CONDITIONS  //
        ///////////////////////////////
        
        //std::cerr << "# SETTING MV PROJECTILE" << std::endl;
        
        SetProjectile();
        
        std::cerr << "# FINISHED SETTING MV PROJECTILE" << std::endl;
        
    }
    
}


#endif
