// IntegratedDistribution.cpp is part of the Dilute-Dense Gluon solver //
// Copyright (C) 2019 Mark Mace //

#ifndef __INTEGRATED_DISTRIBUTION__CPP__
#define __INTEGRATED_DISTRIBUTION__CPP__

namespace Observables{
    
    // SINGLE INTEGRATED SPECTRA //
    void DetermineIntegratedDistribution(CVectorFields *OmS,CVectorFields *OmA){
        
        // TO DO -- MAKE COMMANDLINE //
        const INT NumOfkConfigs=1;
        DOUBLE kIR[NumOfkConfigs]={0.3}; // [GeV] -- NEEDS TO BE SET AT COMPILE TIME //
        DOUBLE kUV[NumOfkConfigs]={3.0}; // [GeV]
        
        // CREATE OUTPUT STREAM FOR INTEGRATED DISTRIBUTION //
        std::ofstream SIOutStream; SIOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"IntegratedSpectrumID",MY_MPI_RNG_SEED,".txt").c_str());
        
        // LOOP OVER LIMITS //
        for(INT kLimitIndex=0;kLimitIndex<NumOfkConfigs;kLimitIndex++){
            
            // MOMENTUM LIMITS //
            DOUBLE kIRloc=kIR[kLimitIndex];
            DOUBLE kUVloc=kUV[kLimitIndex];
            
            std::cerr << "# kIR=" << kIRloc << " kUV=" << kUVloc << std::endl;
            
            // BUFFER FOR ANGULAR INTEGRATION VARIABLES //
            DOUBLE SI=0.0;
            DOUBLE ASI=0.0;
            
            // COMPUTE Psi_n=\int d^2k F_n(k)=\int d^2k SI(nkx,nky)*exp(i n \phi_k) \theta(k-kIR)\theta(kUV-k) //
            FFT2D *F=new FFT2D(Lattice::N[0],Lattice::N[1],5);

            // BIN COUNTER //
            INT Ncount=0;
            
            // COMPUTE F(nkx,nky)=SI(nkx,nky)*exp(i n \phi_k) \theta(k-kIR)\theta(kUV-k)
            for(INT kXIndex=0;kXIndex<Lattice::N[0];kXIndex++){
                for(INT kYIndex=0;kYIndex<Lattice::N[1];kYIndex++){
                    
                    // PHYSICAL MOMENTUM VALUES //
                    DOUBLE kXValue=IndexTok(kXIndex);
                    DOUBLE kYValue=IndexTok(kYIndex);
                    
                    // DEFINE k^2 AND |k| //
                    DOUBLE kSqr=IndexTokSqr(kXIndex,kYIndex);
                    DOUBLE kAbs=sqrt(kSqr);
                    
                    // DEFINE AZIMUTHAL ANGLE //
                    DOUBLE Phik=atan2(kYValue,kXValue);
                    
                    // MOMENTUM STEP FUNCTIONS //
                    if(kAbs>kIRloc && kAbs<=kUVloc){
                        
                        // SINGLE INCLUSIVE BUFFER //
                        COMPLEX SIloc=SymmetricSingleInclusive(kXIndex,kYIndex);
                        
                        COMPLEX AsSIloc;
                        
                        // ASSIGN IF DESIGNATED //
                        if(AS_FLAG==1){
                            AsSIloc=AntiSymmetricSingleInclusive(kXIndex,kYIndex);
                        }
                        else{
                            AsSIloc=COMPLEX(0.0,0.0);
                        }
                        
                        // COMPUTE Phi_n's -- NEGLECT FACTORS OF dphi AND 2*Pi WHICH CANCEL IN RATIOS //
                        SI+=real(SIloc);
                        ASI+=real(AsSIloc);
                        
                        // SET VALUES //
                        F->SetP(kXIndex,kYIndex,0,real(SIloc));
                        F->SetP(kXIndex,kYIndex,1,AsSIloc*exp(1.0*ComplexI*Phik));
                        F->SetP(kXIndex,kYIndex,2,  SIloc*exp(2.0*ComplexI*Phik));
                        F->SetP(kXIndex,kYIndex,3,AsSIloc*exp(3.0*ComplexI*Phik));
                        F->SetP(kXIndex,kYIndex,4,  SIloc*exp(4.0*ComplexI*Phik));

                        // INCREMENT COUNT //
                        Ncount++;
                        
                    } // END k STEP FUNCTION IF STATEMENT //
                    else{
                        
                        F->SetP(kXIndex,kYIndex,0,COMPLEX(0.0));
                        F->SetP(kXIndex,kYIndex,1,COMPLEX(0.0));
                        F->SetP(kXIndex,kYIndex,2,COMPLEX(0.0));
                        F->SetP(kXIndex,kYIndex,3,COMPLEX(0.0));
                        F->SetP(kXIndex,kYIndex,4,COMPLEX(0.0));
                    
                    }
                    
                } // END kYIndex LOOP //
            } // END kXIndex LOOP //
            
            // PERFORM FFT //
            F->ExecutePtoX();
            
            // NORMALIZED VALUES FROM FOURIER TRANSFORM AT x=0 --> ALL p INTEGRATION //
            DOUBLE  Psi0=real(F->GetX(0,0,0)/DOUBLE(Lattice::N[0]*Lattice::N[1]));
            COMPLEX Psi1=F->GetX(0,0,1)/COMPLEX(Lattice::N[0]*Lattice::N[1]);
            COMPLEX Psi2=F->GetX(0,0,2)/COMPLEX(Lattice::N[0]*Lattice::N[1]);
            COMPLEX Psi3=F->GetX(0,0,3)/COMPLEX(Lattice::N[0]*Lattice::N[1]);
            COMPLEX Psi4=F->GetX(0,0,4)/COMPLEX(Lattice::N[0]*Lattice::N[1]);
            
            // OUTPUT //
            //  0 -- kIR, 1 -- kUV, 2 -- Ft0, 3 -- Re(Psi1), 4 -- Im(Psi1), 5 -- Re(Psi2), 6 -- Im(Psi2), 7 -- Re(Psi3), 8 -- Im(Psi3), 9 -- Re(Psi4), 10 -- Im(Psi4), 11 -- SI, 12 -- ASI,
            SIOutStream << kIRloc << " " << kUVloc << " " << Psi0 << " " << real(Psi1) << " " << imag(Psi1) << " " << real(Psi2) << " " << imag(Psi2) << " " << real(Psi3) << " " << imag(Psi3) << " " << real(Psi4) << " " << imag(Psi4) << " " << SI << " " << ASI << std::endl;

            
            // CLEAN-UP //
            delete F;
            
        }// END k LIMIT LOOP //
        
        // CLOSE OUTSTREAM //
        SIOutStream.close();
                
    }

    void DetermineIntegratedDistribution(){
        
        std::cerr << "# CONSTRUCTING INTEGRATED DISTRIBUTION" << std::endl;
        
        // COMPUTE FOURIER TRANSFORM OF OMEGAS AND DETERMINE SINGLE PARTICLE SPECTRA AND STORE TO FILE//
        DetermineIntegratedDistribution(OmegaS::O,OmegaA::O);
        
        std::cerr << "# FINISHED CONSTRUCTING INTEGRATED DISTRIBUTION" << std::endl;
        
    }
    
}


#endif
