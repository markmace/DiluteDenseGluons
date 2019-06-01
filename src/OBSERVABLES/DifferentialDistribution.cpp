// DifferentialDistribution.cpp is part of the Dilute-Dense Gluon solver //
// Copyright (C) 2019 Mark Mace //

#ifndef __DIFFERENTIAL_DISTRIBUTION__CPP__
#define __DIFFERENTIAL_DISTRIBUTION__CPP__

namespace Observables{

    // SINGLE DIFFERENTIAL SPECTRA //
    void DetermineDifferentialDistribution(CVectorFields *OmS,CVectorFields *OmA){

        // CREATE OUTPUT STREAM FOR DIFFERENTIAL DISTRIBUTION //
        std::ofstream SIOutStream; SIOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"DifferentialSpectrumID",MY_MPI_RNG_SEED,".txt").c_str());
        
        INT kMaxIndex=INT((kMax-kMin)/dkStep+0.1);
        
        DOUBLE v0ALL=0.0;
        
        COMPLEX v1ALL={0.0,0.0};
        COMPLEX v2ALL={0.0,0.0};
        COMPLEX v3ALL={0.0,0.0};
        COMPLEX v4ALL={0.0,0.0};
        
        INT NcountALL=0;

        //////////////////////
        // REF MOMENTUM BIN //
        //////////////////////
        
        // LOOP OVER ALL LATTICE MOMENTA FOR ANGULAR INTEGRATION -- IGNORING DOUBLERS //
        for(INT k2XIndex=0;k2XIndex<Lattice::N[0];k2XIndex++){
            if((k2XIndex<Lattice::N[0]/4+1)||(k2XIndex>3*Lattice::N[0]/4-1)){
                for(INT k2YIndex=0;k2YIndex<Lattice::N[1];k2YIndex++){
                    if((k2YIndex<Lattice::N[1]/4+1)||(k2YIndex>3*Lattice::N[0]/4-1)){
                        
                        // MOMENTUM VALUES //
                        DOUBLE k2XValue=IndexTok(k2XIndex);
                        DOUBLE k2YValue=IndexTok(k2YIndex);
                        
                        // DEFINE k^2 AND k //
                        DOUBLE k2Sqr=IndexTokSqr(k2XIndex,k2YIndex);
                        DOUBLE k2Abs=sqrt(k2Sqr);
                        
                        if(k2Abs>=kRefMin && k2Abs<=kRefMax){
                            
                            // SINGLE INCLUSIVE BUFFER //
                            COMPLEX SIloc=SymmetricSingleInclusive(k2XIndex,k2YIndex);
                            
                            COMPLEX AsSIloc;
                            
                            // ASSIGN IF DESIGNATED //
                            if(AS_FLAG==1){
                                AsSIloc=AntiSymmetricSingleInclusive(k2XIndex,k2YIndex);
                            }
                            else{
                                AsSIloc=COMPLEX(0.0,0.0);
                            }
                            
                            // DEFINE AZIMUTHAL ANGLE //
                            DOUBLE Phik2=atan2(k2YValue,k2XValue);
                            
                            if(kDIF_FLAG==1){
                                v0ALL+=real(k2Abs*dkBin*SIloc);
                                v1ALL+=k2Abs*dkBin*AsSIloc*exp(1.0*ComplexI*Phik2);
                                v2ALL+=k2Abs*dkBin*SIloc*exp(2.0*ComplexI*Phik2);
                                v3ALL+=k2Abs*dkBin*AsSIloc*exp(3.0*ComplexI*Phik2);
                                v4ALL+=k2Abs*dkBin*SIloc*exp(4.0*ComplexI*Phik2);
                            }
                            else if(kDIF_FLAG==2){
                                v0ALL+=real(k2Abs*SIloc);
                                v1ALL+=k2Abs*AsSIloc*exp(1.0*ComplexI*Phik2);
                                v2ALL+=k2Abs*SIloc*exp(2.0*ComplexI*Phik2);
                                v3ALL+=k2Abs*AsSIloc*exp(3.0*ComplexI*Phik2);
                                v4ALL+=k2Abs*SIloc*exp(4.0*ComplexI*Phik2);
                            }
                            else if(kDIF_FLAG==3){
                                v0ALL+=real(SIloc);
                                v1ALL+=AsSIloc*exp(1.0*ComplexI*Phik2);
                                v2ALL+=SIloc*exp(2.0*ComplexI*Phik2);
                                v3ALL+=AsSIloc*exp(3.0*ComplexI*Phik2);
                                v4ALL+=SIloc*exp(4.0*ComplexI*Phik2);
                            }
                            else{
                                std::cerr << "# NO K-INTEGRATION SCHEME DENOTED!!!" << std::endl;
                                exit(0);
                            }
                            // INCREMENT COUNT //
                            NcountALL++;
                        }
                    
                    }
                }
            }
        }

        // LOOP THROUGH MOMENTA OF INTEREST //
        for(INT nk=0; nk<=kMaxIndex; nk++){

            // MOMENTUM BUFFER -- GeV //
            DOUBLE K=kMin+dkStep*nk;

            std::cerr << "# K=" << K << " GeV" << std::endl;

            // BUFFER FOR ANGULAR INTEGRATION VARIABLES //
            DOUBLE SI=0.0;
            DOUBLE ASI=0.0;
            DOUBLE v0=0.0;
            
            COMPLEX v1={0.0,0.0};
            COMPLEX v2={0.0,0.0};
            COMPLEX v3={0.0,0.0};
            COMPLEX v4={0.0,0.0};

        
            // BIN COUNTER //
            INT Ncount=0;

            // LOOP OVER ALL LATTICE MOMENTA FOR ANGULAR INTEGRATION -- RESTRICTING FOR CERTAIN KEY BANDS //
            for(INT k2XIndex=0;k2XIndex<Lattice::N[0];k2XIndex++){
                if((k2XIndex<Lattice::N[0]/4+1)||(k2XIndex>3*Lattice::N[0]/4-1)){
                    for(INT k2YIndex=0;k2YIndex<Lattice::N[1];k2YIndex++){
                        if((k2YIndex<Lattice::N[1]/4+1)||(k2YIndex>3*Lattice::N[0]/4-1)){

                            // MOMENTUM VALUES //
                            DOUBLE k2XValue=IndexTok(k2XIndex);
                            DOUBLE k2YValue=IndexTok(k2YIndex);

                            // DEFINE k^2 AND k //
                            DOUBLE k2Sqr=IndexTokSqr(k2XIndex,k2YIndex);
                            DOUBLE k2Abs=sqrt(k2Sqr);
                            
                            // BIN BY MOMENTUM SPACING //
                            if((k2Abs-(K-0.5*dkBin))*(k2Abs-(K+0.5*dkBin))<0.0){
                                
                                // SINGLE INCLUSIVE BUFFER //
                                COMPLEX SIloc=SymmetricSingleInclusive(k2XIndex,k2YIndex);
                                
                                COMPLEX AsSIloc;
                                
                                // ASSIGN IF DESIGNATED //
                                if(AS_FLAG==1){
                                    AsSIloc=AntiSymmetricSingleInclusive(k2XIndex,k2YIndex);
                                }
                                else{
                                    AsSIloc=COMPLEX(0.0,0.0);
                                }
                                
                                // DEFINE AZIMUTHAL ANGLE //
                                DOUBLE Phik2=atan2(k2YValue,k2XValue);


                                if(kINT_FLAG==0){
                                    v0+=real(k2Abs*dkBin*SIloc);
                                    v1+=k2Abs*dkBin*AsSIloc*exp(1.0*ComplexI*Phik2);
                                    v2+=k2Abs*dkBin*SIloc*exp(2.0*ComplexI*Phik2);
                                    v3+=k2Abs*dkBin*AsSIloc*exp(3.0*ComplexI*Phik2);
                                    v4+=k2Abs*dkBin*SIloc*exp(4.0*ComplexI*Phik2);
                                }
                                else if(kINT_FLAG==1){
                                    v0+=real(k2Abs*SIloc);
                                    v1+=k2Abs*AsSIloc*exp(1.0*ComplexI*Phik2);
                                    v2+=k2Abs*SIloc*exp(2.0*ComplexI*Phik2);
                                    v3+=k2Abs*AsSIloc*exp(3.0*ComplexI*Phik2);
                                    v4+=k2Abs*SIloc*exp(4.0*ComplexI*Phik2);
                                }
                                else if(kINT_FLAG==2){
                                    v0+=real(SIloc);
                                    v1+=AsSIloc*exp(1.0*ComplexI*Phik2);
                                    v2+=SIloc*exp(2.0*ComplexI*Phik2);
                                    v3+=AsSIloc*exp(3.0*ComplexI*Phik2);
                                    v4+=SIloc*exp(4.0*ComplexI*Phik2);
                                }
                                else{
                                    std::cerr << "# NO K-INTEGRATION SCHEME DENOTED!!!" << std::endl;
                                    exit(0);
                                }
                                // INCREMENT COUNT //
                                Ncount++;
                                

                            } // END BIN //

                            
                        } // END IF kY STATEMENT //
                    } // END kY LOOP //
                } // END IF kX STATEMENT //
            } // END kX LOOP //

            // OUTPUT //
            SIOutStream << K << " " << v0 << " " << real(v1)/v0 << " " << imag(v1)/v0 << " " << real(v2)/v0 << " " << imag(v2)/v0 << " " << real(v3)/v0 << " " << imag(v3)/v0  << " " << real(v4)/v0 << " " << imag(v4)/v0 << " " << v0ALL << " " << real(v1ALL)/v0ALL << " " << imag(v1ALL)/v0ALL << " " << real(v2ALL)/v0ALL << " " << imag(v2ALL)/v0ALL << " " << real(v3ALL)/v0ALL << " " << imag(v3ALL)/v0ALL  << " " << real(v4ALL)/v0ALL << " " << imag(v4ALL)/v0ALL << " " << SI/Ncount << " " << ASI/Ncount << " " << Ncount << std::endl;
            

        }

        // CLOSE OUTSTREAM //
        SIOutStream.close();

    }

    void DetermineDifferentialDistribution(){

        std::cerr << "## CONSTRUCTING DIFFERENTIAL DISTRIBUTION" << std::endl;

        // COMPUTE FOURIER TRANSFORM OF OMEGAS AND DETERMINE SINGLE PARTICLE SPECTRA AND STORE TO FILE//
        DetermineDifferentialDistribution(OmegaS::O,OmegaA::O);

        std::cerr << "## FINISHED CONSTRUCTING DIFFERENTIAL DISTRIBUTION" << std::endl;

    }

}


#endif
