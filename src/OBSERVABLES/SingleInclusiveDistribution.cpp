// SingleInclusiveDistribution.cpp is part of the Dilute-Dense Gluon solver //
// Copyright (C) 2019 Mark Mace //

#ifndef __SINGLE_INCLUSIVE_DISTRIBUTION__CPP__
#define __SINGLE_INCLUSIVE_DISTRIBUTION__CPP__

namespace Observables{
    
    // ENFORCE PERIODIC BOUNDARIES //
    INT BoundaryIndex(INT i){
        if(i<0){return i+Lattice::N[0];}
        else if(i>Lattice::N[0]-1){return i-Lattice::N[0];}
        else{return i;}
    }
    
    INT SIGN_V(DOUBLE i){
        if(i<0){return -1;}
        else{return 1;}
    }
    
    // GENERAL PURPOSE FUNCTIONS FOR CONVERTING INDICES TO MOMENTA //
    DOUBLE IndexTok(INT Index){
        DOUBLE karg=2.0*PI*Index/(Lattice::N[0]*Lattice::a[0]);
        DOUBLE kValue=1.0/Lattice::a[0]*sin(karg*Lattice::a[0])*invfmtoGeV;
        return kValue;
    }
    
    DOUBLE IndexTokSqr(INT kXIndex,INT kYIndex){
        DOUBLE kXarg=2.0*PI*kXIndex/(Lattice::N[0]*Lattice::a[0]);
        DOUBLE kYarg=2.0*PI*kYIndex/(Lattice::N[1]*Lattice::a[1]);
        return (pow(2.0/Lattice::a[0]*sin(kXarg*Lattice::a[0]/2.0),2)+pow(2.0/Lattice::a[0]*sin(kYarg*Lattice::a[1]/2.0),2))*invfmtoGeV*invfmtoGeV;
    }

    
    // CONSTRUCT OMEGA_S AND OMEGA_A //
    void ConstructOmegas(CVectorFields *OmS,CVectorFields *OmA){
        
        // DETERMINE DERIVATES AND OMEGAS
        #pragma omp for collapse(2)
        for(INT x=0;x<Lattice::N[0];x++){
            for(INT y=0;y<Lattice::N[1];y++){
                
                SU_Nc_FUNDAMENTAL_FORMAT OmegaSBuffer[SUNcGroup::MatrixSize]={SU_Nc_FUNDAMENTAL_FORMAT(0.0),SU_Nc_FUNDAMENTAL_FORMAT(0.0),SU_Nc_FUNDAMENTAL_FORMAT(0.0),SU_Nc_FUNDAMENTAL_FORMAT(0.0),SU_Nc_FUNDAMENTAL_FORMAT(0.0),SU_Nc_FUNDAMENTAL_FORMAT(0.0),SU_Nc_FUNDAMENTAL_FORMAT(0.0),SU_Nc_FUNDAMENTAL_FORMAT(0.0),SU_Nc_FUNDAMENTAL_FORMAT(0.0)};
                SU_Nc_FUNDAMENTAL_FORMAT OmegaABuffer[SUNcGroup::MatrixSize]={SU_Nc_FUNDAMENTAL_FORMAT(0.0),SU_Nc_FUNDAMENTAL_FORMAT(0.0),SU_Nc_FUNDAMENTAL_FORMAT(0.0),SU_Nc_FUNDAMENTAL_FORMAT(0.0),SU_Nc_FUNDAMENTAL_FORMAT(0.0),SU_Nc_FUNDAMENTAL_FORMAT(0.0),SU_Nc_FUNDAMENTAL_FORMAT(0.0),SU_Nc_FUNDAMENTAL_FORMAT(0.0),SU_Nc_FUNDAMENTAL_FORMAT(0.0)};
                
                COMPLEX CHECK=COMPLEX(0.0,0.0);
                
                for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                    
                    // PROJECTILE DERIVATIVES //
                    double dUdxProj=0.5*(ProjSolution::A->Get(BoundaryIndex(x+1),y,0,a)[0]-ProjSolution::A->Get(BoundaryIndex(x-1),y,0,a)[0])/Lattice::a[0];
                    double dUdyProj=0.5*(ProjSolution::A->Get(x,BoundaryIndex(y+1),0,a)[0]-ProjSolution::A->Get(x,BoundaryIndex(y-1),0,a)[0])/Lattice::a[1];
                    
                    SU_Nc_FUNDAMENTAL_FORMAT Wx[SUNcGroup::MatrixSize];
                    SU_Nc_FUNDAMENTAL_FORMAT Wy[SUNcGroup::MatrixSize];
                    
                    // BUFFER //
                    SU_Nc_FUNDAMENTAL_FORMAT Buffer[SUNcGroup::MatrixSize];
                    
                    // LAMBDA BUFFER //
                    SU_Nc_FUNDAMENTAL_FORMAT Generator[SUNcGroup::MatrixSize];
                    COPY_SUNcMatrix(Generator,SUNcAlgebra::GeneratorTa[a]);
                    
                    // +x U^dag t^a U //
                    SU_Nc_FUNDAMENTAL_FORMAT tmpPlus[SUNcGroup::MatrixSize];
                    COPY_SUNcMatrix(tmpPlus,TargetFields::U->Get(BoundaryIndex(x+1),y,0));
                    
                    SUNcGroup::Operations::UU(Generator,tmpPlus,Buffer);
                    
                    SUNcGroup::Operations::DU(TargetFields::U->Get(BoundaryIndex(x+1),y,0),Buffer,tmpPlus);
                    
                    // -x U^dag t^a U //
                    SU_Nc_FUNDAMENTAL_FORMAT tmpMinus[SUNcGroup::MatrixSize];
                    COPY_SUNcMatrix(tmpMinus,TargetFields::U->Get(BoundaryIndex(x-1),y,0));
                    
                    SUNcGroup::Operations::UU(Generator,tmpMinus,Buffer);
                    
                    SUNcGroup::Operations::DU(TargetFields::U->Get(BoundaryIndex(x-1),y,0),Buffer,tmpMinus);
                    
                    for(INT alpha=0;alpha<Nc;alpha++){
                        for(INT beta=0;beta<Nc;beta++){
                            
                            Wx[SUNcGroup::MatrixIndex(alpha,beta)]=DOUBLE(0.5)/Lattice::a[0]*(tmpPlus[SUNcGroup::MatrixIndex(alpha,beta)]-tmpMinus[SUNcGroup::MatrixIndex(alpha,beta)]);
                        }
                    }
                    
                    // +x U^dag t^a U //
                    COPY_SUNcMatrix(tmpPlus,TargetFields::U->Get(x,BoundaryIndex(y+1),0));
                    
                    SUNcGroup::Operations::UU(Generator,tmpPlus,Buffer);
                    
                    SUNcGroup::Operations::DU(TargetFields::U->Get(x,BoundaryIndex(y+1),0),Buffer,tmpPlus);
                    
                    // -x U^dag t^a U //
                    COPY_SUNcMatrix(tmpMinus,TargetFields::U->Get(x,BoundaryIndex(y-1),0));
                    
                    SUNcGroup::Operations::UU(Generator,tmpMinus,Buffer);
                    
                    SUNcGroup::Operations::DU(TargetFields::U->Get(x,BoundaryIndex(y-1),0),Buffer,tmpMinus);
                    
                    for(INT alpha=0;alpha<Nc;alpha++){
                        for(INT beta=0;beta<Nc;beta++){
                            
                            Wy[SUNcGroup::MatrixIndex(alpha,beta)]=DOUBLE(0.5)/Lattice::a[0]*(tmpPlus[SUNcGroup::MatrixIndex(alpha,beta)]-tmpMinus[SUNcGroup::MatrixIndex(alpha,beta)]);
                        }
                    }
                    
                    CHECK+=(OmegaSBuffer[SUNcGroup::MatrixIndex(0,0)]);
                    
                    // SET OMEGA MATRICES //
                    for(INT alpha=0;alpha<Nc;alpha++){
                        for(INT beta=0;beta<Nc;beta++){
                            OmegaSBuffer[SUNcGroup::MatrixIndex(alpha,beta)]+=dUdxProj*Wx[SUNcGroup::MatrixIndex(alpha,beta)]+dUdyProj*Wy[SUNcGroup::MatrixIndex(alpha,beta)];
                            OmegaABuffer[SUNcGroup::MatrixIndex(alpha,beta)]+=dUdxProj*Wy[SUNcGroup::MatrixIndex(alpha,beta)]-dUdyProj*Wx[SUNcGroup::MatrixIndex(alpha,beta)];
                        }
                    }
                    
                } // END VECTORSIZE a LOOP//
                
                COMPLEX OmegaSgenBuffer[SUNcGroup::MatrixSize];
                COMPLEX OmegaAgenBuffer[SUNcGroup::MatrixSize];
                
                SUNcGroup::Operations::GeneratorProjection(OmegaSBuffer,OmegaSgenBuffer);
                SUNcGroup::Operations::GeneratorProjection(OmegaABuffer,OmegaAgenBuffer);
                
                // SET TO VECTOR //
                for(INT a=0;a<SUNcGroup::MatrixSize;a++){
                    
                    OmS->Get(x,y,0,a)[0]=OmegaSgenBuffer[a];
                    OmA->Get(x,y,0,a)[0]=OmegaAgenBuffer[a];
                    
                }
                
            } // END y LOOP //
        } // END x LOOP //
        
    }
    
    void FourierTransformOmegas(CVectorFields *OmS,CVectorFields *OmA){
        
        // INITIALIZE FOURIER TRANSFORM VARIABLES //
        #pragma omp for collapse(2)
        for(INT y=0;y<Lattice::N[1];y++){
            for(INT x=0;x<Lattice::N[0];x++){
                for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                    FourierSpace::OmegaS->SetXc(x,y,a,OmS->Get(x,y,0,a)[0]*Lattice::a[0]*Lattice::a[1]);
                    FourierSpace::OmegaA->SetXc(x,y,a,OmA->Get(x,y,0,a)[0]*Lattice::a[0]*Lattice::a[1]);
                    
                }
            }
        }
        
        // COMPLETE FOURIER TRANSFORMs //
        FourierSpace::OmegaS->ExecuteXtoP();
        FourierSpace::OmegaA->ExecuteXtoP();
    }
    
    COMPLEX SymmetricSingleInclusive(INT kXIndex,INT kYIndex){
        
        COMPLEX SIloc=COMPLEX(0.0,0.0);
        
        // SUM OVER ALEGBRA INDICES OF |OMEGA_S|^2+|OMEGA_A|^2 TO GET SINGLE INCLUSIVE//
        for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
            
            SIloc+=FourierSpace::OmegaS->GetP(BoundaryIndex(Lattice::N[0]-kXIndex),BoundaryIndex(Lattice::N[1]-kYIndex),a)*FourierSpace::OmegaS->GetP(kXIndex,kYIndex,a)+FourierSpace::OmegaA->GetP(BoundaryIndex(Lattice::N[0]-kXIndex),BoundaryIndex(Lattice::N[1]-kYIndex),a)*FourierSpace::OmegaA->GetP(kXIndex,kYIndex,a);
        }
        
        // RE-NORMALIZE 1/(k^2 (2\pi)^3 L^2) //
        DOUBLE kSqr=IndexTokSqr(kXIndex,kYIndex);
        SIloc*=1.0/((kSqr+1e-9)*pow(2.0*PI,3)*Lattice::SizeX*Lattice::SizeY);
        
        return SIloc;
        
    }
    
    // COMPLETE ANTISYMMETRIC SINGLE INCLUSIVE -- BASED ON arXiv:1802.08166 EQN D5 //
    COMPLEX AntiSymmetricSingleInclusive(INT kXIndex, INT kYIndex){
        
        // CONSTANTS //
        // UV CONSTANT TO REMOVE DIVERGENCE //
        double UV=Lattice::a[0]*Lattice::a[1]*1e-5;
        
        DOUBLE kXValue=IndexTok(kXIndex);
        DOUBLE kYValue=IndexTok(kYIndex);
        DOUBLE kSqr=IndexTokSqr(kXIndex,kYIndex);
        
        COMPLEX Integrate=COMPLEX(0.0,0.0);
        
        // q-INTEGRATION //
        for(INT qXIndex=0;qXIndex<Lattice::N[0];qXIndex++){
            for(INT qYIndex=0;qYIndex<Lattice::N[1];qYIndex++){
                
                
                INT kXMinusqXIndex=BoundaryIndex(kXIndex-qXIndex);
                INT kYMinusqYIndex=BoundaryIndex(kYIndex-qYIndex);
                // q_x, q_y //
                DOUBLE qXValue=IndexTok(qXIndex);
                DOUBLE qYValue=IndexTok(qYIndex);
                // q^2 //
                DOUBLE qSqr=IndexTokSqr(qXIndex,qYIndex);
                // (k-q)_i //
                DOUBLE kXMinusqX=IndexTok(kXMinusqXIndex);
                DOUBLE kYMinusqY=IndexTok(kYMinusqYIndex);
                // k.(k-q) //
                DOUBLE kDotkMinusq=kXValue*kXMinusqX+kYValue*kYMinusqY;
                // q.(k-q) //
                DOUBLE qDotkMinusq=qXValue*kXMinusqX+qYValue*kYMinusqY;
                // (k-q)^2
                DOUBLE kMinusqSqr=IndexTokSqr(kXMinusqXIndex,kYMinusqYIndex);
                
                // CROSS PRODUCT ///
                // CHECK THIS //
                DOUBLE kCrossq=kXMinusqX*qYValue-kYMinusqY*qXValue;
                
                COMPLEX Sum=COMPLEX(0.0,0.0);
                
                for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                    for(INT b=0;b<SUNcAlgebra::VectorSize;b++){
                        for(INT c=0;c<SUNcAlgebra::VectorSize;c++){
                            if(SUNcAlgebra::StructureFunctions::f0(a,b,c)*SUNcAlgebra::StructureFunctions::f0(a,b,c)>0){
                                
                                //INTEGRATION -- INHERENTLY GETTING IMAGINARY PART //
                                Sum+=SUNcAlgebra::StructureFunctions::f0(a,b,c)*SIGN_V(kCrossq)*1.0/((qSqr+UV)*(kMinusqSqr+UV))*(
                                                                                                                                 (kSqr*FourierSpace::OmegaA->GetP(qXIndex,qYIndex,a)*FourierSpace::OmegaA->GetP(kXMinusqXIndex,kYMinusqYIndex,b)*FourierSpace::OmegaA->GetP(BoundaryIndex(Lattice::N[0]-kXIndex),BoundaryIndex(Lattice::N[1]-kYIndex),c)
                                                                                                                                  -qDotkMinusq*FourierSpace::OmegaA->GetP(qXIndex,qYIndex,a)*FourierSpace::OmegaA->GetP(kXMinusqXIndex,kYMinusqYIndex,b)*FourierSpace::OmegaA->GetP(BoundaryIndex(Lattice::N[0]-kXIndex),BoundaryIndex(Lattice::N[1]-kYIndex),c)
                                                                                                                                  -qDotkMinusq*FourierSpace::OmegaS->GetP(qXIndex,qYIndex,a)*FourierSpace::OmegaS->GetP(kXMinusqXIndex,kYMinusqYIndex,b)*FourierSpace::OmegaA->GetP(BoundaryIndex(Lattice::N[0]-kXIndex),BoundaryIndex(Lattice::N[1]-kYIndex),c)
                                                                                                                                  +2.0*kDotkMinusq*FourierSpace::OmegaA->GetP(qXIndex,qYIndex,a)*FourierSpace::OmegaS->GetP(kXMinusqXIndex,kYMinusqYIndex,b)*FourierSpace::OmegaS->GetP(BoundaryIndex(Lattice::N[0]-kXIndex),BoundaryIndex(Lattice::N[1]-kYIndex),c))
                                                                                                                                 -(kSqr*FourierSpace::OmegaA->GetP(BoundaryIndex(Lattice::N[0]-qXIndex),BoundaryIndex(Lattice::N[1]-qYIndex),a)*FourierSpace::OmegaA->GetP(BoundaryIndex(Lattice::N[0]-kXMinusqXIndex),BoundaryIndex(Lattice::N[1]-kYMinusqYIndex),b)
                                                                                                                                   *FourierSpace::OmegaA->GetP(kXIndex,kYIndex,c)
                                                                                                                                   -qDotkMinusq*FourierSpace::OmegaA->GetP(BoundaryIndex(Lattice::N[0]-qXIndex),BoundaryIndex(Lattice::N[1]-qYIndex),a)*FourierSpace::OmegaA->GetP(BoundaryIndex(Lattice::N[0]-kXMinusqXIndex),BoundaryIndex(Lattice::N[1]-kYMinusqYIndex),b)*FourierSpace::OmegaA->GetP(kXIndex,kYIndex,c)
                                                                                                                                   -qDotkMinusq*FourierSpace::OmegaS->GetP(BoundaryIndex(Lattice::N[0]-qXIndex),BoundaryIndex(Lattice::N[1]-qYIndex),a)*FourierSpace::OmegaS->GetP(BoundaryIndex(Lattice::N[0]-kXMinusqXIndex),BoundaryIndex(Lattice::N[1]-kYMinusqYIndex),b)*FourierSpace::OmegaA->GetP(kXIndex,kYIndex,c)
                                                                                                                                   +2.0*kDotkMinusq*FourierSpace::OmegaA->GetP(BoundaryIndex(Lattice::N[0]-qXIndex),BoundaryIndex(Lattice::N[1]-qYIndex),a)*FourierSpace::OmegaS->GetP(BoundaryIndex(Lattice::N[0]-kXMinusqXIndex),BoundaryIndex(Lattice::N[1]-kYMinusqYIndex),b)*FourierSpace::OmegaS->GetP(kXIndex,kYIndex,c))
                                                                                                                                 );
                                
                            }// END IF |f|>0 STATEMENT //
                            
                        }// END c LOOP ///
                        //Integrate+=Sum;
                        
                    }// END b LOOP //
                }// END a LOOP //
                Integrate+=Sum;
            } // END qX LOOP //
        } // END qY LOOP //
        
        return ComplexI*0.5*Integrate/(kSqr*std::pow(Lattice::SizeX*Lattice::SizeY,2)*std::pow(2.0*PI,3)); // PER UNIT LENGTH //
        
    }
}

#endif
