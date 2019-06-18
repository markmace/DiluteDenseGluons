// MultiplicityMeasurement.cpp is part of the Dilute-Dense Gluon solver //
// Copyright (C) 2019 Mark Mace //

#ifndef __MULTIPLICITY_MEASREUMENT__CPP__
#define __MULTIPLICITY_MEASREUMENT__CPP__

namespace Observables{
    
    void MeasureMultiplicity(CVectorFields *OmS,CVectorFields *OmA){
        
        // COMMANDLINE OUTPUT //
        std::cerr << "# CALCULATING MULTIPLICITY" << std::endl;
        
        // CREATE OUTPUT STREAM FOR MULTIPLICITY //
        std::ofstream MultOutStream; MultOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"MultiplicityID",MY_MPI_RNG_SEED,".txt").c_str());
        
        // SET SINGLE INCLUSIVE TO FOURIER SCALAR FIELD FOR MULTIPLICITY //
        //#pragma omp for collapse(2)
        for(INT kXIndex=0;kXIndex<Lattice::N[0];kXIndex++){
            for(INT kYIndex=0;kYIndex<Lattice::N[1];kYIndex++){
                
                // DEFINE k^2 //
                DOUBLE kSqr=IndexTokSqr(kXIndex,kYIndex);
                
                // BEGIN OPTION TO REMOVE ZERO MODE EXPLICITLY //
                
                if(kSqr!=0){
                    
                    // SINGLE INCLUSIVE BUFFER //
                    COMPLEX SI={0.0,0.0};
                    
                    // SUM OVER ALEGBRA INDICES OF |OMEGA_S|^2+|OMEGA_A|^2 TO GET SINGLE INCLUSIVE//
                    for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                        
                        SI+=FourierSpace::OmegaS->GetP(BoundaryIndex(Lattice::N[0]-kXIndex),BoundaryIndex(Lattice::N[1]-kYIndex),a)*FourierSpace::OmegaS->GetP(kXIndex,kYIndex,a)+FourierSpace::OmegaA->GetP(BoundaryIndex(Lattice::N[0]-kXIndex),BoundaryIndex(Lattice::N[1]-kYIndex),a)*FourierSpace::OmegaA->GetP(kXIndex,kYIndex,a); // DIMLESS //
                    }
                    
                    // RE-NORMALIZE 1/(k^2 (2\pi)^3) //
                    SI*=1.0/(kSqr*pow(2.0*PI,3)); // GeV^-2 //
                    
                    // SET FOURIER SPACE VARIABLE //
                    FourierSpace::Scalar->SetP(kXIndex,kYIndex,0,SI);
                }
                
                // SET ZERO POINT TO ZERO //
                FourierSpace::Scalar->SetP(0,0,0,0);
                
                // END ZERO MODE REMOVAL OPTION //
                
                // BEGIN OPTION TO USE REGULATOR MASS //
                /*
                // SINGLE INCLUSIVE BUFFER //
                COMPLEX SI={0.0,0.0};
                
                // SUM OVER ALEGBRA INDICES OF |OMEGA_S|^2+|OMEGA_A|^2 TO GET SINGLE INCLUSIVE//
                for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                    
                    SI+=FourierSpace::OmegaS->GetP(BoundaryIndex(Lattice::N[0]-kXIndex),BoundaryIndex(Lattice::N[1]-kYIndex),a)*FourierSpace::OmegaS->GetP(kXIndex,kYIndex,a)+FourierSpace::OmegaA->GetP(BoundaryIndex(Lattice::N[0]-kXIndex),BoundaryIndex(Lattice::N[1]-kYIndex),a)*FourierSpace::OmegaA->GetP(kXIndex,kYIndex,a); // DIMLESS //
                }
                
                // RE-NORMALIZE 1/((k^2+m^2) (2\pi)^3) //
                SI*=1.0/((kSqr+RegMass*RegMass)*pow(2.0*PI,3)); // GeV^-2 //
                
                // SET FOURIER SPACE VARIABLE //
                FourierSpace::Scalar->SetP(kXIndex,kYIndex,0,SI);
                */
                // END REGULATOR MASS OPTION //

                
            }
            
            
        }
        
        // FOURIER TRANSFORM
        FourierSpace::Scalar->ExecutePtoX();
        
        // DEFINE MULTIPLICITY -- RE-NORMALIZE \int d^2p exp(0) dN/d^2p = MULT //
        DOUBLE Multiplicity=real(FourierSpace::Scalar->GetX(0,0,0))*(invfmtoGeV/Lattice::SizeX)*(invfmtoGeV/Lattice::SizeY);
        
        // OUTPUT MULTPLICITY AND IMPACT PARAMETER //
        MultOutStream << Multiplicity << " " << bImpact << std::endl;
        
        std::cerr << "# MULTIPLICTY=" << Multiplicity << " bImpact=" << bImpact << std::endl;
        
        // CLOSE OUTPUT //
        MultOutStream.close();
        
        // COMMANDLINE OUTPUT //
        std::cerr << "# FINISHED CALCULATING MULTIPLICITY" << std::endl;
        
    }
    
    void MeasureMultiplicity(){
        
        //std::cerr << "# MEASURING MULTIPLICITY" << std::endl;
        
        // COMPUTE MULTPLICITY //
        MeasureMultiplicity(OmegaS::O,OmegaA::O);
                
    }
}
#endif
