// FullDistribution.cpp is part of the Dilute-Dense Gluon solver //
// Copyright (C) 2019 Mark Mace //

// NEEDS TO BE CHECKED //

#ifndef __FULL_DISTRIBUTION__CPP__
#define __FULL_DISTRIBUTION__CPP__

namespace Observables{

    // SINGLE DIFFERENTIAL SPECTRA //
    void DetermineFullDistribution(CVectorFields *OmS,CVectorFields *OmA){

        // CREATE OUTPUT STREAM FOR DISTRIBUTION //
        std::ofstream TwoDimDistribution; TwoDimDistribution.open(StringManipulation::StringCast(IO::OutputDirectory,"TwoDimSIDistributionID",MY_MPI_RNG_SEED,".txt").c_str());

        DOUBLE MAXP=2.0*sqrt(2.0)*std::pow(1.0/(Lattice::a[0]*Lattice::a[1]),1.0/2.0)*invfmtoGeV;
        
        INT NumberOfBinsP=MAXP*16;
        
        std::cerr << "# MAXP=" << MAXP << " #NumberOfBinsP " << NumberOfBinsP << std::endl;
        
        Histogram *DistributionHistogram=new Histogram(0.0,MAXP,NumberOfBinsP);

        // SET SINGLE INCLUSIVE TO FOURIER SCALAR FIELD FOR MULTIPLICITY //
        for(INT kXIndex=0;kXIndex<Lattice::N[0];kXIndex++){
            for(INT kYIndex=0;kYIndex<Lattice::N[1];kYIndex++){

                // DEFINE kX,kY,k^2 //
                DOUBLE kXLOC=IndexTok(kXIndex); // GeV //
                DOUBLE kYLOC=IndexTok(kYIndex);
                DOUBLE kSqr=IndexTokSqr(kXIndex,kYIndex); // GeV^2 //
                DOUBLE kAbs=std::sqrt(kSqr); // GeV //
                
                // BEGIN OPTION TO REMOVE ZERO MODE EXPLICITLY //
                
                if(kSqr!=0){

                    // SINGLE INCLUSIVE BUFFER //
                    COMPLEX SI={0.0,0.0};

                    // SUM OVER ALEGBRA INDICES OF |OMEGA_S|^2+|OMEGA_A|^2 TO GET SINGLE INCLUSIVE//
                    for(INT a=0;a<SUNcAlgebra::VectorSize;a++){

                        SI+=FourierSpace::OmegaS->GetP(BoundaryIndex(Lattice::N[0]-kXIndex),BoundaryIndex(Lattice::N[1]-kYIndex),a)*FourierSpace::OmegaS->GetP(kXIndex,kYIndex,a)+FourierSpace::OmegaA->GetP(BoundaryIndex(Lattice::N[0]-kXIndex),BoundaryIndex(Lattice::N[1]-kYIndex),a)*FourierSpace::OmegaA->GetP(kXIndex,kYIndex,a);
                    }

                    // NORMALIZE 1/(k^2)  //
                    SI*=1.0/(kSqr*pow(2.0*PI,3)); // GeV^-2 //

                    // SET FOURIER SPACE VARIABLE //
                    FourierSpace::Scalar->SetP(kXIndex,kYIndex,0,SI);
                    
                    // UPDATE HISTOGRAM //
                    DistributionHistogram->Count(kAbs,real(SI));
                    
                    TwoDimDistribution << kXLOC << " " << kYLOC << " " << real(SI) << std::endl;
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
                    
                    SI+=FourierSpace::OmegaS->GetP(BoundaryIndex(Lattice::N[0]-kXIndex),BoundaryIndex(Lattice::N[1]-kYIndex),a)*FourierSpace::OmegaS->GetP(kXIndex,kYIndex,a)+FourierSpace::OmegaA->GetP(BoundaryIndex(Lattice::N[0]-kXIndex),BoundaryIndex(Lattice::N[1]-kYIndex),a)*FourierSpace::OmegaA->GetP(kXIndex,kYIndex,a);
                }
                
                // NORMALIZE 1/(k^2)  //
                SI*=1.0/((kSqr+RegMass*RegMass)*pow(2.0*PI,3)); // GeV^-2 //

                // SET FOURIER SPACE VARIABLE //
                FourierSpace::Scalar->SetP(kXIndex,kYIndex,0,SI);
                
                // UPDATE HISTOGRAM //
                DistributionHistogram->Count(kAbs,real(SI));
                
                TwoDimDistribution << kXLOC << " " << kYLOC << " " << real(SI) << std::endl;
                */
                // END REGULATOR MASS OPTION //

            }
        }
        
        // SET OUTPUT FILE//
        std::string OutputFileN=StringManipulation::StringCast(IO::OutputDirectory,"SIDistributionID",MY_MPI_RNG_SEED,".txt");
        
        // SET OUTPUT HEADER //
        std::string HeaderMessageN=StringManipulation::StringCast("# kAbs ","dNd2p");
        
        // OUTPUT WILSON LOOP //
        DistributionHistogram->Output(HeaderMessageN,OutputFileN);
        
        // CLEAN-UP //
        delete DistributionHistogram;

        // FOURIER TRANSFORM
        FourierSpace::Scalar->ExecutePtoX();

        // CLOSE OUTPUT //
        TwoDimDistribution.close();

        // COMMANDLINE OUTPUT //
        //std::cerr << "# FINISHED CALCULATING DISTRIBUTION " << std::endl;

    }

    void DetermineFullDistribution(){
        
        std::cerr << "# CONSTRUCTING FULL DISTRIBUTION" << std::endl;

        // COMPUTE FOURIER TRANSFORM OF OMEGAS AND DETERMINE SINGLE PARTICLE SPECTRA AND STORE TO FILE//
        DetermineFullDistribution(OmegaS::O,OmegaA::O);

        std::cerr << "# FINISHED CONSTRUCTING FULL DISTRIBUTION" << std::endl;

    }

}


#endif
