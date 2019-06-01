// SaveConfigurations.cpp is part of the Dilute-Dense Gluon solver //
// Copyright (C) 2019 Mark Mace //

namespace IO{
    // UNDER CONSTRUCTION //
    /*
    void SaveConfiguration(std::string Ufname,std::string Efname,WilsonLines *U,VectorFields *E){
        
        // OUTPUT STREAMS //
        std::ofstream UOutStream,EOutStream;
        
        // OUTPUT FILES //
        std::string UOutputFile=StringManipulation::StringCast(IO::OutputDirectory,Ufname,"ID",RandomNumberGenerator::MySEED,".txt");
        std::string EOutputFile=StringManipulation::StringCast(IO::OutputDirectory,Efname,"ID",RandomNumberGenerator::MySEED,".txt");

        // OPEN FILES //
        UOutStream.open(UOutputFile.c_str()); EOutStream.open(EOutputFile.c_str());
        
        // SET PRECISION //
        UOutStream.precision(OUTPUT_PRECISION); EOutStream.precision(OUTPUT_PRECISION);
        
        // CREATE OUTPUT //
        for(INT z=0;z<U->N[2];z++){
            for(INT y=0;y<U->N[1];y++){
                for(INT x=0;x<U->N[0];x++){
                    for(INT mu=0;mu<Lattice::Dimension;mu++){
                        
                        // OUTPUT GAUGE LINKS //
                        UOutStream << x << " " << y << " " << z << " " << mu << " " << SUNcGroup::IO::MatrixToString(U->Get(x,y,z,mu)) << std::endl;
     
                        // BEGIN SAYANTAN OUTPUT //
                        
                        //COMPLEX UMat[Nc*Nc];
                        //SUNcGroup::Operations::GetMatrix(GLinks::U->Get(x,y,z,mu),UMat);
                        //
                        //UOutStream << std::real(UMat[0]) <<  " " << std::imag(UMat[0]) << std::endl;
                        //UOutStream << std::real(UMat[2]) <<  " " << std::imag(UMat[2]) << std::endl;
                        //UOutStream << std::real(UMat[1]) <<  " " << std::imag(UMat[1]) << std::endl;
                        //UOutStream << std::real(UMat[3]) <<  " " << std::imag(UMat[3]) << std::endl;
                        //// END SAYANTAN OUTPUT //
     
                        // OUTPUT ELECTRIC FIELDS //
                        EOutStream << x << " " << y << " " << z << " " << mu;
                        
                        for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                            EOutStream << " " << E->Get(x,y,z,mu,a)[0];
                        }
                        
                        EOutStream << std::endl;
                        
                    }
                    
                }
            }
        }
        
        // CLOSE OUTPUT STREAM //
        UOutStream.close(); EOutStream.close();
        
    }

    
    void SaveConfigurationU(std::string Ufname,WilsonLines *U){
        
        // OUTPUT STREAMS //
        std::ofstream UOutStream;
        
        // OUTPUT FILES //
        std::string UOutputFile=StringManipulation::StringCast(IO::OutputDirectory,Ufname,"ID",RandomNumberGenerator::MySEED,".txt");
        
        // OPEN FILES //
        UOutStream.open(UOutputFile.c_str());
        
        // SET PRECISION //
        UOutStream.precision(OUTPUT_PRECISION);
        
        // CREATE OUTPUT //
        for(INT z=0;z<U->N[2];z++){
            for(INT y=0;y<U->N[1];y++){
                for(INT x=0;x<U->N[0];x++){
                    for(INT mu=0;mu<Lattice::Dimension;mu++){
                        
                        // OUTPUT GAUGE LINKS //
                        UOutStream << x << " " << y << " " << z << " " << mu << " " << SUNcGroup::IO::MatrixToString(U->Get(x,y,z,mu)) << std::endl;
                        
                    }
                    
                }
            }
        }
        
        // CLOSE OUTPUT STREAM //
        UOutStream.close(); EOutStream.close();
        
    }
    */
    
    void SaveConfigurationE(std::string Efname,VectorFields *E){
        
        // OUTPUT STREAMS //
        std::ofstream EOutStream;
        
        // OUTPUT FILES //
        std::string EOutputFile=StringManipulation::StringCast(IO::OutputDirectory,Efname,"ID",RandomNumberGenerator::MySEED,".txt");
        
        // OPEN FILES //
        EOutStream.open(EOutputFile.c_str());
        
        // SET PRECISION //
        EOutStream.precision(OUTPUT_PRECISION);
        
        // CREATE OUTPUT //
        for(INT y=0;y<E->N[1];y++){
            for(INT x=0;x<E->N[0];x++){
                
                // OUTPUT ELECTRIC FIELDS //
                EOutStream << x << " " << y;
                
                EOutStream << " " << E->Get(x,y,0,0)[0];
                
                EOutStream << std::endl;
                                    
            }
        }
        
        // CLOSE OUTPUT STREAM //
        EOutStream.close();
        
    }

    
    
}
