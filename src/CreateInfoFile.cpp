// CreateInfoFile.cpp is part of the Dilute-Dense Gluon solver //
// Copyright (C) 2019 Mark Mace //

#ifndef __CREATE_INFO_FILE_CPP__
#define __CREATE_INFO_FILE_CPP__

void CreateInfoFile(INT SEED){
    
    // CREATE OUTPUT STREAM FOR INFO FILE //
    std::ofstream InfoStream; InfoStream.open(StringManipulation::StringCast(IO::OutputDirectory,"Info",SEED,".txt").c_str());
    
    InfoStream << "############ BEGIN INFO FILE ############" << std::endl;
    // LATTICE PARAMETERS //
    InfoStream << "#N=" << Lattice::N[0] << " " << Lattice::N[1] << std::endl;
    InfoStream << "#NRap=" << Lattice::NRap << std::endl;
    InfoStream << "#a=" << Lattice::a[0] << " " << Lattice::a[1] << std::endl;
    InfoStream << "#aeta=" << Lattice::aeta << std::endl;
    
    InfoStream << "#RegMass=" << RegMass << std::endl;

    InfoStream << "#kMin=" << kMin << std::endl;
    InfoStream << "#kMax=" << kMax << std::endl;
    InfoStream << "#kRefMin=" << kRefMin << std::endl;
    InfoStream << "#kRefMax=" << kRefMax << std::endl;
    InfoStream << "#dkBin=" << dkBin << std::endl;
    InfoStream << "#dkStep=" << dkStep << std::endl;
    InfoStream << "#kDIF_FLAG=" << kDIF_FLAG << std::endl;
    InfoStream << "#kINT_FLAG=" << kINT_FLAG << std::endl;

    #if IC_FLAG==GIPS_FLAG
    InfoStream << "#InputFile=" << IO::InputFile << std::endl;
    #endif
    
    #if IC_FLAG==MV_FLAG
    InfoStream << "#PROJ_FLAG=" << PROJ_FLAG << std::endl;
    InfoStream << "#g2muP=" << g2muP << " UNIFORM " << std::endl;
    InfoStream << "#g2muT=" << g2muT << " UNIFORM " << std::endl;
    
    #endif

    InfoStream << "############ END INFO FILE ############" << std::endl;
    
    // CLOSE OUTPUT //
    InfoStream.close();
}

#endif
