// CommandLineArguments.cpp is part of the Dilute-Dense Gluon solver //
// Copyright (C) 2019 Mark Mace //

#ifndef __COMMANDLINEARGUMENTS_CPP__
#define __COMMANDLINEARGUMENTS_CPP__

///////////////////////////////////
// PROCESS COMMANDLINE ARGUMENTS //
///////////////////////////////////

void ProcessCommandlineArguments(int argc,char **argv){

    Konfig arguments(argc,argv);

    // SET OUTPUT DIRECTORY //
    char OutDir[512]="OUTPUT";

    arguments.Getval("o",OutDir); // OUTPUT DIRECTORY //

    IO::SetOutputDirectory(OutDir);

    std::cerr << "# OUTPUT DIRECTORY IS " << IO::OutputDirectory <<  std::endl;

    INT NSites=-1;

    arguments.Getval("N",NSites); // SET NUMBER OF SITES IN TRANSVERSE DIRECTION SIMULTANEOUSLY //

    if(NSites>0){
        Lattice::N[0]=NSites;
        Lattice::N[1]=NSites;
        Lattice::Area=NSites*NSites;
        std::cerr << "# LATTICE TRANSVERSE SIZE IS " << Lattice::N[0] << "x" << Lattice::N[1] << std::endl;
    }
    else{
        std::cerr << "# NUMBER OF SITES NOT SPECIFIED -- USING " << Lattice::N[0] << "x" << Lattice::N[1] << std::endl;
    }

    DOUBLE aSpacing=-1.0;

    arguments.Getval("a",aSpacing);  // SET SPACING OF SITES IN TRANSVERSE DIRECTION SIMULTANEOUSLY //

    if(aSpacing>0.0){

        Lattice::a[0]=aSpacing;
        Lattice::a[1]=aSpacing;

        std::cerr << "# LATTICE SPACINGS ARE ax=" << Lattice::a[0] << " ay=" << Lattice::a[1] <<  std::endl;
    }
    else{
        std::cerr << "# LATTICE SPACING NOT SPECIFIED -- USING " << Lattice::a[0] << " " << Lattice::a[1] <<  std::endl;
    }

    Lattice::SizeX=Lattice::a[0]*Lattice::N[0];
    Lattice::SizeY=Lattice::a[1]*Lattice::N[1];

    std::cerr << "# SPATIAL EXTENT IS " << Lattice::SizeX << " x " << Lattice::SizeY << std::endl;

    INT NSitesRap=-1;

    arguments.Getval("NRap",NSitesRap); // NUMBER OF SITES IN THE RAPIDITY DIRECTION //

    if(NSitesRap>0){
        Lattice::NRap=NSitesRap;
        std::cerr << "# NUMBER OF SLICES IN RAPIDITY=" << Lattice::NRap << std::endl;
    }
    else{
        std::cerr << "# NUMBER OF SLICES IN RAPIDITY NOT SPECIFIED -- USING " << Lattice::NRap << std::endl;
    }
    
    // REGULATOR MASS //
    arguments.Getval("RegMass",RegMass);
    std::cerr << "# RegMass=" << RegMass << std::endl;
    
    // MOMENTUM VARIABLES FOR HARMONICS //
    arguments.Getval("kMin",kMin); // MINIMUM DIFFERENTIAL TRIGGER MOMEMNTA //
    arguments.Getval("kMax",kMax); // MAXIMUM DIFFERENTIAL TRIGGER MOMEMNTA //
    arguments.Getval("kRefMin",kRefMin); // MINIMUM DIFFERENTIAL REFERENCE MOMEMNTA //
    arguments.Getval("kRefMax",kRefMax); // MAXIMUM DIFFERENTIAL REFERENCE MOMEMNTA //
    arguments.Getval("dkBin",dkBin); // TRIGGER BIN SIZE IN MOMENTA //
    arguments.Getval("dkStep",dkStep); // TRIGGER STEP SIZE IN MOMENTA //
    arguments.Getval("kDIF_FLAG",kDIF_FLAG); // DETERMINE DIFF HARMONIC CALCULATION : 0 -- OFF, 1 -- \phi dk k , 2 -- \phi k //
    arguments.Getval("kINT_FLAG",kINT_FLAG); // DETERMINE INTEGRATED HARMONIC CALCULATION : 0 -- OFF, 1 -- ON //

    std::cerr << "# kMin=" << kMin << " kMax=" << kMax << " dkStep=" << dkStep << " dkBin=" << dkBin << std::endl;
    std::cerr << "# kRefMin=" << kRefMin << " kRefMax=" << kRefMax << " kDIF_FLAG=" << kDIF_FLAG << " kINT_FLAG=" << kINT_FLAG << std::endl;

    arguments.Getval("AS",AS_FLAG); // CALCULATION OF ANTI-SYMMETRIC SPECTRA (NEEDED FOR ODD HARMONICS) : 0 -- NOT CALCULATED, 1 -- CALCULATED //

    std::cerr << "# MEASUREMENT OF ANTISYMM=" << AS_FLAG << " (0 -- OFF, 1 -- ON)" << std::endl;

    // GLAUBER--IP-SAT SPECIFIC ARGUEMNTS //
    #if IC_FLAG==GIPS_FLAG
    arguments.Getval("OUT",OUTPUT_FLAG); // OUTPUT OF CONFIGURATIONS //
    std::cerr << "# OUTPUT OF CONFIGURATIONS=" << OUTPUT_FLAG << " (0 -- OFF, 1 -- ON, 2 -- ON +RHOS)" << std::endl;
    
    char InFile[512]="input";
    arguments.Getval("if",InFile); // INPUT PARAMETER FILE FOR GLAUBER-IPSAT PORTION OF CODE //
    IO::SetInputFile(InFile);
    std::cerr << "# INPUT FILE IS " << IO::InputFile <<  std::endl;
    #endif
    
    // MV SPECIFIC ARGUEMNTS //
    #if IC_FLAG==MV_FLAG
    arguments.Getval("OUT",OUTPUT_FLAG); // OUTPUT OF CONFIGURATIONS //
    std::cerr << "# OUTPUT OF CONFIGURATIONS=" << OUTPUT_FLAG << " (0 -- OFF, 1 -- ON, 2 -- ON +RHOS)" << std::endl;

    arguments.Getval("PROJ",PROJ_FLAG); // OUTPUT OF CONFIGURATIONS //
    std::cerr << "# PROJ_FLAG=" << PROJ_FLAG << " (0 -- p, 1 -- d, 2 -- INFINITE)" << std::endl;

    // SET COLOR CHARGE DENSITIES SCALES //
    arguments.Getval("g2muP",g2muP); // PROJECTILE COLOR CHARGE DENSITY //
    arguments.Getval("g2muT",g2muT); // TARGET COLOR CHARGE DENSITY //

    std::cerr << "# UNIFORM g2muP=" << g2muP  << " g2muT=" << g2muT << std::endl;
    #endif
}
#endif
