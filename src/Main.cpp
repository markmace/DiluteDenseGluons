///////////////////////////////////////////////////////////////////////////////////////////////////
// PROGRAM TO COMPUTE PARTICLE PRODUCTION IN NUCLEAR COLLISONS IN THE DILUTE DENSE CGC FRAMEWORK //
///////////////////////////////////////////////////////////////////////////////////////////////////
#define IC_FLAG MV_FLAG // OPTIONS : MV MODEL -- MV_FLAG, GLAUBER+IP-SAT -- GIPS_FLAG

#include <iostream>
#include <string>

// MOST BASIC DEFINITIONS //
#include "Definitions.cpp"

// INCLUDE MPI SAMPLING //
#include "MPI/Basic.cpp"

// SET GAUGE GROUP //
#define SU_Nc_FLAG SU3_FLAG

/////////////////////////////////////////
// SU(Nc) GROUP AND ALGEBRA OPERATIONS //
/////////////////////////////////////////

// INCLUDE SU(Nc) ALGEBRA AND GROUP DEFINITION //
#include "SUNc/Definitions.cpp"

// DEFINE 3D GRID //
#include "LATTICE/3DGrid.cpp"

// DEFINE WILSON LINE CLASS //
#include "LATTICE/WilsonLines.cpp"

// DEFINE VECTOR FIELDS CLASSES //
#include "LATTICE/VectorFields.cpp"
#include "LATTICE/ComplexVectorFields.cpp"

// DEFINE SPECIFIC WILSON LINE AND VECTOR FIELD VARIABLES //
#include "LATTICE/Fields.cpp"

// HISTOGRAM //
#include "MISC/HISTOGRAM/Histogram.cpp"
#include "MISC/HISTOGRAM/MultiHistogram.cpp"

///////////////////////////
// SIMULATION PARAMETERS //
///////////////////////////

// DEFINE GLOBAL IMPACT PARAMETER -- SET DYNAMICALLY //
DOUBLE bImpact; // [fm] //

// REGULATOR MASS -- CAN BE SET IN COMMANDLINE //
DOUBLE RegMass=0.1; // [GeV] -- //

// MOMENTUM VARIABLES -- ALL CAN BE SET IN COMMANDLINE  //
DOUBLE kMin=0.5; // DIFFERENTIAL kMin [GeV]
DOUBLE kMax=3.5; // DIFFERENTIAL kMax [GeV]
DOUBLE dkBin=0.5; // MOMENTUM INTEGRATION BIN // ONLY FOR DIFFERENTIAL // [GeV]
DOUBLE dkStep=0.5; // STEP IN MOMENTUM [GeV]
DOUBLE kRefMin=1.0; // MINIMUM REF BIN [GeV]
DOUBLE kRefMax=2.0; // MINIMUM REF BIN [GeV]

INT kDIF_FLAG=2; // DIFFERENTIAL HARMONIC : 0 -- NO, 1 -- e^in\phi dk k , 2 -- e^in\phi k //
INT kINT_FLAG=1; // INTEGRATED HARMONIC : 0 -- N0, 1 -- YES //

// OPTION FOR OUTPUT OF CONFIGURATIONS -- CAN BE SET IN COMMANDLINE //
INT OUTPUT_FLAG=0;

// OPTION FOR CALCULATION OF ASYMMETRIC PART -- NEEDED FOR ODD HARMONICS -- CAN BE SET IN COMMANDLINE //
INT AS_FLAG=0;

// CONSTANT COLOR CHARGE DENSITY AND SIZE -- CAN SET IN COMMANDLINE//
#if IC_FLAG==MV_FLAG
DOUBLE g2muP=0.25; DOUBLE g2muT=1.0;
INT PROJ_FLAG=0; // 0 -- p, 1 -- d, 2 -- infinite //
#endif

// DEFINE RANDOM NUMBER GENERATOR //
#include "MISC/RNG/GSLRandomNumberGenerator.cpp"

///////////////////////////
//  OUTPUT HANDLING      //
///////////////////////////

#include "IO/StringManipulation.cpp"
#include "IO/OutputManagement.cpp"
#include "IO/InputFileManagement.cpp"
#include "IO/SaveConfiguration.cpp"

// GLOBAL RANDOM GENERATOR SEED //
INT GLOBAL_RNG_SEED;
// MPI RANDOM NUMBER SEED //
INT MY_MPI_RNG_SEED;

// INCLUDE INITIALIZATION //
#if IC_FLAG==MV_FLAG
#include "INITIALCONDITIONS/SetMV.cpp"
#endif

#if IC_FLAG==GIPS_FLAG
// INCLUDES IMPACT PARAMETER SAMPLING //
#include "INITIALCONDITIONS/SetGlauberIPSat.cpp"
#endif

// INCLUDE COMMANDLINE PROCESSING //
INT NumberOfConfigurations=1; // CAN BE SET IN COMMANDLINE //
#include "IO/cfile.c"
#include "CommandLineArguments.cpp"

///////////////////////////
//     OBSERVABLES       //
///////////////////////////

// INCLUDE OBSERVABLES
#include "OBSERVABLES/SingleInclusiveDistribution.cpp"
#include "OBSERVABLES/MultiplicityMeasurement.cpp"
#include "OBSERVABLES/DifferentialDistribution.cpp"
#include "OBSERVABLES/IntegratedDistribution.cpp"
#include "OBSERVABLES/FullDistribution.cpp"

// OUTPUT DIRECTORY -- CAN BE SET IN COMMANDLINE //
std::string OutputDirectory;

// CREATE INFO FILE //
#include "CreateInfoFile.cpp"

void Run(int argc,char **argv,INT MPI_RNG_SEED){

    ///////////
    // SETUP //
    ///////////

    MY_MPI_RNG_SEED=MPI_RNG_SEED;

    std::cerr << "# SEED=" << MY_MPI_RNG_SEED << std::endl;
    
    //INITIALIZE RANDOM NUMBER GENERATOR //
    RandomNumberGenerator::Init(MY_MPI_RNG_SEED);

    // INITIALIZE SET GROUP ELEMENTS TO UNITY AND ALGEBRA ELEMENTS TO ZERO //
    Lattice::InitTargetAndProjFields(); ProjSolution::A->SetZero();
    TargetFields::U->SetIdentity(); TempTargetFields::U->SetIdentity();
    #if IC_FLAG==GIPS_FLAG
    Lattice::InitTpp();     Tpp::Targ->SetZero();       Tpp::Proj->SetZero();
    #endif
    Lattice::Initg2mu();    g2mu::Targ->SetZero();      g2mu::Proj->SetZero();

    // INITIALIZE RHO FOURIER SPACE //
    FourierSpace::InitRhos();

    // SET MV INITIAL FIELDS //
    #if IC_FLAG==MV_FLAG
    FourierSpace::InitUtarget();        FourierSpace::InitScalar();
    InitialConditions::SetMVTarget();   InitialConditions::SetMVProjectile();
    FourierSpace::CleanUpUtarget();     FourierSpace::CleanUpScalar();
    #endif
    
    // SET GLAUBER+IP-SAT INITIAL FIELDS //
    #if IC_FLAG==GIPS_FLAG
    InitialConditions::SetGlauberIPSat(MPI_RNG_SEED);
    #endif
    
    // CLEANUP 2mu FEILDS, AND RHO FIELDS //
    Lattice::CleanUpg2mu(); FourierSpace::CleanUpRhos();

    // INITIALIZE OMEGA LATTICE FIELDS, OMEGA AND SCALAR FOURIER SPACE //
    Lattice::InitOmegas();  OmegaS::O->SetZero();   OmegaA::O->SetZero();
    FourierSpace::InitOmegas(); FourierSpace::InitScalar();
    
    // DETERMINE OMEGA_S AND OMEGA_A //
    Observables::ConstructOmegas(OmegaS::O,OmegaA::O);
    Observables::FourierTransformOmegas(OmegaS::O,OmegaA::O);

    // MEASURE MULTIPLICITY //
    // DO NOT REMOVE -- CRITICAL FOR CALCULATING ALL HARMONICS //
    Observables::MeasureMultiplicity();

    // OPTION FOR SINGLE INCLUSIVE DISTRIBUTION OUTPUT
    Observables::DetermineFullDistribution();

    // COMPUTE DIFFERENTIAL AND/OR INTEGRATED HARMONICS FOR A GIVEN SET OF RHOS //
    if(kDIF_FLAG==1 || kDIF_FLAG==2){
        // DIFFERENT SCHEMES : 1 -- \phi dk k , 2 -- \phi k //
        // STANDARD CHOICE IS 1 //
        Observables::DetermineDifferentialDistribution();
    }
    else{
        std::cerr << "# NOT COMPUTING DIFFERENTIAL HARMONICS " << std::endl;
    }
    if(kINT_FLAG==1){
        Observables::DetermineIntegratedDistribution();
    }
    else{
        std::cerr << "# NOT COMPUTING INTEGRATED HARMONICS " << std::endl;
    }
    
    // CLEAN-UP TARGET AND PROJECTILE FIELDS, OMEGAS, SCALAR //
    Lattice::CleanUpTargetAndProjFields(); Lattice::CleanUpOmegas();
    FourierSpace::CleanUpOmegas(); FourierSpace::CleanUpScalar();

    //Lattice::CleanUp();
    RandomNumberGenerator::CleanUp();

}

int main(int argc,char **argv){

    ///////////////////////////////
    //INITIALIZE MPI ENVIRONMENT //
    ///////////////////////////////

    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);

    //SET MPI ID's AND NUMBER OF NODES

    MPI_Comm_rank(MPI_COMM_WORLD,&MPIBasic::ID);
    MPI_Comm_size(MPI_COMM_WORLD,&MPIBasic::NumberOfNodes);

    //////////////////////////////////
    //PROCESS COMMANDLINE ARGUMENTS //
    //////////////////////////////////
    ProcessCommandlineArguments(argc,argv);

    ////////////////////////////////////
    // INITIALIZE OPEN MP ENVIRONMENT //
    ////////////////////////////////////

    //INITIALIZE THREADED FFTW
    int FFTW3_THREAD_STATUS=fftw_init_threads();
    
    //CHECK THREAD COUNT
    std::cerr << "# NUMBER OF THREADS= " << omp_get_max_threads() << " FFTW THREAD STATUS= " << FFTW3_THREAD_STATUS << std::endl;

    if(FFTW3_THREAD_STATUS==1){

        fftw_plan_with_nthreads(omp_get_max_threads());

    }

    //////////////////////////////////
    //PROCESS COMMANDLINE ARGUMENTS //
    //////////////////////////////////

    Konfig arguments(argc,argv);

    // NUMBER OF CONFIGURATIONS RUN IN SERIAL //
    arguments.Getval("nconfs",NumberOfConfigurations);

    std::cerr << "# RUNNING " << NumberOfConfigurations << " CONFIGURATIONS " << std::endl;

    // SAMPLE DIFFERENT CONFIGURATIONS //
    for(INT n=0;n<NumberOfConfigurations;n++){

        //SET GLOBAL RANDOM NUMBER SEED//
        INT GLOBAL_RNG_SEED;

        if(MPIBasic::ID==0){

            GLOBAL_RNG_SEED=time(0);

            // SPECIFY SEED IN COMMANDLINE //
            arguments.Getval("SEED",GLOBAL_RNG_SEED);
            
            // CREATE INFO FILE //
            CreateInfoFile(GLOBAL_RNG_SEED);
            std::cerr << "# INFO FILE CREATED " << std::endl;
            
        }

        // BROADCAST GLOBAL RANDOM SEED //
        MPI_Bcast(&GLOBAL_RNG_SEED, 1, MPI_INT,0,MPI_COMM_WORLD);

        // PERFORM CLASSICAL STATISTICAL SIMULATION //
        Run(argc,argv,GLOBAL_RNG_SEED+MPIBasic::ID);

        // COMMADNLINE NOTIFICATION //
        std::cerr << "# COMPLETED " << GLOBAL_RNG_SEED+MPIBasic::ID << std::endl;

    }

    //SYNCHRONIZE ALL MPI NODES
    MPI_Barrier(MPI_COMM_WORLD);

    //FINALIZE MPI
    MPI_Finalize();

    //////////
    // EXIT //
    //////////

    exit(0);

}
