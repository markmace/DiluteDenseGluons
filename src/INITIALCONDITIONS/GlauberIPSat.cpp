// GlauberIPSat.cpp is adapted from the IP-Glasma solver.
// Copyright (C) 2018 Bjoern Schenke.
// Reproduced with permission by Mark Mace 2019 for Dilute-Dense gluon solver

#include <stdio.h>
#include <string>
#include <cmath>
#include <ctime>
#include <iostream>
#include <complex>
#include <fstream>
#include <vector>

#include "Setup.h"
#include "Init.h"
#include "Random.h"
#include "FFT.h"
#include "Parameters.h"
#include "Matrix.h"
#include "GIPSLattice.h"
#include "Spinor.h"

#define _SECURE_SCL 0
#define _HAS_ITERATOR_DEBUGGING 0
using namespace std;

int readInput(string InputFile,Setup *setup, Parameters *param, string OutDirectory,long int RNG_SEED);

// GLAUBER + IP-Sat INITIALIZATION //
double Initialize(string InputFile,string OutDirectory,long int RNG_SEED, double *PROJ_g2mu2,double *TARG_g2mu2,int OUTPUT_FLAG){
    
    Setup *setup;
    setup = new Setup();
    Random *random;
    random = new Random();
    Parameters *param;
    param = new Parameters();
    
    // READ INPUT PARAMETERS //
    readInput(InputFile,setup,param,OutDirectory,RNG_SEED);
    
    int cells = param->getSize()*param->getSize();
    int Nc2m1 = param->getNc()*param->getNc()-1; // N_c^2-1
    int nn[2];
    int pos;
    double x,y;
    nn[0]=param->getSize();
    nn[1]=param->getSize();
    
    stringstream strup_name;
    strup_name << OutDirectory << "/usedParameters" << RNG_SEED << ".dat";
    string up_name;
    up_name = strup_name.str();
    
    // INITIALIZE INIT CLASS //
    Init *init;
    init = new Init(nn);
    
    // INITIALIZE GROUP CLASS //
    Group *group;
    group = new Group(param->getNc());
    
    // INITIALIZE GLAUBER CLASS //
    Glauber *glauber;
    glauber = new Glauber;

    // INITIALIZE GLAUBER //
    glauber->initGlauber(param->getSigmaNN(),param->getTarget(),param->getProjectile(),param->getb(),100,RNG_SEED);
    
    // RESET SUCCESS INDICATOR //
    param->setSuccess(0);
    
    // INITIALIZE FIRST TRY INDEX //
    int FirstTry=0;
    
    while(param->getSuccess()==0){
        
        // RESET SUCCESS INDICATOR //
        param->setSuccess(0);
        // ALLOCATE LATTICE //
        GIPSLattice *MyLat;
        MyLat=new GIPSLattice(param,param->getNc(),param->getSize());
        
        // SET rnum SEED //
        unsigned long long int rnum;
        if(FirstTry==0){
            rnum=RNG_SEED;
        }
        else{
            rnum=time(0)+RNG_SEED;
        }
        
        // INITIALIZE RNG FOR INITIAL CONDITION METHODS //
        random->init_genrand64(rnum);
        random->gslRandomInit(rnum);
        
        int READFROMFILE = 0;
        init->init(MyLat,group,param,random,glauber,READFROMFILE,OutDirectory,RNG_SEED,OUTPUT_FLAG);
        
        if(param->getSuccess()==0){
            
            std::cerr << "# RESAMPLING COLOR CONFIGURATIONS" << std::endl;
            // DENOTE TO TRY FOR NEW rnum //
            FirstTry++;
            delete MyLat;
            continue;
        }
        
        // CLEAN-UP //
        delete init;
        delete random;
        delete glauber;
        
        // SAVE g2mu //
        std::cerr << "# SETTING FINAL g2mu2 FIELDS " << std::endl;
        for(pos=0;pos<param->getSize()*param->getSize();pos++){
            PROJ_g2mu2[pos]=MyLat->cells[pos]->getg2mu2A(); // DIMLESS //
            TARG_g2mu2[pos]=MyLat->cells[pos]->getg2mu2B(); // DIMLESS //
        }
        
        // SAVE FINAL rNUM TO FILE //
        ofstream IPSatInfo(up_name.c_str(),ios::app);
        IPSatInfo << "Random seed used: " << param->getRandomSeed() << endl;
        IPSatInfo << "rnum used: " << rnum << endl;
        IPSatInfo.close();
        
        // CLEAN-UP //
        delete MyLat;
    }
    // CLEAN-UP //
    delete group;
    delete param;
    delete setup;

    // RETURN USED IMPACT PARAMETER //
    return param->getb();
    
}


int readInput(string file_name,Setup *setup, Parameters *param, string OutDirectory,long int RNG_SEED){
    
    //std::cerr << "# READING IN PARAMETERS FROM FILE ";
    param->setNucleusQsTableFileName(setup->StringFind(file_name,"NucleusQsTableFileName"));
    param->setNucleonPositionsFromFile(setup->IFind(file_name,"nucleonPositionsFromFile"));
    param->setTarget(setup->StringFind(file_name,"Target"));
    param->setProjectile(setup->StringFind(file_name,"Projectile"));
    param->setRunningCoupling(setup->IFind(file_name,"runningCoupling"));
    param->setL(setup->DFind(file_name,"L"));
    param->setBG(setup->DFind(file_name,"BG"));
    param->setBGq(setup->DFind(file_name,"BGq"));
    param->setRnp(setup->DFind(file_name,"rnp"));
    param->setMuZero(setup->DFind(file_name,"muZero"));
    param->setc(setup->DFind(file_name,"c"));
    param->setSize(setup->IFind(file_name,"size"));
    param->setUseFluctuatingx(setup->IFind(file_name,"useFluctuatingx"));
    param->setNc(setup->IFind(file_name,"Nc"));
    param->setSeed(setup->ULLIFind(file_name,"seed"));
    param->setUseSeedList(setup->IFind(file_name,"useSeedList"));
    param->setRoots(setup->DFind(file_name,"roots"));
    param->setg(setup->DFind(file_name,"g"));
    param->setJacobianm(setup->DFind(file_name,"Jacobianm"));
    param->setSigmaNN(setup->DFind(file_name,"SigmaNN"));
    param->setRmax(setup->DFind(file_name,"rmax"));
    param->setbmin(setup->DFind(file_name,"bmin"));
    param->setbmax(setup->DFind(file_name,"bmax"));
    param->setQsmuRatio(setup->DFind(file_name,"QsmuRatio"));
    param->setUsePseudoRapidity(setup->DFind(file_name,"usePseudoRapidity"));
    param->setRapidity(setup->DFind(file_name,"Rapidity"));
    param->setUseNucleus(setup->IFind(file_name,"useNucleus"));
    param->setg2mu(setup->DFind(file_name,"g2mu"));
    param->setRunWithQs(setup->IFind(file_name,"runWith0Min1Avg2MaxQs"));
    param->setRunWithkt(setup->IFind(file_name,"runWithkt"));
    param->setRunWithLocalQs(setup->IFind(file_name,"runWithLocalQs"));
    param->setRunWithThisFactorTimesQs(setup->DFind(file_name,"runWithThisFactorTimesQs"));
    param->setxFromThisFactorTimesQs(setup->DFind(file_name,"xFromThisFactorTimesQs"));
    param->setLinearb(setup->IFind(file_name,"samplebFromLinearDistribution"));
    param->setAverageOverNuclei(setup->IFind(file_name,"averageOverThisManyNuclei"));
    param->setUseTimeForSeed(setup->IFind(file_name,"useTimeForSeed"));
    param->setUseFixedNpart(setup->IFind(file_name,"useFixedNpart"));
    param->setSmearQs(setup->IFind(file_name,"smearQs"));
    param->setSmearingWidth(setup->DFind(file_name,"smearingWidth"));
    param->setGaussianWounding(setup->IFind(file_name,"gaussianWounding"));
    param->setProtonAnisotropy(setup->DFind(file_name,"protonAnisotropy"));
    param->setUseConstituentQuarkProton(setup->DFind(file_name,"useConstituentQuarkProton"));
    std::cerr << "# PARAMETERS READ IN FROM FILE " << file_name << endl;
    
    // write the used parameters into file "usedParameters.dat" as a double check for later
    time_t rawtime = time(0);
    stringstream strup_name;
    strup_name << OutDirectory << "/usedParameters" << RNG_SEED << ".dat";
    string up_name;
    up_name = strup_name.str();
    
    fstream IPSatInfo(up_name.c_str(),ios::out);
    char * timestring = ctime(&rawtime);
    IPSatInfo << "File created on " << timestring << endl;
    IPSatInfo << "##################################################### " << endl;
    IPSatInfo << "Used parameters by Glauber+IP-Sat from IP-Glasma v1.3" << endl;
    IPSatInfo << "##################################################### " << endl;
    IPSatInfo << " " << endl;
    IPSatInfo << "Output by readInput in bGlauberIPSat.cpp: " << endl;
    IPSatInfo << " " << endl;
    IPSatInfo << "Nc " << param->getNc() << endl;
    IPSatInfo << "size " << param->getSize() << endl;
    IPSatInfo << "lattice spacing a " << param->getL()/static_cast<double>(param->getSize()) << " fm " << endl;
    IPSatInfo << "Projectile " << param->getProjectile() << endl;
    IPSatInfo << "Target " << param->getTarget() << endl;
    IPSatInfo << "Gaussian wounding " << param->getGaussianWounding() << endl;
    IPSatInfo << "Using fluctuating x=Qs/root(s) " << param->getUseFluctuatingx() << endl;
    if( param->getRunWithkt()==0)
    IPSatInfo << "Using local Qs to run " << param->getRunWithLocalQs() << endl;
    else
    IPSatInfo << "running alpha_s with k_T" << endl;
    IPSatInfo << "QsmuRatio " << param->getQsmuRatio() << endl;
    IPSatInfo << "smeared mu " << param->getSmearQs() << endl;
    IPSatInfo << "rmax " << param->getRmax() << endl;
    if (param->getSmearQs()==1){
        
        IPSatInfo << "smearing width " << param->getSmearingWidth() << endl;
    }
    IPSatInfo.close();
    
    return 0;
}
