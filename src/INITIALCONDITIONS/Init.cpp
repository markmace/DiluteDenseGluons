// Init.cpp is adapted from the IP-Glasma solver.
// Copyright (C) 2018 Bjoern Schenke.
// Reproduced with permission by Mark Mace 2019 for Dilute-Dense gluon solver

#include "Init.h"

//**************************************************************************
// Init class.

void Init::sampleTA(Parameters *param, Random* random, Glauber* glauber){
    
    ReturnValue rv, rv2;
    std::cerr << "# GLAUBER SAMPLING NUCLEON POSITIONS " << std::endl;
    
    if(param->getNucleonPositionsFromFile()==0){
        
        int A1,A2;
        A1 = static_cast<int>(glauber->nucleusA1())*param->getAverageOverNuclei(); // projectile
        A2 = static_cast<int>(glauber->nucleusA2())*param->getAverageOverNuclei(); // target
        
        std::cerr << "# A1= " << A1 <<  " A2= " << A2 << std::endl;
        if(param->getRnp()>=0.0){

            std::cerr << "# FIXED Rnp= " << param->getRnp() << std::endl;
        }
        
        if((glauber->nucleusA1()==1 || glauber->nucleusA2()==1) && param->getAverageOverNuclei()>1)
        {
            std::cerr << "Averaging not supported for collisions involving protons ... Exiting." << endl;
            exit(1);
        }
        
        if(A1==1)
        {
            rv.x=0.;
            rv.y=0;
            rv.collided=0;
            nucleusA.push_back(rv);
        }
        else if(A1==2) // deuteron
        {
            // IF ASSIGNED IN COMMANDLINE
            if(param->getRnp()>=0.0){

                // SET FIRST NUCLEON AT dnp/2.0 ALONG x-AXIS //
                rv.x=param->getRnp()/2.0;
                rv.y=0.0;
                nucleusA.push_back(rv);
                
                // SET SECOND NUCLEON AT -dnp/2.0 ALONG x-AXIS //
                rv.x=-param->getRnp()/2.0;
                rv.y=0.0;
                rv.collided=0;
                nucleusA.push_back(rv);
            }
            else{
                // SAMPLE DEUTERON DISTANCE VECTOR //
                rv=glauber->SampleTARejection(random,1);
                // SET DISTANCE //
                param->setRnp(sqrt(rv.x*rv.x+rv.y*rv.y));
                // SHIFT TO CoM FOR NEUTRON-PROTON SYSTEM //
                rv.x=rv.x/2.;
                rv.y=rv.y/2.;
                nucleusA.push_back(rv);
                
                // ROTATE BY 180 DEGREES //
                rv.x=-rv.x;
                rv.y=-rv.y;
                rv.collided=0;
                nucleusA.push_back(rv);

            }
            
        }
        else if(A1==3) // He3
        {
            //sample the position in the file
            ifstream fin;
            fin.open("src/MISC/he3.dat");
            
            double dummy;
            double ran2 = random->genrand64_real3();   // sample the position in the file uniformly (13699 events in file)
            int nucleusNumber = static_cast<int>(ran2*13699);
            
            //std::cerr << "# USING PROJ HELIUM-3 CONFIGUARTION = " << nucleusNumber << endl;
            
            // go to the correct line in the file
            fin.seekg(std::ios::beg);
            for(int i=0; i < nucleusNumber; ++i)
            {
                fin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            }
            // am now at the correct line in the file
            
            // start reading one nucleus (3 positions)
            int A=0;
            
            while(A<glauber->nucleusA1())
            {
                if(!fin.eof())
                {
                    fin >> rv.x;
                    fin >> rv.y;
                    fin >> dummy; // don't care about z direction
                    rv.collided=0;
                    nucleusA.push_back(rv);
                    A++;
                    //std::cerr << "### A=" << A << ", x=" << rv.x << ", y=" << rv.y << endl;
                }
            }
            
            
            fin.close();
            
            param->setA1FromFile(A);
            
        }
        else
        {
            for (int i = 0; i < A1; i++) // get all nucleon coordinates
            {
                rv = glauber->SampleTARejection(random,1);
                nucleusA.push_back(rv);
            }
        }
        
        if(A2==1)
        {
            rv2.x=0.;
            rv2.y=0;
            rv2.collided=0;
            nucleusB.push_back(rv2);
        }
        else if(A2==2) // deuteron
        {
            rv = glauber->SampleTARejection(random,2);
            // we sample the neutron proton distance, so distance to the center needs to be divided by 2
            param->setRnp(sqrt(rv.x*rv.x+rv.y*rv.y));
            
            rv.x = rv.x/2.;
            rv.y = rv.y/2.;
            nucleusB.push_back(rv);
            
            // other nucleon is 180 degrees rotated:
            rv.x = -rv.x;
            rv.y = -rv.y;
            rv.collided=0;
            nucleusB.push_back(rv);
            
        }
        else if(A2==3) // He3
        {
            //sample the position in the file
            ifstream fin;
            fin.open("src/MISC/he3.dat");
            
            double dummy;
            double ran2 = random->genrand64_real3();   // sample the position in the file uniformly (13699 events in file)
            int nucleusNumber = static_cast<int>(ran2*13699);
            
            //std::cerr << "# USING TARG HELIUM-3 CONFIGUARTION = " << nucleusNumber << endl;
            
            // go to the correct line in the file
            fin.seekg(std::ios::beg);
            for(int i=0; i < nucleusNumber; ++i)
            {
                fin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            }
            // am now at the correct line in the file
            
            // start reading one nucleus (3 positions)
            int A=0;
            
            while(A<glauber->nucleusA2())
            {
                if(!fin.eof())
                {
                    fin >> rv.x;
                    fin >> rv.y;
                    fin >> dummy; // don't care about z direction
                    rv.collided=0;
                    nucleusB.push_back(rv);
                    A++;
                    //std::cerr << "# A=" << A << ", x=" << rv.x << ", y=" << rv.y << endl;
                }
            }
            
            
            fin.close();
            
            param->setA2FromFile(A);
            
        }
        else
        {
            for (int i = 0; i < A2; i++) // get all nucleon coordinates
            {
                rv2 = glauber->SampleTARejection(random,2);
                //cerr << "N=" << i << ", x=" << rv2.x << ", y=" << rv2.y << endl;
                nucleusB.push_back(rv2);
            }
        }
        std::cerr << "# FINISHED GLAUBER SETTING" << endl;
    }
    else if (param->getNucleonPositionsFromFile()==1)
    {
        std::cerr << "NucleonPositionsFromFile can be 0 (sample nucleons) - no file options yet. you chose " <<
        param->getNucleonPositionsFromFile() << ". Exiting." << endl;
        exit(0);
    }
    else
    {
        std::cerr << "NucleonPositionsFromFile can be 0 (sample nucleons) or 1 or 2 (read from files) - you chose " <<
        param->getNucleonPositionsFromFile() << ". Exiting." << endl;
        exit(0);
    }
}

// SET FOR ONLY PROJECITLE //
void Init::sampleProjectileTA(Parameters *param, Random* random, Glauber* glauber){
    
    ReturnValue rv, rv2;
    std::cerr << "# GLAUBER SAMPLING PROJECTILE NUCLEON POSITIONS " << std::endl;
    
    if(param->getNucleonPositionsFromFile()==0){
        
        int A1;
        A1 = static_cast<int>(glauber->nucleusA1())*param->getAverageOverNuclei(); // projectile
        
        std::cerr << "# A1= " << A1 << std::endl;
        if(param->getRnp()>=0.0){
            
            std::cerr << "# FIXED Rnp= " << param->getRnp() << std::endl;
        }
        
        if((glauber->nucleusA1()==1 || glauber->nucleusA2()==1) && param->getAverageOverNuclei()>1)
        {
            std::cerr << "Averaging not supported for collisions involving protons ... Exiting." << endl;
            exit(1);
        }
        
        if(A1==1)
        {
            rv.x=0.;
            rv.y=0;
            rv.collided=0;
            nucleusA.push_back(rv);
        }
        else if(A1==2) // deuteron
        {
            // IF ASSIGNED IN COMMANDLINE
            if(param->getRnp()>=0.0){
                
                // SET FIRST NUCLEON AT dnp/2.0 ALONG x-AXIS //
                rv.x=param->getRnp()/2.0;
                rv.y=0.0;
                nucleusA.push_back(rv);
                
                // SET SECOND NUCLEON AT -dnp/2.0 ALONG x-AXIS //
                rv.x=-param->getRnp()/2.0;
                rv.y=0.0;
                rv.collided=0;
                nucleusA.push_back(rv);
            }
            else{
                // SAMPLE DEUTERON DISTANCE VECTOR //
                rv=glauber->SampleTARejection(random,1);
                // SET DISTANCE //
                param->setRnp(sqrt(rv.x*rv.x+rv.y*rv.y));
                // SHIFT TO CoM FOR NEUTRON-PROTON SYSTEM //
                rv.x=rv.x/2.;
                rv.y=rv.y/2.;
                nucleusA.push_back(rv);
                
                // ROTATE BY 180 DEGREES //
                rv.x=-rv.x;
                rv.y=-rv.y;
                rv.collided=0;
                nucleusA.push_back(rv);
                
            }
            
        }
        else if(A1==3) // He3
        {
            //sample the position in the file
            ifstream fin;
            fin.open("src/MISC/he3.dat");
            
            double dummy;
            double ran2 = random->genrand64_real3();   // sample the position in the file uniformly (13699 events in file)
            int nucleusNumber = static_cast<int>(ran2*13699);
            
            //std::cerr << "# USING PROJ HELIUM-3 CONFIGUARTION = " << nucleusNumber << endl;
            
            // go to the correct line in the file
            fin.seekg(std::ios::beg);
            for(int i=0; i < nucleusNumber; ++i)
            {
                fin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            }
            // am now at the correct line in the file
            
            // start reading one nucleus (3 positions)
            int A=0;
            
            while(A<glauber->nucleusA1())
            {
                if(!fin.eof())
                {
                    fin >> rv.x;
                    fin >> rv.y;
                    fin >> dummy; // don't care about z direction
                    rv.collided=0;
                    nucleusA.push_back(rv);
                    A++;
                    //std::cerr << "### A=" << A << ", x=" << rv.x << ", y=" << rv.y << endl;
                }
            }
            
            
            fin.close();
            
            param->setA1FromFile(A);
            
        }
        else
        {
            for (int i = 0; i < A1; i++) // get all nucleon coordinates
            {
                rv = glauber->SampleTARejection(random,1);
                nucleusA.push_back(rv);
            }
        }
        
        
        std::cerr << "# FINISHED GLAUBER SETTING PROJECITLE NUCLEUS" << endl;
    }
    else if (param->getNucleonPositionsFromFile()==1)
    {
        std::cerr << "NucleonPositionsFromFile can be 0 (sample nucleons) - no file options yet. you chose " <<
        param->getNucleonPositionsFromFile() << ". Exiting." << endl;
        exit(0);
    }
    else
    {
        std::cerr << "NucleonPositionsFromFile can be 0 (sample nucleons) or 1 or 2 (read from files) - you chose " <<
        param->getNucleonPositionsFromFile() << ". Exiting." << endl;
        exit(0);
    }
}


void Init::readNuclearQs(Parameters *param){
    
    //std::cerr << "# READING Q_s(sum(T_p),y) FROM FILE ";
    // steps in qs0 and Y in the file
    // double y[iymaxNuc];
    // double qs0[ibmax];
    string dummy;
    string T, Qs;
    // open file
    ifstream fin;
    fin.open((param->getNucleusQsTableFileName()).c_str());
    
    //std::cerr << param->getNucleusQsTableFileName() << " ... " ;
    
    if(fin)
    {
        for (int iT=0; iT<iTpmax; iT++)
        {
            for (int iy=0; iy<iymaxNuc; iy++)
            {
                if (!fin.eof())
                {
                    fin >> dummy;
                    fin >> T;
                    Tlist[iT]=atof(T.c_str());
                    fin >> Qs;
                    Qs2Nuclear[iT][iy]=atof(Qs.c_str());
                    //cerr << iT << " " << iy <<  " " << T << " " << iy*deltaYNuc << " " << Qs2Nuclear[iT][iy] << endl;
                }
                else
                {
                    std::cerr << "# CRITICAL ERROR  -- END OF QS FILE REACHED PREMATURELY" << std::endl;
                    exit(1);
                }
            }
        }
        fin.close();
        //std::cerr << " done." << endl;
    }
    else
    {
        std::cerr << "[Init.cpp:readNuclearQs]: File qs2_Adj_Y_qs20_IPSat.dat does not exist. Exiting." << endl;
        exit(1);
    }
}

// Q_s as a function of \sum T_p and y (new in this version of the code - v1.2 and up)
double Init::getNuclearQs2(Parameters *param, Random* random, double T, double y){
    
    double value, fracy, fracT, QsYdown, QsYup;
    int posb, posy, check=0;
    posy = static_cast<int>(floor(y/deltaYNuc+0.0000001));
    
    if (y>iymaxNuc*deltaYNuc){
        
        std::cerr << " [Init:getNuclearQs2]:ERROR: y out of range. Maximum y value is " << iymaxNuc*deltaYNuc << ", you used " << y << ". Exiting." << endl;
        exit(1);
    }
    
    if ( T > Qs2Nuclear[iTpmax-1][iymaxNuc-1] ){
        
        std::cerr << "T=" << T << ", maximal T in table=" << Tlist[iTpmax-1] << std::endl;
        std::cerr << " [Init:getNuclearQs2]:ERROR: out of range. Exiting." << std::endl;
        exit(1);
    }
    
    if ( T < Tlist[0] ){
        
        check = 1;
        return 0.;
    }
    
    for(int iT=0; iT<iTpmax; iT++){
        
        if(T>=Tlist[iT] && T<Tlist[iT+1]){
            
            fracT = (T-Tlist[iT])/(Tlist[iT+1]-Tlist[iT]);
            fracy = (y-static_cast<double>(posy)*deltaYNuc)/deltaYNuc;
            
            QsYdown = (fracT)*(Qs2Nuclear[iT+1][posy])+(1.-fracT)*(Qs2Nuclear[iT][posy]);
            QsYup = (fracT)*(Qs2Nuclear[iT+1][posy+1])+(1.-fracT)*(Qs2Nuclear[iT][posy+1]);
            
            //	  cerr << posy << endl;
            //cerr << Qs2Nuclear[iT+1][posy] << " " <<  QsYdown << " " << QsYup << endl;
            
            value = (fracy*QsYup+(1.-fracy)*QsYdown);//*hbarc*hbarc;
            
            //   	  cerr << "T=" << T << ", lowT=" << Tlist[iT] << ", highT=" << Tlist[iT+1] << endl;
            //   	  cerr << "y=" << y << ", lowy=" << (posy)*deltaYNuc << ", highy=" << (posy+1)*deltaYNuc << endl;
            // 	  cerr << "fracy=" << fracy << endl;
            // 	  cerr << "Qs^2=" << value << endl;
            //   	  cerr << "Qs2Nuclear[iT][posy]=" << Qs2Nuclear[iT][posy] << endl;
            // 	  cerr << "Qs2Nuclear[iT][posy+1]=" << Qs2Nuclear[iT][posy+1] << endl;
            // 	  cerr << "Qs2Nuclear[iT+1][posy]=" << Qs2Nuclear[iT+1][posy] << endl;
            // 	  cerr << "Qs2Nuclear[iT+1][posy+1]=" << Qs2Nuclear[iT+1][posy+1] << endl;
            
            check++;
            continue;
        }
    }
    
    if (check!=1){
        
        cerr << check << ": T=" << T << endl ;
        // cerr << " [Init:getNuclearQs2]:ERROR: something went wrong in determining the value of Qs^2. Exiting." << endl;
        // exit(1);
        std::cerr << " [Init:getNuclearQs2]:ERROR: something went wrong in determining the value of Qs^2. Using maximal T_p" << std::endl;
        value = Tlist[iTpmax-1];
    }
    
    return value;
}


// set g^2\mu^2 as the sum of the individual nucleons' g^2\mu^2, using Q_s(b,y) prop tp g^mu(b,y)
// also compute N_part using Glauber
void Init::setColorChargeDensity(GIPSLattice *lat,Parameters *param,Random *random,Glauber *glauber,string OutDirectory,int RNGSEED,int OUTPUT_FLAG){
    
    int pos,posA,posB;
    int N = param->getSize();
    int A1, A2;
    int check=0;
    double xVal;
    if(param->getNucleonPositionsFromFile()==0)
    {
        A1 = static_cast<int>(glauber->nucleusA1())*param->getAverageOverNuclei();
        A2 = static_cast<int>(glauber->nucleusA2())*param->getAverageOverNuclei();
    }
    else
    {
        A1 = param->getA1FromFile();
        A2 = param->getA2FromFile();
    }
    
    int Npart = 0;
    int Ncoll = 0;
    double g2mu2A, g2mu2B;
    double b = param->getb();
    double x, xm;
    double y, ym;
    double r;
    double L = param->getL();
    double rapidity;
    double P,m;
    if(param->getUsePseudoRapidity()==0){
        rapidity = param->getRapidity();
    }
    else{
        // when using pseudorapidity as input convert to rapidity here. later include Jacobian in multiplicity and energy
        //std::cerr << "# USING PSEUDORAPIDITY " << param->getRapidity() << endl;
        m=param->getJacobianm(); // in GeV
        P=0.13+0.32*pow(param->getRoots()/1000.,0.115); //in GeV
        rapidity = 0.5 * log(sqrt(pow(cosh(param->getRapidity()),2.)+m*m/(P*P))+sinh(param->getRapidity())
                             / (sqrt(pow(cosh(param->getRapidity()),2.)+m*m/(P*P))-sinh(param->getRapidity())));
        //std::cerr << "# CORRESPONDS TO RAPIDITY " << rapidity << endl;
    }
    
    double yIn = rapidity;//param->getRapidity();
    double a = L/N; // lattice spacing in fm
    double dx, dy, dij;
    double d2 = param->getSigmaNN()/(PI*10.);          // in fm^2
    double averageQs = 0.;
    double averageQs2 = 0.;
    double averageQs2Avg = 0.;
    double averageQs2min = 0.;
    double averageQs2min2 = 0.;
    int count = 0;
    int count2 = 0;
    double nucleiInAverage;
    nucleiInAverage = static_cast<double>(param->getAverageOverNuclei());
    
    
    
    //testing distributions
    //   double binsNBD[30];
    //   double binsGauss[30];
    //   double nbd;
    //   double scale=2.;
    //   int bins=30;
    //   int events=2000;
    
    
    double gaussA[A1][3];
    double gaussB[A2][3];
    
    for (int i = 0; i<A1; i++){
        
        for (int iq = 0; iq<3; iq++){
            
            gaussA[i][iq]=1.;
        }
    }
    for (int i = 0; i<A2; i++){
        
        for (int iq = 0; iq<3; iq++){
            
            gaussB[i][iq]=1.;
        }
    }
    
    // let the log fluctuate
    if(param->getSmearQs() == 1){
        if(A1>0){
            for(int i = 0; i<A1; i++){
                for(int iq = 0; iq<3; iq++){
                    
                    gaussA[i][iq] = (exp(random->Gauss(0,param->getSmearingWidth())))/1.13; // dividing by 1.13 restores the same mean Q_s
                    //cerr << i << " " << iq << " " << gaussA[i][iq] << endl;
                    //	  if (gaussA[i]<0)
                    //  gaussA[i]=0.;
                }
            }
        }
        if(A2>0){
            for (int i = 0; i<A2; i++){
                for (int iq = 0; iq<3; iq++){
                    
                    gaussB[i][iq] = (exp(random->Gauss(0,param->getSmearingWidth())))/1.13;
                    
                    //	      cerr << i << " " << iq << " " << gaussB[i][iq] << endl;
                    //if (gaussB[i]<0)
                    //  gaussB[i]=0.;
                }
            }
        }
        
        // CHECK p+p Qs ORDERING FOR DILUTE DENSE //
        if(A1==1 && A2==1){
            // CONSTITUENT QUARK -- USE SUM (NEEDS JUSTIFICATION!) //
            if(param->getUseConstituentQuarkProton()>0){
                
                // SUM Qs //
                double sumA=gaussA[0][0]+gaussA[0][1]+gaussA[0][2];
                double sumB=gaussB[0][0]+gaussB[0][1]+gaussB[0][2];
                
                // TEST //
                if(sumA>sumB){
                    
                    std:cerr << "# FIXING DILUTE-DENSE LIMIT " << std::endl;
                    // SWAP QC Qs //
                    for(int iq=0; iq<3; iq++){
                        double aTemp=gaussA[0][iq];
                        gaussA[0][iq]=gaussB[0][iq];
                        gaussB[0][iq]=aTemp;
                    }
                }
            }
            else{
                if(gaussA[0][0]>gaussA[0][0]){
                    
                    std::cerr << "# FIXING DILUTE-DENSE LIMIT " << std::endl;
                    // SWAP NUCLEON Qs //
                    double aTemp=gaussA[0][0];
                    gaussA[0][0]=gaussB[0][0];
                    gaussB[0][0]=aTemp;
                }
                
            }
        }
    }
    
    param->setQsmuRatioB(param->getQsmuRatio());
    
    // MV INITIAL CONDITIONS -- INFINITE NUCLEUS //
    if(param->getUseNucleus() == 0){
        for(int ix=0; ix<N; ix++){ // loop over all positions
            for(int iy=0; iy<N; iy++){
                
                pos = ix*N+iy;
                lat->cells[pos]->setg2mu2A(param->getg2mu()*param->getg2mu()/param->getg()/param->getg());
                lat->cells[pos]->setg2mu2B(param->getg2mu()*param->getg2mu()/param->getg()/param->getg());
            }
        }
        param->setSuccess(1);
        cerr << "# CONSTANT COLOR CHARGE DENSITY SET" << endl;
        return;
    }
    
    // MV TARGET INITIAL CONDITIONS -- INFINITE NUCLEUS //
    if(param->getUseNucleus() == 2){
        std::cerr << "## CATASTROPHIC FAILURE!!" << std::endl;
        exit(0);
    }
    
    
    for(int ix=0; ix<N; ix++){ // loop over all positions
        for(int iy=0; iy<N; iy++){
            
            pos = ix*N+iy;
            lat->cells[pos]->setg2mu2A(0.);
            lat->cells[pos]->setg2mu2B(0.);
        }
    }
    
    // compute N_part
    // positions are shifted here. not later as in previous versions. bshift below (in init(..)) is zero.
    
    if(A1 == 1 && A2 > 1){
        for(int i = 0; i<A2; i++){
            nucleusB.at(i).x=nucleusB.at(i).x+b;
        }
    }
    else if(A2 == 1 && A1 > 1){
        for(int i = 0; i<A1; i++){
            nucleusA.at(i).x=nucleusA.at(i).x-b;
        }
    }
    else{
        for(int i = 0; i<A1; i++){ // shift the nuclei's position by -b/2 or +b/2 respectively
            nucleusA.at(i).x=nucleusA.at(i).x-b/2.;
        }
        
        for (int i = 0; i<A2; i++){ // shift the nuclei's position by -b/2 or +b/2 respectively
            nucleusB.at(i).x=nucleusB.at(i).x+b/2.;
        }
    }
    
    
    double maxT=0;
    
    double bp2,T,BG;
    BG = param->getBG();
    double BGq = param->getBGq(); // quark size in GeV^-2
    double xi = param->getProtonAnisotropy();
    double phi;
    
    if(xi!=0.0){
        for (int i = 0; i<A1; i++){
            nucleusA.at(i).phi = 2*M_PI*random->genrand64_real2();
        }
        
        for (int i = 0; i<A2; i++){
            nucleusB.at(i).phi = 2*M_PI*random->genrand64_real2();
        }
    }
    
    //  cerr << "BG=" << BG << endl;
    
    
    double xq[A1][param->getUseConstituentQuarkProton()], xq2[A2][param->getUseConstituentQuarkProton()];
    double yq[A1][param->getUseConstituentQuarkProton()], yq2[A2][param->getUseConstituentQuarkProton()];
    double avgxq=0.;
    double avgyq=0.;
    
    if(param->getUseConstituentQuarkProton()>0){
        for (int i=0; i<A1; i++){
            avgxq=0.;
            avgyq=0.;
            for(int iq=0; iq<param->getUseConstituentQuarkProton(); iq++){
                xq[i][iq] = sqrt(BG*hbarc*hbarc)*random->Gauss();
                yq[i][iq] = sqrt(BG*hbarc*hbarc)*random->Gauss();
            }
            for(int iq=0; iq<3; iq++){
                avgxq += xq[i][iq];
                avgyq += yq[i][iq];
            }
            for(int iq=0; iq<3; iq++){
                xq[i][iq] -= avgxq/3.;
                yq[i][iq] -= avgyq/3.;
                //	      cerr << xq[i][iq] << " " << yq[i][iq] << endl;
            }
            
            
            // avgyq=0.;
            // for (int iq=0; iq<3; iq++)
            //   {
            //     avgxq += xq[i][iq];
            //     avgyq += yq[i][iq];
            //   }
            // cerr << avgyq << endl;
        }
    }
    
    if(param->getUseConstituentQuarkProton()>0){
        for(int i=0; i<A2; i++){
            avgxq=0.;
            avgyq=0.;
            for(int iq=0; iq<param->getUseConstituentQuarkProton(); iq++){
                xq2[i][iq] = sqrt(BG*hbarc*hbarc)*random->Gauss();
                yq2[i][iq] = sqrt(BG*hbarc*hbarc)*random->Gauss();
            }
            
            for(int iq=0; iq<3; iq++){
                avgxq += xq2[i][iq];
                avgyq += yq2[i][iq];
            }
            for(int iq=0; iq<3; iq++){
                xq2[i][iq] -= avgxq/3.;
                yq2[i][iq] -= avgyq/3.;
            }
            // avgyq=0.;
            // for (int iq=0; iq<3; iq++)
            //   {
            //     avgxq += xq2[i][iq];
            //     avgyq += yq2[i][iq];
            //   }
            // cerr << avgyq << endl;
        }
        
        
    }
    
    
    //add all T_p's (new in version 1.2)
    for(int ix=0; ix<N; ix++){ // loop over all positions
        
        x = -L/2.+a*ix;
        for(int iy=0; iy<N; iy++){
            y = -L/2.+a*iy;
            
            pos = ix*N+iy;
            
            // nucleus A
            lat->cells[pos]->setTpA(0.);
            for (int i = 0; i<A1; i++){
                xm = nucleusA.at(i).x;
                ym = nucleusA.at(i).y;
                
                if(param->getUseConstituentQuarkProton()>0){
                    T = 0.;
                    for(int iq=0; iq<param->getUseConstituentQuarkProton(); iq++){
                        bp2 = (xm+xq[i][iq]-x)*(xm+xq[i][iq]-x)+(ym+yq[i][iq]-y)*(ym+yq[i][iq]-y);
                        bp2 /= hbarc*hbarc;
                        
                        T += exp(-bp2/(2.*BGq))/(2.*PI*BGq)/(double(param->getUseConstituentQuarkProton()))*gaussA[i][iq]; // I removed the 2/3 here to make it a bit bigger
                        //	      cerr << "A " << i << " " << iq << " " << gaussA[i][iq] << endl;
                    }
                }
                else{
                    phi = nucleusA.at(i).phi;
                    
                    bp2 = (xm-x)*(xm-x)+(ym-y)*(ym-y) + xi*pow((xm-x)*cos(phi) + (ym-y)*sin(phi),2.);
                    bp2 /= hbarc*hbarc;
                    
                    T = sqrt(1+xi)*exp(-bp2/(2.*BG))/(2.*PI*BG)*gaussA[i][0]; // T_p in this cell for the current nucleon
                }
                
                lat->cells[pos]->setTpA(lat->cells[pos]->getTpA()+T/nucleiInAverage); // add up all T_p
                
                maxT=max(lat->cells[pos]->getTpA()+T,maxT);
                
            }
            
            // nucleus B
            lat->cells[pos]->setTpB(0.);
            for(int i = 0; i<A2; i++){
                
                xm = nucleusB.at(i).x;
                ym = nucleusB.at(i).y;
                
                if(param->getUseConstituentQuarkProton()>0){
                    T = 0.;
                    for (int iq=0; iq<param->getUseConstituentQuarkProton(); iq++){
                        bp2 = (xm+xq2[i][iq]-x)*(xm+xq2[i][iq]-x)+(ym+yq2[i][iq]-y)*(ym+yq2[i][iq]-y);
                        bp2 /= hbarc*hbarc;
                        
                        T += exp(-bp2/(2.*BGq))/(2.*PI*BGq)/double(param->getUseConstituentQuarkProton())*gaussB[i][iq];
                        //	      cerr << "B " << i << " " << iq << " " << gaussA[i][iq] << endl;
                        
                    }
                }
                else{
                    phi = nucleusB.at(i).phi;
                    
                    bp2 = (xm-x)*(xm-x)+(ym-y)*(ym-y) + xi*pow((xm-x)*cos(phi) + (ym-y)*sin(phi),2.);
                    bp2 /= hbarc*hbarc;
                    
                    T = sqrt(1+xi)*exp(-bp2/(2.*BG))/(2.*PI*BG)*gaussB[i][0]; // T_p in this cell for the current nucleon
                }
                
                lat->cells[pos]->setTpB(lat->cells[pos]->getTpB()+T/nucleiInAverage); // add up all T_p
                
                maxT=max(lat->cells[pos]->getTpB()+T,maxT);
                
            }
        }
    }
    
    //  cerr << "maximal used T=" << maxT << endl;
//    stringstream strNcoll_name;
//    strNcoll_name << OutDirectory << "/NcollListID" << RNGSEED << ".dat";
//    string Ncoll_name;  Ncoll_name = strNcoll_name.str();
//    ofstream foutNcoll(Ncoll_name.c_str(),ios::out);
    
    if(OUTPUT_FLAG==1 || OUTPUT_FLAG==2){
        stringstream strNucApos_name,strNucBpos_name;
        
        strNucApos_name << OutDirectory << "/NucleusAID" << RNGSEED << ".txt";
        strNucBpos_name << OutDirectory << "/NucleusBID" << RNGSEED << ".txt";
        
        string NucApos_name;  NucApos_name = strNucApos_name.str();
        string NucBpos_name;  NucBpos_name = strNucBpos_name.str();
        
        ofstream ProjPositions(NucApos_name.c_str(),ios::out);
        ofstream TargPositions(NucBpos_name.c_str(),ios::out);
        
        for(int i = 0; i<A1; i++){
            
            ProjPositions << nucleusA.at(i).x << " " << nucleusA.at(i).y << endl;
        }
        
        for(int i = 0; i<A2; i++){
            
            TargPositions << nucleusB.at(i).x << " " << nucleusB.at(i).y << endl;
        }
        
        // CLEAN-UP //
        ProjPositions.close();
        TargPositions.close();
        
    }

    if(param->getGaussianWounding() == 0){
        for (int i = 0; i<A1; i++){
            for (int j = 0 ; j<A2 ;j++){
                
                dx = nucleusB.at(j).x-nucleusA.at(i).x;
                dy = nucleusB.at(j).y-nucleusA.at(i).y;
                dij = dx*dx+dy*dy;
                if(dij<d2){
                    //foutNcoll << (nucleusB.at(j).x+nucleusA.at(i).x)/2. << " " << (nucleusB.at(j).y+nucleusA.at(i).y)/2. << endl;
                    
                    Ncoll++;
                    nucleusB.at(j).collided=1;
                    nucleusA.at(i).collided=1;
                }
            }
        }
    }
    else{
        
        double p;
        double G=0.92;
        double ran;
        
        for (int i = 0; i<A1; i++){
            for (int j = 0 ; j<A2 ;j++){
                
                dx = nucleusB.at(j).x-nucleusA.at(i).x;
                dy = nucleusB.at(j).y-nucleusA.at(i).y;
                dij = dx*dx+dy*dy;
                
                p = G * exp(-G*dij/d2); // Gaussian profile
                
                ran = random->genrand64_real1();
                
                if(ran<p){
                    //foutNcoll << (nucleusB.at(j).x+nucleusA.at(i).x)/2. << " " << (nucleusB.at(j).y+nucleusA.at(i).y)/2. << endl;

                    Ncoll++;
                    nucleusB.at(j).collided=1;
                    nucleusA.at(i).collided=1;
                }
            }
        }
    }

    // CLEAN-UP //
    //foutNcoll.close();
    
    //   for (int i = 0; i<A1; i++)
    //     {
    //       for (int j = 0 ; j<A2 ;j++)
    // 	{
    // 	  dx = nucleusB.at(j).x-nucleusA.at(i).x;
    // 	  dy = nucleusB.at(j).y-nucleusA.at(i).y;
    // 	  dij = dx*dx+dy*dy;
    // 	  if (dij < d2)
    // 	    {
    // 	      nucleusB.at(j).collided=1;
    // 	      nucleusA.at(i).collided=1;
    // 	    }
    // 	}
    //     }
    
    // in p+p assume that they collided in any case
    if ( A1 == 1 && A2 == 1 ){
        
        nucleusB.at(0).collided=1;
        nucleusA.at(0).collided=1;
    }
    
    // stringstream strgmuA_name;
    // strgmuA_name << "gmuA" << param->getMPIRank() << ".dat";
    // string gmuA_name;
    // gmuA_name = strgmuA_name.str();
    
    // ofstream fout(gmuA_name.c_str(),ios::out);
    double outvalue;
    double alphas;
    double Ydeviation = 10000;
    double QsA, QsB, distanceA, distanceB;
    
    QsA = 1;
    QsB = 1;
    
    Npart = 0;
    
    for(int i = 0; i<A1; i++){
        
        if(nucleusA.at(i).collided==1)
            Npart++;
    }
    
    for(int i = 0; i<A2; i++){
        if(nucleusB.at(i).collided==1)
            Npart++;
    }
    
    //   if ( Npart == 0 && param->getUseFixedNpart()==0)
    //     {
    //       cerr << "no collision happened. Exiting." << endl;
    //       exit(1);
    //     }
    
    param->setNpart(Npart);
    
    if(param->getUseFixedNpart()!=0 && Npart!=param->getUseFixedNpart()){
        cerr << "current Npart = " << Npart << endl;
        return;
    }
    
    // get Q_s^2 (and from that g^2mu^2) for a given \sum T_p and Y
    for(int ix=0; ix<N; ix++){ // loop over all positions
        x = -L/2.+a*ix;
        for(int iy=0; iy<N; iy++){
            check = 0;
            y = -L/2.+a*iy;
            Ydeviation = 10000;
            pos = ix*N+iy;
            
            // this version removes noise outside the interaction region
            // by checking whether we are inside a wounded nucleon
            // using 2 times the nucleon size as a cutoff
            // 	  for (int i = 0; i<A1; i++)
            // 	    {
            // 	      xm = nucleusA.at(i).x;
            // 	      ym = nucleusA.at(i).y;
            // 	      r = sqrt((x-xm)*(x-xm)+(y-ym)*(y-ym));
            
            // 	      if(r<sqrt(2.*0.1*param->getSigmaNN()/PI) && nucleusA.at(i).collided==1)
            // 		{
            // 		  check=1;
            // 		}
            // 	    }
            
            // 	  for (int i = 0; i<A2; i++)
            // 	    {
            // 	      xm = nucleusB.at(i).x;
            // 	      ym = nucleusB.at(i).y;
            // 	      r = sqrt((x-xm)*(x-xm)+(y-ym)*(y-ym));
            
            // 	      if(r<sqrt(2.*0.1*param->getSigmaNN()/PI) && nucleusB.at(i).collided==1 && check==1)
            // 		check=2;
            // 	    }
            
            
            // 	  // cut at a radius of ~1.8 fm
            // 	  for (int i = 0; i<A1; i++)
            // 	    {
            // 	      if(lat->cells[pos]->getTpA() > 0.0019 && nucleusA.at(i).collided==1)
            // 		{
            // 		  check=1;
            // 		}
            // 	    }
            
            // 	  for (int i = 0; i<A2; i++)
            // 	    {
            // 	      if(lat->cells[pos]->getTpB() > 0.0019 && nucleusB.at(i).collided==1 && check==1)
            // 		check=2;
            // 	    }
            
            
            // 	  // cut at a radius of ~1.8 fm
            // 	  if(lat->cells[pos]->getTpA() > 0.0019 )
            // 	    {
            // 	      check=1;
            // 	    }
            
            // 	  if(lat->cells[pos]->getTpB() > 0.0019 && check==1)
            // 	    check=2;
            
            
            
            //	  cut proton at a radius of rmax [fm] (about twice the gluonic radius to be generous)
            if(log(2*M_PI*BG*lat->cells[pos]->getTpA())<0.)
                distanceA = sqrt(-2.*BG*log(2*M_PI*BG*lat->cells[pos]->getTpA()))*hbarc;
            else distanceA=0.;
            
            if(log(2*M_PI*BG*lat->cells[pos]->getTpB())<0.)
                distanceB = sqrt(-2.*BG*log(2*M_PI*BG*lat->cells[pos]->getTpB()))*hbarc;
            else distanceB=0.;
            
            // if(distanceA>0.1)
            //   cerr << lat->cells[pos]->getTpA()<< " " << distanceA << " " << param->getRmax() << endl;
            
            if(distanceA < param->getRmax())
            {
                check=1;
            }
            
            if(distanceB < param->getRmax() && check==1)
            {
                check=2;
            }
            
            
            double exponent=5.6; // see 1212.2974 Eq. (17)
            if(check==2)
            {
                if ( param->getUseFluctuatingx() == 1)
                {
                    // iterative loops here to determine the fluctuating Y
                    // _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
                    while (abs(Ydeviation) > 0.001)
                    {
                        if(rapidity>=0)
                            QsA = sqrt(getNuclearQs2(param, random, lat->cells[pos]->getTpA(), abs(rapidity)));
                        else
                        {
                            xVal = QsA*param->getxFromThisFactorTimesQs()/param->getRoots()*exp(yIn);
                            //cerr << " QsA=" << QsA << ", param->getRoots()=" << param->getRoots() << ", exp(yIn)=" << exp(yIn) << endl;
                            if(xVal==0)
                                QsA=0.;
                            else
                                QsA = sqrt(getNuclearQs2(param, random, lat->cells[pos]->getTpA(), 0.))*
                                sqrt(pow((1-xVal)/(1-0.01),exponent)*pow((0.01/xVal),0.2));
                            //cerr << "xVal=" << xVal << endl;
                            //cerr << "QsA=" << QsA << endl;
                        }
                        if(QsA == 0)
                        {
                            Ydeviation = 0;
                            lat->cells[pos]->setg2mu2A(0.);
                        }
                        else
                        {
                            // nucleus A
                            lat->cells[pos]->setg2mu2A(QsA*QsA/param->getQsmuRatio()/param->getQsmuRatio()
                                                       *a*a/hbarc/hbarc/param->getg()); // lattice units? check
                            
                            
                            Ydeviation = rapidity - log(0.01/(QsA*param->getxFromThisFactorTimesQs()/param->getRoots()*exp(yIn)));
                            rapidity = log(0.01/(QsA*param->getxFromThisFactorTimesQs()/param->getRoots()*exp(yIn)));
                        }
                    }
                    if(lat->cells[pos]->getg2mu2A()!=lat->cells[pos]->getg2mu2A())
                    {
                        lat->cells[pos]->setg2mu2A(0.);
                    }
                    
                    // if(ix==N/2 && iy==N/2)
                    //   {
                    //     if(QsA*param->getxFromThisFactorTimesQs()/param->getRoots()*exp(yIn)<=0.01)
                    // 	cerr  << "rapidity_A=" << rapidity << endl;
                    //     else
                    // 	cerr  << "rapidity_A=" << log(0.01/xVal) << endl;
                    // 	cerr  << "xVal=" << xVal << endl;
                    //     cerr  << "Q_sA=" << QsA << endl;
                    //   }
                    
                    
                    
                    Ydeviation = 10000;
                    
                    while (abs(Ydeviation) > 0.001)
                    {
                        if(rapidity>=0)
                            QsB = sqrt(getNuclearQs2(param, random, lat->cells[pos]->getTpB(), abs(rapidity)));
                        else
                        {
                            xVal = QsB*param->getxFromThisFactorTimesQs()/param->getRoots()*exp(-yIn);
                            if(xVal==0)
                                QsB=0.;
                            else
                                QsB = sqrt(getNuclearQs2(param, random, lat->cells[pos]->getTpB(), 0.))*
                                sqrt(pow((1-xVal)/(1-0.01),exponent)*pow((0.01/xVal),0.2));
                        }
                        if(QsB == 0)
                        {
                            Ydeviation = 0;
                            lat->cells[pos]->setg2mu2B(0.);
                        }
                        else
                        {
                            // nucleus B
                            lat->cells[pos]->setg2mu2B(QsB*QsB/param->getQsmuRatioB()/param->getQsmuRatioB()
                                                       *a*a/hbarc/hbarc/param->getg()/param->getg());
                            
                            
                            //      cerr << ix << " " << iy << endl;
                            // 		      cerr << " QsB = " << QsB << " GeV" << endl;
                            // 		      cerr << " x= " << (QsB/2./param->getRoots()) << endl;
                            
                            Ydeviation = rapidity - log(0.01/(QsB*param->getxFromThisFactorTimesQs()/param->getRoots()*exp(-yIn)));
                            rapidity = log(0.01/(QsB*param->getxFromThisFactorTimesQs()/param->getRoots()*exp(-yIn)));
                        }
                    }
                    if(lat->cells[pos]->getg2mu2B()!=lat->cells[pos]->getg2mu2B())
                    {
                        lat->cells[pos]->setg2mu2B(0.);
                    }
                    // if(ix==N/2 && iy==N/2)
                    //   {
                    //     if(QsB*param->getxFromThisFactorTimesQs()/param->getRoots()*exp(-yIn)<=0.01)
                    //       cerr  << "rapidity_B=" << rapidity << endl;
                    //     else
                    //       cerr  << "rapidity_B=" << log(0.01/xVal) << endl;
                    //     cerr  << "Q_sB=" << QsB << endl;
                    //   }
                    
                    
                    // _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
                    // end iterative loops here
                }
                else
                {
                    // nucleus A
                    lat->cells[pos]->setg2mu2A(getNuclearQs2(param, random, lat->cells[pos]->getTpA(), rapidity)/param->getQsmuRatio()/param->getQsmuRatio()
                                               *a*a/hbarc/hbarc/param->getg()/param->getg()); // lattice units? check
                    
                    // nucleus B
                    lat->cells[pos]->setg2mu2B(getNuclearQs2(param, random, lat->cells[pos]->getTpB(), rapidity)/param->getQsmuRatioB()/param->getQsmuRatioB()
                                               *a*a/hbarc/hbarc/param->getg()/param->getg());
                    
                }
            }
        }
    }
    
    
    // output gmu
    count=0;
    count2=0;
    double Tpp=0.;
    
    for(int ix=0; ix<N; ix++) // loop over all positions
    {
        for(int iy=0; iy<N; iy++)
        {
            check = 0;
            pos = ix*N+iy;
            x = -L/2.+a*ix;
            y = -L/2.+a*iy;
            //  outvalue = sqrt(lat->cells[pos]->getg2mu2B())/a*hbarc; // in GeV
            outvalue = lat->cells[pos]->getg2mu2A();
            
            // posA = static_cast<int>(floor((x-b/2.+L/2.)/a+0.00000001))*N+iy;
            // posB = static_cast<int>(floor((x+b/2.+L/2.)/a+0.00000001))*N+iy;
            
            posA = pos;
            posB = pos;
            
            if(posA>0 && posA<(N-1)*N+N-1)
            {
                g2mu2A = lat->cells[posA]->getg2mu2A();
            }
            else
                g2mu2A = 0;
            
            if(posB>0 && posB<(N-1)*N+N-1)
            {
                g2mu2B = lat->cells[posB]->getg2mu2B();
            }
            else
                g2mu2B = 0;
            
            if(g2mu2B>=g2mu2A)
            {
                averageQs2min2 += g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*hbarc*hbarc*param->getg()*param->getg();
            }
            else
            {
                averageQs2min2 += g2mu2B*param->getQsmuRatioB()*param->getQsmuRatioB()/a/a*hbarc*hbarc*param->getg()*param->getg();
            }
            
            for (int i = 0; i<A1; i++)
            {
                xm = nucleusA.at(i).x;
                ym = nucleusA.at(i).y;
                r = sqrt((x-xm)*(x-xm)+(y-ym)*(y-ym));
                
                if(r<sqrt(0.1*param->getSigmaNN()/PI) && nucleusA.at(i).collided==1)
                {
                    check=1;
                }
            }
            
            for (int i = 0; i<A2; i++)
            {
                xm = nucleusB.at(i).x;
                ym = nucleusB.at(i).y;
                r = sqrt((x-xm)*(x-xm)+(y-ym)*(y-ym));
                
                if(r<sqrt(0.1*param->getSigmaNN()/PI) && nucleusB.at(i).collided==1 && check==1)
                    check=2;
            }
            
            if(check==2)
            {
                if(g2mu2B>g2mu2A)
                {
                    averageQs += sqrt(g2mu2B*param->getQsmuRatioB()*param->getQsmuRatioB()/a/a*hbarc*hbarc*param->getg()*param->getg());
                    averageQs2 += g2mu2B*param->getQsmuRatioB()*param->getQsmuRatioB()/a/a*hbarc*hbarc*param->getg()*param->getg();
                    averageQs2min += g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*hbarc*hbarc*param->getg()*param->getg();
                }
                else
                {
                    averageQs += sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*hbarc*hbarc*param->getg()*param->getg());
                    averageQs2 += g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*hbarc*hbarc*param->getg()*param->getg();
                    averageQs2min += g2mu2B*param->getQsmuRatioB()*param->getQsmuRatioB()/a/a*hbarc*hbarc*param->getg()*param->getg();
                }
                averageQs2Avg += (g2mu2B*param->getQsmuRatioB()*param->getQsmuRatioB()+g2mu2A*param->getQsmuRatio()*param->getQsmuRatio())
                /2./a/a*hbarc*hbarc*param->getg()*param->getg();
                count++;
            }
            //else
            //	    fout << x << " " << y << " " << 0. << endl;
            
            // compute T_pp
            
            Tpp += lat->cells[pos]->getTpB()*lat->cells[pos]->getTpA()*a*a/hbarc/hbarc/hbarc/hbarc; // now this quantity is in fm^-2
            // remember: Tp is in GeV^2
        }
        //      fout << endl;
    }
    //  fout.close();
    
    averageQs/=static_cast<double>(count);
    averageQs2/=static_cast<double>(count);
    averageQs2Avg/=static_cast<double>(count);
    averageQs2min/=static_cast<double>(count);
    
    param->setAverageQs(sqrt(averageQs2));
    param->setAverageQsAvg(sqrt(averageQs2Avg));
    param->setAverageQsmin(sqrt(averageQs2min));
    
    param->setTpp(Tpp);
    
    std::cerr << "# N_part=" << Npart << " N_coll=" << Ncoll << endl;
    std::cerr << "# T_pp(" << param->getb() << " fm) = " << Tpp << " 1/fm^2" << endl;
    std::cerr << "# Q_s^2(max) S_T = " << averageQs2*a*a/hbarc/hbarc*static_cast<double>(count) << endl;
    std::cerr << "# Q_s^2(avg) S_T = " << averageQs2Avg*a*a/hbarc/hbarc*static_cast<double>(count) << endl;
    std::cerr << "# Q_s^2(min) S_T = " << averageQs2min2*a*a/hbarc/hbarc << endl;
    
    std::cerr << "# Area = " << a*a*count << " fm^2" << endl;
    
    std::cerr << "# Average Qs(max) = " << param->getAverageQs() << " GeV" << endl;
    std::cerr << "# Average Qs(avg) = " << param->getAverageQsAvg() << " GeV" << endl;
    std::cerr << "# Average Qs(min) = " << param->getAverageQsmin() << " GeV" << endl;
    
    
    // OUTPUT OF SATURATION SCALES //
    stringstream FileName;
    FileName << OutDirectory << "/SaturationID" << RNGSEED << ".txt";
    string FileNameStr;
    FileNameStr=FileName.str();
    
    ofstream SaturationOutput(FileNameStr.c_str(),ios::app);
    
    SaturationOutput << "# N_part -- N_coll -- b[fm] -- Tpp(b)[1/fm^2] -- Q_s^2(max)S_T -- Q_s^2(avg)S_T -- Q_s^2(min)S_T -- Area[fm^2] -- AverageQs(max)[GeV] -- AverageQs(avg)[GeV] -- AverageQs(min)[GeV]" << endl;
    
    SaturationOutput    << Npart << " "
    << Ncoll << " "
    << param->getb() << " "
    << Tpp << " "
    << averageQs2*a*a/hbarc/hbarc*static_cast<double>(count) << " "
    << averageQs2Avg*a*a/hbarc/hbarc*static_cast<double>(count) << " "
    << averageQs2min2*a*a/hbarc/hbarc << " "
    << a*a*count << " "
    << param->getAverageQs() << " "
    << param->getAverageQsAvg() << " "
    << param->getAverageQsmin() << endl;
    
    SaturationOutput.close();
    
    //std::cerr << "# RESULTING Y(Qs(max)*" << param->getxFromThisFactorTimesQs() << ") = " << log(0.01/(param->getAverageQs()*param->getxFromThisFactorTimesQs()/param->getRoots())) << endl;
    //std::cerr << "# RESULTING Y(Qs(avg)*" << param->getxFromThisFactorTimesQs() << ") = "  << log(0.01/(param->getAverageQsAvg()*param->getxFromThisFactorTimesQs()/param->getRoots())) << endl;
    //std::cerr << "# RESULTING Y(Qs(min)*" << param->getxFromThisFactorTimesQs() << ") =  " << log(0.01/(param->getAverageQsmin()*param->getxFromThisFactorTimesQs()/param->getRoots())) << endl;
    
    std::cerr << "# COLOR CHARGE DENSITIES FOR BOTH NUCLEI SET" << endl;
    
    // stringstream strQs2ST_name;
    // strQs2ST_name << "Qs2ST" << param->getMPIRank() << ".dat";
    // string Qs2ST_name;
    // Qs2ST_name = strQs2ST_name.str();
    // ofstream foutQ(Qs2ST_name.c_str(),ios::out);
    // foutQ << b << " " << Npart << " " << averageQs2min2*a*a/hbarc/hbarc << " " << a*a*count << endl;
    // foutQ.close();
    
    
    
    if(param->getRunningCoupling() && param->getRunWithkt()==0)
    {
        if(param->getRunWithQs()==2)
        {
            //std::cerr << "running with " << param->getRunWithThisFactorTimesQs() << " Q_s(max)" << endl;
            alphas = 12.*PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*param->getAverageQs()/0.2)); // 3 flavors
            //std::cerr << "alpha_s(" << param->getRunWithThisFactorTimesQs() << " Qs_max)=" << alphas << endl;
        }
        else if(param->getRunWithQs()==0)
        {
            //std::cerr << "running with " << param->getRunWithThisFactorTimesQs() << " Q_s(min)" << endl;
            alphas = 12.*PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*param->getAverageQsmin()/0.2)); // 3 flavors
            //std::cerr << "alpha_s(" << param->getRunWithThisFactorTimesQs() << " Qs_min)=" << alphas << endl;
        }
        else if(param->getRunWithQs()==1)
        {
            //std::cerr << "running with " << param->getRunWithThisFactorTimesQs() << " <Q_s>" << endl;
            alphas = 12.*PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*param->getAverageQsAvg()/0.2)); // 3 flavors
            //std::cerr << "alpha_s(" << param->getRunWithThisFactorTimesQs() << " <Qs>)=" << alphas << endl;
        }
    }
    else if(param->getRunningCoupling() && param->getRunWithkt()==1)
    {
        //std::cerr << "Multiplicity with running alpha_s(k_T)" << endl;
    }
    else
    {
        //std::cerr << "Using fixed alpha_s" << endl;
        alphas = param->getg()*param->getg()/4./PI;
    }
    
    if(param->getAverageQs() > 0 && param->getAverageQsAvg()>0 && averageQs2>0  && param->getAverageQsmin()>0 && averageQs2Avg>0 && alphas>0 && Npart>=2){
        param->setSuccess(1);
        //std::cerr << "#SUCCESS WITH SETTING Qs " << std::endl;
    }
    param->setalphas(alphas);
    
    stringstream strup_name;
    strup_name << OutDirectory << "/usedParameters" << RNGSEED << ".dat";
    string up_name;
    up_name = strup_name.str();
    
    ofstream IPSatInfo(up_name.c_str(),ios::app);
    IPSatInfo << " " << endl;
    IPSatInfo << " Output by setColorChargeDensity in Init.cpp: " << endl;
    IPSatInfo << " " << endl;
    IPSatInfo << "b = " << b << " fm" << endl;
    IPSatInfo << "Npart = " << Npart << endl;
    IPSatInfo << "Ncoll = " << Ncoll << endl;
    if(param->getRunningCoupling())
    {
        if(param->getRunWithQs()==2)
            IPSatInfo << "<Q_s>(max) = " << param->getAverageQs() << endl;
        else if(param->getRunWithQs()==1)
            IPSatInfo << "<Q_s>(avg) = " << param->getAverageQsAvg() << endl;
        else if(param->getRunWithQs()==0)
            IPSatInfo << "<Q_s>(min) = " << param->getAverageQsmin() << endl;
        IPSatInfo << "alpha_s(" << param->getRunWithThisFactorTimesQs() << " <Q_s>) = " << param->getalphas() << endl;
    }
    else
        IPSatInfo << "using fixed coupling alpha_s=" << param->getalphas() << endl;
    //  IPSatInfo << "Q_s ~ x^-" << param->getxExponent() << endl;
    
    IPSatInfo << " " << endl;
    IPSatInfo << "Rnp=" << param->getRnp() << std::endl;
    IPSatInfo << " " << endl;
    IPSatInfo.close();
    
    //      // output gmu
    //      ofstream foutB("gmuB.dat",ios::out);
    //      for(int ix=0; ix<N; ix++) // loop over all positions
    //        {
    //          for(int iy=0; iy<N; iy++)
    //     	{
    //     	  pos = ix*N+iy;
    //     	  outvalue = sqrt(lat->cells[pos]->getg2mu2B());
    //     	  foutB << ix << " " << iy << " " << outvalue << endl;
    //     	}
    //          foutB << endl;
    //        }
    //      foutB.close();
    //     // output gmu
    //      ofstream foutA("gmuA.dat",ios::out);
    //      for(int ix=0; ix<N; ix++) // loop over all positions
    //        {
    //          for(int iy=0; iy<N; iy++)
    //     	{
    //     	  pos = ix*N+iy;
    //     	  outvalue = sqrt(lat->cells[pos]->getg2mu2A());
    //     	  foutA << ix << " " << iy << " " << outvalue << endl;
    //     	}
    //          foutA << endl;
    //        }
    //      foutA.close();
}


// set g^2\mu^2 as the sum of the individual nucleons' g^2\mu^2, using Q_s(b,y) prop tp g^mu(b,y)
// also compute N_part using Glauber
void Init::setColorChargeDensityMVTarget(GIPSLattice *lat,Parameters *param,Random *random,Glauber *glauber,string OutDirectory,int RNGSEED,int OUTPUT_FLAG){
    
    int pos,posA,posB;
    int N = param->getSize();
    int A1,A2;
    int check=0;
    double xVal;
    if(param->getNucleonPositionsFromFile()==0){
        A1 = static_cast<int>(glauber->nucleusA1())*param->getAverageOverNuclei();
        A2 = static_cast<int>(glauber->nucleusA2())*param->getAverageOverNuclei(); // DUMMY //

    }
    else{
        A1 = param->getA1FromFile();
        A2 = param->getA2FromFile();

    }
    
    int Npart = 0;
    int Ncoll = 0;
    double g2mu2A,g2mu2B;
    double b = param->getb();
    double x, xm;
    double y, ym;
    double r;
    double L = param->getL();
    double rapidity;
    double P,m;
    if(param->getUsePseudoRapidity()==0){
        rapidity = param->getRapidity();
    }
    else{
        // when using pseudorapidity as input convert to rapidity here. later include Jacobian in multiplicity and energy
        //std::cerr << "# USING PSEUDORAPIDITY " << param->getRapidity() << endl;
        m=param->getJacobianm(); // in GeV
        P=0.13+0.32*pow(param->getRoots()/1000.,0.115); //in GeV
        rapidity = 0.5 * log(sqrt(pow(cosh(param->getRapidity()),2.)+m*m/(P*P))+sinh(param->getRapidity())
                             / (sqrt(pow(cosh(param->getRapidity()),2.)+m*m/(P*P))-sinh(param->getRapidity())));
        //std::cerr << "# CORRESPONDS TO RAPIDITY " << rapidity << endl;
    }
    
    double yIn = rapidity;//param->getRapidity();
    double a = L/N; // lattice spacing in fm
    double dx, dy, dij;
    double d2 = param->getSigmaNN()/(PI*10.);          // in fm^2
    double averageQs = 0.;
    double averageQs2 = 0.;
    double averageQs2Avg = 0.;
    double averageQs2min = 0.;
    double averageQs2min2 = 0.;
    int count = 0;
    int count2 = 0;
    double nucleiInAverage;
    nucleiInAverage = static_cast<double>(param->getAverageOverNuclei());
    
    
    double gaussA[A1][3];
    
    for(int i = 0; i<A1; i++){
        for(int iq = 0; iq<3; iq++){
            
            gaussA[i][iq]=1.;
        }
    }
    
    // let the log fluctuate
    if(param->getSmearQs() == 1){
        if(A1>0){
            for(int i = 0; i<A1; i++){
                for(int iq = 0; iq<3; iq++){
                    
                    gaussA[i][iq] = (exp(random->Gauss(0,param->getSmearingWidth())))/1.13; // dividing by 1.13 restores the same mean Q_s
                    //cerr << i << " " << iq << " " << gaussA[i][iq] << endl;
                    //      if (gaussA[i]<0)
                    //  gaussA[i]=0.;
                }
            }
        }
    }
    
    param->setQsmuRatioB(param->getQsmuRatio());

    // MV TARGET INITIAL CONDITIONS -- INFINITE NUCLEUS //
    for(int ix=0; ix<N; ix++){ // loop over all positions
        for(int iy=0; iy<N; iy++){
            pos = ix*N+iy;
            lat->cells[pos]->setg2mu2B(param->getg2mu()*param->getg2mu()/param->getg()/param->getg());
        }
    }
    param->setSuccess(1);
    cerr << "# CONSTANT COLOR CHARGE DENSITY TARGET SET" << endl;
    
    // ERASE PROJECTILE CHARGE DENSITY //
    for(int ix=0; ix<N; ix++){ // loop over all positions
        for(int iy=0; iy<N; iy++){
            
            pos = ix*N+iy;

            lat->cells[pos]->setg2mu2A(0.);
        }
    }

    double maxT=0;
    
    double bp2,T,BG;
    BG = param->getBG();
    double BGq = param->getBGq(); // quark size in GeV^-2
    double xi = param->getProtonAnisotropy();
    double phi;
    
    if(xi!=0.0){
        for (int i = 0; i<A1; i++){
            nucleusA.at(i).phi = 2*M_PI*random->genrand64_real2();
        }
    }
    
    //  cerr << "BG=" << BG << endl;
    
    double xq[A1][param->getUseConstituentQuarkProton()], xq2[A2][param->getUseConstituentQuarkProton()];
    double yq[A1][param->getUseConstituentQuarkProton()], yq2[A2][param->getUseConstituentQuarkProton()];
    double avgxq=0.;
    double avgyq=0.;
    
    if(param->getUseConstituentQuarkProton()>0){
        for (int i=0; i<A1; i++){
            avgxq=0.;
            avgyq=0.;
            for(int iq=0; iq<param->getUseConstituentQuarkProton(); iq++){
                xq[i][iq] = sqrt(BG*hbarc*hbarc)*random->Gauss();
                yq[i][iq] = sqrt(BG*hbarc*hbarc)*random->Gauss();
            }
            for(int iq=0; iq<3; iq++){
                avgxq += xq[i][iq];
                avgyq += yq[i][iq];
            }
            for(int iq=0; iq<3; iq++){
                xq[i][iq] -= avgxq/3.;
                yq[i][iq] -= avgyq/3.;
                //          cerr << xq[i][iq] << " " << yq[i][iq] << endl;
            }
            
            
            // avgyq=0.;
            // for (int iq=0; iq<3; iq++)
            //   {
            //     avgxq += xq[i][iq];
            //     avgyq += yq[i][iq];
            //   }
            // cerr << avgyq << endl;
        }
    }
    
    if(param->getUseConstituentQuarkProton()>0){
        for(int i=0; i<A2; i++){
            avgxq=0.;
            avgyq=0.;
            for(int iq=0; iq<param->getUseConstituentQuarkProton(); iq++){
                xq2[i][iq] = sqrt(BG*hbarc*hbarc)*random->Gauss();
                yq2[i][iq] = sqrt(BG*hbarc*hbarc)*random->Gauss();
            }
            
            for(int iq=0; iq<3; iq++){
                avgxq += xq2[i][iq];
                avgyq += yq2[i][iq];
            }
            for(int iq=0; iq<3; iq++){
                xq2[i][iq] -= avgxq/3.;
                yq2[i][iq] -= avgyq/3.;
            }
            // avgyq=0.;
            // for (int iq=0; iq<3; iq++)
            //   {
            //     avgxq += xq2[i][iq];
            //     avgyq += yq2[i][iq];
            //   }
            // cerr << avgyq << endl;
        }
        
        
    }
    
    //add all T_p's (new in version 1.2)
    for(int ix=0; ix<N; ix++){ // loop over all positions
        
        x = -L/2.+a*ix;
        for(int iy=0; iy<N; iy++){
            y = -L/2.+a*iy;
            
            pos = ix*N+iy;
            
            // nucleus A
            lat->cells[pos]->setTpA(0.);
            for (int i = 0; i<A1; i++){
                xm = nucleusA.at(i).x;
                ym = nucleusA.at(i).y;
                
                if(param->getUseConstituentQuarkProton()>0){
                    T = 0.;
                    for(int iq=0; iq<param->getUseConstituentQuarkProton(); iq++){
                        bp2 = (xm+xq[i][iq]-x)*(xm+xq[i][iq]-x)+(ym+yq[i][iq]-y)*(ym+yq[i][iq]-y);
                        bp2 /= hbarc*hbarc;
                        
                        T += exp(-bp2/(2.*BGq))/(2.*PI*BGq)/(double(param->getUseConstituentQuarkProton()))*gaussA[i][iq]; // I removed the 2/3 here to make it a bit bigger
                        //          cerr << "A " << i << " " << iq << " " << gaussA[i][iq] << endl;
                    }
                }
                else{
                    phi = nucleusA.at(i).phi;
                    
                    bp2 = (xm-x)*(xm-x)+(ym-y)*(ym-y) + xi*pow((xm-x)*cos(phi) + (ym-y)*sin(phi),2.);
                    bp2 /= hbarc*hbarc;
                    
                    T = sqrt(1+xi)*exp(-bp2/(2.*BG))/(2.*PI*BG)*gaussA[i][0]; // T_p in this cell for the current nucleon
                }
                
                lat->cells[pos]->setTpA(lat->cells[pos]->getTpA()+T/nucleiInAverage); // add up all T_p
                
                maxT=max(lat->cells[pos]->getTpA()+T,maxT);
                
            }
            
//            // nucleus B
//            lat->cells[pos]->setTpB(0.);
//            for(int i = 0; i<A2; i++){
//
//                xm = nucleusB.at(i).x;
//                ym = nucleusB.at(i).y;
//
//                if(param->getUseConstituentQuarkProton()>0){
//                    T = 0.;
//                    for (int iq=0; iq<param->getUseConstituentQuarkProton(); iq++){
//                        bp2 = (xm+xq2[i][iq]-x)*(xm+xq2[i][iq]-x)+(ym+yq2[i][iq]-y)*(ym+yq2[i][iq]-y);
//                        bp2 /= hbarc*hbarc;
//
//                        T += exp(-bp2/(2.*BGq))/(2.*PI*BGq)/double(param->getUseConstituentQuarkProton())*gaussB[i][iq];
//                        //          cerr << "B " << i << " " << iq << " " << gaussA[i][iq] << endl;
//
//                    }
//                }
//                else{
//                    phi = nucleusB.at(i).phi;
//
//                    bp2 = (xm-x)*(xm-x)+(ym-y)*(ym-y) + xi*pow((xm-x)*cos(phi) + (ym-y)*sin(phi),2.);
//                    bp2 /= hbarc*hbarc;
//
//                    T = sqrt(1+xi)*exp(-bp2/(2.*BG))/(2.*PI*BG)*gaussB[i][0]; // T_p in this cell for the current nucleon
//                }
//
//                lat->cells[pos]->setTpB(lat->cells[pos]->getTpB()+T/nucleiInAverage); // add up all T_p
//
//                maxT=max(lat->cells[pos]->getTpB()+T,maxT);
//
//            }
        }
    }
    
    //  cerr << "maximal used T=" << maxT << endl;
//    stringstream strNcoll_name;
//    strNcoll_name << OutDirectory << "/NcollListID" << RNGSEED << ".dat";
//    string Ncoll_name;  Ncoll_name = strNcoll_name.str();
//    ofstream foutNcoll(Ncoll_name.c_str(),ios::out);

    if(OUTPUT_FLAG==1 || OUTPUT_FLAG==2){
        stringstream strNucApos_name;
        
        strNucApos_name << OutDirectory << "/NucleusAID" << RNGSEED << ".txt";
        
        string NucApos_name;  NucApos_name = strNucApos_name.str();
        
        ofstream ProjPositions(NucApos_name.c_str(),ios::out);
        
        for(int i = 0; i<A1; i++){
            
            ProjPositions << nucleusA.at(i).x << " " << nucleusA.at(i).y << endl;
        }

        // CLEAN-UP //
        ProjPositions.close();
        
    }
    
    if(param->getGaussianWounding() == 0){
        for(int i = 0; i<A1; i++){
            nucleusA.at(i).collided=1;
        }
    }
    else{
        
        double p;
        double G=0.92;
        double ran;
        
        for(int i = 0; i<A1; i++){
        
            nucleusA.at(i).collided=1;
        }
    }
    
    // CLEAN-UP //
    //foutNcoll.close();
    
    // in p+p assume that they collided in any case
    if(A1 == 1 && A2 == 1){
        
        nucleusB.at(0).collided=1;
        nucleusA.at(0).collided=1;
    }
    
    // stringstream strgmuA_name;
    // strgmuA_name << "gmuA" << param->getMPIRank() << ".dat";
    // string gmuA_name;
    // gmuA_name = strgmuA_name.str();
    
    // ofstream fout(gmuA_name.c_str(),ios::out);
    double outvalue;
    double alphas;
    double Ydeviation = 10000;
    double QsA, QsB, distanceA, distanceB;
    
    QsA = 1;
    QsB = 1;
    
    Npart = 0;
    
    for(int i = 0; i<A1; i++){
        
        if(nucleusA.at(i).collided==1)
            Npart++;
    }
    
    param->setNpart(Npart);
    
    if(param->getUseFixedNpart()!=0 && Npart!=param->getUseFixedNpart()){
        cerr << "current Npart = " << Npart << endl;
        return;
    }
    
    // get Q_s^2 (and from that g^2mu^2) for a given \sum T_p and Y
    for(int ix=0; ix<N; ix++){ // loop over all positions
        x = -L/2.+a*ix;
        for(int iy=0; iy<N; iy++){
            check = 0;
            y = -L/2.+a*iy;
            Ydeviation = 10000;
            pos = ix*N+iy;
            
            // this version removes noise outside the interaction region
            // by checking whether we are inside a wounded nucleon
            // using 2 times the nucleon size as a cutoff
            //       for (int i = 0; i<A1; i++)
            //         {
            //           xm = nucleusA.at(i).x;
            //           ym = nucleusA.at(i).y;
            //           r = sqrt((x-xm)*(x-xm)+(y-ym)*(y-ym));
            
            //           if(r<sqrt(2.*0.1*param->getSigmaNN()/PI) && nucleusA.at(i).collided==1)
            //         {
            //           check=1;
            //         }
            //         }
            
            //       for (int i = 0; i<A2; i++)
            //         {
            //           xm = nucleusB.at(i).x;
            //           ym = nucleusB.at(i).y;
            //           r = sqrt((x-xm)*(x-xm)+(y-ym)*(y-ym));
            
            //           if(r<sqrt(2.*0.1*param->getSigmaNN()/PI) && nucleusB.at(i).collided==1 && check==1)
            //         check=2;
            //         }
            
            
            //       // cut at a radius of ~1.8 fm
            //       for (int i = 0; i<A1; i++)
            //         {
            //           if(lat->cells[pos]->getTpA() > 0.0019 && nucleusA.at(i).collided==1)
            //         {
            //           check=1;
            //         }
            //         }
            
            //       for (int i = 0; i<A2; i++)
            //         {
            //           if(lat->cells[pos]->getTpB() > 0.0019 && nucleusB.at(i).collided==1 && check==1)
            //         check=2;
            //         }
            
            
            //       // cut at a radius of ~1.8 fm
            //       if(lat->cells[pos]->getTpA() > 0.0019 )
            //         {
            //           check=1;
            //         }
            
            //       if(lat->cells[pos]->getTpB() > 0.0019 && check==1)
            //         check=2;
            
            
            
            //      cut proton at a radius of rmax [fm] (about twice the gluonic radius to be generous)
            if(log(2*M_PI*BG*lat->cells[pos]->getTpA())<0.)
                distanceA = sqrt(-2.*BG*log(2*M_PI*BG*lat->cells[pos]->getTpA()))*hbarc;
            else distanceA=0.;
            
            
            // if(distanceA>0.1)
            //   cerr << lat->cells[pos]->getTpA()<< " " << distanceA << " " << param->getRmax() << endl;
            
            if(distanceA < param->getRmax()){

                check=2; // THIS SHOULD ALWAYS BE THE CASE //
            }
            
            
            double exponent=5.6; // see 1212.2974 Eq. (17)
            if(check==2){
                if ( param->getUseFluctuatingx() == 1){
                    
                    // iterative loops here to determine the fluctuating Y
                    // _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
                    while(abs(Ydeviation) > 0.001){
                        
                        if(rapidity>=0)
                            QsA = sqrt(getNuclearQs2(param, random, lat->cells[pos]->getTpA(), abs(rapidity)));
                        else{
                            xVal = QsA*param->getxFromThisFactorTimesQs()/param->getRoots()*exp(yIn);
                            //cerr << " QsA=" << QsA << ", param->getRoots()=" << param->getRoots() << ", exp(yIn)=" << exp(yIn) << endl;
                            if(xVal==0)
                                QsA=0.;
                            else
                                QsA = sqrt(getNuclearQs2(param, random, lat->cells[pos]->getTpA(), 0.))*
                                sqrt(pow((1-xVal)/(1-0.01),exponent)*pow((0.01/xVal),0.2));
                            //cerr << "xVal=" << xVal << endl;
                            //cerr << "QsA=" << QsA << endl;
                        }
                        if(QsA == 0){
                            Ydeviation = 0;
                            lat->cells[pos]->setg2mu2A(0.);
                        }
                        else{
                            // nucleus A
                            lat->cells[pos]->setg2mu2A(QsA*QsA/param->getQsmuRatio()/param->getQsmuRatio()
                                                       *a*a/hbarc/hbarc/param->getg()); // lattice units? check
                            
                            
                            Ydeviation = rapidity - log(0.01/(QsA*param->getxFromThisFactorTimesQs()/param->getRoots()*exp(yIn)));
                            rapidity = log(0.01/(QsA*param->getxFromThisFactorTimesQs()/param->getRoots()*exp(yIn)));
                        }
                    }
                    if(lat->cells[pos]->getg2mu2A()!=lat->cells[pos]->getg2mu2A()){
                        lat->cells[pos]->setg2mu2A(0.);
                    }
                    
                    // if(ix==N/2 && iy==N/2)
                    //   {
                    //     if(QsA*param->getxFromThisFactorTimesQs()/param->getRoots()*exp(yIn)<=0.01)
                    //     cerr  << "rapidity_A=" << rapidity << endl;
                    //     else
                    //     cerr  << "rapidity_A=" << log(0.01/xVal) << endl;
                    //     cerr  << "xVal=" << xVal << endl;
                    //     cerr  << "Q_sA=" << QsA << endl;
                    //   }
                    
                    Ydeviation = 10000;

                }
                else
                {
                    // nucleus A
                    lat->cells[pos]->setg2mu2A(getNuclearQs2(param, random, lat->cells[pos]->getTpA(), rapidity)/param->getQsmuRatio()/param->getQsmuRatio()
                                               *a*a/hbarc/hbarc/param->getg()/param->getg()); // lattice units? check
                    
                }
            }
        }
    }
    
    // output gmu
    count=0;
    count2=0;
    double Tpp=0.;
    
    for(int ix=0; ix<N; ix++){ // loop over all positions
        for(int iy=0; iy<N; iy++){
        
            check = 0;
            pos = ix*N+iy;
            x = -L/2.+a*ix;
            y = -L/2.+a*iy;
            //  outvalue = sqrt(lat->cells[pos]->getg2mu2B())/a*hbarc; // in GeV
            outvalue = lat->cells[pos]->getg2mu2A();
            
            // posA = static_cast<int>(floor((x-b/2.+L/2.)/a+0.00000001))*N+iy;
            // posB = static_cast<int>(floor((x+b/2.+L/2.)/a+0.00000001))*N+iy;
            
            posA = pos;
            posB = pos;
            
            if(posA>0 && posA<(N-1)*N+N-1){
                g2mu2A = lat->cells[posA]->getg2mu2A();
            }
            else{
                g2mu2A = 0;
            }
            
            if(posB>0 && posB<(N-1)*N+N-1){
                g2mu2B = param->getg2mu()*param->getg2mu();//lat->cells[posB]->getg2mu2B();
            }
            else{
                g2mu2B = 0;
            }
            
            if(g2mu2B>=g2mu2A){
                averageQs2min2 += g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*hbarc*hbarc*param->getg()*param->getg();
            }
            else{
                averageQs2min2 += g2mu2B*param->getQsmuRatioB()*param->getQsmuRatioB()/a/a*hbarc*hbarc*param->getg()*param->getg();
            }
            
            for (int i = 0; i<A1; i++){
                
                xm = nucleusA.at(i).x;
                ym = nucleusA.at(i).y;
                r = sqrt((x-xm)*(x-xm)+(y-ym)*(y-ym));
                
                if(r<sqrt(0.1*param->getSigmaNN()/PI) && nucleusA.at(i).collided==1){
                    check=2;
                }
            }
            
//            for (int i = 0; i<A2; i++){
//
//                xm = nucleusB.at(i).x;
//                ym = nucleusB.at(i).y;
//                r = sqrt((x-xm)*(x-xm)+(y-ym)*(y-ym));
//
//                if(r<sqrt(0.1*param->getSigmaNN()/PI) && nucleusB.at(i).collided==1 && check==1)
//                    check=2;
//            }
            
            if(check==2){
                if(g2mu2B>g2mu2A){
                    averageQs += sqrt(g2mu2B*param->getQsmuRatioB()*param->getQsmuRatioB()/a/a*hbarc*hbarc*param->getg()*param->getg());
                    averageQs2 += g2mu2B*param->getQsmuRatioB()*param->getQsmuRatioB()/a/a*hbarc*hbarc*param->getg()*param->getg();
                    averageQs2min += g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*hbarc*hbarc*param->getg()*param->getg();
                }
                else{
                    averageQs += sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*hbarc*hbarc*param->getg()*param->getg());
                    averageQs2 += g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*hbarc*hbarc*param->getg()*param->getg();
                    averageQs2min += g2mu2B*param->getQsmuRatioB()*param->getQsmuRatioB()/a/a*hbarc*hbarc*param->getg()*param->getg();
                }
                averageQs2Avg += (g2mu2B*param->getQsmuRatioB()*param->getQsmuRatioB()+g2mu2A*param->getQsmuRatio()*param->getQsmuRatio())/2./a/a*hbarc*hbarc*param->getg()*param->getg();
                count++;
            }
            //else
            //        fout << x << " " << y << " " << 0. << endl;
            
            // compute T_pp
            
//            Tpp += lat->cells[pos]->getTpB()*lat->cells[pos]->getTpA()*a*a/hbarc/hbarc/hbarc/hbarc; // now this quantity is in fm^-2
            Tpp += 1.0*lat->cells[pos]->getTpA()*a*a/hbarc/hbarc/hbarc/hbarc; // now this quantity is in fm^-2

            // remember: Tp is in GeV^2
        }
        //      fout << endl;
    }
    //  fout.close();
    
    averageQs/=static_cast<double>(count);
    averageQs2/=static_cast<double>(count);
    averageQs2Avg/=static_cast<double>(count);
    averageQs2min/=static_cast<double>(count);
    
    param->setAverageQs(sqrt(averageQs2));
    param->setAverageQsAvg(sqrt(averageQs2Avg));
    param->setAverageQsmin(sqrt(averageQs2min));
    
    param->setTpp(Tpp);
    
    std::cerr << "# N_part=" << Npart << " N_coll=" << Ncoll << endl;
    std::cerr << "# T_pp(" << param->getb() << " fm) = " << Tpp << " 1/fm^2" << endl;
    std::cerr << "# Q_s^2(max) S_T = " << averageQs2*a*a/hbarc/hbarc*static_cast<double>(count) << endl;
    std::cerr << "# Q_s^2(avg) S_T = " << averageQs2Avg*a*a/hbarc/hbarc*static_cast<double>(count) << endl;
    std::cerr << "# Q_s^2(min) S_T = " << averageQs2min2*a*a/hbarc/hbarc << endl;
    
    std::cerr << "# Area = " << a*a*count << " fm^2" << endl;
    
    std::cerr << "# Average Qs(max) = " << param->getAverageQs() << " GeV" << endl;
    std::cerr << "# Average Qs(avg) = " << param->getAverageQsAvg() << " GeV" << endl;
    std::cerr << "# Average Qs(min) = " << param->getAverageQsmin() << " GeV" << endl;
    
    
    // OUTPUT OF SATURATION SCALES //
    stringstream FileName;
    FileName << OutDirectory << "/SaturationID" << RNGSEED << ".txt";
    string FileNameStr;
    FileNameStr=FileName.str();
    
    ofstream SaturationOutput(FileNameStr.c_str(),ios::app);
    
    SaturationOutput << "# N_part -- N_coll -- b[fm] -- Tpp(b)[1/fm^2] -- Q_s^2(max)S_T -- Q_s^2(avg)S_T -- Q_s^2(min)S_T -- Area[fm^2] -- AverageQs(max)[GeV] -- AverageQs(avg)[GeV] -- AverageQs(min)[GeV]" << endl;
    
    SaturationOutput    << Npart << " "
    << Ncoll << " "
    << param->getb() << " "
    << Tpp << " "
    << averageQs2*a*a/hbarc/hbarc*static_cast<double>(count) << " "
    << averageQs2Avg*a*a/hbarc/hbarc*static_cast<double>(count) << " "
    << averageQs2min2*a*a/hbarc/hbarc << " "
    << a*a*count << " "
    << param->getAverageQs() << " "
    << param->getAverageQsAvg() << " "
    << param->getAverageQsmin() << endl;
    
    SaturationOutput.close();
    
    //std::cerr << "# RESULTING Y(Qs(max)*" << param->getxFromThisFactorTimesQs() << ") = " << log(0.01/(param->getAverageQs()*param->getxFromThisFactorTimesQs()/param->getRoots())) << endl;
    //std::cerr << "# RESULTING Y(Qs(avg)*" << param->getxFromThisFactorTimesQs() << ") = "  << log(0.01/(param->getAverageQsAvg()*param->getxFromThisFactorTimesQs()/param->getRoots())) << endl;
    //std::cerr << "# RESULTING Y(Qs(min)*" << param->getxFromThisFactorTimesQs() << ") =  " << log(0.01/(param->getAverageQsmin()*param->getxFromThisFactorTimesQs()/param->getRoots())) << endl;
    
    std::cerr << "# COLOR CHARGE DENSITIES FOR BOTH NUCLEI SET" << endl;
    
    // stringstream strQs2ST_name;
    // strQs2ST_name << "Qs2ST" << param->getMPIRank() << ".dat";
    // string Qs2ST_name;
    // Qs2ST_name = strQs2ST_name.str();
    // ofstream foutQ(Qs2ST_name.c_str(),ios::out);
    // foutQ << b << " " << Npart << " " << averageQs2min2*a*a/hbarc/hbarc << " " << a*a*count << endl;
    // foutQ.close();
    
    if(param->getRunningCoupling() && param->getRunWithkt()==0)
    {
        if(param->getRunWithQs()==2)
        {
            //std::cerr << "running with " << param->getRunWithThisFactorTimesQs() << " Q_s(max)" << endl;
            alphas = 12.*PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*param->getAverageQs()/0.2)); // 3 flavors
            //std::cerr << "alpha_s(" << param->getRunWithThisFactorTimesQs() << " Qs_max)=" << alphas << endl;
        }
        else if(param->getRunWithQs()==0)
        {
            //std::cerr << "running with " << param->getRunWithThisFactorTimesQs() << " Q_s(min)" << endl;
            alphas = 12.*PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*param->getAverageQsmin()/0.2)); // 3 flavors
            //std::cerr << "alpha_s(" << param->getRunWithThisFactorTimesQs() << " Qs_min)=" << alphas << endl;
        }
        else if(param->getRunWithQs()==1)
        {
            //std::cerr << "running with " << param->getRunWithThisFactorTimesQs() << " <Q_s>" << endl;
            alphas = 12.*PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*param->getAverageQsAvg()/0.2)); // 3 flavors
            //std::cerr << "alpha_s(" << param->getRunWithThisFactorTimesQs() << " <Qs>)=" << alphas << endl;
        }
    }
    else if(param->getRunningCoupling() && param->getRunWithkt()==1)
    {
        //std::cerr << "Multiplicity with running alpha_s(k_T)" << endl;
    }
    else
    {
        //std::cerr << "Using fixed alpha_s" << endl;
        alphas = param->getg()*param->getg()/4./PI;
    }
    
    if(param->getAverageQs() > 0 && param->getAverageQsAvg()>0 && averageQs2>0  && param->getAverageQsmin()>0 && averageQs2Avg>0 && alphas>0 && Npart>=2){
        param->setSuccess(1);
        //std::cerr << "#SUCCESS WITH SETTING Qs " << std::endl;
    }
    param->setalphas(alphas);
    
    stringstream strup_name;
    strup_name << OutDirectory << "/usedParameters" << RNGSEED << ".dat";
    string up_name;
    up_name = strup_name.str();
    
    ofstream IPSatInfo(up_name.c_str(),ios::app);
    IPSatInfo << " " << endl;
    IPSatInfo << " Output by setColorChargeDensityMVTarget in Init.cpp: " << endl;
    IPSatInfo << " " << endl;
    IPSatInfo << "b = " << b << " fm" << endl;
    IPSatInfo << "Npart = " << Npart << endl;
    IPSatInfo << "Ncoll = " << Ncoll << endl;
    if(param->getRunningCoupling())
    {
        if(param->getRunWithQs()==2)
            IPSatInfo << "<Q_s>(max) = " << param->getAverageQs() << endl;
        else if(param->getRunWithQs()==1)
            IPSatInfo << "<Q_s>(avg) = " << param->getAverageQsAvg() << endl;
        else if(param->getRunWithQs()==0)
            IPSatInfo << "<Q_s>(min) = " << param->getAverageQsmin() << endl;
        IPSatInfo << "alpha_s(" << param->getRunWithThisFactorTimesQs() << " <Q_s>) = " << param->getalphas() << endl;
    }
    else
        IPSatInfo << "using fixed coupling alpha_s=" << param->getalphas() << endl;
    //  IPSatInfo << "Q_s ~ x^-" << param->getxExponent() << endl;
    
    IPSatInfo << " " << endl;
    IPSatInfo << "Rnp=" << param->getRnp() << std::endl;
    IPSatInfo << " " << endl;
    IPSatInfo.close();
    
    //      // output gmu
    //      ofstream foutB("gmuB.dat",ios::out);
    //      for(int ix=0; ix<N; ix++) // loop over all positions
    //        {
    //          for(int iy=0; iy<N; iy++)
    //         {
    //           pos = ix*N+iy;
    //           outvalue = sqrt(lat->cells[pos]->getg2mu2B());
    //           foutB << ix << " " << iy << " " << outvalue << endl;
    //         }
    //          foutB << endl;
    //        }
    //      foutB.close();
    //     // output gmu
    //      ofstream foutA("gmuA.dat",ios::out);
    //      for(int ix=0; ix<N; ix++) // loop over all positions
    //        {
    //          for(int iy=0; iy<N; iy++)
    //         {
    //           pos = ix*N+iy;
    //           outvalue = sqrt(lat->cells[pos]->getg2mu2A());
    //           foutA << ix << " " << iy << " " << outvalue << endl;
    //         }
    //          foutA << endl;
    //        }
    //      foutA.close();
} // END setColorChargeDensityMVTarget //


void Init::init(GIPSLattice *lat,Group *group,Parameters *param,Random *random,Glauber *glauber,int READFROMFILE,string OutDirectory,int RNGSEED,int OUTPUT_FLAG){
    
    int maxIterations = 100000;
    int N = param->getSize();
    int locNc = param->getNc();
    int bins = param->getSize();
    int ir;
    int count[bins];
    int Nc2m1 = locNc*locNc-1;
    int nn[2];
    int pos, pos1, pos2, pos3, posx, posy, posxm, posym, posxmym;
    int counts, countMe;
    int checkConvergence;
    int alphaCheck;
    int bShift; // number of cells to be shifted by due to impact parameter
    int posU;
    
    double dNc = static_cast<double>(locNc);
    double Fold;
    double Fnew;
    double r;
    double L = param->getL();
    double x;
    double y;
    double a = L/N; // lattice spacing in fm
    std::cerr << "# INITIALIZING FIELDS " << std::endl;
    //param->setRnp(0.); // NOW IN COMMANDLINE //
    
    double b;
    double bmin=param->getbmin();
    double bmax=param->getbmax();
    
    double xb = random->genrand64_real1(); // uniformly distributed random variable
    
    if(param->getUseNucleus() == 0 || param->getUseNucleus() == 2) // use b=0 fm for the constant g^2 mu case
    {
        param->setSuccess(1);
        b=0.;
        std::cerr << "# SETTING b=0 FOR CONSTANT COLOR CHARGE" << std::endl;
    }
    else
    {
        if(param->getLinearb()==1) // use a linear probability distribution for b if we are doing nuclei
        {
            b = sqrt((bmax*bmax-bmin*bmin)*xb+bmin*bmin);
            std::cerr << "# SAMPLING bImpact LINEARLY BETWEEN " << bmin << " and " << bmax << " fm -- SET bImpact=" << b << " fm" << std::endl;
            
        }
        else // use a uniform distribution instead
        {
            b = (bmax-bmin)*xb+bmin;
            std::cerr << "# SAMPLING bImpact UNIFORMLY BETWEEN " << bmin << " and " << bmax << " fm -- SET bImpact=" << b << " fm" << std::endl;
            
        }
    }
    
    param->setb(b);
    
    double Qs2G;
    double temp3;
    double g2mu;
    double trATB[N*N];
    double dr=a;
    double rtrAT2[bins];
    double epsilon;
    double lambda;
    double trATA[N*N];
    double avgEps;
    double avgEpsMag;
    double avgEpsEl;
    
    //  cerr << "m_lat =" << m << endl;
    
    
    // read Q_s^2 from file
    if(param->getUseNucleus() == 1 || param->getUseNucleus() == 2)
    {
        readNuclearQs(param);
    }
    
    
    nn[0]=N;
    nn[1]=N;
    
    // sample nucleon positions
    nucleusA.clear();
    nucleusB.clear();
    
    // to read Wilson lines from file (e.g. after JIMWLK evolution for the 3DGlasma)
    if(READFROMFILE){
        
        //readV(lat, group, param);
        //param->setSuccess(1);
        exit(0);
    }
    // to generate your own Wilson lines
    else{
        
        if(param->getUseNucleus() == 1){
            sampleTA(param, random, glauber);                           // populate the lists nucleusA and nucleusB with position data of the
        }
        
        if(param->getUseNucleus() == 2){
            sampleProjectileTA(param, random, glauber);                           // populate the lists nucleusA with position data of the
        }
        
        // set color charge densities
        if(param->getUseNucleus() == 0 || param->getUseNucleus() == 1){
            setColorChargeDensity(lat, param, random, glauber,OutDirectory,RNGSEED,OUTPUT_FLAG);
        }
        
        if(param->getUseNucleus() == 2){
            setColorChargeDensityMVTarget(lat, param, random, glauber,OutDirectory,RNGSEED,OUTPUT_FLAG);
        }
        
        if(param->getUseNucleus() == 1 && param->getUseFixedNpart()!=0 && param->getNucleonPositionsFromFile()!=1){
            if(param->getNpart()!=param->getUseFixedNpart()){
                while(param->getNpart()!=param->getUseFixedNpart()){
                    
                    std::cerr << "# RESAMPLING TO GET DESIRED Npart=" << param->getUseFixedNpart() << endl;
                    nucleusA.clear();
                    nucleusB.clear();
                    
                    xb = random->genrand64_real1(); // uniformly distributed random variable
                    
                    if(param->getLinearb()==1) // use a linear probability distribution for b if we are doing nuclei
                    {
                        b = sqrt((bmax*bmax-bmin*bmin)*xb+bmin*bmin);
                        std::cerr << "# SAMPLING bImpact LINEARLY BETWEEN " << bmin << " and " << bmax << " fm -- SET bImpact=" << b << " fm" << std::endl;
                        
                    }
                    else // use a uniform distribution instead
                    {
                        b = (bmax-bmin)*xb+bmin;
                        std::cerr << "# SAMPLING bImpact UNIFORMLY BETWEEN " << bmin << " and " << bmax << " fm -- SET bImpact=" << b << " fm" << std::endl;
                    }
                    
                    param->setb(b);
                    //std::cerr << "#SET bImpact=" << b << " fm" << endl;
                    
                    sampleTA(param, random, glauber);                           // populate the lists nucleusA and nucleusB with position data of the
                    setColorChargeDensity(lat, param, random, glauber,OutDirectory,RNGSEED,OUTPUT_FLAG);
                }
            }
            std::cerr << "# USING FIXED Npart=" << param->getNpart() << endl;
        }
        
        if(param->getSuccess()==0){
            
            std::cerr << "## NO COLLISIONS -- RESTARTING" << endl;
            return;
        }
        
    }
    
    // done.
    
    // -----------------------------------------------------------------------------
    // finish
    // -----------------------------------------------------------------------------
}

//void Init::multiplicity(GIPSLattice *lat, Group *group, Parameters *param, Random *random, Glauber *glauber)
//{
//    int N = param->getSize();
//    int pos;
//    double epsilonSum=0.;
//    double L = param->getL();
//    double a = L/N; // lattice spacing in fm
//
//    for(int ix=0; ix<N; ix++)
//    {
//        for(int iy=0; iy<N; iy++)
//        {
//            pos = ix*N+iy;
//            epsilonSum += a*a*lat->cells[pos]->getEpsilon()*hbarc;
//        }
//    }
//    stringstream strtE_name;
//    strtE_name << "totalEnergy.dat";
//    string tE_name;
//    tE_name = strtE_name.str();
//
//    ofstream fout(tE_name.c_str(),ios::out);
//    fout << epsilonSum << endl;
//    fout.close();
//}




