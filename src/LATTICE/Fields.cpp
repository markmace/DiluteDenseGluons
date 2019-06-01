#ifndef __GLINKSEFIELDS__CPP__
#define __GLINKSEFIELDS__CPP__

namespace TargetFields{
    
    WilsonLines *U;
        
    void Init(){
        
        U=new WilsonLines(Lattice::N[0],Lattice::N[1],Lattice::a[0],Lattice::a[1]);
        
    }
    
}

namespace TempTargetFields{
    
    WilsonLines *U;
    
    void Init(){
        
        U=new WilsonLines(Lattice::N[0],Lattice::N[1],Lattice::a[0],Lattice::a[1]);
        
    }
    
}

namespace Tpp{
    
    VectorFields *Proj;
    VectorFields *Targ;

    
    void Init(){
        
        Proj=new VectorFields(Lattice::N[0],Lattice::N[1],Lattice::a[0],Lattice::a[1]);
        Targ=new VectorFields(Lattice::N[0],Lattice::N[1],Lattice::a[0],Lattice::a[1]);
        
    }
    
}

namespace g2mu{
    
    VectorFields *Proj;
    VectorFields *Targ;
    
    
    void Init(){
        
        Proj=new VectorFields(Lattice::N[0],Lattice::N[1],Lattice::a[0],Lattice::a[1]);
        Targ=new VectorFields(Lattice::N[0],Lattice::N[1],Lattice::a[0],Lattice::a[1]);
        
    }
    
}

namespace OmegaS{
    
    CVectorFields *O;
    
    void Init(){
        
        O=new CVectorFields(Lattice::N[0],Lattice::N[1],Lattice::a[0],Lattice::a[1]);
        
    }
    
}

namespace OmegaA{
    
    CVectorFields *O;
    
    void Init(){
        
        O=new CVectorFields(Lattice::N[0],Lattice::N[1],Lattice::a[0],Lattice::a[1]);
        
    }
    
}

namespace ProjSolution{
    
    VectorFields *A;
    
    void Init(){
        
        A=new VectorFields(Lattice::N[0],Lattice::N[1],Lattice::a[0],Lattice::a[1]);
        
    }
    
}

namespace Lattice{
    
    ///////////////////////////
    //INITIALIZATION ROUTINE //
    ///////////////////////////
    
    void InitTargetAndProjFields(){
        
        TargetFields::Init();
        TempTargetFields::Init();
        
        ProjSolution::Init();
        
        Tpp::Init();
        
        g2mu::Init();
        
    }
    
    void InitOmegas(){
        
        OmegaS::Init();
        OmegaA::Init();
    }
    
    void InitTpp(){
        
        Tpp::Init();
        
    }
    
    void Initg2mu(){
        
        g2mu::Init();
        
    }

    
    void InitAll(){
        
        TargetFields::Init();
        TempTargetFields::Init();
        
        OmegaS::Init();
        OmegaA::Init();

        ProjSolution::Init();
        
        Tpp::Init();
        
        g2mu::Init();

    }
    
    // CLEAN-UP //
    void CleanUpTargetAndProjFields(){
        
        delete TargetFields::U;
        delete TempTargetFields::U;
        
        delete ProjSolution::A;

    }
    
    void CleanUpOmegas(){
        
        delete OmegaS::O;
        delete OmegaA::O;
                
    }
    
    void CleanUpg2mu(){

        delete g2mu::Proj;
        delete g2mu::Targ;
        
    }
    
    void CleanUpTpp(){
        
        delete Tpp::Proj;
        delete Tpp::Targ;
        
    }
    
    void CleanUpAll(){
        
        delete TargetFields::U;
        delete TempTargetFields::U;
        
        delete OmegaS::O;
        delete OmegaA::O;
        
        delete ProjSolution::A;
        
        delete Tpp::Proj;
        delete Tpp::Targ;
        
        delete g2mu::Proj;
        delete g2mu::Targ;
        
    }
    
}
#endif
