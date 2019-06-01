#ifndef __FOURIER_SPACE_CPP__
#define __FOURIER_SPACE_CPP__

// INCLUDE TWO DIMENSIONAL FAST FOURIER TRANSFORM //
#include "../FFT/FFT2D.cpp"

namespace FourierSpace{

    FFT2D *RhoT;
    FFT2D *RhoP;
    
    FFT2D *OmegaS;
    FFT2D *OmegaA;
    
    FFT2D *Utarget;
    FFT2D *Scalar;

    // INITIALIZATION //
    void InitRhos(){
        
        RhoT=new FFT2D(TargetFields::U->N[0],TargetFields::U->N[1],SUNcAlgebra::VectorSize);
        RhoP=new FFT2D(TargetFields::U->N[0],TargetFields::U->N[1],SUNcAlgebra::VectorSize);
        
    }
    
    void InitOmegas(){
        
        OmegaS=new FFT2D(TargetFields::U->N[0],TargetFields::U->N[1],SUNcAlgebra::VectorSize);
        OmegaA=new FFT2D(TargetFields::U->N[0],TargetFields::U->N[1],SUNcAlgebra::VectorSize);
        
    }
    
    void InitUtarget(){
        
        Utarget=new FFT2D(TargetFields::U->N[0],TargetFields::U->N[1],SUNcGroup::MatrixSize);
        
    }
    void InitScalar(){

        Scalar=new FFT2D(TargetFields::U->N[0],TargetFields::U->N[1],1);
    }

    void InitAll(){
        
        RhoT=new FFT2D(TargetFields::U->N[0],TargetFields::U->N[1],SUNcAlgebra::VectorSize);
        RhoP=new FFT2D(TargetFields::U->N[0],TargetFields::U->N[1],SUNcAlgebra::VectorSize);
        
        OmegaS=new FFT2D(TargetFields::U->N[0],TargetFields::U->N[1],SUNcAlgebra::VectorSize);
        OmegaA=new FFT2D(TargetFields::U->N[0],TargetFields::U->N[1],SUNcAlgebra::VectorSize);
        
        Utarget=new FFT2D(TargetFields::U->N[0],TargetFields::U->N[1],SUNcGroup::MatrixSize);
        Scalar=new FFT2D(TargetFields::U->N[0],TargetFields::U->N[1],1);

    }
    
    // CLEAN-UP //
    void CleanUpRhos(){
        
        delete RhoT;
        delete RhoP;
        
    }
    
    void CleanUpOmegas(){
        
        delete OmegaS;
        delete OmegaA;
        
        
    }

    void CleanUpScalar(){
        
        delete Scalar;
        
    }
    
    void CleanUpUtarget(){
        
        delete Utarget;
        
    }
    
    void CleanUpAll(){
        
        delete RhoT;
        delete RhoP;
        
        delete OmegaS;
        delete OmegaA;

        delete Utarget;
        delete Scalar;
    
    }
    
}

#endif
