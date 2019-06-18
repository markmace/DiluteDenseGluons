#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace RandomNumberGenerator{
    
    //GSL RANDOM NUMBER GENERATORS AND SEED
    const gsl_rng_type *Type;
    gsl_rng *Generator;
    long int MySEED;
    
    //INITIALIZATION OF RANDOM NUMBER GENERATOR
    void Init(long int SEED){
        
        gsl_rng_env_setup();
        
        //Type=gsl_rng_default;
        Type=gsl_rng_mt19937;
        Generator=gsl_rng_alloc(Type);
        
        MySEED=SEED;
        
        gsl_rng_set(Generator,SEED);
        
    }
    
    // UNIFORM DISTRIBUTED RANDOM NUMBER
    DOUBLE rng(){
        return gsl_rng_uniform(Generator);
    }
    
    //UNIFORM DISTRIBUTED RANDOM NUMBER
    DOUBLE RandomInt(INT MaxNum){
        return gsl_rng_uniform_int(Generator,MaxNum);
    }
    
    //GAUSSIAN DISTRIBUTED RANDOM NUMBER
    DOUBLE Gauss(){
        return gsl_ran_gaussian(Generator,1.0);
    }
    
    //GAUSSIAN DISTRIBUTED RANDOM NUMBER
    DOUBLE Gauss(DOUBLE Amplitude){
        
        return gsl_ran_gaussian(Generator,Amplitude);
    }
    
    //GAUSSIAN DISTRIBUTED RANDOM NUMBER
    std::complex<DOUBLE> ComplexGauss(){
        
        DOUBLE Re=gsl_ran_gaussian(Generator,1.0); DOUBLE Im=gsl_ran_gaussian(Generator,1.0);
        
        return std::complex<DOUBLE>(Re,Im)/D_SQRT2;
    }
    
    
    //SU(Nc) RANDOM MATRIX
    void SUNcMatrix(DOUBLE Amplitude,COMPLEX *V){
        
        SU_Nc_ALGEBRA_FORMAT alpha[SUNcAlgebra::VectorSize];
        
        for(int a=0;a<SUNcAlgebra::VectorSize;a++){
            
            alpha[a]=rng();
        }
        
        SUNcAlgebra::Operations::MatrixIExp(Amplitude,alpha,V);
        
        
    }
    
    void CleanUp(){
    
        gsl_rng_free(Generator);
        
    }
    
    
}
