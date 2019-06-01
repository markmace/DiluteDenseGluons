#ifndef __VECTORFIELDS__CPP__
#define __VECTORFIELDS__CPP__

class VectorFields{
    
public:
    
    INT Area;
    
    const static INT Dimension=2;
    
    DOUBLE a[Dimension];
    
    DOUBLE aSquare;
    
    INT N[Dimension];
    
    SU_Nc_ALGEBRA_FORMAT *E;
    
    INT Index2D(INT x,INT y){
        
        return (MOD(x,N[0])+N[0]*MOD(y,N[1]));
        
    }
    
    INT Index(INT x,INT y,INT mu,INT a){
        
        return a+SUNcAlgebra::VectorSize*(mu+Dimension*Index2D(x,y));
        
    }
    
    SU_Nc_ALGEBRA_FORMAT* Get(INT x,INT y,INT mu,INT a){
        
        return &E[Index(x,y,mu,a)];
    }
    
    void SetZero(){
        
        for(INT y=0;y<=N[1]-1;y++){
            for(INT x=0;x<=N[0]-1;x++){
                for(int mu=0;mu<Dimension;mu++){
                    for(int a=0;a<SUNcAlgebra::VectorSize;a++){
                        
                        this->Get(x,y,mu,a)[0]=DOUBLE(0.0);
                        
                    }
                }
            }
        }
    }
    
    
    VectorFields(INT Nx,INT Ny,DOUBLE ax,DOUBLE ay){
        
        this->Area=Nx*Ny;
        
        this->N[0]=Nx;
        this->N[1]=Ny;
        
        this->a[0]=ax;
        this->a[1]=ay;
        
        this->aSquare=a[0]*a[1];
        
        E=new SU_Nc_ALGEBRA_FORMAT[Dimension*SUNcAlgebra::VectorSize*Area];
        
    }
    
    ~VectorFields(){
        
        delete E;
    }
};

void Copy(VectorFields *ECopy,VectorFields *EOriginal){
    
    for(INT i=0;i<EOriginal->Dimension;i++){
        
        ECopy->a[i]=EOriginal->a[i];
        
        ECopy->N[i]=EOriginal->N[i];
        
    }
    
    ECopy->Area=EOriginal->Area;
    
    ECopy->aSquare=EOriginal->aSquare;
    
    std::memcpy(ECopy->Get(0,0,0,0),EOriginal->Get(0,0,0,0),Lattice::Dimension*SUNcAlgebra::VectorSize*Lattice::Area*sizeof(SU_Nc_ALGEBRA_FORMAT));
    
    
}

#endif

