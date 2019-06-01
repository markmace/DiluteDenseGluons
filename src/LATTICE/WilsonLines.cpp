#ifndef __WILSONLINES__CPP__
#define __WILSONLINES__CPP__

class WilsonLines{
    
public:
    
    INT Area;
    
    const static INT Dimension=2;
    
    DOUBLE a[Dimension];
    
    DOUBLE aSquare;
    
    INT N[Dimension];
    
    SU_Nc_FUNDAMENTAL_FORMAT *U;
    
    INT Index2D(INT x,INT y){
        
        return (MOD(x,N[0])+N[0]*MOD(y,N[1]));
        
    }
    
    INT Index(INT x,INT y,INT mu){
        
        return SUNcGroup::MatrixSize*(mu+Dimension*Index2D(x,y));
    }
    
    void GetPosition(INT Index,INT &x,INT &y){
        
        y=Index/(N[0]*N[1]); x=(Index-y*(N[0]*N[1]))/(N[0]);
        
    }
    
    SU_Nc_FUNDAMENTAL_FORMAT* Get(INT x,INT y,INT mu){
        
        return &(U[Index(x,y,mu)]);
    }
    
    void SetIdentity(){
        
        for(INT y=0;y<=N[1]-1;y++){
            for(INT x=0;x<=N[0]-1;x++){
                for(INT mu=0;mu<Dimension;mu++){
                    
                    COPY_SUNcMatrix(this->Get(x,y,mu),SUNcGroup::UnitMatrix);
                    
                }
            }
        }
    }
    
    
    WilsonLines(INT Nx,INT Ny,DOUBLE ax,DOUBLE ay){
        
        this->Area=Nx*Ny;
        
        this->N[0]=Nx;
        this->N[1]=Ny;
        
        this->a[0]=ax;
        this->a[1]=ay;
        
        this->aSquare=a[0]*a[1];
        
        U=new SU_Nc_FUNDAMENTAL_FORMAT[Dimension*SUNcGroup::MatrixSize*Area];
        
    }
    
    ~WilsonLines(){
        
        delete U;
    }
    
    
};

void Copy(WilsonLines *UCopy,WilsonLines *UOriginal){
    
    for(INT i=0;i<UOriginal->Dimension;i++){
        
        UCopy->a[i]=UOriginal->a[i];
        
        UCopy->N[i]=UOriginal->N[i];
        
    }
    
    UCopy->Area=UOriginal->Area;
    
    UCopy->aSquare=UOriginal->aSquare;
    
    std::memcpy(UCopy->Get(0,0,0),UOriginal->Get(0,0,0),UOriginal->Dimension*SUNcGroup::MatrixSize*UOriginal->Area*sizeof(SU_Nc_FUNDAMENTAL_FORMAT));
    
}
#endif

