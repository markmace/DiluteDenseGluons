// OutputManagement.cpp is part of the Dilute-Dense Gluon solver //
// Copyright (C) 2019 Mark Mace //

namespace  IO {
    
    std::string OutputDirectory;
    
    template<typename GenericArgument>
    
    void SetOutputDirectory(GenericArgument x){
        
        OutputDirectory=StringManipulation::StringCast(x,"/");
        
        //std::cer << "# OUTPUT DIRECTORY IS " << OutputDirectory <<  std::endl;
    }

}
