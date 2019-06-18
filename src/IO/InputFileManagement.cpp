// InputManagement.cpp is part of the Dilute-Dense Gluon solver //
// Copyright (C) 2019 Mark Mace //

namespace  IO {
    
    std::string InputFile;
    
    template<typename GenericArgument>
    
    void SetInputFile(GenericArgument x){
        
        InputFile=StringManipulation::StringCast(x);
        
        //std::cerr << "# INPUT FILE IS " << InputFile <<  std::endl;
    }
    
}
