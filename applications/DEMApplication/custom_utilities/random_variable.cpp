#include "random_variable.h"

namespace Kratos {
    RandomVariable::RandomVariable(){
    }
    RandomVariable::RandomVariable(const Parameters rParameters){
    }

    std::string RandomVariable::Info() const
    {
        std::stringstream buffer;
        buffer << "RandomVariable" ;
        return buffer.str();
    }
    void RandomVariable::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Abstract RandomVariable";
    }
    void RandomVariable::PrintData(std::ostream& rOStream) const {}
}


