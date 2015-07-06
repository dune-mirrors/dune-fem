// vim: set expandtab ts=2 sw=2 sts=2:
#ifndef DUNE_FEM_SOLVERPARAMETER_HH
#define DUNE_FEM_SOLVERPARAMETER_HH

#include <dune/fem/io/parameter.hh>

namespace Dune
{

  namespace Fem
  {


    struct SolverParameter
#ifndef DOXYGEN
    : public LocalParameter< SolverParameter, SolverParameter >
#endif
    {
      protected:
      // key prefix, default is fem.solver (can be overloaded by user)
      const std::string keyPrefix_;

      public:

      SolverParameter ( const std::string keyPrefix = "fem.solver." )
        : keyPrefix_( keyPrefix )
      {}

      virtual bool verbose() const
      {
        return Parameter::getValue< bool >( keyPrefix_ + "verbose", false );
      }

      virtual int errorMeasure() const
      {
        static const std::string errorTypeTable[] = { "absolute", "relative" };
        return Parameter::getEnum( keyPrefix_ + "errormeasure", errorTypeTable, 0 );
      }

    };

  }
}

#endif // #ifndef DUNE_FEM_SOLVERPARAMETER_HH
