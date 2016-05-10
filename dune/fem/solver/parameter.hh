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

      ParameterReader parameter_;

      public:

      explicit SolverParameter ( const ParameterReader &parameter = Parameter::container() )
        : keyPrefix_( "fem.solver." ), parameter_( parameter )
      {}

      explicit SolverParameter ( const std::string keyPrefix, const ParameterReader &parameter = Parameter::container() )
        : keyPrefix_( keyPrefix ), parameter_( parameter )
      {}

      virtual bool verbose() const
      {
        return Parameter::getValue< bool >( keyPrefix_ + "verbose", false );
      }

      virtual int errorMeasure() const
      {
        static const std::string errorTypeTable[] = { "absolute", "relative" };
        return parameter_.getEnum( keyPrefix_ + "errormeasure", errorTypeTable, 0 );
      }

    };

  }
}

#endif // #ifndef DUNE_FEM_SOLVERPARAMETER_HH
