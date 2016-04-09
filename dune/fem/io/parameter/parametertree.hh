#ifndef DUNE_FEM_IO_PARAMETER_PARAMETERTREE_HH
#define DUNE_FEM_IO_PARAMETER_PARAMETERTREE_HH

#include <dune/common/parametertree.hh>

#include <dune/fem/io/parameter/reader.hh>

namespace Dune
{

  namespace Fem
  {

    // parameterReader
    // ---------------

    inline static ParameterReader parameterReader ( const ParameterTree &parameterTree )
    {
      return ParameterReader( [ &parameterTree ] ( const std::string &key, const std::string *defaultValue ) {
            return (parameterTree.hasKey( key ) ? &parameterTree[ key ] : defaultValue);
          } );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_IO_PARAMETER_PARAMETERTREE_HH
