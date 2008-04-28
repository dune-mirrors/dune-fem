#ifndef DUNE_FEM_GRIDHELPER_HH
#define DUNE_FEM_GRIDHELPER_HH

#include <dune/grid/common/capabilities.hh>

#include <dune/fem/misc/metaprogramming.hh>

namespace Dune
{

  namespace Capabilities
  {

    template< class Grid >
    struct hasAllCodimEntities
    {
    private:
      template< unsigned int codim >
      struct Codim
      : public hasEntity< Grid, codim >
      {};
    
    public:
      static const bool v = Loop< MetaAnd, Codim, Grid :: dimension > :: v;
      static const bool value = v;
    };
    
  }
  
}

#endif
