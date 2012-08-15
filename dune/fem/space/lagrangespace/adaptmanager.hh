#ifndef DUNE_FEM_LAGRANGESPACE_ADAPTMANAGER_HH
#define DUNE_FEM_LAGRANGESPACE_ADAPTMANAGER_HH

#include <dune/fem/space/common/localrestrictprolong.hh>
#include <dune/fem/space/lagrangespace/lagrangespace.hh>
#include <dune/fem/space/lagrangespace/restrictprolong.hh>

namespace Dune
{

  namespace Fem 
  {

    template< class FS, class GP, int ord, template< class > class S >
    struct DefaultLocalRestrictProlong< LagrangeDiscreteFunctionSpace< FS, GP, ord, S > >
    : public LagrangeLocalRestrictProlong< typename GP::GridType, ord >
    {
      DefaultLocalRestrictProlong ( const LagrangeDiscreteFunctionSpace< FS, GP, ord, S > & )
      {}
    };

  } // namespace Fem 

} // namespace Dune

#endif // #ifndef DUNE_FEM_LAGRANGESPACE_ADAPTMANAGER_HH
