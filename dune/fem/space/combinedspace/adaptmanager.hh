#ifndef DUNE_FEM_COMBINEDADAPTMANAGER_HH
#define DUNE_FEM_COMBINEDADAPTMANAGER_HH

#include <dune/common/exceptions.hh>
//- local includes  
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>


namespace Dune
{

  namespace Fem
  {

    template< class , class >
    class CombinedDiscreteFunctionSpace; 

    template< class SP1, class SP2 >
    struct DefaultLocalRestrictProlong< CombinedDiscreteFunctionSpace< SP1, SP2 > >
    : public EmptyLocalRestrictProlong< CombinedDiscreteFunctionSpace< SP1, SP2 > >
    {
      DefaultLocalRestrictProlong ( const CombinedDiscreteFunctionSpace< SP1, SP2 > & )
      {}
    };

  } // namespace Fem

} // namespace Dune 

#endif //  #ifndef DUNE_FEM_COMBINEDADAPTMANAGER_HH
