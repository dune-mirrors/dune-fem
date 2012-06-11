#ifndef DUNE_COMBINEDADAPTMANAGERIMP_HH
#define DUNE_COMBINEDADAPTMANAGERIMP_HH

//- local includes  
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/common/restrictprolonginterface.hh>

//- local includes 
#include "combinedspace.hh"


namespace Dune
{

namespace Fem 
{

/** @ingroup RestrictProlongImpl
    @{
**/

#ifdef USE_OLD_COMBINEDSPACE
  /** \brief specialization of RestrictProlongDefault for
      CombinedSpace.
  */
  template < class DiscreteFunctionSpaceImp, 
             int N, 
             DofStoragePolicy policy> 
  struct DefaultLocalRestrictProlong< CombinedSpace<DiscreteFunctionSpaceImp,N,policy> > 
  : public DiscontinuousGalerkinLocalRestrictProlong< CombinedSpace<DiscreteFunctionSpaceImp,N,policy> >
  {
    DefaultLocalRestrictProlong( const CombinedSpace<DiscreteFunctionSpaceImp,N,policy> & space)
    {
      assert( !space.continuous() );
    }
  };
#else 
  template < class DiscreteFunctionSpaceImp, 
             int N, 
             DofStoragePolicy policy> 
  struct DefaultLocalRestrictProlong< CombinedSpace<DiscreteFunctionSpaceImp,N,policy> > 
  : public EmptyLocalRestrictProlong< CombinedSpace<DiscreteFunctionSpaceImp,N,policy> >
  {
    DefaultLocalRestrictProlong( const CombinedSpace<DiscreteFunctionSpaceImp,N,policy> & )
    {}
  };
#endif 

///@}

} // end namespace Fem  

} // end namespace Dune 
#endif
