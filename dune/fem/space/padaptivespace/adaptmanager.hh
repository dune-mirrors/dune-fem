#ifndef DUNE_PADATPTIVESPACE_ADAPTMANAGER_HH
#define DUNE_PADATPTIVESPACE_ADAPTMANAGER_HH

#include <dune/grid/common/capabilities.hh>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>
#include <dune/fem/space/common/localrestrictprolong.hh>
#include <dune/fem/space/dgspace/localrestrictprolong.hh>
#include <dune/fem/space/lagrangespace/lagrangespace.hh>
#include <dune/fem/space/padaptivespace/padaptivespace.hh>
#include <dune/fem/space/padaptivespace/restrictprolong.hh>

namespace Dune
{
  namespace Fem
  {

  // DefaultLocalRestrictProlong
  // ---------------------------

  template< class FS, class GP, int ord, template< class > class S >
  struct DefaultLocalRestrictProlong< Fem::PAdaptiveLagrangeSpace< FS, GP, ord, S > >
  : public PLagrangeLocalRestrictProlong< typename GP::GridType, Fem::PAdaptiveLagrangeSpace< FS, GP, ord, S > >
  {
    DefaultLocalRestrictProlong ( const Fem::PAdaptiveLagrangeSpace< FS, GP, ord, S > &space )
    : PLagrangeLocalRestrictProlong< typename GP::GridType, Fem::PAdaptiveLagrangeSpace< FS, GP, ord, S > >( space )
    {}
  };


  template< class FunctionSpaceImp, class GridPartImp, int polOrd, template< class > class StorageImp >
  struct DefaultLocalRestrictProlong< Fem::PAdaptiveDGSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp > >
  : public DiscontinuousGalerkinLocalRestrictProlong< Fem::PAdaptiveDGSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp >, false > // invert mass matrix or not 
  {
    typedef DiscontinuousGalerkinLocalRestrictProlong< Fem::PAdaptiveDGSpace<
                                                       FunctionSpaceImp, 
                                                       GridPartImp, 
                                                       polOrd, StorageImp >, false > BaseType ;
    DefaultLocalRestrictProlong ( const Fem::PAdaptiveDGSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp > & space ) 
      : BaseType( space )
    {}
  };

  template< class FunctionSpaceImp, class GridPartImp, template< class > class StorageImp >
  struct DefaultLocalRestrictProlong< Fem::PAdaptiveDGSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > >
  : public ConstantLocalRestrictProlong< Fem::PAdaptiveDGSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > >
  {
    DefaultLocalRestrictProlong ( const Fem::PAdaptiveDGSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > & )
    {}
  };



  // pAdaptation
  // -----------

  template <class DF, class Vector, class DFS> 
  void pAdaptation( DF& df, const Vector& polynomialOrders, const DFS &space, const int ) 
  {}

  /** \brief pAdaptation 
      \param df  discrete function to adapt 
      \param polynomialOrders  vector containing polynomial orders for each cell 
      \param space  type of space tp be adapted 
      \param polOrderShift possible shift of polynomial order (i.e. in case of
                           Taylor-Hood put -1 for the pressure) (default = 0)
  */
  template <class DF, class Vector,
            class FS, class GP, int p,
            template< class > class Storage >
  void pAdaptation( DF& df, 
                    const Vector& polynomialOrders, 
                    const Fem::PAdaptiveLagrangeSpace<FS,GP,p,Storage> &space,
                    const int polOrderShift = 0 ) 
  {
    /*
    typedef typename DF :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
    typedef typename GridPartType :: GridType  GridType;
    typedef typename DiscreteFunctionSpaceType :: IteratorType  IteratorType;

    DiscreteFunctionSpaceType& newSpace = const_cast< DiscreteFunctionSpaceType& > (df.space());

    DiscreteFunctionSpaceType oldSpace( df.space().gridPart() );

    typedef DofManager< GridType > DofManagerType;

    typedef typename IteratorType :: Entity  EntityType;

    DofManagerType& dm = DofManagerType :: instance( newSpace.grid() );

    const IteratorType endit = newSpace.end();
    for( IteratorType it = newSpace.begin(); it != endit; ++it ) 
    {
      const EntityType& entity = *it;
      oldSpace.blockMapper().setPolynomOrder( entity, newSpace.blockMapper().polynomOrder( entity ) ); 
    }

    dm.resize();
    dm.compress();

    AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > tmp( "padaptation", oldSpace );

    newSpace.adapt( tmp )
    tmp.assign( df );

    for( IteratorType it = newSpace.begin(); it != endit; ++it ) 
    {
      const EntityType& entity = *it;
      const int polOrder = polynomialOrders[ newSpace.indexSet().index( entity ) ] + polOrderShift ;
      newSpace.blockMapper().setPolynomOrder( entity, polOrder );
    }

    dm.resize();
    dm.compress();

    LagrangeInterpolation< DF > :: interpolateFunction( tmp, df );
    */
  }

  /** \brief pAdaptation 
      \param df  discrete function to adapt 
      \param polynomialOrders  vector containing polynomial orders for each cell 
      \param polOrderShift possible shift of polynomial order (i.e. in case of
        Taylor-Hood put -1 for the pressure) (default = 0)
  */
  template <class DF, class Vector> 
  void pAdaptation( DF& df, 
                    const Vector& polynomialOrders,
                    const int polOrderShift = 0 ) 
  {
    pAdaptation( df, polynomialOrders, df.space(), polOrderShift );
  }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_LAGRANGESPACE_ADAPTMANAGER_HH
