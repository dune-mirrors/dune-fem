#ifndef DUNE_PADATPTIVESPACE_ADAPTMANAGER_HH
#define DUNE_PADATPTIVESPACE_ADAPTMANAGER_HH

#include <dune/grid/common/capabilities.hh>

#include <dune/fem/space/common/restrictprolonginterface.hh>
#include <dune/fem/space/padaptivespace/restrictprolong.hh>

#include <dune/fem/space/lagrangespace/lagrangespace.hh>
#include <dune/fem/space/padaptivespace/padaptivespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>

namespace Dune
{

  /** @ingroup RestrictProlongImpl
   *
   *  \brief Restriction / prolongation operator for Lagrange discrete
   *         function spaces
   */
  template< class DF, class FS, class GP, int ord, template< class > class S >
  class RestrictProlongDefaultImplementation< DF, PAdaptiveLagrangeSpace< FS, GP, ord, S > >
  : public RestrictProlongInterfaceDefault< RestrictProlongTraits
      < RestrictProlongDefaultImplementation< DF, PAdaptiveLagrangeSpace< FS, GP, ord, S > >,
        typename GP::GridType::ctype > >
  {
    typedef RestrictProlongDefaultImplementation
      < DF, PAdaptiveLagrangeSpace< FS, GP, ord, S > >
      ThisType;
    typedef RestrictProlongInterfaceDefault< RestrictProlongTraits< ThisType, typename GP::GridType::ctype > >
      BaseType;

    typedef typename GP::GridType::ctype DomainFieldType;

  public:
    //! type of the discrete function
    typedef DF DiscreteFunctionType;

    //! type of the discrete function space
    typedef PAdaptiveLagrangeSpace< FS, GP, ord, S > DiscreteFunctionSpaceType;

  protected:
    using BaseType::entitiesAreCopies;

  public:
    //! type of the local functions
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

    //! type of the grid
    typedef typename DiscreteFunctionSpaceType::GridType Grid;

  public:
    //! constructor
    explicit
    RestrictProlongDefaultImplementation ( DiscreteFunctionType &discreteFunction )
    : discreteFunction_( discreteFunction ),
      localRestrictProlong_( discreteFunction.space() )
    {}

    //! restrict data to the father
    template< class Entity >
    void restrictLocal ( const Entity &father, const Entity &son, bool initialize ) const
    {
      if( !entitiesAreCopies( discreteFunction_.space().indexSet(), father, son ) )
      {
        LocalFunctionType fatherFunction = discreteFunction_.localFunction( father );
        LocalFunctionType sonFunction = discreteFunction_.localFunction( son );

        localRestrictProlong_.restrictLocal( fatherFunction, sonFunction, initialize );
      }
    }

    //! prolong data to children
    template< class Entity >
    void prolongLocal ( const Entity &father, const Entity &son, bool initialize ) const
    {
      if( !entitiesAreCopies( discreteFunction_.space().indexSet(), father, son ) )
      {
        LocalFunctionType fatherFunction = discreteFunction_.localFunction( father );
        LocalFunctionType sonFunction = discreteFunction_.localFunction( son );

        localRestrictProlong_.prolongLocal( fatherFunction, sonFunction );
      }
    }

    //! add discrete function to communicator 
    template< class Communicator >
    void addToList ( Communicator &comm )
    {
      // for Lagrange spaces this communication is not needed (since
      // data on ghosts is neglected 

      //  comm.addToList( discreteFunction_ );
    }

  private:
    DiscreteFunctionType &discreteFunction_;
    PLagrangeLocalRestrictProlong< Grid, DiscreteFunctionSpaceType > localRestrictProlong_;
  };

  template <class DF, class Vector, class DFS> 
  void pAdaptation( DF& df, const Vector& polynomialOrders, const DFS &space, const int ) 
  {
  }

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
                    const PAdaptiveLagrangeSpace<FS,GP,p,Storage> &space,
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

}

#endif // #ifndef DUNE_LAGRANGESPACE_ADAPTMANAGER_HH
