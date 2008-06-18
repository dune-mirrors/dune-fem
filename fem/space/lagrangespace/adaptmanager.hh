#ifndef DUNE_LAGRANGESPACE_ADAPTMANAGER_HH
#define DUNE_LAGRANGESPACE_ADAPTMANAGER_HH

#include <dune/grid/common/capabilities.hh>

#include <dune/fem/space/common/restrictprolonginterface.hh>

#include "lagrangespace.hh"

namespace Dune
{

  /** @ingroup RestrictProlongImpl
   *
   *  \brief Restriction / prolongation operator for Lagrange discrete
   *         function spaces
   */
  template< class DiscreteFunction, class FunctionSpace, class GridPart,
            int polOrder, template< class > class Storage >
  class RestrictProlongDefaultImplementation
    < DiscreteFunction,
      LagrangeDiscreteFunctionSpace
        < FunctionSpace, GridPart, polOrder, Storage > >
  : public RestrictProlongInterfaceDefault
    < RestrictProlongTraits< RestrictProlongDefaultImplementation
        < DiscreteFunction,
          LagrangeDiscreteFunctionSpace
            < FunctionSpace, GridPart, polOrder, Storage > > > >
  {
  public:
    //! type of the discrete function
    typedef DiscreteFunction DiscreteFunctionType;

    //! type of the discrete function space
    typedef LagrangeDiscreteFunctionSpace
      < FunctionSpace, GridPart, polOrder, Storage >
      DiscreteFunctionSpaceType;

  private:
    typedef RestrictProlongDefaultImplementation
      < DiscreteFunctionType, DiscreteFunctionSpaceType >
      ThisType;
    typedef RestrictProlongInterfaceDefault< RestrictProlongTraits< ThisType > >
      BaseType;

    using BaseType :: checkPersistent;
  public:
    //! field type of the discrete function's domain
    typedef typename DiscreteFunctionType :: DomainFieldType DomainFieldType;
    //! type of the discrete function's domain
    typedef typename DiscreteFunctionType :: DomainType DomainType;
    //! field type of the discrete function's range
    typedef typename DiscreteFunctionType :: RangeFieldType RangeFieldType;
    //! type of the discrete function's range
    typedef typename DiscreteFunctionType :: RangeType RangeType;
    //! type of the local functions
    typedef typename DiscreteFunctionType :: LocalFunctionType
      LocalFunctionType;

    //! type of the grid partition
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
    //! type of the grid
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;
    //! type of the Lagrange point set
    typedef typename DiscreteFunctionSpaceType :: LagrangePointSetType
      LagrangePointSetType;

    enum { dimRange = DiscreteFunctionSpaceType :: dimRange };

    typedef typename LagrangePointSetType
      :: template Codim< 0 > :: SubEntityIteratorType
      EntityDofIteratorType;

  private:
    DiscreteFunctionType &discreteFunction_;
    const DiscreteFunctionSpaceType &discreteFunctionSpace_;
    
  public:
    //! constructor
    inline explicit
    RestrictProlongDefaultImplementation
      ( DiscreteFunctionType &discreteFunction )
    : discreteFunction_( discreteFunction ),
      discreteFunctionSpace_( discreteFunction_.space() )
    {
      // make sure the index set can handle adaptivity
      assert( (Capabilities :: IsUnstructured< GridType > :: v) ?
             ( checkPersistent(discreteFunction_.space().indexSet()) ) : true );
    }

    //! if weight is set, it is assumed that the proportion between father's
    //! and son's volume is constant
    void setFatherChildWeight ( const RangeFieldType &weight ) const
    {
      // We make no use of this information.
    }

    //! restrict data to the father
    template< class EntityType >
    void restrictLocal ( EntityType &father, 
                         EntityType &son,
                         bool initialize ) const
    {
      // make sure the index set can handle adaptivity
      assert( checkPersistent(discreteFunction_.space().indexSet()) );

      typedef typename EntityType :: Geometry GeometryType;

      LocalFunctionType fatherFunction = discreteFunction_.localFunction( father );
      LocalFunctionType sonFunction = discreteFunction_.localFunction( son );

      const LagrangePointSetType &lagrangePointSet
        = discreteFunctionSpace_.lagrangePointSet( father );

      const GeometryType &geometryInFather = son.geometryInFather();

      EntityDofIteratorType it
        = lagrangePointSet.template beginSubEntity< 0 >( 0 );
      const EntityDofIteratorType endit
        = lagrangePointSet.template endSubEntity< 0 >( 0 );
      for( ; it != endit; ++it ) {
        const unsigned int dof = *it;
        const DomainType &pointInFather = lagrangePointSet.point( dof );
        const DomainType pointInSon = geometryInFather.local( pointInFather );
        if( geometryInFather.checkInside( pointInSon ) ) {
          RangeType phi;
          sonFunction.evaluate( pointInSon, phi );
          for( unsigned int coordinate = 0; coordinate < dimRange; ++coordinate )
            fatherFunction[ dimRange * dof + coordinate ] = phi[ coordinate ];
        }
      }
    }

    //! prolong data to children
    template< class EntityType >
    void prolongLocal ( EntityType &father, EntityType &son, bool initialize ) const
    {
      // make sure the index set can handle adaptivity
      assert( checkPersistent(discreteFunction_.space().indexSet()) );

      typedef typename EntityType :: Geometry GeometryType;

      LocalFunctionType fatherFunction = discreteFunction_.localFunction( father );
      LocalFunctionType sonFunction = discreteFunction_.localFunction( son );

      const LagrangePointSetType &lagrangePointSet
        = discreteFunctionSpace_.lagrangePointSet( son );

      const GeometryType &geometryInFather = son.geometryInFather();

      EntityDofIteratorType it
        = lagrangePointSet.template beginSubEntity< 0 >( 0 );
      const EntityDofIteratorType endit
        = lagrangePointSet.template endSubEntity< 0 >( 0 );
      for( ; it != endit; ++it ) {
        const unsigned int dof = *it;
        const DomainType &pointInSon = lagrangePointSet.point( dof );
        const DomainType pointInFather = geometryInFather.global( pointInSon );
        
        RangeType phi;
        fatherFunction.evaluate( pointInFather, phi );
        for( unsigned int coordinate = 0; coordinate < dimRange; ++coordinate )
          sonFunction[ dimRange * dof + coordinate ] = phi[ coordinate ];
      }
    }
  };
}
#endif
