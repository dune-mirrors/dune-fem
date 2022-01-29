#ifndef DUNE_FEM_SPACE_PADAPTIVE_LAGRANGE_HH
#define DUNE_FEM_SPACE_PADAPTIVE_LAGRANGE_HH

#include <dune/fem/common/hybrid.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/lagrange/shapefunctionset.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/space/lagrange/interpolation.hh>

#include "adaptmanager.hh"
#include "declaration.hh"
#include "generic.hh"
#include "mapper.hh"
#include "restrictprolong.hh"


namespace Dune
{

  namespace Fem
  {

    /** \addtogroup PAdaptiveLagrangeSpace
     *
     *  Provides access to base function sets for different element types in
     *  one grid and size of function space and maps from local to global dof
     *  number.
     *
     *  \note This space can only be used with special index sets. If you want
     *  to use the PAdaptiveLagrangeSpace with an index set only
     *  supporting the index set interface you will have to use the
     *  IndexSetWrapper class to provide the required functionality.
     *
     *  \note For adaptive calculations one has to use index sets that are
     *  capable of adaption (i.e. the method adaptive returns true). See also
     *  AdaptiveLeafIndexSet.
     */

    // PAdaptiveLagrangeSpaceTraits
    // ----------------------------

    template< class FunctionSpace, class GridPart, int maxPolOrder, class Storage >
    struct PAdaptiveLagrangeSpaceTraits
    {
      static_assert((maxPolOrder > 0), "LagrangeSpace only defined for maxPolOrder > 0" );

      typedef PAdaptiveLagrangeSpace< FunctionSpace, GridPart, maxPolOrder, Storage > DiscreteFunctionSpaceType;

      typedef FunctionSpace FunctionSpaceType;
      typedef GridPart GridPartType;

      static const int polynomialOrder = maxPolOrder;

      static const bool continuousSpace = true ;
      typedef Hybrid::IndexRange< int, FunctionSpaceType::dimRange > LocalBlockIndices;

      typedef PAdaptiveLagrangeMapper< GridPartType, polynomialOrder > BlockMapperType;

      typedef LagrangePointSet< GridPartType, polynomialOrder > CompiledLocalKeyType;

      static const int codimension = 0;

    private:
      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;

      static const int dimLocal = GridPartType::dimension;
      typedef typename FunctionSpace::ScalarFunctionSpaceType ScalarFunctionSpaceType;
      typedef typename ToNewDimDomainFunctionSpace< ScalarFunctionSpaceType, dimLocal >::Type ShapeFunctionSpaceType;

      typedef LagrangeShapeFunctionInterface< ShapeFunctionSpaceType > ShapeFunctionType;
      typedef SimpleShapeFunctionSet< ShapeFunctionType > SimpleShapeFunctionSetType;

    public:
      typedef SelectCachingShapeFunctionSet< SimpleShapeFunctionSetType, Storage > ScalarShapeFunctionSetType;

      template< int pOrd >
      struct ScalarShapeFunctionSetFactory
      {
        struct Type
        {
          static ScalarShapeFunctionSetType *createObject ( const GeometryType &type )
          {
            typedef LagrangeShapeFunctionFactory< ShapeFunctionSpaceType, maxPolOrder > SimpleShapeFunctionSetFactoryType;
            return new ScalarShapeFunctionSetType( type, SimpleShapeFunctionSetType( SimpleShapeFunctionSetFactoryType( type, pOrd ) ) );
          }

          static void deleteObject ( ScalarShapeFunctionSetType *object ) { delete object; }
        };
      };

      typedef ShapeFunctionSetProxy< ScalarShapeFunctionSetType > ScalarShapeFunctionSetProxyType;
      typedef VectorialShapeFunctionSet< ScalarShapeFunctionSetProxyType, typename FunctionSpaceType::RangeType > ShapeFunctionSetType;

      typedef Dune::Fem::DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType > BasisFunctionSetType;

      template< class DiscreteFunction, class Operation = DFCommunicationOperation::Add >
      struct CommDataHandle
      {
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
        typedef Operation OperationType;
      };
    };



    // PAdaptiveLagrangeSpace
    // ----------------------

    /** \class   PAdaptiveLagrangeSpace
     *
     *  \ingroup PAdaptiveLagrangeSpace
     *
     *  \brief   Lagrange discrete function space
     */
    template< class FunctionSpace, class GridPart, int maxPolOrder, class Storage = CachingStorage >
    class PAdaptiveLagrangeSpace
    : public GenericDiscreteFunctionSpace< PAdaptiveLagrangeSpaceTraits< FunctionSpace, GridPart, maxPolOrder, Storage > >
    {
      typedef PAdaptiveLagrangeSpace< FunctionSpace, GridPart, maxPolOrder, Storage > ThisType;
      typedef GenericDiscreteFunctionSpace< PAdaptiveLagrangeSpaceTraits< FunctionSpace, GridPart, maxPolOrder, Storage > > BaseType;

    public:
      typedef ThisType PAdaptiveLagrangeSpaceType;

      typedef typename BaseType::Traits Traits;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::IntersectionType IntersectionType;
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;
      typedef typename BaseType::EntityType EntityType;

      typedef typename BaseType::CompiledLocalKeyType CompiledLocalKeyType;
      typedef CompiledLocalKeyType LagrangePointSetType;

      typedef LagrangeLocalInterpolation< GridPartType, maxPolOrder, BasisFunctionSetType > InterpolationImplType;
      typedef LocalInterpolationWrapper< ThisType > InterpolationType;

    public:
      using BaseType::blockMapper;
      using BaseType::compiledLocalKey;
      using BaseType::continuous;
      using BaseType::gridPart;
      using BaseType::order;
      using BaseType::basisFunctionSet;

      // default communication interface
      static const InterfaceType defaultInterface = GridPart::indexSetInterfaceType;

      // default communication direction
      static const CommunicationDirection defaultDirection = ForwardCommunication;

      /** \brief constructor
       *
       *  \param[in]  gridPart       grid part for the Lagrange space
       *  \param[in]  order          dynamically set maximal polynomial order between 1 and maxPolOrder
       *  \param[in]  commInterface  communication interface to use (optional)
       *  \param[in]  commDirection  communication direction to use (optional)
       */
      explicit PAdaptiveLagrangeSpace ( GridPartType &gridPart,
                                        const int order,
                                        const InterfaceType commInterface = defaultInterface,
                                        const CommunicationDirection commDirection = defaultDirection )
      : BaseType( gridPart, order, commInterface, commDirection )
      {}

      /** \brief constructor
       *
       *  \param[in]  gridPart       grid part for the Lagrange space
       *  \param[in]  commInterface  communication interface to use (optional)
       *  \param[in]  commDirection  communication direction to use (optional)
       */
      explicit PAdaptiveLagrangeSpace ( GridPartType &gridPart,
                                        const InterfaceType commInterface = defaultInterface,
                                        const CommunicationDirection commDirection = defaultDirection )
      : BaseType( gridPart, maxPolOrder, commInterface, commDirection )
      {}

      // copy constructor needed for p-adaption
      PAdaptiveLagrangeSpace ( const PAdaptiveLagrangeSpace &other )
      : BaseType( other )
      {}

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous (const IntersectionType &intersection) const
      {
        if ( order() > 0 && intersection.conforming())
        {
          return true;
          if (intersection.neighbor())
            return (order(intersection.inside()) == order(intersection.outside()));
          else
            return true;
        }
        return false;
      }

      /** \brief provide access to the Lagrange point set for an entity
       *
       *  \note This method is not part of the DiscreteFunctionSpaceInterface. It
       *        is unique to the LagrangeDiscreteFunctionSpace.
       *
       *  \param[in]  entity  entity the Lagrange point set is requested for
       *
       *  \returns LagrangePointSet
       */
      template< class EntityType >
      const CompiledLocalKeyType &lagrangePointSet ( const EntityType &entity ) const
      {
        return compiledLocalKey( entity.type(), blockMapper().polynomOrder( entity ) );
      }

      InterpolationType interpolation () const
      {
        return InterpolationType( *this );
      }

      [[deprecated]]
      InterpolationImplType interpolation ( const EntityType &entity ) const
      {
        return localInterpolation( entity );
      }

      InterpolationImplType localInterpolation ( const EntityType &entity ) const
      {
        return InterpolationImplType( lagrangePointSet( entity ), basisFunctionSet( entity ) );
      }

    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_PADAPTIVE_LAGRANGE_HH
