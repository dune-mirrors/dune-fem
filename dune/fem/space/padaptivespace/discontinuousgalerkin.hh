#ifndef DUNE_FEM_SPACE_PADAPTIVE_DISCONTINUOUSGALERKIN_HH
#define DUNE_FEM_SPACE_PADAPTIVE_DISCONTINUOUSGALERKIN_HH

#include <dune/fem/common/hybrid.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>

#include "adaptmanager.hh"
#include "declaration.hh"
#include "generic.hh"
#include "lagrange.hh"
#include "mapper.hh"
#include "restrictprolong.hh"

#include <dune/fem/space/discontinuousgalerkin/interpolation.hh>

namespace Dune
{

  namespace Fem
  {

    /** \addtogroup PAdaptiveDGSpace
     *
     *  Provides access to base function sets for different element types in
     *  one grid and size of function space and maps from local to global dof
     *  number.
     *
     *  \note This space can only be used with special index sets. If you want
     *  to use the PAdaptiveDGSpace with an index set only
     *  supporting the index set interface you will have to use the
     *  IndexSetWrapper class to provide the required functionality.
     *
     *  \note For adaptive calculations one has to use index sets that are
     *  capable of adaption (i.e. the method adaptive returns true). See also
     *  AdaptiveLeafIndexSet.
     */

    // PAdaptiveDGSpaceTraits
    // ----------------------

    template< class FunctionSpace, class GridPart, int polOrder, class Storage >
    struct PAdaptiveDGSpaceTraits
    : public PAdaptiveLagrangeSpaceTraits< FunctionSpace, GridPart, polOrder, Storage >
    {
      typedef PAdaptiveDGSpace< FunctionSpace, GridPart, polOrder, Storage > DiscreteFunctionSpaceType;

      static const bool continuousSpace = false ;
      typedef Hybrid::IndexRange< int, FunctionSpace::dimRange > LocalBlockIndices;

      typedef PAdaptiveDGMapper< GridPart, polOrder > BlockMapperType;

      template< class DiscreteFunction,
                class Operation = DFCommunicationOperation :: Copy >
      struct CommDataHandle
      {
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
        typedef Operation OperationType;
      };
    };



    // PAdaptiveDGSpace
    // ----------------

    /** \class   PAdaptiveDGSpace
     *
     *  \ingroup PAdaptiveDGSpace
     *
     *  \brief   adaptive DG discrete function space
     */
    template< class FunctionSpace, class GridPart, int polOrder, class Storage = CachingStorage >
    class PAdaptiveDGSpace
    : public GenericDiscreteFunctionSpace< PAdaptiveDGSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > >
    {
      typedef PAdaptiveDGSpace< FunctionSpace, GridPart, polOrder, Storage > ThisType;
      typedef GenericDiscreteFunctionSpace< PAdaptiveDGSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > > BaseType;

    public:
      typedef ThisType PAdaptiveDGSpaceType;

      typedef typename BaseType::Traits Traits;

      typedef typename BaseType::GridPartType      GridPartType;
      typedef typename BaseType::IntersectionType  IntersectionType;
      typedef typename BaseType::EntityType        EntityType;

      typedef typename BaseType::BasisFunctionSetType  BasisFunctionSetType;
      typedef typename BaseType::CompiledLocalKeyType CompiledLocalKeyType;
      typedef CompiledLocalKeyType LagrangePointSetType;

      typedef DiscontinuousGalerkinLocalL2Projection< GridPartType, BasisFunctionSetType > InterpolationImplType;
      typedef LocalInterpolationWrapper< ThisType > InterpolationType;

    public:
      using BaseType::continuous;
      using BaseType::gridPart;
      using BaseType::blockMapper;
      using BaseType::compiledLocalKey;
      using BaseType::basisFunctionSet;

      // default communication interface
      static const InterfaceType defaultInterface = InteriorBorder_All_Interface;
      // default communication direction
      static const CommunicationDirection defaultDirection = ForwardCommunication;

      /** \brief constructor
       *
       *  \param[in]  gridPart       grid part for the Discontinuous Galerkin space
       *  \param[in]  order          dynamically set maximal polynomial order between 1 and maxPolOrder
       *  \param[in]  commInterface  communication interface to use (optional)
       *  \param[in]  commDirection  communication direction to use (optional)
       */
      explicit PAdaptiveDGSpace ( GridPartType &gridPart,
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
      explicit PAdaptiveDGSpace ( GridPartType &gridPart,
                                  const InterfaceType commInterface = defaultInterface,
                                  const CommunicationDirection commDirection = defaultDirection )
      : BaseType( gridPart, polOrder, commInterface, commDirection )
      {}

      // copy constructor needed for p-adaption
      PAdaptiveDGSpace ( const PAdaptiveDGSpace &other )
      : BaseType( other )
      {}

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      inline bool continuous (const IntersectionType &intersection) const
      {
        return false;
      }

      /** \brief Please doc me. */
      inline const CompiledLocalKeyType &lagrangePointSet( const GeometryType &type, const int order = BaseType::polynomialOrder ) const
      {
        return compiledLocalKey( type, order );
      }

      InterpolationType interpolation() const
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
        return InterpolationImplType( basisFunctionSet( entity ) );
      }

    };

  } // namespace Fem

} // Dune namespace

#endif // #ifndef DUNE_FEM_SPACE_PADAPTIVE_DISCONTINUOUSGALERKIN_HH
