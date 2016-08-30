#ifndef DUNE_FEM_SPACE_COMBINEDSPACE_POWERSPACE_HH
#define DUNE_FEM_SPACE_COMBINEDSPACE_POWERSPACE_HH

#include <algorithm>
#include <memory>
#include <type_traits>

#include <dune/common/math.hh>
#include <dune/common/std/memory.hh>

#include <dune/grid/common/grid.hh>

#include <dune/fem/common/utility.hh>

#include <dune/fem/space/basisfunctionset/vectorial.hh>
#include <dune/fem/space/combinedspace/generic.hh>
#include <dune/fem/space/combinedspace/interpolation.hh>
#include <dune/fem/space/combinedspace/lagrangepointsetexporter.hh>
#include <dune/fem/space/combinedspace/powerlocalrestrictprolong.hh>
#include <dune/fem/space/combinedspace/powermapper.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>

namespace Dune
{

  namespace Fem
  {

    // forward declaration

    template< class DiscreteFunctionSpace, int N >
    class PowerDiscreteFunctionSpace;


    // PowerDiscreteFunctionSpaceTraits
    // --------------------------------

    template< class DiscreteFunctionSpace, int N >
    struct PowerDiscreteFunctionSpaceTraits
    {
      // we need to store pointer to the spaces in the SpaceTuple, since space can not be copied.
      typedef std::unique_ptr< DiscreteFunctionSpace > DiscreteFunctionSpaceTupleType;

      // helper struct to access contained sub spaces
      template< int >
      struct SubDiscreteFunctionSpace
      {
        // type of i-th sub space
        typedef DiscreteFunctionSpace Type;

        // type of i-th sub BlockMapper
        typedef typename Type::BlockMapperType BlockMapperType;

        // access to a const ref of the i-th subspace
        static const Type &subDiscreteFunctionSpace ( const DiscreteFunctionSpaceTupleType &tuple )
        {
          return (*tuple);
        }

        static BlockMapperType &subBlockMapper ( const DiscreteFunctionSpaceTupleType &tuple )
        {
          return subDiscreteFunctionSpace().blockMapper();
        }
      };

      // type of GridPart
      typedef typename DiscreteFunctionSpace::GridPartType GridPartType;
      typedef typename GridPartType::GridType GridType;
      typedef typename GridPartType::IndexSetType IndexSetType;
      typedef typename GridPartType::template Codim< 0 >::IteratorType IteratorType;
      typedef typename IteratorType::Entity EntityType;
      typedef typename GridPartType::IntersectionType IntersectionType;

      // type of this space
      typedef PowerDiscreteFunctionSpace< DiscreteFunctionSpace, N > DiscreteFunctionSpaceType;

      static const int codimension = DiscreteFunctionSpace::Traits::codimension;

      static_assert( DiscreteFunctionSpace::FunctionSpaceType::dimRange == 1,
                     "PowerDiscreteFunctionSpace only works for ContainedSpaces with dimRange = 1" );
      //! implementation of basefunction set
      typedef typename DiscreteFunctionSpace::BasisFunctionSetType ScalarBasisFunctionSetType;
      typedef typename ScalarBasisFunctionSetType::FunctionSpaceType::RangeFieldType RangeFieldType;
      typedef VectorialBasisFunctionSet< ScalarBasisFunctionSetType, FieldVector< RangeFieldType, N >, VerticalDofAlignment > BasisFunctionSetType;

      // type of block mapper
      typedef PowerMapper< GridPartType, typename DiscreteFunctionSpace::BlockMapperType, N > BlockMapperType;

      // in the most general case we will unroll all local blockings
      enum { localBlockSize = DiscreteFunctionSpace::localBlockSize };

      // type functionspace
      typedef typename BasisFunctionSetType::FunctionSpaceType FunctionSpaceType;

      static constexpr int polynomialOrder = DiscreteFunctionSpace::polynomialOrder;

      // type of local Interpolation
      typedef PowerSpaceInterpolation< DiscreteFunctionSpace, N > InterpolationType;

      // review to make it work for all kind of combinations
      template< class DiscreteFunction,
                class Operation =
                  typename DiscreteFunctionSpace::template CommDataHandle< DiscreteFunctionSpace >::OperationType >
      struct CommDataHandle
      {
        //! type of data handle
        typedef typename DiscreteFunctionSpace::
          template CommDataHandle< DiscreteFunction, Operation >::Type Type;
        //! type of operatation to perform on scatter
        typedef typename DiscreteFunctionSpace::
          template CommDataHandle< DiscreteFunction, Operation >::OperationType OperationType;
      };

      // construct new instance of blockMapper
      static BlockMapperType *getBlockMapper ( const DiscreteFunctionSpaceTupleType &spaceTuple )
      {
        return new BlockMapperType( spaceTuple->gridPart(), spaceTuple->blockMapper() );
      }

      // create Tuple of contained subspaces
      static DiscreteFunctionSpaceTupleType createSpaces ( GridPartType &gridPart, InterfaceType commInterface,
                                                           CommunicationDirection commDirection )
      {
        return Std::make_unique< DiscreteFunctionSpace >( gridPart, commInterface, commDirection );
      }

      template< class Entity >
      static BasisFunctionSetType getBasisFunctionSet ( const Entity &entity, const DiscreteFunctionSpaceTupleType &tuple )
      {
        return BasisFunctionSetType( tuple->basisFunctionSet( entity ) );
      }

      static bool continuous ( const DiscreteFunctionSpaceTupleType &tuple )
      {
        return tuple->continuous();
      }

      static bool continuous ( const IntersectionType &intersection, const DiscreteFunctionSpaceTupleType &tuple )
      {
        return tuple->continuous( intersection );
      }
    };



    /** \addtogroup DiscreteFunctionSpace
     *
     *  Provides a DiscreteFunctionSpace combined from arbitrary number of DiscreteFunctionSpaces
     *  of same type into a single \ref Dune::Fem::DiscreteFunctionSpaceInterface ( U_h times V_h times .... ).
     *
     *  \note It is assumed that the each space is build upon the same gridpart
     */

    /** \class   DiscreteFunctionSpace
     *  \ingroup DiscreteFunctionSpace
     *  \brief    discrete function space
     */
    template< class DiscreteFunctionSpace, int N >
    class PowerDiscreteFunctionSpace
      : public GenericCombinedDiscreteFunctionSpace< PowerDiscreteFunctionSpaceTraits< DiscreteFunctionSpace, N > >,
        public CombinedSpaceHelper::LagrangePointSetExporter< DiscreteFunctionSpace >
    {
      typedef PowerDiscreteFunctionSpace< DiscreteFunctionSpace, N > ThisType;
      typedef GenericCombinedDiscreteFunctionSpace< PowerDiscreteFunctionSpaceTraits< DiscreteFunctionSpace, N > > BaseType;
      typedef CombinedSpaceHelper::LagrangePointSetExporter< DiscreteFunctionSpace > LagrangePointSetExporterType;

    public:
      typedef PowerDiscreteFunctionSpaceTraits< DiscreteFunctionSpace, N > Traits;
      //! extract grid informations, it is assumed the both spaces are living on the
      //! same gridPart
      typedef typename Traits::GridPartType GridPartType;

      //! type of contained discrete function space
      typedef DiscreteFunctionSpace ContainedDiscreteFunctionSpaceType;

      //! type of interpolation
      typedef typename Traits::InterpolationType InterpolationType;

      typedef typename Traits::EntityType EntityType;

      /** \brief constructor
       *
       *  \param[in]  gridPart       grid part
       *  \param[in]  commInterface  communication interface to use (optional)
       *  \param[in]  commDirection  communication direction to use (optional)
       */
      PowerDiscreteFunctionSpace ( GridPartType &gridPart,
                                   const InterfaceType commInterface = InteriorBorder_All_Interface,
                                   const CommunicationDirection commDirection = ForwardCommunication )
        : BaseType( gridPart, commInterface, commDirection ),
          LagrangePointSetExporterType( containedSpace() )
      {}

      PowerDiscreteFunctionSpace ( const ThisType& ) = delete;
      ThisType& operator= ( const ThisType& ) = delete;

      //! return tuple of const References to the contained sub spaces
      const ContainedDiscreteFunctionSpaceType &containedSpace () const
      {
        return BaseType::template subDiscreteFunctionSpace< 0 >();
      }

      InterpolationType interpolation ( const EntityType &entity ) const
      {
        return InterpolationType( entity );
      }
    };


    //! specialization of DifferentDiscreteFunctionSpace for PowerDiscreteFunctionSpace
    template< class DiscreteFunctionSpace, int  N, class NewFunctionSpace >
    struct DifferentDiscreteFunctionSpace< PowerDiscreteFunctionSpace< DiscreteFunctionSpace, N >, NewFunctionSpace >
    {
      typedef PowerDiscreteFunctionSpace< DiscreteFunctionSpace, NewFunctionSpace::dimRange > Type;
    };


    // DefaultLocalRestrictProlong
    // ---------------------------

    template< class DiscreteFunctionSpace, int N >
    class DefaultLocalRestrictProlong< PowerDiscreteFunctionSpace< DiscreteFunctionSpace, N > >
      : public PowerLocalRestrictProlong< DiscreteFunctionSpace, N >
    {
      typedef DefaultLocalRestrictProlong< PowerDiscreteFunctionSpace< DiscreteFunctionSpace, N > > ThisType;
      typedef PowerDiscreteFunctionSpace< DiscreteFunctionSpace, N > DiscreteFunctionSpacesType;
      typedef PowerLocalRestrictProlong< DiscreteFunctionSpace, N > BaseType;

    public:
      DefaultLocalRestrictProlong ( const DiscreteFunctionSpacesType &space )
        : BaseType( space.containedSpace() )
      {}

    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_COMBINEDSPACE_POWERSPACE_HH
