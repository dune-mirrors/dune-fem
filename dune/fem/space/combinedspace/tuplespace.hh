#ifndef DUNE_FEM_SPACE_COMBINEDSPACE_TUPLESPACE_HH
#define DUNE_FEM_SPACE_COMBINEDSPACE_TUPLESPACE_HH

#include <algorithm>
#include <memory>
#include <type_traits>
#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/common/math.hh>

#include <dune/grid/common/grid.hh>

#include <dune/fem/common/memory.hh>
#include <dune/fem/common/utility.hh>
#include <dune/fem/space/basisfunctionset/tuple.hh>
#include <dune/fem/space/combinedspace/generic.hh>
#include <dune/fem/space/combinedspace/interpolation.hh>
#include <dune/fem/space/combinedspace/tuplelocalrestrictprolong.hh>
#include <dune/fem/space/combinedspace/tuplemapper.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/space/common/localinterpolation.hh>

namespace Dune
{

  namespace Fem
  {

    // forward declaration
    template< class CombineOp, class ... DiscreteFunctionSpaces >
    class TupleDiscreteFunctionSpaceImpl;

    // TupleDiscreteFunctionSpaceTraits
    // --------------------------------

    template< class CombineOp, class ... DiscreteFunctionSpaces >
    struct TupleDiscreteFunctionSpaceTraits
    {
      static_assert( sizeof ... ( DiscreteFunctionSpaces ) > 0,
                     "You should provide at least one space to the TupleDiscreteFunctionSpace" );

      // we need to store pointer to the spaces in the SpaceTuple, since space can not be copied.
      typedef std::tuple< std::shared_ptr< DiscreteFunctionSpaces > ... > DiscreteFunctionSpaceTupleType;

    public:
      // helper struct to access contained sub spaces
      template< int i >
      struct SubDiscreteFunctionSpace
      {
        // type of i-th sub space
        typedef typename std::tuple_element< i, DiscreteFunctionSpaceTupleType >::type::element_type Type;

        // type of i-th sub BlockMapper
        typedef typename Type::BlockMapperType BlockMapperType;

        // we will unblock all mappers
        typedef NonBlockMapper< BlockMapperType, Type::localBlockSize > NonBlockMapperType;

        // access to a const ref of the i-th subspace
        static const Type &subDiscreteFunctionSpace ( const DiscreteFunctionSpaceTupleType &tuple )
        {
          assert( std::get< i >( tuple ) );
          return *( std::get< i >( tuple ) );
        }

        static BlockMapperType &subBlockMapper ( const DiscreteFunctionSpaceTupleType &tuple )
        {
          return subDiscreteFunctionSpace( tuple ).blockMapper();
        }

        static NonBlockMapperType subNonBlockMapper ( const DiscreteFunctionSpaceTupleType &tuple  )
        {
          return NonBlockMapperType( subDiscreteFunctionSpace( tuple ).blockMapper() );
        }
      };

    public:
      static_assert( Std::are_all_same< typename DiscreteFunctionSpaces::GridPartType::template Codim< 0 >::EntityType ... >::value,
                     "TupleDiscreteFunctionSpace works only for GridPart's with the same entity type" );

      static_assert( Std::are_all_same< std::integral_constant< int, DiscreteFunctionSpaces::Traits::codimension > ... >::value,
                     "TupleDiscreteFunctionSpace for spaces with different codimensions is not supported" );
      static const int codimension = SubDiscreteFunctionSpace< 0 >::Type::Traits::codimension;

      typedef typename SubDiscreteFunctionSpace< 0 >::Type::GridPartType GridPartType;
      typedef typename GridPartType::GridType GridType;
      typedef typename GridPartType::IndexSetType IndexSetType;
      typedef typename GridPartType::template Codim< 0 >::IteratorType IteratorType;
      typedef typename IteratorType::Entity EntityType;
      typedef typename GridPartType::IntersectionType IntersectionType;

      // type of this space
      typedef TupleDiscreteFunctionSpaceImpl< CombineOp, DiscreteFunctionSpaces ... > DiscreteFunctionSpaceType;

      //! implementation of basefunction set
      typedef TupleBasisFunctionSet< CombineOp, typename DiscreteFunctionSpaces::BasisFunctionSetType ... > BasisFunctionSetType;

      // mapper
      typedef TupleMapper< GridPartType, NonBlockMapper< typename DiscreteFunctionSpaces::BlockMapperType, DiscreteFunctionSpaces::localBlockSize > ... > BlockMapperType;

      // in the most general case we will unroll all local blockings
      typedef std::index_sequence< 0 > LocalBlockIndices;

      // type functionspace
      typedef typename BasisFunctionSetType::FunctionSpaceType FunctionSpaceType;

      typedef TupleSpaceInterpolation< CombineOp, DiscreteFunctionSpaces ... > InterpolationImplType;

      // review to make it work for all kind of combinations
      template< class DiscreteFunction, class Operation = DFCommunicationOperation::Copy >
      struct CommDataHandle
      {
        //! type of data handle
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
        //! type of operatation to perform on scatter
        typedef Operation OperationType;
      };

      // construct new instance of blockMapper
      static BlockMapperType *getBlockMapper ( const DiscreteFunctionSpaceTupleType &spaceTuple )
      {
        return getBlockMapper( spaceTuple, std::index_sequence_for< DiscreteFunctionSpaces ... >() );
      }

      // create Tuple of contained subspaces
      static DiscreteFunctionSpaceTupleType createSpaces ( GridPartType &gridPart, InterfaceType commInterface,
                                                           CommunicationDirection commDirection )
      {
        DiscreteFunctionSpaceTupleType tuple;
        Hybrid::forEach( std::make_index_sequence< sizeof ... ( DiscreteFunctionSpaces ) >{},
          [ & ]( auto i )
          {
            typedef typename SubDiscreteFunctionSpace< i >::Type Element;
            std::get< i >( tuple ) = std::make_shared< Element >( gridPart, commInterface, commDirection );
          } );
        return tuple;
      }

      template< class Entity >
      static BasisFunctionSetType getBasisFunctionSet ( const Entity &entity, const DiscreteFunctionSpaceTupleType &tuple )
      {
        return getBasisFunctionSet( entity, tuple, std::index_sequence_for< DiscreteFunctionSpaces ... >() );
      }

      template< class T, class F >
      static T accumulate ( const DiscreteFunctionSpaceTupleType &tuple, T value, F &&f )
      {
        Hybrid::forEach( std::index_sequence_for< DiscreteFunctionSpaces... >{}, [ & ] ( auto &&idx ) {
            const std::size_t i = std::decay_t< decltype( idx ) >::value;
            value = f( value, SubDiscreteFunctionSpace< i >::subDiscreteFunctionSpace( tuple ) );
          } );
        return value;
      }

    protected:
      template< std::size_t ... i >
      static BlockMapperType *getBlockMapper ( const DiscreteFunctionSpaceTupleType &tuple, std::index_sequence< i ... > )
      {
        return new BlockMapperType( SubDiscreteFunctionSpace< 0 >::subDiscreteFunctionSpace( tuple ).gridPart(),
                                    SubDiscreteFunctionSpace< i >::subNonBlockMapper( tuple ) ... );
      }

      template< class Entity, std::size_t ... i >
      static BasisFunctionSetType getBasisFunctionSet ( const Entity &entity, const DiscreteFunctionSpaceTupleType &tuple,
                                                        std::index_sequence< i ... > )
      {
        return BasisFunctionSetType( SubDiscreteFunctionSpace< i >::subDiscreteFunctionSpace( tuple ).basisFunctionSet( entity ) ... );
      }
    };



    /** \addtogroup DiscreteFunctionSpace
     *
     *  Provides a DiscreteFunctionSpace combined from arbitrary number of DiscreteFunctionSpaces
     *  of different types into a single \ref Dune::Fem::DiscreteFunctionSpaceInterface ( U_h times V_h times .... ).
     *
     *  CombineOp describes the way in which the spaces have been combined, options are:
     *  TupleSpaceProduct:    V = V_1 x V_2 x ...
     *  TupleSpaceSummation:  V = V_1 + V_2 + ...
     */

    /** \class   DiscreteFunctionSpace
     *  \ingroup DiscreteFunctionSpace
     *  \brief    discrete function space
     */
    template< class CombineOp, class ... DiscreteFunctionSpaces >
    class TupleDiscreteFunctionSpaceImpl
      : public GenericCombinedDiscreteFunctionSpace< TupleDiscreteFunctionSpaceTraits< CombineOp, DiscreteFunctionSpaces ... > >
    {
      typedef TupleDiscreteFunctionSpaceImpl< CombineOp, DiscreteFunctionSpaces ... > ThisType;
      typedef GenericCombinedDiscreteFunctionSpace< TupleDiscreteFunctionSpaceTraits< CombineOp, DiscreteFunctionSpaces ... > > BaseType;

    public:
      typedef decltype ( std::index_sequence_for< DiscreteFunctionSpaces ... >() ) Sequence;
      typedef typename BaseType::Traits Traits;
      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::EntityType EntityType;

      typedef typename Traits::InterpolationImplType InterpolationImplType;

      typedef LocalInterpolationWrapper< ThisType > InterpolationType;
      typedef typename Traits::DiscreteFunctionSpaceTupleType DiscreteFunctionSpaceTupleType;

      /** \brief constructor
       *
       *  \param[in]  gridPart       reference to the grid part
       *  \param[in]  commInterface  communication interface to use (optional)
       *  \param[in]  commDirection  communication direction to use (optional)
       *
       */
      TupleDiscreteFunctionSpaceImpl ( GridPartType &gridPart,
                                       const InterfaceType commInterface = InteriorBorder_All_Interface,
                                       const CommunicationDirection commDirection = ForwardCommunication )
        : BaseType( gridPart, commInterface, commDirection )
      {}

      /** \brief constructor
       *
       *  \param[in]  spaces        list of move constructable spaces
       *
       *  \note Gridparts, communication interfaces and directions are assumed to be consistent in all spaces.
       *  Otherwise the behaviour of this space is undefined.
       */
      TupleDiscreteFunctionSpaceImpl ( DiscreteFunctionSpaces &&... spaces )
        : BaseType( std::make_tuple( std::make_shared( std::move( spaces ) )... ) )
      {}

      /** \brief constructor
       *
       *  \param[in]  spaces        list of references to spaces
       *
       *  \note Gridparts, communication interfaces and directions are assumed to be consistent in all spaces.
       *  Otherwise the behaviour of this space is undefined.
       */
      TupleDiscreteFunctionSpaceImpl ( const DiscreteFunctionSpaces &... spaces )
        : BaseType( std::make_tuple( referenceToSharedPtr( spaces )... ) )
      {}

      /** \brief constructor
       *
       *  \param[in]  spaces        list of shared pointers to spaces
       *
       *  \note Gridparts, communication interfaces and directions are assumed to be consistent in all spaces.
       *  Otherwise the behaviour of this space is undefined.
       */
      TupleDiscreteFunctionSpaceImpl ( std::shared_ptr< const DiscreteFunctionSpaces >... spaces )
        : BaseType( std::make_tuple( std::move( spaces )... ) )
      {}

      /** \brief constructor
       *
       *  \param[in]  spaceTuple     tuple of unique pointers to the subspaces
       *
       *  \note Gridparts, communication interfaces and directions are assumed to be consistent in all spaces.
       *  Otherwise the behaviour of this space is undefined.
       */
      explicit TupleDiscreteFunctionSpaceImpl ( DiscreteFunctionSpaceTupleType spaceTuple )
        : BaseType( std::move( spaceTuple ) )
      {}

      TupleDiscreteFunctionSpaceImpl ( const ThisType & ) = delete;
      ThisType &operator= ( const ThisType & ) = delete;

      //! return tuple of const References to the contained sub spaces
      std::tuple< const DiscreteFunctionSpaces & ... > spaceTuple () const
      {
        return spaceTuple( std::index_sequence_for< DiscreteFunctionSpaces ... >() );
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
        return localInterpolation( entity, std::index_sequence_for< DiscreteFunctionSpaces ... >() );
      }

    protected:
      template< std::size_t ... i >
      std::tuple< const DiscreteFunctionSpaces & ... > spaceTuple ( std::index_sequence< i ... > ) const
      {
        return std::tuple< const DiscreteFunctionSpaces & ... >( BaseType::template subDiscreteFunctionSpace< i >() ... );
      }

      template< std::size_t ... i >
      InterpolationImplType localInterpolation ( const EntityType &entity, std::index_sequence< i ... > ) const
      {
        return InterpolationImplType( std::get< i >( spaceTuple() ) ..., entity );
      }
    };



    // DifferentDiscreteFunctionSpace
    // ------------------------------

    //! specialization of DifferentDiscreteFunctionSpace for TupleDiscreteFunctionSpace
    template< class CombineOp, class ... DiscreteFunctionSpaces, class NewFunctionSpace >
    struct DifferentDiscreteFunctionSpace< TupleDiscreteFunctionSpaceImpl< CombineOp, DiscreteFunctionSpaces... >, NewFunctionSpace >
    {
      static_assert( (NewFunctionSpace::dimRange % TupleDiscreteFunctionSpaceImpl< CombineOp, DiscreteFunctionSpaces... >::dimRange == 0),
                     "DifferentDiscreteFunctionSpace can only be applied to TupleFunctionSpace, if new dimRange is a multiple of the original one." );

    private:
      static const int factor = (NewFunctionSpace::dimRange / TupleDiscreteFunctionSpaceImpl< CombineOp, DiscreteFunctionSpaces... >::dimRange);

      template< class DiscreteFunctionSpace >
      using NewSubFunctionSpace = typename ToNewDimRangeFunctionSpace< NewFunctionSpace, factor*DiscreteFunctionSpace::dimRange >::Type;

    public:
      typedef TupleDiscreteFunctionSpaceImpl< CombineOp, typename DifferentDiscreteFunctionSpace< DiscreteFunctionSpaces, NewSubFunctionSpace< DiscreteFunctionSpaces > >::Type... > Type;
    };



    // DefaultLocalRestrictProlong
    // ---------------------------

    template< class CombineOp, class ... DiscreteFunctionSpaces >
    class DefaultLocalRestrictProlong< TupleDiscreteFunctionSpaceImpl< CombineOp, DiscreteFunctionSpaces ... > >
      : public TupleLocalRestrictProlong< DiscreteFunctionSpaces ... >
    {
      typedef DefaultLocalRestrictProlong< TupleDiscreteFunctionSpaceImpl< CombineOp, DiscreteFunctionSpaces ... > > ThisType;
      typedef TupleDiscreteFunctionSpaceImpl< CombineOp, DiscreteFunctionSpaces ... > DiscreteFunctionSpacesType;
      typedef TupleLocalRestrictProlong< DiscreteFunctionSpaces ... > BaseType;

    public:
      DefaultLocalRestrictProlong ( const DiscreteFunctionSpacesType &space )
        : BaseType( space.spaceTuple() )
      {}

    };

    // Creating a space V = V_1 x V_2 x ...
    template < class ... DiscreteFunctionSpaces >
    using TupleDiscreteFunctionSpace = TupleDiscreteFunctionSpaceImpl< TupleSpaceProduct, DiscreteFunctionSpaces ... >;

    // Creating a space V = V_1 + V_2 + ...
    template < class ... DiscreteFunctionSpaces >
    using SummationDiscreteFunctionSpace = TupleDiscreteFunctionSpaceImpl< TupleSpaceSummation, DiscreteFunctionSpaces ... >;

    template < class ... DiscreteFunctionSpaces >
    using EnrichedDiscreteFunctionSpace = SummationDiscreteFunctionSpace< DiscreteFunctionSpaces ... >;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_COMBINEDSPACE_TUPLESPACE_HH
