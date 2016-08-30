#ifndef DUNE_FEM_SPACE_COMBINEDSPACE_TUPLESPACE_HH
#define DUNE_FEM_SPACE_COMBINEDSPACE_TUPLESPACE_HH

#include <algorithm>
#include <memory>
#include <type_traits>

#include <dune/common/math.hh>
#include <dune/common/std/memory.hh>

#include <dune/grid/common/grid.hh>

#include <dune/fem/common/utility.hh>
#include <dune/fem/space/basisfunctionset/tuple.hh>
#include <dune/fem/space/combinedspace/generic.hh>
#include <dune/fem/space/combinedspace/interpolation.hh>
#include <dune/fem/space/combinedspace/tuplelocalrestrictprolong.hh>
#include <dune/fem/space/combinedspace/tuplemapper.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>

namespace Dune
{

  namespace Fem
  {

    // forward declaration

    template< class ... DiscreteFunctionSpaces >
    class TupleDiscreteFunctionSpace;


    // TupleDiscreteFunctionSpaceTraits
    // --------------------------------

    template< class ... DiscreteFunctionSpaces >
    struct TupleDiscreteFunctionSpaceTraits
    {
      static_assert( sizeof ... ( DiscreteFunctionSpaces ) > 0,
                     "You should provide at least one space to the TupleDiscreteFunctionSpace" );

      // we need to store pointer to the spaces in the SpaceTuple, since space can not be copied.
      typedef std::tuple< std::unique_ptr< DiscreteFunctionSpaces > ... > DiscreteFunctionSpaceTupleType;

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

    protected:
      // helper struct to create each entry in the tuple
      template< int i >
      struct Constructor
      {
        template< class Tuple, class ... Args >
        static void apply ( Tuple &tuple, Args && ... args )
        {
          typedef typename SubDiscreteFunctionSpace< i >::Type Element;
          std::get< i >( tuple ) = Std::make_unique< Element >( std::forward< Args >( args ) ... );
        }
      };

    public:
      static_assert( Std::are_all_same< std::integral_constant< int, DiscreteFunctionSpaces::Traits::codimension > ... >::value,
                     "TupleDiscreteFunctionSpace for spaces with different codimensions is not supported" );
      static const int codimension = SubDiscreteFunctionSpace< 0 >::Type::Traits::codimension;

      static_assert( Std::are_all_same< typename DiscreteFunctionSpaces::GridPartType ... >::value,
                     "TupleDiscreteFunctionSpace works only for common GridPartTypes" );
      // type of GridPart
      typedef typename SubDiscreteFunctionSpace< 0 >::Type::GridPartType GridPartType;
      typedef typename GridPartType::GridType GridType;
      typedef typename GridPartType::IndexSetType IndexSetType;
      typedef typename GridPartType::template Codim< 0 >::IteratorType IteratorType;
      typedef typename IteratorType::Entity EntityType;
      typedef typename GridPartType::IntersectionType IntersectionType;

      // type of this space
      typedef TupleDiscreteFunctionSpace< DiscreteFunctionSpaces ... > DiscreteFunctionSpaceType;

      //! implementation of basefunction set
      typedef TupleBasisFunctionSet< typename DiscreteFunctionSpaces::BasisFunctionSetType ... > BasisFunctionSetType;

      // mapper
      typedef TupleMapper< GridPartType, NonBlockMapper< typename DiscreteFunctionSpaces::BlockMapperType, DiscreteFunctionSpaces::localBlockSize > ... > BlockMapperType;

      // in the most general case we will unroll all local blockings
      enum { localBlockSize = 1 };

      // type functionspace
      typedef typename BasisFunctionSetType::FunctionSpaceType FunctionSpaceType;

      static constexpr int polynomialOrder = Std::max( (int)DiscreteFunctionSpaces::polynomialOrder ... );

      typedef TupleSpaceInterpolation< DiscreteFunctionSpaces ... > InterpolationType;

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
        ForLoop< Constructor, 0, sizeof ... ( DiscreteFunctionSpaces ) -1 >::apply( tuple, gridPart, commInterface, commDirection );
        return tuple;
      }

      template< class Entity >
      static BasisFunctionSetType getBasisFunctionSet ( const Entity &entity, const DiscreteFunctionSpaceTupleType &tuple )
      {
        return getBasisFunctionSet( entity, tuple, std::index_sequence_for< DiscreteFunctionSpaces ... >() );
      }

      static bool continuous ( const DiscreteFunctionSpaceTupleType &tuple )
      {
        return continuous( tuple, std::index_sequence_for< DiscreteFunctionSpaces ... >() );
      }

      static bool continuous ( const IntersectionType &intersection, const DiscreteFunctionSpaceTupleType &tuple )
      {
        return continuous( tuple, intersection, std::index_sequence_for< DiscreteFunctionSpaces ... >() );
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

      template< std::size_t ... i >
      static bool continuous ( const DiscreteFunctionSpaceTupleType &tuple, std::index_sequence< i ... > )
      {
        return Std::And( SubDiscreteFunctionSpace< i >::subDiscreteFunctionSpace( tuple ).continuous() ... );
      }

      template< std::size_t ... i >
      static bool continuous ( const DiscreteFunctionSpaceTupleType &tuple, const IntersectionType &intersection, std::index_sequence< i ... > )
      {
        return Std::And( SubDiscreteFunctionSpace< i >::subDiscreteFunctionSpace( tuple ).continuous( intersection ) ... );
      }
    };



    /** \addtogroup DiscreteFunctionSpace
     *
     *  Provides a DiscreteFunctionSpace combined from arbitrary number of DiscreteFunctionSpaces
     *  of different types into a single \ref Dune::Fem::DiscreteFunctionSpaceInterface ( U_h times V_h times .... ).
     *
     *  \note It is assumed that the each space is build upon the same gridpart
     */

    /** \class   DiscreteFunctionSpace
     *  \ingroup DiscreteFunctionSpace
     *  \brief    discrete function space
     */
    template< class ... DiscreteFunctionSpaces >
    class TupleDiscreteFunctionSpace
      : public GenericCombinedDiscreteFunctionSpace< TupleDiscreteFunctionSpaceTraits< DiscreteFunctionSpaces ... > >
    {
      typedef TupleDiscreteFunctionSpace< DiscreteFunctionSpaces ... > ThisType;
      typedef GenericCombinedDiscreteFunctionSpace< TupleDiscreteFunctionSpaceTraits< DiscreteFunctionSpaces ... > > BaseType;
      typedef TupleDiscreteFunctionSpaceTraits< DiscreteFunctionSpaces ... > Traits;

    public:
      //! extract grid informations, it is assumed the both spaces are living on the
      //! same gridPart
      typedef typename Traits::GridPartType GridPartType;

      typedef typename Traits::InterpolationType InterpolationType;
      typedef typename Traits::EntityType EntityType;

      /** \brief constructor
       *
       *  \param[in]  gridPart       grid part
       *  \param[in]  commInterface  communication interface to use (optional)
       *  \param[in]  commDirection  communication direction to use (optional)
       */
      TupleDiscreteFunctionSpace ( GridPartType &gridPart,
                                   const InterfaceType commInterface = InteriorBorder_All_Interface,
                                   const CommunicationDirection commDirection = ForwardCommunication )
        : BaseType( gridPart, commInterface, commDirection )
      {}

      TupleDiscreteFunctionSpace ( const ThisType& ) = delete;
      ThisType& operator= ( const ThisType& ) = delete;

      //! return tuple of const References to the contained sub spaces
      std::tuple< const DiscreteFunctionSpaces & ... > spaceTuple () const
      {
        return spaceTuple( std::index_sequence_for< DiscreteFunctionSpaces ... >() );
      }

      InterpolationType interpolation ( const EntityType &entity ) const
      {
        return interpolation( entity, std::index_sequence_for< DiscreteFunctionSpaces ... >() );
      }

    protected:
      template< std::size_t ... i >
      std::tuple< const DiscreteFunctionSpaces & ... > spaceTuple ( std::index_sequence< i ... > ) const
      {
        return std::tuple< const DiscreteFunctionSpaces & ... >( BaseType::template subDiscreteFunctionSpace< i >() ... );
      }

      template< std::size_t ... i >
      InterpolationType interpolation ( const EntityType &entity, std::index_sequence< i ... > ) const
      {
        return InterpolationType( std::get< i >( spaceTuple() ) ..., entity );
      }
    };


    // DefaultLocalRestrictProlong
    // ---------------------------

    template< class ... DiscreteFunctionSpaces >
    class DefaultLocalRestrictProlong< TupleDiscreteFunctionSpace< DiscreteFunctionSpaces ... > >
      : public TupleLocalRestrictProlong< DiscreteFunctionSpaces ... >
    {
      typedef DefaultLocalRestrictProlong< TupleDiscreteFunctionSpace< DiscreteFunctionSpaces ... > > ThisType;
      typedef TupleDiscreteFunctionSpace< DiscreteFunctionSpaces ... > DiscreteFunctionSpacesType;
      typedef TupleLocalRestrictProlong< DiscreteFunctionSpaces ... > BaseType;

    public:
      DefaultLocalRestrictProlong ( const DiscreteFunctionSpacesType &space )
        : BaseType( space.spaceTuple() )
      {}

    };


  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_COMBINEDSPACE_TUPLESPACE_HH
