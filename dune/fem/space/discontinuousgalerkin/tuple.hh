#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_TUPLE_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_TUPLE_HH

#include <tuple>
#include <type_traits>
#include <utility>

#include <dune/common/hybridutilities.hh>

#include <dune/fem/common/hybrid.hh>
#include <dune/fem/common/utility.hh>
#include <dune/fem/function/localfunction/converter.hh>
#include <dune/fem/space/combinedspace/tuplelocalrestrictprolong.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/localinterpolation.hh>
#include <dune/fem/space/discontinuousgalerkin/generic.hh>
#include <dune/fem/space/discontinuousgalerkin/basisfunctionsets.hh>
#include <dune/fem/space/shapefunctionset/tuple.hh>
#include <dune/fem/storage/subvector.hh>


namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class... DFS >
    class TupleDiscontinuousGalerkinSpace;



    // TupleLocalInterpolation
    // -----------------------

    template< class SFS, class... I >
    class TupleLocalInterpolation
    {
      typedef TupleLocalInterpolation< SFS, I... > ThisType;

    public:
      typedef SFS ShapeFunctionSetType;

      static const int dimRange = ShapeFunctionSetType::RangeType::dimension;

    private:
      template< std::size_t i >
      using SubRangeType = typename ShapeFunctionSetType::template SubShapeFunctionSetType< i >::RangeType;

      template< std::size_t... i >
      static constexpr int sumDimSubRange ( std::index_sequence< i... > )
      {
        return Std::sum( SubRangeType< i >::dimension ... );
      }

      static constexpr int sumDimSubRange ( std::index_sequence<> ) { return 0; }

      template< std::size_t i >
      struct RangeConverter
      {
        template< class T >
        FieldVector< T, SubRangeType< i >::dimension > operator() ( const FieldVector< T, dimRange > &in ) const
        {
          FieldVector< T, SubRangeType< i >::dimension > out;
          for( int j = 0; j < SubRangeType< i >::dimension; ++j )
            out[ j ] = in[ j + sumDimSubRange( std::make_index_sequence< i >() ) ];
          return out;
        }

        template< class T, int cols >
        FieldMatrix< T, SubRangeType< i >::dimension, cols > operator() ( const FieldMatrix< T, dimRange, cols > &in ) const
        {
          FieldMatrix< T, SubRangeType< i >::dimension, cols > out;
          for( int j = 0; j < SubRangeType< i >::dimension; ++j )
            out[ j ] = in[ j + sumDimSubRange( std::make_index_sequence< i >() ) ];
          return out;
        }
      };

    public:
      template< class... Args >
      explicit TupleLocalInterpolation ( ShapeFunctionSetType shapeFunctionSet, Args &&... args )
        : shapeFunctionSet_( std::move( shapeFunctionSet ) ),
          interpolations_( std::forward< Args >( args )... )
      {}

      template< class LocalFunction, class LocalDofVector >
      void operator() ( const LocalFunction &lf, LocalDofVector &ldv ) const
      {
        std::size_t offset = 0;
        Hybrid::forEach( std::index_sequence_for< I... >(), [ this, &lf, &ldv, &offset ] ( auto i ) {
            const std::size_t size = shapeFunctionSet_.subShapeFunctionSet( i ).size();
            SubVector< LocalDofVector, OffsetSubMapper > subLdv( ldv, OffsetSubMapper( size, offset ) );
            std::get< i >( interpolations_ ) ( localFunctionConverter( lf, RangeConverter< i >() ), subLdv );
            offset += size;
          } );
      }

      void unbind() {}

    protected:
      ShapeFunctionSetType shapeFunctionSet_;
      std::tuple< I... > interpolations_;
    };



    // TupleDiscontinuousGalerkinSpaceBasisFunctionSets
    // ------------------------------------------------

    template< class... DFS >
    class TupleDiscontinuousGalerkinSpaceBasisFunctionSets
    {
      typedef TupleDiscontinuousGalerkinSpaceBasisFunctionSets< DFS... > ThisType;

      static_assert( sizeof...( DFS ) > 0, "TupleDiscontinuousGalerkinSpace requires at least on space." );

      static_assert( Std::are_all_same< typename DFS::GridPartType... >::value, "TupleDiscontinuousGalerkinSpace only works on spaces with identical GridPart." );

      template< class E, class SFS >
      static std::true_type isDefaultBasisFunctionSet ( DefaultBasisFunctionSet< E, SFS > );

      static std::false_type isDefaultBasisFunctionSet ( ... );

      template< class BFS >
      struct IsDefaultBasisFunctionSet
        : public decltype( isDefaultBasisFunctionSet( std::declval< BFS >() ) )
      {};

      static_assert( Std::And( IsDefaultBasisFunctionSet< typename DFS::BasisFunctionSetType >::value... ), "TupleDiscontinuousGalerkinSpace only works on spaces with DefaultBasisFunctionSets." );

    public:
      template< std::size_t i >
      using SubDiscreteFunctionSpaceType = std::tuple_element_t< i, std::tuple< DFS... > >;

      typedef typename SubDiscreteFunctionSpaceType< 0 >::GridPartType GridPartType;

      typedef TupleShapeFunctionSet< typename DFS::BasisFunctionSetType::ShapeFunctionSetType... > ShapeFunctionSetType;

      static const int codimension = GridPartType::dimension - ShapeFunctionSetType::DomainType::dimension;
      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;

      typedef DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType > BasisFunctionSetType;

      TupleDiscontinuousGalerkinSpaceBasisFunctionSets ( GridPartType &gridPart, InterfaceType commInterface, CommunicationDirection commDirection )
        : subDiscreteFunctionSpaces_( DFS( gridPart, commInterface, commDirection )... )
      {}

      int order () const { return order( std::index_sequence_for< DFS... >() ); }
      int order ( const EntityType &entity ) const { return order( entity, std::index_sequence_for< DFS... >() ); }

      BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const { return BasisFunctionSetType( entity, shapeFunctionSet( entity ) ); }
      ShapeFunctionSetType shapeFunctionSet ( const EntityType &entity ) const { return shapeFunctionSet( entity, std::index_sequence_for< DFS... >() ); }

      template< std::size_t i >
      const SubDiscreteFunctionSpaceType< i > &subDiscreteFunctionSpace ( std::integral_constant< std::size_t, i > = {} ) const
      {
        return std::get< i >( subDiscreteFunctionSpaces_ );
      }

    protected:
      template< class SFS >
      static auto shapeFunctionSet ( const DefaultBasisFunctionSets< GridPartType, SFS > &basisFunctionSets, const EntityType &entity )
      {
        return basisFunctionSets.shapeFunctionSets().shapeFunctionSet( entity.type() );
      }

      template< class BFS >
      static auto shapeFunctionSet ( const BFS &basisFunctionSets, const EntityType &entity )
      {
        return basisFunctionSets.basisFunctionSet( entity ).shapeFunctionSet();
      }

      template< std::size_t... i >
      int order ( std::index_sequence< i... > ) const
      {
        return Std::max( subDiscreteFunctionSpace< i >().order()... );
      }

      template< std::size_t... i >
      int order ( const EntityType &entity, std::index_sequence< i... > ) const
      {
        return Std::max( subDiscreteFunctionSpace< i >().order( entity )... );
      }

      template< std::size_t... i >
      ShapeFunctionSetType shapeFunctionSet ( const EntityType &entity, std::index_sequence< i... > ) const
      {
        return ShapeFunctionSetType( shapeFunctionSet( subDiscreteFunctionSpace< i >().basisFunctionSets(), entity )... );
      }

    private:
      std::tuple< DFS... > subDiscreteFunctionSpaces_;
    };



    // TaylorHoodDiscontinuousGalerkinSpaceTraits
    // ------------------------------------------

    template< class... DFS >
    struct TupleDiscontinuousGalerkinSpaceTraits
    {
      typedef TupleDiscontinuousGalerkinSpace< DFS... > DiscreteFunctionSpaceType;

      typedef TupleDiscontinuousGalerkinSpaceBasisFunctionSets< DFS... > BasisFunctionSetsType;

      typedef typename BasisFunctionSetsType::GridPartType GridPartType;
      typedef typename BasisFunctionSetsType::BasisFunctionSetType BasisFunctionSetType;

      typedef typename BasisFunctionSetType::FunctionSpaceType FunctionSpaceType;

      static const int codimension = BasisFunctionSetsType::codimension;

      typedef Hybrid::CompositeIndexRange< typename DFS::LocalBlockIndices... > LocalBlockIndices;

      static const int polynomialOrder = Std::max( static_cast< int >( DFS::polynomialOrder )... );

      typedef CodimensionMapper< GridPartType, codimension > BlockMapperType;

      template< class DiscreteFunction, class Operation = DFCommunicationOperation::Copy >
      struct CommDataHandle
      {
        typedef Operation OperationType;
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
      };
    };



    // TaylorHoodDiscontinuousGalerkinSpace
    // ------------------------------------

    template< class... DFS >
    class TupleDiscontinuousGalerkinSpace
      : public GenericDiscontinuousGalerkinSpace< TupleDiscontinuousGalerkinSpaceTraits< DFS... > >
    {
      typedef TupleDiscontinuousGalerkinSpace< DFS... > ThisType;
      typedef GenericDiscontinuousGalerkinSpace< TupleDiscontinuousGalerkinSpaceTraits< DFS... > > BaseType;

    public:
      typedef TupleDiscontinuousGalerkinSpaceTraits< DFS... > Traits;

      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;
      typedef typename BaseType::BasisFunctionSetsType BasisFunctionSetsType;
      typedef typename BaseType::EntityType EntityType;
      typedef typename BaseType::GridPartType GridPartType;

      template< std::size_t i >
      using SubDiscreteFunctionSpaceType = typename BasisFunctionSetsType::template SubDiscreteFunctionSpaceType< i >;

      typedef TupleLocalInterpolation< typename BasisFunctionSetsType::ShapeFunctionSetType, typename DFS::InterpolationImplType... > InterpolationImplType;
      typedef LocalInterpolationWrapper< ThisType > InterpolationType;

      using BaseType::basisFunctionSets;

      TupleDiscontinuousGalerkinSpace ( GridPartType &gridPart, InterfaceType commInterface = InteriorBorder_All_Interface, CommunicationDirection commDirection = ForwardCommunication )
        : BaseType( gridPart, BasisFunctionSetsType( gridPart, commInterface, commDirection ), commInterface, commDirection )
      {}

      InterpolationType interpolation () const { return InterpolationType( *this ); }

      [[deprecated]]
      InterpolationImplType interpolation ( const EntityType &entity ) const { return localInterpolation( entity ); }

      InterpolationImplType localInterpolation ( const EntityType &entity ) const { return localInterpolation( entity, std::index_sequence_for< DFS... >() ); }

      template< std::size_t i >
      const SubDiscreteFunctionSpaceType< i > &subDiscreteFunctionSpace ( std::integral_constant< std::size_t, i > = {} ) const
      {
        return basisFunctionSets().subDiscreteFunctionSpace( std::integral_constant< std::size_t, i >() );
      }

    private:
      template< std::size_t... i >
      InterpolationImplType localInterpolation ( const EntityType &entity, std::index_sequence< i... > ) const
      {
        return InterpolationImplType( basisFunctionSets().shapeFunctionSet( entity ), subDiscreteFunctionSpace< i >().localInterpolation( entity )... );
      }
    };



    // DifferentDiscreteFunctionSpace
    // ------------------------------

    template< class... DFS, class NewFunctionSpace >
    struct DifferentDiscreteFunctionSpace< TupleDiscontinuousGalerkinSpace< DFS... >, NewFunctionSpace >
    {
      static_assert( NewFunctionSpace::dimRange % TupleDiscontinuousGalerkinSpace< DFS... >::dimRange == 0,
                     "DifferentDiscreteFunctionSpace can only be applied to TupleDiscontinuousGalerkinSpace, if new dimRange is a multiple of the original one." );

    private:
      static const int factor = (NewFunctionSpace::dimRange / TupleDiscontinuousGalerkinSpace< DFS... >::dimRange);

      template< int dimRange >
      using NewSubFunctionSpace = typename ToNewDimRangeFunctionSpace< NewFunctionSpace, factor*dimRange >::Type;

    public:
      typedef TupleDiscontinuousGalerkinSpace< typename DifferentDiscreteFunctionSpace< DFS, NewSubFunctionSpace< DFS::dimRange > >::Type... > Type;
    };



    // DefaultLocalRestrictProlong
    // ---------------------------

    template< class... DFS >
    class DefaultLocalRestrictProlong< TupleDiscontinuousGalerkinSpace< DFS... > >
      : public TupleLocalRestrictProlong< DFS... >
    {
      typedef DefaultLocalRestrictProlong< TupleDiscontinuousGalerkinSpace< DFS... > > ThisType;
      typedef TupleLocalRestrictProlong< DFS... > BaseType;

    public:
      typedef TupleDiscontinuousGalerkinSpace< DFS... > DiscreteFunctionSpaceType;

      DefaultLocalRestrictProlong ( const DiscreteFunctionSpaceType &space )
        : BaseType( subDiscreteFunctionSpaces( space, std::index_sequence_for< DFS... >() ) )
      {}

    private:
      template< std::size_t... i >
      static std::tuple< const DFS &... > subDiscreteFunctionSpaces ( const DiscreteFunctionSpaceType &space, std::index_sequence< i... > )
      {
        return std::tie( space.subDiscreteFunctionSpace( std::integral_constant< std::size_t, i >() )... );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_TUPLE_HH
