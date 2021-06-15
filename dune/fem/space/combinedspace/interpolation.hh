#ifndef DUNE_FEM_SPACE_COMBINEDSPACE_INTERPOLATION_HH
#define DUNE_FEM_SPACE_COMBINEDSPACE_INTERPOLATION_HH

#include <tuple>
#include <utility>

#include <dune/fem/common/forloop.hh>
#include <dune/fem/function/localfunction/converter.hh>
#include <dune/fem/space/basisfunctionset/tuple.hh>
#include <dune/fem/space/basisfunctionset/vectorial.hh>
#include <dune/fem/storage/subvector.hh>


namespace Dune
{

  namespace Fem
  {

    // PowerSpaceInterpolation
    // ------------------

    template< class Space, int N >
    class PowerSpaceInterpolation
    {
      typedef PowerSpaceInterpolation< Space, N > ThisType;

      struct RangeConverter
      {
        RangeConverter ( std::size_t range ) : range_( range ) {}

        template< class T >
        FieldVector< T, 1 > operator() ( const FieldVector< T, N > &in ) const
        {
          return in[ range_ ];
        }

        template< class T, int j >
        FieldMatrix< T, 1, j > operator() ( const FieldMatrix< T, N, j > &in ) const
        {
          return in[ range_ ];
        }

      protected:
        std::size_t range_;
      };

      template <class DofVector, class DofAlignment>
      struct SubDofVectorWrapper
        : public SubDofVector< DofVector, DofAlignment >
      {
        typedef SubDofVector< DofVector, DofAlignment > BaseType;

        SubDofVectorWrapper( DofVector& dofs, int coordinate, const DofAlignment &dofAlignment )
          : BaseType( dofs, coordinate, dofAlignment )
        {}

        //! do nothing on clear since it's done in apply of this class
        void clear() {}
      };

      // Note: BasisFunctionSetType is VectorialBasisFunctionSet
      typedef typename Space::BasisFunctionSetType::DofAlignmentType DofAlignmentType;

    public:
      typedef typename Space::EntityType EntityType;

      PowerSpaceInterpolation ( const Space &space, const EntityType &entity )
        : interpolation_( space.containedSpace().localInterpolation( entity ) ),
          dofAlignment_( space.basisFunctionSet( entity ).dofAlignment() )
      {}

      template< class LocalFunction, class LocalDofVector >
      void operator () ( const LocalFunction &lv, LocalDofVector &ldv ) const
      {
        apply( lv, ldv );
      }

      template< class LocalFunction, class LocalDofVector >
      void apply ( const LocalFunction &lv, LocalDofVector &ldv ) const
      {
        // clear dofs before something is adedd
        ldv.clear();

        for( std::size_t i = 0; i < N; ++i )
        {
          SubDofVectorWrapper< LocalDofVector, DofAlignmentType > subLdv( ldv, i, dofAlignment_ );
          interpolation_( localFunctionConverter( lv, RangeConverter( i ) ), subLdv );
        }
      }

    protected:
      typename Space::ContainedDiscreteFunctionSpaceType::InterpolationImplType interpolation_;
      DofAlignmentType dofAlignment_;
    };


    // TupleSpaceInterpolation
    // -----------------------
    // CombineOp describes the way in which the spaces have been combined, i.e.
    // Product:    V = V_1 x V_2 x ...
    // Summation:  V = V_1 + V_2 + ...

    template< class CombineOp , class ... Spaces >
    class TupleSpaceInterpolation
    {
      typedef TupleSpaceInterpolation< CombineOp, Spaces ... > ThisType;
      typedef std::tuple< typename Spaces::InterpolationImplType ... > InterpolationTupleType;

      static const int setSize = sizeof ... ( Spaces ) -1;

      template< int >
      struct ProductApply;

      template< int >
      struct SummationApply;

      template< int, class CombOp >
      struct ApplyBase;

      template< int counter >
      struct ApplyBase< counter, TupleSpaceProduct > : public ProductApply< counter >{};

      template< int counter >
      struct ApplyBase< counter, TupleSpaceSummation > : public SummationApply< counter >{};

      template < int counter >
      struct Apply : public ApplyBase< counter, CombineOp > {};

      typedef TupleBasisFunctionSet< CombineOp, typename Spaces::BasisFunctionSetType ... >   BasisFunctionSetType;

    public:

      static_assert( Std::are_all_same< typename Spaces::EntityType ... >::value,
                     "TupleSpaceInterpolation requires Spaces defined over the same grid" );

      typedef typename std::tuple_element< 0, std::tuple< Spaces ... > >::type::EntityType EntityType;

      TupleSpaceInterpolation ( std::tuple< const Spaces & ... > tuple, const EntityType &entity )
        : interpolation_( interpolationTuple( tuple, entity ) ),
          basisFunctionSet_( basisFunctionSetTuple( tuple, entity ) )
      {}

      TupleSpaceInterpolation ( const Spaces & ... spaces, const EntityType &entity )
        : interpolation_( std::make_tuple( spaces.localInterpolation( entity ) ... ) ),
          basisFunctionSet_( std::make_tuple( spaces.basisFunctionSet( entity ) ... )  )
      {}

      template< class LocalFunction, class LocalDofVector >
      void operator() ( const LocalFunction &lf, LocalDofVector &ldv ) const
      {
        Fem::ForLoop< Apply, 0, setSize >::apply( interpolation_, basisFunctionSet_, lf, ldv );
      }

      void unbind() {}

    protected:
      InterpolationTupleType interpolation_;
      BasisFunctionSetType basisFunctionSet_;
    };


    // specialization of TupleSpaceInterpolation::Apply for SpaceProduct
    template< class CombineOp, class ... Spaces >
    template< int i >
    struct TupleSpaceInterpolation< CombineOp, Spaces ... >::ProductApply
    {
      static const int rangeOffset = BasisFunctionSetType::RangeIndices::template offset< i >();
      static const int thisDimRange = BasisFunctionSetType::template SubBasisFunctionSet< i >::type::FunctionSpaceType::dimRange;
      static const int dimRange = BasisFunctionSetType::FunctionSpaceType::dimRange;

      struct RangeConverter
      {
        template< class T >
        FieldVector< T, thisDimRange > operator() ( const FieldVector< T, dimRange > &in ) const
        {
          FieldVector< T, thisDimRange > ret;
          apply( in, ret );
          return ret;
        }

        template< class T, int j >
        FieldMatrix< T, thisDimRange, j > operator() ( const FieldMatrix< T, dimRange, j > &in ) const
        {
          FieldMatrix< T, thisDimRange, j > ret;
          apply( in, ret );
          return ret;
        }

      protected:
        template< class In, class Out >
        void apply ( const In &in, Out &out ) const
        {
          for( std::size_t j = 0; j < thisDimRange; ++j )
            out[ j ] = in[ j + rangeOffset ];
        }
      };

      template< class Tuple, class LocalFunction, class LocalDofVector >
      static void apply ( const Tuple &tuple, const BasisFunctionSetType &basisSet, const LocalFunction &lv, LocalDofVector &ldv )
      {
        SubVector< LocalDofVector, OffsetSubMapper >
          subLdv( ldv, OffsetSubMapper( basisSet.template subBasisFunctionSet< i >().size(), basisSet.offset( i ) ) );
        std::get< i >( tuple ) ( localFunctionConverter( lv, RangeConverter() ), subLdv );
      }
    };


    template< class CombineOp, class ... Spaces >
    template< int i >
    struct TupleSpaceInterpolation< CombineOp, Spaces ... >::SummationApply
    {
      template< class Tuple, class LocalFunction, class LocalDofVector >
      static void apply ( const Tuple &tuple, const BasisFunctionSetType &basisSet, const LocalFunction &lv, LocalDofVector &ldv )
      {
        SubVector< LocalDofVector, OffsetSubMapper >
          subLdv( ldv, OffsetSubMapper( basisSet.template subBasisFunctionSet< i >().size(), basisSet.offset( i ) ) );
        std::get< i >( tuple ) ( lv,  subLdv );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_COMBINEDSPACE_INTERPOLATION_HH
