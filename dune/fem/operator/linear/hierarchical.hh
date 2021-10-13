#ifndef DUNE_FEM_OPERATOR_LINEAR_HIERARCHICAL_HH
#define DUNE_FEM_OPERATOR_LINEAR_HIERARCHICAL_HH

#include <cstddef>

#include <algorithm>
#include <string>
#include <type_traits>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>

#if HAVE_DUNE_ISTL
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>
#endif // #if HAVE_DUNE_ISTL

#include <dune/fem/common/hybrid.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/matrix/functor.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declaration
    // ----------------------------

    template< class DomainFunction, class RangeFunction >
    class HierarchicalLinearOperator;



    namespace Impl
    {

      template< class Dof, class DomainIndices, class RangeIndices >
      struct HierarchicalMatrixChooser;

#if HAVE_DUNE_ISTL
      template< class Dof, int rows, int cols >
      struct HierarchicalMatrixChooser< Dof, Hybrid::IndexRange< int, cols >, Hybrid::IndexRange< int, rows > >
      {
        typedef BCRSMatrix< FieldMatrix< Dof, rows, cols > > Type;
      };

      template< class Dof, class... SD, class... SR >
      struct HierarchicalMatrixChooser< Dof, Hybrid::CompositeIndexRange< SD... >, Hybrid::CompositeIndexRange< SR... > >
      {
      private:
        template< class R >
        using Row = MultiTypeBlockVector< typename HierarchicalMatrixChooser< Dof, SD, R >::Type... >;

      public:
        typedef MultiTypeBlockMatrix< Row< SR >... > Type;
      };
#endif // #if HAVE_DUNE_ISTL

    } // namespace Impl



    // HierarchicalLinearOperator
    // --------------------------

    template< class DomainFunction, class RangeFunction >
    class HierarchicalLinearOperator
      : public AssembledOperator< DomainFunction, RangeFunction >
    {
      typedef HierarchicalLinearOperator< DomainFunction, RangeFunction > ThisType;
      typedef AssembledOperator< DomainFunction, RangeFunction > BaseType;

    public:
      typedef std::common_type_t< typename DomainFunction::DofType, typename RangeFunction::DofType > DofType;

      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;
      typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeSpaceType;

      typedef typename DomainSpaceType::EntityType DomainEntityType;
      typedef typename RangeSpaceType::EntityType RangeEntityType;

      // typedef DomainEntityType ColumnEntityType;
      // typedef RangeEntityType RowEntityType;

      typedef typename Impl::HierarchicalMatrixChooser< DofType, typename DomainSpaceType::LocalBlockIndices, typename RangeSpaceType::LocalBlockIndices >::Type MatrixType;

    private:
      template< class Functor >
      static auto blockFunctor ( Functor &&functor )
      {
        return [ functor ] ( std::pair< std::size_t, std::size_t > local, const auto &global ) {
            local.first *= Hybrid::size( typename RangeSpaceType::LocalBlockIndices() );
            local.second *= Hybrid::size( typename DomainSpaceType::LocalBlockIndices() );
            Hybrid::forEach( typename RangeSpaceType::LocalBlockIndices(), [ functor, local, global ] ( auto &&i ) {
                Hybrid::forEach( typename DomainSpaceType::LocalBlockIndices(), [ functor, local, global, i ] ( auto &&j ) {
                    const auto iGlobal = std::make_pair( static_cast< std::size_t >( global.first ), i );
                    const auto jGlobal = std::make_pair( static_cast< std::size_t >( global.second ), j );
                    functor( std::make_pair( local.first + i, local.second + j ), std::make_pair( iGlobal, jGlobal ) );
                  } );
              } );
          };
      }

    protected:
      MatrixType &matrix () { return matrix_; }

    public:
      HierarchicalLinearOperator ( const std::string &, const DomainSpaceType &domainSpace, const RangeSpaceType &rangeSpace )
        : domainSpace_( domainSpace ), rangeSpace_( rangeSpace )
      {}

      virtual void operator() ( const DomainFunction &u, RangeFunction &w ) const
      {
        w.clear();
        umv( matrix_, u.dofVector().array(), w.dofVector().array() );
        w.communicate();
      }

      void communicate () {}

      const DomainSpaceType &domainSpace () const { return domainSpace_; }
      const RangeSpaceType &rangeSpace () const { return domainSpace_; }

      MatrixType &exportMatrix () const { return matrix_; }

      template< class LocalMatrix >
      void addLocalMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity, const LocalMatrix &localMatrix )
      {
        auto f = blockFunctor( [ this, &localMatrix ] ( auto local, auto global ) {
            ThisType::entry( matrix_, global.first, global.second ) += localMatrix[ local.first ][ local.second ];
          } );
        rangeSpace().blockMapper().mapEach( rangeEntity, makePairFunctor( domainSpace().blockMapper(), domainEntity, std::move( f ) ) );
      }

      template< class LocalMatrix, class Scalar >
      void addScaledLocalMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity, const LocalMatrix &localMatrix, const Scalar &scalar )
      {
        auto f = blockFunctor( [ this, &localMatrix, &scalar ] ( auto local, auto global ) {
            ThisType::entry( matrix_, global.first, global.second ) += scalar * localMatrix[ local.first ][ local.second ];
          } );
        rangeSpace().blockMapper().mapEach( rangeEntity, makePairFunctor( domainSpace().blockMapper(), domainEntity, std::move( f ) ) );
      }

      template< class LocalMatrix >
      void getLocalMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity, LocalMatrix &localMatrix ) const
      {
        auto f = blockFunctor( [ this, &localMatrix ] ( auto local, auto global ) {
            localMatrix[ local.first ][ local.second ] = ThisType::entry( matrix_, global.first, global.second );
          } );
        rangeSpace().blockMapper().mapEach( rangeEntity, makePairFunctor( domainSpace().blockMapper(), domainEntity, std::move( f ) ) );
      }

      template< class LocalMatrix >
      void setLocalMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity, const LocalMatrix &localMatrix )
      {
        auto f = blockFunctor( [ this, &localMatrix ] ( auto local, auto global ) {
            ThisType::entry( matrix_, global.first, global.second ) = localMatrix[ local.first ][ local.second ];
          } );
        rangeSpace().blockMapper().mapEach( rangeEntity, makePairFunctor( domainSpace().blockMapper(), domainEntity, std::move( f ) ) );
      }

      void clear () { clear( matrix_ ); }
      template <class I>
      void unitRow( const I localRow, const double diag = 1.0 )
      {
        DUNE_THROW(NotImplemented,"unitRow not implemented on HierarchicalLinearOperator");
      }

      template< class Stencil >
      void reserve ( const Stencil &stencil )
      {
        reserve( matrix_, stencil );
      }

    private:
#if HAVE_DUNE_ISTL
      template< class... C, class... B, class RangeVector >
      static std::enable_if_t< sizeof...( C ) == sizeof...( B ) > umv ( const MultiTypeBlockVector< C... > &row, const MultiTypeBlockVector< B... > &u, RangeVector &w )
      {
        Hybrid::forEach( std::index_sequence_for< C... >(), [ &row, &u, &w ] ( auto &&j ) {
            ThisType::umv( row[ j ], u[ j ], w );
          } );
      }

      template< class... R, class DomainVector, class... B >
      static std::enable_if_t< sizeof...( R ) == sizeof...( B ) > umv ( const MultiTypeBlockMatrix< R... > &matrix, const DomainVector &u, MultiTypeBlockVector< B... > &w )
      {
        Hybrid::forEach( std::index_sequence_for< R... >(), [ &matrix, &u, &w ] ( auto &&i ) {
            ThisType::umv( matrix[ i ], u, w[ i ] );
          } );
      }

      template< class K, int m, int n, class AM, class AU, class AW >
      static void umv ( const BCRSMatrix< FieldMatrix< K, m, n >, AM > &matrix, const BlockVector< FieldVector< K, n >, AU > &u, BlockVector< FieldVector< K, m >, AW > &w )
      {
        matrix.umv( u, w );
      }

      template< class... C >
      static void clear ( MultiTypeBlockVector< C... > &row )
      {
        Hybrid::forEach( std::index_sequence_for< C... >(), [ &row ] ( auto &&j ) {
            ThisType::clear( row[ j ] );
          } );
      }

      template< class... R >
      static void clear ( MultiTypeBlockMatrix< R... > &matrix )
      {
        Hybrid::forEach( std::index_sequence_for< R... >(), [ &matrix ] ( auto &&i ) {
            ThisType::clear( matrix[ i ] );
          } );
      }

      template< class K, int m, int n, class A >
      static void clear ( BCRSMatrix< FieldMatrix< K, m, n >, A > &matrix )
      {
        for( auto &row : matrix )
          std::fill( row.begin(), row.end(), FieldMatrix< K, m, n >( K( 0 ) ) );
      }

      template< class... C, class I, std::size_t component, class J, J offset, class SJ >
      static decltype( auto ) entry ( const MultiTypeBlockVector< C... > &row, I i, std::pair< std::size_t, Hybrid::CompositeIndex< component, J, offset, SJ > > j )
      {
        return entry( row[ std::integral_constant< std::size_t, component >() ], i, std::make_pair( j.first, j.second.subIndex() ) );
      }

      template< class... C, class I, std::size_t component, class J, J offset, class SJ >
      static decltype( auto ) entry ( MultiTypeBlockVector< C... > &row, I i, std::pair< std::size_t, Hybrid::CompositeIndex< component, J, offset, SJ > > j )
      {
        return entry( row[ std::integral_constant< std::size_t, component >() ], i, std::make_pair( j.first, j.second.subIndex() ) );
      }

      template< class... R, std::size_t component, class I, I offset, class SI, class J >
      static decltype( auto ) entry ( const MultiTypeBlockMatrix< R... > &matrix, std::pair< std::size_t, Hybrid::CompositeIndex< component, I, offset, SI > > i, J j )
      {
        return entry( matrix[ std::integral_constant< std::size_t, component >() ], std::make_pair( i.first, i.second.subIndex() ), j );
      }

      template< class... R, std::size_t component, class I, I offset, class SI, class J >
      static decltype( auto ) entry ( MultiTypeBlockMatrix< R... > &matrix, std::pair< std::size_t, Hybrid::CompositeIndex< component, I, offset, SI > > i, J j )
      {
        return entry( matrix[ std::integral_constant< std::size_t, component >() ], std::make_pair( i.first, i.second.subIndex() ), j );
      }

      template< class K, int m, int n, class A >
      static const K &entry ( const BCRSMatrix< FieldMatrix< K, m, n >, A > &matrix, std::pair< std::size_t, int > i, std::pair< std::size_t, int > j )
      {
        return matrix[ i.first ][ j.first ][ i.second ][ j.second ];
      }

      template< class K, int m, int n, class A >
      static K &entry ( BCRSMatrix< FieldMatrix< K, m, n >, A > &matrix, std::pair< std::size_t, int > i, std::pair< std::size_t, int > j )
      {
        return matrix[ i.first ][ j.first ][ i.second ][ j.second ];
      }

      template< class... C, class Stencil >
      static void reserve ( MultiTypeBlockVector< C... > &row, const Stencil &stencil )
      {
        Hybrid::forEach( std::index_sequence_for< C... >(), [ &row, &stencil ] ( auto &&j ) {
            ThisType::reserve( row[ j ], stencil );
          } );
      }

      template< class... R, class Stencil >
      static void reserve ( MultiTypeBlockMatrix< R... > &matrix, const Stencil &stencil )
      {
        Hybrid::forEach( std::index_sequence_for< R... >(), [ &matrix, &stencil ] ( auto &&i ) {
            ThisType::reserve( matrix[ i ], stencil );
          } );
      }

      template< class K, int m, int n, class A, class Stencil >
      static void reserve ( BCRSMatrix< FieldMatrix< K, m, n >, A > &matrix, const Stencil &stencil )
      {
        // reallocate matrix of correct size (destroys previously allocated matrix)
        if( matrix.buildMode() == matrix.unknown )
          matrix.setBuildMode( matrix.random );
        matrix.setSize( stencil.rows(), stencil.cols() );

        // setup sparsity pattern (using random build mode)
        const auto& globalStencil = stencil.globalStencil();
        const std::size_t nRows = globalStencil.size();
        for( std::size_t row = 0; row<nRows; ++row )
        {
          try {
            const auto& columns = globalStencil.at( row );
            matrix.setrowsize( row, columns.size() );
          }
          catch ( const std::out_of_range& e )
          {
            // if a row is empty then do nothing
            continue ;
          }
        }
        matrix.endrowsizes();

        for( std::size_t row = 0; row<nRows; ++row )
        {
          try {
            const auto& columns = globalStencil.at( row );
            matrix.setIndices( row, columns.begin(), columns.end() );
          }
          catch ( const std::out_of_range& e )
          {
            // if a row is empty then do nothing
            continue ;
          }
        }
        matrix.endindices();
      }
#endif // #if HAVE_DUNE_ISTL

      const DomainSpaceType &domainSpace_;
      const RangeSpaceType &rangeSpace_;
      mutable MatrixType matrix_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_OPERATOR_LINEAR_HIERARCHICAL_HH
