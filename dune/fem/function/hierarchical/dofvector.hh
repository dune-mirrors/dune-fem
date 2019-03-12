#ifndef DUNE_FEM_FUNCTION_HIERARCHICAL_DOFVECTOR_HH
#define DUNE_FEM_FUNCTION_HIERARCHICAL_DOFVECTOR_HH

#include <type_traits>
#include <utility>

#include <dune/common/ftraits.hh>

#if HAVE_DUNE_ISTL
#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockvector.hh>
#endif // #if HAVE_DUNE_ISTL

#include <dune/fem/common/hybrid.hh>
#include <dune/fem/common/utility.hh>
#include <dune/fem/function/blockvectors/defaultblockvectors.hh>

namespace Dune
{

  namespace Fem
  {

    namespace Impl
    {

      template< class T >
      struct BlockIndicesFor;

#if HAVE_DUNE_ISTL
      template< class K, int dim, class A >
      struct BlockIndicesFor< BlockVector< FieldVector< K, dim >, A > >
      {
        typedef Hybrid::IndexRange< int, dim > Type;
      };

      template< class... V >
      struct BlockIndicesFor< MultiTypeBlockVector< V... > >
      {
        typedef Hybrid::CompositeIndexRange< typename BlockIndicesFor< V >::Type... > Type;
      };

      template< class T >
      struct BlockIndicesFor< const T >
      {
        typedef typename BlockIndicesFor< T >::Type Type;
      };
#endif // HAVE_DUNE_ISTL

    } // namespace Impl




    // HierarchicalDofBlock
    // --------------------

    template< class DofContainer >
    struct HierarchicalDofBlock
    {
      friend struct HierarchicalDofBlock< const DofContainer >;

      typedef typename Impl::BlockIndicesFor< DofContainer >::Type BlockIndices;
      static constexpr std::size_t blockSize = Hybrid::size( BlockIndices() );

      HierarchicalDofBlock ( DofContainer &data, std::size_t baseIndex )
        : data_( data ), baseIndex_( baseIndex )
      {}

      HierarchicalDofBlock ( const HierarchicalDofBlock< std::remove_const_t< DofContainer > > &other )
        : data_( other.data_ ), baseIndex_( other.baseIndex_ )
      {}

      template< class Index >
      decltype( auto ) operator[] ( const Index &index ) const
      {
        return access( data_, baseIndex_, index );
      }

      template< class Index >
      decltype( auto ) operator[] ( const Index &index )
      {
        return access( data_, baseIndex_, index );
      }

    private:
#if HAVE_DUNE_ISTL
      template< class... V, std::size_t component, class I, I offset, class SI >
      static decltype( auto ) access ( const MultiTypeBlockVector< V... > &data, std::size_t baseIndex, const Hybrid::CompositeIndex< component, I, offset, SI > &index )
      {
        return access( data[ std::integral_constant< std::size_t, component >() ], baseIndex, index.subIndex() );
      }

      template< class... V, std::size_t component, class I, I offset, class SI >
      static decltype( auto ) access ( MultiTypeBlockVector< V... > &data, std::size_t baseIndex, const Hybrid::CompositeIndex< component, I, offset, SI > &index )
      {
        return access( data[ std::integral_constant< std::size_t, component >() ], baseIndex, index.subIndex() );
      }

      template< class K, int n, class A >
      static const K &access ( const BlockVector< FieldVector< K, n >, A > &data, std::size_t baseIndex, int index )
      {
        return data[ baseIndex ][ index ];
      }

      template< class K, int n, class A >
      static K &access ( BlockVector< FieldVector< K, n >, A > &data, std::size_t baseIndex, int index )
      {
        return data[ baseIndex ][ index ];
      }
#endif // #if HAVE_DUNE_ISTL

      DofContainer &data_;
      std::size_t baseIndex_;
    };


    // SpecialArrayFeature for HierarhicalDofVector
    // --------------------------------------------

    template< class >
    struct SpecialArrayFeatures;


    // HierarchicalDofVector
    // ---------------------

    template< class DC >
    class HierarchicalDofVector
      : public IsBlockVector
    {
      typedef HierarchicalDofVector< DC > ThisType;

    public:
      typedef DC DofContainerType;

      typedef typename FieldTraits< DofContainerType >::field_type FieldType;

      typedef std::size_t SizeType;

      typedef FieldType value_type;
      typedef SizeType size_type;

      typedef typename Impl::BlockIndicesFor< DC >::Type BlockIndices;
      static constexpr std::size_t blockSize = Hybrid::size( BlockIndices() );

      typedef HierarchicalDofBlock< const DofContainerType > ConstDofBlockType;
      typedef HierarchicalDofBlock< DofContainerType > DofBlockType;

      explicit HierarchicalDofVector ( SizeType size ) { resize( size ); }

      ConstDofBlockType operator[] ( SizeType i ) const { return ConstDofBlockType( data_, i ); }
      DofBlockType operator[] ( SizeType i ) { return DofBlockType( data_, i ); }

      ThisType &operator+= ( const ThisType &other ) { data_ += other.data_; return *this; }
      ThisType &operator-= ( const ThisType &other ) { data_ -= other.data_; return *this; }
      ThisType &operator*= ( const FieldType &scalar ) { data_ *= scalar; return *this; }

      FieldType operator* ( const ThisType &other ) const { return data_ * other.data_; }

      void axpy ( const FieldType &scalar, const ThisType &other ) { data_.axpy( scalar, other.data_ ); }

      void clear () { data_ = FieldType( 0 ); }

      /** \brief obtain underlaying DoF storage **/
      const DofContainerType &array () const noexcept { return data_; }
      /** \brief obtain underlaying DoF storage **/
      DofContainerType &array () noexcept { return data_; }

      void reserve ( SizeType size ) { reserve( data_, size ); }
      void resize ( SizeType size ) { resize( data_, size ); }

      SizeType size () const { return size( data_ ); }

      // unimplemented interface methods

      typedef const FieldType *ConstIteratorType;
      typedef FieldType *IteratorType;

      ConstIteratorType begin () const { DUNE_THROW( NotImplemented, "HierarchicalDofVector does not provide iterators" ); }
      IteratorType begin () { DUNE_THROW( NotImplemented, "HierarchicalDofVector does not provide iterators" ); }

      ConstIteratorType end () const { DUNE_THROW( NotImplemented, "HierarchicalDofVector does not provide iterators" ); }
      IteratorType end () { DUNE_THROW( NotImplemented, "HierarchicalDofVector does not provide iterators" ); }

      std::size_t usedMemorySize() const
      {
        return SpecialArrayFeatures< ThisType >::used( *this );
      }

      void setMemoryFactor ( double memFactor )
      {}

      void memMoveBackward ( int length, int oldStartIdx, int newStartIdx )
      {
        SpecialArrayFeatures< ThisType >::memMoveBackward( *this, length, oldStartIdx, newStartIdx );
      }

      void memMoveForward ( int length, int oldStartIdx, int newStartIdx )
      {
        SpecialArrayFeatures< ThisType >::memMoveForward( *this, length, oldStartIdx, newStartIdx );
      }

      void copyContent ( int newIndex, int oldIndex )
      {
        SpecialArrayFeatures< ThisType >::assign( *this, newIndex, oldIndex );
      }

    private:
#if HAVE_DUNE_ISTL
      template< class... V >
      static void reserve ( MultiTypeBlockVector< V... > &data, SizeType size )
      {
        Hybrid::forEach( std::index_sequence_for< V... >(), [ &data, size ] ( auto &&i ) { ThisType::reserve( data[ i ], size ); } );
      }

      template< class B, class A >
      static void reserve ( BlockVector< B, A > &data, SizeType size )
      {
        data.reserve( size );
      }

      template< class... V >
      static void resize ( MultiTypeBlockVector< V... > &data, SizeType size )
      {
        Hybrid::forEach( std::index_sequence_for< V... >(), [ &data, size ] ( auto &&i ) { ThisType::resize( data[ i ], size ); } );
      }

      template< class B, class A >
      static void resize ( BlockVector< B, A > &data, SizeType size )
      {
        data.resize( size );
      }

      template< class... V >
      static SizeType size ( const MultiTypeBlockVector< V... > &data )
      {
        return data[ std::integral_constant< std::size_t, 0 > () ].size();
      }

      template< class B, class A >
      static  SizeType size ( const BlockVector< B, A > &data )
      {
        return data.size();
      }
#endif // HAVE_DUNE_ISTL

      DofContainerType data_;
    };



    // SpecialArrayFeature for HierarhicalDofVector
    // --------------------------------------------

    template< class >
    struct SpecialArrayFeatures;

    template< class DC >
    struct SpecialArrayFeatures< HierarchicalDofVector< DC > >
    {
      static std::size_t used ( const HierarchicalDofVector< DC > &array )
      {
        return used( array.array() );
      }

      static void setMemoryFactor ( HierarchicalDofVector< DC > &array, double memFactor )
      {}

      static void memMoveBackward ( HierarchicalDofVector< DC > &array, int length, int oldStartIdx, int newStartIdx )
      {
        memMoveBackward( array.array(), length, oldStartIdx, newStartIdx );
      }

      static void memMoveForward ( HierarchicalDofVector< DC > &array, int length, int oldStartIdx, int newStartIdx )
      {
        memMoveForward( array.array(), length, oldStartIdx, newStartIdx );
      }

      static void assign ( HierarchicalDofVector< DC > &array, int newIndex, int oldIndex )
      {
        assign( array.array(), newIndex, oldIndex );
      }

    private:
#if HAVE_DUNE_ISTL
      template< class... V >
      static std::size_t used ( const MultiTypeBlockVector< V... > &array )
      {
        std::size_t used( 0 );
        Hybrid::forEach( std::index_sequence_for< V... >(), [ &array, &used ] ( auto &&i ) {
            used += SpecialArrayFeatures< HierarchicalDofVector< DC > >::used( array[ i ] );
          } );
        return used;
      }

      template< class B, class A >
      static std::size_t used ( const BlockVector< B, A > &array )
      {
        return array.size() * sizeof( B );
      }

      template< class... V >
      static void memMoveBackward ( MultiTypeBlockVector< V... > &array, int length, int oldStartIdx, int newStartIdx )
      {
        Hybrid::forEach( std::index_sequence_for< V... >(), [ &array, length, oldStartIdx, newStartIdx ] ( auto &&i ) {
            SpecialArrayFeatures< HierarchicalDofVector< DC > >::memMoveBackward( array[ i ], length, oldStartIdx, newStartIdx );
          } );
      }

      template< class B, class A >
      static void memMoveBackward ( BlockVector< B, A > &array, int length, int oldStartIdx, int newStartIdx )
      {
        for( int oldIdx = oldStartIdx+length-1, newIdx = newStartIdx + length-1; oldIdx >= oldStartIdx; --oldIdx, --newIdx )
          array[ newIdx ] = array[ oldIdx ];
      }

      template< class... V >
      static void memMoveForward ( MultiTypeBlockVector< V... > &array, int length, int oldStartIdx, int newStartIdx )
      {
        Hybrid::forEach( std::index_sequence_for< V... >(), [ &array, length, oldStartIdx, newStartIdx ] ( auto &&i ) {
            SpecialArrayFeatures< HierarchicalDofVector< DC > >::memMoveForward( array[ i ], length, oldStartIdx, newStartIdx );
          } );
      }

      template< class B, class A >
      static void memMoveForward ( BlockVector< B, A > &array, int length, int oldStartIdx, int newStartIdx )
      {
        for( int oldIdx = oldStartIdx, newIdx = newStartIdx; oldIdx < oldStartIdx+length; ++oldIdx, ++newIdx )
          array[ newIdx ] = array[ oldIdx ];
      }

      template< class... V >
      static void assign ( MultiTypeBlockVector< V... > &array, int newIndex, int oldIndex )
      {
        Hybrid::forEach( std::index_sequence_for< V... >(), [ &array, newIndex, oldIndex ] ( auto &&i ) {
            SpecialArrayFeatures< HierarchicalDofVector< DC > >::assign( array[ i ], newIndex, oldIndex );
          } );
      }

      template< class B, class A >
      static void assign ( BlockVector< B, A > &array, int newIndex, int oldIndex )
      {
        array[ newIndex ] = array[ oldIndex ];
      }
#endif // #if HAVE_DUNE_ISTL
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_HIERARCHICAL_DOFVECTOR_HH
