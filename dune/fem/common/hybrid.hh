#ifndef DUNE_FEM_COMMON_HYBRID_HH
#define DUNE_FEM_COMMON_HYBRID_HH

#include <cstddef>

#include <tuple>
#include <type_traits>
#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/common/typetraits.hh>

#include <dune/fem/common/utility.hh>

namespace Dune
{

  namespace Hybrid
  {

    namespace Impl
    {

      template< class Range, class = void >
      struct Index;

      template< class T, T... i >
      struct Index< std::integer_sequence< T, i... > >
      {
        typedef T Type;
      };

      template< class Range >
      struct Index< Range, void_t< typename Range::Index > >
      {
        typedef typename Range::Index Type;
      };


      template< class Index, class = void >
      struct FlatIndex;

      template< class Index >
      struct FlatIndex< Index, std::enable_if_t< std::is_integral< Index >::value > >
      {
        typedef Index Type;
      };

      template< class T, T value >
      struct FlatIndex< std::integral_constant< T, value >, void >
      {
        typedef T Type;
      };

      template< class Index >
      struct FlatIndex< Index, void_t< typename Index::FlatIndex > >
      {
        typedef typename Index::FlatIndex Type;
      };

    } // namespace Impl



    // IndexType
    // ---------

    template< class Range >
    using IndexType = typename Impl::Index< Range >::Type;



    // FlatIndexType
    // -------------

    template< class Index >
    using FlatIndexType = typename Impl::FlatIndex< Index >::Type;



    // IndexRange
    // ----------

    template< class T, T sz >
    struct IndexRange
    {
      typedef T Index;

      static constexpr Index size () { return sz; }
    };



    // CompositeIndexRange
    // -------------------

    template< class... SR >
    struct CompositeIndexRange
    {
      typedef std::common_type_t< IndexType< SR >... > Index;

      template< std::size_t i >
      using SubRange = std::tuple_element_t< i, std::tuple< SR... > >;

      template< std::size_t i >
      static constexpr Index offset ( std::integral_constant< std::size_t, i > = {} )
      {
        return size( std::make_index_sequence< i >() );
      }

      static constexpr Index size () { return size( std::index_sequence_for< SR... >() ); }

    private:
      static constexpr Index size ( std::index_sequence<> ) { return 0; }

      template< std::size_t... i >
      static constexpr Index size ( std::index_sequence< i... > )
      {
        return Std::sum( Hybrid::size( SubRange< i >() )... );
      }
    };



    // forEach for IndexRange
    // ----------------------

    template< class T, T sz, class F >
    inline static void forEach ( IndexRange< T, sz > range, F &&f )
    {
      for( T i = 0; i < sz; ++i )
        f( i );
    }



    // CompositeIndex
    // --------------

    template< std::size_t component, class I, I offset, class SI >
    struct CompositeIndex
    {
      typedef I FlatIndex;
      typedef SI SubIndex;

    private:
      static constexpr std::integral_constant< std::size_t, component > access ( SubIndex subIndex, std::integral_constant< std::size_t, 0 > ) { return {}; }

      static constexpr FlatIndexType< SubIndex > access ( FlatIndexType< SubIndex > subIndex, std::integral_constant< std::size_t, 1 > ) { return subIndex; }

      template< std::size_t _component, class _I, _I _offset, class _SI, std::size_t i >
      static constexpr decltype( auto ) access ( CompositeIndex< _component, _I, _offset, _SI > subIndex, std::integral_constant< std::size_t, i > )
      {
        return subIndex[ std::integral_constant< std::size_t, i-1 >() ];
      }

      template< std::size_t i >
      using Access = decltype( access( std::declval< SubIndex >(), std::integral_constant< std::size_t, i >() ) );

    public:
      explicit constexpr CompositeIndex ( SubIndex subIndex ) : subIndex_( std::move( subIndex ) ) {}

      constexpr operator I () const { return (offset + static_cast< FlatIndexType< SubIndex > >( subIndex() )); }

      template< std::size_t i >
      constexpr auto operator[] ( std::integral_constant< std::size_t, i > ) const
        -> std::enable_if_t< !IsIntegralConstant< Access< i > >::value, Access< i > >
      {
        return access( subIndex(), std::integral_constant< std::size_t, i >() );
      }

      template< std::size_t i >
      constexpr auto operator[] ( std::integral_constant< std::size_t, i > ) const
        -> std::enable_if_t< IsIntegralConstant< Access< i > >::value, std::decay_t< Access< i > > >
      {
        return {};
      }

      const SubIndex &subIndex () const { return subIndex_; }

    private:
      SubIndex subIndex_;
    };



    // forEach for CompositeIndexRange
    // -------------------------------

    template< class... SR, class F >
    inline static void forEach ( CompositeIndexRange< SR... >, F &&f );

    namespace Impl
    {

      template< class Range, class F >
      struct CompositeIndexRangeInnerLoop
      {
        explicit CompositeIndexRangeInnerLoop ( F f )
          : f_( std::move( f ) )
        {}

        template< std::size_t component >
        void operator() ( std::integral_constant< std::size_t, component > )
        {
          typedef IndexType< Range > Index;
          typedef typename Range::template SubRange< component > SubRange;

          Hybrid::forEach( SubRange(), [ this ] ( auto subIndex ) {
              f_( CompositeIndex< component, Index, Range::template offset< component >(), std::decay_t< decltype( subIndex ) > >( std::move( subIndex ) ) );
            } );
        }

      private:
        F f_;
      };

    } // namespace Impl

    template< class... SR, class F >
    inline static void forEach ( CompositeIndexRange< SR... >, F &&f )
    {
      Impl::CompositeIndexRangeInnerLoop< CompositeIndexRange< SR... >, F > innerLoop( std::forward< F >( f ) );
      forEach( std::index_sequence_for< SR... >(), innerLoop );
    }

  } // namespace Hybrid

} // namespace Dune

#endif // #ifndef DUNE_FEM_COMMON_HYBRID_HH
