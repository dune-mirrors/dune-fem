#ifndef DUNE_GEOMETRY_TYPEINDEXSET_HH
#define DUNE_GEOMETRY_TYPEINDEXSET_HH

#include <cstddef>

#include <dune/geometry/type.hh>

namespace Dune
{

  namespace Fem
  {
    namespace hpDG {
    // LocalGeometryTypeIndexSet
    // -------------------------

    /** \class LocalGeometryTypeIndexSet
     *
     *  \brief Please doc me.
     *
     *  \tparam  dim      dimension
     *  \tparam  regular  include geometry type 'none'
     */
    template< int dim, bool regular = false >
    class LocalGeometryTypeIndexSet;

    template< int dim >
    class LocalGeometryTypeIndexSet< dim, true >
    {
      typedef LocalGeometryTypeIndexSet< dim, true > This;

    public:
      inline static constexpr std::size_t size () noexcept
      {
        return (1 << dim) - ((1 << dim) >> 1);
      }

      inline static constexpr std::size_t index ( const GeometryType &type ) noexcept
      {
        return (type.id() >> 1);
      }

      inline static constexpr bool contains ( const GeometryType &type ) noexcept
      {
        return ((type.dim() == dim) && !type.isNone());
      }

      inline static GeometryType type ( std::size_t index ) noexcept
      {
        return GeometryType( static_cast< unsigned int >( index ) << 1, dim );
      }
    };

    template< int dim >
    class LocalGeometryTypeIndexSet< dim, false >
    {
      typedef LocalGeometryTypeIndexSet< dim, false > This;
      typedef LocalGeometryTypeIndexSet< dim, true > RegularTypeIndexSet;

    public:
      inline static constexpr std::size_t size () noexcept
      {
        return (RegularTypeIndexSet::size() + 1);
      }

      inline static constexpr std::size_t index ( const GeometryType &type ) noexcept
      {
        return (type.isNone() ? RegularTypeIndexSet::size() : type.id() >> 1);
      }

      inline static constexpr bool contains ( const GeometryType &type ) noexcept
      {
        return (type.dim() == dim);
      }

      inline static GeometryType type ( std::size_t index ) noexcept
      {
        return (index < RegularTypeIndexSet::size() ? RegularTypeIndexSet::type( index ) : GeometryType( 0, dim, true ));
      }
    };



    // GlobalGeometryTypeIndexSet
    // --------------------------

    /** \class GlobalGeometryTypeIndexSet
     *
     *  \brief Please doc me.
     *
     *  \tparam  maxdim   maximum dimension
     *  \tparam  regular  include geometry type 'none'
     */
    template< int maxdim, bool regular = false >
    class GlobalGeometryTypeIndexSet;

    template< int maxdim >
    class GlobalGeometryTypeIndexSet< maxdim, true >
    {
      typedef GlobalGeometryTypeIndexSet< maxdim, true > This;

    public:
      inline static constexpr std::size_t size () noexcept
      {
        return (1 << maxdim);
      }

      inline static constexpr std::size_t index ( const GeometryType &type ) noexcept
      {
        return ((1 << type.dim()) + type.id()) >> 1;
      }

      inline static constexpr bool contains ( const GeometryType &type ) noexcept
      {
        return ((type.dim() <= maxdim) && !type.isNone());
      }

      inline static GeometryType type ( std::size_t index ) noexcept
      {
        return GeometryType( (index << 1) & ~(1 << dim( index )), dim( index ) );
      }

    private:
      inline static constexpr int dim ( std::size_t index, int d = maxdim )
      {
        return ((d <= 0) || ((index & (1 << (d-1))) != 0) ? d : dim( index, d-1 ));
      }
    };

    template< int maxdim >
    class GlobalGeometryTypeIndexSet< maxdim, false >
    {
      typedef GlobalGeometryTypeIndexSet< maxdim, false > This;
      typedef GlobalGeometryTypeIndexSet< maxdim, true > RegularTypeIndexSet;

    public:
      inline static constexpr std::size_t size () noexcept
      {
        return RegularTypeIndexSet::size() + (maxdim + 1);
      }

      inline static constexpr std::size_t index ( const GeometryType &type ) noexcept
      {
        return (type.isNone() ? RegularTypeIndexSet::size() + type.dim() : RegularTypeIndexSet::index( type ));
      }

      inline static constexpr bool contains ( const GeometryType &type ) noexcept
      {
        return (type.dim() <= maxdim);
      }

      inline static GeometryType type ( std::size_t index ) noexcept
      {
        return (index < RegularTypeIndexSet::size() ? RegularTypeIndexSet::type( index ) : GeometryType( 0, static_cast< int >( index - RegularTypeIndexSet::size() ), true ));
      }
    };



    // SingleGeometryTypeIndexSet
    // --------------------------

    /** \class SingleGeometryTypeIndexSet
     *
     *  \brief Please doc me.
     *
     *  \tparam  topologyId  topology id of the geometry type
     *  \tparam  dim         dimension
     */
    template< unsigned int topologyId, int dim >
    class SingleGeometryTypeIndexSet
    {
      typedef SingleGeometryTypeIndexSet< topologyId, dim > This;

    public:
      inline static constexpr std::size_t size () noexcept
      {
        return 1;
      }

      inline static constexpr std::size_t index ( const GeometryType &type ) noexcept
      {
        return 0;
      }

      inline static constexpr bool contains ( const GeometryType &type ) noexcept
      {
        return (type == GeometryType( topologyId, dim ));
      }

      inline static GeometryType type ( std::size_t index ) noexcept
      {
        return GeometryType( topologyId, dim );
      }
    };



    // GeometryTypes
    // -------------

    template< class TypeIndexSet >
    struct GeometryTypes
    {
      struct Iterator
      {
        bool operator== ( const Iterator &other ) const { return (index == other.index); }
        bool operator!= ( const Iterator &other ) const { return (index != other.index); }

        GeometryType operator* () const { return TypeIndexSet::type( index ); }

        Iterator &operator++ () { ++index; return *this; }

        std::size_t index;
      };

      Iterator begin () const { return Iterator{ 0 }; }
      Iterator end () const { return Iterator{ TypeIndexSet::size() }; }
    };

} // namespace hpDG
} // namespace Fem
} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_TYPEINDEXSET_HH
