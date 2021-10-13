#ifndef DUNE_FEM_GRIDPART_COMMON_SHAREDGEOMETRY_HH
#define DUNE_FEM_GRIDPART_COMMON_SHAREDGEOMETRY_HH

#include <tuple>
#include <type_traits>
#include <utility>

#include <dune/common/fmatrix.hh>

#include <dune/grid/common/geometry.hh>

namespace Dune
{

  // SharedGeometry
  // --------------

  template< class Impl, class Alloc = std::allocator< Impl > >
  class SharedGeometry
  {
    typedef SharedGeometry< Impl, Alloc > This;

    typedef std::pair< Impl, std::size_t > Data;

  public:
    typedef Impl Implementation;
    typedef Alloc Allocator;

    static const int mydimension = Implementation::mydimension;
    static const int coorddimension = Implementation::coorddimension;

    typedef typename Implementation::ctype ctype;

    typedef typename Implementation::LocalCoordinate LocalCoordinate;
    typedef typename Implementation::GlobalCoordinate GlobalCoordinate;

    typedef typename Implementation::JacobianTransposed JacobianTransposed;
    typedef typename Implementation::JacobianInverseTransposed JacobianInverseTransposed;

    SharedGeometry () = default;

    template< class... Args, std::enable_if_t< std::is_constructible< Impl, Args &&... >::value, int > = 0 >
    SharedGeometry ( Args &&... args )
      : data_( construct( std::forward_as_tuple< Args... >( args... ) ) )
    {}

    SharedGeometry ( const This &other )
      : data_( other.data_ ), allocator_( other.allocator_ )
    {
      if (data_)
        ++data_->second;
    }

    SharedGeometry ( This &&other )
      : data_( other.data_ ), allocator_( std::move( other.allocator_ ) )
    {
      other.data_ = nullptr;
    }

    ~SharedGeometry ()
    {
      if( data_ && (--data_->second == 0) )
        destroy();
    }

    This &operator= ( const This &other )
    {
      if( other.data_ )
        ++other.data_->second;
      if( data_ && (--data_->second == 0) )
        destroy();
      data_ = other.data_;
      allocator_ = other.allocator_;
      return *this;
    }

    This &operator= ( This &&other )
    {
      if( data_ && (--data_->second == 0) )
        destroy();
      data_ = other.data_;
      allocator_ = std::move( other.allocator_ );
      other.data_ = nullptr;
      return *this;
    }

    operator bool () const { return static_cast< bool >( data_ ); }

    bool affine () const { return impl().affine(); }
    GeometryType type () const { return impl().type(); }

    int corners () const { return impl().corners(); }
    GlobalCoordinate corner ( int i ) const { return impl().corner( i ); }
    GlobalCoordinate center () const { return impl().center(); }

    GlobalCoordinate global ( const LocalCoordinate &local ) const { return impl().global( local ); }
    LocalCoordinate local ( const GlobalCoordinate &global ) const { return impl().local( global ); }

    ctype integrationElement ( const LocalCoordinate &local ) const { return impl().integrationElement( local ); }
    ctype volume () const { return impl().volume(); }

    JacobianTransposed jacobianTransposed ( const LocalCoordinate &local ) const { return impl().jacobianTransposed( local ); }
    JacobianInverseTransposed jacobianInverseTransposed ( const LocalCoordinate &local ) const { return impl().jacobianInverseTransposed( local ); }

    Allocator allocator () const { return allocator_; }

    const Impl &impl () const { assert( data_ ); return data_->first; }
    Impl &impl () { assert( data_ ); return data_->first; }

  private:
    template< class... Args >
    Data *construct ( std::tuple< Args... > args )
    {
      Data *data = allocator_.allocate( 1 );
      return new (data) Data( std::piecewise_construct, args, std::make_tuple( 1u ) );
    }

    void destroy ()
    {
      data_->~Data();
      allocator_.deallocate( data_, 1 );
    }

    Data *data_ = nullptr;
    typename std::allocator_traits< Allocator >::template rebind_alloc< Data > allocator_;
  };

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_COMMON_SHAREDGEOMETRY_HH
