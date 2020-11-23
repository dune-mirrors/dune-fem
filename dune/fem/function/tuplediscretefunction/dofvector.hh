#ifndef DUNE_FEM_FUNCTION_BLOCKVECTORS_TUPLE_HH
#define DUNE_FEM_FUNCTION_BLOCKVECTORS_TUPLE_HH

#include <limits>
#include <utility>
#include <tuple>

#include <dune/common/hybridutilities.hh>

#include <dune/fem/common/forloop.hh>
#include <dune/fem/common/utility.hh>
#include <dune/fem/function/blockvectors/defaultblockvectors.hh>
#include <dune/fem/function/common/functor.hh>


namespace Dune
{

  namespace Fem
  {


    // TupleDofVector
    // --------------

    template< class ... DofVectors >
    class TupleDofVector
      : public IsBlockVector,
        public std::tuple< DofVectors& ... >
    {
      typedef TupleDofVector< DofVectors ... > ThisType;
      typedef std::tuple< DofVectors & ... > BaseType;

      static_assert( sizeof ... ( DofVectors ) > 0, "TupleDofVector needs at least one DofVector." );

      typedef std::tuple< DofVectors ... > DofVectorTuple;

      typedef decltype ( std::index_sequence_for< DofVectors ... >() ) Sequence;

    public:
      static_assert( Std::are_all_same< typename DofVectors::FieldType ... >::value, "All blocks need to have the same FieldType." );
      typedef typename std::tuple_element< 0, DofVectorTuple >::type::FieldType FieldType;

      struct Iterator;
      struct ConstIterator;

      typedef Iterator IteratorType;
      typedef ConstIterator ConstIteratorType;

      typedef std::size_t SizeType;
      typedef FieldType value_type;
      typedef SizeType size_type;

      static const int blockSize = 1;

      typedef FieldType *DofBlockType;
      typedef const FieldType *ConstDofBlockType;

      typedef Fem::Envelope< DofBlockType > DofBlockPtrType;
      typedef Fem::Envelope< ConstDofBlockType > ConstDofBlockPtrType;

      TupleDofVector ( DofVectors & ... dofVectors )
        : BaseType( std::forward_as_tuple( dofVectors ... ) )
      {}

      // constructor
      TupleDofVector ( BaseType data ) : BaseType( data ) {}
      TupleDofVector ( const ThisType & ) = default;
      TupleDofVector ( ThisType && ) = default;

      const ThisType &operator= ( const ThisType &other )
      {
        assign( other, Sequence() );
        return *this;
      }

      const ThisType &operator+= ( const ThisType &other )
      {
        axpy( 1.0, other, Sequence() );
        return *this;
      }

      const ThisType &operator-= ( const ThisType &other )
      {
        axpy( -1.0, other, Sequence() );
        return *this;
      }

    private:
      template <int i>
      struct Acc
      {
        static void apply(const ThisType &a, const ThisType& b, FieldType&  result )
        {
          result += (a.template subDofVector< i >() * b.template subDofVector<i>());
        }
      };

    public:
      FieldType operator* ( const ThisType &other ) const
      {
        FieldType result( 0 );
        static const int length = std::tuple_size< BaseType >::value;
        Dune::Fem::ForLoop<Acc, 0, length-1>::apply(*this, other, result);
        return result;
      }

      const ThisType &operator*= ( const FieldType &scalar )
      {
        scale( scalar, Sequence() );
        return *this;
      }

      const ThisType &operator/= ( const FieldType &scalar )
      {
        scale( 1.0 / scalar, Sequence() );
        return *this;
      }

      void axpy ( const FieldType &scalar, const ThisType &other )
      {
        axpy( scalar, other, Sequence() );
      }

      void clear () { clear( Sequence() ); }

      SizeType size () const { return size( Sequence() ); }

      // ----------------------------------------------------------------------------

      IteratorType begin () { return IteratorType( *this, 0 ); }
      ConstIteratorType begin () const { return ConstIteratorType( *this, 0 ); }

      IteratorType end () { return IteratorType( *this, size() ); }
      ConstIteratorType end () const { return ConstIteratorType( *this, size() ); }

      DofBlockType operator[] ( std::size_t index )
      {
        return blockAccess( index, std::integral_constant< std::size_t, 0 >() );
      }
      ConstDofBlockType operator[] ( std::size_t index ) const
      {
        return blockAccess( index, std::integral_constant< std::size_t, 0 >() );
      }

      DofBlockPtrType blockPtr ( std::size_t index )
      {
        return DofBlockPtrType( this->operator[]( index ) );
      }
      ConstDofBlockType blockPtr ( std::size_t index ) const
      {
        return ConstDofBlockPtrType( this->operator[]( index ) );
      }

      // ----------------------------------------------------------------------------

      void reserve ( SizeType size ) {}
      void resize ( SizeType size ) {}

      constexpr std::size_t blocks () const { return sizeof ... ( DofVectors ); }

      template< int i >
      const typename std::tuple_element< i, DofVectorTuple >::type &
      subDofVector() const { return std::get< i >( *this ); }

      template< int i >
      typename std::tuple_element< i, DofVectorTuple >::type &
      subDofVector() {
        return std::get< i >( *this );
      }

    protected:
      template< std::size_t ... i >
      void scale ( FieldType scale, std::index_sequence< i ... > )
      {
        std::ignore = std::make_tuple( (std::get< i >( *this ) *= scale, i ) ... );
      }

      template< std::size_t ... i >
      void axpy ( FieldType a, const ThisType &other, std::index_sequence< i ... > )
      {
        std::ignore = std::make_tuple( ( std::get< i >( *this ).axpy( a, std::get< i >( other ) ), i ) ... );
      }

      template< std::size_t ... i >
      void assign ( const ThisType &other, std::index_sequence< i ... > )
      {
        std::ignore = std::make_tuple( ( std::get< i >( *this ) = std::get< i >( other ), i ) ... );
      }

      template< std::size_t ... i >
      SizeType size ( std::index_sequence< i ... > ) const
      {
        return Std::sum( std::get< i >( *this ).size() *
                         std::tuple_element< i, DofVectorTuple >::type::blockSize ... );
      }

      template< std::size_t ... I >
      void clear ( std::index_sequence< I ... > )
      {
        std::ignore = std::make_tuple( ( std::get< I >( *this ).clear(), I ) ... );
      }

      template< std::size_t i >
      FieldType *blockAccess ( std::size_t index, std::integral_constant< std::size_t, i > )
      {
        const std::size_t thisBlockSize = std::tuple_element< i, DofVectorTuple >::type::blockSize;
        std::size_t offset = std::get< i >( *this ).size() * thisBlockSize;
        if( index < offset )
          return &std::get< i >( *this )[ index / thisBlockSize ][ index % thisBlockSize ];
        else
          return blockAccess( index - offset, std::integral_constant< std::size_t, i +1 >() );
      }

      FieldType *blockAccess ( std::size_t index, std::integral_constant< std::size_t, sizeof ... ( DofVectors ) > )
      {
        DUNE_THROW( RangeError, "Index out of range" );
      }

      template< std::size_t i >
      const FieldType *blockAccess ( std::size_t index, std::integral_constant< std::size_t, i > ) const
      {
        const std::size_t thisBlockSize = std::tuple_element< i, DofVectorTuple >::type::blockSize;
        std::size_t offset = std::get< i >( *this ).size() * thisBlockSize;
        if( index < offset )
          return &std::get< i >( *this )[ index / thisBlockSize ][ index % thisBlockSize ];
        else
          return blockAccess( index - offset, std::integral_constant< std::size_t, i +1 >() );
      }

      const FieldType *blockAccess ( std::size_t index, std::integral_constant< std::size_t, sizeof ... ( DofVectors ) > ) const
      {
        DUNE_THROW( RangeError, "Index out of range" );
      }

    };


    template< class ... DofVectors >
    struct TupleDofVector< DofVectors ... >::
    Iterator
    {
      Iterator( TupleDofVector< DofVectors ... > &container, std::size_t it = 0 )
        : container_( container ), iterator_( it ) {}

      Iterator( const Iterator & ) = default;
      Iterator &operator= ( const Iterator & ) = default;

      FieldType &operator* () { return *container_[ iterator_ ]; }
      FieldType *operator-> () { return container_[ iterator_ ]; }

      Iterator &operator++ () { iterator_++; return *this; }
      Iterator operator++ ( int )
      {
        ThisType copy( *this );
        iterator_++;
        return copy;
      }

      bool operator== ( const Iterator &other ) const { return iterator_ == other.iterator_; }
      bool operator!= ( const Iterator &other ) const { return iterator_ != other.iterator_; }

    protected:
      TupleDofVector< DofVectors ... > &container_;
      std::size_t iterator_;
    };

    template< class ... DofVectors >
    struct TupleDofVector< DofVectors ... >::
    ConstIterator
    {
      ConstIterator( const TupleDofVector< DofVectors ... > &container, std::size_t it = 0 )
        : container_( container ), iterator_( it ) {}

      ConstIterator( const ConstIterator & ) = default;

      const FieldType &operator* () const { return *container_[ iterator_ ]; }
      FieldType *operator-> () const { return container_[ iterator_ ]; }

      ConstIterator &operator++ () { iterator_++; return *this; }
      ConstIterator operator++ ( int )
      {
        ThisType copy( *this );
        iterator_++;
        return copy;
      }

      bool operator== ( const ConstIterator &other ) const { return iterator_ == other.iterator_; }
      bool operator!= ( const ConstIterator &other ) const { return iterator_ != other.iterator_; }

    protected:
      const TupleDofVector< DofVectors ... > &container_;
      std::size_t iterator_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_BLOCKVECTORS_TUPLE_HH
