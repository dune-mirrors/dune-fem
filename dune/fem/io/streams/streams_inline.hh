#ifndef DUNE_FEM_STREAMS_INLINE_HH
#define DUNE_FEM_STREAMS_INLINE_HH

#include <vector>
#include <array>

#include <dune/common/fvector.hh>

#include "streams.hh"

namespace Dune
{

  namespace Fem
  {

    template< class Traits >
    inline OutStreamInterface< Traits > &
      operator<< ( OutStreamInterface< Traits > &out,
                   const double value )
    {
      out.writeDouble( value );
      return out;
    }

    template< class Traits >
    inline OutStreamInterface< Traits > &
      operator<< ( OutStreamInterface< Traits > &out,
                   const float value )
    {
      out.writeFloat( value );
      return out;
    }

    template< class Traits >
    inline OutStreamInterface< Traits > &
      operator<< ( OutStreamInterface< Traits > &out,
                   const int value )
    {
      out.writeInt( value );
      return out;
    }

    template< class Traits >
    inline OutStreamInterface< Traits > &
      operator<< ( OutStreamInterface< Traits > &out,
                   const char value )
    {
      out.writeChar( value );
      return out;
    }

    template< class Traits >
    inline OutStreamInterface< Traits > &
      operator<< ( OutStreamInterface< Traits > &out,
                   const bool value )
    {
      out.writeBool( value );
      return out;
    }

    template< class Traits >
    inline OutStreamInterface< Traits > &
      operator<< ( OutStreamInterface< Traits > &out,
                   const std :: string &s )
    {
      out.writeString( s );
      return out;
    }

    template< class Traits >
    inline OutStreamInterface< Traits > &
      operator<< ( OutStreamInterface< Traits > &out,
                   const unsigned int value )
    {
      out.writeUnsignedInt( value );
      return out;
    }

    template< class Traits, class T >
    inline OutStreamInterface< Traits > &
      operator<< ( OutStreamInterface< Traits > &out,
                   const std::complex<T> value )
    {
      out.writeDouble( std::real(value) );
      out.writeDouble( std::imag(value) );
      return out;
    }

    template <class ulongint, class uint64>
    struct SelectUnsignedLongInteger
    {
      // select uint64_t int
      typedef ulongint UnsignedLongIntType;

      template < class Traits >
      static void write( OutStreamInterface< Traits > &out,
                         const UnsignedLongIntType& value )
      {
        // in case uint64_t int and uint64_t are not the same
        // convert long to uint64_t, there will be no information loss
        assert( sizeof(ulongint) <= sizeof(uint64) );
        uint64 value64 = value ;
        out.writeUnsignedInt64( value64 );
      }

      template < class Traits >
      static void read( InStreamInterface< Traits > &in,
                        UnsignedLongIntType& value )
      {
        assert( sizeof(ulongint) <= sizeof(uint64) );
        // always read uint64_t int as uin64_t, since it is always written this way
        uint64 value64;
        in.readUnsignedInt64( value64 );
        value = value64;
      }
    };

    //- in case uint64_t int and uint64_t are the same, do nothing
    template <class ulongint>
    struct SelectUnsignedLongInteger< ulongint, ulongint >
    {
      struct UnsignedLongIntType {};
      template < class Traits >
      static void write( OutStreamInterface< Traits > &out,
                         const UnsignedLongIntType value )
      {
        DUNE_THROW(NotImplemented,"method not implemented");
      }

      template < class Traits >
      static void read( InStreamInterface< Traits > &in,
                        UnsignedLongIntType& value )
      {
        DUNE_THROW(NotImplemented,"method not implemented");
      }
    };

    template< class Traits >
    inline OutStreamInterface< Traits > &
      operator<< ( OutStreamInterface< Traits > &out,
                   const uint64_t value )
    {
      out.writeUnsignedInt64( value );
      return out;
    }

    template< class Traits >
    inline OutStreamInterface< Traits > &
      operator<< ( OutStreamInterface< Traits > &out,
                   const typename SelectUnsignedLongInteger<unsigned long, uint64_t>::UnsignedLongIntType& value )
    {
      SelectUnsignedLongInteger<unsigned long, uint64_t>::write( out, value );
      return out;
    }

    template< class Traits, class T, std::size_t N >
    inline OutStreamInterface< Traits > &
    operator<< ( OutStreamInterface< Traits > &out, const std::array< T, N > &value )
    {
      for( std::size_t i = 0; i < N; ++i )
        out << value[ i ];
      return out;
    }

    template< class Traits, class T, int N >
    inline OutStreamInterface< Traits > &
    operator<< ( OutStreamInterface< Traits > &out, const Dune::FieldVector< T, N > &value )
    {
      for( int i = 0; i < N; ++i )
        out << value[ i ];
      return out;
    }

    template< class Traits, class T, class A >
    inline OutStreamInterface< Traits > &
      operator<< ( OutStreamInterface< Traits > &out,
                   const std::vector< T, A > & value )
    {
      const size_t size = value.size();
      out << size;
      for( size_t i = 0; i < size; ++i )
        out << value[ i ];
      return out;
    }

    template< class Traits >
    inline InStreamInterface< Traits > &
      operator>> ( InStreamInterface< Traits > &in,
                   double &value )
    {
      in.readDouble( value );
      return in;
    }

    template< class Traits >
    inline InStreamInterface< Traits > &
      operator>> ( InStreamInterface< Traits > &in,
                   float &value )
    {
      in.readFloat( value );
      return in;
    }

    template< class Traits >
    inline InStreamInterface< Traits > &
      operator>> ( InStreamInterface< Traits > &in,
                   int &value )
    {
      in.readInt( value );
      return in;
    }

    template< class Traits >
    inline InStreamInterface< Traits > &
      operator>> ( InStreamInterface< Traits > &in,
                   char &value )
    {
      in.readChar( value );
      return in;
    }

    template< class Traits >
    inline InStreamInterface< Traits > &
      operator>> ( InStreamInterface< Traits > &in,
                   bool &value )
    {
      in.readBool( value );
      return in;
    }

    template< class Traits >
    inline InStreamInterface< Traits > &
      operator>> ( InStreamInterface< Traits > &in,
                   std :: string &s )
    {
      in.readString( s );
      return in;
    }

    template< class Traits >
    inline InStreamInterface< Traits > &
      operator>> ( InStreamInterface< Traits > &in,
                   unsigned int &value )
    {
      in.readUnsignedInt( value );
      return in;
    }

    template< class Traits >
    inline InStreamInterface< Traits > &
      operator>> ( InStreamInterface< Traits > &in,
                   uint64_t &value )
    {
      in.readUnsignedInt64( value );
      return in;
    }

    template< class Traits >
    inline InStreamInterface< Traits > &
      operator>> ( InStreamInterface< Traits > &in,
                   typename SelectUnsignedLongInteger<unsigned long, uint64_t>::UnsignedLongIntType& value )
    {
      SelectUnsignedLongInteger<unsigned long, uint64_t>::read( in, value );
      return in;
    }

    template< class Traits, class T, std::size_t N >
    inline InStreamInterface< Traits > &
    operator>> ( InStreamInterface< Traits > &in, std::array< T, N > &value )
    {
      for( std::size_t i = 0; i < N; ++i )
        in >> value[ i ];
      return in;
    }

    template< class Traits, class T, int N >
    inline InStreamInterface< Traits > &
    operator>> ( InStreamInterface< Traits > &in, Dune::FieldVector< T, N > &value )
    {
      for( int i = 0; i < N; ++i )
        in >> value[ i ];
      return in;
    }

    template< class Traits, class T >
    inline InStreamInterface< Traits > &
      operator>> ( InStreamInterface< Traits > &in,
                   std::complex<T> &value )
    {
      T r,i;
      in.readDouble( r );
      in.readDouble( i );
      value = std::complex<T>(r,i);
      return in;
    }

    template< class Traits, class T, class A >
    inline InStreamInterface< Traits > &
      operator>> ( InStreamInterface< Traits > &in,
                   std::vector< T, A > & value )
    {
      size_t size = 0;
      in >> size;
      value.resize( size );
      for( size_t i = 0; i < size; ++i )
        in >> value[ i ];
      return in;
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_STREAMS_INLINE_HH
