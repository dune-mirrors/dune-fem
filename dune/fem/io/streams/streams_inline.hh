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
                   const std::conditional< std::is_same<unsigned long, uint64_t>::value,
                        unsigned long long, // select long long in case long and uint64 are the same
                        unsigned long >::type& value
                 )
    {
      assert( sizeof(value) <= sizeof(uint64_t) );
      uint64_t v = value;
      out.writeUnsignedInt64( v );
      return out;
    }

    template< class Traits >
    inline OutStreamInterface< Traits > &
      operator<< ( OutStreamInterface< Traits > &out,
                   const int64_t value )
    {
      out.writeSignedInt64( value );
      return out;
    }

    template< class Traits >
    inline OutStreamInterface< Traits > &
      operator<< ( OutStreamInterface< Traits > &out,
                   const std::conditional< std::is_same<long, int64_t>::value,
                        long long, // select long long in case long and int64 are the same
                        long >::type& value )
    {
      assert( sizeof(value) <= sizeof(int64_t));
      int64_t v = value;
      out.writeSignedInt64( v );
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
                   std::conditional< std::is_same<unsigned long, uint64_t>::value,
                        unsigned long long,
                        unsigned long >::type& value )
    {
      // in any case, convert to uint64_t
      assert( sizeof(value) <= sizeof(uint64_t));
      uint64_t v;
      in.readUnsignedInt64( v );
      value = v;
      return in;
    }

    template< class Traits >
    inline InStreamInterface< Traits > &
      operator>> ( InStreamInterface< Traits > &in,
                   int64_t &value )
    {
      in.readSignedInt64( value );
      return in;
    }

    template< class Traits >
    inline InStreamInterface< Traits > &
      operator>> ( InStreamInterface< Traits > &in,
                   std::conditional< std::is_same<long, int64_t>::value,
                        long long,
                        long >::type& value )
    {
      // in any case, convert to uint64_t
      assert( sizeof(value) <= sizeof(int64_t));
      int64_t v;
      in.readSignedInt64( v );
      value = v;
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
