#ifndef DUNE_FEM_STREAMS_INLINE_HH
#define DUNE_FEM_STREAMS_INLINE_HH

#include <vector>

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

  template <class ulongint, class uint64, int id >  
  struct SelectUnsignedLongInteger 
  {
    // select uint64_t int 
    typedef ulongint UnsignedLongIntType;

    template < class Traits >
    static void write( OutStreamInterface< Traits > &out,
                       const UnsignedLongIntType& value ) 
    {
      if( sizeof( uint64 ) == 8 ) 
      {
        // in case uint64_t int and uint64_t are not the same 
        // convert long to uint64_t, there will be no information loss
        assert( sizeof(ulongint) <= sizeof(uint64) );
        uint64 value64 = value ;
        out.writeUnsignedInt64( value64 );
      }
      else 
      {
        // in case uint64_t int and uint64_t are not the same 
        // convert long to uint64_t, there will be no information loss
        assert( sizeof(uint64) == 4 );
        unsigned int value32 = value ;
        out.writeUnsignedInt( value32 );
      }
    }

    template < class Traits >
    static void read( InStreamInterface< Traits > &in,
                      UnsignedLongIntType& value ) 
    {
      if( sizeof( uint64 ) == 8 ) 
      {
        assert( sizeof(ulongint) <= sizeof(uint64) );
        // always read uint64_t int as uin64_t, since it is always written this way 
        uint64 value64; 
        in.readUnsignedInt64( value64 );
        value = value64;
      }
      else
      {
        assert( sizeof(uint64) == 4 );
        // always read uint64_t int as uin64_t, since it is always written this way 
        unsigned int value32; 
        in.readUnsignedInt( value32 );
        value = value32;
      }
    }
  };

  //- in case uint64_t int and uint64_t are the same, do nothing 
  template <class ulongint, int id >  
  struct SelectUnsignedLongInteger< ulongint, ulongint, id >
  {
    template <int i> 
    struct DefUnsignedLongIntType {};

    typedef DefUnsignedLongIntType< id >  UnsignedLongIntType;
    
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
                 const typename SelectUnsignedLongInteger<unsigned long, uint64_t, 0>::UnsignedLongIntType& value )
  {
    SelectUnsignedLongInteger<unsigned long, uint64_t, 0>::write( out, value );
    return out;
  }

  template< class Traits >
  inline OutStreamInterface< Traits > &
    operator<< ( OutStreamInterface< Traits > &out,
                 const typename SelectUnsignedLongInteger<size_t, unsigned int, 1>::UnsignedLongIntType& value )
  {
    SelectUnsignedLongInteger<size_t, unsigned int, 1>::write( out, value );
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
                 typename SelectUnsignedLongInteger<unsigned long, uint64_t, 0>::UnsignedLongIntType& value )
  {
    SelectUnsignedLongInteger<unsigned long, uint64_t, 0>::read( in, value );
    return in;
  }

  template< class Traits >
  inline InStreamInterface< Traits > &
    operator>> ( InStreamInterface< Traits > &in,
                 typename SelectUnsignedLongInteger<size_t, unsigned int, 1>::UnsignedLongIntType& value )
  {
    SelectUnsignedLongInteger<size_t, unsigned int, 1>::read( in, value );
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

} // end Namespace Fem   

} // end namespace Dune 

#endif // #ifndef DUNE_FEM_STREAMS_INLINE_HH
