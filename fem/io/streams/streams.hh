#ifndef DUNE_FEM_STREAMS_HH
#define DUNE_FEM_STREAMS_HH

#include <string>

#include <dune/common/bartonnackmanifcheck.hh>

namespace Dune
{

  template< class TraitsImp >
  class OutStreamInterface
  {
  public:
    typedef TraitsImp Traits;

  private:
    typedef OutStreamInterface< Traits > ThisType;

  public:
    typedef ThisType OutStreamInterfaceType;

    typedef typename Traits :: OutStreamType OutStreamType;

  public:
    inline bool operator! () const
    {
      return !valid();
    }

    inline void flush ()
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().flush() );
    }
    
    inline bool valid () const
    {
      CHECK_INTERFACE_IMPLEMENTATON( asImp().valid() );
      return asImp().valid();
    }
   
    inline void writeDouble ( double value )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().writeDouble( value ) );
    }
    
    inline void writeFloat ( float value )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().writeFloat( value ) );
    }
    
    inline void writeInt ( int value )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().writeInt( value ) );
    }

    inline void writeString ( const std :: string &s )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().writeString( s ) );
    }

    inline void writeUnsignedInt ( unsigned int value )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().writeUnsignedInt( value ) );
    }

  protected:
    // Barton-Nackman trick
    inline const OutStreamType &asImp () const
    {
      return static_cast< const OutStreamType & >( *this );
    }

    // Barton-Nackman trick
    inline OutStreamType &asImp ()
    {
      return static_cast< OutStreamType & >( *this );
    }
  };


  
  template< class TraitsImp >
  class InStreamInterface
  {
  public:
    typedef TraitsImp Traits;

  private:
    typedef InStreamInterface< Traits > ThisType;

  public:
    typedef ThisType InStreamInterfaceType;

    typedef typename Traits :: InStreamType InStreamType;

  public:
    inline bool operator! () const
    {
      return !valid();
    }

    inline bool valid () const
    {
      CHECK_INTERFACE_IMPLEMENTATON( asImp().valid() );
      return asImp().valid();
    }
   
    inline void readDouble ( double &value )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().readDouble( value ) );
    }
    
    inline void readFloat ( float &value )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().readFloat( value ) );
    }
    
    inline void readInt ( int &value )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().readInt( value ) );
    }

    inline void readString ( std :: string &s )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().readString( s ) );
    }

    inline void readUnsignedInt ( unsigned int &value )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().readUnsignedInt( value ) );
    }

  protected:
    // Barton-Nackman trick
    inline const InStreamType &asImp () const
    {
      return static_cast< const InStreamType & >( *this );
    }

    // Barton-Nackman trick
    inline InStreamType &asImp ()
    {
      return static_cast< InStreamType & >( *this );
    }
  };
 
}

#include "streams_inline.hh"

#endif
