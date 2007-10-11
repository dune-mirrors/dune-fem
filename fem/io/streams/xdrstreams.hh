#ifndef DUNE_FEM_XDRSTREAMS_HH
#define DUNE_FEM_XDRSTREAMS_HH

#include <rpc/xdr.h>

#include <dune/common/exceptions.hh>

#include <dune/fem/io/streams/streams.hh>

namespace Dune
{

  template< class OutStreamImp >
  struct XDROutStreamTraits
  {
    typedef OutStreamImp OutStreamType;
  };


  
  template< class OutStreamImp >
  class XDRBasicOutStream
  : public OutStreamInterface< XDROutStreamTraits< OutStreamImp > >
  {
  public:
    typedef OutStreamImp OutStreamType;
    typedef XDROutStreamTraits< OutStreamType > Traits;

  private:
    typedef XDRBasicOutStream< OutStreamType > ThisType;
    typedef OutStreamInterface< Traits > BaseType;

    enum { maxStringSize = 2<<18 };
    
  protected:
    XDR xdrs_;
    bool valid_;

  protected:
    inline XDRBasicOutStream ()
    : valid_( true )
    {
    }

  public:
    inline bool valid () const
    {
      return valid_;
    }
 
    inline void writeDouble ( double value )
    {
      valid_ &= (xdr_double( xdrs(), &value ) != 0);
    }

    inline void writeFloat ( float value )
    {
      valid_ &= (xdr_float( xdrs(), &value ) != 0);
    }

    inline void writeInt ( int value )
    {
      valid_ &= (xdr_int( xdrs(), &value ) != 0);
    }

    inline void writeString ( const std :: string &s )
    {
      assert( s.size() < maxStringSize );
      char *cs = s.c_str();
      valid_ &= (xdr_string( xdrs(), &cs, maxStringSize ) != 0);
    }

    inline void writeUnsignedInt ( unsigned int value )
    {
      valid_ &= (xdr_u_int( xdrs(), &value ) != 0);
    }
    
  protected:
    inline XDR *xdrs ()
    {
      return &xdrs_;
    }
  };


  
  template< class InStreamImp >
  struct XDRInStreamTraits
  {
    typedef InStreamImp InStreamType;
  };



  template< class InStreamImp >
  class XDRBasicInStream
  : public InStreamInterface< XDRInStreamTraits< InStreamImp > >
  {
  public:
    typedef InStreamImp InStreamType;
    typedef XDRInStreamTraits< InStreamType > Traits;

  private:
    typedef XDRBasicInStream< InStreamType > ThisType;
    typedef InStreamInterface< Traits > BaseType;

    enum { maxStringSize = 2<<18 };

  protected:
    XDR xdrs_;
    bool valid_;

  protected:
    inline XDRBasicInStream ()
    : valid_( true )
    {
    }

  public:
    inline bool valid () const
    {
      return valid_;
    }
    
    inline void readDouble ( double &value )
    {
      valid_ &= (xdr_double( xdrs(), &value ) != 0);
    }

    inline void readFloat ( float &value )
    {
      valid_ &= (xdr_float( xdrs(), &value ) != 0);
    }

    inline void readInt ( int &value )
    {
      valid_ &= (xdr_int( xdrs(), &value ) != 0);
    }

    inline void readString ( std :: string &s )
    {
      char data[ maxStringSize ];
      char *cs = &(data[ 0 ]);
      xdr_string( xdrs(), &cs, maxStringSize );
      s = data;
    }

    inline void readUnsignedInt ( unsigned int &value )
    {
      valid_ &= (xdr_u_int( xdrs(), &value ) != 0);
    }
    
  protected:
    inline XDR *xdrs ()
    {
      return &xdrs_;
    }
  };



  class XDRFileOutStream
  : public XDRBasicOutStream< XDRFileOutStream >
  {
  private:
    typedef XDRFileOutStream ThisType;
    typedef XDRBasicOutStream< ThisType > BaseType;

  protected:
    using BaseType :: xdrs;

  protected:
    FILE *file_;

  public:
    inline XDRFileOutStream ( const std :: string filename )
    {
      file_ = fopen( filename.c_str(), "wb" );
      if( file_ == 0 )
        DUNE_THROW( IOError, "XDRFileOutStream: Unable to open file." );
      xdrstdio_create( xdrs(), file_, XDR_ENCODE );
    }

    inline ~XDRFileOutStream ()
    {
      xdr_destroy( xdrs() );
      fclose( file_ );
    }
    
    inline void flush ()
    {
      fflush( file_ );
    }
  };


  
  class XDRFileInStream
  : public XDRBasicInStream< XDRFileInStream >
  {
  private:
    typedef XDRFileInStream ThisType;
    typedef XDRBasicInStream< ThisType > BaseType;

  protected:
    using BaseType :: xdrs;

  protected:
    FILE *file_;

  public:
    inline XDRFileInStream ( const std :: string filename )
    {
      file_ = fopen( filename.c_str(), "rb" );
      if( file_ == 0 )
        DUNE_THROW( IOError, "XDRFileInStream: Unable to open file." );
      xdrstdio_create( xdrs(), file_, XDR_DECODE );
    }

    inline ~XDRFileInStream ()
    {
      xdr_destroy( xdrs() );
      fclose( file_ );
    }
  };

}

#endif
