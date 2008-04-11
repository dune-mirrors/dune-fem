#ifndef DUNE_FEM_XDRSTREAMS_HH
#define DUNE_FEM_XDRSTREAMS_HH

#include <cassert>
#include <rpc/xdr.h>

#include <dune/fem/io/streams/streams.hh>

namespace Dune
{

  template< class OutStreamImp >
  struct XDROutStreamTraits
  {
    typedef OutStreamImp OutStreamType;
  };


 
  /** \class XDRBasicOutStream
   *  \ingroup InOutStreams
   *  \brief base implementation for XDR output streams
   *  
   *  This class implements the writing functions for an XDR stream. It must
   *  be associated to a stream by a child class.
   *
   *  The following XDR output streams have been implemented:
   *  -XDRFileOutStream
   */
  template< class OutStreamImp >
  class XDRBasicOutStream
  : public OutStreamInterface< XDROutStreamTraits< OutStreamImp > >
  {
  public:
    //! type of the implementaton (Barton-Nackman)
    typedef OutStreamImp OutStreamType;
    
    //! type of the traits
    typedef XDROutStreamTraits< OutStreamType > Traits;

  private:
    typedef XDRBasicOutStream< OutStreamType > ThisType;
    typedef OutStreamInterface< Traits > BaseType;

    enum { maxStringSize = 2<<18 };
    
  protected:
    XDR xdrs_;

  protected:
    using BaseType :: writeError;

  protected:
    inline XDRBasicOutStream ()
    {
    }

  public:
    /** \copydoc Dune::OutStreamInterface::writeDouble */
    inline void writeDouble ( double value )
    {
      if( xdr_double( xdrs(), &value ) == 0 )
        writeError();
    }

    /** \copydoc Dune::OutStreamInterface::writeFloat */
    inline void writeFloat ( float value )
    {
      if( xdr_float( xdrs(), &value ) == 0 )
        writeError();
    }

    /** \copydoc Dune::OutStreamInterface::writeInt */
    inline void writeInt ( int value )
    {
      if( xdr_int( xdrs(), &value ) == 0 )
        writeError();
    }

    /** \copydoc Dune::OutStreamInterface::writeString */
    inline void writeString ( const std :: string &s )
    {
      assert( s.size() < maxStringSize );
      const char *cs = s.c_str();
      if( xdr_string( xdrs(), (char **)&cs, maxStringSize ) == 0 )
        writeError();
    }

    /** \copydoc Dune::OutStreamInterface::writeUnsignedInt */
    inline void writeUnsignedInt ( unsigned int value )
    {
      if( xdr_u_int( xdrs(), &value ) == 0 )
        writeError();
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


  /** \class XDRBasicInStream
   *  \ingroup InOutStreams
   *  \brief base implementation for XDR input streams
   *  
   *  This class implements the reading functions for an XDR stream. It must
   *  be associated to a stream by a child class.
   *
   *  The following XDR input streams have been implemented:
   *  -XDRFileInStream
   */
  template< class InStreamImp >
  class XDRBasicInStream
  : public InStreamInterface< XDRInStreamTraits< InStreamImp > >
  {
  public:
    //! type of the implementation (Barton-Nackman)
    typedef InStreamImp InStreamType;
    
    //! type of the traits
    typedef XDRInStreamTraits< InStreamType > Traits;

  private:
    typedef XDRBasicInStream< InStreamType > ThisType;
    typedef InStreamInterface< Traits > BaseType;

    enum { maxStringSize = 2<<18 };

  protected:
    XDR xdrs_;

  protected:
    using BaseType :: readError;

  protected:
    inline XDRBasicInStream ()
    {
    }

  public:
    /** \copydoc Dune::InStreamInterface::readDouble */
    inline void readDouble ( double &value )
    {
      if( xdr_double( xdrs(), &value ) == 0 )
        readError();
    }

    /** \copydoc Dune::InStreamInterface::readFloat */
    inline void readFloat ( float &value )
    {
      if( xdr_float( xdrs(), &value ) == 0 )
        readError();
    }

    /** \copydoc Dune::InStreamInterface::readInt */
    inline void readInt ( int &value )
    {
      if( xdr_int( xdrs(), &value ) == 0 )
        readError();
    }

    /** \copydoc Dune::InStreamInterface::readString */
    inline void readString ( std :: string &s )
    {
      char data[ maxStringSize ];
      char *cs = &(data[ 0 ]);
      if( xdr_string( xdrs(), &cs, maxStringSize ) == 0 )
        readError();
      s = data;
    }

    /** \copydoc Dune::InStreamInterface::readUnsignedInt */
    inline void readUnsignedInt ( unsigned int &value )
    {
      if( xdr_u_int( xdrs(), &value ) == 0 )
        readError();
    }
    
  protected:
    inline XDR *xdrs ()
    {
      return &xdrs_;
    }
  };



  /** \class XDRFileOutStream
   *  \ingroup InOutStreams
   *  \brief XDR output stream writing into a file
   *
   *  \newimplementation
   */
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
    /** \brief constructor
     *
     *  \param[in]  filename  name of the file to write to
     */
    inline explicit XDRFileOutStream ( const std :: string filename )
    {
      file_ = fopen( filename.c_str(), "wb" );
      if( file_ == 0 )
        DUNE_THROW( IOError, "XDRFileOutStream: Unable to open file." );
      xdrstdio_create( xdrs(), file_, XDR_ENCODE );
    }

    /** \brief destructor */
    inline ~XDRFileOutStream ()
    {
      xdr_destroy( xdrs() );
      fclose( file_ );
    }
    
    /** \copydoc Dune::OutStreamInterface::flush */
    inline void flush ()
    {
      fflush( file_ );
    }
  };


  
  /** \class XDRFileInStream
   *  \ingroup InOutStreams
   *  \brief XDR output stream reading from a file
   *
   *  \newimplementation
   */
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
    /** \brief constructor
     *
     *  \param[in]  filename  name of the file to read from
     */
    inline explicit XDRFileInStream ( const std :: string filename )
    {
      file_ = fopen( filename.c_str(), "rb" );
      if( file_ == 0 )
        DUNE_THROW( IOError, "XDRFileInStream: Unable to open file." );
      xdrstdio_create( xdrs(), file_, XDR_DECODE );
    }

    /** \brief destructor
     */
    inline ~XDRFileInStream ()
    {
      xdr_destroy( xdrs() );
      fclose( file_ );
    }
  };

}

#endif
