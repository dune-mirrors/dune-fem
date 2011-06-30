#ifndef DUNE_FEM_BINARYSTREAMS_HH
#define DUNE_FEM_BINARYSTREAMS_HH

#include <fstream>

#include <dune/fem/io/streams/streams.hh>

namespace Dune
{

  class BinaryFileOutStream;
  class BinaryFileInStream;

  
  
  struct BinaryFileOutStreamTraits
  {
    typedef BinaryFileOutStream OutStreamType;
  };

  
  /** \class BinaryFileOutStream
   *  \ingroup InOutStreams
   *  \brief output stream writing into a file in binary form
   *
   *  \note This stream directly stores the binary representation of the data.
   *        The binary representation might differ between different machines
   *        (e.g., little endian vs. big endian).
   *
   *  \newimplementation
   */
  class BinaryFileOutStream
  : public OutStreamInterface< BinaryFileOutStreamTraits >
  {
    typedef BinaryFileOutStream ThisType;
    typedef OutStreamInterface< BinaryFileOutStreamTraits > BaseType;

  public:
    //! type of the traits
    typedef BinaryFileOutStreamTraits Traits;
    
  protected:
    using BaseType::writeError;

  public:
    /** \brief constructor
     *
     *  \param[in]  filename  name of a file to write to
     */
    explicit BinaryFileOutStream ( const std::string &filename )
    : stream_( filename.c_str(), std::ios::binary )
    {}

    /** \copydoc Dune::OutStreamInterface::flush */
    void flush ()
    {
      stream_.flush();
    }

    /** \copydoc Dune::OutStreamInterface::writeDouble */
    void writeDouble ( const double value )
    {
      writePrimitive( value );
    }

    /** \copydoc Dune::OutStreamInterface::writeFloat */
    void writeFloat ( const float value )
    {
      writePrimitive( value );
    }

    /** \copydoc Dune::OutStreamInterface::writeInt */
    void writeInt ( const int value )
    {
      writePrimitive( value );
    }

    /** \copydoc Dune::OutStreamInterface::writeChar */
    void writeChar ( const char value )
    {
      writePrimitive( value );
    }

    /** \copydoc Dune::OutStreamInterface::writeChar */
    void writeBool ( const bool value )
    {
      writePrimitive( value );
    }

    /** \copydoc Dune::OutStreamInterface::writeString */
    void writeString ( const std::string &s )
    {
      const unsigned int length = s.length();
      writePrimitive( length );
      for( unsigned int i = 0; i < length; ++i )
        writePrimitive( s[ i ] );
    }

    /** \copydoc Dune::OutStreamInterface::writeUnsignedInt */
    void writeUnsignedInt ( unsigned int value )
    {
      writePrimitive( value );
    }

    /** \copydoc Dune::OutStreamInterface::writeUnsignedLong */
    void writeUnsignedLong ( unsigned long value )
    {
      writePrimitive( value );
    }

  protected:
    bool valid () const
    {
      return stream_.good() | stream_.eof();
    }

    template< class T >
    void writePrimitive ( const T &value )
    {
      union { T value; char bytes[ sizeof( T ) ]; } convert;
      convert.value = value;
      stream_.write( convert.bytes, sizeof( T ) );
      if( !valid () )
        writeError();
    }

  protected:
    std::ofstream stream_;
  };


  
  struct BinaryFileInStreamTraits
  {
    typedef BinaryFileInStream InStreamType;
  };

  

  /** \class BinaryFileInStream
   *  \ingroup InOutStreams
   *  \brief input stream reading from a file in binary form
   *
   *  \note This stream directly stores the binary representation of the data.
   *        The binary representation might differ between different machines
   *        (e.g., little endian vs. big endian).
   *
   *  \newimplementation
   */
  class BinaryFileInStream
  : public InStreamInterface< BinaryFileInStreamTraits >
  {
    typedef BinaryFileInStream ThisType;
    typedef InStreamInterface< BinaryFileInStreamTraits > BaseType;

  public:
    //! type of the traits
    typedef BinaryFileInStreamTraits Traits;
    
  protected:
    using BaseType::readError;

  public:
    /** \brief constructor
     *
     *  \param[in]  filename  name of a file to write to
     */
    BinaryFileInStream ( const std::string &filename )
    : stream_( filename.c_str(), std::ios::binary )
    {}

    /** \copydoc Dune::InStreamInterface::readDouble */
    void readDouble ( double &value )
    {
      readPrimitive( value );
    }

    /** \copydoc Dune::InStreamInterface::readFloat */
    void readFloat ( float &value )
    {
      readPrimitive( value );
    }

    /** \copydoc Dune::InStreamInterface::readInt */
    void readInt ( int &value )
    {
      readPrimitive( value );
    }

    /** \copydoc Dune::InStreamInterface::readChar */
    void readChar ( char &value )
    {
      readPrimitive( value );
    }

    /** \copydoc Dune::InStreamInterface::readBool */
    void readBool ( bool &value )
    {
      readPrimitive( value );
    }

    /** \copydoc Dune::InStreamInterface::readString */
    void readString ( std::string &s )
    {
      unsigned int length;
      readPrimitive( length );
      for( unsigned int i = 0; i < length; ++i )
      {
        char c;
        readPrimitive( c );
        s += c;
      }
    }

    /** \copydoc Dune::InStreamInterface::readUnsignedInt */
    void readUnsignedInt ( unsigned int &value )
    {
      readPrimitive( value );
    }

    /** \copydoc Dune::InStreamInterface::readUnsignedLong */
    void readUnsignedLong ( unsigned long &value )
    {
      readPrimitive( value );
    }

  protected:
    bool valid () const
    {
      return stream_.good() | stream_.eof();
    }

    template< class T >
    void readPrimitive ( T &value )
    {
      union { T value; char bytes[ sizeof( T ) ]; } convert;
      stream_.read( convert.bytes, sizeof( T ) );
      value = convert.value;
      if( !valid() )
        readError();
    }

  protected:
    std::ifstream stream_;
  };

} // end namespace Dune

#endif // #ifndef DUNE_FEM_BINARYSTREAMS_HH
