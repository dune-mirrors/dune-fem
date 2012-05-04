#ifndef DUNE_FEM_BINARYSTREAMS_HH
#define DUNE_FEM_BINARYSTREAMS_HH

#include <fstream>
#include <endian.h>

#include <dune/fem/io/streams/streams.hh>

namespace Dune
{

  namespace Fem 
  {
    struct ByteOrder 
    {
      static const bool swap = 
#if __BYTE_ORDER == __LITTLE_ENDIAN
          true ;
#elif __BYTE_ORDER == __BIG_ENDIAN
          false ;
#endif 

      static inline size_t map( const size_t pos, const size_t size )
      {
        return swap ? ( size - pos - 1 ) : pos ;
      }
    };

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
     *        The binary representation of the stored data is always big endian. 
     *        On read/write the data is converted accordingly on little endian machines.
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
        const size_t tsize = sizeof( T ) ;
        union { T value; char bytes[ tsize ]; } convert;

        // copy  value 
        convert.value = value;

        // if byte order needs a swap 
        if( ByteOrder :: swap ) 
        { 
          union { T value; char bytes[ tsize ]; } convswap;
          convswap.value = convert.value; 
          for( size_t i=0; i<tsize; ++i ) 
          {
            convert.bytes[ i ] = convswap.bytes[ ByteOrder :: map( i, tsize ) ];
          }
        } 

        // write value 
        stream_.write( convert.bytes, tsize );
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
        const size_t tsize = sizeof( T ) ;
        union { T value; char bytes[ tsize ]; } convert;

        // read value 
        stream_.read( convert.bytes, tsize );

        // if byte order needs a swap 
        if( ByteOrder :: swap ) 
        { 
          union { T value; char bytes[ tsize ]; } convswap;
          convswap.value = convert.value; 
          for( size_t i=0; i<tsize; ++i ) 
          {
            convert.bytes[ i ] = convswap.bytes[ ByteOrder :: map( i, tsize ) ];
          }
        } 

        // store result to value 
        value = convert.value;

        if( !valid() )
          readError();
      }

    protected:
      std::ifstream stream_;
    };

  } // end namespace Fem   


  // #if DUNE_FEM_COMPATIBILITY  
  // put this in next version 1.4 

  using Fem :: BinaryFileOutStream ;
  using Fem :: BinaryFileInStream ;
  // #endif // DUNE_FEM_COMPATIBILITY

} // end namespace Dune

#endif // #ifndef DUNE_FEM_BINARYSTREAMS_HH
