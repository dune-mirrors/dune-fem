#ifndef DUNE_FEM_STANDARDSTREAMS_HH
#define DUNE_FEM_STANDARDSTREAMS_HH

#include <fstream>
#include <utility>
#include <memory>

#ifdef SYSTEM_ENDIAN_HEADER 
#include SYSTEM_ENDIAN_HEADER
#endif

#include <dune/fem/io/streams/streams.hh>

namespace Dune
{

  namespace Fem 
  {
    struct ByteOrder 
    {
      static const char defaultEndian = 0;
      static const char order = 
#if __BYTE_ORDER == __LITTLE_ENDIAN
          0 ;
#elif __BYTE_ORDER == __BIG_ENDIAN
          1 ;
#else 
          0 ; // default is zero (in case no endian header was found)
#endif 

      static inline size_t map( const char storedOrder, 
                                const size_t pos, 
                                const size_t size )
      {
        // if byte order differs, swap 
        return ( order == storedOrder ) ? pos : ( size - pos - 1 );
      }

      static inline size_t map( const size_t pos, 
                                const size_t size )
      {
        // if byte order is not big endian, swap 
        return map( 1, pos, size );
      }
    };

    class StandardOutStream; 
    class StandardInStream;
    
    struct StandardOutStreamTraits
    {
      typedef StandardOutStream OutStreamType;
    };
    
    /** \class StandardOutStream
     *  \ingroup InOutStreams
     *  \brief output stream writing into a given std::ostream 
     *
     *  \note This stream directly stores the binary representation of the data.
     *        The binary representation of the stored data is always that of the
     *        current machine. On read the data is converted accordingly on machines 
     *        with different endianess.
     *
     *  \newimplementation
     */
    class StandardOutStream
    : public OutStreamInterface< StandardOutStreamTraits >
    {
      typedef StandardOutStream ThisType;
      typedef OutStreamInterface< StandardOutStreamTraits > BaseType;

    public:
      //! type of the traits
      typedef StandardOutStreamTraits Traits;
      
    protected:
      using BaseType::writeError;

    public:
      /** \brief constructor
       *
       *  \param[in]  stream  std::ostream to write to 
       */
      explicit StandardOutStream ( std::ostream& stream )
      : stream_( stream ) 
      {
        if( ! stream ) 
          DUNE_THROW( Dune::IOError, "Stream not valid!" );
      }

      /** return reference to internal ostream */
      std::ostream& stream() { return stream_; }

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
        return bool( stream_ );
      }

      template< class T >
      void writePrimitive ( const T &value )
      {
        const size_t tsize = sizeof( T );
        union { T value; char bytes[ tsize ]; } convert;

        // copy  value 
        convert.value = value;

        assert( sizeof(char) == 1 ) ;
        // write according to byte order 
        for( size_t i=0; i<tsize; ++i ) 
        {
          stream_.write( &convert.bytes[ ByteOrder :: map( i, tsize ) ], 1 );
        }

        if( !valid () )
          writeError();
      }

    protected:
      std::ostream& stream_;
    };

    
    struct StandardInStreamTraits
    {
      typedef StandardInStream InStreamType;
    };
    

    /** \class StandardInStream
     *  \ingroup InOutStreams
     *  \brief input stream reading from a given std::istream 
     *
     *  \note This stream directly stores the binary representation of the data.
     *        The binary representation of the stored data is always that of the
     *        current machine. On read the data is converted accordingly on machines 
     *        with different endianess.
     *
     *  \newimplementation
     */
    class StandardInStream
    : public InStreamInterface< StandardInStreamTraits >
    {
      typedef StandardInStream ThisType;
      typedef InStreamInterface< StandardInStreamTraits > BaseType;

    public:
      //! type of the traits
      typedef StandardInStreamTraits Traits;
      
    protected:
      using BaseType::readError;

    public:
      /** \brief constructor
       *
       *  \param[in]  stream  std::istream to read from 
       */
      explicit StandardInStream ( std::istream& stream )
      : stream_( stream )
      {
        if( ! valid() ) 
          DUNE_THROW( Dune::IOError, "Stream not valid!" );
      }

      /** return reference to internal istream */
      std::istream& stream() { return stream_; }

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
        // resize string 
        s.resize( length );
        for( unsigned int i = 0; i < length; ++i )
        {
          readPrimitive( s[ i ] );
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

        assert( sizeof(char) == 1 ) ;
        // read according to byte order 
        for( size_t i=0; i<tsize; ++i ) 
        {
          stream_.read( &convert.bytes[ ByteOrder :: map( i, tsize ) ], 1 );
        }

        // store result to value 
        value = convert.value;

        if( !valid() )
          readError();
      }

    protected:
      std::istream& stream_;
      char storedOrder_ ;
    };

  } // end namespace Fem   

  // #if DUNE_FEM_COMPATIBILITY  
  // put this in next version 1.4 

  using Fem :: StandardOutStream ;
  using Fem :: StandardInStream ;
  // #endif // DUNE_FEM_COMPATIBILITY

} // end namespace Dune

#endif // #ifndef DUNE_FEM_BINARYSTREAMS_HH
