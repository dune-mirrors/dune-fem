#ifndef DUNE_FEM_DATASTREAMS_HH
#define DUNE_FEM_DATASTREAMS_HH

#include <cassert>
#include <string>
#include <fstream>

#include <endian.h>

#include <dune/common/exceptions.hh>
#include <dune/fem/io/streams/streams.hh>
#include <dune/fem/io/streams/binarystreams.hh>

#if HAVE_ALUGRID
// inlcude alugrid to have to communicator class from ALUGrid 
#include <dune/grid/alugrid/3d/alugrid.hh>
#endif

namespace Dune
{
  namespace Fem 
  { 

    template< class OutStreamImp >
    struct DataOutStreamTraits 
    {
      typedef OutStreamImp OutStreamType;
    };

   
    /** \class DataOutStream
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
    class DataBasicOutStream
    : public OutStreamInterface< DataOutStreamTraits< OutStreamImp > >
    {
      typedef DataBasicOutStream< OutStreamImp > ThisType;
      typedef OutStreamInterface< DataOutStreamTraits< OutStreamImp > > BaseType;

    public:
      //! type of the implementaton (Barton-Nackman)
      typedef OutStreamImp OutStreamType;
      
      //! type of the traits
      typedef DataOutStreamTraits< OutStreamImp > Traits;

    protected:
      using BaseType::writeError;
      using BaseType::asImp; 

    protected:
      DataBasicOutStream(  )
      {}

    public:
      /** \copydoc Dune::OutStreamInterface::writeDouble */
      void writeDouble ( double value )
      {
        asImp().write( value );
      }

      /** \copydoc Dune::OutStreamInterface::writeFloat */
      void writeFloat ( float value )
      {
        asImp().write( value );
      }

      /** \copydoc Dune::OutStreamInterface::writeInt */
      void writeInt ( int value )
      {
        asImp().write( value );
      }

      /** \copydoc Dune::OutStreamInterface::writeChar */
      void writeChar ( char value )
      {
        asImp().write( value );
      }

      /** \copydoc Dune::OutStreamInterface::writeBool */
      void writeBool ( const bool value )
      {
        asImp().write( value );
      }

      /** \copydoc Dune::OutStreamInterface::writeUnsignedInt */
      void writeUnsignedInt ( unsigned int value )
      {
        asImp().write( value );
      }

      /** \copydoc Dune::OutStreamInterface::writeUnsignedInt64 */
      void writeUnsignedInt64 ( uint64_t value )
      {
        asImp().write( value );
      }
    };


    template< class InStreamImp >
    struct DataInStreamTraits 
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
    class DataBasicInStream
    : public InStreamInterface< DataInStreamTraits< InStreamImp > >
    {
      typedef DataBasicInStream< InStreamImp > ThisType;
      typedef InStreamInterface< DataInStreamTraits< InStreamImp > > BaseType;

    public:
      //! type of the implementation (Barton-Nackman)
      typedef InStreamImp InStreamType;
      
      //! type of the traits
      typedef DataInStreamTraits< InStreamType > Traits;

    protected:
      using BaseType :: asImp;

      DataBasicInStream ()
      {}

    public:
      /** \copydoc Dune::InStreamInterface::readDouble */
      void readDouble ( double &value )
      {
        asImp().read( value );
      }

      /** \copydoc Dune::InStreamInterface::readFloat */
      void readFloat ( float &value )
      {
        asImp().read( value );
      }

      /** \copydoc Dune::InStreamInterface::readInt */
      void readInt ( int &value )
      {
        asImp().read( value );
      }

      /** \copydoc Dune::InStreamInterface::readChar */
      void readChar ( char &value )
      {
        asImp().read( value );
      }

      /** \copydoc Dune::InStreamInterface::readBool */
      void readBool ( bool &value )
      {
        asImp().read( value );
      }

      /** \copydoc Dune::InStreamInterface::readUnsignedInt */
      void readUnsignedInt ( unsigned int &value )
      {
        asImp().read( value );
      }

      /** \copydoc Dune::InStreamInterface::readUnsignedInt64 */
      void readUnsignedInt64 ( uint64_t &value )
      {
        asImp().read( value );
      }
    };


    struct DataObjectStream 
#if HAVE_ALUGRID
      : public ALU3DSPACE ObjectStream 
#endif
    {
#if HAVE_ALUGRID
      typedef ALU3DSPACE ObjectStream BaseType ;
      using BaseType :: _buf ;
      using BaseType :: operator =;
#else 
      DataObjectStream () 
      {
        DUNE_THROW(NotImplemented,"DataOutStream only working with ALUGrid enabled!");
      }

      char* _buf ;
      size_t size() const { return 0; }

      DataObjectStream& operator = (std::pair< char* , int > & osvec)
      {
        return *this;
      }
#endif
      // return pointer to buffer 
      char* buffer () { return _buf; }
      // return pointer to buffer (const)
      const char* buffer () const { return _buf; }
    };


    /** \class XDRFileOutStream
     *  \ingroup InOutStreams
     *  \brief XDR output stream writing into a file
     *
     *  \newimplementation
     */
    class DataOutStream
    : public DataBasicOutStream< DataOutStream >
    {
    private:
      typedef DataOutStream   ThisType;
      typedef DataBasicOutStream< ThisType > BaseType;

    protected:
      DataObjectStream outstream_;
      const std::string filename_;

    public:
      using BaseType :: writeString ;

      /** \brief constructor
       *
       *  \param[in]  filename  name of the file to write to
       */
      explicit DataOutStream ( const std::string &filename )
        : filename_( filename )
      {
#if HAVE_ALUGRID
        outstream_.clear();
#else 
        DUNE_THROW(NotImplemented,"DataOutStream only working with ALUGrid enabled!");
#endif
      }

      /** \copydoc Dune::OutStreamInterface::writeString */
      void writeString( const std::string& s ) 
      {
        const size_t size = s.size();
        // write size 
        write( size );
        // write each character 
        for( size_t i=0; i<size; ++i ) 
          write( s[i] );
      }

      //! write data to stream 
      template <class T> 
      void write( const T& value ) 
      {
#if HAVE_ALUGRID
        const size_t tsize = sizeof( T ) ;
        union { T value; char bytes[ tsize ]; } convert;
        convert.value = value;
        outstream_.reserve( outstream_.size() + tsize );
        for( size_t i=0; i<tsize; ++i ) 
        {
          outstream_.write( convert.bytes[ ByteOrder :: map( i, tsize ) ] );
        }
#else
        DUNE_THROW(NotImplemented,"DataOutStream only working with ALUGrid enabled!");
#endif
      }

      /** \brief return size of internal buffer */
      size_t size() const { return outstream_.size(); }
      
      /** \brief return pointer to buffer  */
      char* buffer() { return outstream_.buffer() ; }
      const char* buffer() const { return outstream_.buffer() ; }

      /** \copydoc Dune::OutStreamInterface::flush */
      inline void flush ()
      {
        std::ofstream file ( filename_.c_str(), std::ios::binary );
        file << size();
        char* buffer = outstream_.buffer();
        file.write( buffer, size() );
        file.flush();
      }
    };


    /** \class DataInStream
     *  \ingroup InOutStreams
     *  \brief stream reading from a char* buffer 
     *
     *  \newimplementation
     */
    class DataInStream
    : public DataBasicInStream< DataInStream >
    {
    private:
      typedef DataInStream   ThisType;
      typedef DataBasicInStream< ThisType > BaseType;

    protected:
      DataObjectStream instream_;

    public:
      /** \brief constructor
       *
       *  \param[in]  filename  name of the file to write to
       *  \param[in]  append    true if stream appends to file filename
       *                        (default = false)
       */
      explicit DataInStream ( const std::string &filename,
                              const size_t pos = 0 )
      {
        std::ifstream file( filename.c_str(), std::ios::binary );
        size_t size = 0 ;
        file >> size ;
        char* buffer = DataObjectStream :: allocateBuffer( size );
        file.read( buffer, size );
        std::pair< char* , int > buff( buffer, int(size) );
        instream_ = buff ;
      }

      /** \copydoc Dune::InStreamInterface::readString */
      void readString( std::string& s ) 
      {
        size_t size = 0;
        read( size );
        // resize the string 
        s.resize( size );
        // read all characters 
        for( size_t i=0; i<size; ++i ) 
          read( s[i] );
      }

      //! write data to stream 
      template <class T> 
      void read( T& value ) 
      {
#if HAVE_ALUGRID
        const size_t tsize = sizeof( T ) ;
        union { T value; char bytes[ tsize ]; } convert;
        for( size_t i=0; i<tsize; ++i ) 
          instream_.read( convert.bytes[ ByteOrder :: map( i, tsize ) ] );
        value = convert.value;
#else
        DUNE_THROW(NotImplemented,"DataOutStream only working with ALUGrid enabled!");
#endif
      }

#if HAVE_ALUGRID
      /** \brief return size of internal buffer */
      size_t size() const { return instream_.size(); }
      
      /** \brief return pointer to buffer  */
      const char* buffer() const { return instream_.buffer() ; }
#endif
    };

  } // end  namespace Fem   

} // end namespace Dune

#endif // #ifndef DUNE_FEM_XDRSTREAMS_HH
