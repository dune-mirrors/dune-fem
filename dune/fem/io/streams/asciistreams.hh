#ifndef DUNE_FEM_ASCIISTREAMS_HH
#define DUNE_FEM_ASCIISTREAMS_HH

#include <iostream>
#include <fstream>

#include <dune/fem/io/streams/streams.hh>

namespace Dune
{

  namespace Fem
  {

    class ASCIIOutStream;
    class ASCIIInStream;

    struct ASCIIOutStreamTraits
    {
      typedef ASCIIOutStream OutStreamType;
    };


    /** \class ASCIIOutStream
     *  \ingroup InOutStreams
     *  \brief output stream writing into an STL output stream using ASCII
     *         encoding
     *
     *  This writes the data into an STL output stream. The data is written in
     *  ASCII format, each basic type on a single line.
     *
     *  \newimplementation
     */
    class ASCIIOutStream
    : public OutStreamInterface< ASCIIOutStreamTraits >
    {
      typedef ASCIIOutStream ThisType;
      typedef OutStreamInterface< ASCIIOutStreamTraits > BaseType;

    public:
      //! type of the traits
      typedef ASCIIOutStreamTraits Traits;

    protected:
      std::ostream &stream_;
      bool mustFreeStream_;

    protected:
      using BaseType::writeError;

    public:
      /** \brief constructor
       *
       *  \param[in]  stream  STL output stream to write to
       */
      explicit ASCIIOutStream ( std::ostream &stream )
      : stream_( stream ),
        mustFreeStream_( false )
      {}

      /** \brief constructor
       *
       *  \param[in]  filename  name of a file to write to
       */
      explicit ASCIIOutStream ( const std::string &filename )
      : stream_( *(new std :: ofstream( filename.c_str() )) ),
        mustFreeStream_( true )
      {}

      /** \brief destructor */
      ~ASCIIOutStream ()
      {
        if( mustFreeStream_ )
          delete &stream_;
      }

      /** \copydoc Dune::Fem::OutStreamInterface::flush */
      void flush ()
      {
        stream_.flush();
      }

      /** \copydoc Dune::Fem::OutStreamInterface::writeDouble */
      void writeDouble ( const double value )
      {
        stream_.setf( std ::ios_base :: scientific, std :: ios_base :: floatfield );
        stream_ .precision( 16 );
        stream_ << value << std :: endl;
        if( !valid () )
          writeError();
      }

      /** \copydoc Dune::Fem::OutStreamInterface::writeFloat */
      void writeFloat ( const float value )
      {
        stream_.setf( std ::ios_base :: scientific, std :: ios_base :: floatfield );
        stream_ .precision( 7 );
        stream_ << value << std :: endl;
        if( !valid() )
          writeError();
      }

      /** \copydoc Dune::Fem::OutStreamInterface::writeInt */
      void writeInt ( const int value )
      {
        stream_ << value << std :: endl;
        if( !valid () )
          writeError();
      }

      /** \copydoc Dune::Fem::OutStreamInterface::writeSignedInt64 */
      void writeSignedInt64 ( int64_t value )
      {
        stream_ << value << std::endl;
        if( !valid () )
          writeError();
      }

      /** \copydoc Dune::Fem::OutStreamInterface::writeChar */
      void writeChar ( const char value )
      {
        // make sure char is written as number
        int val = (int) value;
        writeInt( val );
      }

      /** \copydoc Dune::Fem::OutStreamInterface::writeChar */
      void writeBool ( const bool value )
      {
        std::string val( ( value == true ) ? "true" : "false" );
        writeString( val );
      }

      /** \copydoc Dune::Fem::OutStreamInterface::writeString
       *
       *  \note Strings containing newline characters will not be read back
       *        correctly.
       */
      void writeString ( const std::string &s )
      {
        const unsigned int length = s.length();
        stream_ << length;
        for( unsigned int i = 0; i < length; ++i )
          stream_.put( s[ i ] );
        stream_ << std :: endl;
        if( !valid () )
          writeError();
      }

      /** \copydoc Dune::Fem::OutStreamInterface::writeUnsignedInt */
      void writeUnsignedInt ( unsigned int value )
      {
        stream_ << value << std::endl;
        if( !valid () )
          writeError();
      }

      /** \copydoc Dune::Fem::OutStreamInterface::writeUnsignedInt64 */
      void writeUnsignedInt64 ( uint64_t value )
      {
        stream_ << value << std::endl;
        if( !valid () )
          writeError();
      }

    protected:
      bool valid () const
      {
        return stream_.good() | stream_.eof();
      }
    };



    struct ASCIIInStreamTraits
    {
      typedef ASCIIInStream InStreamType;
    };



    /** \class ASCIIInStream
     *  \ingroup InOutStreams
     *  \brief input stream reading from an STL input stream using ASCII decoding
     *
     *  This writes the data into an STL output stream. The data is written in
     *  ASCII format, each basic type on a single line.
     *
     *  \newimplementation
     */
    class ASCIIInStream
    : public InStreamInterface< ASCIIInStreamTraits >
    {
      typedef ASCIIInStream ThisType;
      typedef InStreamInterface< ASCIIInStreamTraits > BaseType;

    public:
      //! type of the traits
      typedef ASCIIInStreamTraits Traits;

    protected:
      std::istream &stream_;
      bool mustFreeStream_;

    protected:
      using BaseType::readError;

    public:
      /** \brief constructor
       *
       *  \param[in]  stream  STL output stream to write to
       */
      explicit ASCIIInStream ( std::istream &stream )
      : stream_( stream ),
        mustFreeStream_( false )
      {}

      /** \brief constructor
       *
       *  \param[in]  filename  name of a file to write to
       */
      ASCIIInStream ( const std::string &filename )
      : stream_( *(new std :: ifstream( filename.c_str() )) ),
        mustFreeStream_( true )
      {}

      /** \brief destructor */
      ~ASCIIInStream ()
      {
        if( mustFreeStream_ )
          delete &stream_;
      }

      /** \copydoc Dune::Fem::InStreamInterface::readDouble */
      void readDouble ( double &value )
      {
        stream_ >> value;
        if( !valid () )
          readError();
      }

      /** \copydoc Dune::Fem::InStreamInterface::readFloat */
      void readFloat ( float &value )
      {
        stream_ >> value;
        if( !valid () )
          readError();
      }

      /** \copydoc Dune::Fem::InStreamInterface::readInt */
      void readInt ( int &value )
      {
        stream_ >> value;
        if( !valid () )
          readError();
      }

      /** \copydoc Dune::Fem::InStreamInterface::readSignedInt64 */
      void readSignedInt64 (int64_t &value )
      {
        stream_ >> value;
        if( !valid () )
          readError();
      }

      /** \copydoc Dune::Fem::InStreamInterface::readChar */
      void readChar ( char &value )
      {
        int val;
        readInt( val );
        value = (char) val;
      }

      /** \copydoc Dune::Fem::InStreamInterface::readBool */
      void readBool ( bool &value )
      {
        std::string val;
        readString( val );

        if( val == "true" )
          value = true;
        else if ( val == "false" )
          value = false;
        else
          readError();
      }

      /** \copydoc Dune::Fem::InStreamInterface::readString
       *
       *  \note Strings containing newline characters will not be read back
       *        correctly.
       */
      void readString ( std::string &s )
      {
        unsigned int length;
        stream_ >> length;
        for( unsigned int i = 0; i < length; ++i )
          s += stream_.get();
        if( !valid () )
          readError();
      }

      /** \copydoc Dune::Fem::InStreamInterface::readUnsignedInt */
      void readUnsignedInt ( unsigned int &value )
      {
        stream_ >> value;
        if( !valid () )
          readError();
      }

      /** \copydoc Dune::Fem::InStreamInterface::readUnsignedInt64 */
      void readUnsignedInt64 (uint64_t &value )
      {
        stream_ >> value;
        if( !valid () )
          readError();
      }

    protected:
      bool valid () const
      {
        return stream_.good() | stream_.eof();
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ASCIISTREAMS_HH
