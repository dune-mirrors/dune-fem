#ifndef DUNE_FEM_ASCIISTREAMS_HH
#define DUNE_FEM_ASCIISTREAMS_HH

#include <iostream>

#include <dune/fem/io/streams/streams.hh>

namespace Dune
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
  public:
    //! type of the traits
    typedef ASCIIOutStreamTraits Traits;
    
  private:
    typedef ASCIIOutStream ThisType;
    typedef OutStreamInterface< Traits > BaseType;

  protected:
    std :: ostream &stream_;
    bool mustFreeStream_;

  protected:
    using BaseType :: writeError;

  public:
    /** \brief constructor
     *
     *  \param[in]  stream  STL output stream to write to
     */
    inline explicit ASCIIOutStream ( std :: ostream &stream )
    : stream_( stream ),
      mustFreeStream_( false )
    {
    }

    /** \brief constructor
     *
     *  \param[in]  filename  name of a file to write to
     */
    inline explicit ASCIIOutStream ( const std :: string filename )
    : stream_( *(new std :: ofstream( filename.c_str() )) ),
      mustFreeStream_( true )
    {
    }

    /** \brief destructor */
    inline ~ASCIIOutStream ()
    {
      if( mustFreeStream_ )
        delete &stream_;
    }

    /** \copydoc Dune::OutStreamInterface::flush */
    inline void flush ()
    {
      stream_.flush();
    }

    /** \copydoc Dune::OutStreamInterface::writeDouble */
    inline void writeDouble ( double value )
    {
      stream_.setf( std ::ios_base :: scientific, std :: ios_base :: floatfield );
      stream_ .precision( 16 );
      stream_ << value << std :: endl;
      if( !valid () )
        writeError();
    }

    /** \copydoc Dune::OutStreamInterface::writeFloat */
    inline void writeFloat ( float value )
    {
      stream_.setf( std ::ios_base :: scientific, std :: ios_base :: floatfield );
      stream_ .precision( 7 );
      stream_ << value << std :: endl;
      if( !valid() )
        writeError();
    }

    /** \copydoc Dune::OutStreamInterface::writeInt */
    inline void writeInt ( int value )
    {
      stream_ << value << std :: endl;
      if( !valid () )
        writeError();
    }

    /** \copydoc Dune::OutStreamInterface::writeString
     *
     *  \note Strings containing newline characters will not be read back
     *        correctly.
     */
    inline void writeString ( const std :: string &s )
    {
      stream_ << s << std :: endl;
      if( !valid () )
        writeError();
    }

    /** \copydoc Dune::OutStreamInterface::writeUnsignedInt */
    inline void writeUnsignedInt ( unsigned int value )
    {
      stream_ << value << std :: endl;
      if( !valid () )
        writeError();
    }

  protected:
    inline bool valid () const
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
  public:
    //! type of the traits
    typedef ASCIIInStreamTraits Traits;
    
  private:
    typedef ASCIIInStream ThisType;
    typedef InStreamInterface< Traits > BaseType;

  protected:
    std :: istream &stream_;
    bool mustFreeStream_;

  protected:
    using BaseType :: readError;

  public:
    /** \brief constructor
     *
     *  \param[in]  stream  STL output stream to write to
     */
    inline ASCIIInStream ( std :: istream &stream )
    : stream_( stream ),
      mustFreeStream_( false )
    {
    }

    /** \brief constructor
     *
     *  \param[in]  filename  name of a file to write to
     */
    inline ASCIIInStream ( const std :: string filename )
    : stream_( *(new std :: ifstream( filename.c_str() )) ),
      mustFreeStream_( true )
    {
    }

    /** \brief destructor */
    inline ~ASCIIInStream ()
    {
      if( mustFreeStream_ )
        delete &stream_;
    }

    /** \copydoc Dune::InStreamInterface::readDouble */
    inline void readDouble ( double &value )
    {
      stream_ >> value;
      if( !valid () )
        readError();
    }

    /** \copydoc Dune::InStreamInterface::readFloat */
    inline void readFloat ( float &value )
    {
      stream_ >> value;
      if( !valid () )
        readError();
    }

    /** \copydoc Dune::InStreamInterface::readInt */
    inline void readInt ( int &value )
    {
      stream_ >> value;
      if( !valid () )
        readError();
    }

    /** \copydoc Dune::InStreamInterface::readString
     *
     *  \note Strings containing newline characters will not be read back
     *        correctly.
     */
    inline void readString ( std :: string &s )
    {
      stream_ >> s;
      if( !valid () )
        readError();
    }

    /** \copydoc Dune::InStreamInterface::readUnsignedInt */
    inline void readUnsignedInt ( unsigned int value )
    {
      stream_ >> value;
      if( !valid () )
        readError();
    }

  protected:
    inline bool valid () const
    {
      return stream_.good() | stream_.eof();
    }
  };

}

#endif

