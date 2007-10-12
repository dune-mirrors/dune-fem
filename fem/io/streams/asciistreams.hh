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

  

  class ASCIIOutStream
  : public OutStreamInterface< ASCIIOutStreamTraits >
  {
  public:
    typedef ASCIIOutStreamTraits Traits;
    
  private:
    typedef ASCIIOutStream ThisType;
    typedef OutStreamInterface< Traits > BaseType;

  protected:
    std :: ostream &stream_;
    bool mustFreeStream_;

  public:
    inline explicit ASCIIOutStream ( std :: ostream &stream )
    : stream_( stream ),
      mustFreeStream_( false )
    {
    }

    inline explicit ASCIIOutStream ( const std :: string filename )
    : stream_( *(new std :: ofstream( filename.c_str() )) ),
      mustFreeStream_( true )
    {
    }

    inline ~ASCIIOutStream ()
    {
      if( mustFreeStream_ )
        delete &stream_;
    }

    inline void flush ()
    {
      stream_.flush();
    }

    inline bool valid () const
    {
      return stream_.good() | stream_.eof();
    }

    inline void writeDouble ( double value )
    {
      stream_.setf( std ::ios_base :: scientific, std :: ios_base :: floatfield );
      stream_ .precision( 16 );
      stream_ << value << std :: endl;
    }

    inline void writeFloat ( float value )
    {
      stream_.setf( std ::ios_base :: scientific, std :: ios_base :: floatfield );
      stream_ .precision( 7 );
      stream_ << value << std :: endl;
    }

    inline void writeInt ( int value )
    {
      stream_ << value << std :: endl;
    }

    inline void writeString ( const std :: string &s )
    {
      stream_ << s << std :: endl;
    }

    inline void writeUnsignedInt ( unsigned int value )
    {
      stream_ << value << std :: endl;
    }
  };


  
  struct ASCIIInStreamTraits
  {
    typedef ASCIIInStream InStreamType;
  };

  

  class ASCIIInStream
  : public InStreamInterface< ASCIIInStreamTraits >
  {
  public:
    typedef ASCIIInStreamTraits Traits;
    
  private:
    typedef ASCIIInStream ThisType;
    typedef InStreamInterface< Traits > BaseType;

  protected:
    std :: istream &stream_;
    bool mustFreeStream_;

  public:
    inline ASCIIInStream ( std :: istream &stream )
    : stream_( stream ),
      mustFreeStream_( false )
    {
    }

    inline ASCIIInStream ( const std :: string filename )
    : stream_( *(new std :: ifstream( filename.c_str() )) ),
      mustFreeStream_( true )
    {
    }

    inline ~ASCIIInStream ()
    {
      if( mustFreeStream_ )
        delete &stream_;
    }

    inline bool valid () const
    {
      return stream_.good() | stream_.eof();
    }

    inline void readDouble ( double &value )
    {
      stream_ >> value;
    }

    inline void readFloat ( float &value )
    {
      stream_ >> value;
    }

    inline void readInt ( int &value )
    {
      stream_ >> value;
    }

    inline void readString ( std :: string &s )
    {
      stream_ >> s;
    }

    inline void readUnsignedInt ( unsigned int value )
    {
      stream_ >> value;
    }
  };

}

#endif

