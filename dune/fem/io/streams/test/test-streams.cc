#include <config.h>

#include <dune/fem/version.hh>
#include <dune/fem/misc/mpimanager.hh>

#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/streams/asciistreams.hh>
#include <dune/fem/io/streams/binarystreams.hh>
#include <dune/fem/io/streams/sionlibstreams.hh>

using namespace Dune;
using namespace Fem;



struct Data
{
  std :: string my_string;
  unsigned int my_uint;
  uint64_t my_uint64;
  int64_t my_int64;
  unsigned long int my_ulong;
  long int my_long;
  size_t my_sizet;
  double my_double;
  int my_int;
  float my_float;
  char my_char;
  bool my_bool;
};

template< class Traits >
void write ( OutStreamInterface< Traits > &out, const Data &data )
{
  unsigned int versionId = DUNE_MODULE_VERSION_ID(DUNE_FEM);
  out << versionId;

  out << data.my_string
      << data.my_uint
      << data.my_uint64
      << data.my_int64
      << data.my_ulong
      << data.my_long
      << data.my_sizet
      << data.my_double
      << data.my_int
      << data.my_float
      << data.my_char
      << data.my_bool;

  out << DUNE_MODULE_VERSION_ID(DUNE_FEM);
}

template< class Traits >
bool read ( InStreamInterface< Traits > &in, const Data &data )
{
  unsigned int versionId = in.readUnsignedInt();
  if( versionId != DUNE_MODULE_VERSION_ID(DUNE_FEM) )
  {
    std :: cerr << "Incorrect versionId read back: "
                << versionId << " != " << DUNE_MODULE_VERSION_ID(DUNE_FEM)
                << std :: endl;
    return false;
  }

  Data check;
  in >> check.my_string
     >> check.my_uint
     >> check.my_uint64
     >> check.my_int64
     >> check.my_ulong
     >> check.my_long
     >> check.my_sizet
     >> check.my_double
     >> check.my_int >> check.my_float >> check.my_char >> check.my_bool;

  std :: cerr << "Data: " << check.my_string
              << ", " << check.my_uint
              << ", " << check.my_uint64
              << ", " << check.my_int64
              << ", " << check.my_ulong
              << ", " << check.my_long
              << ", " << check.my_sizet
              << ", " << check.my_double
              << ", " << check.my_int
              << ", " << check.my_float
              << ", " << (int)check.my_char
              << ", " << check.my_bool
              << std :: endl;

  versionId = in.readUnsignedInt();
  if( versionId != DUNE_MODULE_VERSION_ID(DUNE_FEM) )
  {
    std :: cerr << "Incorrect versionId read back: "
                << versionId << " != " << DUNE_MODULE_VERSION_ID(DUNE_FEM)
                << std :: endl << std::endl;
    return false;
  }

  std::cerr << std :: endl ;

  bool equal = true;
  equal &= (data.my_string == check.my_string);
  equal &= (data.my_uint   == check.my_uint);
  equal &= (data.my_uint64 == check.my_uint64);
  equal &= (data.my_int64  == check.my_int64);
  equal &= (data.my_ulong  == check.my_ulong);
  equal &= (data.my_long   == check.my_long);
  equal &= (data.my_sizet  == check.my_sizet);
  equal &= (data.my_double == check.my_double);
  equal &= (data.my_int    == check.my_int);
  equal &= (data.my_float  == check.my_float);
  equal &= (data.my_char   == check.my_char );
  equal &= (data.my_bool   == check.my_bool);
  return equal;
}

int main ( int argc, char** argv )
{
  MPIManager::initialize( argc, argv );

  try
  {
    const bool writeStreams = ( argc > 1 ) ? false : true ;

    Data data;
    data.my_string = "Hello, World!";
    data.my_uint   = 42;
    data.my_uint64 = -4 ; // this results in 18446744073709551612
    data.my_int64  = -4 ; // this results in -4
    data.my_ulong  = uint32_t(-4) ; // this results in 4294967292
    data.my_long   = -long(uint32_t(-4)); // this results in -4294967292
    data.my_sizet  = uint32_t(-16) ; // this results in 4294967280
    data.my_double = 1.2345678901234;
    data.my_int    = -767;
    data.my_float  = 1.23456;
    data.my_char   = 123;
    data.my_bool   = true;

    bool failed = false ;
    std::stringstream filestr;
    std::cout << "Path: "<< Parameter :: commonOutputPath() << std::endl;
    filestr << Parameter :: commonOutputPath() << "/test." << MPIManager :: rank() << ".";
    {
      std::string filename( filestr.str() + "ascii" );
      std :: cerr << "Checking ASCII streams..." << std :: endl;
      if( writeStreams )
      {
        Fem :: ASCIIOutStream aout( filename.c_str() );
        write( aout, data );
        aout.flush();
      }
      Fem :: ASCIIInStream ain( filename.c_str() );
      if( !read( ain, data ) )
        failed = true ;
    }

    {
      std::string filename( filestr.str() + "binary" );
      std :: cerr << "Checking Binary streams..." << std :: endl;
      if( writeStreams )
      {
        Fem :: BinaryFileOutStream bout( filename.c_str() );
        write( bout, data );
        bout.flush();
      }
      Fem :: BinaryFileInStream bin( filename.c_str() );
      if( !read( bin, data ) )
        failed = true ;
    }

#if HAVE_SIONLIB
    {
      std :: cerr << "Checking SIONlib streams..." << std :: endl;
      std::stringstream file;
      file << Parameter::commonOutputPath() << "/test.sion." << MPIManager :: size() ;
      std::string filename( file.str() );
      if( writeStreams )
      {
        Fem :: SIONlibOutStream sionout( filename.c_str() );
        write( sionout, data );
        sionout.flush();
      }
      // check parallel read
      {
        Fem :: SIONlibInStream sionin( filename.c_str() );
        if( ! read( sionin, data ) )
          failed = true ;
      }

      if( MPIManager :: size () == 1 )
      {
        int size = (argc > 2) ? atoi( argv[ 2 ] ) : MPIManager :: size();
        std::cout << "Check serial read for " << size << " procs" << std::endl;
        std::stringstream file;
        file <<  Parameter::commonOutputPath() << "/test.sion." << size ;
        // check serial read
        for( int rank=0; rank<size; ++rank )
        {
          Fem :: SIONlibInStream sionin( filename.c_str(), rank );
          if( ! read( sionin, data ) )
            failed = true ;
        }
      }

    }
#endif

    if( failed ) return 1;

  }
  catch( const Exception& e )
  {
    std :: cerr << e.what() << std :: endl;
    return 1;
  }
  catch( ... )
  {
    return 1;
  }
  return 0;
}
