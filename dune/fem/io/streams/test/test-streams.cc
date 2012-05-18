#include <config.h>

#include <dune/fem/version.hh>
#include <dune/fem/misc/mpimanager.hh>

#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/streams/asciistreams.hh>
#include <dune/fem/io/streams/binarystreams.hh>
#include <dune/fem/io/streams/datastreams.hh>
#include <dune/fem/io/streams/xdrstreams.hh>
#include <dune/fem/io/streams/sionlibstreams.hh>

using namespace Dune;

struct Data
{
  std :: string my_string;
  unsigned int my_uint;
  unsigned long my_ulong;
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
      << data.my_uint << data.my_ulong << data.my_double
      << data.my_int << data.my_float << data.my_char << data.my_bool;

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
     >> check.my_uint >> check.my_ulong >> check.my_double
     >> check.my_int >> check.my_float >> check.my_char >> check.my_bool;

  std :: cerr << "Data: " << check.my_string 
              << ", " << check.my_uint
              << ", " << check.my_ulong
              << ", " << check.my_double << ", " << check.my_int
              << ", " << check.my_float 
              << ", " << (int)check.my_char
              << ", " << check.my_bool
              << std :: endl;

  versionId = in.readUnsignedInt();
  if( versionId != DUNE_MODULE_VERSION_ID(DUNE_FEM) )
  {
    std :: cerr << "Incorrect versionId read back: "
                << versionId << " != " << DUNE_MODULE_VERSION_ID(DUNE_FEM)
                << std :: endl;
    return false;
  }

  bool equal = true;
  equal &= (data.my_string == check.my_string);
  equal &= (data.my_uint == check.my_uint);
  equal &= (data.my_ulong == check.my_ulong);
  equal &= (data.my_double == check.my_double);
  equal &= (data.my_int == check.my_int);
  equal &= (data.my_float == check.my_float);
  equal &= (data.my_char  == check.my_char );
  equal &= (data.my_bool  == check.my_bool);
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
    data.my_uint = 42;
    data.my_ulong = -4 ; // this results in 18446744073709551612 
    data.my_double = 1.2345678901234;
    data.my_int = -767;
    data.my_float = 1.23456;
    data.my_char  = 123;
    data.my_bool  = true;
    
    std::stringstream filestr;
    filestr << Dune :: Parameter :: commonOutputPath() << "/test." << MPIManager :: rank() << ".";
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
        return 1;
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
        return 1;
    }

    {
      std::string filename( filestr.str() + "xdr" );
      std :: cerr << "Checking XDR streams..." << std :: endl;
      if( writeStreams ) 
      {
        Fem :: XDRFileOutStream xout( filename.c_str() );
        write( xout, data );
        xout.flush();
      }
      Fem :: XDRFileInStream xin( filename.c_str() );
      if( !read( xin, data ) )
        return 1;
    }

#if HAVE_ALUGRID
    {
      std::string filename( filestr.str() + "data" );
      std :: cerr << "Checking Data streams..." << std :: endl;
      Fem :: DataOutStream dout( filename.c_str() );
      write( dout, data );
      dout.flush();
      Fem :: DataInStream din( filename.c_str() );
      if( ! read( din, data ) )
        return 1;  
    }
#endif

#if HAVE_SIONLIB
    {
      std :: cerr << "Checking SIONlib streams..." << std :: endl;
      std::stringstream file;
      file << Dune :: Parameter::commonOutputPath() << "/test.sion." << MPIManager :: size() ;
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
          return 1;  
      }

      if( MPIManager :: size () == 1 ) 
      {
        int size = (argc > 2) ? atoi( argv[ 2 ] ) : MPIManager :: size();
        std::cout << "Check serial read for " << size << " procs" << std::endl;
        std::stringstream file;
        file << Dune :: Parameter::commonOutputPath() << "/test.sion." << size ;
        // check serial read
        for( int rank=0; rank<size; ++rank ) 
        {
          Fem :: SIONlibInStream sionin( filename.c_str(), rank );
          if( ! read( sionin, data ) )
            return 1;  
        }
      }

    }
#endif

  }
  catch( Exception e )
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
