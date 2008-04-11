#include <config.h>

#include <dune/fem/version.hh>
#include <dune/fem/io/streams/asciistreams.hh>
#include <dune/fem/io/streams/xdrstreams.hh>

using namespace Dune;

struct Data
{
  std :: string my_string;
  unsigned int my_uint;
  double my_double;
  int my_int;
  float my_float;
};

template< class Traits >
void write ( OutStreamInterface< Traits > &out, const Data &data )
{
  out << DuneFEM :: versionId();
  out << DuneFEM :: version();
  
  out << data.my_string << data.my_uint << data.my_double
      << data.my_int << data.my_float;

  out << DuneFEM :: versionId();
}

template< class Traits >
bool read ( InStreamInterface< Traits > &in, const Data &data )
{
  unsigned int versionId = in.readUnsignedInt();
  if( versionId != DuneFEM :: versionId() )
  {
    std :: cerr << "Incorrect versionId read back: "
                << versionId << " != " << DuneFEM :: versionId()
                << std :: endl;
    return false;
  }

  std :: string version;
  in >> version;
  if( version != DuneFEM :: version() )
  {
    std :: cerr << "Incorrect version read back: "
                << version << " != " << DuneFEM :: version()
                << std :: endl;
    return false;
  }

  Data check;
  in >> check.my_string >> check.my_uint >> check.my_double
     >> check.my_int >> check.my_float;

  std :: cerr << "Data: " << check.my_string << ", " << check.my_uint
              << ", " << check.my_double << ", " << check.my_int
              << ", " << check.my_float << std :: endl;

  if( in.readUnsignedInt() != DuneFEM :: versionId() )
  {
    std :: cerr << "Incorrect versionId read back: "
                << versionId << " != " << DuneFEM :: versionId()
                << std :: endl;
    return false;
  }

  bool equal = true;
  equal &= (data.my_string == check.my_string);
  equal &= (data.my_uint == check.my_uint);
  equal &= (data.my_double == check.my_double);
  equal &= (data.my_int == check.my_int);
  equal &= (data.my_float == check.my_float);
  return equal;
}

int main ()
{
  try
  {
    Data data;
    data.my_string = "Hello, World!";
    data.my_uint = 42;
    data.my_double = 1.2345678901234;
    data.my_int = -767;
    data.my_float = 1.23456;
    
    std :: cerr << "Checking ASCII streams..." << std :: endl;
    ASCIIOutStream aout( "test.ascii" );
    write( aout, data );
    aout.flush();
    ASCIIInStream ain( "test.ascii" );
    if( !read( ain, data ) )
      return 1;

    std :: cerr << "Checking XDR streams..." << std :: endl;
    XDRFileOutStream xout( "test.xdr" );
    write( xout, data );
    xout.flush();
    XDRFileInStream xin( "test.xdr" );
    if( !read( xin, data ) )
      return 1;
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
