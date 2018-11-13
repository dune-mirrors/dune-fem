# include "config.h"

// C++ includes
#include <string>

// dune-common includes
#include <dune/common/exceptions.hh>

// dune-fem includes
#include <dune/fem/io/parameter.hh>

#include <dune/fem/io/parameter/subreader.hh>


int main(int argc, char** argv)
try
{
  Dune::Fem::MPIManager::initialize( argc, argv );

  auto& parameter = Dune::Fem::Parameter::container();
  parameter.append( argc, argv );

  auto subParameter = subParameterReader( "level0." , parameter );
  auto subSubParameter = subParameterReader( "level0.", subParameterReader( "level1." , parameter ));

  std::string default_( "default" );

  auto value0 = parameter.getValue< std::string >( "level0.level1.value", default_ );
  auto value1 = subParameter.getValue< std::string >( "level1.value", default_ );
  auto value2 = subSubParameter.getValue< std::string >( "value", default_ );

  if( value0 == default_ )
    DUNE_THROW( Dune::Exception, "Value read is default." );

  if( value0 != value1 || value0 != value2 )
    DUNE_THROW( Dune::Exception, "Values read differ. [" << value0 << ", " << value1 << ", " << value2 << " ]" );

  return 0;
}
catch( const std::exception &e )
{
  std::cerr << "Exception: " << e.what() << std::endl;
  return 1;
}
