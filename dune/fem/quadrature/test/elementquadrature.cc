#include <config.h>

#include <dune/common/exceptions.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/quadrature/elementquadrature.hh>
#include <dune/fem/quadrature/dunequadratures.hh>

template< int dim >
struct FakeGridPart
{
  typedef double ctype;
  static const int dimension = dim;

  template< int codim >
  struct Codim
  {
    typedef int EntityType;
  };
};


template< int dim, bool useDuneQuad >
void checkQuadraturePoints ( Dune::GeometryType type, int order )
{
  const char* duneQuad = useDuneQuad ? "DUNE " : " " ;
  std::cout << ">>> Checking " << duneQuad << "ElementQuadrature for type " << type << " and order " << order << "..." << std::endl;

  typedef typename std::conditional< useDuneQuad,
      Dune::Fem::ElementQuadrature< FakeGridPart< dim >, 0, Dune::Fem::DuneQuadratureTraits >,
      Dune::Fem::ElementQuadrature< FakeGridPart< dim >, 0 > > :: type  Quadrature;

  Quadrature quadrature( type, order );
  const auto refElement = Dune::referenceElement< double, dim >( type );
  for( std::size_t qp = 0; qp < quadrature.nop(); ++qp )
  {
    if( !refElement.checkInside( Dune::Fem::coordinate( quadrature[ qp ] ) ) )
      DUNE_THROW( Dune::RangeError, "Quadrature Point " << Dune::Fem::coordinate( quadrature[ qp ] ) << " not within reference element." );
  }
}


int main ( int argc, char **argv )
try
{
  Dune::Fem::MPIManager::initialize( argc, argv );

  for( int order = 0; order < 12; ++order )
    checkQuadraturePoints< 1, false >( Dune::GeometryTypes::line, order );
  for( int order = 0; order < 12; ++order )
    checkQuadraturePoints< 2, false >( Dune::GeometryTypes::quadrilateral, order );
  for( int order = 0; order < 11; ++order )
    checkQuadraturePoints< 2, false >( Dune::GeometryTypes::triangle, order );
  for( int order = 0; order < 12; ++order )
    checkQuadraturePoints< 3, false >( Dune::GeometryTypes::hexahedron, order );
  for( int order = 0; order < 9; ++order )
    checkQuadraturePoints< 3, false >( Dune::GeometryTypes::tetrahedron, order );
  for( int order = 0; order < 11; ++order )
    checkQuadraturePoints< 3, false >( Dune::GeometryTypes::prism, order );
  for( int order = 0; order < 3; ++order )
    checkQuadraturePoints< 3, false >( Dune::GeometryTypes::pyramid, order );

  for( int order = 0; order < 12; ++order )
    checkQuadraturePoints< 1, true >( Dune::GeometryTypes::line, order );
  for( int order = 0; order < 12; ++order )
    checkQuadraturePoints< 2, true >( Dune::GeometryTypes::quadrilateral, order );
  for( int order = 0; order < 11; ++order )
    checkQuadraturePoints< 2, true >( Dune::GeometryTypes::triangle, order );
  for( int order = 0; order < 12; ++order )
    checkQuadraturePoints< 3, true >( Dune::GeometryTypes::hexahedron, order );
  for( int order = 0; order < 9; ++order )
    checkQuadraturePoints< 3, true >( Dune::GeometryTypes::tetrahedron, order );
  for( int order = 0; order < 11; ++order )
    checkQuadraturePoints< 3, true >( Dune::GeometryTypes::prism, order );
  for( int order = 0; order < 3; ++order )
    checkQuadraturePoints< 3, true >( Dune::GeometryTypes::pyramid, order );

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
