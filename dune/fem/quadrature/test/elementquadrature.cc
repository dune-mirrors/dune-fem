#include <config.h>

#include <dune/common/exceptions.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/quadrature/elementquadrature.hh>

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


template< int dim >
void checkQuadraturePoints ( Dune::GeometryType type, int order )
{
  std::cout << ">>> Checking ElementQuadrature for type " << type << " and order " << order << "..." << std::endl;

  Dune::Fem::ElementQuadrature< FakeGridPart< dim >, 0 > quadrature( type, order );
  const Dune::ReferenceElement< double, dim > &refElement = Dune::ReferenceElements< double, dim >::general( type );
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

  using Dune::GeometryType;
  for( int order = 0; order < 12; ++order )
    checkQuadraturePoints< 1 >( GeometryType( GeometryType::simplex, 1 ), order );
  for( int order = 0; order < 12; ++order )
    checkQuadraturePoints< 2 >( GeometryType( GeometryType::cube, 2 ), order );
  for( int order = 0; order < 19; ++order )
    checkQuadraturePoints< 2 >( GeometryType( GeometryType::simplex, 2 ), order );
  for( int order = 0; order < 12; ++order )
    checkQuadraturePoints< 3 >( GeometryType( GeometryType::cube, 3 ), order );
  for( int order = 0; order < 12; ++order )
    checkQuadraturePoints< 3 >( GeometryType( GeometryType::simplex, 3 ), order );
  for( int order = 0; order < 12; ++order )
    checkQuadraturePoints< 3 >( GeometryType( GeometryType::prism, 3 ), order );
  for( int order = 0; order < 12; ++order )
    checkQuadraturePoints< 3 >( GeometryType( GeometryType::pyramid, 3 ), order );

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
