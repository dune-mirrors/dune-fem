#include <config.h>

#include <dune/common/exceptions.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/quadrature/elementquadrature.hh>
#include <dune/fem/quadrature/dunequadratures.hh>
#include <dune/fem/quadrature/interpolationquadrature.hh>

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


template <class ctype, int dim>
ctype analyticalSolution (Dune::GeometryType t, int p, int direction )
{
  using Dune::GeometryType;
  ctype exact=0;

  if (t.isCube())
  {
    exact=1.0/(p+1);
    return exact;
  }

  if (t.isSimplex())
  {
    /* 1/(prod(k=1..dim,(p+k)) */
    exact = ctype( 1 );
    for( int k = 1; k <= dim; ++k )
      exact *= p+k;
    exact = ctype( 1 ) / exact;
    return exact;
  }

  if (t.isPrism())
  {
    const int pdim = (dim > 0 ? dim-1 : 0);
    if( direction < dim-1 )
    {
      Dune::GeometryType nt = Dune::GeometryTypes::simplex( dim-1 );
      if( dim > 0 )
        exact = analyticalSolution< ctype, pdim >( nt, p, direction );
      else
        exact = ctype( 1 );
    }
    else
      exact = ctype( 1 ) / ctype( Dune::factorial( pdim ) * (p+1));
    return exact;
  }

  if (t.isPyramid())
  {
    switch( direction )
    {
    case 0 :
    case 1 :
      exact=1.0/((p+3)*(p+1));
      break;
    case 2 :
      exact=2.0/((p+1)*(p+2)*(p+3));
      break;
    };
    return exact;
  }

  DUNE_THROW(Dune::NotImplemented, __func__ << " for " << t);
  return exact;
}


template<class Quadrature>
bool checkIntegration(const Quadrature &quad)
{
  typedef typename Quadrature::RealType ctype;
  const unsigned int dim = Quadrature::dimension;
  const unsigned int p = quad.order();
  const Dune::GeometryType& t = quad.geometry();
  typename Quadrature::CoordinateType integral(0);
  for (const auto& qp : quad )
  {
    // pos of integration point
    const auto& x = Dune::Fem::coordinate( qp );
    const ctype weight = qp.weight();

    for (unsigned int d=0; d<dim; d++)
      integral[d] += weight*std::pow(x[d],double(p));
  }

  ctype maxRelativeError = 0;
  for(unsigned int d=0; d<dim; d++)
  {
    ctype exact = analyticalSolution<ctype,dim>(t,p,d);
    ctype relativeError = std::abs(integral[d]-exact) /
                          (std::abs(integral[d])+std::abs(exact));
    if (relativeError > maxRelativeError)
      maxRelativeError = relativeError;
  }
  ctype epsilon = std::pow(2.0,double(p))*p*std::numeric_limits<double>::epsilon();
  if (p==0)
    epsilon = 2.0*std::numeric_limits<double>::epsilon();
  if (maxRelativeError > epsilon) {
    std::cerr << "Error: Quadrature for " << t << " and order=" << p << " failed" << std::endl;
    for (unsigned int d=0; d<dim; d++)
    {
      ctype exact = analyticalSolution<ctype,dim>(t,p,d);
      ctype relativeError = std::abs(integral[d]-exact) /
                            (std::abs(integral[d])+std::abs(exact));
      std::cerr << "       relative error " << relativeError << " in direction " << d << " (exact = " << exact << " numerical = " << integral[d] << ")" << std::endl;
    }
    return false;
  }
  return true;
}


template< int dim, template <class,int> class QuadTraits = Dune::Fem::DefaultQuadratureTraits >
void checkQuadraturePoints ( Dune::GeometryType type, int order )
{
  std::cout << ">>> Checking ElementQuadrature for type " << type << " and order " << order << "..." << std::endl;

  typedef typename Dune::Fem::ElementQuadrature< FakeGridPart< dim >, 0, QuadTraits > Quadrature;

  Quadrature quadrature( type, order );
  const auto refElement = Dune::referenceElement< double, dim >( type );
  for( std::size_t qp = 0; qp < quadrature.nop(); ++qp )
  {
    if( !refElement.checkInside( Dune::Fem::coordinate( quadrature[ qp ] ) ) )
      DUNE_THROW( Dune::RangeError, "Quadrature Point " << Dune::Fem::coordinate( quadrature[ qp ] ) << " not within reference element." );
  }

  double sum = 0 ;
  for( std::size_t qp = 0; qp < quadrature.nop(); ++qp )
  {
    sum += quadrature.weight( qp );
  }
  if( std::abs( sum - refElement.volume() ) > 1e-8 )
      DUNE_THROW( Dune::RangeError, "Weights don't sum up to 1");

  checkIntegration( quadrature );
}


int main ( int argc, char **argv )
try
{
  Dune::Fem::MPIManager::initialize( argc, argv );

  std::cout << "****  Checking DefaultQuadratureTraits  ****" << std::endl;
  for( int order = 0; order < 12; ++order )
    checkQuadraturePoints< 1 >( Dune::GeometryTypes::line, order );
  for( int order = 0; order < 12; ++order )
    checkQuadraturePoints< 2 >( Dune::GeometryTypes::quadrilateral, order );
  for( int order = 0; order < 11; ++order )
    checkQuadraturePoints< 2 >( Dune::GeometryTypes::triangle, order );
  for( int order = 0; order < 12; ++order )
    checkQuadraturePoints< 3 >( Dune::GeometryTypes::hexahedron, order );
  for( int order = 0; order < 9; ++order )
    checkQuadraturePoints< 3 >( Dune::GeometryTypes::tetrahedron, order );
  for( int order = 0; order < 11; ++order )
    checkQuadraturePoints< 3 >( Dune::GeometryTypes::prism, order );
  for( int order = 0; order < 3; ++order )
    checkQuadraturePoints< 3 >( Dune::GeometryTypes::pyramid, order );

  std::cout << "****  Checking DuneQuadratureTraits  ****" << std::endl;
  for( int order = 0; order < 12; ++order )
    checkQuadraturePoints< 1, Dune::Fem::DuneQuadratureTraits >( Dune::GeometryTypes::line, order );
  for( int order = 0; order < 12; ++order )
    checkQuadraturePoints< 2, Dune::Fem::DuneQuadratureTraits >( Dune::GeometryTypes::quadrilateral, order );
  for( int order = 0; order < 11; ++order )
    checkQuadraturePoints< 2, Dune::Fem::DuneQuadratureTraits >( Dune::GeometryTypes::triangle, order );
  for( int order = 0; order < 12; ++order )
    checkQuadraturePoints< 3, Dune::Fem::DuneQuadratureTraits >( Dune::GeometryTypes::hexahedron, order );
  for( int order = 0; order < 9; ++order )
    checkQuadraturePoints< 3, Dune::Fem::DuneQuadratureTraits >( Dune::GeometryTypes::tetrahedron, order );
  for( int order = 0; order < 11; ++order )
    checkQuadraturePoints< 3, Dune::Fem::DuneQuadratureTraits >( Dune::GeometryTypes::prism, order );
  for( int order = 0; order < 3; ++order )
    checkQuadraturePoints< 3, Dune::Fem::DuneQuadratureTraits >( Dune::GeometryTypes::pyramid, order );

#if HAVE_DUNE_LOCALFUNCTIONS
  std::cout << "****  Checking GaussLobattoQuadratureTraits  ****" << std::endl;
  for( int order = 0; order < 12; ++order )
    checkQuadraturePoints< 2, Dune::Fem::GaussLobattoQuadratureTraits >( Dune::GeometryTypes::quadrilateral, order );
  for( int order = 0; order < 12; ++order )
    checkQuadraturePoints< 3, Dune::Fem::GaussLobattoQuadratureTraits >( Dune::GeometryTypes::hexahedron, order );

  std::cout << "****  Checking GaussLobattoQuadratureTraits  ****" << std::endl;
  for( int order = 0; order < 12; ++order )
    checkQuadraturePoints< 2, Dune::Fem::GaussLegendreQuadratureTraits >( Dune::GeometryTypes::quadrilateral, order );
  for( int order = 0; order < 12; ++order )
    checkQuadraturePoints< 3, Dune::Fem::GaussLegendreQuadratureTraits >( Dune::GeometryTypes::hexahedron, order );
#endif
  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
