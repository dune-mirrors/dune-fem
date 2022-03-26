// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "config.h"

#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <dune/fem/misc/mpimanager.hh>

#include <dune/fem/quadrature/elementquadrature.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/quadrature/intersectionquadrature.hh>
#include <dune/fem/gridpart/leafgridpart.hh>

#include <dune/alugrid/grid.hh>
#include <dune/alugrid/dgf.hh>

static const int dim = 2;
using GridType = Dune::ALUGrid<dim, dim, Dune::cube, Dune::nonconforming>;

template <class Grid>
void checkPeriodic( Grid& grid )
{
  typedef Dune::Fem::LeafGridPart< Grid > GridPartType;
  GridPartType gp( grid );

  for(const auto& elem : Dune::elements( gp ) )
  {
    for( const auto& intersection : Dune::intersections( gp, elem ))
    {
      // check periodic boundary
      if( intersection.neighbor() && intersection.boundary() )
      {
        typedef Dune::Fem::CachingQuadrature< GridPartType, 1 > FaceQuadratureType;

        // use IntersectionQuadrature to create appropriate face quadratures
        typedef Dune::Fem::IntersectionQuadrature< FaceQuadratureType, true > IntersectionQuadratureType;
        typedef typename IntersectionQuadratureType :: FaceQuadratureType QuadratureImp;

        // create intersection quadrature (no neighbor check here)
        IntersectionQuadratureType interQuad( gp, intersection, 4, true );

        // get appropriate references
        const QuadratureImp &faceQuadInner = interQuad.inside();
        const QuadratureImp &faceQuadOuter = interQuad.outside();

        const auto& geomInside = intersection.inside().geometry();
        const auto& geomOutside = intersection.outside().geometry();

        const int nop = faceQuadInner.nop();
        assert( nop == int(faceQuadOuter.nop()) );
        for( int qp = 0; qp<nop; ++qp )
        {
          auto gIn  = geomInside.global ( faceQuadInner.point(qp ) );
          auto gOut = geomOutside.global ( faceQuadOuter.point(qp ) );

          auto diff = gIn - gOut;

          if( std::abs(diff.two_norm() - 1.0 ) > 1e-12 )
          {
            std::cout << gIn << " " << gOut << " " << diff << std::endl;
            DUNE_THROW(Dune::InvalidStateException,"Centers of periodic boundary differ");
          }
        }
      }
    }
  }
}

std::shared_ptr< GridType > createGrid( const double length, const int cells )
{
  std::stringstream file;
  file << "DGF" << std::endl;
  file << "Interval" << std::endl;
  for( int i=0; i<dim; ++i )
    file << "0 ";
  file << std::endl;
  for( int i=0; i<dim; ++i )
    file << length << " ";
  file << std::endl;
  for( int i=0; i<dim; ++i )
    file << cells << " ";
  file << std::endl;
  file << "#" << std::endl;
  file << "PERIODICFACETRANSFORMATION" << std::endl;
  file << length << " 0, 0 " << length << " + " << length << " 0" << std::endl;
  file << length << " 0, 0 " << length << " + " << "0 " << length << std::endl;
  file << "#" << std::endl;
  file << "GridParameter" << std::endl;
  file << "overlap 0" << std::endl;
  file << "#" << std::endl;

  std::shared_ptr< GridType > ptr;
  ptr.reset( Dune::GridPtr< GridType > (file).release() );
  return ptr;
}

int main( int argc, char** argv )
try
{
  Dune::Fem::MPIManager :: initialize( argc, argv );

  std::shared_ptr< GridType > gridPtr;

  if( argc > 1 )
    gridPtr.reset( Dune::GridPtr< GridType > (argv[1]).release() );
  else
    gridPtr = createGrid( 1.0, 4 );

  GridType& grid = *gridPtr;
  grid.loadBalance();

  checkPeriodic( grid );

  for( int i=0; i<3; ++i )
  {
    grid.globalRefine( 1 );
    checkPeriodic( grid );
  }

  return 0;
}
catch ( Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
catch (std::exception &e) {
  std::cerr << e.what() << std::endl;
  return 1;
}
catch ( ... )
{
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
