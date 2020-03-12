#ifndef DUNE_FEM_GRIDPART_COMMON_SIMPLEGEOMETRY_HH
#define DUNE_FEM_GRIDPART_COMMON_SIMPLEGEOMETRY_HH

#include <type_traits>
#include <utility>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/geometry.hh>

#include <dune/fem/common/fmatrixcol.hh>

namespace Dune
{

  // SimpleGeometry
  // --------------

  template< class BasicGeometry >
  struct SimpleGeometry
    : public BasicGeometry
  {
    typedef typename BasicGeometry::ctype ctype;

    using BasicGeometry::mydimension;
    using BasicGeometry::coorddimension;

    using BasicGeometry::global;
    using BasicGeometry::jacobianTransposed;
    using BasicGeometry::quadrature;
    using BasicGeometry::type;

    typedef FieldVector< ctype, mydimension > LocalCoordinate;
    typedef FieldVector< ctype, coorddimension > GlobalCoordinate;

    typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;
    typedef FieldMatrix< ctype, coorddimension, mydimension > JacobianInverseTransposed;
    typedef JacobianInverseTransposed Jacobian;

    // Helper class to compute a matrix pseudo inverse
    typedef Impl::FieldMatrixHelper< ctype > MatrixHelper;

    template< class... Args, std::enable_if_t< std::is_constructible< BasicGeometry, Args &&... >::value, int > = 0 >
    SimpleGeometry ( Args &&...args )
      : BasicGeometry( std::forward< Args >( args )... )
    {}

    int corners () const { return referenceElement().size( mydimension ); }
    GlobalCoordinate corner ( int i ) const { return global( referenceElement().position( i, mydimension ) ); }

    bool affine () const { return false; }

    LocalCoordinate local ( const GlobalCoordinate &global ) const
    {
      const ctype tolerance = 1e-12; // use something better here e.g. Traits::tolerance();
      LocalCoordinate local = referenceElement().position( 0, 0 );
      LocalCoordinate dlocal;
      do
      {
        // Newton's method: DF^n dx^n = F^n, x^{n+1} -= dx^n
        const GlobalCoordinate dglobal = this->global( local ) - global;
        MatrixHelper::template xTRightInvA< mydimension, coorddimension >( jacobianTransposed( local ), dglobal, dlocal );
        local -= dlocal;
        assert( referenceElement().checkInside( local ) );
      }
      while( dlocal.two_norm2() > tolerance );
      return local;
    }

    ctype integrationElement ( const LocalCoordinate &local ) const
    {
      return MatrixHelper::template sqrtDetAAT< mydimension, coorddimension >( jacobianTransposed( local ) );
    }

    JacobianInverseTransposed jacobianInverseTransposed ( const LocalCoordinate &local ) const
    {
      JacobianInverseTransposed jacInverseTransposed( 0 );
      MatrixHelper::template rightInvA< mydimension, coorddimension >( jacobianTransposed( local ), jacInverseTransposed );
      return jacInverseTransposed;
    }

    GlobalCoordinate center () const
    {
      GlobalCoordinate center( 0 );
      ctype volume( 0 );
      for( const auto &qp : quadrature( 0 ) )
      {
        const ctype weight = qp.weight() * integrationElement( qp.position() );
        center.axpy( weight, global( qp.position() ) );
        volume += weight;
      }
      return center /= volume;
    }

    ctype volume () const
    {
      ctype volume( 0 );
      for( const auto &qp : quadrature( 0 ) )
        volume += qp.weight() * integrationElement( qp.position() );
      return volume;
    }

    auto referenceElement () const { return Dune::referenceElement<double>( type(), Dune::Dim<mydimension>() ); }
  };

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_COMMON_SIMPLEGEOMETRY_HH
