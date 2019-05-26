#ifndef DUNE_FEM_GRIDPART_COMMON_COMPOSITEGEOMETRY_HH
#define DUNE_FEM_GRIDPART_COMMON_COMPOSITEGEOMETRY_HH

#include <type_traits>
#include <utility>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>

#include <dune/grid/common/geometry.hh>

#include <dune/fem/common/fmatrixcol.hh>

namespace Dune
{

  // CompositeGeometry
  // -----------------

  template< class Geometry, class Embedding >
  struct CompositeGeometry
  {
    static const int mydimension = Embedding::mydimension;
    static const int coorddimension = Geometry::coorddimension;

    typedef typename Geometry::ctype ctype;
    typedef FieldVector< ctype, mydimension > LocalCoordinate;
    typedef FieldVector< ctype, coorddimension > GlobalCoordinate;

    typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;
    typedef FieldMatrix< ctype, coorddimension, mydimension > JacobianInverseTransposed;
    typedef JacobianInverseTransposed Jacobian;

    // Helper class to compute a matrix pseudo inverse
    typedef Impl::FieldMatrixHelper< ctype > MatrixHelper;

    CompositeGeometry ( Geometry geometry, Embedding embedding, int order )
      : geometry_( std::move( geometry ) ), embedding_( std::move( embedding ) ), order_( order )
    {}

    GeometryType type () const { return embedding_.type(); }

    int corners () const { return embedding_.corners(); }
    GlobalCoordinate corner( int i ) const { return geometry_.global( embedding_.corner( i ) ); }

    bool affine () const { return geometry_.affine() && embedding_.affine(); }

    GlobalCoordinate global ( const LocalCoordinate &local ) const { return geometry_.global( embedding_.global( local ) ); }
    LocalCoordinate local ( const GlobalCoordinate &global ) const { return embedding_.local( geometry_.local( global ) ); }

    JacobianTransposed jacobianTransposed ( const LocalCoordinate &local ) const
    {
      const FieldMatrix< ctype, mydimension, Embedding::coorddimension > jacEmbedding( embedding_.jacobianTransposed( local ) );
      const auto jacGeometry = geometry_.jacobianTransposed( embedding_.global( local ) );

      JacobianTransposed jacTransposed( 0 );
      for( int i = 0; i < mydimension; ++i )
        jacGeometry.mtv( jacEmbedding[ i ], jacTransposed[ i ] );
      return jacTransposed;
    }

    JacobianInverseTransposed jacobianInverseTransposed ( const LocalCoordinate &local ) const
    {
      JacobianInverseTransposed jacInverseTransposed( 0 );
      MatrixHelper::template rightInvA< mydimension, coorddimension >( jacobianTransposed( local ), jacInverseTransposed );
      return jacInverseTransposed;
    }

    ctype integrationElement ( const LocalCoordinate &local ) const
    {
      return MatrixHelper::template sqrtDetAAT< mydimension, coorddimension >( jacobianTransposed( local ) );
    }

    GlobalCoordinate center () const
    {
      GlobalCoordinate center( 0 );
      ctype volume( 0 );
      for( const auto &qp : QuadratureRules< ctype, mydimension >::rule( type(), order_+1 ) )
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
      for( const auto &qp : QuadratureRules< ctype, mydimension >::rule( type(), order_ ) )
        volume += qp.weight() * integrationElement( qp.position() );
      return volume;
    }

  private:
    Geometry geometry_;
    Embedding embedding_;
    int order_;
  };

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_GEOMETRY_HH
