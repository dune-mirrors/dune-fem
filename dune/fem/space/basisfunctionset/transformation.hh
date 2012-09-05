#ifndef DUNE_FEM_BASEFUNCTIONSET_TRANSFORMATION_HH
#define DUNE_FEM_BASEFUNCTIONSET_TRANSFORMATION_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

namespace Dune
{

  namespace Fem
  {

    // JacobianTransformation
    // ----------------------

    template< class Geometry >
    struct JacobianTransformation
    {
      typedef typename Geometry::LocalCoordinate LocalCoordinate;
      typedef typename Geometry::Jacobian GeometryJacobianInverseTransposed;

      JacobianTransformation ( const Geometry &geometry, const LocalCoordinate &x )
      : gjit_( geometry.jacobianInverseTransposed( x ) )
      {}

      template< class K, int ROWS >
      void operator() ( const FieldMatrix< K, ROWS, Geometry::mydimension > &a,
                        FieldMatrix< K, ROWS, Geometry::coorddimension > &b )
      {
        for( int r = 0; r < ROWS; ++r )
          gjit_.mv( a[ r ], b[ r ] );
      }

    private:
      const GeometryJacobianInverseTransposed &gjit_;
    };



    // HessianTransformation
    // ---------------------

    template< class Geometry >
    struct HessianTransformation
    {
      typedef typename Geometry::LocalCoordinate LocalCoordinate;
      typedef typename Geometry::Jacobian GeometryJacobianInverseTransposed;

      HessianTransformation ( const Geometry &geometry, const LocalCoordinate &x )
      : gjit_( geometry.jacobianInverseTransposed( x ) )
      {
        if( !geometry.affine() )
          DUNE_THROW( NotImplemented, "HessianTransformation not implemented for non-affine geometries." );
      }

      template< class K, int SIZE >
      void operator() ( const FieldVector< FieldMatrix< K, Geometry::mydimension, Geometry::mydimension >, SIZE > &a,
                        FieldVector< FieldMatrix< K, Geometry::coorddimension, Geometry::coorddimension >, SIZE > &b )
      {
        for( int r = 0; r < SIZE; ++r )
        {
          b[ r ] = K( 0 );
          for( int k = 0; k < Geometry::mydimension; ++k )
          {
            FieldVector< K, Geometry::coorddimension > c;
            gjit_.mv( a[ r ][ k ], c );
            for( int j = 0; j < Geometry::coorddimension; ++j )
              b[ r ][ j ].axpy( gjit_[ j ][ k ], c );
          }
        }
      }

    private:
      const GeometryJacobianInverseTransposed &gjit_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_BASEFUNCTIONSET_TRANSFORMATION_HH
