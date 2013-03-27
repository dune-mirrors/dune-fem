#ifndef DUNE_FEM_BASEFUNCTIONSET_TRANSFORMATION_HH
#define DUNE_FEM_BASEFUNCTIONSET_TRANSFORMATION_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

namespace Dune
{

  namespace Fem
  {

    // jacobianTransformation
    // ----------------------


    template< class GeometryJacobianInverseTransposed, class K, int ROWS >
    void jacobianTransformation ( const GeometryJacobianInverseTransposed &gjit,
                                  const FieldMatrix< K, ROWS, GeometryJacobianInverseTransposed::cols > &a,
                                  FieldMatrix< K, ROWS, GeometryJacobianInverseTransposed::rows > &b )
    {
      for( int r = 0; r < ROWS; ++r )
        gjit.mv( a[ r ], b[ r ] );
    }



    // JacobianTransformation
    // ----------------------

    template< class Geometry >
    struct JacobianTransformation
    {
      typedef typename Geometry::LocalCoordinate LocalCoordinate;
      typedef typename Geometry::JacobianInverseTransposed GeometryJacobianInverseTransposed;

      JacobianTransformation ( const Geometry &geometry, const LocalCoordinate &x )
      : gjit_( geometry.jacobianInverseTransposed( x ) )
      {}

      template< class A, class B >
      void operator() ( const A &a, B &b ) const
      {
        jacobianTransformation( gjit_, a, b );
      }

    private:
      const GeometryJacobianInverseTransposed &gjit_;
    };



    // hessianTransformation
    // ---------------------

    template< class GeometryJacobianInverseTransposed, class K, int SIZE >
    void hessianTransformation ( const GeometryJacobianInverseTransposed &gjit,
                                 const FieldVector< FieldMatrix< K, GeometryJacobianInverseTransposed::cols, GeometryJacobianInverseTransposed::cols >, SIZE > &a,
                                 FieldVector< FieldMatrix< K, GeometryJacobianInverseTransposed::rows, GeometryJacobianInverseTransposed::rows >, SIZE > &b )
    {
      for( int r = 0; r < SIZE; ++r )
      {
        b[ r ] = K( 0 );
        for( int k = 0; k < GeometryJacobianInverseTransposed::cols; ++k )
        {
          FieldVector< K, GeometryJacobianInverseTransposed::rows > c;
          gjit.mv( a[ r ][ k ], c );
          for( int j = 0; j < GeometryJacobianInverseTransposed::rows; ++j )
            b[ r ][ j ].axpy( gjit[ j ][ k ], c );
        }
      }
    }



    // HessianTransformation
    // ---------------------

    template< class Geometry >
    struct HessianTransformation
    {
      typedef typename Geometry::LocalCoordinate LocalCoordinate;
      typedef typename Geometry::JacobianInverseTransposed GeometryJacobianInverseTransposed;

      HessianTransformation ( const Geometry &geometry, const LocalCoordinate &x )
      : gjit_( geometry.jacobianInverseTransposed( x ) )
      {
        if( !geometry.affine() )
          DUNE_THROW( NotImplemented, "HessianTransformation not implemented for non-affine geometries." );
      }

      template< class A, class B >
      void operator() ( const A &a, B &b ) const
      {
        hessianTransformation( gjit_, a, b );
      }

    private:
      const GeometryJacobianInverseTransposed &gjit_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_BASEFUNCTIONSET_TRANSFORMATION_HH
