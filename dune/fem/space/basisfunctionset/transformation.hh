#ifndef DUNE_FEM_SPACE_BASISFUNCTIONSET_TRANSFORMATION_HH
#define DUNE_FEM_SPACE_BASISFUNCTIONSET_TRANSFORMATION_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/fem/common/explicitfieldvector.hh>
#include <dune/fem/common/fmatrixcol.hh>

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
      GeometryJacobianInverseTransposed gjit_;
    };



    // hessianTransformation
    // ---------------------

    template< class GeometryJacobianInverseTransposed, class K, int SIZE >
    void hessianTransformation ( const GeometryJacobianInverseTransposed &gjit,
                                 const ExplicitFieldVector< FieldMatrix< K, GeometryJacobianInverseTransposed::cols, GeometryJacobianInverseTransposed::cols >, SIZE > &a,
                                 ExplicitFieldVector< FieldMatrix< K, GeometryJacobianInverseTransposed::rows, GeometryJacobianInverseTransposed::rows >, SIZE > &b )
    {
      const int dimLocal = GeometryJacobianInverseTransposed::cols;
      const int dimGlobal = GeometryJacobianInverseTransposed::rows;

      for( int r = 0; r < SIZE; ++r )
      {
        // c = J^{-T} a_r^T
        FieldMatrix< K, dimLocal, dimGlobal > c;
        for( int i = 0; i < dimLocal; ++i )
          gjit.mv( a[ r ][ i ], c[ i ] );

        // b_r = J^{-T} c
        for( int i = 0; i < dimGlobal; ++i )
        {
          FieldMatrixColumn< const FieldMatrix< K, dimLocal, dimGlobal > > ci( c, i );
          FieldMatrixColumn< FieldMatrix< K, dimGlobal, dimGlobal > > bi( b[ r ], i );
          gjit.mv( ci, bi );
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
      GeometryJacobianInverseTransposed gjit_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_BASISFUNCTIONSET_TRANSFORMATION_HH
