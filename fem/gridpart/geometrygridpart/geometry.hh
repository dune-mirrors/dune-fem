// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_GEOMETRY_HH
#define DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_GEOMETRY_HH

#include <type_traits>

#include <dune/common/fmatrix.hh>

#include <dune/geometry/affinegeometry.hh>

#include <dune/grid/common/geometry.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/quadrature/intersectionquadrature.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< int, int, class >
    class GeometryGridPartGeometry;



    // GeometryGridPartGeometryTraits
    // ------------------------------

    template< int mydim, class GridFamily >
    struct GeometryGridPartGeometryTraits
    {
      typedef typename std::remove_const< GridFamily >::type::Traits::HostGridPartType HostGridPartType;

      static const int dimension = HostGridPartType::dimension;
      static const int mydimension = mydim;
      static const int codimension = dimension - mydimension;
      static const int dimensionworld = GridFamily::dimensionworld;

      typedef typename HostGridPartType::template Codim< codimension >::GeometryType HostGeometryType;
    };



    // GeometryGridPartBasicGeometry
    // -----------------------------

    template< int mydim, int cdim, class GridFamily, class = void >
    struct GeometryGridPartBasicGeometry
    {
      typedef GeometryGridPartGeometryTraits< mydim, GridFamily > Traits;
      typedef typename Traits::HostGeometryType HostGeometryType;

      static const int dimension = HostGeometryType::dimension;
      static const int mydimension = HostGeometryType::mydimension;
      static const int coorddimension = Traits::dimensionworld;

      typedef typename HostGeometryType::ctype ctype;
      typedef FieldVector< ctype, mydimension > LocalVector;
      typedef FieldVector< ctype, coorddimension > GlobalVector;

      typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;
      typedef FieldMatrix< ctype, coorddimension, mydimension > JacobianInverseTransposed;

      typedef typename GridFamily::GridFunctionType GridFunctionType;

      GeometryType type () const
      {
        DUNE_THROW( NotImplemented, "GeometryGridPart does not implement Geometry for codimension " << (GridFamily::dimension-mydim) );
      }

      bool affine () const { return false; }

      int corners () const
      {
        DUNE_THROW( NotImplemented, "GeometryGridPart does not implement Geometry for codimension " << (GridFamily::dimension-mydim) );
      }

      GlobalVector corner ( const int i ) const
      {
        DUNE_THROW( NotImplemented, "GeometryGridPart does not implement Geometry for codimension " << (GridFamily::dimension-mydim) );
      }

      GlobalVector center () const
      {
        DUNE_THROW( NotImplemented, "GeometryGridPart does not implement Geometry for codimension " << (GridFamily::dimension-mydim) );
      }

      GlobalVector global ( const LocalVector &local ) const
      {
        DUNE_THROW( NotImplemented, "GeometryGridPart does not implement Geometry for codimension " << (GridFamily::dimension-mydim) );
      }

      LocalVector local ( const GlobalVector &global ) const
      {
        DUNE_THROW( NotImplemented, "GeometryGridPart does not implement Geometry for codimension " << (GridFamily::dimension-mydim) );
      }

      ctype integrationElement ( const LocalVector &local ) const
      {
        DUNE_THROW( NotImplemented, "GeometryGridPart does not implement Geometry for codimension " << (GridFamily::dimension-mydim) );
      }

      ctype volume () const
      {
        DUNE_THROW( NotImplemented, "GeometryGridPart does not implement Geometry for codimension " << (GridFamily::dimension-mydim) );
      }

      FieldMatrix< ctype, mydimension, coorddimension > jacobianTransposed ( const LocalVector &local ) const
      {
        DUNE_THROW( NotImplemented, "GeometryGridPart does not implement Geometry for codimension " << (GridFamily::dimension-mydim) );
      }

      FieldMatrix< ctype, coorddimension, mydimension > jacobianInverseTransposed ( const LocalVector &local ) const
      {
        DUNE_THROW( NotImplemented, "GeometryGridPart does not implement Geometry for codimension " << (GridFamily::dimension-mydim) );
      }
    };

    template< int mydim, int cdim, class GridFamily >
    struct GeometryGridPartBasicGeometry< mydim, cdim, GridFamily, std::enable_if_t< mydim == GridFamily::dimension-1 > >
    {
      typedef GeometryGridPartGeometryTraits< mydim, GridFamily > Traits;
      typedef typename Traits::HostGeometryType HostGeometryType;

      static const int dimension = HostGeometryType::dimension;
      static const int mydimension = HostGeometryType::mydimension;
      static const int coorddimension = Traits::dimensionworld;
      // static const int dimensionworld = HostGeometryType::dimensionworld;

      typedef typename HostGeometryType::ctype ctype;
      typedef FieldVector< ctype, mydimension > LocalVector;
      typedef FieldVector< ctype, coorddimension > GlobalVector;

      typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;
      typedef FieldMatrix< ctype, coorddimension, mydimension > JacobianInverseTransposed;
      typedef JacobianInverseTransposed Jacobian;

      typedef typename GridFamily::GridFunctionType GridFunctionType;

      typedef typename Traits::HostGridPartType HostGridPartType;
      typedef Dune::Fem::CachingQuadrature< HostGridPartType, 0 > ElementQuadratureType;

      typedef Dune::AffineGeometry< ctype, mydim, GridFamily::dimension > AffineGeometryType;
      typedef typename HostGridPartType::IntersectionType HostIntersectionType;

      // Helper class to compute a matrix pseudo inverse
      typedef Impl::FieldMatrixHelper< ctype > MatrixHelper;

      GeometryGridPartBasicGeometry ( const GeometryGridPartBasicGeometry< GridFamily::dimension, cdim, GridFamily > elementGeo,
                                      const GridFunctionType *gridFunction,
                                      const HostIntersectionType *hostIntersection,
                                      const AffineGeometryType &affineGeometry )
        : elementGeo_( elementGeo ),
          gridFunction_( *gridFunction ),
          hostIntersection_( hostIntersection ),
          affineGeo_( affineGeometry )
      {}

      GeometryType type () const { return affineGeo_.type(); }
      bool affine () const { return false; }

      int corners () const { return affineGeo_.corners(); }
      GlobalVector corner ( const int i ) const { return elementGeo_.global( affineGeo_.corner( i ) ); }
      GlobalVector center () const { return elementGeo_.global( affineGeo_.center() ); }

      GlobalVector global ( const LocalVector &local ) const { return elementGeo_.global( affineGeo_.global( local ) ); }

      LocalVector local ( const GlobalVector &global ) const
      {
        const ctype tolerance = 1e-12; // use something better here e.g. Traits::tolerance();
        LocalVector x( 0 ); // do something better e.g. the center refElement().position( 0, 0 );
        LocalVector dx;
        do
        {
          // Newton's method: DF^n dx^n = F^n, x^{n+1} -= dx^n
          const GlobalVector dglobal = (*this).global( x ) - global;
          MatrixHelper::template xTRightInvA< mydimension, coorddimension >( jacobianTransposed( x ), dglobal, dx );
          x -= dx;
          /// assert( refElement().checkInside( x ) );
        } while( dx.two_norm2() > tolerance );
        return x;
      }

      ctype integrationElement ( const LocalVector &local ) const
      {
        FieldMatrix< ctype, mydimension, GridFamily::dimension > gradFeT( affineGeo_.jacobianTransposed( local ) );
        FieldMatrix< ctype, GridFamily::dimension, coorddimension > gradFT( elementGeo_.jacobianTransposed( affineGeo_.global( local )) );

        FieldMatrix< ctype, mydimension, coorddimension > gradFeTGradFT;
        //   = jacobianTransposed(local);
        Dune::FMatrixHelp::multMatrix( gradFeT, gradFT, gradFeTGradFT );
        return MatrixHelper::template sqrtDetAAT< mydimension, coorddimension >( gradFeTGradFT );
      }

      ctype volume () const
      {
        IntersectionQuadrature< CachingQuadrature< HostGridPartType, 1 >, true >  faceQuadInside( gridFunction_.gridPart(), *hostIntersection_, 2*gridFunction_.space().order() + 1 );

        ctype vol = 0.0;
        const int numQuadraturePoints = faceQuadInside.nop();
        for( int qp = 0; qp < numQuadraturePoints; ++qp )
        {
          const double weight = faceQuadInside.weight( qp );
          vol += weight * integrationElement( faceQuadInside.localPoint( qp ) );
        }

        return vol;
      }

      const FieldMatrix< ctype, mydimension, coorddimension > &
      jacobianTransposed ( const LocalVector &local ) const
      {
        FieldMatrix< ctype, mydimension, GridFamily::dimension > gradFeT( affineGeo_.jacobianTransposed( local ) );
        FieldMatrix< ctype, GridFamily::dimension, coorddimension > gradFT( elementGeo_.jacobianTransposed( affineGeo_.global( local )) );

        Dune::FMatrixHelp::multMatrix( gradFeT, gradFT, jacTransposed_ );

        return jacTransposed_;
      }

      const FieldMatrix< ctype, coorddimension, mydimension > &
      jacobianInverseTransposed ( const LocalVector &local ) const
      {
        MatrixHelper::template rightInvA< mydimension, coorddimension >( jacobianTransposed( local ), jacInverseTransposed_ );
        return jacInverseTransposed_;
      }

    private:
      const GeometryGridPartBasicGeometry< GridFamily::dimension, cdim, GridFamily > elementGeo_;
      const GridFunctionType gridFunction_;
      const HostIntersectionType *hostIntersection_;
      const AffineGeometryType affineGeo_;
      mutable FieldMatrix< ctype, mydimension, coorddimension > jacTransposed_;
      mutable FieldMatrix< ctype, coorddimension, mydimension > jacInverseTransposed_;
    };

    template< int mydim, int cdim, class GridFamily >
    struct GeometryGridPartBasicGeometry< mydim, cdim, GridFamily, std::enable_if_t< mydim == GridFamily::dimension > >
    {
      typedef GeometryGridPartGeometryTraits< mydim, GridFamily > Traits;

      typedef typename Traits::HostGeometryType HostGeometryType;

      static const int dimension = HostGeometryType::dimension;


      static const int mydimension = HostGeometryType::mydimension;
      static const int coorddimension = Traits::dimensionworld;

      typedef typename HostGeometryType::ctype ctype;
      typedef FieldVector< ctype, mydimension > LocalVector;
      typedef FieldVector< ctype, coorddimension > GlobalVector;

      typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;
      typedef FieldMatrix< ctype, coorddimension, mydimension > JacobianInverseTransposed;
      typedef JacobianInverseTransposed Jacobian;

      typedef typename GridFamily::GridFunctionType GridFunctionType;
      typedef typename GridFunctionType::LocalFunctionType LocalFunctionType;

      typedef typename Traits::HostGridPartType HostGridPartType;
      typedef Dune::Fem::CachingQuadrature< HostGridPartType, 0 > ElementQuadratureType;
      typedef typename HostGridPartType::template Codim< 0 >::EntityType HostEntityType;

      typedef typename HostGridPartType::IntersectionType HostIntersectionType;

      // Helper class to compute a matrix pseudo inverse
      typedef Impl::FieldMatrixHelper< ctype > MatrixHelper;

      GeometryGridPartBasicGeometry ( const HostGeometryType &hostGeometry, const HostEntityType &hostEntity, const GridFunctionType *gridFunction )
        : hostGeometry_( hostGeometry ),
          localFunction_( *gridFunction )
      {
        localFunction_.init( hostEntity );
      }

      GeometryGridPartBasicGeometry ( const HostEntityType &hostEntity, const GridFunctionType *gridFunction )
        : hostGeometry_( hostEntity.geometry() ),
          localFunction_( *gridFunction )
      {
        localFunction_.init( hostEntity );
      }

      GeometryType type () const { return hostGeometry_.type(); }
      bool affine () const { return false; }

      int corners () const { return hostGeometry_.corners(); }

      GlobalVector corner ( const int i ) const
      {
        GlobalVector ret;
        localFunction_.evaluate( hostGeometry_.local( hostGeometry_.corner( i ) ), ret );
        return ret;
      }

      GlobalVector center () const
      {
        GlobalVector ret;
        localFunction_.evaluate( hostGeometry_.local( hostGeometry_.center() ), ret );
        return ret;
      }

      GlobalVector global ( const LocalVector &local ) const
      {
        GlobalVector ret;
        localFunction_.evaluate( local, ret );
        return ret;
      }

      LocalVector local ( const GlobalVector &global ) const
      {
        std::cout << "Need to implement local() method" << std::endl;
        return hostGeometry_.local( global );
      }

      ctype integrationElement ( const LocalVector &local ) const
      {
        FieldMatrix< ctype, coorddimension, coorddimension > gradPhi;
        FieldMatrix< ctype, mydimension, coorddimension > gradFT( hostGeometry_.jacobianTransposed( local ) );
        FieldMatrix< ctype, coorddimension, mydimension > gradF;

        for( int i = 0; i != coorddimension; ++i )
          for( int j = 0; j != mydimension; ++j )
            gradF[ i ][ j ] = gradFT[ j ][ i ];

        FieldMatrix< ctype, coorddimension, mydimension > gradPhiGradF;
        FieldMatrix< ctype, mydimension, mydimension > retMatrix;

        localFunction_.jacobian( local, gradPhi );

        Dune::FMatrixHelp::multMatrix( gradPhi, gradF, gradPhiGradF );

        Dune::FMatrixHelp::multTransposedMatrix( gradPhiGradF, retMatrix );

        return std::sqrt( retMatrix.determinant() );
      }

      ctype volume () const
      {
        const HostEntityType &entity = localFunction_.entity();
        ElementQuadratureType quad( entity, 2*localFunction_.order() + 1 );

        ctype vol = 0.0;
        const int numQuadraturePoints = quad.nop();
        for( int qp = 0; qp < numQuadraturePoints; ++qp )
        {
          const double weight = quad.weight( qp );
          vol += weight * integrationElement( quad.point( qp ) );
        }

        return vol;
      }

      const FieldMatrix< ctype, mydimension, coorddimension > &
      jacobianTransposed ( const LocalVector &local ) const
      {
        FieldMatrix< ctype, mydimension, coorddimension > gradFT( hostGeometry_.jacobianTransposed( local ) );
        FieldMatrix< ctype, coorddimension, coorddimension > gradPhi, gradPhiT;

        localFunction_.jacobian( local, gradPhi );

        for( int i = 0; i != coorddimension; ++i )
          for( int j = 0; j != coorddimension; ++j )
            gradPhiT[ i ][ j ] = gradPhi[ j ][ i ];

        Dune::FMatrixHelp::multMatrix( gradFT, gradPhiT, jacTransposed_ );

        return jacTransposed_;
      }

      const FieldMatrix< ctype, coorddimension, mydimension > &
      jacobianInverseTransposed ( const LocalVector &local ) const
      {
        MatrixHelper::template rightInvA< mydimension, coorddimension >( jacobianTransposed( local ), jacInverseTransposed_ );

        return jacInverseTransposed_;
      }

    private:
      HostGeometryType hostGeometry_;
      LocalFunctionType localFunction_;
      mutable FieldMatrix< ctype, mydimension, coorddimension > jacTransposed_;
      mutable FieldMatrix< ctype, coorddimension, mydimension > jacInverseTransposed_;
    };



    // GeometryGridPartGeometry
    // ------------------------

    template< int mydim, int cdim, class GridFamily >
    struct GeometryGridPartGeometry
      : public GeometryGridPartBasicGeometry< mydim, cdim, GridFamily >
    {
      typedef GeometryGridPartBasicGeometry< mydim, cdim, GridFamily > Base;
      typedef typename Base::HostGeometryType HostGeometryType;
      typedef typename HostGeometryType::ctype ctype;

      typedef typename Base::GridFunctionType GridFunctionType;

      typedef typename std::remove_const< GridFamily >::type::Traits::HostGridPartType HostGridPartType;

      typedef typename HostGridPartType::template Codim< 0 >::EntityType HostEntityType;
      typedef typename HostGridPartType::IntersectionType HostIntersectionType;

      typedef Dune::AffineGeometry< ctype, mydim, GridFamily::dimension > AffineGeometryType;

      GeometryGridPartGeometry () = default;

      GeometryGridPartGeometry ( const GeometryGridPartBasicGeometry< GridFamily::dimension, cdim, GridFamily > elementGeo,
                                 const GridFunctionType *gridFunction,
                                 const HostIntersectionType *hostIntersection,
                                 const AffineGeometryType &affineGeometry )
        : Base( elementGeo, gridFunction, hostIntersection, affineGeometry )
      {}

      GeometryGridPartGeometry ( const HostGeometryType &hostGeometry,
                                 const HostEntityType &hostEntity,
                                 const GridFunctionType *gridFunction )
        : Base( hostGeometry, hostEntity, gridFunction )
      {}

      GeometryGridPartGeometry ( const HostEntityType &hostEntity, const GridFunctionType *gridFunction )
        : Base( hostEntity, gridFunction )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_GEOMETRY_HH
