#ifndef DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_INTERSECTION_HH
#define DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_INTERSECTION_HH

#include <dune/common/version.hh>

#include <dune/fem/gridpart/geometrygridpart/geometry.hh>
#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/genericgeometry/matrixhelper.hh>

namespace Dune
{

  namespace Fem
  {

    // GeometryGridPartIntersection
    // ----------------------------

    template< class GridFamily >
    class GeometryGridPartIntersection
    {
      typedef typename remove_const< GridFamily >::type::Traits Traits;

    public:
      typedef typename remove_const< GridFamily >::type::ctype ctype;

      static const int dimension = remove_const< GridFamily >::type::dimension;
      static const int dimensionworld = remove_const< GridFamily >::type::dimensionworld;

      typedef typename Traits::template Codim< 0 >::Entity Entity;
#if ! DUNE_VERSION_NEWER( DUNE_GRID, 3, 0 )
      typedef typename Traits::template Codim< 0 >::EntityPointer EntityPointer;
#endif // #if ! DUNE_VERSION_NEWER( DUNE_GRID, 3, 0 )
      typedef typename Traits::template Codim< 0 >::Geometry ElementGeometry;
      typedef typename Traits::template Codim< 1 >::Geometry Geometry;
      typedef typename Traits::template Codim< 1 >::LocalGeometry LocalGeometry;

      typedef typename Traits::GridFunctionType GridFunctionType;
      typedef GeometryGridPartGeometry<1,dimensionworld,GridFamily> GeometryGridPartGeometryType;
      typedef Dune::AffineGeometry< ctype, 1, GridFamily::dimension > AffineGeometryType;

    private:
#if ! DUNE_VERSION_NEWER( DUNE_GRID, 3, 0 )
      typedef typename EntityPointer::Implementation EntityPointerImplType;
#endif // #if ! DUNE_VERSION_NEWER( DUNE_GRID, 3, 0 )
      typedef typename ElementGeometry::Implementation ElementGeometryImplType;
      typedef typename Geometry::Implementation GeometryImplType;

      typedef typename Traits::HostGridPartType HostGridPartType;
      typedef typename HostGridPartType::IntersectionType HostIntersectionType;

    public:
      GeometryGridPartIntersection ( const GridFunctionType &gridFunction, const ElementGeometry &insideGeo )
      : hostIntersection_( 0 ),
        gridFunction_( &gridFunction ),
        insideGeo_( insideGeo.impl() )
      {}

      operator bool () const { return bool( hostIntersection_ ); }

#if DUNE_VERSION_NEWER( DUNE_GRID, 3, 0 )
      Entity inside () const { return EntityImplType( hostIntersection().inside(), gridFunction() ); }
      Entity outside () const { return EntityImplType( hostIntersection().outside(), gridFunction() ); }
#else // #if DUNE_VERSION_NEWER( DUNE_GRID, 3, 0 )
      EntityPointer inside () const
      {
        return EntityPointerImplType( hostIntersection().inside(), gridFunction() );
      }

      EntityPointer outside () const
      {
        return EntityPointerImplType( hostIntersection().outside(), gridFunction() );
      }
#endif // #else // #if DUNE_VERSION_NEWER( DUNE_GRID, 3, 0 )

      bool boundary () const
      {
        return hostIntersection().boundary();
      }

      bool conforming () const
      {
        return hostIntersection().conforming();
      }

      int twistInSelf() const
      {
        return hostIntersection().impl().twistInSelf();
      }

      int twistInNeighbor() const
      {
        return hostIntersection().impl().twistInNeighbor();
      }

      bool neighbor () const
      {
        return hostIntersection().neighbor();
      }

      int boundaryId () const
      {
        return hostIntersection().boundaryId();
      }

      size_t boundarySegmentIndex () const
      {
        return hostIntersection().boundarySegmentIndex();
      }

      LocalGeometry geometryInInside () const
      {
        return LocalGeometry( hostIntersection().geometryInInside() );
      }

      LocalGeometry geometryInOutside () const
      {
        return LocalGeometry( hostIntersection().geometryInOutside() );
      }

      Geometry geometry () const
      {

        //! affineGeometry == hostIntersection().geometryInInside()
        AffineGeometryType affineGeometry ( type(), hostIntersection().geometryInInside().corner(0),
              hostIntersection().geometryInInside().jacobianTransposed(typename
                LocalGeometry::LocalCoordinate(0) ) );
        return Geometry( GeometryGridPartGeometryType(
              insideGeo_,
              gridFunction_, hostIntersection_, affineGeometry) );
      }

      GeometryType type () const
      {
        return hostIntersection().type();
      }

      int indexInInside () const
      {
        return hostIntersection().indexInInside();
      }

      int indexInOutside () const
      {
        return hostIntersection().indexInOutside();
      }

      FieldVector< ctype, dimensionworld >
      integrationOuterNormal ( const FieldVector< ctype, dimension-1 > &local ) const
      {
        const ReferenceElement< ctype, dimension > &refElement
          = ReferenceElements< ctype, dimension>::general( insideGeo_.type() );

        FieldVector< ctype, dimensionworld > normal;
        FieldVector< ctype, dimension > x( geometryInInside().global( local ) );

        const FieldMatrix< ctype, dimensionworld, dimension > &jit = insideGeo_.jacobianInverseTransposed( x );
        const FieldVector< ctype, dimension > &refNormal = refElement.volumeOuterNormal( indexInInside() );
        jit.mv( refNormal, normal );
        normal *= geometry().integrationElement(local)/normal.two_norm();
        // double det = std::sqrt( GeometryGridPartGeometryType::MatrixHelper::template detATA<dimensionworld,dimension>( jit ) );
        // normal *= ctype( 1 ) / sqrt(det);
        return normal;

#if 0
//! Alternative implementation (gives the same result)
FieldVector< ctype, dimension > x( geometryInInside().global( local ) );
FieldVector< ctype, dimensionworld > normal, tau, nu;
FieldMatrix< ctype, 1, dimensionworld >  tauT = geometry().jacobianTransposed(local);

for (int i = 0; i != dimensionworld; ++i)
  tau[i] = tauT[0][i];

tau /= tau.two_norm();

FieldMatrix< ctype, dimensionworld, dimensionworld >  localFnGrad, localFnGradT;
gridFunction().localFunction(*hostIntersection().inside()).jacobian(x, localFnGrad);

for (int i = 0; i != dimensionworld; ++i)
         for (int j = 0; j != dimensionworld; ++j )
           localFnGradT[i][j] = localFnGrad[j][i];

crossProduct(localFnGradT[0], localFnGradT[1], nu );
nu /= nu.two_norm();

/*
nu = insideGeo_.global(x);
nu /= nu.two_norm();
*/

crossProduct(nu, tau, normal);

// get conormal in same direction as discrete conormal
if (normal*hostIntersection().integrationOuterNormal(local) < 0)
  normal *= -1;

// normal /= normal.two_norm();
// normal *= GeometryImplType(insideGeo_, gridFunction_, hostIntersection_, affineGeometry).integrationElement(local);
normal *= geometry().integrationElement(local)/normal.two_norm();
       // return hostIntersection().integrationOuterNormal(local);
#endif
       return normal;
      }

  void crossProduct(const FieldVector< ctype, dimensionworld > &vec1, const FieldVector< ctype, dimensionworld > &vec2, FieldVector< ctype, dimensionworld > &ret) const
  {
    assert( dimensionworld == 3 );

    ret[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
    ret[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
    ret[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
  }

      FieldVector< ctype, dimensionworld >
      outerNormal ( const FieldVector< ctype, dimension-1 > &local ) const
      {

        // return hostIntersection().outerNormal( local );


        const Dune::ReferenceElement< ctype, dimension > &refElement
          = Dune::ReferenceElements< ctype, dimension >::general( insideGeo_.type() );


        FieldVector< ctype, dimension > x( geometryInInside().global( local ) );
        typedef typename ElementGeometryImplType::JacobianInverseTransposed JacobianInverseTransposed;
        const JacobianInverseTransposed &jit = insideGeo_.jacobianInverseTransposed( x );
        const FieldVector< ctype, dimension > &refNormal = refElement.integrationOuterNormal( indexInInside() );

        FieldVector< ctype, dimensionworld > normal;
        jit.mv( refNormal, normal );
        return normal;

      }

      FieldVector< ctype, dimensionworld >
      unitOuterNormal ( const FieldVector< ctype, dimension-1 > &local ) const
      {
        FieldVector< ctype, dimensionworld > normal = outerNormal( local );
        normal *= (ctype( 1 ) / normal.two_norm());
        return normal;

      }

      FieldVector< ctype, dimensionworld > centerUnitOuterNormal () const
      {
        const ReferenceElement< ctype, dimension-1 > &refFace
          = ReferenceElements< ctype, dimension-1 >::general( type() );
        return unitOuterNormal( refFace.position( 0, 0 ) );
      }

      const HostIntersectionType &hostIntersection () const
      {
        assert( *this );
        return *hostIntersection_;
      }

      void invalidate ()
      {
        hostIntersection_ = 0;
      }

      void initialize ( const HostIntersectionType &hostIntersection )
      {
        assert( !(*this) );
        hostIntersection_ = &hostIntersection;
      }

      const GridFunctionType &gridFunction () const
      {
        assert( gridFunction_ );
        return *gridFunction_;
      }

    private:
      const HostIntersectionType *hostIntersection_;
      const GridFunctionType *gridFunction_;
      ElementGeometryImplType insideGeo_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_INTERSECTION_HH
