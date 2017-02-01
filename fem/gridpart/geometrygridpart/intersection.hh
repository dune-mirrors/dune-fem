#ifndef DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_INTERSECTION_HH
#define DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_INTERSECTION_HH

#include <type_traits>

#include <dune/common/version.hh>
#include <dune/fem/misc/compatibility.hh>

#include <dune/geometry/affinegeometry.hh>

namespace Dune
{

  namespace Fem
  {

    // GeometryGridPartIntersection
    // ----------------------------

    template< class GridFamily >
    class GeometryGridPartIntersection
    {
      typedef GeometryGridPartIntersection<GridFamily> This;
      typedef typename std::remove_const< GridFamily >::type::Traits Traits;

    public:
      typedef typename std::remove_const< GridFamily >::type::ctype ctype;

      static const int dimension = std::remove_const< GridFamily >::type::dimension;
      static const int dimensionworld = std::remove_const< GridFamily >::type::dimensionworld;

      typedef typename Traits::template Codim< 0 >::Entity Entity;
#if ! DUNE_VERSION_NEWER( DUNE_GRID, 2, 4 )
      typedef typename Traits::template Codim< 0 >::EntityPointer EntityPointer;
#endif // #if ! DUNE_VERSION_NEWER( DUNE_GRID, 2, 4 )
      typedef typename Traits::template Codim< 0 >::Geometry ElementGeometry;
      typedef typename Traits::template Codim< 1 >::Geometry Geometry;
      typedef typename Traits::template Codim< 1 >::LocalGeometry LocalGeometry;

      typedef typename Traits::GridFunctionType GridFunctionType;

      typedef FieldVector< ctype, dimensionworld > GlobalCoordinate;
      typedef FieldVector< ctype, dimension-1 > LocalCoordinate;

    private:
#if ! DUNE_VERSION_NEWER( DUNE_GRID, 2, 4 )
      typedef typename EntityPointer::Implementation EntityPointerImplType;
#endif // #if ! DUNE_VERSION_NEWER( DUNE_GRID, 2, 4 )
      typedef typename Geometry::Implementation GeometryImplType;

      typedef typename Traits::HostGridPartType HostGridPartType;
      typedef typename HostGridPartType::IntersectionType HostIntersectionType;

    public:
      GeometryGridPartIntersection ( const GridFunctionType &gridFunction, const ElementGeometry &insideGeo )
        : gridFunction_( &gridFunction ), insideGeo_( insideGeo.impl() )
      {}

      operator bool () const { return bool( hostIntersection_ ); }

#if DUNE_VERSION_NEWER( DUNE_GRID, 2, 4 )
      Entity inside () const { return typename Entity::Implementation( gridFunction(), hostIntersection().inside() ); }
      Entity outside () const { return typename Entity::Implementation( gridFunction(), hostIntersection().outside() ); }
#else // #if DUNE_VERSION_NEWER( DUNE_GRID, 2, 4 )
      EntityPointer inside () const
      {
        return EntityPointerImplType( hostIntersection().inside(), gridFunction() );
      }

      EntityPointer outside () const
      {
        return EntityPointerImplType( hostIntersection().outside(), gridFunction() );
      }
#endif // #else // #if DUNE_VERSION_NEWER( DUNE_GRID, 2, 4 )

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
        return hostIntersection().geometryInInside();
      }

      LocalGeometry geometryInOutside () const
      {
        return hostIntersection().geometryInOutside();
      }

      Geometry geometry () const
      {
        typedef typename Geometry::Implementation Impl;
        return Geometry( Impl( insideGeo_, geometryInInside(), 2*insideGeo_.localFunction().order()+1 ) );
      }

      bool equals(const This &other) const
      {
        return hostIntersection() == other.hostIntersection();
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

      GlobalCoordinate integrationOuterNormal ( const LocalCoordinate &local ) const
      {
        const auto &refElement = ReferenceElements< ctype, dimension >::general( insideGeo_.type() );
        const auto &refNormal = refElement.integrationOuterNormal( indexInInside() );

        const auto jit = insideGeo_.jacobianInverseTransposed( geometryInInside().global( local ) );

        GlobalCoordinate normal;
        jit.mv( refNormal, normal );
        // double det = std::sqrt( GeometryGridPartGeometryType::MatrixHelper::template detATA<dimensionworld,dimension>( jit ) );
        // return normal *= ctype( 1 ) / sqrt(det);
        return normal *= geometry().integrationElement( local ) / normal.two_norm();

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
       return normal;
#endif
      }

#if 0
  void crossProduct(const FieldVector< ctype, dimensionworld > &vec1, const FieldVector< ctype, dimensionworld > &vec2, FieldVector< ctype, dimensionworld > &ret) const
  {
    assert( dimensionworld == 3 );

    ret[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
    ret[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
    ret[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
  }
#endif

      GlobalCoordinate outerNormal ( const LocalCoordinate &local ) const
      {
        const auto &refElement = Dune::ReferenceElements< ctype, dimension >::general( insideGeo_.type() );
        const auto &refNormal = refElement.integrationOuterNormal( indexInInside() );

        GlobalCoordinate normal;
        insideGeo_.jacobianInverseTransposed( geometryInInside().global( local ) ).mv( refNormal, normal );
        return normal;
      }

      GlobalCoordinate unitOuterNormal ( const LocalCoordinate &local ) const
      {
        GlobalCoordinate normal = outerNormal( local );
        return normal *= (ctype( 1 ) / normal.two_norm());
      }

      GlobalCoordinate centerUnitOuterNormal () const
      {
        const auto &refElement = Dune::ReferenceElements< ctype, dimension >::general( insideGeo_.type() );
        const auto &refNormal = refElement.integrationOuterNormal( indexInInside() );

        GlobalCoordinate normal;
        insideGeo_.jacobianInverseTransposed( geometryInInside().center() ).mv( refNormal, normal );
        return normal *= (ctype( 1 ) / normal.two_norm());
      }

      const HostIntersectionType &hostIntersection () const
      {
        assert( *this );
        return *hostIntersection_;
      }

      void invalidate () { hostIntersection_ = nullptr; }

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
      const HostIntersectionType *hostIntersection_ = nullptr;
      const GridFunctionType *gridFunction_ = nullptr;
      typename ElementGeometry::Implementation insideGeo_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_INTERSECTION_HH
