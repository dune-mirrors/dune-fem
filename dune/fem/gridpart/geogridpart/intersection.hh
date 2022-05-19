#ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_INTERSECTION_HH
#define DUNE_FEM_GRIDPART_GEOGRIDPART_INTERSECTION_HH

#include <type_traits>
#include <utility>

#include <dune/fem/gridpart/geogridpart/cornerstorage.hh>

namespace Dune
{

  namespace Fem
  {

    // GeoIntersection
    // --------------

    template< class GridFamily >
    class GeoIntersection
    {
      typedef typename std::remove_const< GridFamily >::type::Traits Traits;

    public:
      typedef typename std::remove_const< GridFamily >::type::ctype ctype;

      static const int dimension = std::remove_const< GridFamily >::type::dimension;
      static const int dimensionworld = std::remove_const< GridFamily >::type::dimensionworld;

      typedef typename Traits::template Codim< 0 >::Entity Entity;
      typedef typename Traits::template Codim< 0 >::Geometry ElementGeometry;
      typedef typename Traits::template Codim< 1 >::Geometry Geometry;
      typedef typename Traits::template Codim< 1 >::LocalGeometry LocalGeometry;

      typedef typename Traits::CoordFunctionType CoordFunctionType;

    private:
      typedef typename Entity::Implementation EntityImplType;
      typedef typename ElementGeometry::Implementation ElementGeometryImplType;
      typedef typename Geometry::Implementation GeometryImplType;

      typedef typename Traits::HostGridPartType HostGridPartType;
      typedef typename HostGridPartType::IntersectionType HostIntersectionType;

      typedef GeoIntersectionCoordVector< GridFamily > CoordVectorType;

    public:
      GeoIntersection ( const CoordFunctionType &coordFunction, const ElementGeometry &insideGeo, HostIntersectionType hostIntersection )
      : coordFunction_( &coordFunction ),
        insideGeo_( insideGeo.impl() ),
        hostIntersection_( std::move( hostIntersection ) )
      {}

      GeoIntersection ()
      : coordFunction_( nullptr ),
        insideGeo_(),
        hostIntersection_()
      {}

      Entity inside () const
      {
        return Entity( EntityImplType( coordFunction(), hostIntersection().inside() ) );
      }

      Entity outside () const
      {
        return Entity( EntityImplType( coordFunction(), hostIntersection().outside() ) );
      }

      bool boundary () const
      {
        return hostIntersection().boundary();
      }

      bool conforming () const
      {
        return hostIntersection().conforming();
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
        const LocalGeometry &localGeo = geometryInInside();
        CoordVectorType coords( insideGeo_, localGeo );
        return Geometry( GeometryImplType( type(), coords ) );
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
        auto refElement = referenceElement< ctype, dimension>( insideGeo_.type() );

        FieldVector< ctype, dimension > x( geometryInInside().global( local ) );
        typedef typename ElementGeometryImplType::JacobianInverseTransposed JacobianInverseTransposed;
        const JacobianInverseTransposed &jit = insideGeo_.jacobianInverseTransposed( x );
        const FieldVector< ctype, dimension > &refNormal = refElement.integrationOuterNormal( indexInInside() );

        FieldVector< ctype, dimensionworld > normal;
        jit.mv( refNormal, normal );
        normal *= ctype( 1 ) / jit.det();
        return normal;
      }

      FieldVector< ctype, dimensionworld >
      outerNormal ( const FieldVector< ctype, dimension-1 > &local ) const
      {
        auto refElement = referenceElement< ctype, dimension>( insideGeo_.type() );

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
        auto refFace = referenceElement< ctype, dimension-1 >( type() );
        return unitOuterNormal( refFace.position( 0, 0 ) );
      }

      const CoordFunctionType &coordFunction () const
      {
        assert( coordFunction_ );
        return *coordFunction_;
      }

      const HostIntersectionType &hostIntersection () const
      {
        return hostIntersection_;
      }

    private:
      const CoordFunctionType *coordFunction_ = nullptr;
      ElementGeometryImplType insideGeo_;
      HostIntersectionType hostIntersection_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_INTERSECTION_HH
