#ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_INTERSECTION_HH
#define DUNE_FEM_GRIDPART_GEOGRIDPART_INTERSECTION_HH

#include <dune/fem/gridpart/geogridpart/cornerstorage.hh>
#include <dune/fem/gridpart/geogridpart/entitypointer.hh>

namespace Dune
{

  namespace Fem
  {

    // GeoIntersection
    // --------------

    template< class GridFamily >
    class GeoIntersection
    {
      typedef typename remove_const< GridFamily >::type::Traits Traits;

    public:
      typedef typename remove_const< GridFamily >::type::ctype ctype;
      
      static const int dimension = remove_const< GridFamily >::type::dimension;
      static const int dimensionworld = remove_const< GridFamily >::type::dimensionworld;

      typedef typename Traits::template Codim< 0 >::Entity Entity;
      typedef typename Traits::template Codim< 0 >::EntityPointer EntityPointer;
      typedef typename Traits::template Codim< 0 >::Geometry ElementGeometry;
      typedef typename Traits::template Codim< 1 >::Geometry Geometry;
      typedef typename Traits::template Codim< 1 >::LocalGeometry LocalGeometry;

      typedef typename Traits::CoordFunctionType CoordFunctionType;

    private:
      typedef typename EntityPointer::Implementation EntityPointerImpl;
      typedef typename Geometry::Implementation GeometryImpl;

      typedef typename Traits::HostGridPartType HostGridPartType;
      typedef typename HostGridPartType::IntersectionType HostIntersectionType;

      typedef GeoIntersectionCoordVector< GridFamily > CoordVectorType;

    public:
      GeoIntersection ( const CoordFunctionType &coordFunction, const ElementGeometry &insideGeo )
      : coordFunction_( &coordFunction ),
        insideGeo_( insideGeo.impl() ),
        hostIntersection_( 0 ),
        geo_( GeometryImpl() )
      {}

      GeoIntersection ( const GeoIntersection &other )
      : coordFunction_( other.coordFunction_ ),
        insideGeo_( other.insideGeo_.impl() ),
        hostIntersection_( other.hostIntersection_ ),
        geo_( other.geo_.impl() )
      {}

      const GeoIntersection &operator= ( const GeoIntersection &other )
      {
        coordFunction_ = other.coordFunction_;
        insideGeo_.impl() = other.insideGeo_.impl();
        invalidate();
        return *this;
      }

      operator bool () const { return bool( hostIntersection_ ); }

      EntityPointer inside () const
      {
        return EntityPointerImpl( coordFunction(), hostIntersection().inside() );
      }
      
      EntityPointer outside () const
      {
        return EntityPointerImpl( coordFunction(), hostIntersection().outside() );
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
          
      const LocalGeometry &geometryInInside () const
      {
        return hostIntersection().geometryInInside();
      }
      
      const LocalGeometry &geometryInOutside () const
      {
        return hostIntersection().geometryInOutside();
      }
     
      const Geometry &geometry () const
      {
        GeometryImpl &geo = geo_.impl();
        if( !geo )
        {
          const LocalGeometry &localGeo = geometryInInside();
          CoordVectorType coords( insideGeometry(), localGeo );
          geo = GeometryImpl( type(), coords );
        }
        return geo_;
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
        const ElementGeometry &geo = insideGeometry();

        const GenericReferenceElement< ctype, dimension > &refElement
          = GenericReferenceElements< ctype, dimension>::general( geo.type() );

        FieldVector< ctype, dimension > x( geometryInInside().global( local ) );
        typedef typename ElementGeometry::Implementation::JacobianInverseTransposed JacobianInverseTransposed;
        const JacobianInverseTransposed &jit = geo.impl().jacobianInverseTransposed( x );
        const FieldVector< ctype, dimension > &refNormal = refElement.volumeOuterNormal( indexInInside() );

        FieldVector< ctype, dimensionworld > normal;
        jit.mv( refNormal, normal );
        normal *= ctype( 1 ) / jit.det();
        return normal;
      }
      
      FieldVector< ctype, dimensionworld >
      outerNormal ( const FieldVector< ctype, dimension-1 > &local ) const
      {
        const ElementGeometry &geo = insideGeometry();

        const GenericReferenceElement< ctype, dimension > &refElement
          = GenericReferenceElements< ctype, dimension>::general( geo.type() );

        FieldVector< ctype, dimension > x( geometryInInside().global( local ) );
        typedef typename ElementGeometry::Implementation::JacobianInverseTransposed JacobianInverseTransposed;
        const JacobianInverseTransposed &jit = geo.impl().jacobianInverseTransposed( x );
        const FieldVector< ctype, dimension > &refNormal = refElement.volumeOuterNormal( indexInInside() );

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
        const GenericReferenceElement< ctype, dimension-1 > &refFace
          = GenericReferenceElements< ctype, dimension-1 >::general( type() );
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
        geo_.impl() = GeometryImpl();
      }

      void initialize ( const HostIntersectionType &hostIntersection )
      {
        assert( !(*this) );
        hostIntersection_ = &hostIntersection;
      }

      const CoordFunctionType &coordFunction () const
      {
        assert( coordFunction_ );
        return *coordFunction_;
      }

    private:
      const ElementGeometry &insideGeometry () const { return insideGeo_; }

      const CoordFunctionType *coordFunction_;
      ElementGeometry insideGeo_;
      const HostIntersectionType *hostIntersection_;
      mutable Geometry geo_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_INTERSECTION_HH
