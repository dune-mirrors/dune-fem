#ifndef DUNE_FEM_GRIDPART_IDGRIDPART_INTERSECTION_HH
#define DUNE_FEM_GRIDPART_IDGRIDPART_INTERSECTION_HH

#include <type_traits>
#include <utility>

#include <dune/fem/gridpart/idgridpart/geometry.hh>

namespace Dune
{

  namespace Fem
  {

    // IdIntersection
    // --------------

    template< class GridFamily >
    class IdIntersection
    {
      typedef typename std::remove_const< GridFamily >::type::Traits Traits;

      typedef typename Traits::HostGridPartType HostGridPartType;

    public:
      typedef typename std::remove_const< GridFamily >::type::ctype ctype;

      static const int dimension = std::remove_const< GridFamily >::type::dimension;
      static const int dimensionworld = std::remove_const< GridFamily >::type::dimensionworld;

      typedef typename Traits::template Codim< 0 >::Entity Entity;
      typedef typename Traits::template Codim< 1 >::Geometry Geometry;
      typedef typename Traits::template Codim< 1 >::LocalGeometry LocalGeometry;

      typedef typename Traits::ExtraData  ExtraData;

    private:
      typedef typename Entity::Implementation EntityImpl;

      typedef typename HostGridPartType::IntersectionType HostIntersectionType;

    public:
      IdIntersection () = default;

      IdIntersection ( ExtraData data, HostIntersectionType hostIntersection )
      : data_( std::move( data ) ),
        hostIntersection_( std::move( hostIntersection ) )
      {}

      Entity inside () const
      {
        return Entity( EntityImpl( data(), hostIntersection().inside() ) );
      }

      Entity outside () const
      {
        return Entity( EntityImpl( data(), hostIntersection().outside() ) );
      }

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
        return Geometry( hostIntersection().geometry() );
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
        return hostIntersection().integrationOuterNormal( local );
      }

      FieldVector< ctype, dimensionworld >
      outerNormal ( const FieldVector< ctype, dimension-1 > &local ) const
      {
        return hostIntersection().outerNormal( local );
      }

      FieldVector< ctype, dimensionworld >
      unitOuterNormal ( const FieldVector< ctype, dimension-1 > &local ) const
      {
        return hostIntersection().unitOuterNormal( local );
      }

      FieldVector< ctype, dimensionworld > centerUnitOuterNormal () const
      {
        return hostIntersection().centerUnitOuterNormal();
      }

      const ExtraData &data () const { return data_; }

      const HostIntersectionType &hostIntersection () const
      {
        return hostIntersection_;
      }

    protected:
      ExtraData data_;
      HostIntersectionType hostIntersection_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_IDGRIDPART_INTERSECTION_HH
