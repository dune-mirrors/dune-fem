#ifndef DUNE_FEM_GRIDPART_IDGRIDPART_INTERSECTION_HH
#define DUNE_FEM_GRIDPART_IDGRIDPART_INTERSECTION_HH

#include <dune/fem/gridpart/idgridpart/entitypointer.hh>
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
      typedef typename remove_const< GridFamily >::type::Traits Traits;

      typedef typename Traits::HostGridPartType HostGridPartType;

    public:
      typedef typename remove_const< GridFamily >::type::ctype ctype;
      
      static const int dimension = remove_const< GridFamily >::type::dimension;
      static const int dimensionworld = remove_const< GridFamily >::type::dimensionworld;

      typedef typename Traits::template Codim< 0 >::Entity Entity;
      typedef typename Traits::template Codim< 0 >::EntityPointer EntityPointer;
      typedef typename Traits::template Codim< 1 >::Geometry Geometry;
      typedef typename Traits::template Codim< 1 >::LocalGeometry LocalGeometry;

    private:
      typedef typename EntityPointer::Implementation EntityPointerImpl;

      typedef typename HostGridPartType::IntersectionType HostIntersectionType;

    public:
      IdIntersection ()
      : hostIntersection_( 0 )
      {}

      explicit IdIntersection ( const HostIntersectionType &hostIntersection )
      : hostIntersection_( &hostIntersection )
      {}

      operator bool () const { return bool( hostIntersection_ ); }

      EntityPointer inside () const
      {
        return EntityPointerImpl( hostIntersection().inside() );
      }
      
      EntityPointer outside () const
      {
        return EntityPointerImpl( hostIntersection().outside() );
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

      const HostIntersectionType &hostIntersection () const
      {
        assert( *this );
        return *hostIntersection_;
      }

    private:
      const HostIntersectionType *hostIntersection_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_IDGRIDPART_INTERSECTION_HH
