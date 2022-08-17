#ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_ENTITY_HH
#define DUNE_FEM_GRIDPART_GEOGRIDPART_ENTITY_HH

#include <type_traits>
#include <utility>

#include <dune/grid/common/entity.hh>
#include <dune/grid/common/gridenums.hh>

#include <dune/fem/gridpart/common/defaultgridpartentity.hh>
#include <dune/fem/gridpart/geogridpart/cornerstorage.hh>

namespace Dune
{

  namespace Fem
  {

    // GeoEntity
    // ---------

    template< int codim, int dim, class GridFamily >
    class GeoEntity
    : public DefaultGridPartEntity < codim, dim, GridFamily >
    {
      typedef typename std::remove_const< GridFamily >::type::Traits Traits;

    public:
      static const int codimension = codim;
      static const int dimension = std::remove_const< GridFamily >::type::dimension;
      static const int mydimension = dimension - codimension;
      static const int dimensionworld = std::remove_const< GridFamily >::type::dimensionworld;

      typedef typename std::remove_const< GridFamily >::type::ctype ctype;

      typedef typename Traits::template Codim< codimension >::EntitySeed EntitySeed;
      typedef typename Traits::template Codim< codimension >::Geometry Geometry;

      typedef typename Traits::HostGridPartType HostGridPartType;

    private:
      typedef typename Traits::CoordFunctionType CoordFunctionType;

      typedef typename Geometry::Implementation GeometryImplType;

      typedef GeoCoordVector< mydimension, GridFamily > CoordVectorType;

    public:
      typedef typename HostGridPartType::template Codim< codimension >::EntityType HostEntityType;

      GeoEntity () = default;

      GeoEntity ( const CoordFunctionType &coordFunction, HostEntityType hostEntity )
      : coordFunction_( &coordFunction ),
        hostEntity_( std::move( hostEntity ) )
      {}

      GeometryType type () const
      {
        return hostEntity().type();
      }

      PartitionType partitionType () const
      {
        return hostEntity().partitionType();
      }

      Geometry geometry () const
      {
        if( !geo_ )
        {
          CoordVectorType coords( coordFunction(), hostEntity() );
          geo_ = GeometryImplType( type(), coords );
        }
        return Geometry( geo_ );
      }

      EntitySeed seed () const { return EntitySeed( hostEntity().seed() ); }

      bool equals ( const GeoEntity &rhs ) const
      {
        return hostEntity() == rhs.hostEntity();
      }

      const CoordFunctionType &coordFunction () const
      {
        assert( coordFunction_ );
        return *coordFunction_;
      }

      const HostEntityType &hostEntity () const
      {
        return hostEntity_;
      }

      unsigned int subEntities ( unsigned int cdim ) const { return hostEntity().subEntities( cdim ); }

    private:
      const CoordFunctionType *coordFunction_ = nullptr;
      HostEntityType hostEntity_;

      mutable GeometryImplType geo_;
    };



    // GeoEntity for codimension 0
    // ---------------------------

    template< int dim, class GridFamily >
    class GeoEntity< 0, dim, GridFamily >
    : public DefaultGridPartEntity < 0, dim, GridFamily >
    {
      typedef typename std::remove_const< GridFamily >::type::Traits Traits;

    public:
      static const int codimension = 0;
      static const int dimension = std::remove_const< GridFamily >::type::dimension;
      static const int mydimension = dimension - codimension;
      static const int dimensionworld = std::remove_const< GridFamily >::type::dimensionworld;

      typedef typename std::remove_const< GridFamily >::type::ctype ctype;

      typedef typename Traits::template Codim< codimension >::EntitySeed EntitySeed;
      typedef typename Traits::template Codim< codimension >::Geometry Geometry;
      typedef typename Traits::template Codim< codimension >::LocalGeometry LocalGeometry;

      typedef typename Traits::HierarchicIterator HierarchicIterator;
      typedef typename Traits::LeafIntersectionIterator LeafIntersectionIterator;
      typedef typename Traits::LevelIntersectionIterator LevelIntersectionIterator;

      typedef typename Traits::HostGridPartType HostGridPartType;
    private:

      typedef typename Traits::CoordFunctionType CoordFunctionType;

      typedef typename Geometry::Implementation GeometryImplType;

      typedef GeoCoordVector< mydimension, GridFamily > CoordVectorType;

    public:
      typedef typename HostGridPartType::template Codim< codimension >::EntityType HostEntityType;

      GeoEntity () = default;

      GeoEntity ( const CoordFunctionType &coordFunction, HostEntityType hostEntity )
      : coordFunction_( &coordFunction ),
        hostEntity_( std::move( hostEntity ) )
      {}

      template< class LocalFunction >
      GeoEntity ( const GeoEntity &other, const LocalFunction &localCoordFunction )
      : coordFunction_( other.coordFunction_ ),
        hostEntity_( other.hostEntity_ )
      {
        GeoLocalCoordVector< mydimension, GridFamily, LocalFunction > coords( localCoordFunction );
        geo_ = GeometryImplType( type(), coords );
      }

      GeometryType type () const
      {
        return hostEntity().type();
      }

      PartitionType partitionType () const
      {
        return hostEntity().partitionType();
      }

      Geometry geometry () const
      {
        if( !geo_ )
        {
          CoordVectorType coords( coordFunction(), hostEntity() );
          geo_ = GeometryImplType( type(), coords );
        }
        return Geometry( geo_ );
      }

      EntitySeed seed () const { return EntitySeed( hostEntity().seed() ); }

      template< int codim >
      int count () const
      {
        return hostEntity().template count< codim >();
      }

      unsigned int subEntities ( unsigned int codim ) const { return hostEntity().subEntities( codim ); }

      template< int codim >
      typename Traits::template Codim< codim >::Entity
      subEntity ( int i ) const
      {
        typedef typename Traits::template Codim< codim >::Entity::Implementation EntityImpl;
        return EntityImpl( *coordFunction_, hostEntity().template subEntity< codim >( i ) );
      }

      bool hasBoundaryIntersections () const
      {
        return hostEntity().hasBoundaryIntersections();
      }

      bool equals ( const GeoEntity &rhs ) const
      {
        return hostEntity() == rhs.hostEntity();
      }

      const CoordFunctionType &coordFunction () const
      {
        assert( coordFunction_ );
        return *coordFunction_;
      }

      const HostEntityType &hostEntity () const
      {
        return hostEntity_;
      }

      void setHostEntity ( const HostEntityType &hostEntity )
      {
        hostEntity_ = &hostEntity;
      }

    private:
      const CoordFunctionType *coordFunction_ = nullptr;
      HostEntityType hostEntity_;

      mutable GeometryImplType geo_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_ENTITY_HH
