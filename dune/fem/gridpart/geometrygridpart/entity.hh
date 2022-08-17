#ifndef DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_ENTITY_HH
#define DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_ENTITY_HH

#include <type_traits>
#include <utility>

#include <dune/grid/common/entity.hh>
#include <dune/grid/common/gridenums.hh>

#include <dune/fem/gridpart/common/defaultgridpartentity.hh>

namespace Dune
{

  namespace Fem
  {

    // GeometryGridPartEntity
    // ----------------------

    template< int codim, int dim, class GridFamily >
    class GeometryGridPartEntity
      : public DefaultGridPartEntity < codim, dim, GridFamily >
    {
      typedef typename std::remove_const< GridFamily >::type::Traits Traits;
      typedef typename Traits::GridFunctionType GridFunctionType;

      //typedef typename GridFamily::GridType GridType;
      //typedef typename GridType::template Codim<codim>::Entity GridEntityType;

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
      typedef typename Geometry::Implementation GeometryImplType;

    public:
      typedef typename HostGridPartType::template Codim< codimension >::EntityType HostEntityType;

      GeometryGridPartEntity () = default;

      GeometryGridPartEntity ( const GridFunctionType &gridFunction, HostEntityType hostEntity )
        : hostEntity_( std::move( hostEntity ) ), gridFunction_( &gridFunction )
      {}

      GeometryType type () const
      {
        return hostEntity().type();
      }

      PartitionType partitionType () const
      {
        return hostEntity().partitionType();
      }

      unsigned int subEntities ( unsigned int c ) const { return hostEntity().subEntities( c ); }

      Geometry geometry () const
      {
        DUNE_THROW( NotImplemented, "GeometryGridPart only implements the geometry for entities of codimension 0." );
      }

      EntitySeed seed () const { return EntitySeed( hostEntity().seed() ); }

      bool equals ( const GeometryGridPartEntity &rhs ) const
      {
        return hostEntity() == rhs.hostEntity();
      }

      const HostEntityType &hostEntity () const
      {
        return hostEntity_;
      }

      const GridFunctionType &gridFunction () const
      {
        assert( gridFunction_ );
        return *gridFunction_;
      }

    private:
      HostEntityType hostEntity_;
      const GridFunctionType *gridFunction_ = nullptr;
    };



    // GeometryGridPartEntity for codimension 0
    // ----------------------------------------

    template< int dim, class GridFamily >
    class GeometryGridPartEntity< 0, dim, GridFamily >
      : public DefaultGridPartEntity < 0, dim, GridFamily >
    {
      typedef typename std::remove_const< GridFamily >::type::Traits Traits;
      typedef typename Traits::GridFunctionType GridFunctionType;

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
      typedef typename HostGridPartType::GridType GridType;

    public:
      typedef typename HostGridPartType::template Codim< codimension >::EntityType HostEntityType;

      GeometryGridPartEntity () = default;

      GeometryGridPartEntity ( const GridFunctionType &gridFunction, HostEntityType hostEntity )
        : hostEntity_( std::move( hostEntity ) ), gridFunction_( &gridFunction )
      {}

      template< class LocalFunction >
      GeometryGridPartEntity ( const GeometryGridPartEntity &other, const LocalFunction &localGridFunction )
        : hostEntity_( other.hostEntity_ ), gridFunction_( other.gridFunction_ )
      {
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
        typedef typename Geometry::Implementation Impl;
        Impl impl( gridFunction() );
        impl.impl().bind( hostEntity() );
        return Geometry( impl );
      }

      EntitySeed seed () const { return EntitySeed( hostEntity().seed() ); }

      template< int codim >
      int count () const
      {
        return hostEntity().template count< codim >();
      }

      unsigned int subEntities ( unsigned int codim ) const { return hostEntity().subEntities( codim ); }

      template< int codim >
      typename Traits::template Codim< codim >::Entity subEntity ( int i ) const
      {
        typedef typename Traits::template Codim< codim >::Entity::Implementation EntityImpl;
        return EntityImpl( *gridFunction_, hostEntity().template subEntity< codim >( i ) );
      }

      bool hasBoundaryIntersections () const
      {
        return hostEntity().hasBoundaryIntersections();
      }

      bool equals ( const GeometryGridPartEntity &rhs ) const
      {
        return hostEntity() == rhs.hostEntity();
      }

      const HostEntityType &hostEntity () const
      {
        return hostEntity_;
      }

      void setHostEntity ( const HostEntityType &hostEntity )
      {
        hostEntity_ = &hostEntity;
      }

      const GridFunctionType &gridFunction () const
      {
        assert( gridFunction_ );
        return *gridFunction_;
      }
    private:
      HostEntityType hostEntity_;
      const GridFunctionType *gridFunction_ = nullptr;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_ENTITY_HH
