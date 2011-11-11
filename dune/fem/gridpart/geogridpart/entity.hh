#ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_ENTITY_HH
#define DUNE_FEM_GRIDPART_GEOGRIDPART_ENTITY_HH

#include <dune/grid/common/entity.hh>
#include <dune/grid/common/gridenums.hh>

#include <dune/fem/gridpart/geogridpart/cornerstorage.hh>

namespace Dune
{

  namespace Fem
  {

    // GeoEntity
    // ---------

    template< int codim, int dim, class GridFamily >
    class GeoEntity
    {
      typedef typename remove_const< GridFamily >::type::Traits Traits;

    public:
      static const int codimension = codim;
      static const int dimension = remove_const< GridFamily >::type::dimension;
      static const int mydimension = dimension - codimension;
      static const int dimensionworld = remove_const< GridFamily >::type::dimensionworld;

      typedef typename remove_const< GridFamily >::type::ctype ctype;

      typedef typename Traits::template Codim< codimension >::EntitySeed EntitySeed;
      typedef typename Traits::template Codim< codimension >::Geometry Geometry;

    private:
      typedef typename Traits::HostGridPartType HostGridPartType;
      typedef typename Traits::CoordFunctionType CoordFunctionType;

      typedef typename Geometry::Implementation GeometryImpl;

      typedef GeoCoordVector< mydimension, GridFamily > CoordVectorType;

    public:
      typedef typename HostGridPartType::template Codim< codimension >::EntityType HostEntityType;
      typedef typename HostGridPartType::template Codim< codimension >::EntityPointerType HostEntityPointerType;

      explicit GeoEntity ( const CoordFunctionType &coordFunction ) 
      : coordFunction_( &coordFunction ),
        hostEntity_( 0 ),
        geo_( GeometryImpl() )
      {}

      GeoEntity ( const CoordFunctionType &coordFunction, const HostEntityType &hostEntity )
      : coordFunction_( &coordFunction ),
        hostEntity_( &hostEntity ),
        geo_( GeometryImpl() )
      {}

      GeoEntity ( const GeoEntity &other )
      : coordFunction_( other.coordFunction_ ),
        hostEntity_( other.hostEntity_ ),
        geo_( other.geo_.impl() )
      {}

      const GeoEntity &operator= ( const GeoEntity &other )
      {
        coordFunction_ = other.coordFunction_;
        hostEntity_ = other.hostEntity_;
        geo_.impl() = other.geo_.impl();
        return *this;
      }

      operator bool () const { return bool( hostEntity_ ); }

      GeometryType type () const
      {
        return hostEntity().type();
      }

      int level () const
      {
        return hostEntity().level();
      }
      
      PartitionType partitionType () const
      {
        return hostEntity().partitionType();
      }

      const Geometry &geometry () const
      {
        GeometryImpl &geo = geo_.impl();
        if( !geo )
        {
          CoordVectorType coords( coordFunction(), hostEntity() );
          geo = GeometryImpl( type(), coords );
        }
        return geo_;
      }

      EntitySeed seed () const { return EntitySeed( hostEntity().seed() ); }

      const HostEntityType &hostEntity () const
      {
        assert( *this );
        return *hostEntity_;
      }

      const CoordFunctionType &coordFunction () const
      {
        assert( coordFunction_ );
        return *coordFunction_;
      }

    private:
      const CoordFunctionType *coordFunction_;
      const HostEntityType *hostEntity_;
      mutable Geometry geo_;
    };



    // GeoEntity for codimension 0
    // ---------------------------

    template< int dim, class GridFamily >
    class GeoEntity< 0, dim, GridFamily >
    {
      typedef typename remove_const< GridFamily >::type::Traits Traits;

    public:
      static const int codimension = 0;
      static const int dimension = remove_const< GridFamily >::type::dimension;
      static const int mydimension = dimension - codimension;
      static const int dimensionworld = remove_const< GridFamily >::type::dimensionworld;

      typedef typename remove_const< GridFamily >::type::ctype ctype;

      typedef typename Traits::template Codim< codimension >::EntitySeed EntitySeed;
      typedef typename Traits::template Codim< codimension >::Geometry Geometry;
      typedef typename Traits::template Codim< codimension >::LocalGeometry LocalGeometry;
      typedef typename Traits::template Codim< codimension >::EntityPointer EntityPointer;

      typedef typename Traits::HierarchicIterator HierarchicIterator;
      typedef typename Traits::LeafIntersectionIterator LeafIntersectionIterator;
      typedef typename Traits::LevelIntersectionIterator LevelIntersectionIterator;

    private:
      typedef typename Traits::HostGridPartType HostGridPartType;
      typedef typename Traits::CoordFunctionType CoordFunctionType;

      typedef Dune::GenericGeometry::Geometry< mydimension, dimensionworld, const GridFamily > GeometryImpl;

      typedef GeoCoordVector< mydimension, GridFamily > CoordVectorType;

    public:
      typedef typename HostGridPartType::template Codim< codimension >::EntityType HostEntityType;
      typedef typename HostGridPartType::template Codim< codimension >::EntityPointerType HostEntityPointerType;

      explicit GeoEntity ( const CoordFunctionType &coordFunction ) 
      : coordFunction_( &coordFunction ),
        hostEntity_( 0 ),
        geo_( GeometryImpl() )
      {}

      GeoEntity ( const CoordFunctionType &coordFunction, const HostEntityType &hostEntity )
      : coordFunction_( &coordFunction ),
        hostEntity_( &hostEntity ),
        geo_( GeometryImpl() )
      {}

      GeoEntity ( const GeoEntity &other )
      : coordFunction_( other.coordFunction_ ),
        hostEntity_( other.hostEntity_ ),
        geo_( other.geo_.impl() )
      {}

      const GeoEntity &operator= ( const GeoEntity &other )
      {
        coordFunction_ = other.coordFunction_;
        hostEntity_ = other.hostEntity_;
        geo_.impl() = other.geo_.impl();
        return *this;
      }

      operator bool () const { return bool( hostEntity_ ); }

      GeometryType type () const
      {
        return hostEntity().type();
      }

      int level () const
      {
        DUNE_THROW( InvalidStateException, "Trying to access level intersection iterator from a grid part." );
      }
      
      PartitionType partitionType () const
      {
        return hostEntity().partitionType();
      }

      const Geometry &geometry () const
      {
        GeometryImpl &geo = geo_.impl();
        if( !geo )
        {
          CoordVectorType coords( coordFunction(), hostEntity() );
          geo = GeometryImpl( type(), coords );
        }
        return geo_;
      }

      EntitySeed seed () const { return EntitySeed( hostEntity().seed() ); }

      template< int codim >
      int count () const
      {
        return hostEntity().template count< codim >();
      }
      
      template< int codim >
      typename Traits::template Codim< codim >::EntityPointer
      subEntity ( int i ) const
      {
        typedef typename Traits::template Codim< codim >::EntityPointerImpl EntityPointerImpl;
        return EntityPointerImpl( hostEntity().template subEntity< codim >( i ) );
      }

      LevelIntersectionIterator ilevelbegin () const
      {
        DUNE_THROW( InvalidStateException, "Trying to access level intersection iterator from a grid part." );
      }
      
      LevelIntersectionIterator ilevelend () const
      {
        DUNE_THROW( InvalidStateException, "Trying to access level intersection iterator from a grid part." );
      }
      
      LeafIntersectionIterator ileafbegin () const
      {
        DUNE_THROW( InvalidStateException, "Trying to access leaf intersection iterator from a grid part." );
      }
      
      LeafIntersectionIterator ileafend () const
      {
        DUNE_THROW( InvalidStateException, "Trying to access leaf intersection iterator from a grid part." );
      }

      bool hasBoundaryIntersections () const
      {
        return hostEntity().hasBoundaryIntersections();
      }

      bool isLeaf () const
      {
        DUNE_THROW( InvalidStateException, "Trying to access hierarchy information from a grid part." );
      }
   
      EntityPointer father () const
      {
        DUNE_THROW( InvalidStateException, "Trying to access hierarchy information from a grid part." );
      }

      bool hasFather () const
      {
        DUNE_THROW( InvalidStateException, "Trying to access hierarchy information from a grid part." );
      }
        
      const LocalGeometry &geometryInFather () const
      {
        DUNE_THROW( InvalidStateException, "Trying to access hierarchy information from a grid part." );
      }
   
      HierarchicIterator hbegin ( int maxLevel ) const
      {
        DUNE_THROW( InvalidStateException, "Trying to access hierarchy information from a grid part." );
      }
      
      HierarchicIterator hend ( int maxLevel ) const
      {
        DUNE_THROW( InvalidStateException, "Trying to access hierarchy information from a grid part." );
      }

      bool isRegular () const
      {
        DUNE_THROW( InvalidStateException, "Trying to access hierarchy information from a grid part." );
      }

      bool isNew () const
      {
        DUNE_THROW( InvalidStateException, "Trying to access hierarchy information from a grid part." );
      }
            
      bool mightVanish () const
      {
        DUNE_THROW( InvalidStateException, "Trying to access hierarchy information from a grid part." );
      }

      const HostEntityType &hostEntity () const
      {
        assert( *this );
        return *hostEntity_;
      }

      const CoordFunctionType &coordFunction () const
      {
        assert( coordFunction_ );
        return *coordFunction_;
      }

    private:
      const CoordFunctionType *coordFunction_;
      const HostEntityType *hostEntity_;
      mutable Geometry geo_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_ENTITY_HH
