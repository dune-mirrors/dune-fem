#ifndef DUNE_FEM_GRIDPART_IDGRIDPART_ENTITY_HH
#define DUNE_FEM_GRIDPART_IDGRIDPART_ENTITY_HH

#include <dune/grid/common/entity.hh>
#include <dune/grid/common/gridenums.hh>

#include <dune/fem/gridpart/common/defaultgridpartentity.hh>
#include <dune/fem/gridpart/idgridpart/geometry.hh>


namespace Dune
{

  namespace Fem
  {

    // IdEntity
    // --------

    template< int codim, int dim, class GridFamily >
    class IdEntity 
    : public DefaultGridPartEntity < codim, dim, GridFamily > 
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

    public:
      typedef typename HostGridPartType::template Codim< codimension >::EntityType HostEntityType;
      typedef typename HostGridPartType::template Codim< codimension >::EntityPointerType HostEntityPointerType;

      IdEntity () 
      : hostEntity_( 0 )
      {}

      explicit IdEntity ( const HostEntityType &hostEntity )
      : hostEntity_( &hostEntity )
      {}

      operator bool () const { return bool( hostEntity_ ); }

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
        return Geometry( hostEntity().geometry() );
      }

      EntitySeed seed () const { return EntitySeed( hostEntity().seed() ); }

      const HostEntityType &hostEntity () const
      {
        assert( *this );
        return *hostEntity_;
      }

    private:
      const HostEntityType *hostEntity_;
    };



    // IdEntity for codimension 0
    // --------------------------

    template< int dim, class GridFamily >
    class IdEntity< 0, dim, GridFamily > 
    : public DefaultGridPartEntity < 0, dim, GridFamily > 
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

    public:
      typedef typename HostGridPartType::template Codim< codimension >::EntityType HostEntityType;
      typedef typename HostGridPartType::template Codim< codimension >::EntityPointerType HostEntityPointerType;

    public:
      IdEntity ()
      : hostEntity_( 0 )
      {}

      explicit IdEntity ( const HostEntityType &hostEntity )
      : hostEntity_( &hostEntity )
      {}

      operator bool () const { return bool( hostEntity_ ); }

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
        return Geometry( hostEntity().geometry() );
      }

      EntitySeed seed () const { return EntitySeed( hostEntity().seed() ); }

      template< int codim >
      int count () const
      {
        return hostEntity().template count< codim >();
      }
      
      unsigned int subEntities ( unsigned int codim ) const { return hostEntity().subEntities(codim); }

      template< int codim >
      typename Traits::template Codim< codim >::EntityPointer
      subEntity ( int i ) const
      {
        typedef typename Traits::template Codim< codim >::EntityPointerImpl EntityPointerImpl;
        return EntityPointerImpl( hostEntity().template subEntity< codim >( i ) );
      }

      bool hasBoundaryIntersections () const
      {
        return hostEntity().hasBoundaryIntersections();
      }

      const HostEntityType &hostEntity () const
      {
        assert( *this );
        return *hostEntity_;
      }

    private:
      const HostEntityType *hostEntity_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_IDGRIDPART_ENTITY_HH
