#ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_ENTITYPOINTER_HH
#define DUNE_FEM_GRIDPART_GEOGRIDPART_ENTITYPOINTER_HH

#include <dune/grid/common/entitypointer.hh>

namespace Dune
{

  namespace Fem
  {

    // External Forward Declarations
    // -----------------------------

    template< int, int, class >
    class GeoEntity;



    // GeoEntityPointerTraits
    // ----------------------
    
    template< int codim, class GridFamily >
    struct GeoEntityPointerTraits
    {
      typedef GeoEntityPointerTraits< codim, const GridFamily > BaseTraits;

      typedef typename remove_const< GridFamily >::type::ctype ctype;

      static const int dimension = remove_const< GridFamily >::type::dimension;
      static const int codimension = codim;

      typedef GeoEntity< codimension, dimension, const GridFamily > EntityImpl;
      typedef Dune::Entity< codimension, dimension, const GridFamily, GeoEntity > Entity;

      typedef typename remove_const< GridFamily >::type::Traits::CoordFunctionType CoordFunctionType;
      typedef typename remove_const< GridFamily >::type::Traits::HostGridPartType HostGridPartType;

      typedef typename HostGridPartType::template Codim< codim >::EntityType HostEntityType;
      typedef typename HostGridPartType::template Codim< codim >::EntityPointerType HostEntityPointerType;
      typedef HostEntityPointerType HostIteratorType;
    };



    // GeoEntityPointer
    // ----------------

    template< class Traits >
    class GeoEntityPointer
    {
      typedef GeoEntityPointer< Traits > ThisType;

      friend class GeoEntityPointer< typename Traits::BaseTraits >;

    public:
      static const int dimension = Traits::dimension;
      static const int codimension = Traits::codimension;
   
      typedef typename Traits::Entity Entity;

      typedef GeoEntityPointer< typename Traits::BaseTraits > EntityPointerImp;

      typedef typename Traits::CoordFunctionType CoordFunctionType;

    protected:
      typedef typename Traits::HostEntityPointerType HostEntityPointerType;
      typedef typename Traits::HostIteratorType HostIteratorType;

    private:
      typedef typename Traits::EntityImpl EntityImpl;

    public:
      GeoEntityPointer ( const CoordFunctionType &coordFunction, const HostIteratorType &hostIterator )
      : entity_( EntityImpl( coordFunction ) ),
        hostIterator_( hostIterator )
      {}

      GeoEntityPointer ( const EntityImpl &entity )
      : entity_( EntityImpl( entity.coordFunction() ) ),
        hostIterator_( entity.hostEntity() )
      {}

      GeoEntityPointer ( const ThisType &other )
      : entity_( EntityImpl( other.entity_.impl().coordFunction() ) ),
        hostIterator_( other.hostIterator_ )
      {}

      template< class T >
      explicit GeoEntityPointer ( const GeoEntityPointer< T > &other )
      : entity_( EntityImpl( other.entity_.impl().coordFunction() ) ),
        hostIterator_( other.hostIterator_ )
      {}
      
      ThisType &operator= ( const ThisType &other )
      {
        entity_.impl() = EntityImpl( other.entity_.impl().coordFunction() );
        hostIterator_ = other.hostIterator_;
        return *this;
      }

      operator const EntityPointerImp & () const
      {
        return reinterpret_cast< const EntityPointerImp & >( *this );
      }

      template< class T >
      bool equals ( const GeoEntityPointer< T > &other ) const
      {
        return (hostIterator() == other.hostIterator());
      }
      
      Entity &dereference () const
      {
        if( !entity_.impl() )
          entity_.impl() = EntityImpl( entity_.impl().coordFunction(), *hostIterator() );
        return entity_;
      }
      
      int level () const
      {
        return hostIterator().level();
      }

      const HostIteratorType &hostIterator() const
      {
        return hostIterator_;
      }

    protected:
      void releaseEntity ()
      {
        entity_.impl() = EntityImpl( entity_.impl().coordFunction() );
      }

    private:
      mutable Entity entity_;

    protected:
      HostIteratorType hostIterator_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_ENTITYPOINTER_HH
