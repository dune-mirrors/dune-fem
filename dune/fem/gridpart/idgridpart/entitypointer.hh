#ifndef DUNE_FEM_GRIDPART_IDGRIDPART_ENTITYPOINTER_HH
#define DUNE_FEM_GRIDPART_IDGRIDPART_ENTITYPOINTER_HH

#include <dune/grid/common/entitypointer.hh>

namespace Dune
{

  namespace Fem
  {

    // External Forward Declarations
    // -----------------------------

    template< int, int, class >
    class IdEntity;



    // IdEntityPointerTraits
    // ---------------------

    template< int codim, class GridFamily >
    struct IdEntityPointerTraits
    {
      typedef IdEntityPointerTraits< codim, const GridFamily > BaseTraits;

      typedef typename remove_const< GridFamily >::type::ctype ctype;

      static const int dimension = remove_const< GridFamily >::type::dimension;
      static const int codimension = codim;

      typedef IdEntity< codimension, dimension, const GridFamily > EntityImpl;
      typedef Dune::Entity< codimension, dimension, const GridFamily, IdEntity > Entity;

      typedef typename remove_const< GridFamily >::type::Traits::HostGridPartType HostGridPartType;

      typedef typename HostGridPartType::template Codim< codim >::EntityType HostEntityType;
      typedef typename HostGridPartType::template Codim< codim >::EntityPointerType HostEntityPointerType;
      typedef HostEntityPointerType HostIteratorType;
    };



    // IdEntityPointer
    // ---------------

    template< class Traits >
    class IdEntityPointer
    {
      typedef IdEntityPointer< Traits > ThisType;

      friend class IdEntityPointer< typename Traits::BaseTraits >;

    public:
      static const int dimension = Traits::dimension;
      static const int codimension = Traits::codimension;

      typedef typename Traits::Entity Entity;

      typedef IdEntityPointer< typename Traits::BaseTraits > EntityPointerImp;

    protected:
      typedef typename Traits::HostEntityPointerType HostEntityPointerType;
      typedef typename Traits::HostIteratorType HostIteratorType;

    private:
      typedef typename Traits::EntityImpl EntityImpl;

    public:
      IdEntityPointer ( const HostIteratorType &hostIterator )
      : entity_( EntityImpl() ),
        hostIterator_( hostIterator )
      {}

      IdEntityPointer ( const EntityImpl &entity )
      : entity_( EntityImpl() ),
        hostIterator_( entity.hostEntity() )
      {}

      IdEntityPointer ( const ThisType &other )
      : entity_( EntityImpl() ),
        hostIterator_( other.hostIterator_ )
      {}

      template< class T >
      explicit IdEntityPointer ( const IdEntityPointer< T > &other )
      : entity_( EntityImpl() ),
        hostIterator_( other.hostIterator_ )
      {}

      ThisType &operator= ( const ThisType &other )
      {
        releaseEntity();
        hostIterator_ = other.hostIterator_;
        return *this;
      }

      operator const EntityPointerImp & () const
      {
        return reinterpret_cast< const EntityPointerImp & >( *this );
      }

      template< class T >
      bool equals ( const IdEntityPointer< T > &other ) const
      {
        return (hostIterator() == other.hostIterator());
      }

      Entity &dereference () const
      {
        if( !entity_.impl() )
          entity_.impl() = EntityImpl( *hostIterator() );
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
        entity_.impl() = EntityImpl();
      }

    private:
      mutable Entity entity_;

    protected:
      HostIteratorType hostIterator_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_IDGRIDPART_ENTITYPOINTER_HH
