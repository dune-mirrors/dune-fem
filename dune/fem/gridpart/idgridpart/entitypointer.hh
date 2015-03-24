#ifndef DUNE_FEM_GRIDPART_IDGRIDPART_ENTITYPOINTER_HH
#define DUNE_FEM_GRIDPART_IDGRIDPART_ENTITYPOINTER_HH

#include <type_traits>
#include <utility>

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

      typedef typename GridFamily::Traits::ExtraData  ExtraData;
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

      typedef typename Traits::ExtraData  ExtraData;

    protected:
      typedef typename Traits::HostEntityPointerType HostEntityPointerType;
      typedef typename Traits::HostIteratorType HostIteratorType;

    private:
      typedef typename Traits::EntityImpl EntityImpl;

    public:
      IdEntityPointer ( ExtraData data, HostIteratorType hostIterator )
      : data_( std::move( data ) ),
        hostIterator_( std::move( hostIterator ) )
      {}

      IdEntityPointer ( const EntityImpl &entity )
      : data_( entity.data() ),
        hostIterator_( entity.hostEntity() )
      {}

      template< class T >
      explicit IdEntityPointer ( const IdEntityPointer< T > &other )
      : data_( other.data() ),
        hostIterator_( other.hostIterator_ )
      {}

      IdEntityPointer ( const ThisType & ) = default;

      IdEntityPointer ( ThisType && ) = default;

      ThisType &operator= ( const ThisType & ) = default;

      ThisType &operator= ( ThisType && ) = default;

      template< class T >
      bool equals ( const IdEntityPointer< T > &other ) const
      {
        return hostIterator() == other.hostIterator();
      }

      Entity dereference () const
      {
        return EntityImpl( data(), *hostIterator() );
      }

      int level () const
      {
        return hostIterator().level();
      }

      const ExtraData &data () const { return data_; }

      const HostIteratorType &hostIterator() const
      {
        return hostIterator_;
      }

    private:
      ExtraData data_;

    protected:
      HostIteratorType hostIterator_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_IDGRIDPART_ENTITYPOINTER_HH
