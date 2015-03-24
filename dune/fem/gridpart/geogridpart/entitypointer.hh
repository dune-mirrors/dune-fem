#ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_ENTITYPOINTER_HH
#define DUNE_FEM_GRIDPART_GEOGRIDPART_ENTITYPOINTER_HH

#include <type_traits>
#include <utility>

#include <dune/grid/common/entitypointer.hh>
#include <dune/fem/function/localfunction/localfunction.hh>

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
      GeoEntityPointer () = default;

      GeoEntityPointer ( const CoordFunctionType &coordFunction, HostIteratorType hostIterator )
      : coordFunction_( &coordFunction ),
        hostIterator_( std::move( hostIterator ) )
      {}

      GeoEntityPointer ( const EntityImpl &entity )
      : coordFunction_( &entity.coordFunction() ),
        hostIterator_( entity.hostEntity() )
      {}

      template< class T >
      explicit GeoEntityPointer ( const GeoEntityPointer< T > &other )
      : coordFunction_( other.coordFunction_ ),
        hostIterator_( other.hostIterator_ )
      {}

      GeoEntityPointer ( const ThisType & ) = default;

      GeoEntityPointer ( ThisType && ) = default;

      ThisType &operator= ( const ThisType & ) = default;

      ThisType &operator= ( ThisType && ) = default;

      template< class T >
      bool equals ( const GeoEntityPointer< T > &other ) const
      {
        return hostIterator() == other.hostIterator();
      }

      Entity dereference () const
      {
        return EntityImpl( coordFunction(), *hostIterator() );
      }

      int level () const
      {
        return hostIterator().level();
      }

      const CoordFunctionType &coordFunction () const
      {
        assert( coordFunction_ );
        return *coordFunction_;
      }

      const HostIteratorType &hostIterator() const
      {
        return hostIterator_;
      }

    private:
      const CoordFunctionType *coordFunction_ = nullptr;

    protected:
      HostIteratorType hostIterator_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_ENTITYPOINTER_HH
