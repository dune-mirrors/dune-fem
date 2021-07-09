#ifndef DUNE_FEM_SPACE_LAGRANGE_STORAGE_HH
#define DUNE_FEM_SPACE_LAGRANGE_STORAGE_HH

// local includes
#include "dofmappercode.hh"


namespace Dune
{

  namespace Fem
  {

    // LagrangeMapperSingletonKey
    // --------------------------

    template< class GridPart, class LagrangePointSetContainer >
    class LagrangeMapperSingletonKey
    {
      typedef LagrangeMapperSingletonKey< GridPart, LagrangePointSetContainer > ThisType;

    public:
      typedef GridPart GridPartType;
      typedef LagrangePointSetContainer LagrangePointSetContainerType;

      LagrangeMapperSingletonKey ( const GridPartType &gridPart,
                                   const LagrangePointSetContainerType &pointSet,
                                   const int polOrd )
      : gridPart_( gridPart ),
        pointSet_( pointSet ),
        polOrd_( polOrd )
      {
      }

      bool operator== ( const ThisType &other ) const
      {
        return ((&indexSet() == &other.indexSet()) && (polOrd_ == other.polOrd_));
      }

      bool operator!= ( const ThisType &other ) const
      {
        return !( *this == other );
      }

      const GridPartType &gridPart () const
      {
        return gridPart_;
      }

      const typename GridPartType::IndexSetType &indexSet () const
      {
        return gridPart().indexSet();
      }

      const LagrangePointSetContainerType &pointSet () const {
        return pointSet_;
      }

    private:
      const GridPartType &gridPart_;
      const LagrangePointSetContainerType &pointSet_;
      const int polOrd_;
    };



    // LagrangeMapperSingletonFactory
    // ------------------------------

    template< class Key, class Object >
    struct LagrangeMapperSingletonFactory
    {
      typedef typename Key::LagrangePointSetContainerType LagrangePointSetContainerType;

      static Object *createObject ( const Key &key )
      {
        LagrangeDofMapperCodeFactory< LagrangePointSetContainerType > codeFactory( key.pointSet() );
        return new Object( key.gridPart(), codeFactory );
      }

      static void deleteObject ( Object *obj )
      {
        delete obj;
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LAGRANGE_STORAGE_HH
