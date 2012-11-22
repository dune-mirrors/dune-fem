#ifndef DUNE_FEM_DEFAULTINDEXSETS_HH
#define DUNE_FEM_DEFAULTINDEXSETS_HH

//- system includes 
#include <vector>
#include <rpc/rpc.h>

//- Dune includes 
#include <dune/common/misc.hh>
#include <dune/common/typetraits.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/common/adaptcallback.hh> // for compatibility only
#include <dune/fem/gridpart/emptyindexset.hh>

#include <dune/fem/misc/capabilities.hh>
#include <dune/fem/misc/mpimanager.hh>

/** \file
 *  \author Robert Kloefkorn
 */

namespace Dune
{

  namespace Fem 
  {

    //! Wraps the interface methods of indexsets and adds the addiotnal needed
    //! functions 
    template< class IndexSetImp >
    class IndexSetWrapper
    : public Fem :: EmptyIndexSet 
    {
      typedef IndexSetWrapper< IndexSetImp > ThisType;
      typedef Fem :: EmptyIndexSet  BaseType ;

    public:
      typedef typename IndexSetImp::IndexType IndexType;

      //! store const reference to set 
      IndexSetWrapper ( const IndexSetImp &set )
      : BaseType(),
        set_( set )
      {}
      
      //! copy constructor
      IndexSetWrapper ( const ThisType &other )
      : BaseType(),
        set_( other.set_ )
      {}

      //! return persistent status 
      bool persistent () const { return false; }
     
      //! return size of set for codim  
      IndexType size ( GeometryType type ) const   
      {
        return set_.size(type);
      }

      //! return size of grid entities per level and codim 
      IndexType size ( int codim ) const   
      {
        return set_.size(codim);
      }

      //! return index of en 
      template< class Entity >
      IndexType index ( const Entity &entity ) const
      {
        return set_.index( entity );
      }

      template< class Entity >
      IndexType subIndex ( const Entity &entity, int num, unsigned int codim ) const
      {
        return set_.subIndex( entity, num, codim );
      }

      //! wrap geomTypes method of set 
      const std::vector< GeometryType > &geomTypes ( const int codim ) const 
      {
        return set_.geomTypes(codim); 
      }

      //! returns true if this set provides an index for given entity
      template< class Entity >
      bool contains ( const Entity &entity ) const
      {
        return set_.contains( entity );
      }

    private: 
      const IndexSetImp &set_;
    };



    //! Wraps LevelIndexSet for use with LagrangeFunctionSpace 
    template< class GridType >
    class WrappedLevelIndexSet
    : public IndexSetWrapper< typename GridType::Traits::LevelIndexSet >
    {
      typedef WrappedLevelIndexSet< GridType > ThisType;
      typedef IndexSetWrapper< typename GridType::Traits::LevelIndexSet > BaseType;

      // my type, to be revised 
      enum { myType = 1 };

    public:
      //! number of codimensions 
      enum { ncodim = GridType::dimension + 1 };

      //! Constructor getting grid and level for Index Set 
      WrappedLevelIndexSet ( const GridType &grid, const int level )
      : BaseType( grid.levelIndexSet( level ) )
      {}
     
      //! return type of index set (for input/output)
      static int type() { return myType; }
    };



    //! compile time chooser for hierarchic or leaf index set
    /** \deprecated */
    template< class Grid >
    class HierarchicIndexSetSelector
    {
      template< bool >
      struct HierarchicIndexSetGetter
      {
        typedef typename Grid::HierarchicIndexSet IndexSet;

        static const IndexSet &indexSet ( const Grid &grid )
        {
          return grid.hierarchicIndexSet();
        }
      };

      template< bool >
      struct LeafIndexSetGetter
      {
        typedef typename Grid::LeafIndexSet IndexSet;

        static const IndexSet &indexSet ( const Grid &grid )
        {
          if( Dune::Fem::MPIManager::rank() == 0 )
            std::cerr << "Warning: Grid does not provide a HierarchicIndexSet, using LeafIndexSet instead." << std::endl;
          return grid.leafIndexSet();
        }
      };

      static const bool hasHierarchicIndexSet = Dune::Capabilities::hasHierarchicIndexSet< Grid >::v;
      typedef typename SelectType< hasHierarchicIndexSet, HierarchicIndexSetGetter< true >, LeafIndexSetGetter< false > >::Type IndexSetGetter;
      
    public: 
      //! \brief type of HierarchicIndexSet, default is LeafIndexSet
      typedef typename IndexSetGetter::IndexSet HierarchicIndexSet;
     
      //! \brief return reference to hierarchic index set 
      static const HierarchicIndexSet &hierarchicIndexSet ( const Grid &grid )
      { 
        return IndexSetGetter::indexSet( grid );
      }
    };



    //! Wraps HierarchicIndex Sets of AlbertaGrid and ALUGrid 
    /** \deprecated */
    template< class GridType >
    class WrappedHierarchicIndexSet
    : public IndexSetWrapper< typename HierarchicIndexSetSelector< GridType >::HierarchicIndexSet >
    {
      typedef WrappedHierarchicIndexSet< GridType > ThisType;
      typedef IndexSetWrapper< typename HierarchicIndexSetSelector< GridType >::HierarchicIndexSet > BaseType;

      // my type, to be revised 
      enum { myType = 0 };

      // type of hierarchic index set selector
      typedef HierarchicIndexSetSelector< GridType > SelectorType;
      
    public:
      //! number of codimensions 
      enum { ncodim = GridType::dimension + 1 };

      //! constructor 
      WrappedHierarchicIndexSet ( const GridType &grid, const int level =-1 )
      : BaseType( SelectorType::hierarchicIndexSet( grid ) )
      {}
         
      //! return type (for Grape In/Output)
      static int type ()
      {
        return myType;
      }

      //! returns reference to singleton
      static ThisType &instance ( const GridType &grid )
      {
        static ThisType set( grid );
        std::cerr << "Warning: WrappedHierarchicIndexSet::instance( grid ) can only handle one grid." << std::endl;
        return set;
      }
      //! return persistent status 
      bool persistent () const { return true; }
    };



    //! Wraps LeafIndexSet of Dune Grids for use with LagrangeFunctionSpace 
    template< class GridType >
    class WrappedLeafIndexSet
    : public IndexSetWrapper< typename GridType::Traits::LeafIndexSet >
    {
      typedef WrappedLeafIndexSet< GridType > ThisType;
      typedef IndexSetWrapper< typename GridType::Traits::LeafIndexSet > BaseType;

      // my type, to be revised 
      enum { myType = 5 };

    public:
      //! number of codimensions 
      enum { ncodim = GridType::dimension + 1 };

      //! constructor 
      WrappedLeafIndexSet ( const GridType &grid )
      : BaseType( grid.leafIndexSet() )
      {}

      WrappedLeafIndexSet ( const GridType &grid, const int level ) DUNE_DEPRECATED
      : BaseType( grid.leafIndexSet() )
      {}

      //! constructor taking grid part 
      template< class GridPartType >
      WrappedLeafIndexSet ( const GridPartType &gridPart )
      : BaseType( gridPart.grid().leafIndexSet() )
      {}

      //! return type (for Grape In/Output)
      static int type() { return myType; }
    };

  } // namespace Fem

} // namespace Dune 

#endif // #ifndef DUNE_FEM_DEFAULTINDEXSETS_HH
