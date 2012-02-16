#ifndef DUNE_FEM_GRIDPART_COMMON_GRIDPART_HH
#define DUNE_FEM_GRIDPART_COMMON_GRIDPART_HH

//- dune-common includes
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/deprecated.hh>

//- dune-grid includes
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/entity.hh>
#include <dune/grid/common/grid.hh>

//- dune-fem includes
#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/gridpart/common/gridpartview.hh>

namespace Dune
{

  /** \addtogroup GridPart
   *
   * Grid parts allow to define a view on a given DUNE grid, treating the
   * underlying grid as a container for entities.
   *
   * All parts of the dune-fem package rely on grid parts to access the entities
   * of the grid. For example, discrete functions are defined on the set of
   * entities accesseable by the given GridPart implementation using the iterator
   * and index set provided by the GridPart.
   *
   * \section GridPart Interface and available Implementations
   * 
   * The interface for a GridPart is implemented by the class template
   * GridPartInterface. Basically, a GridPart provides the following
   * functionality:
   * - The underlying grid can be accessed through the grid method.
   * - The indexSet method provides a suitable dune-fem index set for the grid
   *   part.
   * - Pairs of begin / end methods provide iterators over the entities of a
   *   given codimension belonging to the grid part.
   * - A pair of ibegin / iend methods provide suitable intersection iterators
   *   for a given entity of codimension 0.
   * - For parallel computations, a suitable communicate method is provided.
   * .
   *
   * The following grid parts have been implemented:
   * - LeafGridPart: A view of the leaf grid,
   * - LevelGridPart: A view of a given grid level,
   * - FilteredGridPart: A view filtering another grid part.
   *
   * \todo Implement a grid part for a given grid view (Suggestion: use the
   *       name GridPart).
   */



  /** 
   * @addtogroup GridPart
   *
   * @{ 
   */

  //! \brief Interface for the GridPart classes
  //! A GridPart class allows to access only a specific subset of a grid's
  //! entities. A GridPart implementation provides the corresponding index set
  //! and a begin/end iterator pair for accessing those entities, the
  //! corresponding intersection iterators and a appropriate communication
  //! method. 
  //! GridParts are used to parametrize spaces (see DiscreteFunctionSpaceDefault [in dune-fem]).
  template< class GridPartTraits >
  class GridPartInterface
  {
    typedef GridPartInterface< GridPartTraits > ThisType;

  public:
    //! \brief Type of the Traits
    typedef GridPartTraits Traits;

    //! \brief Type of the implementation
    typedef typename Traits::GridPartType GridPartType;
   
    //! \brief type of Grid implementation
    typedef typename Traits::GridType GridType;

    //! \brief Index set implementation
    typedef typename Traits::IndexSetType IndexSetType;

    //! \brief Maximum Partition type, the index set provides indices for
    static const PartitionIteratorType indexSetPartitionType
      = Traits::indexSetPartitionType;
    static const InterfaceType indexSetInterfaceType
      = Traits::indexSetInterfaceType;

    //! \brief type of IntersectionIterator
    typedef typename Traits::IntersectionIteratorType IntersectionIteratorType;

    //! \brief type of Intersection
    typedef typename IntersectionIteratorType::Intersection IntersectionType;

    /** \brief is true if grid on this view only has conforming intersections;
     *         use the grid part capability 
     *         \code         
Dune::Fem::GridPartCapabilities::isConforming< GridPartType >::v
     *         \endcode
     *         instead.
     */
    static const bool conforming DUNE_DEPRECATED = Traits::conforming;

    typedef GridView< Fem::GridPartViewTraits< GridPartType > > GridViewType;

    typedef typename GridType::ctype ctype;

    static const int dimension = GridType::dimension;
    static const int dimensionworld = GridType::dimensionworld;

    template< int codim >
    struct Codim
    {
      typedef typename Traits::template Codim< codim >::GeometryType       GeometryType;
      typedef typename Traits::template Codim< codim >::LocalGeometryType  LocalGeometryType;

      typedef typename Traits::template Codim< codim >::EntityPointerType  EntityPointerType;
      typedef typename Traits::template Codim< codim >::EntityType         EntityType;

      template< PartitionIteratorType pitype >
      struct Partition
      {
        typedef typename Traits::template Codim< codim >::template Partition< pitype >::IteratorType
          IteratorType;
      };

      typedef typename Partition< InteriorBorder_Partition >::IteratorType IteratorType;
    };
    
  public:
    //! \brief Returns const reference to the underlying grid
    const GridType &grid () const
    { 
      CHECK_INTERFACE_IMPLEMENTATION((asImp().grid()));
      return asImp().grid(); 
    }
    //! \brief Returns reference to the underlying grid
    GridType &grid ()
    { 
      CHECK_INTERFACE_IMPLEMENTATION((asImp().grid()));
      return asImp().grid(); 
    }

    GridViewType gridView () const
    {
      typedef typename GridViewType::GridViewImp Impl;
      return GridViewType( Impl( asImp() ) );
    }
    
    //! \brief Returns reference to index set of the underlying grid
    const IndexSetType& indexSet() const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION((asImp().indexSet()));
      return asImp().indexSet(); 
    }

    /** \brief obtain begin iterator for the interior-border partition
     *
     *  \tparam  codim  codimension for which the iterator is requested
     */
    template< int codim >
    typename Codim< codim >::IteratorType
    begin () const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION( (asImp().template begin< codim >()) );
      return asImp().template begin< codim >();
    }

    /** \brief obtain begin iterator for the given partition
     *
     *  \tparam  codim   codimension for which the iterator is requested
     *  \tparam  pitype  requested partition iterator type
     */
    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim >::template Partition< pitype >::IteratorType
    begin () const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION( (asImp().template begin< codim, pitype >()) );
      return asImp().template begin< codim, pitype >(); 
    }

    /** \brief obtain end iterator for the interior-border partition
     *
     *  \tparam  codim  codimension for which the iterator is requested
     */
    template< int codim >
    typename Codim< codim >::IteratorType
    end () const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION( (asImp().template end< codim >()) );
      return asImp().template end< codim >();
    }

    /** \brief obtain end iterator for the given partition
     *
     *  \tparam  codim   codimension for which the iterator is requested
     *  \tparam  pitype  requested partition iterator type
     */
    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim >::template Partition< pitype >::IteratorType
    end () const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION( (asImp().template end< codim, pitype >()) );
      return asImp().template end< codim, pitype >();
    }

    //! \brief Level of the grid part
    int level () const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION((asImp().level()));
      return asImp().level(); 
    }

    //! \brief ibegin of corresponding intersection iterator for given entity
    IntersectionIteratorType
    ibegin ( const typename Codim< 0 >::EntityType &entity ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION( (asImp().ibegin( entity )) );
      return asImp().ibegin( entity ); 
    }
    
    //! \brief iend of corresponding intersection iterator for given entity
    IntersectionIteratorType iend ( const typename Codim< 0 >::EntityType &entity ) const 
    {
      CHECK_INTERFACE_IMPLEMENTATION( (asImp().iend( entity )) );
      return asImp().iend( entity ); 
    }

    int boundaryId ( const IntersectionType &intersection ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().boundaryId( intersection ) );
      return asImp().boundaryId( intersection );
    }

    //! \brief corresponding communication method for grid part
    template< class DataHandleImp, class DataType >
    void communicate ( CommDataHandleIF< DataHandleImp, DataType > &data,
                       InterfaceType iftype, CommunicationDirection dir ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( (asImp().communicate( data, iftype, dir )) );
    }

    /*! \brief convert the grid's entity to a grid part entity
        Usually the parameter is GridType :: Codim< codim > :: Entity  
        and the return is Codim< codim > :: EntityType. 
        In general these types are the same, but for overloaded entities on grid parts
        this can differ. 
      */
    template <class Entity> 
    const Entity& convert( const Entity& entity ) const 
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().convert( entity ) );
      return asImp().convert( entity );
    }

  protected: 
    //! do not create explicit instances of this class 
    GridPartInterface () {}  

  private:
    GridPartType &asImp () { return static_cast< GridPartType & >( *this ); }
    const GridPartType &asImp () const { return static_cast< const GridPartType & >( *this ); }
  };



  //! \brief Default implementation for the GridPart classes
  template< class GridPartTraits >
  class GridPartDefault
  : public GridPartInterface< GridPartTraits >
  {
    typedef GridPartDefault< GridPartTraits > ThisType;

  public:
    //! \brief Type of the Traits
    typedef GridPartTraits Traits;
    //! \brief Grid implementation
    typedef typename Traits::GridType GridType;
    //! \brief Index set implementation
    typedef typename Traits::IndexSetType IndexSetType;

  protected:
    GridType &grid_;

  protected:
    //! constructor
    GridPartDefault ( GridType &grid )
    : grid_( grid )
    {}

    GridPartDefault ( const ThisType &other )
    : grid_( other.grid_ )
    {}

    ~GridPartDefault ()
    {}

  public:  
    //! Returns const reference to the underlying grid
    const GridType &grid () const { return grid_; }

    //! Returns reference to the underlying grid
    GridType &grid () { return grid_; }

    /* \brief \copydoc GridPartInterface::convert 
       
       The default implementation does nothing but return the same entity 
     */
    template <class Entity> 
    const Entity& convert( const Entity& entity ) const 
    {
      return entity;
    }

  };


#undef CHECK_INTERFACE_IMPLEMENTATION
#undef CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
  /** @} */

  namespace Fem
  {

    template< class Entity >
    struct GridEntityAccess;

    template< int codim, int dim, class Grid, template< int, int, class > class EntityImpl >
    struct GridEntityAccess< Dune::Entity< codim, dim, Grid, EntityImpl > >
    {
      typedef Dune::Entity< codim, dim, Grid, EntityImpl > EntityType;
      typedef Dune::Entity< codim, dim, Grid, EntityImpl > GridEntityType;

      static const GridEntityType &gridEntity ( const EntityType &entity )
      {
        return entity;
      }
    };

    template< class Entity >
    const typename GridEntityAccess< Entity >::GridEntityType &
    gridEntity ( const Entity &entity )
    {
      return GridEntityAccess< Entity >::gridEntity( entity );
    }

  } // end namespace Fem

} // end namespace Dune

#endif // #define DUNE_FEM_GRIDPART_COMMON_GRIDPART_HH
