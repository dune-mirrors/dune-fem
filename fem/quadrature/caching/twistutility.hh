#ifndef DUNE_TWISTUTILITY_HH
#define DUNE_TWISTUTILITY_HH

#include <cassert>

#include <dune/common/interfaces.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/version.hh>
#include <dune/common/geometrytype.hh>

#include <dune/grid/common/capabilities.hh>

#if HAVE_DUNE_GEOGRID
#include <dune/grid/utility/hostgridaccess.hh>
#endif

namespace Dune
{

  /** \brief Utility to get twist from IntersectionIterator, 
      if provided by grid (i.e. AlbertaGrid, ALUGrid)
      otherwise return default values (correct for YASP/SGRID).
      
      The twist (t) of a face is defined in the following way:
      - sign(t) gives information on the relationship between the
        orientation of the intersection geometry and the geometry
        of the corresponding codim 1 entity of the inside/outside
        entity:
        - sign(t)>=0: same orientation
        - sign(t)<0:  opposite orientation

      - The value of the twist gives information on the local numbering
        of the corners of the corresponding geometries. This value
        is only correctly defined for conforming grids, i.e.,
        the intersection is identical to an codim 1 entity of inside/outside.
        In this case we have the following definition:
        - sign(t)>=0: corner[0] of inside/outside face is equal to
                      corner[t] of intersection.
        - sign(t)<0:  corner[0] of inside/outside face is equal to
                      corner[t'] of intersection with t' = abs(t)+1.
  */
  template< class Grid > 
  class TwistUtility
  {
    dune_static_assert( (!Conversion< Grid, HasHierarchicIndexSet > :: exists),
                        "The default TwistUtility is only for SGrid, "
			"YaspGrid, UGGrid and OneDGrid." );

  public:
    typedef Grid GridType;

  public:
    //! \brief constructor taking grid reference 
    TwistUtility(const GridType& grid) {}

    //! \brief return 0 for inner face 
    template <class IntersectionIterator> 
    static inline int twistInSelf(const GridType &, const IntersectionIterator&)
    {
      return 0;
    }
    
    //! \brief return 0 for inner face 
    template <class IntersectionIterator> 
    int twistInSelf(const IntersectionIterator& it) const {
      return 0;
    }
    
    //! \brief return 0 for outer face 
    template <class IntersectionIterator> 
    int twistInNeighbor(const IntersectionIterator& it) const {
      return 0;
    }

    //! \brief return 0 for outer face 
    template <class IntersectionIterator> 
    static inline int twistInNeighbor(const GridType &, const IntersectionIterator&) 
    {
      return 0;
    }

    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    bool DUNE_DEPRECATED conforming (const IntersectionIterator& it) const { return true; }
    
    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    static bool DUNE_DEPRECATED
    conforming (const GridType &, const IntersectionIterator&) { return true; }

    /** \brief return geometry type of inside or outside entity */
    template <class Intersection>  
    static inline GeometryType
    elementGeometry(const Intersection& intersection, 
                    const bool inside) 
    {
      if( Capabilities::IsUnstructured<GridType>::v ) 
      {
        return (inside) ? intersection.inside()->type() :  
                          intersection.outside()->type();  
      }
      else 
      {
#ifndef NDEBUG
        GeometryType geoType( GeometryType::cube, GridType :: dimension );
        GeometryType realType = (inside) ? intersection.inside()->type() :
                                           intersection.outside()->type();
        assert ( realType == geoType );
#endif
        return GeometryType( GeometryType::cube, GridType :: dimension );
      }
    }
  };



  // Specialization for AlbertaGrid
  // ------------------------------
  
#ifdef ENABLE_ALBERTA
  template< int dim, int dimW >
  class AlbertaGrid;
  
  /** \brief Specialization of TwistUtility for AlbertaGrid. 
  */
  template< int dim, int dimW >
  struct TwistUtility< AlbertaGrid< dim, dimW > >
  {
    typedef AlbertaGrid<dim, dimW> GridType;
    typedef typename GridType::Traits::LeafIntersectionIterator  LeafIntersectionIterator;
    typedef typename LeafIntersectionIterator::Intersection  LeafIntersection;
    typedef typename GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
    typedef typename LevelIntersectionIterator::Intersection LevelIntersection;

    static const int dimension = GridType::dimension;

  private:
    const GridType &grid_;

  public:
    //! \brief constructor taking grid reference 
    TwistUtility ( const GridType &grid )
    : grid_( grid )
    {}

    //! \brief return twist for inner face 
    static int twistInSelf ( const GridType &grid, const LeafIntersection &it )
    {
      const int map3d[ 6 ] = {-2, -3, -1, 0, 2, 1};
      const int twist = grid.getTwistInInside( it );
      return (dimension == 3 ? map3d[ twist + 3 ] : twist);
    }
    
    //! \brief return twist for inner face 
    int twistInSelf ( const LeafIntersection &it ) const
    {
      return twistInSelf( grid_, it );
    }
    
    //! \brief return twist for outer face 
    static int twistInNeighbor ( const GridType &grid, const LeafIntersection &it )
    {
      const int map3d[ 6 ] = {-2, -3, -1, 0, 2, 1};
      const int twist = grid.getTwistInOutside( it );
      return (dimension == 3 ? map3d[ twist + 3 ] : twist);
    }

    //! \brief return twist for outer face 
    int twistInNeighbor ( const LeafIntersection &it ) const
    {
      return twistInNeighbor( grid_, it );
    }
    
    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    bool DUNE_DEPRECATED conforming (const IntersectionIterator& it) const
    {
      return true;
    }
    
    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    static bool DUNE_DEPRECATED
    conforming (const GridType & grid, const IntersectionIterator& it)
    {
      return true;
    }
    
    /** \brief return element geometry type of inside or outside entity 
    */
    template <class Intersection>  
    static inline GeometryType
    elementGeometry(const Intersection& intersection, 
                    const bool inside)
    {
      return GeometryType( GeometryType::simplex, dimension );
    }
  };
#endif



  // Specialization for ALUGrid
  // --------------------------

#ifdef ENABLE_ALUGRID
  template< int dim, int dimW >
  class ALUSimplexGrid;

  template< int dim, int dimW >
  class ALUCubeGrid;

  template< int dim, int dimW >
  class ALUConformGrid;

  /** \brief Specialization of TwistUtility for ALUGridSimplex. 
  */
  template< int dim >
  struct TwistUtility< ALUSimplexGrid< dim, dim > >
  {
    typedef ALUSimplexGrid< dim, dim > GridType;

    typedef typename GridType::Traits::LeafIntersectionIterator
      LeafIntersectionIterator;
    typedef typename LeafIntersectionIterator::Intersection LeafIntersection;

    typedef typename GridType::Traits::LevelIntersectionIterator
      LevelIntersectionIterator;
    typedef typename LevelIntersectionIterator::Intersection LevelIntersection;

    static const int dimension = GridType::dimension;

  public:
    //! \brief constructor taking grid reference 
    TwistUtility(const GridType& grid) :
      grid_(grid)
    {}

    //! \brief return twist for inner face 
    template <class IntersectionIterator> 
    static inline int twistInSelf(const GridType & grid, 
                           const IntersectionIterator& it)
    {
      return grid.getRealIntersectionIterator(it).twistInSelf();
    }
    
    //! \brief return twist for inner face 
    int twistInSelf(const LeafIntersection& it) const {
      return grid_.getRealIntersectionIterator(it).twistInSelf();
    }
    
    //! \brief return twist for inner face 
    int twistInSelf(const LevelIntersection& it) const {
      return grid_.getRealIntersectionIterator(it).twistInSelf();
    }

    //! \brief return twist for outer face 
    int twistInNeighbor(const LeafIntersection& it) const {
      return grid_.getRealIntersectionIterator(it).twistInNeighbor();
    }
    
    //! \brief return twist for outer face 
    int twistInNeighbor(const LevelIntersection& it) const {
      return grid_.getRealIntersectionIterator(it).twistInNeighbor();
    }

    //! \brief return twist for outer face 
    template <class IntersectionIterator>
    static inline int twistInNeighbor(const GridType & grid, 
                               const IntersectionIterator& it) 
    {
      return grid.getRealIntersectionIterator(it).twistInNeighbor();
    }
    
    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    bool DUNE_DEPRECATED conforming (const IntersectionIterator& it) const
    { 
      return grid_.getRealIntersectionIterator(it).conforming(); 
    }

    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    static bool DUNE_DEPRECATED
    conforming (const GridType & grid, const IntersectionIterator& it)
    { 
      return grid.getRealIntersectionIterator(it).conforming(); 
    }
    
    /** \brief return element geometry type of inside or outside entity 
    */
    template <class Intersection>  
    static inline GeometryType
    elementGeometry(const Intersection& intersection, 
                    const bool inside)
    {
      return GeometryType( GeometryType::simplex, dimension );
    }
  private:
    TwistUtility(const TwistUtility&);
    TwistUtility& operator=(const TwistUtility&);
    
  private:
    const GridType& grid_; 
  };

  /** \brief Specialization of TwistUtility for ALUGridSimplex. 
  */
  template< int dim >
  struct TwistUtility< ALUCubeGrid< dim, dim > >
  {
    typedef ALUCubeGrid< dim, dim > GridType;
    typedef typename GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename LeafIntersectionIterator::Intersection LeafIntersection;
    typedef typename GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
    typedef typename LevelIntersectionIterator::Intersection LevelIntersection;

  public:
    //! \brief constructor taking grid reference 
    TwistUtility(const GridType& grid) :
      grid_(grid)
    {}

    //! \brief return twist for inner face 
    template <class IntersectionIterator> 
    static inline int twistInSelf(const GridType & grid, 
                           const IntersectionIterator& it)
    {
      return grid.getRealIntersectionIterator(it).twistInSelf();
    }
    
    //! \brief return twist for inner face 
    int twistInSelf(const LeafIntersection& it) const {
      return grid_.getRealIntersectionIterator(it).twistInSelf();
    }
    
    //! \brief return twist for inner face 
    int twistInSelf(const LevelIntersection& it) const {
      return grid_.getRealIntersectionIterator(it).twistInSelf();
    }

    //! \brief return twist for outer face 
    int twistInNeighbor(const LeafIntersection& it) const {
      return grid_.getRealIntersectionIterator(it).twistInNeighbor();
    }
    
    //! \brief return twist for outer face 
    template <class IntersectionIterator>
    static inline int twistInNeighbor(const GridType & grid, 
                               const IntersectionIterator& it) 
    {
      return grid.getRealIntersectionIterator(it).twistInNeighbor();
    }
    
    //! \brief return twist for outer face 
    int twistInNeighbor(const LevelIntersection& it) const {
      return grid_.getRealIntersectionIterator(it).twistInNeighbor();
    }

    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    bool DUNE_DEPRECATED conforming (const IntersectionIterator& it) const
    { 
      return grid_.getRealIntersectionIterator(it).conforming(); 
    }

    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    static bool DUNE_DEPRECATED
    conforming (const GridType & grid, const IntersectionIterator& it)
    { 
      return grid.getRealIntersectionIterator(it).conforming(); 
    }
    
    /** \brief return element geometry type of inside or outside entity 
    */
    template <class Intersection>  
    static inline GeometryType
    elementGeometry(const Intersection& intersection, 
                    const bool inside)
    {
      return GeometryType( GeometryType::cube, GridType :: dimension );
    }
  private:
    TwistUtility(const TwistUtility&);
    TwistUtility& operator=(const TwistUtility&);
    
  private:
    const GridType& grid_; 
  };

  template< int dim >
  struct TwistUtility< ALUConformGrid< dim, dim > >
  {
    typedef ALUConformGrid< dim, dim > GridType;
    typedef typename GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename LeafIntersectionIterator::Intersection LeafIntersection;
    typedef typename GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
    typedef typename LevelIntersectionIterator::Intersection LevelIntersection;

  public:
    //! \brief constructor taking grid reference 
    TwistUtility(const GridType& grid) : grid_(grid) {}

    //! \brief return twist for inner face 
    static inline int twistInSelf(const GridType &, const LeafIntersection&)
    {
      return 0;
    }
    
    //! \brief return twist for inner face 
    int twistInSelf(const LeafIntersection& it) const {
      return 0;
    }
    
    //! \brief return twist for inner face 
    int twistInSelf(const LevelIntersection& it) const {
      return 0;
    }

    //! \brief return twist for outer face 
    static inline int twistInNeighbor(const GridType &, const LeafIntersection&)
    {
      return 1;
    }
    
    //! \brief return twist for outer face 
    int twistInNeighbor(const LeafIntersection& it) const {
      return 1;
    }
    
    //! \brief return twist for outer face 
    int twistInNeighbor(const LevelIntersection& it) const {
      return 1;
    }

    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    bool DUNE_DEPRECATED conforming (const IntersectionIterator& it) const
    { 
      return grid_.getRealIntersectionIterator(it).conforming(); 
    }
    
    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    static bool DUNE_DEPRECATED
    conforming (const GridType & grid, const IntersectionIterator& it)
    { 
      return grid.getRealIntersectionIterator(it).conforming(); 
    }
    
    /** \brief return element geometry type of inside or outside entity 
    */
    template <class Intersection>  
    static inline GeometryType
    elementGeometry(const Intersection& intersection, 
                    const bool inside)
    {
      return GeometryType( GeometryType::simplex, GridType :: dimension );
    }
    
  private:
    TwistUtility(const TwistUtility&);
    TwistUtility& operator=(const TwistUtility&);

  private:
    const GridType& grid_; 
  };
#endif



  // Specialization for UGGrid
  // -------------------------

#ifdef ENABLE_UG
  template< int dim >
  class UGGrid;

  template< int dim >
  struct TwistUtility< UGGrid< dim > >
  {
    typedef UGGrid< dim > GridType;

    typedef typename GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename LeafIntersectionIterator::Intersection LeafIntersection;
    typedef typename GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
    typedef typename LevelIntersectionIterator::Intersection LevelIntersection;

    explicit TwistUtility ( const GridType &grid )
    {}

    static int twistInSelf ( const GridType &grid, const LeafIntersection &it );
    int twistInSelf( const LeafIntersection &it ) const;
    static int twistInSelf ( const GridType &grid, const LevelIntersection &it );
    int twistInSelf( const LevelIntersection &it ) const;

    static int twistInNeighbor ( const GridType &grid, const LeafIntersection &it );
    int twistInNeighbor( const LeafIntersection &it ) const;
    static int twistInNeighbor ( const GridType &grid, const LevelIntersection &it );
    int twistInNeighbor( const LevelIntersection &it ) const;

    static bool DUNE_DEPRECATED conforming ( const GridType &grid, const LeafIntersection &it );
    bool DUNE_DEPRECATED conforming ( const LeafIntersection &it ) const;
    static bool DUNE_DEPRECATED conforming ( const GridType &grid, const LevelIntersection &it );
    bool DUNE_DEPRECATED conforming ( const LevelIntersection &it ) const;

    template< class Intersection >
    static GeometryType
    elementGeometry ( const Intersection &intersection, const bool inside )
    {
      return (inside ? intersection.inside()->type() : intersection.outside()->type());
    }
  };
#endif // #ifdef ENABLE_UG



  // Specialization for GeoGrid
  // --------------------------

#if HAVE_DUNE_GEOGRID
  template< class HostGrid, class CoordFunction >
  class GeometryGrid;

  template< class HostGrid, class CoordFunction >
  struct TwistUtility< GeometryGrid< HostGrid, CoordFunction > >
  {
    typedef GeometryGrid< HostGrid, CoordFunction > GridType;
    typedef typename GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename LeafIntersectionIterator::Intersection LeafIntersection;
    typedef typename GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
    typedef typename LevelIntersectionIterator::Intersection LevelIntersection;

  private:
    typedef TwistUtility< HostGrid > HostTwistUtility;
    typedef Dune::HostGridAccess< GridType > HostGridAccess;

    const GridType &grid_;

  public:
    //! \brief constructor
    TwistUtility ( const GridType &grid )
    : grid_( grid )
    {}

    //! \brief return twist for inner face
    static int twistInSelf ( const GridType &grid,
                             const LeafIntersection &intersection )
    {
      typedef typename HostGridAccess :: HostLeafIntersection HostIntersection;
      const HostIntersection &hostIntersection
        = HostGridAccess :: getIntersection( intersection );
      return HostTwistUtility :: twistInSelf( grid.hostGrid(), hostIntersection );
    }

    //! \brief return twist for inner face
    static int twistInSelf ( const GridType &grid,
                             const LevelIntersection &intersection )
    {
      typedef typename HostGridAccess :: HostLevelIntersection HostIntersection;
      const HostIntersection &hostIntersection
        = HostGridAccess :: getIntersection( intersection );
      return HostTwistUtility :: twistInSelf( grid.hostGrid(), hostIntersection );
    }

    //! \brief return twist for inner face
    template< class G, template< class > class II, template< class > class I >
    static int DUNE_DEPRECATED
    twistInSelf ( const GridType &grid,
                  const IntersectionIterator< G, II, I > &intersectionIterator )
    {
      return twistInSelf( grid, *intersectionIterator );
    }

    //! \brief return twist for inner face
    template< class Intersection >
    int twistInSelf ( const Intersection &intersection ) const
    {
      return twistInSelf( grid_, intersection );
    }


    //! \brief return twist for outer face
    static int twistInNeighbor ( const GridType &grid,
                                 const LeafIntersection &intersection )
    {
      typedef typename HostGridAccess :: HostLeafIntersection HostIntersection;
      const HostIntersection &hostIntersection
        = HostGridAccess :: getIntersection( intersection );
      return HostTwistUtility :: twistInNeighbor( grid.hostGrid(), hostIntersection );
    }

    //! \brief return twist for outer face
    static int twistInNeighbor ( const GridType &grid,
                                 const LevelIntersection &intersection )
    {
      typedef typename HostGridAccess :: HostLevelIntersection HostIntersection;
      const HostIntersection &hostIntersection
        = HostGridAccess :: getIntersection( intersection );
      return HostTwistUtility :: twistInNeighbor( grid.hostGrid(), hostIntersection );
    }

    //! \brief return twist for inner face
    template< class G, template< class > class II, template< class > class I >
    static int DUNE_DEPRECATED
    twistInNeighbor ( const GridType &grid,
                      const IntersectionIterator< G, II, I > &intersectionIterator )
    {
      return twistInNeighbor( grid, *intersectionIterator );
    }

    //! \brief return twist for outer face
    template< class Intersection >
    int twistInNeighbor ( const Intersection &intersection ) const
    {
      return twistInNeighbor( grid_, intersection );
    }


    //! \brief return true if intersection is conform
    static bool DUNE_DEPRECATED
    conforming ( const GridType &grid, const LeafIntersection &intersection )
    {
      typedef typename HostGridAccess :: HostLeafIntersection HostIntersection;
      const HostIntersection &hostIntersection
        = HostGridAccess :: getIntersection( intersection );
      return HostTwistUtility :: conforming( grid.hostGrid(), hostIntersection );
    }

    //! \brief return true if intersection is conform
    static bool DUNE_DEPRECATED
    conforming ( const GridType &grid, const LevelIntersection &intersection )
    {
      typedef typename HostGridAccess :: HostLevelIntersection HostIntersection;
      const HostIntersection &hostIntersection
        = HostGridAccess :: getIntersection( intersection );
      return HostTwistUtility :: conforming( grid.hostGrid(), hostIntersection );
    }

    //! \brief return true if intersection is conform
    template< class G, template< class > class II, template< class > class I >
    static bool DUNE_DEPRECATED
    conforming ( const GridType &grid,
                 const IntersectionIterator< G, II, I > &intersectionIterator )
    {
      return conforming( grid, *intersectionIterator );
    }

    //! \brief return true if intersection is conform
    template< class Intersection >
    bool DUNE_DEPRECATED conforming ( const Intersection &intersection ) const
    {
      return conforming( grid_, intersection );
    }


    /** \brief return element geometry type of inside or outside entity */
    static GeometryType
    elementGeometry ( const LeafIntersection &intersection, bool inside )
    {
      typedef typename HostGridAccess :: HostLeafIntersection HostIntersection;
      const HostIntersection &hostIntersection
        = HostGridAccess :: getIntersection( intersection );
      return HostTwistUtility :: elementGeometry( hostIntersection, inside );
    }

    /** \brief return element geometry type of inside or outside entity */
    static GeometryType
    elementGeometry ( const LevelIntersection &intersection, bool inside )
    {
      typedef typename HostGridAccess :: HostLevelIntersection HostIntersection;
      const HostIntersection &hostIntersection
        = HostGridAccess :: getIntersection( intersection );
      return HostTwistUtility :: elementGeometry( hostIntersection, inside );
    }

    /** \brief return element geometry type of inside or outside entity */
    template< class G, template< class > class II, template< class > class I >
    static GeometryType DUNE_DEPRECATED
    elementGeometry ( const IntersectionIterator< G, II, I > &intersectionIterator,
                      bool inside )
    {
      return elementGeometry( *intersectionIterator, inside );
    }

  private:
    TwistUtility( const TwistUtility & );
    TwistUtility &operator=( const TwistUtility & );
  };
#endif // #if HAVE_DUNE_GEOGRID
  
} // end namespace Dune 

#endif // #ifndef DUNE_TWISTUTILITY_HH
