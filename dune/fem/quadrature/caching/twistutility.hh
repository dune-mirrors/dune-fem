#ifndef DUNE_TWISTUTILITY_HH
#define DUNE_TWISTUTILITY_HH

#include <cassert>

#include <dune/common/static_assert.hh>
#include <dune/common/version.hh>
#if HAVE_DUNE_GEOMETRY
#include <dune/geometry/type.hh>
#include <dune/geometry/genericgeometry/topologytypes.hh>
#else
#include <dune/common/geometrytype.hh>
#include <dune/grid/genericgeometry/topologytypes.hh>
#endif

#include <dune/grid/alugrid/common/interfaces.hh>

#if HAVE_DUNE_GEOGRID
#include <dune/grid/utility/hostgridaccess.hh>
#endif

// this also includes the forward declarations 
#include <dune/fem/misc/capabilities.hh>

namespace Dune
{

  /** \brief TwistFreeTwistUtility provides the default implementation for twistfree grid
       such as Cartesian grids. 
  */
  template< class Grid > 
  struct TwistFreeTwistUtility
  {
    typedef Grid GridType;

    //! \brief return 0 for inner face 
    template< class Intersection >
    static int twistInSelf ( const GridType &, const Intersection & )
    {
      return 0;
    }
    
    //! \brief return 0 for outer face 
    template< class Intersection >
    static int twistInNeighbor ( const GridType &, const Intersection & )
    {
      return 0;
    }

    /** \brief return geometry type of inside or outside entity */
    template< class Intersection >
    static GeometryType elementGeometry ( const Intersection &intersection, const bool inside )
    {
      typedef Capabilities::hasSingleGeometryType< GridType > hasSingleGeometryType;
      if( hasSingleGeometryType::v )
        return GeometryType( hasSingleGeometryType::topologyId, GridType::dimension );
      else
        return (inside ? intersection.inside()->type() : intersection.outside()->type());
    }
  };



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

      \note This class needs to be explicitly specialized for each grid! 
  */
  template< class Grid > 
  struct TwistUtility;
  //{
  //  dune_static_assert( false, "TwistUtility not specialized, please specialize TwistUtility!");
  //};

  // Specialization for YaspGrid
  // ------------------------------
  
  template< int dimw  > 
  struct TwistUtility< YaspGrid< dimw > > 
    : public TwistFreeTwistUtility< YaspGrid< dimw > >
  {
  };

  
  // Specialization for SGrid
  // ------------------------------
  
  template< int dim, int dimworld, class ctype > 
  struct TwistUtility< SGrid<dim, dimworld, ctype> > 
    : public TwistFreeTwistUtility< SGrid<dim, dimworld, ctype> >
  {
  };


  // Specialization for OneDGrid
  // ------------------------------
  
  template <> 
  struct TwistUtility< OneDGrid > 
    : public TwistFreeTwistUtility< OneDGrid >
  {
  };


  // Specialization for AlbertaGrid
  // ------------------------------
  
#if HAVE_ALBERTA
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

  public:
    //! \brief return twist for inner face 
    static int twistInSelf ( const GridType &grid, const LeafIntersection &it )
    {
      return grid.getTwistInInside( it );
    }
    
    //! \brief return twist for outer face 
    static int twistInNeighbor ( const GridType &grid, const LeafIntersection &it )
    {
      return grid.getTwistInOutside( it );
    }

    /** \brief return element geometry type of inside or outside entity 
    */
    template <class Intersection>  
    static inline GeometryType
    elementGeometry(const Intersection& intersection, 
                    const bool inside)
    {
      return GeometryType( GenericGeometry::SimplexTopology< dimension >::type::id,
                           dimension );
    }
  };
#endif



  // Specialization for ALUGrid
  // --------------------------

#if HAVE_ALUGRID
  /** \brief Specialization of TwistUtility for ALUGridSimplex. 
  */
  template< int dim, int dimw >
  struct TwistUtility< ALUSimplexGrid< dim, dimw > >
  {
    typedef ALUSimplexGrid< dim, dimw > GridType;

    typedef typename GridType::Traits::LeafIntersectionIterator
      LeafIntersectionIterator;
    typedef typename LeafIntersectionIterator::Intersection LeafIntersection;

    typedef typename GridType::Traits::LevelIntersectionIterator
      LevelIntersectionIterator;
    typedef typename LevelIntersectionIterator::Intersection LevelIntersection;

    static const int dimension = GridType::dimension;

  public:
    //! \brief return twist for inner face 
    template <class Intersection> 
    static inline int twistInSelf(const GridType & grid, 
                                  const Intersection& intersection)
    {
      return grid.getRealIntersection( intersection ).twistInInside();
    }
    
    //! \brief return twist for outer face 
    template <class Intersection>
    static inline int twistInNeighbor(const GridType & grid, 
                                      const Intersection& intersection) 
    {
      return grid.getRealIntersection( intersection ).twistInOutside();
    }
    
    /** \brief return element geometry type of inside or outside entity 
    */
    template <class Intersection>  
    static inline GeometryType
    elementGeometry(const Intersection& intersection, 
                    const bool inside)
    {
      return GeometryType( GenericGeometry::SimplexTopology< dim >::type::id,
                           dim );
    }
  private:
    TwistUtility(const TwistUtility&);
    TwistUtility& operator=(const TwistUtility&);
  };

  /** \brief Specialization of TwistUtility for ALUSimplexGrid. 
  */
  template< int dim, int dimw >
  struct TwistUtility< ALUCubeGrid< dim, dimw > >
  {
    typedef ALUCubeGrid< dim, dimw > GridType;
    typedef typename GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename LeafIntersectionIterator::Intersection LeafIntersection;
    typedef typename GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
    typedef typename LevelIntersectionIterator::Intersection LevelIntersection;

  public:
    //! \brief return twist for inner face 
    template <class Intersection> 
    static inline int twistInSelf(const GridType & grid, 
                                  const Intersection& intersection)
    {
      return grid.getRealIntersection( intersection ).twistInInside();
    }
    
    //! \brief return twist for outer face 
    template <class Intersection>
    static inline int twistInNeighbor(const GridType & grid, 
                                      const Intersection& intersection) 
    {
      return grid.getRealIntersection( intersection ).twistInOutside();
    }
    
    /** \brief return element geometry type of inside or outside entity 
    */
    template <class Intersection>  
    static inline GeometryType
    elementGeometry(const Intersection& intersection, 
                    const bool inside)
    {
      return GeometryType( GenericGeometry::CubeTopology< dim >::type::id,
                           dim );
    }

  private:
    TwistUtility(const TwistUtility&);
    TwistUtility& operator=(const TwistUtility&);
  };

  /** \brief Specialization of TwistUtility for ALUConformGrid. 
  */
  template< int dim, int dimw >
  struct TwistUtility< ALUConformGrid< dim, dimw > >
  {
    typedef ALUConformGrid< dim, dimw > GridType;
    typedef typename GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename LeafIntersectionIterator::Intersection LeafIntersection;
    typedef typename GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
    typedef typename LevelIntersectionIterator::Intersection LevelIntersection;

  public:
    //! \brief return twist for inner face 
    static inline int twistInSelf(const GridType & grid, const LeafIntersection& intersection)
    {
      return grid.getRealIntersection( intersection ).twistInInside();
    }
    
    //! \brief return twist for outer face 
    static inline int twistInNeighbor(const GridType &grid, const LeafIntersection& intersection )
    {
      return grid.getRealIntersection( intersection ).twistInOutside();
    }
    
    /** \brief return element geometry type of inside or outside entity 
    */
    template <class Intersection>  
    static inline GeometryType
    elementGeometry(const Intersection& intersection, 
                    const bool inside)
    {
      return GeometryType( GenericGeometry::SimplexTopology< dim >::type::id,
                           dim );
    }
    
  private:
    TwistUtility(const TwistUtility&);
    TwistUtility& operator=(const TwistUtility&);
  };


  /** \brief Specialization of TwistUtility for ALUGrid. 
  */
  template< int dim, int dimw, ALUGridElementType elType, ALUGridRefinementType refineType, class Comm >
  struct TwistUtility< ALUGrid< dim, dimw, elType, refineType, Comm > >
  {
    typedef ALUGrid< dim, dimw, elType, refineType, Comm > GridType;
    typedef typename GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename LeafIntersectionIterator::Intersection LeafIntersection;
    typedef typename GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
    typedef typename LevelIntersectionIterator::Intersection LevelIntersection;

  public:
    //! \brief return twist for inner face 
    static inline int twistInSelf(const GridType & grid, const LeafIntersection& intersection)
    {
      return grid.getRealIntersection( intersection ).twistInInside();
    }
    
    //! \brief return twist for outer face 
    static inline int twistInNeighbor(const GridType &grid, const LeafIntersection& intersection )
    {
      return grid.getRealIntersection( intersection ).twistInOutside();
    }
    
    /** \brief return element geometry type of inside or outside entity 
    */
    template <class Intersection>  
    static inline GeometryType
    elementGeometry(const Intersection& intersection, 
                    const bool inside)
    {
      return GeometryType( Capabilities :: hasSingleGeometryType< GridType > :: topologyId, 
                           dim );
    }
    
  private:
    TwistUtility(const TwistUtility&);
    TwistUtility& operator=(const TwistUtility&);
  };
#endif



  // Specialization for UGGrid
  // -------------------------

#if HAVE_UG
  template< int dim >
  struct TwistUtility< UGGrid< dim > >
  {
    typedef UGGrid< dim > GridType;

    typedef typename GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename LeafIntersectionIterator::Intersection LeafIntersection;
    typedef typename GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
    typedef typename LevelIntersectionIterator::Intersection LevelIntersection;

    static int twistInSelf ( const GridType &grid, const LeafIntersection &it );
    static int twistInSelf ( const GridType &grid, const LevelIntersection &it );

    static int twistInNeighbor ( const GridType &grid, const LeafIntersection &it );
    static int twistInNeighbor ( const GridType &grid, const LevelIntersection &it );

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

  public:
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

  private:
    TwistUtility( const TwistUtility & );
    TwistUtility &operator=( const TwistUtility & );
  };
#endif // #if HAVE_DUNE_GEOGRID
  
} // end namespace Dune 

#endif // #ifndef DUNE_TWISTUTILITY_HH
