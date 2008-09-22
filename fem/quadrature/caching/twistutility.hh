#ifndef DUNE_TWISTUTILITY_HH
#define DUNE_TWISTUTILITY_HH

#include <cassert>

#include <dune/common/interfaces.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/geometrytype.hh>

#include <dune/grid/common/capabilities.hh>

#ifdef ENABLE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif

//#ifdef ENABLE_ALBERTA
//#include <dune/grid/albertagrid.hh>
//#endif

#ifdef ENABLE_UG
#include <dune/grid/uggrid.hh>
#endif

#if HAVE_DUNE_GEOGRID
#include <dune/grid/utility/hostgridaccess.hh>
#endif

//#include <dune/grid/sgrid.hh>
//#include <dune/grid/yaspgrid.hh>

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
    bool conforming (const IntersectionIterator& it) const { return true; }
    
    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    static inline bool conforming (const GridType &, const IntersectionIterator&) { return true; }

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
  
#ifdef ENABLE_ALBERTA
  template< int dim, int dimW >
  class AlbertaGrid;
  
  /** \brief Specialization of TwistUtility for AlbertaGrid. 
  */
  template <int dim, int dimW>
  class TwistUtility<AlbertaGrid<dim, dimW> >
  {
  public:
    typedef AlbertaGrid<dim, dimW> GridType;
    typedef typename GridType::Traits::LeafIntersectionIterator  LeafIntersectionIterator;
    typedef typename LeafIntersectionIterator::Intersection  LeafIntersection;
    typedef typename GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
    typedef typename LevelIntersectionIterator::Intersection LevelIntersection;
  public:
    //! \brief constructor taking grid reference 
    TwistUtility(const GridType& grid) :
      grid_(grid)
    {}

    //! \brief return twist for inner face 
    static inline int twistInSelf(const GridType & grid, 
                           const LeafIntersection& it) {
      return grid.getRealIntersectionIterator(it).twistInSelf();
    }
    
    //! \brief return twist for inner face 
    int twistInSelf(const LeafIntersection& it) const {
      return grid_.getRealIntersectionIterator(it).twistInSelf();
    }
    
    //int twistInSelf(const LevelIntersectionIterator& it) const {
    //  return grid_.getRealIntersectionIterator(it).twistInSelf();
    //}

    //! \brief return twist for outer face 
    int twistInNeighbor(const LeafIntersection& it) const {
      return grid_.getRealIntersectionIterator(it).twistInNeighbor();
    }
    
    //! \brief return twist for outer face 
    static inline int twistInNeighbor(const GridType & grid, 
                               const LeafIntersection& it)
    {
      return grid.getRealIntersectionIterator(it).twistInNeighbor();
    }
    
    //int twistInNeighbor(const LevelIntersectionIterator& it) const {
    //  return grid_.getRealIntersectionIterator(it).twistInNeighbor();
    //}

    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    bool conforming (const IntersectionIterator& it) const { return true; }
    
    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    static inline bool conforming (const GridType & grid, 
                            const IntersectionIterator& it) { return true; }
    
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
    const GridType& grid_;
  };
#endif

#ifdef ENABLE_ALUGRID
  /** \brief Specialization of TwistUtility for ALUGridSimplex. 
  */
  template <>
  class TwistUtility<ALUSimplexGrid<3,3>  >
  {
  public:
    typedef ALUSimplexGrid<3,3> GridType;
    typedef GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
    typedef LeafIntersectionIterator::Intersection LeafIntersection;
    typedef GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
    typedef LevelIntersectionIterator::Intersection LevelIntersection;
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
    bool conforming (const IntersectionIterator& it) const 
    { 
      return grid_.getRealIntersectionIterator(it).conforming(); 
    }

    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    static inline bool conforming (const GridType & grid, 
                     const IntersectionIterator& it)
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

  /** \brief Specialization of TwistUtility for ALUGridSimplex. 
  */
  template <>
  class TwistUtility<ALUCubeGrid<3,3>  >
  {
  public:
    typedef ALUCubeGrid<3,3> GridType;
    typedef GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
    typedef LeafIntersectionIterator::Intersection LeafIntersection;
    typedef GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
    typedef LevelIntersectionIterator::Intersection LevelIntersection;
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
    bool conforming (const IntersectionIterator& it) const 
    { 
      return grid_.getRealIntersectionIterator(it).conforming(); 
    }

    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    static inline bool conforming (const GridType & grid, 
                            const IntersectionIterator& it)
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
  
  /** \brief Specialization of TwistUtility for ALUGridSimplex. 
  */
  template <>
  class TwistUtility<ALUSimplexGrid<2,2>  >
  {
  public:
    typedef ALUSimplexGrid<2, 2> GridType;
    typedef GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
    typedef LeafIntersectionIterator::Intersection LeafIntersection;
    typedef GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
    typedef LevelIntersectionIterator::Intersection LevelIntersection;
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
    bool conforming (const IntersectionIterator& it) const 
    { 
      return grid_.getRealIntersectionIterator(it).conforming(); 
    }
    
    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    static inline bool conforming (const GridType & grid, 
                     const IntersectionIterator& it)
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
  
  template <>
  class TwistUtility<ALUConformGrid<2,2>  >
  {
  public:
    typedef ALUConformGrid<2, 2> GridType;
    typedef GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
    typedef LeafIntersectionIterator::Intersection LeafIntersection;
    typedef GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
    typedef LevelIntersectionIterator::Intersection LevelIntersection;
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
    bool conforming (const IntersectionIterator& it) const 
    { 
      return grid_.getRealIntersectionIterator(it).conforming(); 
    }
    
    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    static inline bool conforming (const GridType & grid, 
                     const IntersectionIterator& it)
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

#ifdef ENABLE_UG
  //! 2d twist utility for UGGrid.
  template <> 
  class TwistUtility< UGGrid<2> > 
  {
    enum { dim = 2 };
    typedef UGGrid<dim> GridImp;
    // this default implementation only is for SGrid, YaspGrid, UGGrid
    // and OneDGrid. 
  public:
    typedef GridImp GridType;
  public:
    //! \brief constructor taking grid reference 
    TwistUtility(const GridType& grid) {}

    //! \brief return 0 for inner face 
    template <class Intersection> 
    static inline int twistInSelf(const GridType &, const Intersection& it)
    {
      // for simplex twist is 0 
      // for cube twist is 1 for side 0 and 3 
      // for 1 and 2 is 0 
      return (it.inside()->type().isSimplex()) ? 0 : 
        (it.numberInSelf() == 1 || it.numberInSelf() == 2) ? 0 : 1;
    }    
    //! \brief return 0 for inner face 
    template <class Intersection> 
    int twistInSelf(const Intersection& it) const 
    {
      // for simplex twist is 0 
      // for cube twist is 1 for side 0 and 3 
      // for 1 and 2 is 0 
      return (it.inside()->type().isSimplex()) ? 0 : 
        (it.numberInSelf() == 1 || it.numberInSelf() == 2) ? 0 : 1;
    }
    
    //! \brief return 0 for outer face 
    template <class Intersection> 
    int twistInNeighbor(const Intersection& it) const {
      assert( it.neighbor() );
      return (it.outside()->type().isSimplex()) ? 1 : 
        (it.numberInNeighbor() == 1 || it.numberInNeighbor() == 2) ? 1 : 0;
    }

    //! \brief return 0 for outer face 
    template <class Intersection> 
    static inline int twistInNeighbor(const GridType &, const Intersection& it) 
    {
      assert( it.neighbor() );
      return (it.outside()->type().isSimplex()) ? 1 : 
        (it.numberInNeighbor() == 1 || it.numberInNeighbor() == 2) ? 1 : 0;
    }

    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    bool conforming (const IntersectionIterator& it) const 
    { 
      return (it.neighbor()) ? 
        (it.inside()->level() == it.outside()->level()) : true; 
    }
    
    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    static inline bool conforming (const GridType &, const IntersectionIterator& it)
    { 
      return (it.neighbor()) ? 
        (it.inside()->level() == it.outside()->level()) : true; 
    }
    
    /** \brief return geometry type of inside or outside entity */
    template <class Intersection>  
    static inline GeometryType
    elementGeometry(const Intersection& intersection, 
                    const bool inside) 
    {
      return (inside) ? intersection.inside()->type() :  
                        intersection.outside()->type();  
    }
  };

  //! 3d twist utility for UGGrid.
  template <> 
  class TwistUtility< UGGrid<3> > 
  {
    enum { dim = 3 };
    typedef UGGrid<dim> GridImp;
    // this default implementation only is for SGrid, YaspGrid, UGGrid
    // and OneDGrid. 
  public:
    typedef GridImp GridType;

    struct CubeTwists 
    {
      static inline int twistInNeighbor(const int face)
      {
        static const int twistInNeigh[6] = { 0, -2, -2,  0,  0, -1 };
        assert( face >= 0 && face < 6 );
        return twistInNeigh[face];
      }
      
      static inline int twistInSelf(const int face)
      {
        static const int twistInSelf[6]  = {-2,  0,  0, -2, -1,  0 };
        assert( face >= 0 && face < 6 );
        return twistInSelf[face];
      }
    };
    
    template <class IndexSetType,
              class EntityType>
    static inline int calculateSimplexTwistInNeighbor
                                    (const IndexSetType& set, 
                                     const EntityType& en,
                                     const int inSelf,
                                     const EntityType& nb,
                                     const int inNeigh) 
    {
      typedef typename EntityType :: ctype ctype;
      const ReferenceElement< ctype, dim > & enRef = 
        ReferenceElements< ctype, dim >::general(en.geometry().type());

      const ReferenceElement< ctype, dim > & nbRef = 
        ReferenceElements< ctype, dim >::general(nb.geometry().type());
        
      // number of vertices of face 
      const int numVert = enRef.size (inSelf,1, dim);
      int enVx [4];
      int nbVx [4];
      
      int faceMap[4] = { 0, 1, 2, 3};

      bool allRight = true;
      for(int i=0; i<numVert; ++i ) 
      {
        enVx[i] = set.template subIndex<dim>(en, enRef.subEntity(inSelf,1,i,dim)); 
        nbVx[i] = set.template subIndex<dim>(nb, nbRef.subEntity(inNeigh,1,i,dim)); 
        if( enVx[i] != nbVx[i] ) allRight = false;
      }

      if( !allRight )
      {
        for(int i=0; i<numVert; ++i)
        {
          if(enVx[i] != nbVx[i])
          {
            for(int k=1; k<numVert; ++k)
            {
              int newvx = (i+k) % numVert;
              if( enVx[i] == nbVx[newvx] ) faceMap[i] = newvx;
            }
          }
        }
      }

      // return twist 

      if (faceMap[1] == (faceMap[0]+1) % numVert) 
      {
        return faceMap[0];
      }
      else 
      {
        int twst = faceMap[1] - numVert;
        if( numVert == 3 ) 
        {
          // same bug as in Alberta (check reference elements)
          if( twst == -3 ) return -2;
          else if ( twst == -2 ) return -3; 
          else return twst; 
        }
        else 
          return twst; 
      }
    }

  public:
    //! \brief constructor taking grid reference 
    TwistUtility(const GridType& grid) {}

    //! \brief return 0 for inner face 
    template <class Intersection> 
    static inline int twistInSelf(const GridType &, const Intersection& it)
    {
      if( it.inside()->type().isSimplex() )
      {
        // inside twist for simplices is zero 
        return 0;
      }
      else 
      {
        assert( it.inside()->type().isCube() );
        return CubeTwists :: twistInSelf(it.numberInSelf());
      }
    }    

    //! \brief return 0 for inner face 
    template <class IntersectionIterator> 
    int twistInSelf(const IntersectionIterator& it) const 
    {
      DUNE_THROW(NotImplemented,"not implemented because grid is missing!");
      abort();
      return 0; 
    }
    
    //! \brief return 0 for outer face 
    template <class IntersectionIterator> 
    int twistInNeighbor(const IntersectionIterator& it) const 
    {
      DUNE_THROW(NotImplemented,"not implemented because grid is missing!");
      abort();
      return 0;
    }

    //! \brief return 0 for outer face 
    template <class Intersection> 
    static inline int twistInNeighbor(const GridType & grid, const Intersection& it) 
    {
      assert( it.neighbor () );
      if( it.outside()->type().isSimplex() )
      {
        return calculateSimplexTwistInNeighbor( grid.leafIndexSet(),
                                           * it.inside(),
                                           it.numberInSelf(),
                                           * it. outside(),
                                           it.numberInNeighbor());
      }
      else 
      {
        assert( it.outside()->type().isCube() );
        return CubeTwists :: twistInNeighbor(it.numberInNeighbor());
      }
    }

    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    bool conforming (const IntersectionIterator& it) const 
    { 
      return (it.neighbor()) ? 
        (it.inside()->level() == it.outside()->level()) : true; 
    }
    
    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    static inline bool conforming (const GridType &, const IntersectionIterator& it)
    { 
      return (it.neighbor()) ? 
        (it.inside()->level() == it.outside()->level()) : true; 
    }
    
    /** \brief return geometry type of inside or outside entity */
    template <class Intersection>  
    static inline GeometryType
    elementGeometry(const Intersection& intersection, 
                    const bool inside) 
    {
      return (inside) ? intersection.inside()->type() :  
                        intersection.outside()->type();  
    }
  };
#endif


#if HAVE_DUNE_GEOGRID
  template< class HostGrid, class CoordFunction >
  class GeometryGrid;

  template< class HostGrid, class CoordFunction >
  class TwistUtility< GeometryGrid< HostGrid, CoordFunction > >
  {
  public:
    typedef GeometryGrid< HostGrid, CoordFunction > GridType;
    typedef typename GridType :: Traits :: LeafIntersectionIterator
      LeafIntersectionIterator;
    typedef typename LeafIntersectionIterator :: Intersection LeafIntersection;
    typedef typename GridType :: Traits :: LevelIntersectionIterator
      LevelIntersectionIterator;
    typedef typename LevelIntersectionIterator :: Intersection LevelIntersection;

  private:
    typedef TwistUtility< HostGrid > HostTwistUtility;
    typedef Dune :: HostGridAccess< GridType > HostGridAccess;

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
    static bool conforming ( const GridType &grid,
                             const LeafIntersection &intersection )
    {
      typedef typename HostGridAccess :: HostLeafIntersection HostIntersection;
      const HostIntersection &hostIntersection
        = HostGridAccess :: getIntersection( intersection );
      return HostTwistUtility :: conforming( grid.hostGrid(), hostIntersection );
    }

    //! \brief return true if intersection is conform
    static bool conforming ( const GridType &grid,
                             const LevelIntersection &intersection )
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
    bool conforming ( const Intersection &intersection ) const
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
#endif
  
} // end namespace Dune 

#endif
