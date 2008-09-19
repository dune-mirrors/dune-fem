#ifndef DUNE_GRIDWIDTH_HH
#define DUNE_GRIDWIDTH_HH

//- Dune includes 
#include <dune/common/fvector.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/fem/storage/singletonlist.hh>
#include <dune/fem/space/common/allgeomtypes.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/common/dofmanager.hh>

namespace Dune {

template <class GridPartImp>
class GridWidthProvider;

//! utility functions for calculating the maximum grid width 
struct GridWidth 
{
  // for convenience  
  template <class GridPartType, class EntityType> 
  static inline double maxEdgeWidth(const GridPartType& gridPart,
                                    const EntityType &en)
  {
    typedef typename GridPartType :: GridType GridType;

    // get geo infor for elements 
    typedef AllGeomTypes< typename GridPartType :: IndexSetType, GridType > GeomInfoType;
    GeomInfoType geoInfo( gridPart.indexSet() );

    // get geo infor for faces  
    typedef GeometryInformation< GridType , 1 > FaceGeometryInformationType;
    FaceGeometryInformationType faceGeoInfo( geoInfo.geomTypes(1) );

    typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;
    IntersectionIteratorType begin = gridPart.ibegin(en);
    const IntersectionIteratorType end = gridPart.iend(en);

    return maxEdgeWidth(begin, end, gridPart, geoInfo, faceGeoInfo, en);
  }
  
  template <class IntersectionIteratorType, class ElemGeoInfoType, class FaceGeoInfoType, class EntityType> 
  static inline double maxEdgeWidth(IntersectionIteratorType& it, 
                                    const IntersectionIteratorType& endit,
                                    const ElemGeoInfoType& geoInfo,
                                    const FaceGeoInfoType& faceGeoInfo, 
                                    const EntityType &en)
  {
    typedef typename EntityType :: Geometry Geometry;
    enum { dim = EntityType::dimension };

    const Geometry& geo = en.geometry();
    const double elemVol = geo.volume() * geoInfo.referenceVolume( geo.type() );
    
    double faceVol = 1e10;
    int numberInSelf = -1;
    double currVol = -1e10;
    double refFaceVol = -1e10;
    
    for( ; it != endit; ++it)
    {
      typedef typename IntersectionIteratorType :: Intersection IntersectionType;
      const IntersectionType& inter = *it;
      
      typedef typename IntersectionIteratorType :: Geometry LocalGeometryType;
      const LocalGeometryType& interGeo = inter.intersectionGlobal();

      // calculate face Volume also for non-conforming intersections  
      if( numberInSelf != inter.numberInSelf() )
      {
        if (numberInSelf >= 0) 
        {
          faceVol = std::min( faceVol, currVol * refFaceVol );
        }
        
        refFaceVol = faceGeoInfo.referenceVolume( interGeo.type() );
        numberInSelf = inter.numberInSelf();
        currVol = 0.0;
      }

      currVol += interGeo.volume();
    }

    return elemVol/faceVol;
  }

  template <class GridPartType> 
  static inline double calcGridWidth (const GridPartType & gridPart)
  {     
    double maxwidth = 0.0;
    typedef typename GridPartType :: GridType GridType;
    typedef typename GridPartType :: template Codim<0> :: IteratorType IteratorType; 

    // get geo infor for elements 
    typedef AllGeomTypes< typename GridPartType :: IndexSetType, GridType > GeomInfoType;
    GeomInfoType geoInfo( gridPart.indexSet() );

    // get geo infor for faces  
    typedef GeometryInformation< GridType , 1 > FaceGeometryInformationType;
    FaceGeometryInformationType faceGeoInfo( geoInfo.geomTypes(1) );

    // unstructured case 
    if( Capabilities::IsUnstructured<GridType>::v )
    {
      IteratorType endit = gridPart.template end<0> (); 
      for(IteratorType it = gridPart.template begin<0> (); 
          it != endit; ++it )
      {
        const double w = maxEdgeWidth(gridPart, geoInfo, faceGeoInfo, *it);
        if(w > maxwidth) maxwidth = w;
      }
    }
    else 
    {
      // here we only need to check one element 
      IteratorType it = gridPart.template begin<0> (); 
      if( it != gridPart.template end<0> () )
      {
        const double w = maxEdgeWidth(gridPart, geoInfo, faceGeoInfo, *it);
        if(w > maxwidth) maxwidth = w;
      }
    }
    return maxwidth;
  }   

  template <class EntityType, class GeomInfoType, class FaceGeomInfoType> 
  static inline double  
  leafEdgeWidth(const EntityType &en,
                const GeomInfoType& geoInfo,
                const FaceGeomInfoType& faceGeoInfo)
  {
    typedef typename EntityType :: LeafIntersectionIterator IntersectionIteratorType;
    const IntersectionIteratorType endit = en.ileafend();
    IntersectionIteratorType it = en.ileafbegin();

    return maxEdgeWidth(it, endit, geoInfo, faceGeoInfo, en);
  }

  template <class GridType, class GeomInfoType, class FaceGeomInfoType> 
  static inline double calcLeafWidth (
      const GridType & grid, 
      const GeomInfoType& geoInfo,
      const FaceGeomInfoType& faceGeoInfo)
  {     
    double width = 1e308;
    typedef typename GridType::template Codim<0> :: LeafIterator IteratorType; 
    
    // unstructured case 
    if( Capabilities::IsUnstructured<GridType>::v )
    {
      IteratorType endit = grid.template leafend<0> (); 
      for(IteratorType it = grid.template leafbegin<0> (); 
          it != endit; ++it )
      {
        width = std::min( width , leafEdgeWidth(*it, geoInfo, faceGeoInfo) );
      }
    }
    else 
    {
      // here we only need to check one element 
      IteratorType it = grid.template leafbegin<0> (); 
      if( it != grid.template leafend<0> () )
      {
        width = std::min( width , leafEdgeWidth(*it, geoInfo, faceGeoInfo) );
      }
    }
    return width;
  }   
};

//! utility functions for calculating the maximum grid width 
template <class GridType>
class GridWidthProvider 
{
  typedef GridWidthProvider < GridType > ThisType;
public:
  //! type of singleton provider 
  typedef SingletonList< const GridType* , ThisType > ProviderType;

protected:
  typedef DofManager<GridType> DofManagerType;
  typedef DofManagerFactory< DofManagerType > DMFactoryType;

  typedef HierarchicGridPart< GridType > GridPartType; 
  
  // get geo infor for elements 
  typedef AllGeomTypes< typename GridPartType :: IndexSetType, GridType > GeomInfoType;

  // get geo infor for faces  
  typedef GeometryInformation< GridType , 1 > FaceGeometryInformationType;

  const GridType& grid_;
  const DofManagerType& dm_;

  GridPartType gridPart_;

  const GeomInfoType geoInfo_;
  const FaceGeometryInformationType faceGeoInfo_;
  
  mutable double width_;
  mutable int sequence_;

  GridWidthProvider( const ThisType& );
public:
  //! constructor taking grid part 
  GridWidthProvider(const GridType* grid) 
    : grid_( *grid )
    , dm_( DMFactoryType::getDofManager( grid_ ))
    , gridPart_( const_cast<GridType& > (grid_) )
    , geoInfo_( gridPart_.indexSet() )
    , faceGeoInfo_( geoInfo_.geomTypes(1) ) 
    , width_(-1.0)
    , sequence_(-1)
  {
  }

  //! return characteristic grid width 
  double gridWidth () const 
  {
    calcWidths();
    return width_;
  }
  
  //! return characteristic grid width 
  double minGridWidth () const DUNE_DEPRECATED
  {
    return gridWidth();
  }
  
  //! return characteristic grid width 
  double maxGridWidth () const DUNE_DEPRECATED
  {
    return gridWidth();
  }

protected:  
  void calcWidths() const 
  {
    if( dm_.sequence() != sequence_ )
    {
      width_ = GridWidth::calcLeafWidth( grid_ , geoInfo_, faceGeoInfo_ );
      assert( width_ > 0 );
      sequence_  = dm_.sequence();
    }
  }
};

} // end namespace Dune 
#endif

