#ifndef DUNE_GRIDWIDTH_HH
#define DUNE_GRIDWIDTH_HH

//- Dune includes 
#include <dune/common/fvector.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/fem/space/common/singletonlist.hh>
#include <dune/fem/space/common/dofmanager.hh>

namespace Dune {

template <class GridPartImp>
class GridWidthProvider;

//! utility functions for calculating the maximum grid width 
struct GridWidth 
{
  template <class GridPartType, class EntityType> 
  static inline double maxEdgeWidth(const GridPartType& gridPart, const EntityType &en)
  {
    typedef typename GridPartType::IntersectionIteratorType
      IntersectionIteratorType;
    enum { dim = EntityType::dimension };

    Dune::FieldVector<double, dim-1> tmp1(0.5);
    
    double vol = en.geometry().volume();
    double fak = 2.0;
    if(dim == 3) fak = 1.0;

    double minFace = 1e308;
    IntersectionIteratorType endit = gridPart.iend(en);
    for(IntersectionIteratorType it = gridPart.ibegin(en);
        it != endit; ++it)
    {
      double face = it.integrationOuterNormal(tmp1).two_norm();
      minFace = std::min(minFace,face);
    }
    return fak*vol/minFace;
  }

  template <class GridPartType> 
  static inline double calcGridWidth (const GridPartType & gridPart)
  {     
    double maxwidth = 0.0;
    typedef typename GridPartType::template Codim<0> :: IteratorType IteratorType; 
    
    // unstructured case 
    if( Capabilities::IsUnstructured<GridType>::v )
    {
      IteratorType endit = gridPart.template end<0> (); 
      for(IteratorType it = gridPart.template begin<0> (); 
          it != endit; ++it )
      {
        double w = maxEdgeWidth(gridPart,*it);
        if(w > maxwidth) maxwidth = w;
      }
    }
    else 
    {
      // here we only need to check one element 
      IteratorType it = gridPart.template begin<0> (); 
      if( it != gridPart.template end<0> () )
      {
        double w = maxEdgeWidth(gridPart,*it);
        if(w > maxwidth) maxwidth = w;
      }
    }
    return maxwidth;
  }   

  template <class EntityType> 
  static inline double maxLeafEdgeWidth(const EntityType &en)
  {
    typedef typename EntityType :: LeafIntersectionIterator IntersectionIteratorType;
    enum { dim = EntityType::dimension };

    Dune::FieldVector<double, dim-1> tmp1(0.5);
    
    double vol = en.geometry().volume();
    double fak = 2.0;
    if(dim == 3) fak = 1.0;

    double minFace = 1e308;
    IntersectionIteratorType endit = en.ileafend();
    for(IntersectionIteratorType it = en.ileafbegin();
        it != endit; ++it)
    {
      double face = it.integrationOuterNormal(tmp1).two_norm();
      minFace = std::min(minFace,face);
    }
    return fak*vol/minFace;
  }

  template <class GridType> 
  static inline double calcLeafWidth (const GridType & grid)
  {     
    double maxwidth = 0.0;
    typedef typename GridType::template Codim<0> :: LeafIterator IteratorType; 
    
    // unstructured case 
    if( Capabilities::IsUnstructured<GridType>::v )
    {
      IteratorType endit = grid.template leafend<0> (); 
      for(IteratorType it = grid.template leafbegin<0> (); 
          it != endit; ++it )
      {
        double w = maxLeafEdgeWidth(*it);
        if(w > maxwidth) maxwidth = w;
      }
    }
    else 
    {
      // here we only need to check one element 
      IteratorType it = grid.template leafbegin<0> (); 
      if( it != grid.template leafend<0> () )
      {
        double w = maxLeafEdgeWidth(*it);
        if(w > maxwidth) maxwidth = w;
      }
    }
    return maxwidth;
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

private:
  typedef DofManager<GridType> DofManagerType;
  typedef DofManagerFactory< DofManagerType > DMFactoryType;

  const GridType& grid_;
  const DofManagerType& dm_;
  mutable double gridWidth_;
  mutable int sequence_;

  GridWidthProvider( const ThisType& );
public:
  //! constructor taking grid part 
  GridWidthProvider(const GridType& grid) 
    : grid_( grid )
    , dm_( DMFactoryType::getDofManager( grid_ ))
    , gridWidth_(-1.0)
    , sequence_(-1)
  {}

  //! return characteristic grid width 
  double gridWidth () const 
  {
    if( dm_.sequence() != sequence_ )
    {
      gridWidth_ = GridWidth::calcLeafWidth( grid_ );
      sequence_  = dm_.sequence();
    }
    return gridWidth_;
  }
};

} // end namespace Dune 
#endif

