#ifndef DUNE_GRIDWIDTH_HH
#define DUNE_GRIDWIDTH_HH

//- Dune includes 
#include <dune/common/fvector.hh>
#include <dune/grid/common/capabilities.hh>

namespace Dune {

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

};

} // end namespace Dune 
#endif

