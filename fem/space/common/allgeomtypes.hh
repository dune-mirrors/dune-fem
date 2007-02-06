#ifndef DUNE_ALLGEOMTYPES_HH
#define DUNE_ALLGEOMTYPES_HH

//- system includes 
#include <vector>

//- Dune includes 
#include <dune/common/geometrytype.hh>
#if HAVE_UG 
#include <dune/grid/uggrid.hh>
#endif

namespace Dune {

  /** \brief default implementation uses method geomTypes of given index
       set. 
   */
  template <class IndexSetImp, class GridImp> 
  class AllGeomTypes 
  {
    typedef IndexSetImp IndexSetType;
    const IndexSetType & indexSet_;
  public:
    //! constructor storing index set reference 
    AllGeomTypes(const IndexSetType & indexSet) 
      : indexSet_(indexSet) 
    {}                  
    
    //! returns vector with geometry tpyes this index set has indices for
    const std::vector <GeometryType> & geomTypes (int codim) const
    {
      return indexSet_.geomTypes(codim);
    }

    //! UGGrid might have different geom types 
    static bool multipleGeomTypes () { return false; }
  };

#if HAVE_UG 
  /** \brief specialisation fir UGGrid, because geomTypes method of index
      sets not usable in this case. 
  */
  template <class IndexSetImp, int dimworld> 
  class AllGeomTypes<IndexSetImp, UGGrid<dimworld> >
  {
    typedef IndexSetImp IndexSetType;
    enum { ncodim = dimworld + 1 };
    
    std::vector< GeometryType > geomTypes_[ncodim]; 
  public:
    //! constructor storing index set reference 
    AllGeomTypes(const IndexSetType & indexSet) 
    {
      // UGGrid only working for dimension 2 or 3 
      assert( (dimworld == 2) || (dimworld == 3) ); 
      // vertices 
      {
        geomTypes_[dimworld].push_back( GeometryType(GeometryType::simplex,0)); 
      }
      
      {
        if( dimworld == 2 )
        {
          // elements 
          geomTypes_[0].push_back( GeometryType(GeometryType::simplex,2)); 
          geomTypes_[0].push_back( GeometryType(GeometryType::cube,2)); 

          // faces 
          geomTypes_[1].push_back( GeometryType(GeometryType::cube,1));
        }

        if( dimworld == 3 )
        {
          // elements 
          geomTypes_[0].push_back( GeometryType(GeometryType::simplex,3)); 
          geomTypes_[0].push_back( GeometryType(GeometryType::cube,3)); 
          geomTypes_[0].push_back( GeometryType(GeometryType::prism,3)); 
          geomTypes_[0].push_back( GeometryType(GeometryType::pyramid,3)); 

          // faces 
          geomTypes_[1].push_back( GeometryType(GeometryType::simplex,2)); 
          geomTypes_[1].push_back( GeometryType(GeometryType::cube,2)); 

          // edges 
          geomTypes_[2].push_back( GeometryType(GeometryType::cube,1));
        }
      }
    }                  
    
    //! returns vector with geometry tpyes this index set has indices for
    const std::vector <GeometryType> & geomTypes (int codim) const
    {
      assert( codim >= 0 );
      assert( codim < ncodim );
      return geomTypes_[codim];
    }

    //! UGGrid might have different geom types 
    static bool multipleGeomTypes () { return true; }
  };
#endif

} // end namespace Dune 
#endif
