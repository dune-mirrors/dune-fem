#ifndef DUNE_FILTEREDGRID_HH
#define DUNE_FILTEREDGRID_HH

//- System includes
#include <vector>
#include <cassert>

//- Dune includes
#include <dune/common/interfaces.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/referenceelements.hh>
#include <dune/grid/common/defaultindexsets.hh>
#include <dune/grid/common/sizecache.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/grid/utility/grapedataioformattypes.hh>

//***************************************************************************
//! A FilteredGridPart is a subset of a GridPart and a GridPart itself 
//! (without iterators for codim \neq 0). FilteredGridPart will work with
//! LeafGridPart and HierarchicGridPart but not with a LevelGridPart. If 
//! you wish to use FilteredGridPart together with a LevelGridPart you are
//! asked to write a specialization of FilteredGridPart, which should be 
//! very easy - imO only the constructor has to be changed (see below).
//! The codim 0 entities that belong to the FilteredGrid are defined by a 
//! filter class. 
//! Again: Be careful, yet we have only iterators for codim 0 entities on the 
//! FilteredGridPart!
//
// @author Christoph Gersbacher
//***************************************************************************

namespace Dune {

  // forward declarations
  template <class FilterImp, class GridImp>
  struct DefaultFilterTraits;
  template <class FilterTraits>
  class FilterInterface;
  //template <class GridType>
  //class TrueFilter;
  template <class GridType>
  class RadialFilter;
  template <class GridPartImp, class FilterImp, PartitionIteratorType pitype>
  class FilteredGridPart;

//***************************************************************************
// 
// DefaultFilterTraits
//
//***************************************************************************

  //! Type definitions 
  template <class FilterImp, class GridImp>
  struct DefaultFilterTraits {
    typedef FilterImp FilterType;
    typedef GridImp GridType;
    typedef typename GridType::template Codim<0>::Entity EntityCodim0Type;   
    typedef typename GridType::template Codim<0>::EntityPointer EntityPointerCodim0Type;

    template <int cd>
    struct Codim {
      typedef typename GridType::template Codim<cd>::Entity EntityType;
      typedef typename GridType::template Codim<cd>::EntityPointer EntityPointerType;
    };
  }; // end DefaultFilterTraits

//***************************************************************************
// 
// FilteredInterface
//
//***************************************************************************
  
  template <class FilterTraits>
  class FilterInterface
  {
  public:
    //! \brief Type of the filter implementation
    typedef typename FilterTraits::FilterType FilterType;  
    //! \brief type of Grid implementation
    typedef typename FilterTraits::GridType GridType;
    //! \brief type of Entity with codim=0 
    typedef typename FilterTraits::EntityCodim0Type EntityCodim0Type; 
    //! \brief type of EntityPointer with codim=0 
    typedef typename FilterTraits::EntityPointerCodim0Type EntityPointerCodim0Type;
      
    inline int has0Entity(EntityCodim0Type & e) const {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().has0Entitiy(e)));
      return asImp().has0Entity(e);
    }
    inline int has0Entity(EntityPointerCodim0Type & e) const {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().has0Entitiy(e)));
      return asImp().has0Entity(e);
    }
    inline bool boundary() const {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().boundary()));
      return asImp().boundary();
    }
    inline int boundaryId() const {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().boundaryId()));
      return asImp().boundaryId();
    }
     inline bool neighbor() const {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().neighbor()));
      return asImp().neighbor();
    }
  protected: 
    //! do not create explict instances of this class 
    FilterInterface () {}  

  private:
    // Barton-Nackman 
    FilterType& asImp() { 
      return static_cast<GridType&>(*this); 
    }
    
    // const Barton-Nackman 
    const FilterType& asImp() const { 
      return static_cast<const GridType&>(*this);
    }  
  
  }; // end FilteredInterface
  
//***************************************************************************
// 
// RadialFilter
//
//***************************************************************************

  template <class GridType>
  class RadialFilter : 
    public FilterInterface<DefaultFilterTraits<RadialFilter<GridType>, GridType> >
  {

  public:
    typedef DefaultFilterTraits<RadialFilter<GridType>, GridType> Traits;
    typedef FilterInterface<DefaultFilterTraits<RadialFilter<GridType>,GridType> > BaseType;
    typedef typename BaseType::FilterType FilterType;
    typedef typename BaseType::EntityCodim0Type EntityCodim0Type;
    typedef typename BaseType::EntityPointerCodim0Type EntityPointerCodim0Type;

  private:
    typedef typename GridType::Traits::template Codim<0>::Geometry GeometryType;
    enum{dim = GeometryType::dimension};
    typedef typename Dune::FieldVector<typename GridType::ctype, dim> FieldVectorType;

  public:
    RadialFilter(FieldVectorType center, double radius) : center_(center), radius_(radius){ }

    inline int has0Entity(EntityPointerCodim0Type & e) const {
      return has0Entity(*e);
    }
    
    // a very expensive computation...
    inline int has0Entity(EntityCodim0Type & e) const {
      const GeometryType& geom = e.geometry();
      FieldVectorType weigh(0);
      for(int i = 0; i < geom.corners(); ++i)
        weigh += geom[i];
      weigh /= geom.corners();
      weigh -= center_;
      double dist = weigh.two_norm();
      return (dist >= radius_) ? 0 : 1;
    }

    inline bool boundary() const {
      return true;
    }
    inline int boundaryId() const {
      return 1;
    }
    inline bool neighbor() const {
      return true;
    }

  private:
    FieldVectorType center_;
    double radius_;
  }; // end RadialFilter

//***************************************************************************
// 
// FilteredGridPart
//
//***************************************************************************

  template <class GridPartImp, class FilterImp, PartitionIteratorType pitype = Interior_Partition>
  class FilteredGridPart :
    public GridPartImp
  {

  private:
    // forward declaration of IteratorWrappers
    template <class GridPartType, int cd, class IteratorType>
    struct IteratorWrapper;
    template <class GridPartType, class IteratorType>
    struct IntersectionIteratorWrapper;
    // type of this
    typedef FilteredGridPart<GridPartImp, FilterImp, pitype> ThisType;
    // type of underlying GridPart
    typedef GridPartImp GridPartType;
    // type of filter
    typedef FilterImp FilterType;
    // the codim 0 entities type
    typedef typename GridType::template Codim<0>::Entity EntityCodim0Type;
    // the original IntersectionIteratorType
    typedef typename GridPartImp::IntersectionIteratorType IntersectionIteratorImpType;

  public:
    //- Public typedefs and enums    
    //! Grid implementation type
    typedef typename GridPartImp::GridType GridType;
    //! The leaf index set of the grid implementation
    typedef typename GridPartImp::IndexSetType IndexSetType;
    //! The corresponding IntersectionIterator 
    typedef IntersectionIteratorWrapper<GridPartType, IntersectionIteratorImpType> IntersectionIteratorType;
    
    //! Struct providing types of the iterators on codimension cd
    template <int cd>
    struct Codim {
    private:
      typedef typename GridPartImp::template Codim<cd>::IteratorType IteratorImpType;
    public:
      typedef IteratorWrapper<GridPartType, cd, IteratorImpType> IteratorType;
    };

  public:
    //- Public methods
    //! Constructor
    FilteredGridPart(const GridType& grid, FilterType &filter) :
      GridPartImp(grid),   // This is only legal for LeafGridPart and HierarchicGridParts!
                           // A constructor LevelGridPart(GridType grid) exists, but will  
                           // return a LevelGridPart<GridType, maxLevel>.
      filter_(filter),
      maxlevel_(0)
    {
      updateStatus();
    }

    //! Begin iterator on the leaf level
    template <int cd>
    typename ThisType::template Codim<cd>::IteratorType begin() const {
      //typedef typename ThisType::template Codim<cd>::IteratorType IteratorType;
      return typename ThisType::template Codim<cd>::IteratorType(this, filter_, false);
    }

    //! End iterator on the leaf level
    template <int cd>
    typename ThisType::template Codim<cd>::IteratorType end() const {
      //typedef typename ThisType::template Codim<cd>::IteratorType IteratorType;
      return typename ThisType::template Codim<cd>::IteratorType(this, filter_, true);
    }

    //! ibegin of corresponding intersection iterator for given entity
    IntersectionIteratorType ibegin(const EntityCodim0Type & en) const 
    {
      return typename ThisType::IntersectionIteratorType(this, filter_, en, false);
    }
    
    //! iend of corresponding intersection iterator for given entity
    IntersectionIteratorType iend(const EntityCodim0Type & en) const 
    {
      return typename ThisType::IntersectionIteratorType(this, filter_, en, true);
    }

    //! Returns maxlevel of the grid
    int level() const { return maxlevel_; }

    //! corresponding communication method for this grid part
    template <class DataHandleImp,class DataType>
    void communicate(CommDataHandleIF<DataHandleImp,DataType> & data, 
                     InterfaceType iftype, CommunicationDirection dir) const 
    {
      this->grid().communicate(data,iftype,dir);
    }

  private:   
    inline void updateStatus()
    {
      calcMaxlevel();
    }
    inline void calcMaxlevel()
    {
      maxlevel_ = GridPartImp::level();
    }

  protected:
    FilterType & filter() const { return filter_; }

  private: 
    FilterType & filter_;
    int maxlevel_;

  private:
    //**********************************************************************
    // IteratorWrapper
    //**********************************************************************

    template <class GridPartType, int cd, class IteratorType>
    struct IteratorWrapper : public IteratorType {
     
      public:
      inline IteratorWrapper(const GridPartType * gridPart, FilterType & filter,  bool endIter):
        IteratorType(endIter?gridPart->template end<cd>():gridPart->template begin<cd>()),
        gridPart_(gridPart),
        filter_(filter),          
      	endIter_(gridPart->template end<cd>())
        { 
	        if (!endIter)
	          if(!filter_.has0Entity(*this))
	            operator++();           
        }

        inline IteratorWrapper & operator++()
        {	  
          do {
            IteratorType::operator++();
      	    if (*this==endIter_) break;
          } while(!filter_.has0Entity(*this));
          return *this;
        }
        
      private:
        const GridPartType * gridPart_;        
        FilterType & filter_;
        const IteratorType endIter_;
    }; // end IteratorWrapper

    //**********************************************************************
    // IntersectionIteratorWrapper
    //**********************************************************************

    template <class GridPartType, class IteratorType>
    struct IntersectionIteratorWrapper : public IteratorType {

      // type of codim 0 entity      
      // typedef typename GridPartType::EntityCodim0Type EntityCodim0Type;
      typedef typename GridPartType::GridType::template Codim<0>::Entity EntityCodim0Type;

      struct neighborInfo {    
        neighborInfo() : boundary_(false), boundaryId_(-1), neighbor_(false) { }
        neighborInfo(const neighborInfo & org) : boundary_(org.boundary_), 
          boundaryId_(org.boundaryId_), neighbor_(org.neighbor_){ }
        neighborInfo & operator = (const neighborInfo & org) {
          boundary_ = org.boundary_;
          boundaryId_ = org.boundaryId_; 
          neighbor_ = org.neighbor_;        
          return *this;
        }        
        bool boundary_;
        int boundaryId_; 
        bool neighbor_;        
      } nInfo;

      public:
      inline IntersectionIteratorWrapper(const GridPartType * gridPart, FilterType & filter, 
                                         EntityCodim0Type &en, bool endIter):
        IteratorType(endIter?gridPart->template iend(en):gridPart->template ibegin(en)),
        gridPart_(gridPart),
        filter_(filter),          
       	endIter_(gridPart->template iend(en))
        { 
	        //if (!endIter)
	        //  if(!filter_.has0Entity(*this))
       	  //    operator++();           
          //return *this;
        }
       
        inline IntersectionIteratorWrapper & operator++()
        {         
          if (*this!=endIter_)
            IteratorType::operator++();
          if (*this!=endIter_) {    
            if (IteratorType::neighbor()) {
              EntityCodim0Type& neigh = *(this->outside());
              if (filter_.has0Entity(neigh)) {
                nInfo.boundary = false;
	        nInfo.boundaryId = 0;
       	        nInfo.neighbor_ = true;	  
	      }
              else{
                nInfo.boundary = filter_.intersectionBoundary();
       	        nInfo.boundaryId = filter_.intersectionBoundaryId();
	        nInfo.neighbor_ = filter_.intersectionNeighbor();	  
              }
            }
            else {
              nInfo.boundary = filter_.intersectionBoundary();
     	      nInfo.boundaryId = filter_.intersectionBoundaryId();
	      nInfo.neighbor_ = false;	  
	    }
      	  }          
          return *this;
        }
          
        inline bool boundary() { return nInfo.boundary_; }
        inline int boundaryId() { return nInfo.boundaryId_; }
        inline bool neighbor(){ return nInfo.neighbor_; }

      private:
        const GridPartType * gridPart_;        
        FilterType & filter_;
        const IteratorType endIter_;        
    }; // end IntersectionIteratorWrapper
   
  }; // end FilteredGridPart

}  // end namespace Dune

#endif
