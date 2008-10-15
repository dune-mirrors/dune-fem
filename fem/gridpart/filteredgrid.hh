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

#include <dune/fem/gridpart/gridpart.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/grid/utility/grapedataioformattypes.hh>

namespace Dune
{

  // forward declarations
  template <class FilterImp, class GridImp>
  struct DefaultFilterTraits;
  template <class FilterTraits>
  class FilterInterface;
  //template <class GridType>
  //class TrueFilter;
  template <class GridType>
  class RadialFilter;

  template< class GridPartImp, class FilterImp >
  class FilteredGridPart;

//***************************************************************************
// 
// DefaultFilterTraits
//
//***************************************************************************

  //! Type definitions 
  template <class FilterImp, class GridPartImp>
  struct DefaultFilterTraits 
  {
    typedef FilterImp FilterType;
    typedef GridPartImp GridPartType; 
    typedef typename GridPartType :: GridType GridType;
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
/** @ingroup FilterGridPart 
 @brief
 Interface class for filter to use with a Dune::FilteredGridPart
 **/
  template <class FilterTraits>
  class FilterInterface
  {
  public:
    //! \brief Type of the filter implementation
    typedef typename FilterTraits::FilterType FilterType;  
    //! \brief type of Grid implementation
    typedef typename FilterTraits::GridType GridType;
    //! \brief type of original grid part 
    typedef typename FilterTraits::GridPartType GridPartType;
    //! \brief type of Entity with codim=0 
    typedef typename FilterTraits::EntityCodim0Type EntityCodim0Type; 
    //! \brief type of EntityPointer with codim=0 
    typedef typename FilterTraits::EntityPointerCodim0Type EntityPointerCodim0Type;
      
    //! returns true if the given entity of the pointer in the domain 
    inline bool has0Entity(EntityCodim0Type & e) const {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().has0Entity(e)));
      return asImp().has0Entity(e);
    }
    
    //! returns true if the given entity is in the domain 
    inline bool has0Entity(EntityPointerCodim0Type & e) const {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().has0Entity(e)));
      return asImp().has0Entity(e);
    }
    
    //! returns true if an intersection is interior 
    //! (allows boundarys within a given domain)
    template<class IntersectionIteratorType>
    inline bool interiorIntersection(const IntersectionIteratorType & it) const {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().interiorIntersection(it)));
      return asImp().interiorIntersection(it);
    }

    //! returns true if an intersection is a boundary intersection 
    template<class IntersectionIteratorType>
    inline bool intersectionBoundary(IntersectionIteratorType & it) const {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().intersectionBoundary<IntersectionIteratorType>(it)));
      return asImp().intersectionBoundary<IntersectionIteratorType>(it);
    }
    
    //! returns the boundary id for an intersection 
    template<class IntersectionIteratorType>
    inline int intersectionBoundaryId(IntersectionIteratorType & it) const {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().intersectionBoundaryId<IntersectionIteratorType>(it)));
      return asImp().intersectionBoundaryId<IntersectionIteratorType>(it);
    }

    //! returns true if for an intersection a neighbor exsits 
    template<class IntersectionIteratorType>
    inline bool intersectionNeighbor(IntersectionIteratorType & it) const {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().intersectionNeighbor<IntersectionIteratorType>(it)));
      return asImp().intersectionNeighbor<IntersectionIteratorType>(it);
    }

    //! return object instance of filter for object creation  
    inline static FilterType createObject(const GridPartType& gridPart)
    {
      CHECK_INTERFACE_IMPLEMENTATION( FilterType :: createObject(gridPart));
      return FilterType :: createObject(gridPart);
    }

  protected: 
    //! do not create explict instances of this class 
    FilterInterface () {}  

  protected:
    // Barton-Nackman 
    FilterType& asImp() { return static_cast<FilterType&>(*this); }
    
    // const Barton-Nackman 
    const FilterType& asImp() const { return static_cast<const FilterType&>(*this); }  
  }; // end FilteredInterface

  template <class FilterTraits>
  class FilterDefaultImplementation : public FilterInterface<FilterTraits>
  {
    typedef FilterInterface<FilterTraits> BaseType;
    using BaseType :: asImp;
    using BaseType :: has0Entity;
  public:
    //! \brief Type of the filter implementation
    typedef typename FilterTraits::FilterType FilterType;  
    //! \brief type of Grid implementation
    typedef typename FilterTraits::GridType GridType;
    //! \brief type of Entity with codim=0 
    typedef typename FilterTraits::EntityCodim0Type EntityCodim0Type; 
    //! \brief type of EntityPointer with codim=0 
    typedef typename FilterTraits::EntityPointerCodim0Type EntityPointerCodim0Type;
    //! \brief type of original grid part 
    typedef typename FilterTraits::GridPartType GridPartType;
      
    //! default implementation returns hasEntity0 from neighbor
    template<class IntersectionIteratorType>
    inline bool interiorIntersection(IntersectionIteratorType & it) const 
    {
      return asImp().has0Entity( it.outside() );
    }

    //! \brief default createObject method calling FilterType(gridPart) 
    inline static FilterType createObject(const GridPartType& gridPart)
    {
      return FilterType(gridPart);      
    }
  };
  
//***************************************************************************
// 
// Example: RadialFilter
//
//***************************************************************************

  /**! \brief RadialFilter is an example filter class cutting out a
       circle with given center and radius 
  */
  template <class GridPartType>
  class RadialFilter : 
    public FilterDefaultImplementation<
      DefaultFilterTraits<RadialFilter<GridPartType>, GridPartType> >
  {
  public:
    typedef typename GridPartType :: GridType GridType;
    typedef DefaultFilterTraits<RadialFilter<GridPartType>, GridType> Traits;
    typedef FilterDefaultImplementation<
      DefaultFilterTraits<RadialFilter<GridPartType>,GridPartType> > BaseType;
    typedef typename BaseType::FilterType FilterType;
    typedef typename BaseType::EntityCodim0Type EntityCodim0Type;
    typedef typename BaseType::EntityPointerCodim0Type EntityPointerCodim0Type;

  private:
    typedef RadialFilter<GridPartType> ThisType;

    typedef typename GridType::Traits::template Codim<0>::Geometry GeometryType;
    enum{dim = GeometryType::dimension};
    typedef typename Dune::FieldVector<typename GridType::ctype, dim> FieldVectorType;

  public:
    //! constructor defining center and radius of circle to cut out 
    RadialFilter(const FieldVectorType& center, 
                 const double radius) 
      : center_(center), radius_(radius)
    {}

    //! copy constructor 
    RadialFilter(const ThisType& org) 
      : center_(org.center_), radius_(org.radius_)
    {}

    //! check whether bary center is inside of circle 
    inline bool has0Entity(EntityPointerCodim0Type & e) const {
      return has0Entity(*e);
    }
    
    //! check whether bary center is inside of circle 
    inline bool has0Entity(EntityCodim0Type & e) const 
    {
      const GeometryType& geom = e.geometry();
      FieldVectorType weigh(0);
      for(int i = 0; i < geom.corners(); ++i)
        weigh += geom[i];
      weigh /= geom.corners();
      weigh -= center_;
      double dist = weigh.two_norm();
      return (dist <= radius_);
    }
    
    //! return what boundary id we have in case of boundary intersection 
    //! which is either it.boundary == true or has0Entity (it.ouside()) == false 
    //! so here true is a good choice 
    template <class IntersectionIteratorType>
    inline bool intersectionBoundary(IntersectionIteratorType & it) const 
    {
      return true;
    }
    //! return what boundary id we have in case of boundary intersection 
    //! which is either it.boundary == true or has0Entity (it.ouside()) == false 
    template <class IntersectionIteratorType>
    inline int intersectionBoundaryId(IntersectionIteratorType & it) const {
      return 1;
    }

    //! if has0Entity is true then we have an interior entity 
    template <class IntersectionIteratorType>
    inline bool intersectionNeighbor(IntersectionIteratorType & it) const {
      return true;
    }

    static ThisType createObject(const GridPartType& gridPart) 
    {
      DUNE_THROW(NotImplemented,"Method createObject not implemented");
      FieldVectorType center(0);
      double radius = 0;
      return ThisType(center, radius);
    }

  private:
    const FieldVectorType center_;
    const double radius_;
  }; // end RadialFilter

//***************************************************************************
// 
// FilteredGridPart
//
/** @addtogroup FilterGridPart 
 
 A FilteredGridPart is a subset of a GridPart and a GridPart itself 
 (without iterators for codim \f$\neq\f$ 0). FilteredGridPart will work with
 LeafGridPart and HierarchicGridPart but not with a LevelGridPart. If 
 you wish to use FilteredGridPart together with a LevelGridPart you are
 asked to write a specialization of FilteredGridPart, which should be 
 very easy - imO only the constructor has to be changed (see below).
 The codim 0 entities that belong to the FilteredGrid are defined by a 
 filter class. 
 On a codim 0 entitiy there is a method 
   hasBoundaryIntersection().
 This method will not work correctly since the entity is not wrapped. 
 Again: Be careful, yet we have only iterators for codim 0 entities on the 
 FilteredGridPart!
**/

/** @ingroup FilterGridPart
 @brief
 A FilteredGridPart allows to extract a set of entities from a grid
 satisfying a given constrainted defined through a filter class.
**/ 
  template< class GridPartImp, class FilterImp >
  class FilteredGridPart
  : public GridPartImp
  {
    typedef FilteredGridPart< GridPartImp, FilterImp > ThisType;

    // forward declaration of IteratorWrappers
    template< class GridPart, int codim, class Iterator >
    struct IteratorWrapper;
    template< class GridPart, class Iterator >
    struct IntersectionIteratorWrapper;

    // the original IntersectionIteratorType
    typedef typename GridPartImp::IntersectionIteratorType IntersectionIteratorImpType;

  public:
    //! type of filter
    typedef FilterImp FilterType;

    //- Public typedefs and enums    
    //! Grid implementation type
    typedef typename GridPartImp::GridType GridType;

    //! type of underlying grid part 
    //typedef FilteredGridPart<GridPartImp, FilterImp > GridPartType;

    //! The index set of the gridpart implementation
    typedef typename GridPartImp::IndexSetType IndexSetType;
    
    //! The corresponding IntersectionIterator 
    typedef IntersectionIteratorWrapper<GridPartImp, IntersectionIteratorImpType> IntersectionIteratorType;
    // the codim 0 entities type
    typedef typename GridType::template Codim<0>::Entity EntityCodim0Type;

    //! Struct providing types of the iterators on codimension cd
    template< int codim >
    struct Codim
    {
      template< PartitionIteratorType pitype >
      struct Partition
      {
      private:
        typedef typename GridPartImp :: template Codim< codim >
          :: template Partition< pitype > :: IteratorType IteratorImpType;

      public:
        typedef IteratorWrapper< GridPartImp, codim, IteratorImpType > IteratorType;
      };

      typedef typename Partition< InteriorBorder_Partition > :: IteratorType
        IteratorType;
    };

  public:
    //- Public methods
    //! Constructor
    FilteredGridPart(GridType& grid, const FilterType& filter) :
      GridPartImp(grid),   // This is only legal for LeafGridPart and HierarchicGridParts!
                           // A constructor LevelGridPart(GridType grid) exists, but will  
                           // return a LevelGridPart<GridType, maxLevel>.
      filter_(filter),
      maxlevel_(0)
    {
      updateStatus();
    }

    //! Copy Constructor
    FilteredGridPart(const FilteredGridPart& other) :
      GridPartImp(other), 
      filter_(other.filter_),
      maxlevel_(other.maxlevel_)
    {
    }

    //! constructor only taking grid 
    FilteredGridPart(GridType& grid) :
      GridPartImp(grid),   // This is only legal for LeafGridPart and HierarchicGridParts!
                           // A constructor LevelGridPart(GridType grid) exists, but will  
                           // return a LevelGridPart<GridType, maxLevel>.
      filter_( FilterType :: createObject(static_cast<const GridPartImp&> (*this) )), 
      maxlevel_(0)
    {
      updateStatus();
    }

    //! Begin iterator on the leaf level
    template< int codim >
    typename Codim< codim > :: IteratorType
    begin () const
    {
      return begin< codim, InteriorBorder_Partition >();
    }

    //! Begin iterator on the leaf level
    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim > :: template Partition< pitype > :: IteratorType
    begin () const
    {
      typedef typename Codim< codim > :: template Partition< pitype > :: IteratorType IteratorType;
      return IteratorType( this, &filter_, GridPartImp :: template begin< codim, pitype >() );
    }

    //! Begin iterator on the leaf level
    template< int codim >
    typename Codim< codim > :: IteratorType
    end () const
    {
      return end< codim, InteriorBorder_Partition >();
    }

    //! End iterator on the leaf level
    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim > :: template Partition< pitype > :: IteratorType
    end () const
    {
      typedef typename Codim< codim > :: template Partition< pitype > :: IteratorType IteratorType;
      return IteratorType( this, &filter_, GridPartImp :: template end< codim, pitype >() );
    }

    //! ibegin of corresponding intersection iterator for given entity
    inline IntersectionIteratorType ibegin(const EntityCodim0Type & en) const 
    {
      return typename ThisType::IntersectionIteratorType(this, &filter_, en );
    }
    
    //! iend of corresponding intersection iterator for given entity
    inline IntersectionIteratorType iend(const EntityCodim0Type & en) const 
    {
      return typename ThisType::IntersectionIteratorType(this, &filter_, GridPartImp::iend(en) );
    }

    //! Returns maxlevel of the grid
    int level() const 
    { 
       DUNE_THROW(Dune::NotImplemented,"Method FilteredGridPart::level() not implemented");
      return maxlevel_; 
    }

    //! corresponding communication method for this grid part
    template <class DataHandleImp,class DataType>
    void communicate(CommDataHandleIF<DataHandleImp,DataType> & data, 
                     InterfaceType iftype, CommunicationDirection dir) const 
    {
      this->grid().communicate(data,iftype,dir);
    }

    //! return reference to filter 
    const FilterType & filter() const { return filter_; }
    
  private:   
    inline void updateStatus()
    {
      calcMaxlevel();
    }
    inline void calcMaxlevel()
    {
      maxlevel_ = GridPartImp::level();
    }

  private: 
    FilterType filter_;
    int maxlevel_;

  private:
    //**********************************************************************
    // IntersectionIteratorWrapper
    //**********************************************************************

    template <class GridPartType, class IteratorType>
    class IntersectionIteratorWrapper : public IteratorType 
    {
      // type of codim 0 entity      
      typedef typename GridPartType::GridType::template Codim<0>::EntityPointer EntityPointerCodim0Type;
      typedef typename GridPartType::GridType::template Codim<0>::Entity EntityCodim0Type;

      typedef IntersectionIteratorWrapper<GridPartType,IteratorType>  ThisType;
      
    protected:
      class neighborInfo 
      {
        public:
        inline neighborInfo() 
          : boundaryId_(-1), boundary_(false), neighbor_(false) 
        {}
        
        inline neighborInfo(const neighborInfo & org) 
          : boundaryId_(org.boundaryId_), 
            boundary_(org.boundary_), 
            neighbor_(org.neighbor_)
        {}
        
        inline neighborInfo & operator = (const neighborInfo & org) 
        {
          boundary_   = org.boundary_;
          boundaryId_ = org.boundaryId_; 
          neighbor_   = org.neighbor_;        
          return *this;
        }        

        int boundaryId_; 
        bool boundary_;
        bool neighbor_;        
      } nInfo_;

     public:
      //! constructor 
      inline IntersectionIteratorWrapper(const GridPartType* gridPart, 
                                         const FilterType* filter, 
                                         const EntityCodim0Type & en)
        : IteratorType(gridPart->ibegin(en)),
          nInfo_(),
          gridPart_(gridPart),
          filter_(filter),          
          endIter_(gridPart->iend(en))
        { 
          assert( *this != endIter_ );
          writeNeighborInfo();
        }
        
      //! constructor creating end iterator 
      inline IntersectionIteratorWrapper(const GridPartType * gridPart, 
                                         const FilterType* filter, 
                                         const IteratorType& endIter)
        : IteratorType(endIter),
          nInfo_(),
          gridPart_(gridPart),
          filter_(filter),          
          endIter_(endIter)
      { 
      }
        
      //! copy constructor 
      inline IntersectionIteratorWrapper(const IntersectionIteratorWrapper& other) 
        : IteratorType(other),
          nInfo_(other.nInfo_),
          gridPart_(other.gridPart_),
          filter_(other.filter_),          
          endIter_(other.endIter_)
      { 
      }
        
      //! assignment operator 
      inline IntersectionIteratorWrapper& operator = (const IntersectionIteratorWrapper& other) 
      {
        IteratorType :: operator = (other);
        nInfo_    = other.nInfo_; 
        gridPart_ = other.gridPart_;
        filter_   = other.filter_;
        endIter_  = other.endIter_;
        return *this;
      }
        
      protected:
        //! write information for current intersection 
        inline void writeNeighborInfo() 
        {
          if ( IteratorType::neighbor() ) 
          { 
            if ( filter_->interiorIntersection( asBase() ) )
            {
              nInfo_.boundary_   = false;
              nInfo_.boundaryId_ = 0;
              nInfo_.neighbor_   = true;
            }
            else 
            {
              // otherwise get boundary information from filter 
              nInfo_.boundary_   = filter_->intersectionBoundary( asBase() );
              nInfo_.boundaryId_ = filter_->intersectionBoundaryId( asBase() );
              nInfo_.neighbor_   = filter_->intersectionNeighbor( asBase() );
            }
          }
          else 
          {
            // for real boundary get boundary from filter 
            nInfo_.boundary_   = true;
            nInfo_.boundaryId_ = filter_->intersectionBoundaryId( asBase() );
            nInfo_.neighbor_   = false;
          }    
        }

      public:
        //! increment intersection iterator 
        inline IntersectionIteratorWrapper & operator++()
        { 
          // if iterator in-valid , od nothing 
          if (*this == endIter_) return *this; 
          
          // increment real iterator 
          IteratorType::operator++();
          if( *this == endIter_ ) return *this; 
            
          // write new infos 
          writeNeighborInfo();

          return *this;
        }
          
        //! overloaded boundary method 
        inline bool boundary() const  { return nInfo_.boundary_; }
        //! overloaded boundaryId method 
        inline int boundaryId() const { return nInfo_.boundaryId_; }
        //! overloaded neighbor method 
        inline bool neighbor() const { return nInfo_.neighbor_; }

        //! type of Intersection 
        typedef ThisType Intersection;
        //! dereference operator 
        inline const Intersection& operator *() const { return *this; }
        //! de-pointer operator 
        inline const Intersection* operator ->() const { return this; }
      protected:
        //! return reference to base class 
        inline IteratorType & asBase() { return static_cast<IteratorType &>(*this); }
        
        const GridPartType* gridPart_;        
        const FilterType* filter_;
        const IteratorType endIter_;        
    }; // end IntersectionIteratorWrapper
  }; // end FilteredGridPart



  // FilteredGridPart :: IteratorWrapper
  // -----------------------------------

  template< class GridPartImp, class FilterImp >
  template< class GridPart, int codim, class Iterator >
  class FilteredGridPart< GridPartImp, FilterImp > :: IteratorWrapper
  : public Iterator
  {
    typedef IteratorWrapper< GridPart, codim, Iterator > ThisType;
    typedef Iterator BaseType;

    const GridPart *gridPart_;        
    const FilterType *filter_;
    Iterator endIter_;

  public:
    IteratorWrapper( const GridPart *gridPart,
                     const FilterType* filter,
                     const Iterator &iterator,
                     const Iterator &endIterator )
    : BaseType( iterator ),
      gridPart_( gridPart ),
      filter_( filter ),
      endIter_( endIterator )
    {
      while( (*this != endIter_) && (!filter_->has0Entity( *this )) )
        BaseType :: operator++();
    }

#if 0
    //! copy constructor 
    IteratorWrapper( const ThisType &other )
    : IteratorType( other ),
      gridPart_( other.gridPart_ ),
      filter_( other.filter_ ),
      endIter_( other.endIter_ )
    {}

    //! assignment operator
    ThisType &operator= ( const ThisType &other )
    {
      BaseType :: operator=( other );
      gridPart_ = other.gridPart_;
      filter_ = other.filter_;
      endIter_ = other.endIter_;
      return *this;
    }
#endif

    //! overloaded increment 
    ThisType &operator++ ()
    {
      do
        BaseType :: operator++();
      while( (*this != endIter_) && (!filter_->has0Entity( *this )) );
      return *this;
    }
  }; // end IteratorWrapper

}  // end namespace Dune

#endif
