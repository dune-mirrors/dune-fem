#ifndef DUNE_FILTEREDGRID_HH
#define DUNE_FILTEREDGRID_HH

//- System includes
#include <vector>
#include <cassert>

//- Dune includes

#include <dune/grid/common/grid.hh>
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/grid/common/sizecache.hh>

#include <dune/fem/gridpart/gridpart.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/utility/grapedataioformattypes.hh>

#include <dune/fem/gridpart/adaptiveleafindexset.hh>

namespace Dune
{

  // Forward Declarations
  // --------------------

  template <class FilterImp, class GridImp>
  struct DefaultFilterTraits;

  template< class FilterTraits >
  class FilterInterface;

  template< class FilterTraits >
  class FilterDefaultImplementation;

  //template <class GridType>
  //class TrueFilter;
  template <class GridType>
  class RadialFilter;

  template< class GridPartImp, class FilterImp, bool useFilteredIndexSet = false >
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



  // FilterInterface
  // ---------------

  /** \ingroup FilterGridPart
   *  \brief   Interface class for filter to use with a Dune::FilteredGridPart
   */
  template< class FilterTraits >
  class FilterInterface
  {
    typedef FilterInterface< FilterTraits > ThisType;

    friend class FilterDefaultImplementation< FilterTraits >;

  public:
    //! \brief Type of the filter implementation
    typedef typename FilterTraits :: FilterType FilterType;
    //! \brief type of original grid part 
    typedef typename FilterTraits :: GridPartType GridPartType;

    //! \brief type of Grid implementation
    typedef typename GridPartType :: GridType GridType;
    //! \brief type of Entity with codim=0
    typedef typename GridType :: template Codim< 0 > :: Entity EntityCodim0Type;
    //! \brief type of EntityPointer with codim=0 
    typedef typename GridType :: template Codim< 0 > :: EntityPointer EntityPointerCodim0Type;
      
    //! \brief type of Elements (i.e., Entities with codim=0)
    typedef typename GridType :: template Codim< 0 > :: Entity ElementType;
    //! \brief type of Element Pointers (i.e. Entity Pointers with codim=0)
    typedef typename GridType :: template Codim< 0 > :: EntityPointer ElementPointerType;

  private: 
    FilterInterface ()
    {}

    FilterInterface ( const ThisType & );
    ThisType &operator= ( const ThisType & );

  public:
    //! returns true if the given entity is in the domain 
    bool has0Entity ( const ElementPointerType &elementPtr ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().has0Entity( elementPtr ) );
      return asImp().has0Entity( elementPtr );
    }
    
    //! returns true if the given entity of the pointer in the domain 
    template <class EntityType>
    bool has0Entity ( const EntityType &element ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().has0Entity( element ) );
      return asImp().has0Entity( element );
    }
    
    //! returns true if an intersection is interior 
    //! (allows boundarys within a given domain)
    template< class Intersection >
    bool interiorIntersection ( const Intersection &intersection ) const
    {
      return asImp().interiorIntersection( intersection );
    }

    //! returns true if an intersection is a boundary intersection 
    template< class Intersection >
    bool intersectionBoundary( const Intersection &intersection ) const
    {
      return asImp().intersectionBoundary( intersection );
    }
    
    //! returns the boundary id for an intersection 
    template< class Intersection >
    int intersectionBoundaryId ( const Intersection &intersection ) const
    {
      return asImp().intersectionBoundaryId( intersection );
    }

    //! returns true if for an intersection a neighbor exsits 
    template< class Intersection >
    bool intersectionNeighbor ( const Intersection &intersection ) const
    {
      return asImp().intersectionNeighbor( intersection );
    }

    //! return object instance of filter for object creation  
    static FilterType createObject( const GridPartType &gridPart )
    {
      return FilterType :: createObject( gridPart );
    }

  protected:
    FilterType &asImp ()
    {
      return static_cast< FilterType & >( *this );
    }
    
    const FilterType &asImp () const
    {
      return static_cast< const FilterType & >( *this );
    }
  };



  // FilterDefaultImplementation
  // ---------------------------

  template< class FilterTraits >
  class FilterDefaultImplementation
  : public FilterInterface< FilterTraits >
  {
    typedef FilterDefaultImplementation< FilterTraits > ThisType;
    typedef FilterInterface< FilterTraits > BaseType;

  public:
    //! \brief Type of the filter implementation
    typedef typename BaseType :: FilterType FilterType;
    //! \brief type of original grid part 
    typedef typename BaseType :: GridPartType GridPartType;
    

    //! \brief type of Elements (i.e., Entities with codim=0)
    typedef typename BaseType :: ElementType ElementType;
    //! \brief type of Element Pointers (i.e. Entity Pointers with codim=0)
    typedef typename BaseType :: ElementPointerType ElementPointerType;
      
  protected:
    FilterDefaultImplementation ()
    {}

  private:
    FilterDefaultImplementation ( const ThisType & );
    ThisType &operator= ( const ThisType & );

  public:
    //! default implementation returns hasEntity0 from neighbor
    template< class Intersection >
    bool interiorIntersection( const Intersection &intersection ) const
    {
      const ElementPointerType outside = intersection.outside();
      return asImp().has0Entity( outside );
    }

    //! \brief default createObject method calling FilterType(gridPart) 
    static FilterType createObject( const GridPartType &gridPart )
    {
      return FilterType( gridPart );
    }

  private:
    bool has0Entity ( const ElementType &element ) const;
    
    bool has0Entity ( const ElementPointerType &elementPtr ) const;
    
    template< class Intersection >
    bool intersectionBoundary( const Intersection &intersection ) const;
    
    //! returns the boundary id for an intersection 
    template< class Intersection >
    int intersectionBoundaryId ( const Intersection &intersection ) const;

    template< class Intersection >
    bool intersectionNeighbor ( const Intersection &intersection ) const;

  protected:
    using BaseType :: asImp;
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
    inline bool has0Entity( const EntityPointerCodim0Type & e) const {
      return has0Entity(*e);
    }
    
    //! check whether bary center is inside of circle 
    inline bool has0Entity( const EntityCodim0Type & e) const 
    {
      const GeometryType& geom = e.geometry();
      FieldVectorType weigh(0);
      for(int i = 0; i < geom.corners(); ++i)
        weigh += geom.corner( i );
      weigh /= geom.corners();
      weigh -= center_;
      double dist = weigh.two_norm();
      return (dist <= radius_);
    }
    
    //! return what boundary id we have in case of boundary intersection 
    //! which is either it.boundary == true or has0Entity (it.ouside()) == false 
    //! so here true is a good choice 
    template <class IntersectionIteratorType>
    inline bool intersectionBoundary(const IntersectionIteratorType & it) const 
    {
      return true;
    }
    //! return what boundary id we have in case of boundary intersection 
    //! which is either it.boundary == true or has0Entity (it.ouside()) == false 
    template <class IntersectionIteratorType>
    inline int intersectionBoundaryId(const IntersectionIteratorType & it) const {
      return 1;
    }

    //! if has0Entity is true then we have an interior entity 
    template <class IntersectionIteratorType>
    inline bool intersectionNeighbor(const IntersectionIteratorType & it) const {
      return true;
    }

    static ThisType createObject(const GridPartType& gridPart) 
    {
      std::cout << "Warning, creating standard readial filter! " <<
        std::endl;
      //DUNE_THROW(NotImplemented,"Method createObject not implemented");
      FieldVectorType center(0);
      double radius = 0.5;
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
  template< class GridPartImp, class FilterImp, bool useFilteredIndexSet > 
  class FilteredGridPart
    : public GridPartImp
  {
    typedef FilteredGridPart< GridPartImp, FilterImp, useFilteredIndexSet > ThisType;

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

    template <int dummy, bool useNewIndexSet > 
    struct IndexSetSpecialization
    {
      typedef AdaptiveLeafIndexSet< GridPartImp > IndexSetType;
      static IndexSetType* create(const ThisType& gridPart) 
      {
        return new IndexSetType( gridPart );
      }

      template <class IndexSetPtr>
      static inline const IndexSetType& 
      indexSet(const GridPartImp& gridPart, const IndexSetPtr* idxSetPtr )
      {
        assert( idxSetPtr );
        return *idxSetPtr;
      }
    };

    //! when index set from gridpartimp is used return 0 
    template <int dummy>
    struct IndexSetSpecialization<dummy, false > 
    {
      typedef typename GridPartImp :: IndexSetType IndexSetType;

      static IndexSetType* create(const GridPartImp&) 
      {
        return 0;
      }
      template <class IndexSetPtr>
      static inline const IndexSetType& 
      indexSet(const GridPartImp& gridPart, const IndexSetPtr* )
      {
        return gridPart.indexSet();
      }
    };

    //! The index set use in this gridpart 
    typedef typename IndexSetSpecialization<0, useFilteredIndexSet> :: IndexSetType IndexSetType;
    
    //! The corresponding IntersectionIterator 
    typedef IntersectionIteratorWrapper<GridPartImp, IntersectionIteratorImpType> IntersectionIteratorType;

    typedef typename IntersectionIteratorType::Intersection IntersectionType;

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
      maxlevel_(0),
      indexSetPtr_( 0 )
    {
      updateStatus();
      indexSetPtr_ = IndexSetSpecialization<0, useFilteredIndexSet > :: create( *this );
    }

    //! Copy Constructor
    FilteredGridPart(const FilteredGridPart& other) :
      GridPartImp(other), 
      filter_(other.filter_),
      maxlevel_(other.maxlevel_),
      indexSetPtr_( IndexSetSpecialization<0, useFilteredIndexSet > :: create( *this ))
    {
    }

    //! constructor only taking grid 
    FilteredGridPart(GridType& grid) :
      GridPartImp(grid),   // This is only legal for LeafGridPart and HierarchicGridParts!
                           // A constructor LevelGridPart(GridType grid) exists, but will  
                           // return a LevelGridPart<GridType, maxLevel>.
      filter_( FilterType :: createObject(static_cast<const GridPartImp&> (*this) )), 
      maxlevel_(0),
      indexSetPtr_( 0 )
    {
      updateStatus();
      indexSetPtr_ = IndexSetSpecialization<0, useFilteredIndexSet > :: create( *this );
    }

    //! destructor 
    ~FilteredGridPart()
    {
      delete indexSetPtr_; 
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
      return IteratorType( this, &filter_, 
                           GridPartImp :: template begin< codim, pitype >() ,
                           GridPartImp :: template end< codim, pitype >() );
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
      return IteratorType( this, &filter_, 
                           GridPartImp :: template end< codim, pitype >() );
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

    int boundaryId ( const IntersectionType &intersection ) const
    {
      return intersection.boundaryId();
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

    //! return index set of this grid part 
    //! if IndexSetType is from GridPartImp the original index set is returned 
    const IndexSetType& indexSet() const 
    {
      return IndexSetSpecialization<0, useFilteredIndexSet > :: indexSet( *this, indexSetPtr_ );
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
    FilterType filter_;
    int maxlevel_;
    const IndexSetType* indexSetPtr_; 

  private:
    //**********************************************************************
    // IntersectionIteratorWrapper
    //**********************************************************************

    template< class GridPartType, class IteratorType >
    class IntersectionIteratorWrapper
    : public IteratorType 
    {
      // type of codim 0 entity      
      typedef typename GridPartType::GridType::template Codim<0>::EntityPointer EntityPointerCodim0Type;
      typedef typename GridPartType::GridType::template Codim<0>::Entity EntityCodim0Type;

      typedef IntersectionIteratorWrapper<GridPartType,IteratorType>  ThisType;
      
      typedef typename IteratorType::Intersection RealIntersection;

    public:
      typedef typename RealIntersection::ctype ctype;

      static const int dimension = RealIntersection::dimension;
      static const int dimensionworld = RealIntersection::dimensionworld;

      typedef FieldVector< ctype, dimensionworld > NormalVector;
      typedef FieldVector< ctype, dimension-1 > LocalVector;

      typedef typename RealIntersection::Entity Entity;
      typedef typename RealIntersection::EntityPointer EntityPointer;
      typedef typename RealIntersection::Geometry Geometry;
      typedef typename RealIntersection::LocalGeometry LocalGeometry;

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
          if ( asBase()->neighbor() ) 
          { 
            if ( filter_->interiorIntersection( *asBase() ) )
            {
              nInfo_.boundary_   = false;
              nInfo_.boundaryId_ = 0;
              nInfo_.neighbor_   = true;
            }
            else 
            {
              // otherwise get boundary information from filter 
              nInfo_.boundary_   = filter_->intersectionBoundary( *asBase() );
              nInfo_.boundaryId_ = filter_->intersectionBoundaryId( *asBase() );
              nInfo_.neighbor_   = filter_->intersectionNeighbor( *asBase() );
            }
          }
          else 
          {
            // for real boundary get boundary from filter 
            nInfo_.boundary_   = true;
            nInfo_.boundaryId_ = filter_->intersectionBoundaryId( *asBase() );
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
        bool boundary () const
        {
          return nInfo_.boundary_;
        }

        //! overloaded boundaryId method 
        int boundaryId () const
        {
          return nInfo_.boundaryId_;
        }

        //! overloaded neighbor method 
        bool neighbor () const
        {
          return nInfo_.neighbor_;
        }

        EntityPointer inside () const
        {
          return asBase()->inside();
        }

        EntityPointer outside () const
        {
          return asBase()->outside();
        }

        bool conforming () const
        {
          return asBase()->conforming();
        }

        const LocalGeometry &geometryInInside () const
        {
          return asBase()->geometryInInside();
        }

        const LocalGeometry &geometryInOutside () const
        {
          return asBase()->geometryInOutside();
        }

        const Geometry &geometry () const
        {
          return asBase()->geometry();
        }

        GeometryType type () const
        {
          return asBase()->type();
        }

        int indexInInside () const
        {
          return asBase()->indexInInside();
        }

        int indexInOutside () const
        {
          return asBase()->indexInOutside();
        }

        NormalVector outerNormal ( const LocalVector &local ) const
        {
          return asBase()->outerNormal( local );
        }

        NormalVector integrationOuterNormal ( const LocalVector &local ) const
        {
          return asBase()->integrationOuterNormal( local );
        }

        NormalVector unitOuterNormal( const LocalVector &local ) const
        {
          return asBase()->unitOuterNormal( local );
        }

        //! type of Intersection 
        typedef ThisType Intersection;
        //! dereference operator 
        inline const Intersection& operator *() const { return *this; }
        //! de-pointer operator 
        inline const Intersection* operator ->() const { return this; }
      protected:
        //! return reference to base class 
        inline IteratorType & asBase() { return static_cast<IteratorType &>(*this); }
        //! return reference to base class 
        inline const IteratorType & asBase() const { return static_cast<const IteratorType &>(*this); }
        
        const GridPartType* gridPart_;        
        const FilterType* filter_;
        const IteratorType endIter_;        
    }; // end IntersectionIteratorWrapper
  }; // end FilteredGridPart

  // FilteredGridPart :: IteratorWrapper
  // -----------------------------------

  template< class GridPartImp, class FilterImp, bool useFilteredIndexSet >
  template< class GridPart, int codim, class Iterator >
  class FilteredGridPart< GridPartImp, FilterImp, useFilteredIndexSet > :: IteratorWrapper
  : public Iterator
  {
    typedef IteratorWrapper< GridPart, codim, Iterator > ThisType;
    typedef Iterator BaseType;

    const GridPart *gridPart_;        
    const FilterType *filter_;
    Iterator endIter_;

  public:
    //! constructor creating begin iterator 
    IteratorWrapper( const GridPart *gridPart,
                     const FilterType* filter,
                     const Iterator &iterator,
                     const Iterator &endIterator )
    : BaseType( iterator ),
      gridPart_( gridPart ),
      filter_( filter ),
      endIter_( endIterator )
    {
      typedef typename BaseType :: Entity Entity ;
      while( (*this != endIter_) && (!filter_->has0Entity(
              static_cast<const Entity&> (*(*this)) )) )
        BaseType :: operator++();
    }

    //! constructor creating end iterator 
    IteratorWrapper( const GridPart *gridPart,
                     const FilterType* filter,
                     const Iterator &endIterator )
    : BaseType( endIterator ),
      gridPart_( gridPart ),
      filter_( filter ),
      endIter_( endIterator )
    {
    }

    //! overloaded increment 
    ThisType &operator++ ()
    {
      typedef typename BaseType :: Entity Entity ;
      do
        BaseType :: operator++();
      while( (*this != endIter_) && (!filter_->has0Entity( static_cast<const Entity&> (*(*this) ))) );
      return *this;
    }
  }; // end IteratorWrapper

}  // end namespace Dune

#endif
