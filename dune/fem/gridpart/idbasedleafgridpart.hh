#ifndef DUNE_IDBASEDLEAFGRIDPART_HH
#define DUNE_IDBASEDLEAFGRIDPART_HH

//- Dune includes 
#include <dune/fem/gridpart/gridpart.hh>

//- local includes 
#include "idbasedleafindexset.hh"

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class GridImp >
  struct IdBasedLeafGridPart;



  // Traits for IdBasedLeafGridPart
  // ------------------------------

  template< class GridImp >
  struct IdBasedLeafGridPartTraits
  {
    typedef GridImp GridType;
    typedef IdBasedLeafGridPart< GridImp > GridPartType;

    typedef IdBasedLeafIndexSet< GridImp > IndexSetType;
    static const PartitionIteratorType indexSetPartitionType = All_Partition;

    typedef typename GridType::template Codim<0>::Entity::
      LeafIntersectionIterator IntersectionIteratorType;
    
    template< int codim >
    struct Codim
    {
      template< PartitionIteratorType pitype >
      struct Partition
      {
        typedef typename GridType :: template Codim< codim >
          :: template Partition< pitype > :: LeafIterator
          IteratorType;
      };
    };

    static const bool conforming = Capabilities::isLeafwiseConforming<GridType>::v;
  };



  // IdBasedLeafGridPart
  // -------------------

  /** \class IdBasedLeafGridPart
   *  \ingroup AdaptiveLeafGP
   *  \brief index set based on the grid's local ids
   *
   *  Special implementation of the AdaptiveLeafGridPart with an underlying
   *  index set is only defined for entities with codimension 0 for use with
   *  the Dune::DiscontinuousGalerkinSpace of FiniteVolumeSpaces.
   *
   *  The underlying \ref IdBasedLeafIndexSet "index set" is a singleton for
   *  each different grid.
   *
   *  \note: This index set provides only indices for codimension 0.
   *  \note: The indices are stored in maps using the grid's local ids as the
   *         key.
   *
   *  \todo add support for higher codimensions
   */
  template< class GridImp > 
  class IdBasedLeafGridPart
  : public GridPartDefault< IdBasedLeafGridPartTraits< GridImp > > 
  {
    typedef IdBasedLeafGridPart< GridImp > ThisType;
    typedef GridPartDefault< IdBasedLeafGridPartTraits< GridImp > > BaseType;

  public:
    //- Public typedefs and enums
    //! Type definitions
    typedef IdBasedLeafGridPartTraits< GridImp > Traits;

    //! Grid implementation type
    typedef typename Traits::GridType GridType;
    //! The leaf index set of the grid implementation
    typedef typename Traits::IndexSetType IndexSetType;
    
    //! The corresponding IntersectionIterator 
    typedef typename Traits::IntersectionIteratorType IntersectionIteratorType ;
    
    //! Struct providing types of the leaf iterators on codimension cd
    template< int codim >
    struct Codim
    : public BaseType :: template Codim< codim >
    {};

  protected:
    // singleton provider 
    typedef SingletonList<const GridType* , IndexSetType > IndexSetProviderType;  

    // type of entity of codim 0
    typedef typename GridType::template Codim<0>::Entity EntityCodim0Type;

  public:
    //- Public methods
    //! Constructor
    IdBasedLeafGridPart ( GridType &grid )
    : BaseType( grid, IndexSetProviderType :: getObject( &grid ) )
    {}
    //! copy Constructor
    IdBasedLeafGridPart ( const ThisType &other )
    : BaseType( other.grid_, IndexSetProviderType :: getObject( &(other.grid()) ) )
    {}

    /** \brief Destrcutor removeing index set, if only one reference left, index set
        removed.  */
    ~IdBasedLeafGridPart() 
    { 
      IndexSetProviderType::removeObject(this->indexSet());
    }

    //! Begin iterator on the leaf level
    template< int codim >
    typename Codim< codim > :: IteratorType
    begin () const
    {
      return BaseType :: template begin< codim >();
    }

    //! Begin iterator on the leaf level
    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim > :: template Partition< pitype > :: IteratorType
    begin () const
    {
      return (*this).grid().template leafbegin< codim, pitype >();
    }

    //! Begin iterator on the leaf level
    template< int codim >
    typename Codim< codim > :: IteratorType
    end () const
    {
      return BaseType :: template end< codim >();
    }

    //! End iterator on the leaf level
    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim > :: template Partition< pitype > :: IteratorType
    end () const
    {
      return (*this).grid().template leafend< codim, pitype >();
    }

    //! ibegin of corresponding intersection iterator for given entity
    IntersectionIteratorType ibegin(const EntityCodim0Type & en) const 
    {
      return en.ileafbegin();
    }
    
    //! iend of corresponding intersection iterator for given entity
    IntersectionIteratorType iend(const EntityCodim0Type & en) const 
    {
      return en.ileafend();
    }

    //! Returns maxlevel of the grid
    int level() const { return this->grid().maxLevel(); }

    //! corresponding communication method for this grid part
    template <class DataHandleImp,class DataType>
    void communicate(CommDataHandleIF<DataHandleImp,DataType> & data, 
                     InterfaceType iftype, CommunicationDirection dir) const 
    {
      this->grid().communicate(data,iftype,dir);
    }
  };

} // end namespace Dune 

#endif
