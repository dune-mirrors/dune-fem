#ifndef DUNE_COMBINEDDOFSTORAGE_HH
#define DUNE_COMBINEDDOFSTORAGE_HH

//- local includes 
#include <dune/fem/space/common/dofstorage.hh>
#include <dune/fem/storage/subarray.hh>

namespace Dune
{

  //! Utility class that helps in the transformation between dofs in the
  //! combined space and its enclosed spaces
  template< class ContainedMapper, DofStoragePolicy policy >
  class CombinedDofConversionUtility;


  //! does the same as DofConversionUtility<PointBased>, just other
  //! construtor 
  template< class ContainedMapper >
  class CombinedDofConversionUtility< ContainedMapper, PointBased >
  : public DofConversionUtility< PointBased >
  {
  public:
    typedef ContainedMapper ContainedMapperType;

  private:
    typedef DofConversionUtility< PointBased > BaseType;

  public:
    inline CombinedDofConversionUtility ( const ContainedMapperType &,
                                          int numComponents )
    : BaseType( numComponents )
    {}
  };

  //! Specialisation for VariableBased approach
  template< class ContainedMapper >
  class CombinedDofConversionUtility< ContainedMapper, VariableBased >
  {
  public:
    typedef ContainedMapper ContainedMapperType;

  protected:
    const ContainedMapperType &mapper_;

  public:
    /** \brief constructor
     *  
     *  \param[in]  mapper  mapper of the contained space
     *  \param[in]  size    number of global DoFs per component
     */
    inline CombinedDofConversionUtility ( const ContainedMapperType &mapper,
                                          int size )
    : mapper_( mapper )
    {}

    //! Find out what type of policy this is.
    inline static DofStoragePolicy policy ()
    {
      return VariableBased;
    }

    //! Set new size after adaptation.
    inline void newSize ( int size )
    {}

    //! Component which the actual base function index gives a contribution
    //! \return is in range {0, dimRange-1}
    int component ( int combinedIndex ) const
    { 
      return combinedIndex / containedSize();
    }

    //! Number of the (scalar) base function belonging to base function index
    int containedDof ( int combinedIndex ) const
    {
      return combinedIndex % containedSize();
    }

    //! Reverse operation of containedDof, component
    //! i == combinedDof(containedDof(i), component(i))
    int combinedDof ( int containedIndex, int component ) const
    {
      return containedIndex + (component * containedSize());
    }

  protected:
    inline int containedSize () const
    {
      return mapper_.size();
    }
  };


  
  template< class CombinedSpace >
  class CombinedSubMapper
  : public IndexMapperInterface< CombinedSubMapper< CombinedSpace > >
  {
  public:
    //- Typedefs and enums
    typedef CombinedSpace CombinedSpaceType;
    
  private:
    typedef CombinedSubMapper< CombinedSpaceType > ThisType;

  public:
    typedef typename CombinedSpaceType :: ContainedDiscreteFunctionSpaceType
      ContainedDiscreteFunctionSpaceType;
    typedef typename ContainedDiscreteFunctionSpaceType :: MapperType
      ContainedMapperType;
    typedef typename CombinedSpaceType :: DofConversionType DofConversionType;

  public:
    //- Public methods
    CombinedSubMapper ( const CombinedSpaceType &spc,
                        unsigned int component )
    : mapper_( spc.containedSpace().mapper() ),
      component_( component ),
      utilGlobal_( mapper_,
                  spc.myPolicy() == PointBased ? 
                  spc.numComponents() :
                  spc.size()/spc.numComponents())
    {
      assert(component_<CombinedSpaceType::dimRange);
    }

    CombinedSubMapper(const ThisType& other) :
      mapper_(other.mapper_),
      component_(other.component_),
      utilGlobal_(other.utilGlobal_) { 
      assert(component_<CombinedSpaceType::dimRange);
    }

    //! Total number of degrees of freedom
    inline unsigned int size () const {
      return mapper_.size();
    }
    inline unsigned int range () const {
      return size()*CombinedSpaceType::dimRange;
    }
    inline const unsigned int operator[] ( unsigned int index ) const {
      utilGlobal_.newSize(mapper_.size());
      return utilGlobal_.combinedDof(index, component_);
    }

  private:
    ThisType& operator=(const ThisType& other);
    //- Data members
    const ContainedMapperType& mapper_;
    unsigned int component_;
    mutable DofConversionType utilGlobal_;
  };

} // end namespace Dune

#endif
