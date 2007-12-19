#ifndef DUNE_COMBINEDDOFSTORAGE_HH
#define DUNE_COMBINEDDOFSTORAGE_HH

//- local includes 
#include <dune/fem/space/common/dofstorage.hh>
#include <dune/fem/storage/subarray.hh>

namespace Dune {

  //! Utility class that helps in the transformation between dofs in the
  //! combined space and its enclosed spaces
  template <class DiscreteFunctionSpaceImp, DofStoragePolicy p>
  class CombinedDofConversionUtility; 


  //! does the same as DofConversionUtility<PointBased>, just other
  //! construtor 
  template <class DiscreteFunctionSpaceImp>
  class CombinedDofConversionUtility<DiscreteFunctionSpaceImp,PointBased> 
  : public DofConversionUtility<PointBased> 
  {
  public:
    CombinedDofConversionUtility(const DiscreteFunctionSpaceImp & , 
                                 int numComponents) : 
      DofConversionUtility<PointBased>(numComponents) {}
  };

  //! Specialisation for VariableBased approach
  template <class DiscreteFunctionSpaceImp>
  class CombinedDofConversionUtility<DiscreteFunctionSpaceImp,VariableBased> {
  public:
    //! Constructor
    //! \param spc the discrete function space
    //! \param size Number of global dofs per component.
    CombinedDofConversionUtility(const DiscreteFunctionSpaceImp & spc, int size) :
      spc_(spc)
    {}

    //! Find out what type of policy this is.
    static DofStoragePolicy policy() { return VariableBased; }

    //! Set new size after adaptation.
    void newSize(int size) {}

    //! Component which the actual base function index gives a contribution
    //! \return is in range {0, dimRange-1}
    int component(int combinedIndex) const { 
      return combinedIndex / orgSize(); 
    }

    //! Number of the (scalar) base function belonging to base function index
    int containedDof(int combinedIndex) const {
      return combinedIndex % orgSize();
    }

    //! Reverse operation of containedDof, component
    //! i == combinedDof(containedDof(i), component(i))
    int combinedDof(int containedIndex, int component) const {
      return containedIndex + (component * orgSize());
    }

  private:
    int orgSize () const { return spc_.size(); }
    const DiscreteFunctionSpaceImp & spc_; 
  };


  template <class CombinedSpaceImp>
  class CombinedSubMapper : 
    public IndexMapperInterface<CombinedSubMapper< CombinedSpaceImp > >
  {
  public:
    //- Typedefs and enums
    typedef CombinedSpaceImp CombinedSpaceType;
    
    typedef CombinedSubMapper<CombinedSpaceType> ThisType;

    typedef typename CombinedSpaceType::ContainedDiscreteFunctionSpaceType 
    ContainedDiscreteFunctionSpaceType; 
    typedef typename ContainedDiscreteFunctionSpaceType::MapperType 
    ContainedMapperType;
    typedef typename CombinedSpaceType::DofConversionType DofConversionType;

  public:
    //- Public methods
    CombinedSubMapper(const CombinedSpaceType& spc,
              unsigned int component) :
      mapper_(spc.containedSpace().mapper()),
      component_(component),
      utilGlobal_(spc.containedSpace(),
                  spc.myPolicy() == PointBased ? 
                  spc.numComponents() :
                  spc.size()/spc.numComponents())
    {
      assert(component_<CombinedSpaceType::DimRange);
    }
    CombinedSubMapper(const ThisType& other) :
      mapper_(other.mapper_),
      component_(other.component_),
      utilGlobal_(other.utilGlobal_) { 
      assert(component_<CombinedSpaceType::DimRange);
    }

    //! Total number of degrees of freedom
    inline unsigned int size () const {
      return mapper_.size();
    }
    inline unsigned int range () const {
      return size()*CombinedSpaceType::DimRange;
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
