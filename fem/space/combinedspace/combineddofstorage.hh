#ifndef DUNE_COMBINEDDOFSTORAGE_HH
#define DUNE_COMBINEDDOFSTORAGE_HH

//- local includes 
#include <dune/fem/space/common/dofstorage.hh>

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

} // end namespace Dune

#endif
