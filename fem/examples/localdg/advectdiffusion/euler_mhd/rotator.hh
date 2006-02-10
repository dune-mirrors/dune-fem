#ifndef ADI_ROTATOR_HH
#define ADI_ROTATOR_HH

//- Dune includes
#include <dune/common/misc.hh>

//- Local includes

using namespace Dune;

namespace Adi {

  template <class FunctionSpaceT>
  class FieldRotator {
  public:
    //- Global typedefs
    typedef typename FunctionSpaceT::DomainType NormalType;
    typedef typename FunctionSpaceT::RangeType ValueType;

    //! Constructor
    //! The vector components to be rotated must be consecutive in the
    //! vector of unknows.
    //! \param startIdx Specifies first component of vector component
    FieldRotator(int startIdx) : 
      idx_(startIdx) {
      assert(startIdx < ValueType::size-NormalType::size+1);
    }

    //! Rotate data from basic coordinate system into normal coordinate system
    void rotateForth(ValueType& res, 
                     const NormalType& n) const {
      rotateForth(res, n, Int2Type<NormalType::size>());
    }

    //! Rotate data from normal coordinate system into basic coordinate system
    void rotateBack(ValueType& res, 
                    const NormalType& n) const {
      rotateBack(res, n, Int2Type<NormalType::size>());
    }

  private:
    // Local methods
    void rotateForth(ValueType& res, 
                     const NormalType& n,
                     Int2Type<1>) const;
    void rotateForth(ValueType& res, 
                     const NormalType& n,
                     Int2Type<2>) const;
    void rotateForth(ValueType& res, 
                     const NormalType& n,
                     Int2Type<3>) const;
    void rotateBack(ValueType& res, 
                    const NormalType& n,
                    Int2Type<1>) const;
    void rotateBack(ValueType& res, 
                    const NormalType& n,
                    Int2Type<2>) const;
    void rotateBack(ValueType& res, 
                    const NormalType& n,
                    Int2Type<3>) const;

    
    const int idx_;
    const static double eps_;
  };

} // end namespace Adi

#include "rotator.cc"

#endif
