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
    typedef typename FunctionSpaceT::Domain NormalType;
    typedef typename FunctionSpaceT::Range ValueType;

    //! Constructor
    //! The vector components to be rotated must be consecutive in the
    //! vector of unknows.
    //! \param startIdx Specifies first component of vector component
    FieldRotator(int startIdx) : 
      idx_(startIdx) {
      assert(startIdx < ValueType::size-NormalType::size+1);
    }

    //! Rotate data from basic coordinate system into normal coordinate system
    void rotateForth(const ValueType& arg, 
                     ValueType& res, 
                     const NormalType& n) const {
      rotateForth(arg, res, n, Int2Type<NormalType::size>());
    }

    //! Rotate data from normal coordinate system into basic coordinate system
    void rotateBack(const ValueType& arg, 
                    ValueType& res, 
                    const NormalType& n) const {
      rotateBack(arg, res, n, Int2Type<NormalType::size>());
    }

  private:
    // Local methods
    void rotateForth(const ValueType& arg, 
                     ValueType& res, 
                     const NormalType& n,
                     Int2Type<2>) const;
    void rotateForth(const ValueType& arg, 
                     ValueType& res, 
                     const NormalType& n,
                     Int2Type<3>) const;
    void rotateBack(const ValueType& arg, 
                    ValueType& res, 
                    const NormalType& n,
                    Int2Type<2>) const;
    void rotateBack(const ValueType& arg, 
                    ValueType& res, 
                    const NormalType& n,
                    Int2Type<3>) const;

    
    const int idx_;
    const static double eps_;
  };

} // end namespace Adi

#include "rotator.cc"

#endif
