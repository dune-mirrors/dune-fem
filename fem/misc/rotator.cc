#include "rotator.hh"

namespace Adi {
  
  template <class FunctionSpaceT>
  const double FieldRotator<FunctionSpaceT>::eps_ = 1.0e-14;

  template <class FunctionSpaceT>
  inline void FieldRotator<FunctionSpaceT>::
  rotateForth(const ValueType& arg, 
              ValueType& res, 
              const NormalType& n,
              Int2Type<2>) const {
    res = arg;
    res[idx_] = n[0]*arg[idx_] + n[1]*arg[idx_+1];
    res[idx_+1] = -n[1]*arg[idx_] + n[0]*arg[idx_+1];
  }
  
  template <class FunctionSpaceT>
  inline void FieldRotator<FunctionSpaceT>::
  rotateForth(const ValueType& arg, 
              ValueType& res, 
              const NormalType& n,
              Int2Type<3>) const {
    res = arg;
    
    double d = std::sqrt(n[0]*n[0]+n[1]*n[1]);

    if (d > 1.0e-8) {
      double d_1 = 1.0/d;
      res[idx_]   =   n[0] * arg[idx_]
                    + n[1] * arg[idx_+1]
                    + n[2] * arg[idx_+2];
      res[idx_+1] = - n[1] * d_1 * arg[idx_]
                    + n[0] * d_1 * arg[idx_+1];
      res[idx_+2] = - n[0] * n[2] * d_1 * arg[idx_]
                    - n[1] * n[2] * d_1 * arg[idx_+1]
                    + d                 * arg[idx_+2];
    } else {
      res[idx_]   =   n[2] * arg[idx_+2];
      res[idx_+1] =          arg[idx_+1];
      res[idx_+2] = - n[2] * arg[idx_];
    }
  }

  template <class FunctionSpaceT>
  inline void FieldRotator<FunctionSpaceT>::
  rotateBack(const ValueType& arg, 
              ValueType& res, 
              const NormalType& n,
              Int2Type<2>) const {
    res = arg;

    res[idx_] = n[0]*arg[idx_] - n[1]*arg[idx_+1];
    res[idx_+1] = n[1]*arg[idx_] + n[0]*arg[idx_+1];
  }

  template <class FunctionSpaceT>
  inline void FieldRotator<FunctionSpaceT>::
  rotateBack(const ValueType& arg, 
              ValueType& res, 
              const NormalType& n,
              Int2Type<3>) const {
    res = arg;

    double d = std::sqrt(n[0]*n[0]+n[1]*n[1]);

    if (d > 1.0e-8) {
      double d_1 = 1.0/d;
      res[idx_]   =   n[0]              * arg[idx_]
                    - n[1] * d_1        * arg[idx_+1]
                    - n[0] * n[2] * d_1 * arg[idx_+2];
      res[idx_+1] =   n[1]              * arg[idx_]
                    + n[0] * d_1        * arg[idx_+1]
                    - n[1] * n[2] * d_1 * arg[idx_+2];
      res[idx_+2] =   n[2]              * arg[idx_]
                    + d                 * arg[idx_+2];
    } else {
      res[idx_]   = - n[2] * arg[idx_+2];
      res[idx_+1] =          arg[idx_+1];
      res[idx_+2] =   n[2] * arg[idx_];
    }
    
  }

} // end namespace Adi
