#include "rotator.hh"

namespace Adi {
  
  template <class FunctionSpaceT>
  const double FieldRotator<FunctionSpaceT>::eps_ = 1.0e-14;

  template <class FunctionSpaceT>
  inline void FieldRotator<FunctionSpaceT>::
  rotateForth(ValueType& res, 
              const NormalType& n,
              Int2Type<1>) const {
    // res = arg ;
    res[idx_] = res[idx_] * n[0];
  }

  template <class FunctionSpaceT>
  inline void FieldRotator<FunctionSpaceT>::
  rotateForth(ValueType& res, 
              const NormalType& n,
              Int2Type<2>) const {
    // res = arg;
    double a[2] = {res[idx_],res[idx_+1]};
    res[idx_] = n[0]*a[0] + n[1]*a[1];
    res[idx_+1] = -n[1]*a[0] + n[0]*a[1];
  }
  
  template <class FunctionSpaceT>
  inline void FieldRotator<FunctionSpaceT>::
  rotateForth(ValueType& res, 
              const NormalType& n,
              Int2Type<3>) const {
    // res = arg;

    double a[3]={res[idx_],res[idx_+1],res[idx_+2]};
 
    double d = std::sqrt(n[0]*n[0]+n[1]*n[1]);

    if (d > 1.0e-8) {
      double d_1 = 1.0/d;
      res[idx_]   =   n[0] * a[0]
                    + n[1] * a[1]
                    + n[2] * a[2];
      res[idx_+1] = - n[1] * d_1 * a[0]
                    + n[0] * d_1 * a[1];
      res[idx_+2] = - n[0] * n[2] * d_1 * a[0]
                    - n[1] * n[2] * d_1 * a[1]
                    + d                 * a[2];
    } else {
      res[idx_]   =   n[2] * a[2];
      res[idx_+1] =          a[1];
      res[idx_+2] = - n[2] * a[0];
    }
  }

  template <class FunctionSpaceT>
  inline void FieldRotator<FunctionSpaceT>::
  rotateBack(ValueType& res, 
	     const NormalType& n,
	     Int2Type<1>) const {
    // res = arg;

    res[idx_] = res[idx_] * n[0]; 
  }

  template <class FunctionSpaceT>
  inline void FieldRotator<FunctionSpaceT>::
  rotateBack(ValueType& res, 
	     const NormalType& n,
	     Int2Type<2>) const {
    // res = arg;

    double a[2] = {res[idx_],res[idx_+1]};
    res[idx_] = n[0]*a[0] - n[1]*a[1];
    res[idx_+1] = n[1]*a[0] + n[0]*a[1];
  }

  template <class FunctionSpaceT>
  inline void FieldRotator<FunctionSpaceT>::
  rotateBack(ValueType& res, 
	     const NormalType& n,
	     Int2Type<3>) const {
    // res = arg;
    double a[3]={res[idx_],res[idx_+1],res[idx_+2]};

    double d = std::sqrt(n[0]*n[0]+n[1]*n[1]);

    if (d > 1.0e-8) {
      double d_1 = 1.0/d;
      res[idx_]   =   n[0]              * a[0]
                    - n[1] * d_1        * a[1]
                    - n[0] * n[2] * d_1 * a[2];
      res[idx_+1] =   n[1]              * a[0]
                    + n[0] * d_1        * a[1]
                    - n[1] * n[2] * d_1 * a[2];
      res[idx_+2] =   n[2]              * a[0]
                    + d                 * a[2];
    } else {
      res[idx_]   = - n[2] * a[2];
      res[idx_+1] =          a[1];
      res[idx_+2] =   n[2] * a[0];
    }
    
  }

} // end namespace Adi
