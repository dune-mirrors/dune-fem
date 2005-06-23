#ifndef ADI_FINITE_VOLUME_HH
#define ADI_FINITE_VOLUME_HH

// Include Dune headers
#include <dune/common/fvector.hh>

using namespace Dune;

namespace Adi {

  //! Base class describing interface of numerical flux function
  template <
    typename FunctionSpaceImp,
    typename NumericalFluxImp
    >
  class NumericalFlux {
  public:
    typedef FieldVector<double, FunctionSpaceImp::DimDomain> NormalType;
    typedef FieldVector<double, FunctionSpaceImp::DimRange> ValueType;
    typedef FieldMatrix<double, 
                        FunctionSpaceImp::DimRange,
                        FunctionSpaceImp::DimDomain> JacobianType;

    NumericalFlux() :
    stdnorm_(0.0) {}

    virtual ~NumericalFlux() {}

    //! Evaluation of numerical flux function
    double operator() (const ValueType& valsEntity,
                       const ValueType& valsNeigh,
                       const NormalType& normal,
                       ValueType& result) const {
      return asImp(valsEntity, valsNeigh, normal, result);
    }

    //! Evaluation of corresponding analytical flux function
    void analyticalFlux(const ValueType& vals,
                        JacobianType& result) const {
      asImp(vals, result);
    }

  protected:
    void normalizeN(const NormalType& normal) const {
      stdnorm_ = normal;
      stdnorm_ /= normal.two_norm();
    }

    mutable NormalType stdnorm_;
  private:
    NumericalFluxImp& asImp() { 
      return static_cast<NumericalFluxImp&>(*this); 
    }

    const NumericalFluxImp& asImp() const {
      return static_cast<const NumericalFluxImp&>(*this);
    }

  }; // end class NumericalFlux



  //! Implementation of a linear advection flux function using upwind flux
  template <class FunctionSpaceImp>
  class LinearAdvectionUpwindFlux : 
    public NumericalFlux<FunctionSpaceImp,
                         LinearAdvectionUpwindFlux<FunctionSpaceImp> >
  {
    typedef NumericalFlux<FunctionSpaceImp,
                          LinearAdvectionUpwindFlux<FunctionSpaceImp> > BaseType;
  public:
    typedef typename BaseType::NormalType NormalType;
    typedef typename BaseType::ValueType ValueType;
    typedef typename BaseType::JacobianType JacobianType;

    LinearAdvectionUpwindFlux(const NormalType& velocity) :
      vel_(velocity) {}

    ~LinearAdvectionUpwindFlux() {}

    double operator() (const ValueType& valsEntity,
                       const ValueType& valsNeigh,
                       const NormalType& normal,
                       ValueType& result) const
    {
      normalizeN(normal);
      const double project = vel_*this->stdnorm_;
      result = 0.0;
      result += (project > 0) ? valsEntity : valsNeigh;
      result *= project;
      return project;
    }

    void analyticalFlux(const ValueType& vals,
                        JacobianType& result) const {
      for (int i = 0; i < ValueType::size; ++i) {
        result[i] = vel_;
        result[i] *= vals[i];
     }
    }

  private:
    const NormalType vel_;
  }; // end class LinearAdvectionUpwindFlux

} // end namespace Adi

#endif
