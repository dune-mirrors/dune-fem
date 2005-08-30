#ifndef DUNE_FINITEVOLUME_HH
#define DUNE_FINITEVOLUME_HH

// Include Dune headers
#include <dune/common/fvector.hh>

namespace Dune {

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
                       const NormalType& position,
                       const NormalType& normal,
                       ValueType& result) const {
      return asImp(valsEntity, valsNeigh, position, normal, result);
    }

    //! Evaluation of corresponding analytical flux function
    void analyticalFlux(const ValueType& vals,
                        const NormalType& position,
                        JacobianType& result) const {
      asImp(vals, position, result);
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
  class ConstantLinearAdvectionUpwindFlux : 
    public NumericalFlux<FunctionSpaceImp,
                         ConstantLinearAdvectionUpwindFlux<FunctionSpaceImp> >
  {
    typedef NumericalFlux<FunctionSpaceImp,
                          ConstantLinearAdvectionUpwindFlux<FunctionSpaceImp> > BaseType;
  public:
    typedef typename BaseType::NormalType NormalType;
    typedef typename BaseType::ValueType ValueType;
    typedef typename BaseType::JacobianType JacobianType;

    ConstantLinearAdvectionUpwindFlux(const NormalType& velocity) :
      vel_(velocity) {}

    ~ConstantLinearAdvectionUpwindFlux() {}

    double operator() (const ValueType& valsEntity,
                       const ValueType& valsNeigh,
                       const NormalType& position,
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
                        const NormalType& position,
                        JacobianType& result) const {
      for (int i = 0; i < ValueType::size; ++i) {
        result[i] = vel_;
        result[i] *= vals[i];
     }
    }

  private:
    const NormalType vel_;
  }; // end class ConstantLinearAdvectionUpwindFlux

  template <class FunctionSpaceImp, class VelocityImp>
  class LinearAdvectionUpwindFlux : 
    public NumericalFlux<FunctionSpaceImp, LinearAdvectionUpwindFlux<FunctionSpaceImp, VelocityImp> >
  {
    typedef NumericalFlux<FunctionSpaceImp,
                          LinearAdvectionUpwindFlux<FunctionSpaceImp, VelocityImp> > BaseType;
  public:
    typedef typename BaseType::NormalType NormalType;
    typedef typename BaseType::ValueType ValueType;
    typedef typename BaseType::JacobianType JacobianType;

    LinearAdvectionUpwindFlux(const VelocityImp& velo) :
      vel_(velo) {}

    ~LinearAdvectionUpwindFlux() {}
    
    double operator() (const ValueType& valsEntity,
                       const ValueType& valsNeigh,
                       const NormalType& position,
                       const NormalType& normal,
                       ValueType& result) const 
    {
      normalizeN(normal);
      const double project = vel_(position)*this->stdnorm_;
      result = 0.0;
      if (project > 0) {
        result += valsEntity;
      }
      else if (project < 0) {
        result += valsNeigh;
      }
      result *= project;
      return project;
    }

    void analyticalFlux(const ValueType& vals,
                        const NormalType& position,
                        JacobianType& result) const {
      for (int i = 0; i < ValueType::size; ++i) {
        result[i] = vel_(position);
        result[i] *= vals[i];
      }
    }
    
  public:
    VelocityImp vel_;
  }; // end class LinearUpwindFlux

} // end namespace Dune

#endif
