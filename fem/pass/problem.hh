#ifndef DUNE_PROBLEM_HH
#define DUNE_PROBLEM_HH

#include <dune/common/fvector.hh>

namespace Dune {

  /**
   * @brief Interface class for problem definition in the LDG context.
   *
   * The problem interface prescribes the methods needed by the implementation
   * of local DG passes. Users need to derive from either this class or the
   * ProblemDefault class.
   */
  template <class ProblemTraits>
  class ProblemInterface 
  {
  public:
    typedef ProblemTraits Traits;
    typedef typename Traits::ProblemType ProblemType;
    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef typename Traits::GridType GridType;
    typedef typename GridType::Traits::IntersectionIterator IntersectionIterator;
    typedef typename GridType::template Codim<0>::Entity EntityType;

    //typedef FieldVector<typename DiscreteFunctionSpaceType::DomainFieldType,
    //                    DiscreteFunctionSpaceType::DimDomain-1> FaceDomainType;

  public:
    //! Returns true if problem has a flux contribution.
    //! If you program this method to return true, make sure to implement
    //! numericalFlux and analyticalFlux as well.
    bool hasFlux() const { return asImp().hasFlux(); }

    //! Returns true if problem has a source term.
    //! If you program this method to return true, make sure to implement
    //! method source well.
    bool hasSource() const { return asImp().hasSource(); }

    //! Computes the numerical flux at a cell interface.
    //! ... More to come ...
    template <class ArgumentTuple, class FaceDomainType>
    double numericalFlux(IntersectionIterator& it,
                         double time, const FaceDomainType& x,
                         const ArgumentTuple& uLeft, 
                         const ArgumentTuple& uRight,
                         RangeType& gLeft,
                         RangeType& gRight)
    { return asImp().numericalFlux(); }

    template <class ArgumentTuple, class FaceDomainType>
    double boundaryFlux(IntersectionIterator& it,
                        double time, const FaceDomainType& x,
                        const ArgumentTuple& uLeft,
                        RangeType& gLeft)
    { return asImp().boundaryFlux(); }

    //! Computes the analytical flux of the problem.
    //! ... More to come ...
    template <class ArgumentTuple>
    void analyticalFlux(EntityType& en,
                        double time, const DomainType& x,
                        const ArgumentTuple& u, JacobianRangeType& f) 
    { asImp().analyticalFlux(); }

    //! Implements the source term of the problem.
    template <class ArgumentTuple, class JacobianTuple>
    void source(EntityType& en, 
                double time, const DomainType& x,
                const ArgumentTuple& u, 
                const JacobianTuple& jac, 
                RangeType& s) 
    { asImp().source(); }

    // Need second source function with grad u as additional argument 

  private:
    ProblemType& asImp() { return static_cast<ProblemType&>(*this); }
    const ProblemType& asImp() const { return static_cast<const ProblemType&>(*this); }
  };
  
  //! Default implementation of the ProblemInterface where methods for 
  //! the fluxes and the source term do nothing, so that the user needn't
  //! implement them if not needed.
  template <class ProblemTraits>
  class ProblemDefault : 
    public ProblemInterface<ProblemTraits> {
  public:
    typedef ProblemTraits Traits;
    typedef typename Traits::ProblemType ProblemType;
    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef typename Traits::GridType GridType;
    typedef typename GridType::template Codim<0>::Entity EntityType;
    typedef typename GridType::Traits::IntersectionIterator IntersectionIterator;

    //typedef FieldVector<typename DiscreteFunctionSpaceType::DomainFieldType,
    //                    DiscreteFunctionSpaceType::DimDomain-1> FaceDomainType;
  public:
    //! Empty implementation that fails if problem claims to have a flux
    //! contribution.
    template <class ArgumentTuple, class FaceDomainType>
    double numericalFlux(IntersectionIterator& it,
                         double time, const FaceDomainType& x,
                         const ArgumentTuple& uLeft, 
                         const ArgumentTuple& uRight,
                         RangeType& gLeft,
                         RangeType& gRight)
    { 
      assert(!this->hasFlux()); 
      gLeft = 0.0;
      gRight = 0.0;
      return 0.0;
    }

    template <class ArgumentTuple, class FaceDomainType>
    double boundaryFlux(IntersectionIterator& it,
                        double time, const FaceDomainType& x,
                        const ArgumentTuple& uLeft,
                        RangeType& gLeft)
    {
      assert(!this->hasFlux());
      gLeft = 0.0;
      return 0.0;
    }

    //! Empty implementation that fails if problem claims to have a flux
    //! contribution.
    template <class ArgumentTuple>
    void analyticalFlux(EntityType& en,
                        double time, const DomainType& x,
                        const ArgumentTuple& u, JacobianRangeType& f)
    { 
      assert(!this->hasFlux()); 
      f = 0.0;
    }

    //! Empty implementation that fails if problem claims to have a source 
    //! term.
    template <class ArgumentTuple, class JacobianTuple>
    void source(EntityType& en, 
                double time, const DomainType& x,
                const ArgumentTuple& u, const JacobianTuple& jac, 
                RangeType& s)
    { 
      assert(!this->hasSource()); 
      s = 0.0;
    }
  };

}  // end namespace Dune
#endif
